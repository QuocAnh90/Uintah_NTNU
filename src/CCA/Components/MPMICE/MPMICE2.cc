/*
 * The MIT License
 *
 * Copyright (c) 1997-2020 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

 // MPMICE2.cc
#include <CCA/Components/MPMICE/MPMICE2.h>
#include <CCA/Components/MPMICE/Core/MPMICELabel.h>

#include <CCA/Components/ICE/AMRICE.h>
#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/ICE/CustomBCs/BoundaryCond.h>
#include <CCA/Components/ICE/EOS/EquationOfState.h>
#include <CCA/Components/ICE/ICE.h>
#include <CCA/Components/MPM/Materials/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/RigidMPM.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/ShellMPM.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>

#include <CCA/Ports/Output.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SolverInterface.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/AMR_CoarsenRefine.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/Utils.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Math/MiscMath.h>

#include <cfloat>
#include <cstdio>
#include <Core/Util/DebugStream.h>

#include <iomanip>
#include <errno.h>
#include <fenv.h>


using namespace Uintah;
using namespace std;
//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "MPMICE_NORMAL_COUT:+,MPMICE_DOING_COUT".....
//  MPMICE_NORMAL_COUT:  dumps out during problemSetup 
//  MPMICE_DOING_COUT:   dumps when tasks are scheduled and performed
//  default is OFF


static DebugStream cout_norm("MPMICE2_NORMAL_COUT", false);
static DebugStream cout_doing("MPMICE2_DOING_COUT", false);
static DebugStream ds_EqPress("DBG_EqPress", false);

MPMICE2::MPMICE2(const ProcessorGroup* myworld,
    const MaterialManagerP materialManager,
    MPMType2 mpmtype, const bool doAMR)
    : ApplicationCommon(myworld, materialManager)
{
    MIlb = scinew MPMICELabel();

    d_rigidMPM = false;
    d_testForNegTemps_mpm = true;

    switch (mpmtype) {

    case RIGID_MPMICE2:
        d_mpm = scinew RigidMPM(myworld, m_materialManager);
        d_rigidMPM = true;
        break;
    case SHELL_MPMICE2:
        d_mpm = scinew ShellMPM(myworld, m_materialManager);
        break;

    default:
        d_mpm = scinew SerialMPM(myworld, m_materialManager);
    }

    // Don't do AMRICE with MPMICE for now...
    if (doAMR) {
        d_ice = scinew AMRICE(myworld, m_materialManager);
    }
    else {
        d_ice = scinew ICE(myworld, m_materialManager);
    }

    d_exchModel = d_ice->d_exchModel;

    Ilb = d_ice->lb;
    Mlb = d_mpm->lb;

    d_SMALL_NUM = d_ice->d_SMALL_NUM;
    d_TINY_RHO = 1.e-12;  // Note, within MPMICE2, d_TINY_RHO is only applied
                           // to MPM materials, for which its value is hardcoded,
                           // unlike the situation for ice materials

    d_switchCriteria = 0;
}
//______________________________________________________________________
//
MPMICE2::~MPMICE2()
{
    d_mpm->releaseComponents();
    d_ice->releaseComponents();

    delete MIlb;
    delete d_mpm;
    delete d_ice;

    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->releaseComponents();
            delete am;
        }
    }
}

//__________________________________
//    For recomputing timesteps
double MPMICE2::recomputeDelT(const double delT)
{
    return delT / 2.0;
}
//______________________________________________________________________
//
void MPMICE2::problemSetup(const ProblemSpecP& prob_spec,
    const ProblemSpecP& restart_prob_spec,
    GridP& grid)
{
    cout_doing << "Doing MPMICE2::problemSetup " << endl;

    ProblemSpecP phys_cons_ps = prob_spec->findBlock("PhysicalConstants");
    if (phys_cons_ps) {
        phys_cons_ps->require("reference_pressure", d_ref_press);
        phys_cons_ps->require("gravity", d_gravity);
    }
    else {
        throw ProblemSetupException(
            "\n Could not find the <PhysicalConstants> section in the input file.  This section contains <gravity> and <reference pressure> \n"
            " This pressure is used during the problem intialization and when\n"
            " the pressure gradient is interpolated to the MPM particles \n"
            " you must have it for all MPMICE and multimaterial ICE problems\n",
            __FILE__, __LINE__);
    }

    //__________________________________
    //  M P M
    d_mpm->setComponents(this);
    dynamic_cast<ApplicationCommon*>(d_mpm)->problemSetup(prob_spec);

    d_mpm->setWithICE();
    d_mpm->setWithMPMICE2();
    d_ice->setWithMPMICE2();
    d_mpm->problemSetup(prob_spec, restart_prob_spec, grid);

    d_8or27 = d_mpm->flags->d_8or27;
    if (d_8or27 == 8) {
        NGN = 1;
    }
    else {
        NGN = 2;
    }

    Ghost::GhostType gp;
    int ngc_p;
    d_mpm->getParticleGhostLayer(gp, ngc_p);

    //__________________________________
    //  I C E
    d_ice->setComponents(this);
    dynamic_cast<ApplicationCommon*>(d_ice)->problemSetup(prob_spec);

    d_ice->setWithMPM();

    // Communicate the particle ghost from MPM to ICE. Used only in the
    // HEChem/Unsteady_Burn model.
    d_ice->setParticleGhostLayer(gp, ngc_p);

    if (d_rigidMPM) {
        d_ice->setWithRigidMPM();
    }

    d_switchCriteria = dynamic_cast<SwitchingCriteria*>
        (getPort("switch_criteria"));

    if (d_switchCriteria) {
        d_switchCriteria->problemSetup(prob_spec, restart_prob_spec, m_materialManager);
    }

    d_ice->problemSetup(prob_spec, restart_prob_spec, grid);

    ProblemSpecP mpm_ps = 0;
    mpm_ps = prob_spec->findBlock("MPM");

    if (!mpm_ps) {
        mpm_ps = restart_prob_spec->findBlock("MPM");
    }
    mpm_ps->get("testForNegTemps_mpm", d_testForNegTemps_mpm);

    //__________________________________
    //  bulletproofing
    if (isAMR() && !isLockstepAMR()) {
        ostringstream msg;
        msg << "\n ERROR: You must add \n"
            << " <useLockStep> true </useLockStep> \n"
            << " inside of the <AMR> section for MPMICE2 and AMR. \n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }

    if (cout_norm.active()) {
        cout_norm << "Done with problemSetup \t\t\t MPMICE2" << endl;
        cout_norm << "--------------------------------\n" << endl;
    }

    //__________________________________
    //  Set up data analysis modules
    d_analysisModules = AnalysisModuleFactory::create(d_myworld,
        m_materialManager,
        prob_spec);

    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->setComponents(dynamic_cast<ApplicationInterface*>(this));
            am->problemSetup(prob_spec, restart_prob_spec, grid,
                d_mpm->d_particleState, d_mpm->d_particleState_preReloc);
        }
    }
}

//______________________________________________________________________
//
void MPMICE2::outputProblemSpec(ProblemSpecP& root_ps)
{
    d_mpm->outputProblemSpec(root_ps);
    d_ice->outputProblemSpec(root_ps);

    // Global flags required by MPMICE2
    ProblemSpecP mpm_ps = root_ps->findBlock("MPM");
    mpm_ps->appendElement("testForNegTemps_mpm", d_testForNegTemps_mpm);

    //__________________________________
    //  output data analysis modules
    if (d_analysisModules.size() != 0) {

        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;

            am->outputProblemSpec(root_ps);
        }
    }
}

//______________________________________________________________________
//
void MPMICE2::scheduleInitialize(const LevelP& level,
    SchedulerP& sched)
{
    printSchedule(level, cout_doing, "MPMICE2::scheduleInitialize");

    d_mpm->scheduleInitialize(level, sched);
    d_ice->scheduleInitialize(level, sched);

    //__________________________________
    //  What isn't initialized in either ice or mpm
    Task* t = scinew Task("MPMICE2::actuallyInitialize",
        this, &MPMICE2::actuallyInitialize);

    // because there is ICE material inside the MPM material so there is no necessary to 
    // create CCLabel for MPM material in MPMICE2!!!

    // Get the material subsets
    const MaterialSubset* ice_matls = m_materialManager->allMaterials("ICE")->getUnion();
    //MaterialSubset* press_matl = d_ice->d_press_matl;

    // This is compute in d_ice->actuallyInitalize(...), and it is needed in 
    //  MPMICE2's actuallyInitialize()
    t->requires(Task::NewDW, Ilb->vol_frac_CCLabel, ice_matls, Ghost::None, 0);

    //t->requires(Task::NewDW, Ilb->Porosity_CCLabel, press_matl, Ghost::None, 0);

    if (d_switchCriteria) {
        d_switchCriteria->scheduleInitialize(level, sched);
    }

    //__________________________________
    // dataAnalysis 
    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->scheduleInitialize(sched, level);
        }
    }

    sched->addTask(t, level->eachPatch(), m_materialManager->allMaterials());

}

//______________________________________________________________________
//       A C T U A L   S T E P S :
//______________________________________________________________________
void MPMICE2::actuallyInitialize(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse*,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    new_dw->get(timeStep, VarLabel::find(timeStep_name));

    //bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing actuallyInitialize ");
        //__________________________________
        //output material indices
        if (patch->getID() == 0) {
            cout << "Materials Indicies:   MPM [" << *(m_materialManager->allMaterials("MPM")) << "] "
                << "ICE[" << *(m_materialManager->allMaterials("ICE")) << "]" << endl;

            cout << "Material Names:";
            unsigned int numAllMatls = m_materialManager->getNumMatls();
            for (unsigned int m = 0; m < numAllMatls; m++) {
                Material* matl = m_materialManager->getMaterial(m);
                cout << " " << matl->getDWIndex() << ") " << matl->getName();
            }
            cout << "\n";
        }

        // Sum variable for testing that the volume fractions sum to Porosity
        CCVariable<double> vol_frac_sum;
        new_dw->allocateTemporary(vol_frac_sum, patch);
        vol_frac_sum.initialize(0.0);

        CCVariable<double> errorThresholdTop, errorThresholdBottom;
        new_dw->allocateTemporary(errorThresholdTop, patch);
        new_dw->allocateTemporary(errorThresholdBottom, patch);
        errorThresholdTop.initialize(0.0);
        errorThresholdBottom.initialize(0.0);

        //constCCVariable<double> Porosity_CC;
        //new_dw->get(Porosity_CC, Ilb->Porosity_CCLabel, 0, patch, Ghost::None, 0);

        //___________________________________
        //   B U L L E T  P R O O F I N G
        // Verify volume fractions sum to 1.0
        // Loop through ICE materials to get their contribution to volume fraction
        unsigned int numICE_matls = m_materialManager->getNumMatls("ICE");
        for (unsigned int m = 0; m < numICE_matls; m++) {
            constCCVariable<double> vol_frac;
            ICEMaterial* ice_matl = (ICEMaterial*)m_materialManager->getMaterial("ICE", m);
            int indx = ice_matl->getDWIndex();

            // Get the Volume Fraction computed in ICE's actuallyInitialize(...)
            new_dw->get(vol_frac, Ilb->vol_frac_CCLabel, indx, patch, Ghost::None, 0);

            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                vol_frac_sum[c] += vol_frac[c];
                //errorThresholdTop[c] = Porosity_CC[c] + 1.0e-10;
                //errorThresholdBottom[c] = Porosity_CC[c] - 1.0e-10;
            }
        }  // num_ICE_matls loop


        // Loop of MPM Materials
        unsigned int numMPM_matls = m_materialManager->getNumMatls("MPM");
        for (unsigned int m = 0; m < numMPM_matls; m++) {

            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();

            CCVariable<double> heatFlux;
            new_dw->allocateAndPut(heatFlux, Mlb->heatRate_CCLabel, indx, patch);
            heatFlux.initialize(0.0);
        } // num_MPM_matls loop


        // Tempoerary turn it off
        /*
        for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            Point pt = patch->getCellPosition(c);

            if (!(vol_frac_sum[c] <= errorThresholdTop[c] && vol_frac_sum[c] >= errorThresholdBottom[c])) {

                ostringstream warn;
                warn << "ERROR MPMICE2::actuallyInitialize cell " << *iter << " position" << pt << "\n\n"
                    << "volume fraction (" << std::setprecision(13) << vol_frac_sum[*iter] << ") does not sum to 1 +- 1e-10.\n"
                    << "Verify that this region of the domain contains at least 1 geometry object.  If you're using the optional\n"
                    << "'volumeFraction' tags verify that they're correctly specified.\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
        } // cell iterator for volume fraction
        */
    } // Patch loop
}

//______________________________________________________________________
//
void MPMICE2::scheduleRestartInitialize(const LevelP& level,
    SchedulerP& sched)
{
    printSchedule(level, cout_doing, "MPMICE2::scheduleInitialize");

    d_mpm->scheduleRestartInitialize(level, sched);
    d_ice->scheduleRestartInitialize(level, sched);

    //__________________________________
    // dataAnalysis 
    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->scheduleRestartInitialize(sched, level);
        }
    }
}

//______________________________________________________________________
//
void MPMICE2::restartInitialize()
{
    if (cout_doing.active())
        cout_doing << "Doing restartInitialize \t\t\t MPMICE2" << endl;

    d_mpm->restartInitialize();
    d_ice->restartInitialize();

    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->restartInitialize();
        }
    }
}

//______________________________________________________________________
//
void MPMICE2::scheduleComputeStableTimeStep(const LevelP& level,
    SchedulerP& sched)
{
    // Schedule computing the ICE stable timestep
    d_ice->scheduleComputeStableTimeStep(level, sched);
    // MPM stable timestep is a by product of the CM
}

//______________________________________________________________________
//
void
MPMICE2::scheduleTimeAdvance(const LevelP& inlevel, SchedulerP& sched)
{
    // Only do scheduling on level 0 for lockstep AMR
    if (inlevel->getIndex() > 0 && isLockstepAMR())
        return;

    // If we have a finer level, then assume that we are doing multilevel MPMICE2
    // Otherwise, it is plain-ole MPMICE2
    do_mlmpmice2 = false;
    if (inlevel->hasFinerLevel()) {
        do_mlmpmice2 = true;
    }
    const LevelP& mpm_level = do_mlmpmice2 ? inlevel->getGrid()->getLevel(inlevel->getGrid()->numLevels() - 1) : inlevel;

    const PatchSet* mpm_patches = mpm_level->eachPatch();
    const MaterialSet* ice_matls = m_materialManager->allMaterials("ICE");
    const MaterialSet* mpm_matls = m_materialManager->allMaterials("MPM");
    const MaterialSet* all_matls = m_materialManager->allMaterials();
    MaterialSubset* press_matl = d_ice->d_press_matl;
    MaterialSubset* one_matl = d_ice->d_press_matl;

    const MaterialSubset* ice_matls_sub = ice_matls->getUnion();
    const MaterialSubset* mpm_matls_sub = mpm_matls->getUnion();
    cout_doing << "---------------------------------------------------------Level ";
    if (do_mlmpmice2) {
        cout_doing << inlevel->getIndex() << " (ICE) " << mpm_level->getIndex() << " (MPM)" << endl;;
    }
    else {
        cout_doing << inlevel->getIndex() << endl;
    }

    //__________________________________
    // Scheduling
    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
        const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
        d_ice->scheduleComputeThermoTransportProperties(sched, ice_level, ice_matls);

        d_ice->scheduleMaxMach_on_Lodi_BC_Faces(sched, ice_level, ice_matls);
    }

    d_mpm->scheduleApplyExternalLoads(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleComputeCurrentParticleSize(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleInterpolateParticlesToGrid(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleComputeInternalForce(sched, mpm_patches, mpm_matls);
    // Move internal force here because we need nodal stress but we need to set
    // p.pressure = 0 because the stress is subtracted p.pressure in MPMICE, not in
    // MPMICE2

    d_mpm->scheduleComputeHeatExchange(sched, mpm_patches, mpm_matls);
    if (d_mpm->flags->d_computeNormals) {
        d_mpm->scheduleComputeNormals(sched, mpm_patches, mpm_matls);
    }
    d_mpm->scheduleExMomInterpolated(sched, mpm_patches, mpm_matls);

    // schedule the interpolation of mass and volume to the cell centers
    scheduleInterpolateNCToCC_0(sched, mpm_patches, one_matl, press_matl,
        mpm_matls);

    // do coarsens in reverse order, and before the other tasks
    if (do_mlmpmice2) {
        for (int l = inlevel->getGrid()->numLevels() - 2; l >= 0; l--) {
            const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
            const PatchSet* ice_patches = ice_level->eachPatch();

            scheduleCoarsenCC_0(sched, ice_patches, mpm_matls);
            scheduleCoarsenNCMass(sched, ice_patches, mpm_matls);
        }
    }

    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
        const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
        const PatchSet* ice_patches = ice_level->eachPatch();

        // Only compute the ICE  pressure
        //d_ice->scheduleComputePressure(sched, ice_patches, press_matl,
        //    ice_matls);
        //scheduleComputePressure(sched, ice_patches, press_matl, ice_matls);

        scheduleComputePressure(sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            press_matl,
            all_matls);

        // Currently there is no need for face centered temp
        /*
        d_ice->scheduleComputeTempFC(sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            all_matls);
        */

        d_ice->scheduleComputeModelSources(sched, ice_level, all_matls);

        d_ice->scheduleUpdateVolumeFraction(sched, ice_level, press_matl,
            all_matls);

        // Replace the ComputeVel_FC by two functions
        /*
        d_ice->scheduleComputeVel_FC(sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            press_matl,
            all_matls);
            */

        scheduleComputeVelICE_FC(sched, ice_patches, ice_matls_sub,
            press_matl,
            ice_matls);

        scheduleComputeVelMPM_FC(sched, ice_patches, mpm_matls_sub,
            press_matl,
            mpm_matls);

        /*
        scheduleComputeVel_FC(sched, ice_patches, ice_matls_sub, mpm_matls_sub,
            press_matl,
            all_matls);
        */

        d_ice->d_exchModel->sched_PreExchangeTasks(sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            all_matls);

        d_ice->d_exchModel->sched_AddExch_VelFC(sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            all_matls,
            d_ice->d_BC_globalVars,
            false);
    }

    if (d_ice->d_impICE) {        //  I M P L I C I T, won't work with AMR yet
      // we should use the AMR multi-level pressure solve
        for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
            const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
            const PatchSet* ice_patches = ice_level->eachPatch();

            d_ice->scheduleSetupRHS(sched, ice_patches, one_matl,
                all_matls,
                false,
                "computes");
            d_ice->scheduleCompute_maxRHS(sched, ice_level, one_matl,
                all_matls);

            d_ice->scheduleImplicitPressureSolve(sched, ice_level, ice_patches,
                one_matl,
                press_matl,
                ice_matls_sub,
                mpm_matls_sub,
                all_matls);

            d_ice->scheduleComputeDel_P(sched, ice_level, ice_patches,
                one_matl,
                press_matl,
                all_matls);
        }
    }                           //  IMPLICIT AND EXPLICIT

    if (!(d_ice->d_impICE)) {       //  E X P L I C I T  (was not tested)
        for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
            const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
            const PatchSet* ice_patches = ice_level->eachPatch();

            d_ice->scheduleComputeDelPressAndUpdatePressCC(
                sched, ice_patches, press_matl,
                ice_matls_sub,
                mpm_matls_sub,
                ice_matls);
        }
    }  

    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
        const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
        const PatchSet* ice_patches = ice_level->eachPatch();

        d_ice->scheduleComputePressFC(sched, ice_patches, press_matl,
            all_matls);

        d_ice->scheduleVelTau_CC(sched, ice_patches, ice_matls);

        d_ice->scheduleViscousShearStress(sched, ice_patches, ice_matls);

        d_ice->scheduleAccumulateMomentumSourceSinks(
            sched, ice_patches, press_matl,
            ice_matls_sub,
            mpm_matls_sub,
            all_matls);

        d_ice->scheduleAccumulateEnergySourceSinks(sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            press_matl,
            all_matls);
    }

    // This one can be used for vizualization of pore pressure in particles if (vizualization)
    //if (!d_rigidMPM) {
        //    if(do_mlmpmice2){
        //      scheduleRefinePressCC(                  sched, mpm_patches, press_matl,
        //                                                                  mpm_matls);
        //    }
        
        scheduleInterpolatePressCCToPressNC(sched, mpm_patches, press_matl,
            mpm_matls);
        scheduleInterpolatePAndGradP(sched, mpm_patches, press_matl,
            one_matl,
            mpm_matls_sub,
            mpm_matls);
            
    //}

    d_mpm->scheduleComputeInternalHeatRate(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleComputeNodalHeatFlux(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleSolveHeatEquations(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleComputeAndIntegrateAcceleration(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleIntegrateTemperatureRate(sched, mpm_patches, mpm_matls);

    scheduleComputeLagrangianValuesMPM(sched, mpm_patches, one_matl,
        mpm_matls);

    // do coarsens in reverse order, and before the other tasks
    if (do_mlmpmice2) {
        for (int l = inlevel->getGrid()->numLevels() - 2; l >= 0; l--) {
            const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
            const PatchSet* ice_patches = ice_level->eachPatch();
            scheduleCoarsenLagrangianValuesMPM(sched, ice_patches, mpm_matls);
        }
    }

    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
        const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
        const PatchSet* ice_patches = ice_level->eachPatch();

        d_ice->scheduleComputeLagrangianValues(sched, ice_patches, press_matl, ice_matls);

        d_ice->d_exchModel->sched_AddExch_Vel_Temp_CC(
            sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            all_matls,
            d_ice->d_BC_globalVars);

        scheduleComputeLagrangianSpecificVolume(
            sched, ice_patches, ice_matls_sub,
            mpm_matls_sub,
            press_matl,
            all_matls);

        d_ice->scheduleComputeLagrangian_Transported_Vars(
            sched, ice_patches, ice_matls);

    }

    scheduleComputeCCVelAndTempRates(sched, mpm_patches, mpm_matls);

    //  if(do_mlmpmice2){
    //    scheduleRefineCC(                         sched, mpm_patches, mpm_matls);
    //  }

    scheduleInterpolateCCToNC(sched, mpm_patches, mpm_matls);

    d_mpm->scheduleExMomIntegrated(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleSetGridBoundaryConditions(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleInterpolateToParticlesAndUpdate(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleComputeParticleGradients(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleComputeStressTensor(sched, mpm_patches, mpm_matls);
    d_mpm->scheduleFinalParticleUpdate(sched, mpm_patches, mpm_matls);

    if (d_mpm->flags->d_computeScaleFactor) {
        d_mpm->scheduleComputeParticleScaleFactor(sched, mpm_patches, mpm_matls);
    }

    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
        const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
        const PatchSet* ice_patches = ice_level->eachPatch();

        d_ice->scheduleAdvectAndAdvanceInTime(sched, ice_patches, ice_matls_sub,
            ice_matls);

        d_ice->scheduleConservedtoPrimitive_Vars(sched, ice_patches, ice_matls_sub, press_matl,
            ice_matls, "afterAdvection");
    }
} // end scheduleTimeAdvance()

// Start the schedule
//______________________________________________________________________
//
void MPMICE2::scheduleInterpolateNCToCC_0(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* one_matl,
    const MaterialSubset* press_matl,
    const MaterialSet* mpm_matls)
{
    if (d_mpm->flags->doMPMOnLevel(getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {

        printSchedule(patches, cout_doing, "MPMICE2::scheduleInterpolateNCToCC_0");

        /* interpolateNCToCC */
        Task* t = scinew Task("MPMICE2::interpolateNCToCC_0",
            this, &MPMICE2::interpolateNCToCC_0);
        const MaterialSubset* mss = mpm_matls->getUnion();
        t->requires(Task::NewDW, Mlb->gMassLabel, Ghost::AroundCells, 1);
        t->requires(Task::NewDW, Mlb->gStressForSavingLabel, Ghost::AroundCells, 1);
        t->requires(Task::NewDW, Mlb->gVolumeLabel, Ghost::AroundCells, 1);
        t->requires(Task::NewDW, Mlb->gVolumeFractionLabel, Ghost::AroundCells, 1);
        t->requires(Task::NewDW, Mlb->gVelocityBCLabel, Ghost::AroundCells, 1);
        t->requires(Task::NewDW, Mlb->gTemperatureLabel, Ghost::AroundCells, 1);
        t->requires(Task::OldDW, Mlb->NC_CCweightLabel, one_matl,
            Ghost::AroundCells, 1);
        //t->requires(Task::OldDW, Ilb->sp_vol_CCLabel, Ghost::None, 0);
        //t->requires(Task::OldDW, MIlb->temp_CCLabel, Ghost::None, 0);

        t->requires(Task::OldDW, Ilb->timeStepLabel);

        t->computes(MIlb->cMassLabel);
        //t->computes(MIlb->cVolumeLabel);
        //t->computes(Ilb->Porosity_CCLabel);
        t->computes(Ilb->Porosity_CCLabel, press_matl);
        t->computes(MIlb->cPermeabilityLabel);
        t->computes(MIlb->vel_CCLabel);
        t->computes(MIlb->temp_CCLabel);
        t->computes(Ilb->rho_CCLabel, mss);
        t->computes(MIlb->stress_CCLabel);
        t->computes(Ilb->VolumeFraction_CCLabel, mss);

        t->computes(MIlb->stressX_CCLabel);
        t->computes(MIlb->stressY_CCLabel);
        t->computes(MIlb->stressZ_CCLabel);

        sched->addTask(t, patches, mpm_matls);
    }
}

//______________________________________________________________________
//
void MPMICE2::interpolateNCToCC_0(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, VarLabel::find(timeStep_name));

    bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing interpolateNCToCC_0");

        unsigned int numMatls = m_materialManager->getNumMatls("MPM");
        Vector dx = patch->dCell();
        double cell_vol = dx.x() * dx.y() * dx.z();
        Ghost::GhostType  gac = Ghost::AroundCells;

        constNCVariable<double> NC_CCweight;
        old_dw->get(NC_CCweight, Mlb->NC_CCweightLabel, 0, patch, gac, 1);

        CCVariable<double> Porosity_CC;
        new_dw->allocateAndPut(Porosity_CC, Ilb->Porosity_CCLabel, 0, patch);
        Porosity_CC.initialize(1.0);

        for (unsigned int m = 0; m < numMatls; m++) {
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();
            // Create arrays for the grid data
            //constNCVariable<double> gSp_vol, gPorosity;
            constNCVariable<double> gmass, gvolume, gtemperature, gVolumeFraction;
            constNCVariable<Vector> gvelocity;
            CCVariable<double> cmass, Temp_CC, Permeability_CC, VolumeFraction_CC;
            CCVariable<double> rho_CC;
            CCVariable<Vector> vel_CC;
            constNCVariable<Matrix3> gstress;
            CCVariable<Matrix3> stress_CC;

            new_dw->allocateAndPut(cmass, MIlb->cMassLabel, indx, patch);
            //new_dw->allocateAndPut(Volume_CC, MIlb->cVolumeLabel, indx, patch);
            //new_dw->allocateAndPut(Porosity_CC, Ilb->Porosity_CCLabel, indx, patch);    
            new_dw->allocateAndPut(VolumeFraction_CC, Ilb->VolumeFraction_CCLabel, indx, patch);
            VolumeFraction_CC.initialize(0.0);
            new_dw->allocateAndPut(Permeability_CC, MIlb->cPermeabilityLabel, indx, patch);
            new_dw->allocateAndPut(vel_CC, MIlb->vel_CCLabel, indx, patch);
            new_dw->allocateAndPut(Temp_CC, MIlb->temp_CCLabel, indx, patch);
            new_dw->allocateAndPut(rho_CC, Ilb->rho_CCLabel, indx, patch);
            new_dw->allocateAndPut(stress_CC, MIlb->stress_CCLabel, indx, patch);

            // Hacking
            CCVariable<double> stressX_CC, stressY_CC, stressZ_CC;
            new_dw->allocateAndPut(stressX_CC, MIlb->stressX_CCLabel, indx, patch);
            new_dw->allocateAndPut(stressY_CC, MIlb->stressY_CCLabel, indx, patch);
            new_dw->allocateAndPut(stressZ_CC, MIlb->stressZ_CCLabel, indx, patch);
            stressX_CC.initialize(0.0);
            stressY_CC.initialize(0.0);
            stressZ_CC.initialize(0.0);

            double very_small_mass = d_TINY_RHO * cell_vol;
            cmass.initialize(very_small_mass);
            stress_CC.initialize(Matrix3(0.0));

            new_dw->get(gmass, Mlb->gMassLabel, indx, patch, gac, 1);
            new_dw->get(gstress, Mlb->gStressForSavingLabel, indx, patch, gac, 1);
            new_dw->get(gvolume, Mlb->gVolumeLabel, indx, patch, gac, 1);
            new_dw->get(gVolumeFraction, Mlb->gVolumeFractionLabel, indx, patch, gac, 1);
            new_dw->get(gvelocity, Mlb->gVelocityBCLabel, indx, patch, gac, 1);
            new_dw->get(gtemperature, Mlb->gTemperatureLabel, indx, patch, gac, 1);
            //new_dw->get(gSp_vol, Mlb->gSp_volLabel, indx, patch, gac, 1);
            //old_dw->get(sp_vol_CC_ice, Ilb->sp_vol_CCLabel, indx, patch, gn, 0);
            //old_dw->get(Temp_CC_ice, MIlb->temp_CCLabel, indx, patch, gn, 0);
            IntVector nodeIdx[8];
            double IniPermeability = mpm_matl->getInitialPermeability();

                  //__________________________________
                  //  compute CC Variables
            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                patch->findNodesFromCell(*iter, nodeIdx);

                // double Vol_CC_mpm = 0.0;
                //double Porosity_CC_mpm = 0.0;
                double Volumefraction_CC_mpm = 0.0;
                double Temp_CC_mpm = 0.0;
                Vector vel_CC_mpm = Vector(0.0, 0.0, 0.0);
                Matrix3 stress_CC_mpm = Matrix3(0.0);

                for (int in = 0; in < 8; in++) {
                    double NC_CCw_mass = NC_CCweight[nodeIdx[in]] * gmass[nodeIdx[in]];
                    cmass[c] += NC_CCw_mass;
                  
                    Volumefraction_CC_mpm += gVolumeFraction[nodeIdx[in]] * NC_CCw_mass;
                    vel_CC_mpm += gvelocity[nodeIdx[in]] * NC_CCw_mass;
                    Temp_CC_mpm += gtemperature[nodeIdx[in]] * NC_CCw_mass;
                    stress_CC_mpm += gstress[nodeIdx[in]] * NC_CCw_mass;

                    // Average in a cell rather than interpolation to avoid diffusion
                    //Volumefraction_CC_mpm += gVolumeFraction[nodeIdx[in]] * 0.125;
                    //vel_CC_mpm += gvelocity[nodeIdx[in]] * 0.125;
                    //Temp_CC_mpm += gtemperature[nodeIdx[in]] * 0.125;
                    //stress_CC_mpm += gstress[nodeIdx[in]] * 0.125;
                }
                           
                double inv_cmass = 1.0 / cmass[c];
                vel_CC_mpm *= inv_cmass;
                Temp_CC_mpm *= inv_cmass;
                Volumefraction_CC_mpm *= inv_cmass;
                stress_CC_mpm *= inv_cmass;
                
                // Locate on cell
                vel_CC[c] = vel_CC_mpm;
                stress_CC[c] = stress_CC_mpm;
                Temp_CC[c] = Temp_CC_mpm;
                VolumeFraction_CC[c] = Volumefraction_CC_mpm;

                //cerr << " Porosity_CC" << Porosity_CC[c] << endl;

                // Global porosity variable (press_matl) 1 - sum(volumrfractionMPM)
                if (cmass[c] > very_small_mass) { // For cell with particles (only MPM material)
                    Porosity_CC[c] -= Volumefraction_CC_mpm;
                }            

                //if (Porosity_CC[c] > 1) { // Cap the porosity
                //    Porosity_CC[c] = 1;
                //}

                //cerr << " VolumeFraction_CC[c]" << VolumeFraction_CC[c] << endl;

                rho_CC[c] = cmass[c] / cell_vol;

                // Hacking here because stuggling to deal with stress_CC when computing face centered velocity
                stressX_CC[c] = stress_CC[c](0, 0);
                stressY_CC[c] = stress_CC[c](1, 1);
                stressZ_CC[c] = stress_CC[c](2, 2);

                // Calculate cell centered permeability (temporary set constant but should be porosity dependent)
                Permeability_CC[c] = IniPermeability;
            } // cell
                          
            //  Set BC's
            setBC(Temp_CC, "Temperature", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(rho_CC, "Density", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(vel_CC, "Velocity", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            //  Set if symmetric Boundary conditions
            setBC(cmass, "set_if_sym_BC", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(Porosity_CC, "set_if_sym_BC", patch, m_materialManager, 0, new_dw, isNotInitialTimeStep);

            //---- B U L L E T   P R O O F I N G------
            // ignore BP if recompute time step has already been requested
            IntVector neg_cell;
            ostringstream warn;
            bool rts = new_dw->recomputeTimeStep();

            int L = getLevel(patches)->getIndex();
            if (d_testForNegTemps_mpm) {
                if (!areAllValuesPositive(Temp_CC, neg_cell) && !rts) {
                    warn << "ERROR MPMICE2:(" << L << "):interpolateNCToCC_0, mat " << indx
                        << " cell "
                        << neg_cell << " Temp_CC " << Temp_CC[neg_cell] << "\n ";
                    throw InvalidValue(warn.str(), __FILE__, __LINE__);
                }
            }
            if (!areAllValuesPositive(rho_CC, neg_cell) && !rts) {
                warn << "ERROR MPMICE2:(" << L << "):interpolateNCToCC_0, mat " << indx
                    << " cell " << neg_cell << " rho_CC " << rho_CC[neg_cell] << "\n ";
                throw InvalidValue(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesPositive(Porosity_CC, neg_cell) && !rts) {
                warn << "ERROR MPMICE2:(" << L << "):interpolateNCToCC_0, mat " << indx
                    << " cell " << neg_cell << " Porosity_CC " << Porosity_CC[neg_cell] << "\n ";
                throw InvalidValue(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesPositive(Permeability_CC, neg_cell) && !rts) {
                warn << "ERROR MPMICE2:(" << L << "):interpolateNCToCC_0, mat " << indx
                    << " cell " << neg_cell << " Permeability_CC " << Permeability_CC[neg_cell] << "\n ";
                throw InvalidValue(warn.str(), __FILE__, __LINE__);
            }

            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;

                //cerr << "interpolateNCToCC_0 " << endl;
                //cerr << "cmass " << cmass[c] << endl;
                //cerr << "vel_CC " << vel_CC[c] << endl;
                //cerr << "stress_CC " << stress_CC[c] << endl;
                //cerr << "Temp_CC " << Temp_CC[c] << endl;
                //cerr << "Porosity_CC " << Porosity_CC[c] << endl;
                //cerr << "rho_CC " << rho_CC[c] << endl;
            }

        }  //materials
    }  //patches
}

/*_____________________________________________________________________
 Function~  MPMICE2::scheduleComputePressure--
 Note:  Handles both Rate and Equilibration form of solution technique
 This task is a bit strange as the index of materials is important for 
 the boundary condition so we also need MPM material index even though
 MPMICE2 does not update any MPM's variables here
_____________________________________________________________________*/
void MPMICE2::scheduleComputePressure(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* ice_matls,
    const MaterialSubset* mpm_matls,
    const MaterialSubset* press_matl,
    const MaterialSet* all_matls)
{
    Task* t = nullptr;

    printSchedule(patches, cout_doing, "MPMICE2::scheduleComputeEquilibrationPressure");

    if (m_materialManager->getNumMatls("ICE") == 1) {
        t = scinew Task("MPMICE2::computeEquilPressure_1_matl",
            this, &MPMICE2::computeEquilPressure_1_matl);
    }
    else {
        t = scinew Task("MPMICE2::computeEquilibrationPressure",
            this, &MPMICE2::computeEquilibrationPressure);
    }
    
    Ghost::GhostType  gn = Ghost::None;

    t->requires(Task::OldDW, Ilb->timeStepLabel);
    t->requires(Task::OldDW, Ilb->delTLabel, getLevel(patches));
    t->requires(Task::NewDW, Ilb->Porosity_CCLabel, press_matl, gn);
    t->requires(Task::OldDW, Ilb->press_CCLabel, press_matl, gn);

    // I C E   
    t->requires(Task::OldDW, Ilb->temp_CCLabel, ice_matls, gn);
    t->requires(Task::OldDW, Ilb->rho_CCLabel, ice_matls, gn);
    t->requires(Task::OldDW, Ilb->sp_vol_CCLabel, ice_matls, gn);
    t->requires(Task::NewDW, Ilb->specific_heatLabel, ice_matls, gn);
    t->requires(Task::NewDW, Ilb->gammaLabel, ice_matls, gn);   

    // M P M
    t->requires(Task::NewDW, MIlb->temp_CCLabel, mpm_matls, gn);
    t->requires(Task::NewDW, Ilb->rho_CCLabel, mpm_matls, gn);

    computesRequires_CustomBCs(t, "EqPress", Ilb, ice_matls,
        d_ice->d_BC_globalVars);

    //  A L L _ M A T L S
    t->computes(Ilb->f_theta_CCLabel);
    t->computes(Ilb->compressibilityLabel);
    t->computes(Ilb->speedSound_CCLabel);  // Value of MPM will not be used !!
    t->computes(Ilb->vol_frac_CCLabel);    
    t->computes(Ilb->sumKappaLabel, press_matl);
    //t->computes(Ilb->TMV_CCLabel, press_matl);
    t->computes(Ilb->press_equil_CCLabel, press_matl);
    t->computes(Ilb->sum_imp_delPLabel, press_matl);  // needed by implicit ICE
    t->computes(Ilb->sp_vol_CCLabel);
    t->computes(Ilb->rho_micro_CCLabel);
    t->computes(Ilb->rho_CCLabel, ice_matls);
    //t->computes(Ilb->rho1_CCLabel);

    sched->addTask(t, patches, all_matls);
}

/* _____________________________________________________________________
 Function~  MPMICE2::computeEquilPressure_1_matl--
 Purpose~   Simple EOS evaluation
_____________________________________________________________________*/
void MPMICE2::computeEquilPressure_1_matl(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset* matls,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, Ilb->timeStepLabel);

    bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);

        printTask(patches, patch, cout_doing, "Doing MPMICE2::computeEquilPressure_1_matl");

        unsigned int numICEMatls = m_materialManager->getNumMatls("ICE");
        unsigned int numMPMMatls = m_materialManager->getNumMatls("MPM");
        unsigned int numALLMatls = numICEMatls + numMPMMatls;

        //rho_micro also need MPM mateiral index
        //std::vector<double> dp_drho(numALLMatls), dp_de(numALLMatls);
        std::vector<CCVariable<double> > vol_frac(numALLMatls);
        std::vector<CCVariable<double> > rho_micro(numALLMatls);
        std::vector<CCVariable<double> > rho_CC_new(numALLMatls);
        //std::vector<CCVariable<double> > rho1_CC(numALLMatls);
        std::vector<CCVariable<double> > sp_vol_new(numALLMatls);
        std::vector<CCVariable<double> > speedSound(numALLMatls);
        std::vector<CCVariable<double> > speedSound_new(numALLMatls);
        std::vector<CCVariable<double> > f_theta(numALLMatls);
        std::vector<CCVariable<double> > kappa(numALLMatls);
        std::vector<constCCVariable<double> > Temp(numALLMatls);
        std::vector<constCCVariable<double> > rho_CC(numALLMatls);
        std::vector<constCCVariable<double> > sp_vol_CC(numALLMatls);
        std::vector<constCCVariable<double> > cv(numALLMatls);
        std::vector<constCCVariable<double> > gamma(numALLMatls);
        std::vector<constCCVariable<double> > placeHolder(0);

        std::vector<MPMMaterial*> mpm_matl(numALLMatls);
        std::vector<ICEMaterial*> ice_matl(numALLMatls);
        for (unsigned int m = 0; m < numALLMatls; m++) {
            Material* matl = m_materialManager->getMaterial(m);
            ice_matl[m] = dynamic_cast<ICEMaterial*>(matl);
            mpm_matl[m] = dynamic_cast<MPMMaterial*>(matl);
        }

        //__________________________________ 
        Ghost::GhostType  gn = Ghost::None;
        constCCVariable<double> Porosity_CC;
        new_dw->get(Porosity_CC, Ilb->Porosity_CCLabel, 0, patch, gn, 0);
        CCVariable<double> press_eq, sumKappa, sum_imp_delP;
        new_dw->allocateAndPut(press_eq, Ilb->press_equil_CCLabel, 0, patch);
        new_dw->allocateAndPut(sumKappa, Ilb->sumKappaLabel, 0, patch);
        new_dw->allocateAndPut(sum_imp_delP, Ilb->sum_imp_delPLabel, 0, patch);
        sum_imp_delP.initialize(0.0);

        for (unsigned int m = 0; m < numALLMatls; m++) {
            Ghost::GhostType  gn = Ghost::None;
            Material* matl = m_materialManager->getMaterial(m);
            int indx = matl->getDWIndex();
            if (ice_matl[m]) {                    // I C E

                old_dw->get(Temp[m], Ilb->temp_CCLabel, indx, patch, gn, 0);
                old_dw->get(rho_CC[m], Ilb->rho_CCLabel, indx, patch, gn, 0);
                old_dw->get(sp_vol_CC[m], Ilb->sp_vol_CCLabel, indx, patch, gn, 0);
                new_dw->get(cv[m], Ilb->specific_heatLabel, indx, patch, gn, 0);
                new_dw->get(gamma[m], Ilb->gammaLabel, indx, patch, gn, 0);
                new_dw->allocateAndPut(rho_CC_new[m], Ilb->rho_CCLabel, indx, patch);
                        
            }
            if (mpm_matl[m]) {                    // M P M
                new_dw->get(Temp[m], MIlb->temp_CCLabel, indx, patch, gn, 0);
                new_dw->get(rho_CC[m], Ilb->rho_CCLabel, indx, patch, gn, 0);
            }
            new_dw->allocateAndPut(kappa[m], Ilb->compressibilityLabel, indx, patch);
            new_dw->allocateAndPut(sp_vol_new[m], Ilb->sp_vol_CCLabel, indx, patch);
            new_dw->allocateAndPut(speedSound_new[m], Ilb->speedSound_CCLabel,
                indx, patch);
            new_dw->allocateAndPut(vol_frac[m], Ilb->vol_frac_CCLabel, indx, patch);
            //new_dw->allocateAndPut(rho1_CC[m], Ilb->rho1_CCLabel, indx, patch);
            new_dw->allocateAndPut(rho_micro[m], Ilb->rho_micro_CCLabel, indx, patch);
            //new_dw->allocateTemporary(rho_micro[m], patch);
            new_dw->allocateAndPut(f_theta[m], Ilb->f_theta_CCLabel, indx, patch);
        }

        //______________________________________________________________________
        //  Main loop
        for (unsigned int m = 0; m < numALLMatls; m++) {
            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;

                if (ice_matl[m]) {                // I C E
                    //vol_frac[m][c] = 1;
                    vol_frac[m][c] = Porosity_CC[c];
                    rho_micro[m][c] = rho_CC_new[m][c] = rho_CC[m][c];
                    //rho1_CC[m][c] = rho_CC_new[m][c] * Porosity_CC[c];   
                    sp_vol_new[m][c] = 1.0 / rho_CC[m][c];
                    double dp_drho, dp_de, c_2;

                    //__________________________________
                    // evaluate EOS
                    ice_matl[m]->getEOS()->computePressEOS(rho_micro[m][c], gamma[m][c],
                        cv[m][c], Temp[m][c], press_eq[c],
                        dp_drho, dp_de);

                    c_2 = dp_drho + dp_de * press_eq[c] / (rho_micro[m][c] * rho_micro[m][c]);
                    speedSound_new[m][c] = sqrt(c_2);
                    
                    //  compute f_theta  
                    kappa[m][c] = sp_vol_new[m][c] / (speedSound_new[m][c] * speedSound_new[m][c]);
                    //sumKappa[c] = kappa[m][c];
                    //f_theta[m][c] = 1.0;
                }

                else if (mpm_matl[m]) {                //  M P M  I need rho_micro index for setBC
                    rho_micro[m][c] = mpm_matl[m]->getInitialDensity();
                    //rho1_CC[m][c] = rho_micro[m][c] * (1 - Porosity_CC[c]);
                    //speedSound_new[m][c] = 30;
                    //vol_frac[m][c] = 0;
                    vol_frac[m][c] = (1 - Porosity_CC[c]);
                    sp_vol_new[m][c] = 1.0 / rho_CC[m][c];
                    //kappa[m][c] = 0.00001;

                    double dp_drho, c_2;
                    //double press_eq_MPM;

                    // Pointwise computation of thermodynamic quantities from Murnaghan
                    // Parameter
                    double P0 = 101325;
                    double n = 7.4;
                    double K = 39e-11;
                    double rho0 = 2000;

                    if (rho_micro[m][c] >= rho0) {
                        //press_eq_MPM = P0 + (1. / (n * K)) * (pow(rho_micro[m][c] / rho0, n) - 1.);
                        dp_drho = (1. / (K * rho0)) * pow((rho_micro[m][c] / rho0), n - 1.);
                    }
                    else {
                        //press_eq_MPM = P0 * pow(rho_micro[m][c] / rho0, (1. / (K * P0)));
                        dp_drho = (1. / (K * rho0)) * pow(rho_micro[m][c] / rho0, (1. / (K * P0) - 1.));
                    }

                    c_2 = dp_drho;
                    speedSound[m][c] = sqrt(c_2);

                    //  compute f_theta  
                    kappa[m][c] = sp_vol_new[m][c] / (speedSound[m][c] * speedSound[m][c]);
                    //sumKappa[c] = kappa[c];
                    //f_theta[m][c] = 1.0;                
                }

                sumKappa[c] += vol_frac[m][c] * kappa[m][c];
                f_theta[m][c] = vol_frac[m][c] * kappa[m][c] / sumKappa[c];
                
            }
        }

        //__________________________________
        // - apply Boundary conditions
        customBC_localVars* BC_localVars = scinew customBC_localVars();

        preprocess_CustomBCs("EqPress", old_dw, new_dw, Ilb, patch,
            999, d_ice->d_BC_globalVars, BC_localVars);

        setBC(press_eq, rho_micro, placeHolder, d_ice->d_surroundingMatl_indx,
            "rho_micro", "Pressure", patch, m_materialManager, 0, new_dw,
            d_ice->d_BC_globalVars, BC_localVars, isNotInitialTimeStep/*, d_with_mpmice2*/);

        delete_CustomBCs(d_ice->d_BC_globalVars, BC_localVars);

    }  // patch loop
}

/* _____________________________________________________________________
 Function~  MPMICE2::computeEquilibrationPressure--
 Purpose~   Find the equilibration pressure
 Reference: Flow of Interpenetrating Material Phases, J. Comp, Phys
               18, 440-464, 1975, see the equilibration section

 Steps
 ----------------
    - Compute rho_micro_CC, SpeedSound, vol_frac
    For each cell
    _ WHILE LOOP(convergence, max_iterations)
        - compute the pressure and dp_drho from the EOS of each material.
        - Compute delta Pressure
        - Compute delta volume fraction and update the
          volume fraction and the celldensity.
        - Test for convergence of delta pressure and delta volume fraction
    - END WHILE LOOP
    - bulletproofing
    end
Note:  The nomenclature follows the reference.
_____________________________________________________________________*/
void MPMICE2::computeEquilibrationPressure(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)

{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, Ilb->timeStepLabel);

    bool isNotInitialTimeStep = (timeStep > 0);

    const Level* level = getLevel(patches);
    int L_indx = level->getIndex();

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing MPMICE2::computeEquilibrationPressure");

        double    converg_coeff = 10e7;
        double    convergence_crit = converg_coeff * DBL_EPSILON;
        double    sum = 0., tmp;

        unsigned int numICEMatls = m_materialManager->getNumMatls("ICE");
        unsigned int numMPMMatls = m_materialManager->getNumMatls("MPM");
        unsigned int numALLMatls = numICEMatls + numMPMMatls;

        //unsigned int       numMatls = m_materialManager->getNumMatls("ICE");
        static int n_passes;
        n_passes++;
        
        //rho_micro also need MPM mateiral index
        std::vector<double> press_eos(numALLMatls);
        std::vector<double> dp_drho(numALLMatls), dp_de(numALLMatls);
        std::vector<CCVariable<double> > vol_frac(numALLMatls);
        std::vector<CCVariable<double> > rho_micro(numALLMatls);
        std::vector<CCVariable<double> > rho_CC_new(numALLMatls);
        //std::vector<CCVariable<double> > rho1_CC(numALLMatls);
        std::vector<CCVariable<double> > sp_vol_new(numALLMatls);
        std::vector<CCVariable<double> > speedSound(numALLMatls);
        std::vector<CCVariable<double> > speedSound_new(numALLMatls);
        std::vector<CCVariable<double> > f_theta(numALLMatls);
        std::vector<CCVariable<double> > kappa(numALLMatls);
        std::vector<constCCVariable<double> > Temp(numALLMatls);
        std::vector<constCCVariable<double> > rho_CC(numALLMatls);
        std::vector<constCCVariable<double> > sp_vol_CC(numALLMatls);
        std::vector<constCCVariable<double> > cv(numALLMatls);
        std::vector<constCCVariable<double> > gamma(numALLMatls);
        std::vector<constCCVariable<double> > placeHolder(0);

        CCVariable<int> n_iters_equil_press;
        constCCVariable<double> press;
        constCCVariable<double> Porosity_CC;
        CCVariable<double> press_new, sumKappa, sum_imp_delP;
        Ghost::GhostType  gn = Ghost::None;

        //__________________________________ 
        old_dw->get(press, Ilb->press_CCLabel, 0, patch, gn, 0);
        new_dw->get(Porosity_CC, Ilb->Porosity_CCLabel, 0, patch, gn, 0);
        new_dw->allocateAndPut(press_new, Ilb->press_equil_CCLabel, 0, patch);
        new_dw->allocateAndPut(sumKappa, Ilb->sumKappaLabel, 0, patch);
        new_dw->allocateAndPut(sum_imp_delP, Ilb->sum_imp_delPLabel, 0, patch);

        sum_imp_delP.initialize(0.0); //-- initialize for implicit pressure

        std::vector<MPMMaterial*> mpm_matl(numALLMatls);
        std::vector<ICEMaterial*> ice_matl(numALLMatls);
        for (unsigned int m = 0; m < numALLMatls; m++) {
            Material* matl = m_materialManager->getMaterial(m);
            ice_matl[m] = dynamic_cast<ICEMaterial*>(matl);
            mpm_matl[m] = dynamic_cast<MPMMaterial*>(matl);
        }   

        for (unsigned int m = 0; m < numALLMatls; m++) {
            Material* matl = m_materialManager->getMaterial(m);
            int indx = matl->getDWIndex();
            if (ice_matl[m]) {                    // I C E

                old_dw->get(Temp[m], Ilb->temp_CCLabel, indx, patch, gn, 0);
                old_dw->get(rho_CC[m], Ilb->rho_CCLabel, indx, patch, gn, 0);
                old_dw->get(sp_vol_CC[m], Ilb->sp_vol_CCLabel, indx, patch, gn, 0);
                new_dw->get(cv[m], Ilb->specific_heatLabel, indx, patch, gn, 0);
                new_dw->get(gamma[m], Ilb->gammaLabel, indx, patch, gn, 0);               
                new_dw->allocateAndPut(rho_CC_new[m], Ilb->rho_CCLabel, indx, patch);              
                
                
            }
            if (mpm_matl[m]) {                    // M P M
                new_dw->get(Temp[m], MIlb->temp_CCLabel, indx, patch, gn, 0);
                new_dw->get(rho_CC[m], Ilb->rho_CCLabel, indx, patch, gn, 0);
                //new_dw->get(sp_vol_CC[m], Ilb->sp_vol_CCLabel, indx, patch, gn, 0);
            }
            new_dw->allocateAndPut(sp_vol_new[m], Ilb->sp_vol_CCLabel, indx, patch);
            new_dw->allocateAndPut(speedSound_new[m], Ilb->speedSound_CCLabel,
                indx, patch);
            new_dw->allocateAndPut(vol_frac[m], Ilb->vol_frac_CCLabel, indx, patch);
            //new_dw->allocateAndPut(rho1_CC[m], Ilb->rho1_CCLabel, indx, patch);
            new_dw->allocateAndPut(rho_micro[m], Ilb->rho_micro_CCLabel, indx, patch);
            new_dw->allocateAndPut(kappa[m], Ilb->compressibilityLabel, indx, patch);
            new_dw->allocateAndPut(f_theta[m], Ilb->f_theta_CCLabel, indx, patch);
            //new_dw->allocateTemporary(rho_micro[m], patch);           
        }       
        press_new.copyData(press);
    
        //__________________________________
        // Compute rho_micro, volfrac
        for (unsigned int m = 0; m < numALLMatls; m++) {
            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;

                if (ice_matl[m]) {                // I C E
                    rho_micro[m][c] = 1.0 / sp_vol_CC[m][c];
                    //vol_frac[m][c] = Porosity_CC[c] * rho_CC[m][c] * sp_vol_CC[m][c];
                    vol_frac[m][c] = rho_CC[m][c] * sp_vol_CC[m][c];
                }
            
                else if (mpm_matl[m]) {                //  M P M  I need rho_micro index for setBC
                    rho_micro[m][c] = mpm_matl[m]->getInitialDensity();
                    vol_frac[m][c] = (1 - Porosity_CC[c]);
                    //vol_frac[m][c] = 0;
                }
            }
        }

        // Debug 
        for (unsigned int m = 0; m < numALLMatls; m++) {
            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;

                //cerr << "Cell vol frac of m " << m << " of cell " << c << " " << vol_frac[m][c] << endl;

                if (ice_matl[m]) {                // I C E
                   //cerr << "Cell vol frac of m " << m << " of cell " << c << " " << vol_frac[m][c] << endl;
                   //cerr << "Cell rho_CC of m " << m << " of cell " << c << " " << rho_CC[m][c] << endl;
                   //cerr << "Cell sp_vol_CC of m " << m << " of cell " << c << " " << sp_vol_CC[m][c] << endl;
                }
            }
        }
      
        //______________________________________________________________________
        // Done with preliminary calcs, now loop over every cell
        int count, test_max_iter = 0;
        for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            double delPress = 0.;
            bool converged = false;
            count = 0;
            vector<EqPress_dbg> dbgEqPress;

            while (count < d_ice->d_max_iter_equilibration && converged == false) {
                count++;
               
                //cerr << "interation" << count << endl;

                //__________________________________
                // evaluate press_eos at cell i,j,k               
                for (unsigned int m = 0; m < numALLMatls; m++) {
                    if (ice_matl[m]) {    // ICE
                        ice_matl[m]->getEOS()->computePressEOS(rho_micro[m][c], gamma[m][c],
                            cv[m][c], Temp[m][c], press_eos[m], dp_drho[m], dp_de[m]);
                    }
                }
                   
                //__________________________________
                // - compute delPress
                // - update press_CC     
                double A = 0., B = 0., C = 0.;
                for (unsigned int m = 0; m < numALLMatls; m++) {
                    if (ice_matl[m]) {
                        double Q = press_new[c] - press_eos[m];
                        double div_y = (vol_frac[m][c] * vol_frac[m][c])
                            / (dp_drho[m] * rho_CC[m][c] + d_SMALL_NUM);
                        A += vol_frac[m][c];
                        B += Q * div_y;
                        C += div_y;
                    }
                }
                //double vol_frac_not_close_packed = 1.0;
                double vol_frac_not_close_packed = Porosity_CC[c];
                delPress = (A - vol_frac_not_close_packed - B) / C;

                press_new[c] += delPress;

                //cerr << "delPress" << delPress << endl;

                //cerr << "press_new[c]" << press_new[c] << endl;

                //__________________________________
                // backout rho_micro_CC at this new pressure
                for (unsigned int m = 0; m < numALLMatls; m++) {
                    if (ice_matl[m]) {
                    rho_micro[m][c] =
                        ice_matl[m]->getEOS()->computeRhoMicro(press_new[c], gamma[m][c],
                            cv[m][c], Temp[m][c], rho_micro[m][c]);

                    //cerr << "rho_micro computed of m " << m << " of cell " << c << " " << rho_micro[m][c] << endl;

                    double div = 1. / rho_micro[m][c];

                    // - updated volume fractions
                    //vol_frac[m][c] = Porosity_CC[c] * rho_CC[m][c] * div;
                    vol_frac[m][c] = rho_CC[m][c] * div;
                    }

                    if (mpm_matl[m]) {
                        rho_micro[m][c] = mpm_matl[m]->getInitialDensity();
                        //vol_frac[m][c] = 0;
                        vol_frac[m][c] = 1 - Porosity_CC[c];
                    }               
                }
                              
                //__________________________________
                // - Test for convergence
                //  If sum of vol_frac_CC ~= vol_frac_not_close_packed then converged
                sum = 0.0;
                for (unsigned int m = 0; m < numALLMatls; m++) {
                    //if (ice_matl[m]) {
                        sum += vol_frac[m][c];    

                        //cerr << "vol frac of material " << m << "at cell " << c << "is " << vol_frac[m][c] << endl;
                   // }
                }
                              
                //cerr << "sum - 1 = " << sum - 1.0 << endl;
                //cerr << "convergence_crit" << convergence_crit << endl;

                if (fabs(sum - 1) < convergence_crit) {
                    converged = true;
                    //__________________________________
                    // Find the speed of sound based on converged solution
                    for (unsigned int m = 0; m < numALLMatls; m++) {

                        if (ice_matl[m]) {
                        ice_matl[m]->getEOS()->computePressEOS(rho_micro[m][c], gamma[m][c],
                            cv[m][c], Temp[m][c],
                            press_eos[m], dp_drho[m], dp_de[m]);

                        tmp = dp_drho[m]
                            + dp_de[m] * press_eos[m] / (rho_micro[m][c] * rho_micro[m][c]);
                        speedSound_new[m][c] = sqrt(tmp);
                        }

                        if (mpm_matl[m]) {
                            //speedSound_new[m][c] = 30;

                            double dp_drho, c_2;
                            //double press_eq_MPM;  // not used

                            // Pointwise computation of thermodynamic quantities from Murnaghan
                            // Parameter
                            double P0 = 101325;
                            double n = 7.4;
                            double K = 39e-11;
                            double rho0 = 2000;

                            if (rho_micro[m][c] >= rho0) {
                                //press_eq_MPM = P0 + (1. / (n * K)) * (pow(rho_micro[m][c] / rho0, n) - 1.);
                                dp_drho = (1. / (K * rho0)) * pow((rho_micro[m][c] / rho0), n - 1.);
                            }
                            else {
                                //press_eq_MPM = P0 * pow(rho_micro[m][c] / rho0, (1. / (K * P0)));
                                dp_drho = (1. / (K * rho0)) * pow(rho_micro[m][c] / rho0, (1. / (K * P0) - 1.));
                            }

                            c_2 = dp_drho;
                            speedSound_new[m][c] = sqrt(c_2);
                        }
                    }
                }
                              
                // Save iteration data for output in case of crash
                if (ds_EqPress.active()) {
                    EqPress_dbg dbg;
                    dbg.delPress = delPress;
                    dbg.press_new = press_new[c];
                    dbg.sumVolFrac = sum;
                    dbg.count = count;

                    for (unsigned int m = 0; m < numALLMatls; m++) {
                        EqPress_dbgMatl dmatl;
                        dmatl.press_eos = press_eos[m];
                        dmatl.volFrac = vol_frac[m][c];
                        dmatl.rhoMicro = rho_micro[m][c];
                        dmatl.rho_CC = rho_CC[m][c];
                        dmatl.temp_CC = Temp[m][c];
                        dmatl.mat = m;
                        dbg.matl.push_back(dmatl);
                    }
                    dbgEqPress.push_back(dbg);
                }
            }   // end of converged           

            test_max_iter = std::max(test_max_iter, count);

            //__________________________________
            //      BULLET PROOFING
            // ignore BP if a recompute time step has already been requested
            bool rts = new_dw->recomputeTimeStep();

            string message;
            bool allTestsPassed = true;
            if (test_max_iter == d_ice->d_max_iter_equilibration && !rts) {
                allTestsPassed = false;
                message += "Max. iterations reached ";
            }

            for (unsigned int m = 0; m < numALLMatls; m++) {
                //if (ice_matl[m]) {
                    if ((vol_frac[m][c] > 0.0) || (vol_frac[m][c] < 1.0)) {
                        message += " ( vol_frac[m][c] > 0.0 ) ||( vol_frac[m][c] < 1.0) ";
                    }
                //}
            }

            if (fabs(sum - 1.0) > convergence_crit && !rts) {
                allTestsPassed = false;
                message += " sum (volumeFractions) != 1 ";
            }

            if (press_new[c] < 0.0 && !rts) {
                allTestsPassed = false;
                message += " Computed pressure is < 0 ";
            }

            for (unsigned int m = 0; m < numALLMatls; m++) {
                if (ice_matl[m]) {
                    if ((rho_micro[m][c] < 0.0 || vol_frac[m][c] < 0.0) && !rts) {
                        allTestsPassed = false;
                        message += " rho_micro < 0 || vol_frac < 0";
                    }
                }
            }

            if (allTestsPassed != true) {  // throw an exception of there's a problem
                Point pt = patch->getCellPosition(c);

                ostringstream warn;
                warn << "\nICE::ComputeEquilibrationPressure: Cell " << c << " position: " << pt << ", L-" << L_indx << "\n"
                    << message
                    << "\nThis usually means that something much deeper has gone wrong with the simulation. "
                    << "\nCompute equilibration pressure task is rarely the problem. "
                    << "For more debugging information set the environmental variable:  \n"
                    << "   SCI_DEBUG DBG_EqPress:+\n\n";

                warn << "INPUTS: \n";

                for (unsigned int m = 0; m < numALLMatls; m++) {
                    if (ice_matl[m]) {
                        warn << "\n matl: " << m << "\n"
                            << "   rho_CC:     " << rho_CC[m][c] << "\n"
                            << "   Temperature:   " << Temp[m][c] << "\n";
                    }

                }
                if (ds_EqPress.active()) {
                    warn << "\nDetails on iterations " << endl;
                    vector<EqPress_dbg>::iterator dbg_iter;

                    for (dbg_iter = dbgEqPress.begin(); dbg_iter != dbgEqPress.end(); dbg_iter++) {
                        EqPress_dbg& d = *dbg_iter;
                        warn << "Iteration:   " << d.count
                            << "  press_new:   " << d.press_new
                            << "  sumVolFrac:  " << d.sumVolFrac
                            << "  delPress:    " << d.delPress << "\n";

                        for (unsigned int m = 0; m < numALLMatls; m++) {
                            if (ice_matl[m]) {
                                warn << "  matl: " << d.matl[m].mat
                                    << "  press_eos:  " << d.matl[m].press_eos
                                    << "  volFrac:    " << d.matl[m].volFrac
                                    << "  rhoMicro:   " << d.matl[m].rhoMicro
                                    << "  rho_CC:     " << d.matl[m].rho_CC
                                    << "  Temp:       " << d.matl[m].temp_CC << "\n";
                            }
                        }
                    }
                }
                throw InvalidValue(warn.str(), __FILE__, __LINE__);
            }
        } // end of cell interator
        

        //for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
            //IntVector c = *iter;
            //cerr << "press_new[c]" << press_new[c] << endl;
        //}


        cout_norm << "max. iterations in any cell " << test_max_iter <<
            " on patch " << patch->getID() << endl;

        //__________________________________
        // carry rho_cc forward 
        // MPMICE computes rho_CC_new
        // therefore need the machinery here
        for (unsigned int m = 0; m < numALLMatls; m++) {
            if (ice_matl[m]) {
                rho_CC_new[m].copyData(rho_CC[m]);             
            }

            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
             
                if (ice_matl[m]) {
                    //cerr << "ComputePressure " << endl;
                    //cerr << "ICE rho_CC_new " << rho_CC_new[m][c] << endl;
                    //cerr << "ICE rho_micro " << rho_micro[m][c] << endl;
                    //cerr << "ICE rho1_CC " << rho1_CC[m][c] << endl;
                    //cerr << "Porosity_CC in Pressure " << Porosity_CC[c] << endl;
                }
                if (mpm_matl[m]) {
                   //rho1_CC[m][c] = rho_micro[m][c] * (1 - Porosity_CC[c]);
                    //cerr << "MPM rho_micro " << rho_micro[m][c] << endl;
                    //cerr << "MPM rho1_CC " << rho1_CC[m][c] << endl;
                }               
            }
        }
    
        //__________________________________
        // - update Boundary conditions
        customBC_localVars* BC_localVars = scinew customBC_localVars();

        preprocess_CustomBCs("EqPress", old_dw, new_dw, Ilb, patch,
            999,d_ice->d_BC_globalVars, BC_localVars);

        setBC(press_new, rho_micro, placeHolder, d_ice->d_surroundingMatl_indx,
            "rho_micro", "Pressure", patch, m_materialManager, 0, new_dw,
            d_ice->d_BC_globalVars, BC_localVars, isNotInitialTimeStep/*, 0  d_ice->d_with_mpmice2*/);

        delete_CustomBCs(d_ice->d_BC_globalVars, BC_localVars);

        //__________________________________
        // compute sp_vol_CC
        // - Set BCs on rhoMicro. using press_CC 
        // - backout sp_vol_new 
        for (unsigned int m = 0; m < numALLMatls; m++) { //for both MPM and ICE
            
                for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                    IntVector c = *iter;
                    sp_vol_new[m][c] = 1.0 / rho_micro[m][c];        
                }
                
                Material* matl = m_materialManager->getMaterial(m);
                int indx = matl->getDWIndex();

                setSpecificVolBC(sp_vol_new[m], "SpecificVol", false, rho_CC[m], vol_frac[m],
                    patch, m_materialManager, indx);               
        }

        //__________________________________
        //  compute f_theta  
        for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            sumKappa[c] = 0.0;
            for (unsigned int m = 0; m < numALLMatls; m++) { // for ICE only
                //if (ice_matl[m]) {
                  //  kappa[m][c] = sp_vol_new[m][c] / (speedSound_new[m][c] * speedSound_new[m][c]);
                 //   sumKappa[c] += vol_frac[m][c] * kappa[m][c];
                    //cerr << "kappa " << kappa[m][c] << endl;
                    //cerr << "sumKappa " << sumKappa[c] << endl;
               // }

               // if (mpm_matl[m]) {
               //     kappa[m][c] = 0.000001;
               // }

                kappa[m][c] = sp_vol_new[m][c] / (speedSound_new[m][c] * speedSound_new[m][c]);
                sumKappa[c] += vol_frac[m][c] * kappa[m][c];

            }

            for (unsigned int m = 0; m < numALLMatls; m++) { // for ICE only
                //if (ice_matl[m]) {
                    f_theta[m][c] = vol_frac[m][c] * kappa[m][c] / sumKappa[c];
                    //cerr << "ICE f_theta " << f_theta[m][c] << endl;
                //}
            }
        }
    }  // patch loop
}

//This task is optional for multi-levelMPMICE2
//
void MPMICE2::scheduleCoarsenCC_0(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* mpm_matls)
{
    printSchedule(patches, cout_doing, "MPMICE2::scheduleCoarsenCC_0");

    bool modifies = false;

    scheduleCoarsenVariableCC(sched, patches, mpm_matls, MIlb->cMassLabel,
        1.9531e-15, modifies, "sum");

    scheduleCoarsenVariableCC(sched, patches, mpm_matls, MIlb->temp_CCLabel,
        0., modifies, "massWeighted");

    scheduleCoarsenVariableCC(sched, patches, mpm_matls, MIlb->vel_CCLabel,
        Vector(0, 0, 0), modifies, "massWeighted");

    scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->sp_vol_CCLabel,
        0.8479864471, modifies, "massWeighted");

    scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->rho_CCLabel,
        1.e-12, modifies, "std");
}

//______________________________________________________________________
//
void MPMICE2::scheduleCoarsenNCMass(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* mpm_matls)
{
    printSchedule(patches, cout_doing, "MPMICE2::scheduleCoarsenNCMass");

    bool modifies = false;

    scheduleCoarsenVariableNC(sched, patches, mpm_matls, Mlb->gMassLabel,
        1.e-200, modifies, "sum");
}

//______________________________________________________________________
// 

void MPMICE2::scheduleComputeVelICE_FC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* ice_matls,
    const MaterialSubset* press_matl,
    const MaterialSet* all_matls)
{
    int levelIndex = getLevel(patches)->getIndex();
    Task* t = 0;

    cout_doing << d_myworld->myRank() << " MPMICE2::scheduleComputeVelICE_FC"
        << "\t\t\t\t\tL-" << levelIndex << endl;

    t = scinew Task("MPMICE2::computeVelICE_FC",
        this, &MPMICE2::computeVelICE_FC);

    Ghost::GhostType  gac = Ghost::AroundCells;
    Task::MaterialDomainSpec oims = Task::OutOfDomain;  //outside of ice matlSet.
    t->requires(Task::OldDW, Ilb->delTLabel, getLevel(patches));
    t->requires(Task::NewDW, Ilb->press_equil_CCLabel, press_matl, oims, gac, 1);
    t->requires(Task::NewDW, Ilb->sp_vol_CCLabel, ice_matls, gac, 1);
    t->requires(Task::OldDW, Ilb->rho_CCLabel, ice_matls, gac, 1);
    t->requires(Task::OldDW, Ilb->vel_CCLabel, ice_matls, gac, 1);

    t->computes(Ilb->uvel_FCLabel, ice_matls);
    t->computes(Ilb->vvel_FCLabel, ice_matls);
    t->computes(Ilb->wvel_FCLabel, ice_matls);
    t->computes(Ilb->grad_P_XFCLabel);
    t->computes(Ilb->grad_P_YFCLabel);
    t->computes(Ilb->grad_P_ZFCLabel);
    sched->addTask(t, patches, all_matls);
}

//______________________________________________________________________
//                       
void MPMICE2::computeVelICE_FC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    const Level* level = getLevel(patches);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing MPMICE2::computeVelICE_FCVel");

        unsigned int numICE_Matls = m_materialManager->getNumMatls("ICE");
        Vector dx = patch->dCell();
        Vector gravity = getGravity();

        cerr << "gravity " << gravity << endl;

        constCCVariable<double> press_CC;
        Ghost::GhostType  gac = Ghost::AroundCells;
        new_dw->get(press_CC, Ilb->press_equil_CCLabel, 0, patch, gac, 1);

        delt_vartype delT;
        old_dw->get(delT, Ilb->delTLabel, level);

        // Compute the face centered velocities
        for (unsigned int m = 0; m < numICE_Matls; m++) {

            ICEMaterial* ice_matl = (ICEMaterial*)m_materialManager->getMaterial("ICE", m);
            int indx = ice_matl->getDWIndex();

            constCCVariable<double> rho_CC, sp_vol_CC;
            constCCVariable<Vector> vel_CC;

            new_dw->get(rho_CC, Ilb->rho_CCLabel, indx, patch, gac, 1);
            old_dw->get(vel_CC, Ilb->vel_CCLabel, indx, patch, gac, 1);
            new_dw->get(sp_vol_CC, Ilb->sp_vol_CCLabel, indx, patch, gac, 1);

            SFCXVariable<double> uvel_FC, grad_P_XFC;
            SFCYVariable<double> vvel_FC, grad_P_YFC;
            SFCZVariable<double> wvel_FC, grad_P_ZFC;

            new_dw->allocateAndPut(uvel_FC, Ilb->uvel_FCLabel, indx, patch);
            new_dw->allocateAndPut(vvel_FC, Ilb->vvel_FCLabel, indx, patch);
            new_dw->allocateAndPut(wvel_FC, Ilb->wvel_FCLabel, indx, patch);
            // debugging variables
            new_dw->allocateAndPut(grad_P_XFC, Ilb->grad_P_XFCLabel, indx, patch);
            new_dw->allocateAndPut(grad_P_YFC, Ilb->grad_P_YFCLabel, indx, patch);
            new_dw->allocateAndPut(grad_P_ZFC, Ilb->grad_P_ZFCLabel, indx, patch);

            IntVector lowIndex(patch->getExtraSFCXLowIndex());
            uvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCXHighIndex());
            vvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCYHighIndex());
            wvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCZHighIndex());

            grad_P_XFC.initialize(0.0);
            grad_P_YFC.initialize(0.0);
            grad_P_ZFC.initialize(0.0);

            vector<IntVector> adj_offset(3);
            adj_offset[0] = IntVector(-1, 0, 0);    // X faces
            adj_offset[1] = IntVector(0, -1, 0);    // Y faces
            adj_offset[2] = IntVector(0, 0, -1);   // Z faces     

            CellIterator XFC_iterator = patch->getSFCXIterator();
            CellIterator YFC_iterator = patch->getSFCYIterator();
            CellIterator ZFC_iterator = patch->getSFCZIterator();

            //__________________________________
            //  Compute vel_FC for each face
            computeVelICEFace<SFCXVariable<double> >(0, XFC_iterator,
                adj_offset[0], dx[0], delT, gravity[0],
                rho_CC, sp_vol_CC, vel_CC, press_CC,
                uvel_FC, grad_P_XFC);

            computeVelICEFace<SFCYVariable<double> >(1, YFC_iterator,
                adj_offset[1], dx[1], delT, gravity[1],
                rho_CC, sp_vol_CC, vel_CC, press_CC,
                vvel_FC, grad_P_YFC);

            computeVelICEFace<SFCZVariable<double> >(2, ZFC_iterator,
                adj_offset[2], dx[2], delT, gravity[2],
                rho_CC, sp_vol_CC, vel_CC, press_CC,
                wvel_FC, grad_P_ZFC);

            //__________________________________
            // (*)vel_FC BC are updated in 
            // ICE::addExchangeContributionToFCVel()

        } // matls loop
    }  // patch loop
}

/* _____________________________________________________________________
 Function~  MPMICE2::computeFaceCenteredVelocities--
 Purpose~   compute the face centered velocities minus the exchange
            contribution.
_____________________________________________________________________*/
template<class T> void MPMICE2::computeVelICEFace(int dir,
    CellIterator it,
    IntVector adj_offset,
    double dx,
    double delT, double gravity,
    constCCVariable<double>& rho_CC,
    constCCVariable<double>& sp_vol_CC,
    constCCVariable<Vector>& vel_CC,
    constCCVariable<double>& press_CC,
    T& vel_FC,
    T& grad_P_FC)
{
    double inv_dx = 1.0 / dx;

    for (; !it.done(); it++) {
        IntVector R = *it;
        IntVector L = R + adj_offset;

        double rho1_FC = rho_CC[L] + rho_CC[R];

        #if SCI_ASSERTION_LEVEL >=2
                if (rho_FC <= 0.0) {
                    cout << d_myworld->myRank() << " rho_fc <= 0: " << rho_FC << " with L= " << L << " ("
                        << rho_CC[L] << ") R= " << R << " (" << rho_CC[R] << ")\n";
                }
        #endif
        ASSERT(rho_FC > 0.0);

        //__________________________________
        // interpolation to the face
        double term1 = (rho_CC[L] * vel_CC[L][dir] +
            rho_CC[R] * vel_CC[R][dir]) / (rho1_FC);
        //__________________________________
        // pressure term           
        double sp_vol_brack = 2. * (sp_vol_CC[L] * sp_vol_CC[R]) /
            (sp_vol_CC[L] + sp_vol_CC[R]);

        grad_P_FC[R] = (press_CC[R] - press_CC[L]) * inv_dx;
        double term2 = delT * sp_vol_brack * grad_P_FC[R];

        //__________________________________
        // gravity term
        double term3 = delT * gravity;

        vel_FC[R] = term1 - term2 + term3;

        //cerr << "vel_FC at " << R << " is " << vel_FC[R] << endl;
        //cerr << "term1 " << term1 << endl;
        cerr << "term2 ice " << term2 << endl;
       // cerr << "term3 " << term3 << endl;
        //cerr << "delT " << delT << endl;
       // cerr << "gravity " << gravity << endl;
    }
}

void MPMICE2::scheduleComputeVelMPM_FC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* mpm_matls,
    const MaterialSubset* press_matl,
    const MaterialSet* all_matls)
{
    int levelIndex = getLevel(patches)->getIndex();
    Task* t = 0;

    cout_doing << d_myworld->myRank() << " MPMICE2::scheduleComputeVelMPM_FC"
        << "\t\t\t\t\tL-" << levelIndex << endl;

    t = scinew Task("MPMICE2::computeVelMPM_FC",
        this, &MPMICE2::computeVelMPM_FC);

    Ghost::GhostType  gac = Ghost::AroundCells;
    Task::MaterialDomainSpec oims = Task::OutOfDomain;  //outside of ice matlSet.
    //const MaterialSubset* mss = mpm_matls->getUnion();

    t->requires(Task::OldDW, Ilb->delTLabel, getLevel(patches));
    t->requires(Task::NewDW, Ilb->press_equil_CCLabel, press_matl, oims, gac, 1);
    t->requires(Task::NewDW, Ilb->rho_CCLabel, mpm_matls, gac, 1);
    t->requires(Task::NewDW, MIlb->vel_CCLabel, mpm_matls, gac, 1);
    //t->requires(Task::NewDW, MIlb->stress_CCLabel, mpm_matls, gac, 1);

    // Hacking
    t->requires(Task::NewDW, MIlb->stressX_CCLabel, mpm_matls, gac, 1);
    t->requires(Task::NewDW, MIlb->stressY_CCLabel, mpm_matls, gac, 1);
    t->requires(Task::NewDW, MIlb->stressZ_CCLabel, mpm_matls, gac, 1);

    t->computes(Ilb->uvel_FCLabel, mpm_matls);
    t->computes(Ilb->vvel_FCLabel, mpm_matls);
    t->computes(Ilb->wvel_FCLabel, mpm_matls);
    t->computes(Ilb->grad_P_XFCLabel);
    t->computes(Ilb->grad_P_YFCLabel);
    t->computes(Ilb->grad_P_ZFCLabel);
    t->computes(MIlb->grad_stress_XFCLabel);
    t->computes(MIlb->grad_stress_YFCLabel);
    t->computes(MIlb->grad_stress_ZFCLabel);
    sched->addTask(t, patches, all_matls);
}

//______________________________________________________________________
//                       
void MPMICE2::computeVelMPM_FC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    const Level* level = getLevel(patches);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing MPMICE2::computeVelMPM_FCVel");

        unsigned int numMPM_matls = m_materialManager->getNumMatls("MPM");
        Vector dx = patch->dCell();
        Vector gravity = getGravity();

        constCCVariable<double> press_CC;
        Ghost::GhostType  gac = Ghost::AroundCells;
        new_dw->get(press_CC, Ilb->press_equil_CCLabel, 0, patch, gac, 1);

        delt_vartype delT;
        old_dw->get(delT, Ilb->delTLabel, level);

        // Compute the face centered velocities
        for (unsigned int m = 0; m < numMPM_matls; m++) {

            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();

            double rho_init = mpm_matl->getInitialDensity();

            constCCVariable<double> rho_CC, sp_vol_CC;
            constCCVariable<Vector> vel_CC;
            //constCCVariable<Matrix3> stress_CC;

            new_dw->get(rho_CC, Ilb->rho_CCLabel, indx, patch, gac, 1);
            new_dw->get(vel_CC, MIlb->vel_CCLabel, indx, patch, gac, 1);
            //new_dw->get(stress_CC, MIlb->stress_CCLabel, indx, patch, gac, 1);

            constCCVariable<double> stressX_CC, stressY_CC, stressZ_CC;
            new_dw->get(stressX_CC, MIlb->stressX_CCLabel, indx, patch, gac, 1);
            new_dw->get(stressY_CC, MIlb->stressY_CCLabel, indx, patch, gac, 1);
            new_dw->get(stressZ_CC, MIlb->stressZ_CCLabel, indx, patch, gac, 1);

            SFCXVariable<double> uvel_FC, grad_stress_XFC, grad_P_XFC;
            SFCYVariable<double> vvel_FC, grad_stress_YFC, grad_P_YFC;
            SFCZVariable<double> wvel_FC, grad_stress_ZFC, grad_P_ZFC;

            new_dw->allocateAndPut(uvel_FC, Ilb->uvel_FCLabel, indx, patch);
            new_dw->allocateAndPut(vvel_FC, Ilb->vvel_FCLabel, indx, patch);
            new_dw->allocateAndPut(wvel_FC, Ilb->wvel_FCLabel, indx, patch);

            // debugging variables
            new_dw->allocateAndPut(grad_stress_XFC, MIlb->grad_stress_XFCLabel, indx, patch);
            new_dw->allocateAndPut(grad_stress_YFC, MIlb->grad_stress_YFCLabel, indx, patch);
            new_dw->allocateAndPut(grad_stress_ZFC, MIlb->grad_stress_ZFCLabel, indx, patch);

            new_dw->allocateAndPut(grad_P_XFC, Ilb->grad_P_XFCLabel, indx, patch);
            new_dw->allocateAndPut(grad_P_YFC, Ilb->grad_P_YFCLabel, indx, patch);
            new_dw->allocateAndPut(grad_P_ZFC, Ilb->grad_P_ZFCLabel, indx, patch);

            IntVector lowIndex(patch->getExtraSFCXLowIndex());
            uvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCXHighIndex());
            vvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCYHighIndex());
            wvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCZHighIndex());

            grad_stress_XFC.initialize(0.0);
            grad_stress_YFC.initialize(0.0);
            grad_stress_ZFC.initialize(0.0);

            grad_P_XFC.initialize(0.0);
            grad_P_YFC.initialize(0.0);
            grad_P_ZFC.initialize(0.0);

            vector<IntVector> adj_offset(3);
            adj_offset[0] = IntVector(-1, 0, 0);    // X faces
            adj_offset[1] = IntVector(0, -1, 0);    // Y faces
            adj_offset[2] = IntVector(0, 0, -1);   // Z faces     

            CellIterator XFC_iterator = patch->getSFCXIterator();
            CellIterator YFC_iterator = patch->getSFCYIterator();
            CellIterator ZFC_iterator = patch->getSFCZIterator();

            //__________________________________
            //  Compute vel_FC for each face
            computeVelMPMFace<SFCXVariable<double> >(0, XFC_iterator,
                adj_offset[0], dx[0], delT, gravity[0], rho_init,
                rho_CC, stressX_CC, vel_CC, press_CC,
                uvel_FC, grad_P_XFC, grad_stress_XFC);

            computeVelMPMFace<SFCYVariable<double> >(1, YFC_iterator,
                adj_offset[1], dx[1], delT, gravity[1], rho_init,
                rho_CC, stressY_CC, vel_CC, press_CC,
                vvel_FC, grad_P_YFC, grad_stress_YFC);

            computeVelMPMFace<SFCZVariable<double> >(2, ZFC_iterator,
                adj_offset[2], dx[2], delT, gravity[2], rho_init,
                rho_CC, stressZ_CC, vel_CC, press_CC,
                wvel_FC, grad_P_ZFC, grad_stress_ZFC);

            //__________________________________
            // (*)vel_FC BC are updated in 
            // ICE::addExchangeContributionToFCVel()

        } // matls loop
    }  // patch loop
}


/*___________________________________________________________________
 Function~  MPMICE2::computeFaceCenteredVelocities--
 Purpose~   compute the face centered velocities minus the exchange
            contribution.
_____________________________________________________________________*/
template<class T> void MPMICE2::computeVelMPMFace(int dir,
    CellIterator it,
    IntVector adj_offset,
    double dx,
    double delT, double gravity, double rho_init,
    constCCVariable<double>& rho_CC,
    constCCVariable<double>& stress_CC,
    constCCVariable<Vector>& vel_CC,
    constCCVariable<double>& press_CC,
    T& vel_FC,
    T& grad_P_FC,
    T& grad_stress_FC)
{
    double inv_dx = 1.0 / dx;

    for (; !it.done(); it++) {
        IntVector R = *it;
        IntVector L = R + adj_offset;

        double rho_FC = rho_CC[L] + rho_CC[R];

#if SCI_ASSERTION_LEVEL >=2
        if (rho_FC <= 0.0) {
            cout << d_myworld->myRank() << " rho_fc <= 0: " << rho_FC << " with L= " << L << " ("
                << rho_CC[L] << ") R= " << R << " (" << rho_CC[R] << ")\n";
        }
#endif
        ASSERT(rho_FC > 0.0);

        //__________________________________
        // interpolation to the face
        double term1 = (rho_CC[L] * vel_CC[L][dir] +
            rho_CC[R] * vel_CC[R][dir]) / (rho_FC);
        //__________________________________
        // pressure term           

        grad_stress_FC[R] = (stress_CC[R] - stress_CC[L]) * inv_dx;
        double term2 = delT * grad_stress_FC[R] / rho_FC;

        grad_P_FC[R] = (press_CC[R] - press_CC[L]) * inv_dx;
        double term3 = delT * grad_P_FC[R] / rho_init;

        //__________________________________
        // gravity term
        double term4 = delT * gravity;

        vel_FC[R] = term1 + term2 - term3 + term4;

        //cerr << "MPM term 1 " << term1 << endl;
        //cerr << "MPM term 2 " << term2 << endl;
        cerr << "MPM term 3 " << term3 << endl;
        //cerr << "MPM term 4 " << term4 << endl;
    }
}

//______________________________________________________________________
//
void MPMICE2::scheduleInterpolatePressCCToPressNC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* press_matl,
    const MaterialSet* matls)
{
    const Level* level = getLevel(patches);
    int L_indx = level->getIndex();
    if (!d_mpm->flags->doMPMOnLevel(L_indx, level->getGrid()->numLevels()))
        return;

    printSchedule(patches, cout_doing, "MPMICE2::scheduleInterpolatePressCCToPressNC");

    Task* t = scinew Task("MPMICE2::interpolatePressCCToPressNC",
        this, &MPMICE2::interpolatePressCCToPressNC);

    Ghost::GhostType  gac = Ghost::AroundCells;
    t->requires(Task::NewDW, Ilb->press_CCLabel, press_matl, gac, 1);
    t->computes(MIlb->press_NCLabel, press_matl);

    sched->addTask(t, patches, matls);
}

//______________________________________________________________________
//
void MPMICE2::interpolatePressCCToPressNC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse*,
    DataWarehouse* new_dw)
{
    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing interpolatePressCCToPressNC");

        constCCVariable<double> pressCC;
        NCVariable<double> pressNC;

        Ghost::GhostType  gac = Ghost::AroundCells;
        new_dw->get(pressCC, Ilb->press_CCLabel, 0, patch, gac, 1);
        new_dw->allocateAndPut(pressNC, MIlb->press_NCLabel, 0, patch);
        pressNC.initialize(0.0);

        IntVector cIdx[8];
        // Interpolate CC pressure to nodes
        for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
            patch->findCellsFromNode(*iter, cIdx);
            for (int in = 0; in < 8; in++) {
                pressNC[*iter] += .125 * pressCC[cIdx[in]];
            }
        }

        // Apply grid boundary conditions to the pressure before storing the data
        string inter_type = d_mpm->flags->d_interpolator_type;
        MPMBoundCond bc;
        bc.setBoundaryCondition(patch, 0, "Pressure", pressNC, inter_type);
    }
}

//______________________________________________________________________
//
void MPMICE2::scheduleInterpolatePAndGradP(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* press_matl,
    const MaterialSubset* /*one_matl*/,
    const MaterialSubset* mpm_matl,
    const MaterialSet* all_matls)
{
    if (!d_mpm->flags->doMPMOnLevel(getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels()))
        return;

    printSchedule(patches, cout_doing, "MPMICE2::scheduleInterpolatePAndGradP");

    Task* t = scinew Task("MPMICE2::interpolatePAndGradP",
        this, &MPMICE2::interpolatePAndGradP);
    Ghost::GhostType  gac = Ghost::AroundCells;

    t->requires(Task::NewDW, MIlb->press_NCLabel, press_matl, gac, NGN);
    t->requires(Task::NewDW, MIlb->cMassLabel, mpm_matl, gac, 1);
    t->requires(Task::OldDW, Mlb->pXLabel, mpm_matl, Ghost::None);
    t->requires(Task::NewDW, Mlb->pCurSizeLabel, mpm_matl, Ghost::None);

    t->computes(Mlb->pPressureLabel, mpm_matl);
    sched->addTask(t, patches, all_matls);
}

//______________________________________________________________________
//
void MPMICE2::interpolatePAndGradP(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing interpolatePressureToParticles");

        ParticleInterpolator* interpolator = d_mpm->flags->d_interpolator->clone(patch);
        vector<IntVector> ni(interpolator->size());
        vector<double> S(interpolator->size());

        //double p_ref = d_ice->getRefPress();
        constNCVariable<double>   pressNC;
        Ghost::GhostType  gac = Ghost::AroundCells;
        new_dw->get(pressNC, MIlb->press_NCLabel, 0, patch, gac, NGN);

        for (unsigned int m = 0; m < m_materialManager->getNumMatls("MPM"); m++) {
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();

            ParticleSubset* pset = old_dw->getParticleSubset(indx, patch);
            ParticleVariable<double> pPressure;
            constParticleVariable<Point> px;
            constParticleVariable<Matrix3> psize;
            constParticleVariable<Matrix3> deformationGradient;
            new_dw->get(psize, Mlb->pCurSizeLabel, pset);
            old_dw->get(px, Mlb->pXLabel, pset);
            new_dw->allocateAndPut(pPressure, Mlb->pPressureLabel, pset);

            //__________________________________
            // Interpolate NC pressure to particles
            for (ParticleSubset::iterator iter = pset->begin();
                iter != pset->end(); iter++) {
                particleIndex idx = *iter;
                double press = 0.;

                // Get the node indices that surround the cell
                int NN = interpolator->findCellAndWeights(px[idx], ni, S, psize[idx]);

                for (int k = 0; k < NN; k++) {
                    press += pressNC[ni[k]] * S[k];
                }
                pPressure[idx] = press;               
                //pPressure[idx] = press - p_ref;
            }
        }  // numMPMMatls
        delete interpolator;
    } //patches
}

//______________________________________________________________________
//
void MPMICE2::scheduleComputeLagrangianValuesMPM(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* one_matl,
    const MaterialSet* mpm_matls)
{
    if (d_mpm->flags->doMPMOnLevel(getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {

        printSchedule(patches, cout_doing, "MPMICE2::scheduleComputeLagrangianValuesMPM");

        Task* t = scinew Task("MPMICE2::computeLagrangianValuesMPM",
            this, &MPMICE2::computeLagrangianValuesMPM);

        const MaterialSubset* mss = mpm_matls->getUnion();
        Ghost::GhostType  gac = Ghost::AroundCells;
        Ghost::GhostType  gn = Ghost::None;
        t->requires(Task::NewDW, Mlb->gVelocityStarLabel, mss, gac, 1);
        t->requires(Task::NewDW, Mlb->gMassLabel, gac, 1);
        t->requires(Task::NewDW, Mlb->gTemperatureStarLabel, gac, 1);
        t->requires(Task::OldDW, Mlb->NC_CCweightLabel, one_matl, gac, 1);
        t->requires(Task::NewDW, MIlb->cMassLabel, gn);
        t->requires(Task::NewDW, Ilb->mom_source_CCLabel, gn);

        t->requires(Task::OldDW, Ilb->timeStepLabel);

        if (d_ice->d_models.size() > 0 && !do_mlmpmice2) {
            t->requires(Task::NewDW, Ilb->modelMass_srcLabel, gn);
            t->requires(Task::NewDW, Ilb->modelMom_srcLabel, gn);
            t->requires(Task::NewDW, Ilb->modelEng_srcLabel, gn);
        }

        t->computes(Ilb->mass_L_CCLabel);
        t->computes(Ilb->mom_L_CCLabel);
        t->computes(Ilb->int_eng_L_CCLabel);

        sched->addTask(t, patches, mpm_matls);
    }
}

//______________________________________________________________________
//
void MPMICE2::computeLagrangianValuesMPM(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, VarLabel::find(timeStep_name));

    bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing computeLagrangianValuesMPM");

        unsigned int numMatls = m_materialManager->getNumMatls("MPM");
        Vector dx = patch->dCell();
        double cellVol = dx.x() * dx.y() * dx.z();
        double very_small_mass = d_TINY_RHO * cellVol;
        Ghost::GhostType  gn = Ghost::None;
        Ghost::GhostType  gac = Ghost::AroundCells;

        constNCVariable<double> NC_CCweight;
        old_dw->get(NC_CCweight, Mlb->NC_CCweightLabel, 0, patch, gac, 1);
        for (unsigned int m = 0; m < numMatls; m++) {
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();

            // Create arrays for the grid data
            constNCVariable<double> gmass, gvolume, gtempstar;
            constNCVariable<Vector> gvelocity;
            CCVariable<Vector> cmomentum;
            CCVariable<double> int_eng_L, mass_L;
            constCCVariable<double> cmass;
            constCCVariable<Vector> mom_source;
            new_dw->get(gmass, Mlb->gMassLabel, indx, patch, gac, 1);
            new_dw->get(gvelocity, Mlb->gVelocityStarLabel, indx, patch, gac, 1);
            new_dw->get(gtempstar, Mlb->gTemperatureStarLabel, indx, patch, gac, 1);
            new_dw->get(cmass, MIlb->cMassLabel, indx, patch, gn, 0);
            new_dw->get(mom_source, Ilb->mom_source_CCLabel, indx, patch, gn, 0);

            new_dw->allocateAndPut(mass_L, Ilb->mass_L_CCLabel, indx, patch);
            new_dw->allocateAndPut(cmomentum, Ilb->mom_L_CCLabel, indx, patch);
            new_dw->allocateAndPut(int_eng_L, Ilb->int_eng_L_CCLabel, indx, patch);

            cmomentum.initialize(Vector(0.0, 0.0, 0.0));
            int_eng_L.initialize(0.);
            mass_L.initialize(0.);
            double cv = mpm_matl->getSpecificHeat();

            IntVector nodeIdx[8];

            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
                iter++) {
                IntVector c = *iter;
                patch->findNodesFromCell(c, nodeIdx);
                double int_eng_L_mpm = 0.0;
                Vector cmomentum_mpm = Vector(0.0, 0.0, 0.0);

                for (int in = 0; in < 8; in++) {
                    double NC_CCw_mass = NC_CCweight[nodeIdx[in]] * gmass[nodeIdx[in]];
                    cmomentum_mpm += gvelocity[nodeIdx[in]] * NC_CCw_mass;
                    int_eng_L_mpm += gtempstar[nodeIdx[in]] * cv * NC_CCw_mass;
                }

                // Add pore water pressure term  (1-n) * gradP
                cmomentum_mpm -= mom_source[c];
                cmomentum[c] = cmomentum_mpm;
                int_eng_L[c] = int_eng_L_mpm;
                //cerr << "MPM cmomentum" << cmomentum[c] << endl;
                //cerr << "MPM int_eng_L" << int_eng_L[c] << endl;
            }
            //__________________________________
            //  NO REACTION
            if (d_ice->d_models.size() == 0 || do_mlmpmice2) {
                for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
                    iter++) {
                    IntVector c = *iter;
                    mass_L[c] = cmass[c];
                    //cerr << "MPM mass_L" << mass_L[c] << endl;
                }
            }
            //__________________________________
            //   M O D E L   B A S E D   E X C H A N G E (NOT YET CONSIDERED IN MPMICE2)
            // The reaction can't completely eliminate 
            //  all the mass, momentum and internal E.
            // If it does then we'll get erroneous vel,
            // and temps in CCMomExchange.  If the mass
            // goes to min_mass then cmomentum and int_eng_L
            // need to be scaled by min_mass to avoid inf temp and vel_CC
            // in 
            if (d_ice->d_models.size() > 0 && !do_mlmpmice2) {
                constCCVariable<double> modelMass_src;
                constCCVariable<double> modelEng_src;
                constCCVariable<Vector> modelMom_src;
                new_dw->get(modelMass_src, Ilb->modelMass_srcLabel, indx, patch, gn, 0);
                new_dw->get(modelMom_src, Ilb->modelMom_srcLabel, indx, patch, gn, 0);
                new_dw->get(modelEng_src, Ilb->modelEng_srcLabel, indx, patch, gn, 0);

                for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
                    iter++) {
                    IntVector c = *iter;
                    //  must have a minimum mass
                    double min_mass = very_small_mass;
                    double inv_cmass = 1.0 / cmass[c];
                    mass_L[c] = std::max((cmass[c] + modelMass_src[c]), min_mass);

                    //  must have a minimum momentum 
                    for (int dir = 0; dir < 3; dir++) {  //loop over all three directons
                        double min_mom_L = min_mass * cmomentum[c][dir] * inv_cmass;
                        double mom_L_tmp = cmomentum[c][dir] + modelMom_src[c][dir];

                        // Preserve the original sign on momemtum     
                        // Use d_SMALL_NUMs to avoid nans when mom_L_temp = 0.0
                        double plus_minus_one = (mom_L_tmp + d_SMALL_NUM) /
                            (fabs(mom_L_tmp + d_SMALL_NUM));

                        mom_L_tmp = (mom_L_tmp / mass_L[c]) * (cmass[c] + modelMass_src[c]);

                        cmomentum[c][dir] = plus_minus_one *
                            std::max(fabs(mom_L_tmp), fabs(min_mom_L));
                    }
                    // must have a minimum int_eng   
                    double min_int_eng = min_mass * int_eng_L[c] * inv_cmass;
                    int_eng_L[c] = (int_eng_L[c] / mass_L[c]) * (cmass[c] + modelMass_src[c]);
                    int_eng_L[c] = std::max((int_eng_L[c] + modelEng_src[c]), min_int_eng);
                }
            }  // if(model.size() >0)      

             //__________________________________
             //  Set Boundary conditions
            setBC(cmomentum, "set_if_sym_BC", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(int_eng_L, "set_if_sym_BC", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);

            //---- B U L L E T   P R O O F I N G------
            // ignore BP if recompute time step has already been requested
            IntVector neg_cell;
            ostringstream warn;
            bool rts = new_dw->recomputeTimeStep();

            if (d_testForNegTemps_mpm) {
                if (!areAllValuesPositive(int_eng_L, neg_cell) && !rts) {
                    int L = getLevel(patches)->getIndex();
                    warn << "ERROR MPMICE2:(" << L << "):computeLagrangianValuesMPM, mat "
                        << indx << " cell "
                        << neg_cell << " int_eng_L_CC " << int_eng_L[neg_cell] << "\n ";
                    throw InvalidValue(warn.str(), __FILE__, __LINE__);
                }
            }
        }  //numMatls
    }  //patches
}

/* _____________________________________________________________________
 Function~  MPMICE2:: scheduleComputeLagrangianSpecificVolume--
_____________________________________________________________________*/
void MPMICE2::scheduleComputeLagrangianSpecificVolume(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* ice_matls,
    const MaterialSubset* mpm_matls,
    const MaterialSubset* press_matl,
    const MaterialSet* matls)
{
    int levelIndex = getLevel(patches)->getIndex();
    Task* t = 0;
    cout_doing << d_myworld->myRank() << " MPMICE2::scheduleComputeLagrangianSpecificVolume"
        << "\t\t\tL-" << levelIndex << endl;
    t = scinew Task("MPMICE2::computeLagrangianSpecificVolume",
        this, &MPMICE2::computeLagrangianSpecificVolume);

    Ghost::GhostType  gn = Ghost::None;
    Ghost::GhostType  gac = Ghost::AroundCells;
    Task::MaterialDomainSpec oims = Task::OutOfDomain;  //outside of ice matlSet.

    t->requires(Task::NewDW, Ilb->uvel_FCMELabel, gac, 2);
    t->requires(Task::NewDW, Ilb->vvel_FCMELabel, gac, 2);
    t->requires(Task::NewDW, Ilb->wvel_FCMELabel, gac, 2);

    t->requires(Task::OldDW, Mlb->delTLabel, getLevel(patches));
    t->requires(Task::NewDW, Ilb->rho_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->sp_vol_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->Tdot_CCLabel, gn);

    t->requires(Task::NewDW, Ilb->f_theta_CCLabel, ice_matls, gn);
    t->requires(Task::NewDW, Ilb->compressibilityLabel, gn);
    t->requires(Task::NewDW, Ilb->vol_frac_CCLabel, gac, 1);

    t->requires(Task::OldDW, Ilb->temp_CCLabel, ice_matls, gn);
    t->requires(Task::NewDW, Ilb->specific_heatLabel, ice_matls, gn);
    if (mpm_matls){
        t->requires(Task::NewDW, Ilb->temp_CCLabel, mpm_matls, gn);
        t->requires(Task::NewDW, Ilb->VolumeFraction_CCLabel, mpm_matls, gn);
    }

    t->requires(Task::NewDW, Ilb->delP_DilatateLabel, press_matl, oims, gn);
    t->requires(Task::NewDW, Ilb->press_CCLabel, press_matl, oims, gn);
    t->requires(Task::NewDW, Ilb->sumKappaLabel, press_matl, oims, gn);
    //t->requires(Task::NewDW, Ilb->Porosity_CCLabel, press_matl, oims, gn);

    if (d_ice->d_models.size() > 0) {
        t->requires(Task::NewDW, Ilb->modelVol_srcLabel, gn);
    }

    t->computes(Ilb->sp_vol_L_CCLabel);
    t->computes(Ilb->sp_vol_src_CCLabel);

    //t->computes(Ilb->vol_fracX_FCLabel);
    //t->computes(Ilb->vol_fracY_FCLabel);
    //t->computes(Ilb->vol_fracZ_FCLabel);

    t->computes(VarLabel::find(abortTimeStep_name));
    t->computes(VarLabel::find(recomputeTimeStep_name));

    sched->addTask(t, patches, matls);
}

/* _____________________________________________________________________
 Function~  MPMICE2::computeLagrangianSpecificVolume--
 _____________________________________________________________________  */
void MPMICE2::computeLagrangianSpecificVolume(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset* /*matls*/,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    const Level* level = getLevel(patches);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);

        printTask(patches, patch, cout_doing, "Doing MPMICE2::computeLagrangianSpecificVolume");

        delt_vartype delT;
        old_dw->get(delT, Mlb->delTLabel, level);

        Advector* advector = d_ice->d_advector->clone(new_dw, patch, isRegridTimeStep());

        unsigned int numICEMatls = m_materialManager->getNumMatls("ICE");
        unsigned int numMPMMatls = m_materialManager->getNumMatls("MPM");
        unsigned int numALLMatls = numICEMatls + numMPMMatls;

        Vector  dx = patch->dCell();
        double vol = dx.x() * dx.y() * dx.z();
        Ghost::GhostType  gn = Ghost::None;
        Ghost::GhostType  gac = Ghost::AroundCells;

        std::vector<constCCVariable<double> > Tdot(numALLMatls);
        std::vector<constCCVariable<double> > vol_frac(numALLMatls);
        std::vector<constCCVariable<double> > Temp_CC(numALLMatls);
        std::vector<constCCVariable<double> > VolumeFraction_CC(numALLMatls);
        std::vector<CCVariable<double> > alpha(numALLMatls);
        constCCVariable<double> rho_CC, f_theta, sp_vol_CC, cv;
        constCCVariable<double> delP, P;
        //constCCVariable<double> Porosity_CC;
        CCVariable<double> sum_therm_exp;
        constCCVariable<double>sumKappa;
        
        new_dw->allocateTemporary(sum_therm_exp, patch);
        new_dw->get(delP, Ilb->delP_DilatateLabel, 0, patch, gn, 0);
        new_dw->get(P, Ilb->press_CCLabel, 0, patch, gn, 0);
        //new_dw->get(Porosity_CC, Ilb->Porosity_CCLabel, 0, patch, gn, 0);
        sum_therm_exp.initialize(0.);
        new_dw->get(sumKappa, Ilb->sumKappaLabel, 0, patch, gn, 0);

        // Aditional line
        CCVariable<double> termICE, termMPM, q_advectedICE, q_advectedMPM, vol_fracICE;
        new_dw->allocateTemporary(termICE, patch);
        new_dw->allocateTemporary(termMPM, patch);
        new_dw->allocateTemporary(q_advectedICE, patch);
        new_dw->allocateTemporary(q_advectedMPM, patch);
        new_dw->allocateTemporary(vol_fracICE, patch);
        termICE.initialize(0.);
        termMPM.initialize(0.);
        q_advectedICE.initialize(0.);
        q_advectedMPM.initialize(0.);
        vol_fracICE.initialize(0.);

        CCVariable<double> q_advected;
        new_dw->allocateTemporary(q_advected, patch);
        q_advected.initialize(0.);
     
        for (unsigned int m = 0; m < numALLMatls; m++) {
            Material* matl = m_materialManager->getMaterial(m);
            MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
            ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
            int indx = matl->getDWIndex();

            new_dw->get(Tdot[m], Ilb->Tdot_CCLabel, indx, patch, gn, 0);
            new_dw->get(vol_frac[m], Ilb->vol_frac_CCLabel, indx, patch, gac, 1);
            new_dw->allocateTemporary(alpha[m], patch);

            if (ice_matl) {
                old_dw->get(Temp_CC[m], Ilb->temp_CCLabel, indx, patch, gn, 0);
            }
            if (mpm_matl) {
                new_dw->get(Temp_CC[m], Ilb->temp_CCLabel, indx, patch, gn, 0);
                new_dw->get(VolumeFraction_CC[m], Ilb->VolumeFraction_CCLabel, indx, patch, gn, 0);
            }
        }
        
        //__________________________________
        // Sum of thermal expansion and avection of porosity
        for (unsigned int m = 0; m < numALLMatls; m++) {
            Material* matl = m_materialManager->getMaterial(m);
            ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
            //MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);

            int indx = matl->getDWIndex();

            constSFCXVariable<double> uvel_FC;
            constSFCYVariable<double> vvel_FC;
            constSFCZVariable<double> wvel_FC;
            constCCVariable<double> mass_L;
            //constCCVariable<double> kappa;

            Ghost::GhostType  gac = Ghost::AroundCells;
            new_dw->get(uvel_FC, Ilb->uvel_FCMELabel, indx, patch, gac, 2);
            new_dw->get(vvel_FC, Ilb->vvel_FCMELabel, indx, patch, gac, 2);
            new_dw->get(wvel_FC, Ilb->wvel_FCMELabel, indx, patch, gac, 2);
            new_dw->get(mass_L, Ilb->mass_L_CCLabel, indx, patch, gac, 2);
            //new_dw->get(kappa, Ilb->compressibilityLabel, indx, patch, gac, 2);

            /*
            SFCXVariable<double> vol_fracX_FC;
            SFCYVariable<double> vol_fracY_FC;
            SFCZVariable<double> vol_fracZ_FC;

            new_dw->allocateAndPut(vol_fracX_FC, Ilb->vol_fracX_FCLabel, indx, patch);
            new_dw->allocateAndPut(vol_fracY_FC, Ilb->vol_fracY_FCLabel, indx, patch);
            new_dw->allocateAndPut(vol_fracZ_FC, Ilb->vol_fracZ_FCLabel, indx, patch);


            // lowIndex is the same for all vel_FC
            IntVector lowIndex(patch->getExtraSFCXLowIndex());
            double nan = getNan();
            vol_fracX_FC.initialize(nan, lowIndex, patch->getExtraSFCXHighIndex());
            vol_fracY_FC.initialize(nan, lowIndex, patch->getExtraSFCYHighIndex());
            vol_fracZ_FC.initialize(nan, lowIndex, patch->getExtraSFCZHighIndex());
            */

            if (ice_matl) {
                new_dw->get(sp_vol_CC, Ilb->sp_vol_CCLabel, indx, patch, gn, 0);             
                new_dw->get(cv, Ilb->specific_heatLabel, indx, patch, gn, 0);

                for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                    IntVector c = *iter;
                    alpha[m][c] =
                        ice_matl->getEOS()->getAlpha(Temp_CC[m][c], sp_vol_CC[c], P[c], cv[c]);
                    sum_therm_exp[c] += vol_frac[m][c] * alpha[m][c] * Tdot[m][c];
                }
            }

            /*
            // Calculate the divergence of velocity
            //__________________________________            
            // Advection preprocessing
            bool bulletProof_test = true;
            advectVarBasket* varBasket = scinew advectVarBasket();
            varBasket->new_dw = new_dw;
            varBasket->old_dw = old_dw;
            varBasket->indx = indx;
            varBasket->patch = patch;
            varBasket->level = level;
            varBasket->lb = Ilb;
            varBasket->doRefluxing = isAMR(); // always reflux with amr
            varBasket->is_Q_massSpecific = false;
            varBasket->useCompatibleFluxes = d_ice->d_useCompatibleFluxes;
            varBasket->AMR_subCycleProgressVar = 0;       // for lockstep it's always 0

            advector->inFluxOutFluxVolume(uvel_FC, vvel_FC, wvel_FC, delT, patch, indx,
                bulletProof_test, new_dw, varBasket);

            //__________________________________
            //   advect vol_frac * Porosity
            if (ice_matl) {
                for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                    IntVector c = *iter;
                    //vol_fracICE[c] = Porosity_CC[c] * vol_frac[m][c];
                    vol_fracICE[c] = vol_frac[m][c];
                }
                //advector->advectQ(vol_fracICE, patch, q_advectedICE, varBasket,
                 //   vol_fracX_FC, vol_fracY_FC, vol_fracZ_FC, new_dw);
                
                advector->advectQ(vol_fracICE, mass_L, q_advectedICE, varBasket);
            }

            if (mpm_matl) {         
                //advector->advectQ(vol_fracMPM, patch, q_advectedMPM, varBasket,
                 //   vol_fracX_FC, vol_fracY_FC, vol_fracZ_FC, new_dw);

                //advector->advectQ(VolumeFraction_CC[m], mass_L, q_advectedMPM, varBasket);
                advector->advectQ(vol_frac[m], mass_L, q_advectedMPM, varBasket);
            }

            delete varBasket;

            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
               
                double inv_sumKappa = 1.0 / sumKappa[c];
                termICE[c] -= q_advectedICE[c] * inv_sumKappa; 
                termMPM[c] -= q_advectedMPM[c] * vol_frac[m][c];
            }
            */

            /* // This is not work
            //__________________________________
            // Advection preprocessing
            // - divide vol_frac_cc/vol
            bool bulletProof_test = true;
            advectVarBasket* varBasket = scinew advectVarBasket();

            advector->inFluxOutFluxVolume(uvel_FC, vvel_FC, wvel_FC, delT, patch, indx,
                bulletProof_test, new_dw, varBasket);
            //__________________________________
            //   advect vol_frac
            varBasket->doRefluxing = false;  // don't need to reflux here
            advector->advectQ(vol_frac, patch, q_advected, varBasket,
                vol_fracX_FC, vol_fracY_FC, vol_fracZ_FC, new_dw);

            delete varBasket;

            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                //term2[c] -= q_advected[c];

                double inv_sumKappa = 1.0 / sumKappa[c];
                if (ice_matl) {
                    termICE[c] -= q_advected[c] * inv_sumKappa;
                }
                if (mpm_matl) {
                    termMPM[c] -= q_advected[c] * vol_frac[m][c];
                }
            }*/
        }
        

        delete advector;
        //double very_small_mass = d_TINY_RHO * vol;

        //__________________________________ 
        for (unsigned int m = 0; m < numICEMatls; m++) {
            //Material* matl = m_materialManager->getMaterial(m);
            // indx = matl->getDWIndex();

            ICEMaterial* ice_matl = (ICEMaterial*)m_materialManager->getMaterial("ICE", m);
            int indx = ice_matl->getDWIndex();
            
            constCCVariable<double> mass_L;
            Ghost::GhostType  gac = Ghost::AroundCells;
            new_dw->get(mass_L, Ilb->mass_L_CCLabel, indx, patch, gac, 2);

            constCCVariable<double> kappa;
            new_dw->get(kappa, Ilb->compressibilityLabel, indx, patch, gac, 2);

            CCVariable<double> sp_vol_L, sp_vol_src;
            new_dw->allocateAndPut(sp_vol_L, Ilb->sp_vol_L_CCLabel, indx, patch);
            new_dw->allocateAndPut(sp_vol_src, Ilb->sp_vol_src_CCLabel, indx, patch);
            sp_vol_src.initialize(0.);
            double tiny_rho = 1.e-15;
            tiny_rho = ice_matl->getTinyRho();

            new_dw->get(sp_vol_CC, Ilb->sp_vol_CCLabel, indx, patch, gn, 0);
            new_dw->get(rho_CC, Ilb->rho_CCLabel, indx, patch, gn, 0);
            new_dw->get(f_theta, Ilb->f_theta_CCLabel, indx, patch, gn, 0);
                    
            //__________________________________
            //  compute sp_vol_L * mass
            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                 sp_vol_L[c] = (rho_CC[c] * vol) * sp_vol_CC[c]; // mass * sp_vol_CC
                 //sp_vol_L[c] = (Porosity_CC[c] * rho_CC[c] * vol) * sp_vol_CC[c]; // mass * sp_vol_CC
            }

            //__________________________________
            //   Contributions from models
            constCCVariable<double> Modelsp_vol_src;
            if (d_ice->d_models.size() > 0) {
                new_dw->get(Modelsp_vol_src, Ilb->modelVol_srcLabel, indx, patch, gn, 0);
                for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                    IntVector c = *iter;
                    sp_vol_L[c] += Modelsp_vol_src[c];
                }
            }

            //__________________________________
            //  add the sources to sp_vol_L
            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                //__________________________________
                //  term1
                //double term1 = vol * (termICE[c] + termMPM[c]);
                //double term1 = vol * (termMPM[c]);
                //double term1 = 0;
                //double inv_sumKappa = 1.0 / sumKappa[c];
                //double termICEtotal = -termICE[c] * vol_frac[m][c] * kappa[c];
                // double termMPMtotal = -termMPM[c];
                //double term1 = vol * (termICEtotal + termMPMtotal);

                // term 1 
                double term1 = -vol_frac[indx][c] * kappa[c] * vol * delP[c];

                //  term2
                double term2 = delT * vol *
                    (vol_frac[indx][c] * alpha[indx][c] * Tdot[indx][c] -
                        f_theta[c] * sum_therm_exp[c]);

                // This is actually mass * sp_vol              
                //double src = 0;

                //if (mass_L[c] > very_small_mass)
                //{
                double src = term1 + term2;

               // cerr << "vol_frac[m][c] of m " << indx << "at cell " << c << " is " << vol_frac[indx][c] << endl;
                //cerr << "kappa 1 of m " << indx << "at cell " << c << " is " << kappa[c] << endl;
                //cerr << "term 1 of m " << indx << "at cell " << c << " is " << term1 << endl;
                cerr << "delP " << delP[c] << endl;
                //}

                //cerr << "material " << m << endl;
                //cerr << "mass_L[c] " << mass_L[c] << endl;
                //cerr << "very_small_mass " << very_small_mass << endl;
                //cerr << "termICEtotal " << termICEtotal << endl;
                //cerr << "termMPM[c] " << termMPM[c] << endl;
                //cerr << "term1 " << term1 << endl;
                //cerr << "ICE before sp_vol_L of m " << m << " " << sp_vol_L[c] << endl;

                sp_vol_L[c] += src;
                sp_vol_src[c] = src / (rho_CC[c] * vol);
                //sp_vol_src[c] = src / (Porosity_CC[c] * rho_CC[c] * vol);

                //cerr << "ICE src at cell "  << c << " " << src << endl;
                //cerr << "ICE sp_vol_src[c] at cell " << c << " " << sp_vol_src[c] << endl;
                //cerr << " "  << endl;
            }

            if (d_ice->d_clampSpecificVolume) {
                for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                    IntVector c = *iter;
                    /*`==========TESTING==========*/
                    sp_vol_L[c] = max(sp_vol_L[c], tiny_rho * vol * sp_vol_CC[c]);
                    //sp_vol_L[c] = max(sp_vol_L[c], Porosity_CC[c] * tiny_rho * vol * sp_vol_CC[c]);
                    /*==========TESTING==========`*/
                }
            }

            // Debug
            for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                //cerr << "ICE after sp_vol_L of m " << m << " " << sp_vol_L[c] << endl;
                //cerr << "ICE sp_vol_src" << sp_vol_src[c] << endl;
            }

            //__________________________________
            // Apply boundary conditions
            setSpecificVolBC(sp_vol_L, "SpecificVol", true, rho_CC, vol_frac[m],
                patch, m_materialManager, indx);


            //____ B U L L E T   P R O O F I N G----
            // ignore BP if recompute time step has already been requested
            IntVector neg_cell;
            bool rts = new_dw->recomputeTimeStep();

            if (!areAllValuesPositive(sp_vol_L, neg_cell) && !rts) {
                cout << "\nMPMICE2:WARNING......Negative specific Volume" << endl;
                cout << "cell              " << neg_cell << " level " << level->getIndex() << endl;
                cout << "matl              " << indx << endl;
                cout << "sum_thermal_exp   " << sum_therm_exp[neg_cell] << endl;
                cout << "sp_vol_src        " << sp_vol_src[neg_cell] << endl;
                cout << "mass sp_vol_L     " << sp_vol_L[neg_cell] << endl;
                cout << "mass sp_vol_L_old "
                    << (rho_CC[neg_cell] * vol * sp_vol_CC[neg_cell]) << endl;
                cout << "-----------------------------------" << endl;
                //        ostringstream warn;
                //        int L = level->getIndex();
                //        warn<<"ERROR ICE:("<<L<<"):computeLagrangianSpecificVolumeRF, mat "<<indx
                //            << " cell " <<neg_cell << " sp_vol_L is negative\n";
                //        throw InvalidValue(warn.str(), __FILE__, __LINE__);

                new_dw->put(bool_or_vartype(true), VarLabel::find(abortTimeStep_name));
                new_dw->put(bool_or_vartype(true), VarLabel::find(recomputeTimeStep_name));
            }
        }  // end numALLMatl loop
    }  // patch loop
}

//______________________________________________________________________
// Exact the same with MPMICE
//
void MPMICE2::scheduleComputeCCVelAndTempRates(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* mpm_matls)
{
    printSchedule(patches, cout_doing, "MPMICE2::scheduleComputeCCVelAndTempRates");

    Task* t = scinew Task("MPMICE2::computeCCVelAndTempRates",
        this, &MPMICE2::computeCCVelAndTempRates);

    Ghost::GhostType  gn = Ghost::None;

    t->requires(Task::OldDW, Ilb->timeStepLabel);
    t->requires(Task::OldDW, Ilb->delTLabel, getLevel(patches));
    t->requires(Task::NewDW, Ilb->mass_L_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->mom_L_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->int_eng_L_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->mom_L_ME_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->eng_L_ME_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->int_eng_source_CCLabel, gn);
    t->requires(Task::NewDW, Ilb->mom_source_CCLabel, gn);
    t->requires(Task::OldDW, Mlb->heatRate_CCLabel, gn);

    t->computes(Ilb->dTdt_CCLabel);
    t->computes(Ilb->dVdt_CCLabel);
    t->computes(Mlb->heatRate_CCLabel);

    sched->addTask(t, patches, mpm_matls);
}

//______________________________________________________________________
//
void MPMICE2::computeCCVelAndTempRates(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, VarLabel::find(timeStep_name));

    bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing computeCCVelAndTempRates");

        //__________________________________
        // This is where I interpolate the CC 
        // changes to NCs for the MPMMatls
        unsigned int numMPMMatls = m_materialManager->getNumMatls("MPM");

        delt_vartype delT;
        old_dw->get(delT, Ilb->delTLabel);

        for (unsigned int m = 0; m < numMPMMatls; m++) {
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();
            CCVariable<double> dTdt_CC, heatRate;
            CCVariable<Vector> dVdt_CC;

            constCCVariable<double> mass_L_CC, old_heatRate;
            constCCVariable<Vector> mom_L_ME_CC, old_mom_L_CC, mom_source;
            constCCVariable<double> eng_L_ME_CC, old_int_eng_L_CC, int_eng_src;

            double cv = mpm_matl->getSpecificHeat();

            Ghost::GhostType  gn = Ghost::None;
            new_dw->get(old_mom_L_CC, Ilb->mom_L_CCLabel, indx, patch, gn, 0);
            new_dw->get(old_int_eng_L_CC, Ilb->int_eng_L_CCLabel, indx, patch, gn, 0);
            new_dw->get(mass_L_CC, Ilb->mass_L_CCLabel, indx, patch, gn, 0);
            new_dw->get(mom_L_ME_CC, Ilb->mom_L_ME_CCLabel, indx, patch, gn, 0);
            new_dw->get(eng_L_ME_CC, Ilb->eng_L_ME_CCLabel, indx, patch, gn, 0);
            old_dw->get(old_heatRate, Mlb->heatRate_CCLabel, indx, patch, gn, 0);
            new_dw->get(mom_source, Ilb->mom_source_CCLabel, indx, patch, gn, 0);
            new_dw->get(int_eng_src, Ilb->int_eng_source_CCLabel, indx, patch, gn, 0);

            new_dw->allocateAndPut(dTdt_CC, Ilb->dTdt_CCLabel, indx, patch);
            new_dw->allocateAndPut(dVdt_CC, Ilb->dVdt_CCLabel, indx, patch);
            new_dw->allocateAndPut(heatRate, Mlb->heatRate_CCLabel, indx, patch);

            dTdt_CC.initialize(0.0);
            dVdt_CC.initialize(Vector(0.0));
            //__________________________________
            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                if (!d_rigidMPM) {
                    dVdt_CC[c] = (mom_L_ME_CC[c] - (old_mom_L_CC[c] - mom_source[c]))
                        / (mass_L_CC[c] * delT);
                }
                dTdt_CC[c] = (eng_L_ME_CC[c] - (old_int_eng_L_CC[c] - int_eng_src[c]))
                    / (mass_L_CC[c] * cv * delT);
                double heatRte = (eng_L_ME_CC[c] - old_int_eng_L_CC[c]) / delT;
                heatRate[c] = .05 * heatRte + .95 * old_heatRate[c];

               // cerr << "dTdt_CC " << dTdt_CC[c] << endl;
               //cerr << "dVdt_CC " << dVdt_CC[c] << endl;
                //cerr << "heatRate " << heatRate[c] << endl;
            }

            for (CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                if (!d_rigidMPM) {
                    dVdt_CC[c] = (mom_L_ME_CC[c] - (old_mom_L_CC[c] - mom_source[c]))
                        / (mass_L_CC[c] * delT);
                }
                dTdt_CC[c] = (eng_L_ME_CC[c] - (old_int_eng_L_CC[c] - int_eng_src[c]))
                    / (mass_L_CC[c] * cv * delT);
                double heatRte = (eng_L_ME_CC[c] - old_int_eng_L_CC[c]) / delT;
                heatRate[c] = .05 * heatRte + .95 * old_heatRate[c];

                // cerr << "dTdt_CC " << dTdt_CC[c] << endl;
                //cerr << "dVdt_CC " << dVdt_CC[c] << endl;
                //cerr << "heatRate " << heatRate[c] << endl;
            }
            setBC(dTdt_CC, "set_if_sym_BC", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            setBC(dVdt_CC, "set_if_sym_BC", patch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
        }
    }  //patches
}

//______________________________________________________________________
//
void MPMICE2::scheduleInterpolateCCToNC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* mpm_matls)
{
    if (!d_mpm->flags->doMPMOnLevel(getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels()))
        return;

    printSchedule(patches, cout_doing, "MPMICE2::scheduleInterpolateCCToNC");

    Task* t = scinew Task("MPMICE2::interpolateCCToNC",
        this, &MPMICE2::interpolateCCToNC);
    const MaterialSubset* mss = mpm_matls->getUnion();
    Ghost::GhostType  gan = Ghost::AroundNodes;
    Ghost::GhostType  gac = Ghost::AroundCells;

    t->requires(Task::OldDW, Ilb->delTLabel, getLevel(patches));
    t->requires(Task::NewDW, Ilb->dVdt_CCLabel, gan, 1);
    t->requires(Task::NewDW, Ilb->dTdt_CCLabel, gan, 1);

    if (d_ice->d_models.size() > 0) {
        t->requires(Task::NewDW, MIlb->cMassLabel, gac, 1);
        t->requires(Task::NewDW, Ilb->modelMass_srcLabel, gac, 1);
    }

    t->modifies(Mlb->gVelocityStarLabel, mss);
    t->modifies(Mlb->gAccelerationLabel, mss);
    t->computes(Mlb->massBurnFractionLabel, mss);
    t->computes(Mlb->dTdt_NCLabel);

    sched->addTask(t, patches, mpm_matls);
}

//______________________________________________________________________
//
void MPMICE2::interpolateCCToNC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing interpolateCCToNC");

        //__________________________________
        // This is where I interpolate the CC 
        // changes to NCs for the MPMMatls
        unsigned int numMPMMatls = m_materialManager->getNumMatls("MPM");

        delt_vartype delT;
        old_dw->get(delT, Ilb->delTLabel);

        for (unsigned int m = 0; m < numMPMMatls; m++) {
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int indx = mpm_matl->getDWIndex();
            NCVariable<Vector> gacceleration, gvelocity;
            NCVariable<double> dTdt_NC, massBurnFraction;

            constCCVariable<double> dTdt_CC;
            constCCVariable<Vector> dVdt_CC;

            new_dw->getModifiable(gvelocity, Mlb->gVelocityStarLabel, indx, patch);
            new_dw->getModifiable(gacceleration, Mlb->gAccelerationLabel, indx, patch);

            Ghost::GhostType  gan = Ghost::AroundNodes;
            new_dw->get(dTdt_CC, Ilb->dTdt_CCLabel, indx, patch, gan, 1);
            new_dw->get(dVdt_CC, Ilb->dVdt_CCLabel, indx, patch, gan, 1);

            new_dw->allocateAndPut(massBurnFraction,
                Mlb->massBurnFractionLabel, indx, patch);
            new_dw->allocateAndPut(dTdt_NC, Mlb->dTdt_NCLabel, indx, patch);

            dTdt_NC.initialize(0.0);
            massBurnFraction.initialize(0.);
            IntVector cIdx[8];
            //__________________________________
            //  
            for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
                patch->findCellsFromNode(*iter, cIdx);
                for (int in = 0; in < 8; in++) {
                    gvelocity[*iter] += dVdt_CC[cIdx[in]] * delT * .125;
                    gacceleration[*iter] += dVdt_CC[cIdx[in]] * .125;
                    dTdt_NC[*iter] += dTdt_CC[cIdx[in]] * .125;
                }               
            }

            for (NodeIterator iter = patch->getNodeIterator();
                !iter.done(); iter++) {
                IntVector c = *iter;
                //cerr << "gvelocity " << gvelocity[c] << endl;
                //cerr << "gacceleration of material " << indx  << " at node " << c <<  " is " << gacceleration[c] << endl;
                //cerr << "dTdt_NC " << dTdt_NC[c] << endl;
            }

            //__________________________________
            //  inter-material phase transformation
            if (d_ice->d_models.size() > 0) {
                constCCVariable<double> modelMass_src, mass_CC;
                Ghost::GhostType  gac = Ghost::AroundCells;
                new_dw->get(modelMass_src, Ilb->modelMass_srcLabel, indx, patch, gac, 1);
                new_dw->get(mass_CC, MIlb->cMassLabel, indx, patch, gac, 1);

                for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
                    patch->findCellsFromNode(*iter, cIdx);
                    for (int in = 0; in < 8; in++) {
                        massBurnFraction[*iter] +=
                            (fabs(modelMass_src[cIdx[in]]) / mass_CC[cIdx[in]]) * .125;

                    }
                }
            }  // if(models >0 )
        }  // mpmMatls
    }  //patches
}

// Optional function //
//______________________________________________________________________
//
void MPMICE2::scheduleRefinePressCC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSubset* press_matl,
    const MaterialSet* matls)
{
    printSchedule(patches, cout_doing, "MPMICE2::scheduleRefinePressCC");

    MaterialSet* press_matls = scinew MaterialSet();
    press_matls->add(0);
    press_matls->addReference();

    scheduleRefineVariableCC(sched, patches, press_matls, Ilb->press_CCLabel);
    if (press_matls->removeReference())
        delete press_matls;
}

//if (do_mlmpmice2)//

//______________________________________________________________________
//
void MPMICE2::scheduleCoarsenLagrangianValuesMPM(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* mpm_matls)
{
    printSchedule(patches, cout_doing, "MPMICE2:scheduleCoarsenLagrangianValues mpm_matls");

    scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->rho_CCLabel,
        1e-12, true, "std"); // modifies
    scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->mass_L_CCLabel,
        1.9531e-15, false, "sum");
    scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->mom_L_CCLabel,
        Vector(0, 0, 0), false, "sum");
    scheduleCoarsenVariableCC(sched, patches, mpm_matls, Ilb->int_eng_L_CCLabel,
        0.0, false, "sum");
}

//______________________________________________________________________
//
void MPMICE2::scheduleRefineCC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* mpm_matls)
{
    if (!d_mpm->flags->doMPMOnLevel(getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels()))
        return;

    printSchedule(patches, cout_doing, "MPMICE2::scheduleRefineCC");
    scheduleRefineVariableCC(sched, patches, mpm_matls, Ilb->dTdt_CCLabel);
    scheduleRefineVariableCC(sched, patches, mpm_matls, Ilb->dVdt_CCLabel);
}

void MPMICE2::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
    if (d_switchCriteria) {
        d_switchCriteria->scheduleSwitchTest(level, sched);
    }
}

//______________________________________________________________________
void MPMICE2::scheduleRefineInterface(const LevelP& fineLevel,
    SchedulerP& scheduler,
    bool needOld, bool needNew)
{
    d_ice->scheduleRefineInterface(fineLevel, scheduler, needOld, needNew);
    d_mpm->scheduleRefineInterface(fineLevel, scheduler, needOld, needNew);

    if (fineLevel->getIndex() > 0 && scheduler->isCopyDataTimestep() &&
        d_mpm->flags->doMPMOnLevel(fineLevel->getIndex(),
            fineLevel->getGrid()->numLevels())) {
        cout_doing << d_myworld->myRank()
            << " MPMICE2::scheduleRefineInterface \t\t\tL-"
            << fineLevel->getIndex() << endl;

        Task* task = scinew Task("MPMICE2::refineCoarseFineInterface",
            this, &MPMICE2::refineCoarseFineInterface);

        const MaterialSet* all_matls = m_materialManager->allMaterials();
        const MaterialSubset* one_matl = d_ice->d_press_matl;

        task->modifies(Mlb->NC_CCweightLabel, one_matl);

        scheduler->addTask(task, fineLevel->eachPatch(), all_matls);
    }
}


void MPMICE2::refineCoarseFineInterface(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset*,
    DataWarehouse* fine_old_dw,
    DataWarehouse* fine_new_dw)
{
    // This isn't actually refining anything, it is simply reinitializing
    // NC_CCweight after regridding on all levels finer than 0 because
    // copyData doesn't copy extra cell data.
    const Level* level = getLevel(patches);
    if (level->getIndex() > 0) {
        cout_doing << d_myworld->myRank()
            << " Doing refineCoarseFineInterface" << "\t\t\t MPMICE2 L-"
            << level->getIndex() << " Patches: " << *patches << endl;

        for (int p = 0; p < patches->size(); p++) {
            const Patch* patch = patches->get(p);
            //__________________________________
            //NC_CCweight
            NCVariable<double> NC_CCweight;
            fine_new_dw->getModifiable(NC_CCweight, Mlb->NC_CCweightLabel, 0, patch);
            //__________________________________
            // - Initialize NC_CCweight = 0.125
            // - Find the walls with symmetry BC and double NC_CCweight
            NC_CCweight.initialize(0.125);
            vector<Patch::FaceType>::const_iterator iter;
            vector<Patch::FaceType> bf;
            patch->getBoundaryFaces(bf);

            for (iter = bf.begin(); iter != bf.end(); ++iter) {
                Patch::FaceType face = *iter;
                int mat_id = 0;
                if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {

                    for (CellIterator iter = patch->getFaceIterator(face, Patch::FaceNodes);
                        !iter.done(); iter++) {
                        NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
                    } // cell iterator
                } // if symmetry
            } // for patch faces
        } // for patches
    } // if level
}
//______________________________________________________________________
void MPMICE2::scheduleRefine(const PatchSet* patches,
    SchedulerP& sched)
{
    d_ice->scheduleRefine(patches, sched);
    d_mpm->scheduleRefine(patches, sched);

    printSchedule(patches, cout_doing, "MPMICE2::scheduleRefine");

    Task* task = scinew Task("MPMICE2::refine", this, &MPMICE2::refine);

    task->requires(Task::OldDW, Ilb->timeStepLabel);

    task->computes(Mlb->heatRate_CCLabel);
    task->computes(Ilb->sp_vol_CCLabel);
    task->computes(MIlb->vel_CCLabel);
    task->computes(Ilb->temp_CCLabel);

    sched->addTask(task, patches, m_materialManager->allMaterials("MPM"));
}

//______________________________________________________________________    
void MPMICE2::scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched)
{
    d_ice->scheduleCoarsen(coarseLevel, sched);
    d_mpm->scheduleCoarsen(coarseLevel, sched);
}

//______________________________________________________________________
void MPMICE2::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
    SchedulerP& sched)
{
    d_ice->scheduleInitialErrorEstimate(coarseLevel, sched);
    d_mpm->scheduleInitialErrorEstimate(coarseLevel, sched);
}

//______________________________________________________________________
void MPMICE2::scheduleErrorEstimate(const LevelP& coarseLevel,
    SchedulerP& sched)
{
    d_ice->scheduleErrorEstimate(coarseLevel, sched);
    d_mpm->scheduleErrorEstimate(coarseLevel, sched);
}

//______________________________________________________________________
void MPMICE2::scheduleRefineVariableCC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* matls,
    const VarLabel* variable)
{
    ostringstream taskName;
    taskName << "MPMICE2::refineVariable(" << variable->getName() << ")";
    Task* t;

    // the sgis don't like accepting a templated function over a function call for some reason...
    void (MPMICE2:: * func)(const ProcessorGroup*, const PatchSubset*, const MaterialSubset*,
        DataWarehouse*, DataWarehouse*, const VarLabel*);

    switch (variable->typeDescription()->getSubType()->getType()) {
    case TypeDescription::double_type:
        func = &MPMICE2::refineVariableCC<double>;
        t = scinew Task(taskName.str().c_str(), this, func, variable);
        break;
    case TypeDescription::Vector:
        func = &MPMICE2::refineVariableCC<Vector>;
        t = scinew Task(taskName.str().c_str(), this, func, variable);
        break;
    default:
        throw InternalError("Unknown variable type for refine", __FILE__, __LINE__);
    }

    Ghost::GhostType  gac = Ghost::AroundCells;
    t->requires(Task::NewDW, variable, 0, Task::CoarseLevel, 0, Task::NormalDomain, gac, 1);
    t->computes(variable);
    sched->addTask(t, patches, matls);
}

//______________________________________________________________________
template<typename T>
void MPMICE2::scheduleCoarsenVariableCC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* matls,
    const VarLabel* variable,
    T defaultValue,
    bool modifies,
    const string& coarsenMethod)
{
    // The SGI compiler does't like accepting a templated function over
    // a function call for some reason...  We use this hack to force it
    // to figure out the correct type of the function.
    void (MPMICE2:: * func)(const ProcessorGroup*, const PatchSubset*,
        const MaterialSubset*, DataWarehouse*, DataWarehouse*,
        const VarLabel*, T, bool, string);
    func = &MPMICE2::coarsenVariableCC<T>;
    ostringstream taskName;

    taskName << "MPMICE2::coarsenVariableCC(" << variable->getName()
        << (modifies ? " modified" : "") << ")";

    Task* t = scinew Task(taskName.str().c_str(), this, func,
        variable, defaultValue, modifies, coarsenMethod);

    Ghost::GhostType  gn = Ghost::None;
    Task::MaterialDomainSpec ND = Task::NormalDomain;

    t->requires(Task::OldDW, Ilb->timeStepLabel);

    t->requires(Task::NewDW, variable, 0, Task::FineLevel, 0, ND, gn, 0);

    if (coarsenMethod == "massWeighted") {
        t->requires(Task::NewDW, MIlb->cMassLabel, 0, Task::FineLevel, 0, ND, gn, 0);
    }

    if (modifies) {
        t->modifies(variable);
    }
    else {
        t->computes(variable);
    }
    sched->addTask(t, patches, matls);
}


//______________________________________________________________________
template<typename T>
void MPMICE2::scheduleCoarsenVariableNC(SchedulerP& sched,
    const PatchSet* patches,
    const MaterialSet* matls,
    const VarLabel* variable,
    T defaultValue,
    bool modifies,
    string coarsenMethod)
{
    // The SGI compiler does't like accepting a templated function over
    // a function call for some reason...  We use this hack to force it
    // to figure out the correct type of the function.
    void (MPMICE2:: * func)(const ProcessorGroup*, const PatchSubset*,
        const MaterialSubset*, DataWarehouse*, DataWarehouse*,
        const VarLabel*, T, bool, string);
    func = &MPMICE2::coarsenVariableNC<T>;
    ostringstream taskName;

    taskName << "MPMICE2::coarsenVariableNC(" << variable->getName()
        << (modifies ? " modified" : "") << ")";

    Task* t = scinew Task(taskName.str().c_str(), this, func,
        variable, defaultValue, modifies, coarsenMethod);

    //Ghost::GhostType  gn = Ghost::None;
    Ghost::GhostType  gan = Ghost::AroundNodes;
    Task::MaterialDomainSpec ND = Task::NormalDomain;

    const LevelP fineLevel = getLevel(patches)->getFinerLevel();
    IntVector refineRatio(fineLevel->getRefinementRatio());
    int ghost = max(refineRatio.x(), refineRatio.y());
    ghost = max(ghost, refineRatio.z());

    t->requires(Task::NewDW, variable, 0, Task::FineLevel, 0, ND, gan, ghost);

    if (modifies) {
        t->modifies(variable);
    }
    else {
        t->computes(variable);
    }
    sched->addTask(t, patches, matls);
}

//______________________________________________________________________
void
MPMICE2::refine(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset* /*matls*/,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw)
{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, VarLabel::find(timeStep_name));

    bool isNotInitialTimeStep = (timeStep > 0);

    for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        printTask(patches, patch, cout_doing, "Doing refine");

        unsigned int numMPMMatls = m_materialManager->getNumMatls("MPM");

        for (unsigned int m = 0; m < numMPMMatls; m++) {
            MPMMaterial* mpm_matl = (MPMMaterial*)m_materialManager->getMaterial("MPM", m);
            int dwi = mpm_matl->getDWIndex();

            cout_doing << d_myworld->myRank() << " Doing refine on patch "
                << patch->getID() << " material # = " << dwi << endl;

            // for now, create 0 heat flux
            CCVariable<double> heatFlux;
            new_dw->allocateAndPut(heatFlux, Mlb->heatRate_CCLabel, dwi, patch);
            heatFlux.initialize(0.0);

            CCVariable<double> rho_micro, sp_vol_CC, rho_CC, Temp_CC, vol_frac_CC;
            CCVariable<Vector> vel_CC;

            new_dw->allocateTemporary(rho_micro, patch);
            new_dw->allocateTemporary(rho_CC, patch);
            new_dw->allocateTemporary(vol_frac_CC, patch);

            new_dw->allocateAndPut(sp_vol_CC, Ilb->sp_vol_CCLabel, dwi, patch);
            new_dw->allocateAndPut(Temp_CC, MIlb->temp_CCLabel, dwi, patch);
            new_dw->allocateAndPut(vel_CC, MIlb->vel_CCLabel, dwi, patch);

            mpm_matl->initializeDummyCCVariables(rho_micro, rho_CC,
                Temp_CC, vel_CC,
                vol_frac_CC, patch);
            //__________________________________
            //  Set boundary conditions                                     
            setBC(rho_micro, "Density", patch, m_materialManager, dwi, new_dw, isNotInitialTimeStep);
            setBC(Temp_CC, "Temperature", patch, m_materialManager, dwi, new_dw, isNotInitialTimeStep);
            setBC(vel_CC, "Velocity", patch, m_materialManager, dwi, new_dw, isNotInitialTimeStep);

            for (CellIterator iter = patch->getExtraCellIterator();
                !iter.done(); iter++) {
                sp_vol_CC[*iter] = 1.0 / rho_micro[*iter];
            }

            //__________________________________
            //    B U L L E T   P R O O F I N G
            IntVector neg_cell;
            ostringstream warn;
            if (!areAllValuesPositive(rho_CC, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);

                warn << "ERROR MPMICE2::actuallyInitialize, mat " << dwi << " cell: "
                    << neg_cell << ", position: " << pt << ", rho_CC is negative\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesPositive(Temp_CC, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);

                warn << "ERROR MPMICE2::actuallyInitialize, mat " << dwi << " cell "
                    << neg_cell << ", position: " << pt << ", Temp_CC is negative\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
            if (!areAllValuesPositive(sp_vol_CC, neg_cell)) {
                Point pt = patch->getCellPosition(neg_cell);

                warn << "ERROR MPMICE2::actuallyInitialize, mat " << dwi << " cell "
                    << neg_cell << ", position: " << pt << ", sp_vol_CC is negative\n";
                throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
            }
        }  //mpmMatls
    }  //patches
}

//______________________________________________________________________
//
template<typename T>
void MPMICE2::refineVariableCC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset* matls,
    DataWarehouse*,
    DataWarehouse* new_dw,
    const VarLabel* variable)
{
    const Level* fineLevel = getLevel(patches);
    const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();
    IntVector refineRatio(fineLevel->getRefinementRatio());

    for (int p = 0; p < patches->size(); p++) {
        const Patch* finePatch = patches->get(p);
        ostringstream message;
        message << "Doing refineVariableCC (" << variable->getName() << ")\t\t\t";
        printTask(patches, finePatch, cout_doing, message.str());

        // region of fine space that will correspond to the coarse we need to get
        IntVector cl, ch, fl, fh;
        IntVector bl(0, 0, 0);  // boundary layer cells
        int nGhostCells = 1;
        bool returnExclusiveRange = true;

        getCoarseLevelRange(finePatch, coarseLevel, cl, ch, fl, fh, bl,
            nGhostCells, returnExclusiveRange);

        for (int m = 0; m < matls->size(); m++) {
            int indx = matls->get(m);

            CCVariable<T> fine_q_CC;
            new_dw->allocateAndPut(fine_q_CC, variable, indx, finePatch);

            constCCVariable<T> coarse_q_CC;

            new_dw->getRegion(coarse_q_CC, variable, indx, coarseLevel, cl, ch, false);

            // Only interpolate over the intersection of the fine and coarse patches
            // coarse cell 
      //      linearInterpolation<T>(coarse_q_CC, coarseLevel, fineLevel,
      //                             refineRatio, lo, hi, fine_q_CC);

            piecewiseConstantInterpolation<T>(coarse_q_CC, fineLevel,
                fl, fh, fine_q_CC);
        }
    }
}



//__________________________________
//
template<typename T>
void MPMICE2::coarsenDriver_stdNC(IntVector cl,
    IntVector ch,
    IntVector refinementRatio,
    double ratio,
    const Level* coarseLevel,
    constNCVariable<T>& fine_q_NC,
    NCVariable<T>& coarse_q_NC)
{
    T zero(0.0);
    // iterate over coarse level cells
    const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
    Vector DX = coarseLevel->dCell();
    IntVector range(refinementRatio.x() / 2, refinementRatio.y() / 2, refinementRatio.z() / 2);

    IntVector varLow = fine_q_NC.getLowIndex();
    IntVector varHigh = fine_q_NC.getHighIndex();

    for (NodeIterator iter(cl, ch); !iter.done(); iter++) {
        IntVector c = *iter;
        IntVector fineNode = coarseLevel->mapNodeToFiner(c);
        Point coarseLoc = coarseLevel->getNodePosition(c);

        IntVector start = Max(fineNode - range, varLow);
        IntVector end = Min(fineNode + range, varHigh);

        // for each coarse level cell iterate over the fine level cells
        T q_NC_tmp(zero);

        for (NodeIterator inner(start, end); !inner.done(); inner++) {
            IntVector fc(*inner);
            Point fineLoc = fineLevel->getNodePosition(fc);
            Vector C2F = fineLoc - coarseLoc;
            Vector Vweight = C2F / DX;
            double weight = (1. - fabs(Vweight.x())) *
                (1. - fabs(Vweight.y())) *
                (1. - fabs(Vweight.z()));
            q_NC_tmp += fine_q_NC[fc] * weight;
        }
        coarse_q_NC[c] = q_NC_tmp;
    }
}

//______________________________________________________________________
//
template<typename T>
void MPMICE2::coarsenVariableCC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset* matls,
    DataWarehouse* old_dw,
    DataWarehouse* new_dw,
    const VarLabel* variable,
    T defaultValue,
    bool modifies,
    string coarsenMethod)
{
    timeStep_vartype timeStep;
    old_dw->get(timeStep, Ilb->timeStepLabel);

    bool isNotInitialTimeStep = (timeStep > 0);

    const Level* coarseLevel = getLevel(patches);
    const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();

    IntVector refineRatio(fineLevel->getRefinementRatio());
    double ratio = 1. / (refineRatio.x() * refineRatio.y() * refineRatio.z());

    for (int p = 0; p < patches->size(); p++) {
        const Patch* coarsePatch = patches->get(p);
        ostringstream message;
        message << "Doing CoarsenVariableCC (" << variable->getName() << ")\t\t\t";
        printTask(patches, coarsePatch, cout_doing, message.str());

        for (int m = 0; m < matls->size(); m++) {
            int indx = matls->get(m);

            CCVariable<T> coarse_q_CC;
            if (modifies) {
                new_dw->getModifiable(coarse_q_CC, variable, indx, coarsePatch);
            }
            else {
                new_dw->allocateAndPut(coarse_q_CC, variable, indx, coarsePatch);
            }
            coarse_q_CC.initialize(defaultValue);

            Level::selectType finePatches;
            coarsePatch->getFineLevelPatches(finePatches);
            for (unsigned int i = 0; i < finePatches.size(); i++) {
                const Patch* finePatch = finePatches[i];

                IntVector cl, ch, fl, fh;
                getFineLevelRange(coarsePatch, finePatch, cl, ch, fl, fh);
                if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
                    continue;
                }

                constCCVariable<T> fine_q_CC;
                new_dw->getRegion(fine_q_CC, variable, indx, fineLevel, fl, fh, false);

                //__________________________________
                //  call the coarsening function
                ASSERT((coarsenMethod == "std" || coarsenMethod == "sum"
                    || coarsenMethod == "massWeighted"));
                if (coarsenMethod == "std") {
                    coarsenDriver_std(cl, ch, fl, fh, refineRatio, ratio, coarseLevel,
                        fine_q_CC, coarse_q_CC);
                }
                if (coarsenMethod == "sum") {
                    ratio = 1.0;
                    coarsenDriver_std(cl, ch, fl, fh, refineRatio, ratio, coarseLevel,
                        fine_q_CC, coarse_q_CC);
                }
                if (coarsenMethod == "massWeighted") {
                    constCCVariable<double> cMass;
                    new_dw->getRegion(cMass, MIlb->cMassLabel, indx, fineLevel, fl, fh, false);

                    coarsenDriver_massWeighted(cl, ch, fl, fh, refineRatio, coarseLevel,
                        cMass, fine_q_CC, coarse_q_CC);
                }
            }  // fine patches
            // Set BCs on coarsened data.  This sucks--Steve
            if (variable->getName() == "temp_CC") {
                setBC(coarse_q_CC, "Temperature", coarsePatch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            }
            else if (variable->getName() == "rho_CC") {
                setBC(coarse_q_CC, "Density", coarsePatch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            }
            else if (variable->getName() == "vel_CC") {
                setBC(coarse_q_CC, "Velocity", coarsePatch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            }
            else if (variable->getName() == "c.mass" ||
                variable->getName() == "sp_vol_CC" ||
                variable->getName() == "mom_L_CC" ||
                variable->getName() == "int_eng_L_CC") {
                setBC(coarse_q_CC, "set_if_sym_BC", coarsePatch, m_materialManager, indx, new_dw, isNotInitialTimeStep);
            }
        }  // matls
    }  // coarse level
}

//______________________________________________________________________
//
template<typename T>
void MPMICE2::coarsenVariableNC(const ProcessorGroup*,
    const PatchSubset* patches,
    const MaterialSubset* matls,
    DataWarehouse*,
    DataWarehouse* new_dw,
    const VarLabel* variable,
    T defaultValue,
    bool modifies,
    string coarsenMethod)
{
    const Level* coarseLevel = getLevel(patches);
    const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();

    IntVector refineRatio(fineLevel->getRefinementRatio());
    double ratio = 1. / (refineRatio.x() * refineRatio.y() * refineRatio.z());

    for (int p = 0; p < patches->size(); p++) {
        const Patch* coarsePatch = patches->get(p);
        ostringstream message;
        message << "Doing CoarsenVariableNC (" << variable->getName() << ")\t\t\t";
        printTask(patches, coarsePatch, cout_doing, message.str());

        for (int m = 0; m < matls->size(); m++) {
            int indx = matls->get(m);

            NCVariable<T> coarse_q_NC;
            if (modifies) {
                new_dw->getModifiable(coarse_q_NC, variable, indx, coarsePatch);
            }
            else {
                new_dw->allocateAndPut(coarse_q_NC, variable, indx, coarsePatch);
            }
            coarse_q_NC.initialize(defaultValue);

            Level::selectType finePatches;
            coarsePatch->getFineLevelPatches(finePatches);
            for (unsigned int i = 0; i < finePatches.size(); i++) {
                const Patch* finePatch = finePatches[i];

                IntVector cl, ch, fl, fh;

                IntVector padding(refineRatio.x() / 2, refineRatio.y() / 2, refineRatio.z() / 2);
                getFineLevelRangeNodes(coarsePatch, finePatch, cl, ch, fl, fh, padding);


                if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
                    continue;
                }

                constNCVariable<T> fine_q_NC;
                new_dw->getRegion(fine_q_NC, variable, indx, fineLevel, fl, fh, false);

                //__________________________________
                //  call the coarsening function
                ASSERT(coarsenMethod == "sum");
                if (coarsenMethod == "sum") {
                    coarsenDriver_stdNC(cl, ch, refineRatio, ratio, coarseLevel,
                        fine_q_NC, coarse_q_NC);
                }
            }  // fine patches
        }  // matls
    }  // coarse level
}

/* _____________________________________________________________________
MPMICE2::scheduleFinalizeTimestep--
This task called at the very bottom of the timestep,
after scheduleTimeAdvance and the scheduleCoarsen.
This is scheduled on every level.
_____________________________________________________________________*/
void
MPMICE2::scheduleFinalizeTimestep(const LevelP& level, SchedulerP& sched)
{
    cout_doing << "----------------------------" << endl;
    cout_doing << d_myworld->myRank() << " MPMICE2::scheduleFinalizeTimestep\t\t\t\tL-" << level->getIndex() << endl;

    const PatchSet* ice_patches = level->eachPatch();
    const MaterialSet* ice_matls = m_materialManager->allMaterials("ICE");
    const MaterialSet* all_matls = m_materialManager->allMaterials();
    const MaterialSet* mpm_matls = m_materialManager->allMaterials("MPM");
    const MaterialSubset* ice_matls_sub = ice_matls->getUnion();
    const MaterialSubset* press_matl = d_ice->d_press_matl;

    d_ice->scheduleConservedtoPrimitive_Vars(sched, ice_patches, ice_matls_sub, press_matl,
        ice_matls,
        "finalizeTimestep");

    d_ice->scheduleTestConservation(sched, ice_patches, ice_matls_sub,
        all_matls);

    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->scheduleDoAnalysis_preReloc(sched, level);
        }
    }

    // only do on finest level until we get AMR MPM
    if (level->getIndex() == level->getGrid()->numLevels() - 1)
        sched->scheduleParticleRelocation(level,
            Mlb->pXLabel_preReloc,
            d_mpm->d_particleState_preReloc,
            Mlb->pXLabel,
            d_mpm->d_particleState,
            Mlb->pParticleIDLabel, mpm_matls);

    //__________________________________
    //  on the fly analysis
    if (d_analysisModules.size() != 0) {
        vector<AnalysisModule*>::iterator iter;
        for (iter = d_analysisModules.begin();
            iter != d_analysisModules.end(); iter++) {
            AnalysisModule* am = *iter;
            am->scheduleDoAnalysis(sched, level);
        }
    }

    cout_doing << "---------------------------------------------------------" << endl;
}