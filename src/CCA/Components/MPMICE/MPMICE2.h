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

#ifndef UINTAH_HOMEBREW_MPMICE2_H
#define UINTAH_HOMEBREW_MPMICE2_H

 // NOTE: SOMETHING IS FUBAR IN THE DEFINITION OF fflux AS IT RELATES TO
 // swapbytes. AS SUCH, Advector.h MUST BE CALLED BEFORE ApplicationCommon.h

#include <CCA/Components/ICE/Advection/Advector.h>
#include <CCA/Components/Application/ApplicationCommon.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/Models/MultiMatlExchange/ExchangeModel.h>
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <CCA/Components/MPM/Materials/Contact/Contact.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/RigidMPM.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>
#include <CCA/Ports/SwitchingCriteria.h>
#include <Core/Geometry/Vector.h>

namespace Uintah {
    class ICE;
    class ICELabel;
    class MPMLabel;
    class MPMICELabel;
    class Output;

    /**************************************

    CLASS
       MPMICE2

       Short description...

    GENERAL INFORMATION

       MPMICE2.h

       Steven G. Parker
       Department of Computer Science
       University of Utah

       Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


    KEYWORDS
       MPMICE2

    DESCRIPTION
       Long description...

    WARNING

    ****************************************/

    enum MPMType2 { STAND_MPMICE2 = 0, RIGID_MPMICE2, SHELL_MPMICE2, FRACTURE_MPMICE2 };

    class MPMICE2 : public ApplicationCommon {

    public:
        MPMICE2(const ProcessorGroup* myworld,
            const MaterialManagerP materialManager,
            MPMType2 type, const bool doAMR = false);

        virtual ~MPMICE2();

        virtual double recomputeDelT(const double delT);

        virtual void problemSetup(const ProblemSpecP& params,
            const ProblemSpecP& restart_prob_spec,
            GridP& grid);

        virtual void outputProblemSpec(ProblemSpecP& ps);

        virtual void scheduleInitialize(const LevelP& level,
            SchedulerP&);

        virtual void scheduleRestartInitialize(const LevelP& level,
            SchedulerP& sched);

        virtual void restartInitialize();

        virtual void scheduleComputeStableTimeStep(const LevelP& level,
            SchedulerP&);

        // scheduleTimeAdvance version called by the AMR simulation controller.
        virtual void scheduleTimeAdvance(const LevelP& level,
            SchedulerP&);

        virtual void scheduleFinalizeTimestep(const LevelP& level,
            SchedulerP&);

        void scheduleInterpolateNCToCC_0(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void interpolateNCToCC_0(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void scheduleCoarsenCC_0(SchedulerP&,
            const PatchSet*,
            const MaterialSet*);

        void scheduleCoarsenNCMass(SchedulerP&,
            const PatchSet*,
            const MaterialSet*);       
        
        void scheduleComputePressure(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void computeEquilibrationPressure(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse*,
            DataWarehouse*);

        void computeEquilPressure_1_matl(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);


        void scheduleComputeVelICE_FC(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void computeVelICE_FC(const ProcessorGroup*,
            const PatchSubset*,
            const MaterialSubset*,
            DataWarehouse*,
            DataWarehouse*);
            
            /*
        void scheduleComputeVel_FC(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void computeVel_FC(const ProcessorGroup*,
            const PatchSubset*,
            const MaterialSubset*,
            DataWarehouse*,
            DataWarehouse*);
            */

        template<class T> void computeVelICEFace(int dir, CellIterator it,
            IntVector adj_offset,
            double dx,
            double delT,
            double gravity,
            constCCVariable<double>& rho_CC,
            constCCVariable<double>& sp_vol_CC,
            constCCVariable<Vector>& vel_CC,
            constCCVariable<double>& press_CC,
            T& vel_FC,
            T& gradP_FC);
        
        void scheduleComputeVelMPM_FC(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void computeVelMPM_FC(const ProcessorGroup*,
            const PatchSubset*,
            const MaterialSubset*,
            DataWarehouse*,
            DataWarehouse*);
        
        template<class T> void computeVelMPMFace(int dir, CellIterator it,
            IntVector adj_offset,
            double dx,
            double delT,
            double gravity,
            double rho_init,
            constCCVariable<double>& rho_CC,
            constCCVariable<double>& stress_CC,
            constCCVariable<Vector>& vel_CC,
            constCCVariable<double>& press_CC,
            T& vel_FC,
            T& grad_P_FC,
            T& grad_stress_FC);

        void scheduleInterpolatePressCCToPressNC(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSet*);

        void interpolatePressCCToPressNC(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void scheduleInterpolatePAndGradP(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void interpolatePAndGradP(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void scheduleComputeLagrangianValuesMPM(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSet*);

        void computeLagrangianValuesMPM(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void scheduleCoarsenLagrangianValuesMPM(SchedulerP&,
            const PatchSet*,
            const MaterialSet*);


        void scheduleComputeLagrangianSpecificVolume(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSubset*,
            const MaterialSet*);

        void computeLagrangianSpecificVolume(const ProcessorGroup*,
            const PatchSubset* patches,
            const MaterialSubset* matls,
            DataWarehouse*,
            DataWarehouse*);

        void scheduleComputeCCVelAndTempRates(SchedulerP&, const PatchSet*,
            const MaterialSet*);

        void computeCCVelAndTempRates(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void scheduleInterpolateCCToNC(SchedulerP&, const PatchSet*,
            const MaterialSet*);

        void interpolateCCToNC(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        //______________________________________________________________________
        //       A C T U A L   S T E P S : 
        void actuallyInitialize(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse*,
            DataWarehouse* new_dw);


        void actuallyInitializeAddedMPMMaterial(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse*,
            DataWarehouse* new_dw);

        void interpolateCCToNCRefined(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void HEChemistry(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

        void printData(int indx,
            const Patch* patch,
            int   include_EC,
            const std::string& message1,
            const std::string& message2,
            const NCVariable<double>& q_NC);

        void printNCVector(int indx,
            const Patch* patch, int include_EC,
            const std::string& message1,
            const std::string& message2,
            int     component,
            const NCVariable<Vector>& q_NC);

        void scheduleRefinePressCC(SchedulerP&,
            const PatchSet*,
            const MaterialSubset*,
            const MaterialSet*);

        void scheduleRefineCC(SchedulerP&, const PatchSet*,
            const MaterialSet*);

        //__________________________________
        //    R A T E   F O R M                   
        void computeRateFormPressure(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse*,
            DataWarehouse*);

        virtual void scheduleSwitchTest(const LevelP& level, SchedulerP& sched);

        //__________________________________
        //   AMR
        virtual void scheduleRefineInterface(const LevelP& fineLevel,
            SchedulerP& scheduler,
            bool needOld, bool needNew);

        virtual void scheduleRefine(const PatchSet* patches,
            SchedulerP& sched);

        virtual void scheduleCoarsen(const LevelP& coarseLevel,
            SchedulerP& sched);

        void refine(const ProcessorGroup*,
            const PatchSubset* patches,
            const MaterialSubset* matls,
            DataWarehouse*,
            DataWarehouse* new_dw);

        virtual void scheduleInitialErrorEstimate(const LevelP& coarseLevel,
            SchedulerP& sched);

        virtual void scheduleErrorEstimate(const LevelP& coarseLevel,
            SchedulerP& sched);

        void scheduleRefineVariableCC(SchedulerP& sched,
            const PatchSet* patches,
            const MaterialSet* matls,
            const VarLabel* variable);

        template<typename T>
        void scheduleCoarsenVariableCC(SchedulerP& sched,
            const PatchSet* patches,
            const MaterialSet* matls,
            const VarLabel* variable,
            T defaultValue,
            bool modifies,
            const std::string& coarsenMethod);

        template<typename T>
        void scheduleCoarsenVariableNC(SchedulerP& sched,
            const PatchSet* patches,
            const MaterialSet* matls,
            const VarLabel* variable,
            T defaultValue,
            bool modifies,
            std::string coarsenMethod);

        template<typename T>
        void refineVariableCC(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw,
            const VarLabel* variable);

        //
        template<typename T>
        void coarsenDriver_stdNC(IntVector cl,
            IntVector ch,
            IntVector refinementRatio,
            double ratio,
            const Level* coarseLevel,
            constNCVariable<T>& fine_q_NC,
            NCVariable<T>& coarse_q_NC);

        template<typename T>
        void coarsenVariableCC(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw,
            const VarLabel* variable,
            T defaultValue,
            bool modifies,
            std::string coarsenMethod);

        template<typename T>
        void coarsenVariableNC(const ProcessorGroup*,
            const PatchSubset* patch,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw,
            const VarLabel* variable,
            T defaultValue,
            bool modifies,
            std::string coarsenMethod);

        void refineCoarseFineInterface(const ProcessorGroup*,
            const PatchSubset* patches,
            const MaterialSubset*,
            DataWarehouse* fine_old_dw,
            DataWarehouse* fine_new_dw);

        Vector getGravity() const {
            return d_gravity;
        }

        double getRefPress() const {
            return d_ref_press;
        }

    private:

        enum bctype {
            NONE = 0,
            FIXED,
            SYMMETRY,
            NEIGHBOR
        };

        //Advector* d_advector;

    protected:
        MPMICE2(const MPMICE2&);
        MPMICE2& operator=(const MPMICE2&);

        MPMLabel* Mlb;
        ICELabel* Ilb;
        MPMICELabel* MIlb;
        ExchangeModel* d_exchModel;

        bool             d_rigidMPM;
        SerialMPM* d_mpm;
        ICE* d_ice;
        int              d_8or27;
        int              NGN;
        bool             d_testForNegTemps_mpm;
        bool             do_mlmpmice2;

        int              pbx_matl_num;
        MaterialSubset* pbx_matl;

        std::vector<AnalysisModule*>  d_analysisModules;

        SwitchingCriteria* d_switchCriteria;

        std::vector<MPMPhysicalBC*> d_physicalBCs;
        double d_SMALL_NUM;
        double d_TINY_RHO;
        Vector d_gravity;
        double d_ref_press;
    };

} // End namespace Uintah

#endif
