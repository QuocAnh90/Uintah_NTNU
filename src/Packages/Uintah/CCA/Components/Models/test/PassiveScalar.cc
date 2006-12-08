
#include <Packages/Uintah/CCA/Components/ICE/ICEMaterial.h>
#include <Packages/Uintah/CCA/Components/ICE/ConservationTest.h>
#include <Packages/Uintah/CCA/Components/ICE/BoundaryCond.h>
#include <Packages/Uintah/CCA/Components/ICE/Diffusion.h>
#include <Packages/Uintah/CCA/Components/Models/test/PassiveScalar.h>
#include <Packages/Uintah/CCA/Components/Regridder/PerPatchVars.h>
#include <Packages/Uintah/CCA/Ports/Scheduler.h>
#include <Packages/Uintah/Core/Exceptions/ParameterNotFound.h>
#include <Packages/Uintah/Core/Exceptions/ProblemSetupException.h>
#include <Packages/Uintah/Core/Exceptions/InvalidValue.h>
#include <Packages/Uintah/Core/GeometryPiece/GeometryPieceFactory.h>
#include <Packages/Uintah/Core/GeometryPiece/UnionGeometryPiece.h>
#include <Packages/Uintah/Core/Grid/Box.h>
#include <Packages/Uintah/Core/Grid/Variables/CellIterator.h>
#include <Packages/Uintah/Core/Grid/Variables/PerPatch.h>
#include <Packages/Uintah/Core/Grid/SimulationState.h>
#include <Packages/Uintah/Core/Grid/Variables/VarTypes.h>
#include <Packages/Uintah/Core/Exceptions/ParameterNotFound.h>
#include <Packages/Uintah/Core/Parallel/ProcessorGroup.h>
#include <Packages/Uintah/Core/Exceptions/InvalidValue.h>

#include <Core/Containers/StaticArray.h>
#include <Core/Math/MiscMath.h>
#include <iostream>
#include <Core/Util/DebugStream.h>
#include <stdio.h>

using namespace Uintah;
using namespace std;
//__________________________________
//  To turn on the output
//  setenv SCI_DEBUG "MODELS_DOING_COUT:+,PASSIVE_SCALAR_DBG_COUT:+"
//  PASSIVE_SCALAR_DBG:  dumps out during problemSetup 
static DebugStream cout_doing("MODELS_DOING_COUT", false);
static DebugStream cout_dbg("PASSIVE_SCALAR_DBG_COUT", false);
//______________________________________________________________________              
PassiveScalar::PassiveScalar(const ProcessorGroup* myworld, 
                            ProblemSpecP& params,
                            const bool doAMR)
  : ModelInterface(myworld), params(params)
{
  d_doAMR = doAMR;
  d_matl_set = 0;
  lb  = scinew ICELabel();
  Slb = scinew PassiveScalarLabel();
}

//__________________________________
PassiveScalar::~PassiveScalar()
{
  if(d_matl_set && d_matl_set->removeReference()) {
    delete d_matl_set;
  }
  
  VarLabel::destroy(d_scalar->scalar_CCLabel);
  VarLabel::destroy(d_scalar->scalar_source_CCLabel);
  VarLabel::destroy(d_scalar->diffusionCoefLabel);
  VarLabel::destroy(d_scalar->mag_grad_scalarLabel);
  VarLabel::destroy(Slb->lastProbeDumpTimeLabel);
  VarLabel::destroy(Slb->sum_scalar_fLabel);

  delete lb;
  delete Slb;
  
  for(vector<Region*>::iterator iter = d_scalar->regions.begin();
                                iter != d_scalar->regions.end(); iter++){
    Region* region = *iter;
    delete region;
  }
}

//__________________________________
PassiveScalar::Region::Region(GeometryPieceP piece, ProblemSpecP& ps)
  : piece(piece)
{
  ps->require("scalar", initialScalar);
  ps->getWithDefault("sinusoidalInitialize", sinusoidalInitialize, false);
  if(sinusoidalInitialize){
    ps->getWithDefault("freq",freq,IntVector(0,0,0));
  }
  
  ps->getWithDefault("linearInitialize", linearInitialize, false);
  if(linearInitialize){
    ps->getWithDefault("slope",slope,Vector(0,0,0));
  }
  ps->getWithDefault("cubicInitialize",     cubicInitialize,     false);
  ps->getWithDefault("quadraticInitialize", quadraticInitialize, false);
  if(quadraticInitialize){
    ps->getWithDefault("coeff",coeff,Vector(0,0,0));
  }
  
  uniformInitialize = true;
  if(sinusoidalInitialize || linearInitialize || quadraticInitialize || cubicInitialize){
    uniformInitialize = false;
  }
}
//______________________________________________________________________
//     P R O B L E M   S E T U P
void PassiveScalar::problemSetup(GridP&, SimulationStateP& in_state,
                        ModelSetup* setup)
{
  cout_doing << "Doing problemSetup \t\t\t\tPASSIVE_SCALAR" << endl;
  d_sharedState = in_state;
  d_matl = d_sharedState->parseAndLookupMaterial(params, "material");

  vector<int> m(1);
  m[0] = d_matl->getDWIndex();
  d_matl_set = new MaterialSet();
  d_matl_set->addAll(m);
  d_matl_set->addReference();
  d_matl_sub = d_matl_set->getUnion();
  
  //__________________________________
  // - create Label names
  // - register the scalar to be transported
  d_scalar = new Scalar();
  d_scalar->index = 0;
  d_scalar->name  = "f";
  
  const TypeDescription* td_CCdouble = CCVariable<double>::getTypeDescription();
  //const TypeDescription* td_CCVector = CCVariable<Vector>::getTypeDescription();
    
  d_scalar->scalar_CCLabel =     VarLabel::create("scalar-f",       td_CCdouble);
  d_scalar->diffusionCoefLabel = VarLabel::create("scalar-diffCoef",td_CCdouble);
  d_scalar->scalar_source_CCLabel = 
                                 VarLabel::create("scalar-f_src",   td_CCdouble);
  d_scalar->mag_grad_scalarLabel = 
                               VarLabel::create("mag_grad_scalar-f",td_CCdouble);                                 
  
                               
                                 
                                 
  Slb->lastProbeDumpTimeLabel =  VarLabel::create("lastProbeDumpTime", 
                                            max_vartype::getTypeDescription());
  Slb->sum_scalar_fLabel      =  VarLabel::create("sum_scalar_f", 
                                            sum_vartype::getTypeDescription());
  
  d_modelComputesThermoTransportProps = true;
  
  setup->registerTransportedVariable(d_matl_set->getSubset(0),
                                     d_scalar->scalar_CCLabel,
                                     d_scalar->scalar_source_CCLabel);  

  //__________________________________
  //  register the AMRrefluxing variables                               
  if(d_doAMR){
    setup->registerAMR_RefluxVariable(d_matl_set->getSubset(0),
                                      d_scalar->scalar_CCLabel);
  }
  //__________________________________
  // Read in the constants for the scalar
   ProblemSpecP child = params->findBlock("scalar");
   if (!child){
     throw ProblemSetupException("PassiveScalar: Couldn't find scalar tag", __FILE__, __LINE__);    
   }

   child->getWithDefault("test_conservation", d_test_conservation, false);

   ProblemSpecP const_ps = child->findBlock("constants");
   if(!const_ps) {
     throw ProblemSetupException("PassiveScalar: Couldn't find constants tag", __FILE__, __LINE__);
   }
       
   const_ps->getWithDefault("initialize_diffusion_knob",       
                            d_scalar->initialize_diffusion_knob,   0);
                            
   const_ps->getWithDefault("diffusivity",  d_scalar->diff_coeff, 0.0);
   
   const_ps->getWithDefault("AMR_Refinement_Criteria", d_scalar->refineCriteria,1e100);

  //__________________________________
  //  Read in the geometry objects for the scalar
  for (ProblemSpecP geom_obj_ps = child->findBlock("geom_object");
    geom_obj_ps != 0;
    geom_obj_ps = geom_obj_ps->findNextBlock("geom_object") ) {
    vector<GeometryPieceP> pieces;
    GeometryPieceFactory::create(geom_obj_ps, pieces);

    GeometryPieceP mainpiece;
    if(pieces.size() == 0){
     throw ParameterNotFound("No piece specified in geom_object", __FILE__, __LINE__);
    } else if(pieces.size() > 1){
     mainpiece = scinew UnionGeometryPiece(pieces);
    } else {
     mainpiece = pieces[0];
    }

    d_scalar->regions.push_back(scinew Region(mainpiece, geom_obj_ps));
  }
  if(d_scalar->regions.size() == 0) {
    throw ProblemSetupException("Variable: scalar-f does not have any initial value regions",
                                __FILE__, __LINE__);
  }

  //__________________________________
  //  Read in probe locations for the scalar field
  d_usingProbePts = false;
  ProblemSpecP probe_ps = child->findBlock("probePoints");
  if (probe_ps) {
    probe_ps->require("probeSamplingFreq", d_probeFreq);
     
    Vector location = Vector(0,0,0);
    map<string,string> attr;                    
    for (ProblemSpecP prob_spec = probe_ps->findBlock("location"); prob_spec != 0; 
                      prob_spec = prob_spec->findNextBlock("location")) {
                      
      prob_spec->get(location);
      prob_spec->getAttributes(attr);
      string name = attr["name"];
      
      d_probePts.push_back(location);
      d_probePtsNames.push_back(name);
      d_usingProbePts = true;
    }
  } 
}
//__________________________________
//
void PassiveScalar::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP model_ps = ps->appendChild("Model");
  model_ps->setAttribute("type","PassiveScalar");

  model_ps->appendElement("material",d_matl->getName());
  ProblemSpecP scalar_ps = model_ps->appendChild("scalar");
  scalar_ps->setAttribute("name","f");
  scalar_ps->appendElement("test_conservation",d_test_conservation);


  ProblemSpecP const_ps = scalar_ps->appendChild("constants");
  const_ps->appendElement("initialize_diffusion_knob",
                          d_scalar->initialize_diffusion_knob);
  const_ps->appendElement("diffusivity",d_scalar->diff_coeff);
  const_ps->appendElement("AMR_Refinement_Criteria",d_scalar->refineCriteria);

  vector<Region*>::const_iterator iter;
  for ( iter = d_scalar->regions.begin(); iter != d_scalar->regions.end(); iter++) {
    ProblemSpecP geom_ps = scalar_ps->appendChild("geom_object");
    (*iter)->piece->outputProblemSpec(geom_ps);
    geom_ps->appendElement("scalar",(*iter)->initialScalar);
  }

  if (d_usingProbePts) {
    throw ProblemSetupException("Can't output tag",__FILE__,__LINE__);
    ProblemSpecP probe_ps = scalar_ps->appendChild("probePoints");
    probe_ps->appendElement("probeSamplingFreq",d_probeFreq);
    
    for (unsigned int i = 0; i < d_probePts.size(); i++) {
      probe_ps->appendElement("location",d_probePts[i]);
      probe_ps->setAttribute("name",d_probePtsNames[i]);
    }
  }

}


//______________________________________________________________________
//      S C H E D U L E   I N I T I A L I Z E
void PassiveScalar::scheduleInitialize(SchedulerP& sched,
                                   const LevelP& level,
                                   const ModelInfo*)
{
  cout_doing << "PassiveScalar::scheduleInitialize " << endl;
  Task* t = scinew Task("PassiveScalar::initialize", 
                  this, &PassiveScalar::initialize);
  
  t->computes(d_scalar->scalar_CCLabel);
  t->computes(Slb->lastProbeDumpTimeLabel);
  
  sched->addTask(t, level->eachPatch(), d_matl_set);
}
//______________________________________________________________________
//       I N I T I A L I Z E
void PassiveScalar::initialize(const ProcessorGroup*, 
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  cout_doing << "Doing Initialize \t\t\t\t\tPASSIVE_SCALAR" << endl;
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    int indx = d_matl->getDWIndex();
    
    CCVariable<double>  f;
    new_dw->allocateAndPut(f, d_scalar->scalar_CCLabel, indx, patch);

    f.initialize(0);
    
    //__________________________________
    //  Uniform initialization scalar field in a region
    for(vector<Region*>::iterator iter = d_scalar->regions.begin();
                                  iter != d_scalar->regions.end(); iter++){
      Region* region = *iter;
      
      if(region->uniformInitialize){
      
        for(CellIterator iter = patch->getExtraCellIterator();
            !iter.done(); iter++){
          IntVector c = *iter;
          Point p = patch->cellPosition(c);            
          if(region->piece->inside(p)) {
            f[c] = region->initialScalar;
          }
        } // Over cells
      }
      
      if(region->linearInitialize){
      
        for(CellIterator iter = patch->getExtraCellIterator();
            !iter.done(); iter++){
          IntVector c = *iter;
          Point p = patch->cellPosition(c);            
          if(region->piece->inside(p)) {
            f[c] = region->initialScalar;
          }
        } // Over cells
      }

      
      //__________________________________
      // Sinusoidal & linear initialization
      if(!region->uniformInitialize){
        IntVector freq = region->freq;
        // bulletproofing
        if(region->sinusoidalInitialize && freq.x()==0 && freq.y()==0 && freq.z()==0){
          throw ProblemSetupException("PassiveScalar: you need to specify a <freq> whenever you use sinusoidalInitialize", __FILE__, __LINE__);
        }
        
        Vector slope = region->slope;
        if(region->linearInitialize && slope.x()==0 && slope.y()==0 && slope.z()==0){
          throw ProblemSetupException("PassiveScalar: you need to specify a <slope> whenever you use linearInitialize", __FILE__, __LINE__);
        }
        
        Vector coeff = region->coeff;
        cout << " coeff " << coeff << endl;
        if(region->quadraticInitialize && coeff.x()==0 && coeff.y()==0 && coeff.z()==0){
          throw ProblemSetupException("PassiveScalar: you need to specify a <coeff> whenever you use quadraticInitialize", __FILE__, __LINE__);
        }

        Point lo = region->piece->getBoundingBox().lower();
        Point hi = region->piece->getBoundingBox().upper();
        Vector dist = hi.asVector() - lo.asVector();
        
        for(CellIterator iter = patch->getExtraCellIterator();
            !iter.done(); iter++){
          IntVector c = *iter;
          Point p = patch->cellPosition(c);            
          if(region->piece->inside(p)) {
            // normalized distance
            Vector d = (p.asVector() - lo.asVector() )/dist;
            
            if(region->sinusoidalInitialize){
              f[c] = sin( 2.0 * freq.x() * d.x() * M_PI) + sin( 2.0 * freq.y() * d.y() * M_PI)  + sin( 2.0 * freq.z() * d.z() * M_PI);
            }
            if(region->linearInitialize){
              f[c] = slope.x() * d.x() + slope.y() * d.y() + slope.z() * d.z(); 
            }
            if(region->cubicInitialize){    
              if(d.x() <= 0.5)
                f[c] = -1.3333333*pow(d.x(),3)  + pow(d.x(),2);
              else{
                f[c] = -1.3333333*pow( (1.0 - d.x()),3) + pow( (1.0 - d.x()),2);
              } 
            }            
            if(region->quadraticInitialize){
              f[c] = coeff.x() * d.x() * d.x() 
                   + coeff.y() * d.y() * d.y() 
                   + coeff.z() * d.z() * d.z();
            }
          }
        }
      }  // sinusoidal Initialize  
    } // regions
    setBC(f,"scalar-f", patch, d_sharedState,indx, new_dw);
     
    //__________________________________
    //  Dump out a header for the probe point files
    new_dw->put(max_vartype(0.0), Slb->lastProbeDumpTimeLabel);
    if (d_usingProbePts){
      FILE *fp;
      IntVector cell;
      string udaDir = d_dataArchiver->getOutputLocation();
      
        for (unsigned int i =0 ; i < d_probePts.size(); i++) {
          if(patch->findCell(Point(d_probePts[i]),cell) ) {
            string filename=udaDir + "/" + d_probePtsNames[i].c_str() + ".dat";
            fp = fopen(filename.c_str(), "a");
            fprintf(fp, "%% Time Scalar Field at [%e, %e, %e], at cell [%i, %i, %i]\n", 
                    d_probePts[i].x(),d_probePts[i].y(), d_probePts[i].z(),
                    cell.x(), cell.y(), cell.z() );
            fclose(fp);
        }
      }  // loop over probes
    }  // if using probe points
  }  // patches
}

//______________________________________________________________________     
void PassiveScalar::scheduleModifyThermoTransportProperties(SchedulerP& sched,
                                                   const LevelP& level,
                                                   const MaterialSet* /*ice_matls*/)
{
  cout_doing << "PASSIVE_SCALAR::scheduleModifyThermoTransportProperties" << endl;
  Task* t = scinew Task("PassiveScalar::modifyThermoTransportProperties", 
                   this,&PassiveScalar::modifyThermoTransportProperties);
  t->computes(d_scalar->diffusionCoefLabel);
  sched->addTask(t, level->eachPatch(), d_matl_set);
}
//______________________________________________________________________
void PassiveScalar::modifyThermoTransportProperties(const ProcessorGroup*, 
                                                const PatchSubset* patches,
                                                const MaterialSubset*,
                                                DataWarehouse* /*old_dw*/,
                                                DataWarehouse* new_dw)
{ 
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    cout_doing << "Doing modifyThermoTransportProperties on patch "
               <<patch->getID()<< "\t PassiveScalar" << endl;
   
    int indx = d_matl->getDWIndex();
    CCVariable<double> diffusionCoeff;
    new_dw->allocateAndPut(diffusionCoeff, 
                           d_scalar->diffusionCoefLabel,indx, patch);  
        
    diffusionCoeff.initialize(d_scalar->diff_coeff);
  }
} 

//______________________________________________________________________
void PassiveScalar::computeSpecificHeat(CCVariable<double>& ,
                                    const Patch* ,
                                    DataWarehouse* ,
                                    const int )
{ 
  //none
} 

//______________________________________________________________________
void PassiveScalar::scheduleComputeModelSources(SchedulerP& sched,
                                            const LevelP& level,
                                            const ModelInfo* mi)
{
  cout_doing << "PASSIVE_SCALAR::scheduleComputeModelSources " << endl;
  Task* t = scinew Task("PassiveScalar::computeModelSources", 
                   this,&PassiveScalar::computeModelSources, mi);
                     
  Ghost::GhostType  gac = Ghost::AroundCells;
  t->requires(Task::NewDW, d_scalar->diffusionCoefLabel, gac,1);
  t->requires(Task::OldDW, d_scalar->scalar_CCLabel,     gac,1); 
  t->modifies(d_scalar->scalar_source_CCLabel);

  //  if dumping out probePts
  if (d_usingProbePts){
    t->requires(Task::OldDW, Slb->lastProbeDumpTimeLabel);
    t->computes(Slb->lastProbeDumpTimeLabel);
  }

  sched->addTask(t, level->eachPatch(), d_matl_set);
}

//______________________________________________________________________
void PassiveScalar::computeModelSources(const ProcessorGroup*, 
                                    const PatchSubset* patches,
                                    const MaterialSubset* /*matls*/,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw,
                                    const ModelInfo* mi)
{
  const Level* level = getLevel(patches);
  delt_vartype delT;
  old_dw->get(delT, mi->delT_Label, level);     
  Ghost::GhostType gac = Ghost::AroundCells;
  
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    cout_doing << "Doing computeModelSources... on patch "<<patch->getID()
               << "\t\t\t PassiveScalar" << endl;
   
    constCCVariable<double> f_old, diff_coeff;
    CCVariable<double>  f_src;
    
    int indx = d_matl->getDWIndex();
    old_dw->get(f_old,      d_scalar->scalar_CCLabel,      indx, patch, gac,1);
    new_dw->get(diff_coeff, d_scalar->diffusionCoefLabel,  indx, patch, gac,1);
    new_dw->allocateAndPut(f_src, d_scalar->scalar_source_CCLabel,
                                                           indx, patch);
    f_src.initialize(0.0);
                                                            
    //__________________________________
    //  scalar diffusion
    double diff_coeff_test = d_scalar->diff_coeff;
    if(diff_coeff_test != 0.0){ 
      bool use_vol_frac = false; // don't include vol_frac in diffusion calc.
      SFCXVariable<double> placeHolderX;
      SFCYVariable<double> placeHolderY;
      SFCZVariable<double> placeHolderZ;

      scalarDiffusionOperator(new_dw, patch, use_vol_frac, f_old,
                              placeHolderX, placeHolderY, placeHolderZ,
                              f_src, diff_coeff, delT);
    }

    //__________________________________
    //  dump out the probe points
    if (d_usingProbePts ) {
      
      max_vartype lastDumpTime;
      old_dw->get(lastDumpTime, Slb->lastProbeDumpTimeLabel);
      double oldProbeDumpTime = lastDumpTime;
      
      double time = d_dataArchiver->getCurrentTime();
      double nextDumpTime = oldProbeDumpTime + 1.0/d_probeFreq;
      
      if (time >= nextDumpTime){        // is it time to dump the points
        FILE *fp;
        string udaDir = d_dataArchiver->getOutputLocation();
        IntVector cell_indx;
        
        // loop through all the points and dump if that patch contains them
        for (unsigned int i =0 ; i < d_probePts.size(); i++) {
          if(patch->findCell(Point(d_probePts[i]),cell_indx) ) {
            string filename=udaDir + "/" + d_probePtsNames[i].c_str() + ".dat";
            fp = fopen(filename.c_str(), "a");
            fprintf(fp, "%16.15E  %16.15E\n",time, f_old[cell_indx]);
            fclose(fp);
          }
        }
        oldProbeDumpTime = time;
      }  // time to dump
      new_dw->put(max_vartype(oldProbeDumpTime), Slb->lastProbeDumpTimeLabel);
    } // if(probePts)  
  }
}
//__________________________________      
void PassiveScalar::scheduleComputeStableTimestep(SchedulerP&,
                                      const LevelP&,
                                      const ModelInfo*)
{
  // None necessary...
}

//______________________________________________________________________
void PassiveScalar::scheduleTestConservation(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const ModelInfo* mi)
{
  const Level* level = getLevel(patches);
  int L = level->getIndex();
  
  if(d_test_conservation && L == 0){
    cout_doing << "PASSIVESCALAR::scheduleTestConservation " << endl;
    Task* t = scinew Task("PassiveScalar::testConservation", 
                     this,&PassiveScalar::testConservation, mi);

    Ghost::GhostType  gn = Ghost::None;
    // compute sum(scalar_f * mass)
    t->requires(Task::NewDW, d_scalar->scalar_CCLabel, gn,0); 
    t->requires(Task::NewDW, mi->density_CCLabel,      gn,0);
    t->requires(Task::NewDW, lb->uvel_FCMELabel,       gn,0); 
    t->requires(Task::NewDW, lb->vvel_FCMELabel,       gn,0); 
    t->requires(Task::NewDW, lb->wvel_FCMELabel,       gn,0); 
    t->computes(Slb->sum_scalar_fLabel);

    sched->addTask(t, patches, d_matl_set);
  }
}

//______________________________________________________________________
void PassiveScalar::testConservation(const ProcessorGroup*, 
                                     const PatchSubset* patches,
                                     const MaterialSubset* /*matls*/,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw,
                                     const ModelInfo* mi)
{
  const Level* level = getLevel(patches);
  delt_vartype delT;
  old_dw->get(delT, mi->delT_Label, level);     
  Ghost::GhostType gn = Ghost::None; 
  
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    cout_doing << "Doing testConservation on patch "<<patch->getID()
               << "\t\t\t PassiveScalar" << endl;
               
    //__________________________________
    //  conservation of f test
    constCCVariable<double> rho_CC, f;
    constSFCXVariable<double> uvel_FC;
    constSFCYVariable<double> vvel_FC;
    constSFCZVariable<double> wvel_FC;
    int indx = d_matl->getDWIndex();
    new_dw->get(f,       d_scalar->scalar_CCLabel,indx,patch,gn,0);
    new_dw->get(rho_CC,  mi->density_CCLabel,     indx,patch,gn,0); 
    new_dw->get(uvel_FC, lb->uvel_FCMELabel,      indx,patch,gn,0); 
    new_dw->get(vvel_FC, lb->vvel_FCMELabel,      indx,patch,gn,0); 
    new_dw->get(wvel_FC, lb->wvel_FCMELabel,      indx,patch,gn,0); 
    Vector dx = patch->dCell();
    double cellVol = dx.x()*dx.y()*dx.z();

    CCVariable<double> q_CC;
    new_dw->allocateTemporary(q_CC, patch);

    for(CellIterator iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      q_CC[c] = rho_CC[c]*cellVol*f[c];
    }

    double sum_mass_f;
    conservationTest<double>(patch, delT, q_CC, uvel_FC, vvel_FC, wvel_FC, 
                             sum_mass_f);

    new_dw->put(sum_vartype(sum_mass_f), Slb->sum_scalar_fLabel);
  }
}

//______________________________________________________________________
//
void PassiveScalar::scheduleErrorEstimate(const LevelP& coarseLevel,
                                          SchedulerP& sched)
{
  cout_doing << "PassiveScalar::scheduleErrorEstimate \t\t\tL-" 
             << coarseLevel->getIndex() << '\n';
  
  Task* t = scinew Task("PassiveScalar::errorEstimate", 
                  this, &PassiveScalar::errorEstimate, false);  
  
  Ghost::GhostType  gac  = Ghost::AroundCells; 
  
  t->requires(Task::NewDW, d_scalar->scalar_CCLabel,  d_matl_sub, gac,1);
  
  t->computes(d_scalar->mag_grad_scalarLabel, d_matl_sub);
  t->modifies(d_sharedState->get_refineFlag_label(),      d_sharedState->refineFlagMaterials());
  t->modifies(d_sharedState->get_refinePatchFlag_label(), d_sharedState->refineFlagMaterials());
 
  // define the material set of 0 and whatever the passive scalar index is 
  // don't add matl 0 twice 
  MaterialSet* matl_set;
  vector<int> m;
  m.push_back(0);
  if(d_matl->getDWIndex() != 0){
    m.push_back(d_matl->getDWIndex());
  }
  matl_set = new MaterialSet();
  matl_set->addAll(m);
  matl_set->addReference();
    
  sched->addTask(t, coarseLevel->eachPatch(), matl_set); 
}
/*_____________________________________________________________________
 Function~  PassiveScalar::errorEstimate--
______________________________________________________________________*/
void PassiveScalar::errorEstimate(const ProcessorGroup*,
			             const PatchSubset* patches,
			             const MaterialSubset*,
			             DataWarehouse*,
			             DataWarehouse* new_dw,
                                  bool)
{
  cout_doing << "Doing errorEstimate \t\t\t\t\t PassiveScalar"<< endl;
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    
    Ghost::GhostType  gac  = Ghost::AroundCells;
    const VarLabel* refineFlagLabel = d_sharedState->get_refineFlag_label();
    const VarLabel* refinePatchLabel= d_sharedState->get_refinePatchFlag_label();
    
    CCVariable<int> refineFlag;
    new_dw->getModifiable(refineFlag, refineFlagLabel, 0, patch);      

    PerPatch<PatchFlagP> refinePatchFlag;
    new_dw->get(refinePatchFlag, refinePatchLabel, 0, patch);

    int indx = d_matl->getDWIndex();
    constCCVariable<double> f;
    CCVariable<double> mag_grad_f;
    
    new_dw->get(f,                     d_scalar->scalar_CCLabel,  indx ,patch,gac,1);
    new_dw->allocateAndPut(mag_grad_f, d_scalar->mag_grad_scalarLabel, 
                           indx,patch);
    mag_grad_f.initialize(0.0);
    
    //__________________________________
    // compute gradient
    Vector dx = patch->dCell(); 
    
    for(CellIterator iter = patch->getCellIterator();!iter.done();iter++){
      IntVector c = *iter;
      Vector grad_f;
      for(int dir = 0; dir <3; dir ++ ) {
        IntVector r = c;
        IntVector l = c; 
        double inv_dx = 0.5 /dx[dir];
        r[dir] += 1;
        l[dir] -= 1;
        grad_f[dir] = (f[r] - f[l])*inv_dx;
      }
      mag_grad_f[c] = grad_f.length();
    }
    //__________________________________
    // set refinement flag
    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();
    for(CellIterator iter = patch->getCellIterator();!iter.done();iter++){
      IntVector c = *iter;
      if( mag_grad_f[c] > d_scalar->refineCriteria){
        refineFlag[c] = true;
        refinePatch->set();
      }
    }
  }  // patches
}
