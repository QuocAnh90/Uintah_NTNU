/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMDiffusionLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>
using namespace Uintah;
using namespace std;

MPMLabel::MPMLabel()
{

  // Time Step
  timeStepLabel =
    VarLabel::create(timeStep_name, timeStep_vartype::getTypeDescription() );

  // Simulation Time
  simulationTimeLabel =
    VarLabel::create(simTime_name, simTime_vartype::getTypeDescription() );

  // delta t
  VarLabel* nonconstDelt =
    VarLabel::create(delT_name, delt_vartype::getTypeDescription() );
  nonconstDelt->allowMultipleComputes();
  delTLabel = nonconstDelt;


  diffusion = scinew MPMDiffusionLabel();

  // Heat flux from fire

  fAllParticlesUpdated = VarLabel::create("finalUpdateFlag",
                                          PerPatch<bool>::getTypeDescription());

  heatRate_CCLabel = 
    VarLabel::create("heatRate_CC",  CCVariable<double>::getTypeDescription());

  // Particle Variables

  //non PermanentParticleState
  pPressureLabel  = VarLabel::create( "p.pressure",
                        ParticleVariable<double>::getTypeDescription() );
  
  pScratchVecLabel  = VarLabel::create( "p.scratchvec",
                        ParticleVariable<Vector>::getTypeDescription() );
  
  pScratchLabel  = VarLabel::create( "p.scratch",
                        ParticleVariable<double>::getTypeDescription() );
  
  // for visualization only
  pScaleFactorLabel = VarLabel::create( "p.scalefactor",
                        ParticleVariable<Matrix3>::getTypeDescription());
  
  pScaleFactorLabel_preReloc = VarLabel::create( "p.scalefactor+",
                        ParticleVariable<Matrix3>::getTypeDescription());
  
  // for thermal stress
  pTempCurrentLabel = VarLabel::create( "p.tempCurrent",
                        ParticleVariable<double>::getTypeDescription() ); 
    
  pXXLabel  = VarLabel::create( "p.xx",
                        ParticleVariable<Point>::getTypeDescription() );
                     
  p_qLabel  = VarLabel::create( "p.q",
                        ParticleVariable<double>::getTypeDescription() );
  
  p_qLabel_preReloc  = VarLabel::create( "p.q+",
                        ParticleVariable<double>::getTypeDescription() );
  
  pColorLabel  = VarLabel::create( "p.color",
                        ParticleVariable<double>::getTypeDescription() );
                     
  pColorLabel_preReloc  = VarLabel::create( "p.color+",
                        ParticleVariable<double>::getTypeDescription() );
                     
  pPartitionUnityLabel   = VarLabel::create( "p.partitionUnity",
                        ParticleVariable<double>::getTypeDescription() );
  
  //PermanentParticleState
  pDeformationMeasureLabel = VarLabel::create("p.deformationMeasure",
                        ParticleVariable<Matrix3>::getTypeDescription());

  pVelGradLabel = VarLabel::create("p.velGrad",
                        ParticleVariable<Matrix3>::getTypeDescription());

  pStressLabel = VarLabel::create( "p.stress",
                        ParticleVariable<Matrix3>::getTypeDescription() );
  
  pVolumeLabel = VarLabel::create( "p.volume",
                        ParticleVariable<double>::getTypeDescription());
  
  pVolumeDeformedLabel = VarLabel::create( "p.volumedeformed",
                        ParticleVariable<double>::getTypeDescription());
  
  pMassLabel = VarLabel::create( "p.mass",
                        ParticleVariable<double>::getTypeDescription() );
  
  pVelocityLabel = VarLabel::create( "p.velocity", 
                        ParticleVariable<Vector>::getTypeDescription() );
  
  pVelocitySSPlusLabel = VarLabel::create( "p.velocitySSPlus", 
                        ParticleVariable<Vector>::getTypeDescription() );
  
  pExternalForceLabel = VarLabel::create( "p.externalforce",
                        ParticleVariable<Vector>::getTypeDescription() );

  pExternalForceCorner1Label = VarLabel::create( "p.externalforcecorner1",
                        ParticleVariable<Point>::getTypeDescription() );

  pExternalForceCorner2Label = VarLabel::create( "p.externalforcecorner2",
                        ParticleVariable<Point>::getTypeDescription() );

  pExternalForceCorner3Label = VarLabel::create( "p.externalforcecorner3",
                        ParticleVariable<Point>::getTypeDescription() );

  pExternalForceCorner4Label = VarLabel::create( "p.externalforcecorner4",
                        ParticleVariable<Point>::getTypeDescription() );
  
  pXLabel = VarLabel::create("p.x",
                             ParticleVariable<Point>::getTypeDescription(),
                             IntVector(0,0,0), VarLabel::PositionVariable);
  
  pTemperatureLabel = VarLabel::create( "p.temperature",
                        ParticleVariable<double>::getTypeDescription() );
  
  // for thermal stress
  pTempPreviousLabel = VarLabel::create( "p.tempPrevious",
                        ParticleVariable<double>::getTypeDescription() ); 
                     
  pdTdtLabel = VarLabel::create( "p.dTdt",
                        ParticleVariable<double>::getTypeDescription() );
  
  pExternalHeatRateLabel = VarLabel::create( "p.externalHeatRate",
                        ParticleVariable<double>::getTypeDescription() );

  pExternalHeatFluxLabel = VarLabel::create( "p.externalHeatFlux",
                        ParticleVariable<double>::getTypeDescription() );
  
  pSurfLabel = VarLabel::create( "p.surface",
                        ParticleVariable<double>::getTypeDescription() );

  pSurfGradLabel = VarLabel::create( "p.surfgrad",
                        ParticleVariable<Vector>::getTypeDescription() );

  pSurfGradLabel_preReloc = VarLabel::create( "p.surfgrad+",
                        ParticleVariable<Vector>::getTypeDescription() );

  pLastLevelLabel = VarLabel::create( "p.lastlevel",
                        ParticleVariable<int>::getTypeDescription() );

  pParticleIDLabel = VarLabel::create("p.particleID",
                        ParticleVariable<long64>::getTypeDescription() );

  pSizeLabel = VarLabel::create( "p.size",
                        ParticleVariable<Matrix3>::getTypeDescription());
  
  pSizeLabel_preReloc = VarLabel::create( "p.size+",
                        ParticleVariable<Matrix3>::getTypeDescription());
  
  pCurSizeLabel = VarLabel::create( "p.currentsize",
                        ParticleVariable<Matrix3>::getTypeDescription());
  
  pLocalizedMPMLabel  = VarLabel::create( "p.localizedMPM",
                        ParticleVariable<int>::getTypeDescription() );
  
  pLocalizedMPMLabel_preReloc  = VarLabel::create( "p.localizedMPM+",
                        ParticleVariable<int>::getTypeDescription() );
  
  pRefinedLabel       = VarLabel::create( "p.refinedMPM",
                        ParticleVariable<int>::getTypeDescription() );
  
  pRefinedLabel_preReloc  = VarLabel::create( "p.refinedMPM+",
                        ParticleVariable<int>::getTypeDescription() );

  pESPotential = VarLabel::create( "p.esPotential",
      ParticleVariable<double>::getTypeDescription());

  pESGradPotential = VarLabel::create( "p.esGradPotential",
                     ParticleVariable<Vector>::getTypeDescription() );

  pFiberDirLabel  = VarLabel::create( "p.fiberdir",
                        ParticleVariable<Vector>::getTypeDescription() );

  pFiberDirLabel_preReloc  = VarLabel::create( "p.fiberdir+",
                        ParticleVariable<Vector>::getTypeDescription() );

  // Particle Variables 
  pDeformationMeasureLabel_preReloc = VarLabel::create("p.deformationMeasure+",
                        ParticleVariable<Matrix3>::getTypeDescription());
  
  pVelGradLabel_preReloc = VarLabel::create("p.velGrad+",
                        ParticleVariable<Matrix3>::getTypeDescription());

  pStressLabel_preReloc = VarLabel::create( "p.stress+",
                        ParticleVariable<Matrix3>::getTypeDescription() );

  pStressVizualLabel = VarLabel::create("p.stressVizual",
      ParticleVariable<Matrix3>::getTypeDescription());

  pStressVizualLabel_preReloc = VarLabel::create("p.stressVizual+",
      ParticleVariable<Matrix3>::getTypeDescription());

  pVolumeLabel_preReloc = VarLabel::create( "p.volume+",
                        ParticleVariable<double>::getTypeDescription());
  
  pMassLabel_preReloc = VarLabel::create( "p.mass+",
                        ParticleVariable<double>::getTypeDescription() );
  
  pVelocityLabel_preReloc = VarLabel::create( "p.velocity+", 
                        ParticleVariable<Vector>::getTypeDescription() );
  
  pExtForceLabel_preReloc = VarLabel::create( "p.externalforce+",
                        ParticleVariable<Vector>::getTypeDescription() );
  
  pXLabel_preReloc = VarLabel::create( "p.x+",
                        ParticleVariable<Point>::getTypeDescription(),
                        IntVector(0,0,0),
                        VarLabel::PositionVariable);
  
  pTemperatureLabel_preReloc = VarLabel::create( "p.temperature+",
                        ParticleVariable<double>::getTypeDescription() );
    
  // for thermal stress
  pTempPreviousLabel_preReloc = VarLabel::create( "p.tempPrevious+",
                        ParticleVariable<double>::getTypeDescription() ); 

  pdTdtLabel_preReloc = VarLabel::create( "p.dTdt+",
                        ParticleVariable<double>::getTypeDescription() );

  pExternalHeatRateLabel_preReloc = VarLabel::create( "p.externalHeatRate+",
                        ParticleVariable<double>::getTypeDescription() );

  pExternalHeatFluxLabel_preReloc = VarLabel::create( "p.externalHeatFlux+",
                        ParticleVariable<double>::getTypeDescription() );
  
  pSurfLabel_preReloc = VarLabel::create( "p.surface+",
                        ParticleVariable<double>::getTypeDescription() );

  pLastLevelLabel_preReloc = VarLabel::create( "p.lastlevel+",
                        ParticleVariable<int>::getTypeDescription() );

  pParticleIDLabel_preReloc = VarLabel::create("p.particleID+",
                        ParticleVariable<long64>::getTypeDescription() );

  pTemperatureGradientLabel = VarLabel::create( "p.temperatureGradient",
                        ParticleVariable<Vector>::getTypeDescription() );

  pTemperatureGradientLabel_preReloc =VarLabel::create("p.temperatureGradient+",
                        ParticleVariable<Vector>::getTypeDescription() );

  pPosChargeLabel = VarLabel::create( "p.posCharge",
      ParticleVariable<double>::getTypeDescription() );

  pPosChargeLabel_preReloc = VarLabel::create( "p.posCharge+",
      ParticleVariable<double>::getTypeDescription() );

  pNegChargeLabel = VarLabel::create( "p.negCharge",
      ParticleVariable<double>::getTypeDescription() );

  pNegChargeLabel_preReloc = VarLabel::create( "p.negCharge+",
      ParticleVariable<double>::getTypeDescription() );

  pPosChargeGradLabel = VarLabel::create( "p.posChargeGrad",
      ParticleVariable<Vector>::getTypeDescription() );

  pPosChargeGradLabel_preReloc =VarLabel::create("p.posChargeGrad+",
      ParticleVariable<Vector>::getTypeDescription() );

  pNegChargeGradLabel = VarLabel::create( "p.negChargeGrad",
          ParticleVariable<Vector>::getTypeDescription() );

  pNegChargeGradLabel_preReloc =VarLabel::create("p.negChargeGrad+",
          ParticleVariable<Vector>::getTypeDescription() );

  pPermittivityLabel = VarLabel::create( "p.permittivity",
      ParticleVariable<double>::getTypeDescription() );

  pPermittivityLabel_preReloc = VarLabel::create( "p.permittivity+",
        ParticleVariable<double>::getTypeDescription() );

  pPosChargeFluxLabel = VarLabel::create("p.posChargeFlux",
        ParticleVariable<Vector>::getTypeDescription() );

  pPosChargeFluxLabel_preReloc = VarLabel::create("p.posChargeFlux+",
        ParticleVariable<Vector>::getTypeDescription() );

  pNegChargeFluxLabel = VarLabel::create("p.negChargeFlux",
        ParticleVariable<Vector>::getTypeDescription() );

  pNegChargeFluxLabel_preReloc = VarLabel::create("p.negChargeFlux+",
        ParticleVariable<Vector>::getTypeDescription() );

  // Node Centered Variables
  
  gAccelerationLabel = VarLabel::create( "g.acceleration",
                        NCVariable<Vector>::getTypeDescription() );

  gColorLabel = VarLabel::create( "g.color",
                        NCVariable<double>::getTypeDescription() );
  
  gMassLabel = VarLabel::create( "g.mass",
                        NCVariable<double>::getTypeDescription() );
  
  gMassAllLabel = VarLabel::create( "g.massall",
                        NCVariable<double>::getTypeDescription() );
  
  gMassF0Label = VarLabel::create( "g.massF0",
                        NCVariable<double>::getTypeDescription() );

  gMassF1Label = VarLabel::create( "g.massF1",
                        NCVariable<double>::getTypeDescription() );

  gVelocityLabel = VarLabel::create( "g.velocity",
                        NCVariable<Vector>::getTypeDescription() );
  
  gVelocityF0Label = VarLabel::create( "g.velocityF0",
                        NCVariable<Vector>::getTypeDescription() );
  
  gVelocityF1Label = VarLabel::create( "g.velocityF1",
                        NCVariable<Vector>::getTypeDescription() );
  
  gVelocityStarF0Label = VarLabel::create( "g.velocity_starF0",
                        NCVariable<Vector>::getTypeDescription() );
  
  gVelocityStarF1Label = VarLabel::create( "g.velocity_starF1",
                        NCVariable<Vector>::getTypeDescription() );
  
  gAccelerationF0Label = VarLabel::create( "g.accelerationF0",
                        NCVariable<Vector>::getTypeDescription() );
  
  gAccelerationF1Label = VarLabel::create( "g.accelerationF1",
                        NCVariable<Vector>::getTypeDescription() );
  
  gInternalForceF0Label = VarLabel::create( "g.internalforceF0",
                        NCVariable<Vector>::getTypeDescription() );
  
  gInternalForceF1Label = VarLabel::create( "g.internalforceF1",
                        NCVariable<Vector>::getTypeDescription() );
  
  gExternalForceF0Label = VarLabel::create( "g.externalforceF0",
                        NCVariable<Vector>::getTypeDescription() );
  
  gExternalForceF1Label = VarLabel::create( "g.externalforceF1",
                        NCVariable<Vector>::getTypeDescription() );
  
  gVelocityBCLabel = VarLabel::create( "g.velocityBC",
                        NCVariable<Vector>::getTypeDescription() );
  
  gVelSPSSPLabel = VarLabel::create( "g.velocitySPLusSSPlus",
                        NCVariable<Vector>::getTypeDescription() );
  
  gMatlProminenceLabel = VarLabel::create( "g.matlProminence",
                        NCVariable<double>::getTypeDescription() );
  
  gAlphaMaterialLabel = VarLabel::create( "g.alphaMaterial",
                        NCVariable<int>::getTypeDescription() );
  
  gNormAlphaToBetaLabel = VarLabel::create( "g.normAlphaToBeta",
                        NCVariable<Vector>::getTypeDescription() );
  
  gPositionLabel = VarLabel::create( "g.position",
                        NCVariable<Point>::getTypeDescription() );
  
  gPositionF0Label = VarLabel::create( "g.positionF0",
                        NCVariable<Point>::getTypeDescription() );
  
  gPositionF1Label = VarLabel::create( "g.positionF1",
                        NCVariable<Point>::getTypeDescription() );
  
  gExternalForceLabel = VarLabel::create( "g.externalforce",
                        NCVariable<Vector>::getTypeDescription() );
  
  gInternalForceLabel = VarLabel::create( "g.internalforce",
                        NCVariable<Vector>::getTypeDescription() );
  
  gContactLabel       = VarLabel::create( "g.contact",
                        NCVariable<int>::getTypeDescription() );
  
  gVelocityStarLabel = VarLabel::create( "g.velocity_star",
                        NCVariable<Vector>::getTypeDescription() );
  
  gTemperatureLabel = VarLabel::create("g.temperature",
                        NCVariable<double>::getTypeDescription());

  gTemperatureNoBCLabel = VarLabel::create("g.temperaturenobc",
                        NCVariable<double>::getTypeDescription());

  gTemperatureStarLabel = VarLabel::create("g.temperatureStar",
                        NCVariable<double>::getTypeDescription());

  gTemperatureRateLabel = VarLabel::create("g.temperatureRate",
                        NCVariable<double>::getTypeDescription());

  gdTdtLabel = VarLabel::create("g.dTdt",
                        NCVariable<double>::getTypeDescription());

  gHeatFluxLabel = VarLabel::create("g.HeatFlux",
                        NCVariable<Vector>::getTypeDescription());

  gExternalHeatRateLabel = VarLabel::create("g.externalHeatRate",
                        NCVariable<double>::getTypeDescription());

  gExternalHeatFluxLabel = VarLabel::create("g.externalHeatFlux",
                       NCVariable<double>::getTypeDescription());

  NC_CCweightLabel     = VarLabel::create("NC_CCweight",
                     NCVariable<double>::getTypeDescription());

  gThermalContactTemperatureRateLabel = 
     VarLabel::create("g.thermalContactTemperatureRate",
     NCVariable<double>::getTypeDescription());

  gNormTractionLabel = VarLabel::create( "g.normtraction",
                   NCVariable<double>::getTypeDescription() );

  gNormTractionF0Label = VarLabel::create( "g.normtractionF0",
                   NCVariable<double>::getTypeDescription() );

  gNormTractionF1Label = VarLabel::create( "g.normtractionF1",
                   NCVariable<double>::getTypeDescription() );

  gSurfNormLabel = VarLabel::create( "g.surfnorm",
                   NCVariable<Vector>::getTypeDescription() );

  gSurfNormF0Label = VarLabel::create( "g.surfnormF0",
                   NCVariable<Vector>::getTypeDescription() );

  gSurfNormF1Label = VarLabel::create( "g.surfnormF1",
                   NCVariable<Vector>::getTypeDescription() );

  gSurfLabel     = VarLabel::create( "g.surf",
                   NCVariable<double>::getTypeDescription() );

  gSurfGradLabel = VarLabel::create( "g.surfgrad",
                   NCVariable<Vector>::getTypeDescription() );

  gStressLabel   = VarLabel::create( "g.stress",
                   NCVariable<Matrix3>::getTypeDescription() );

  gStressF0Label   = VarLabel::create( "g.stressF0",
                   NCVariable<Matrix3>::getTypeDescription() );

  gStressF1Label   = VarLabel::create( "g.stressF1",
                   NCVariable<Matrix3>::getTypeDescription() );

  gStressForSavingLabel   = VarLabel::create( "g.stressFS",
                   NCVariable<Matrix3>::getTypeDescription() );

  gVolumeLabel     = VarLabel::create("g.volume",
                        NCVariable<double>::getTypeDescription());

  gVolumeF0Label   = VarLabel::create("g.volumeF0",
                        NCVariable<double>::getTypeDescription());

  gVolumeF1Label   = VarLabel::create("g.volumeF1",
                        NCVariable<double>::getTypeDescription());

  gZOILabel     = VarLabel::create("g.zoi",
                        NCVariable<Stencil7>::getTypeDescription());

  MPMRefineCellLabel  = VarLabel::create( "MPMRefineCell",
                     CCVariable<double>::getTypeDescription() );

  cVolumeLabel  = VarLabel::create( "c.volume",
                     CCVariable<double>::getTypeDescription() );

  numLocInCellLabel  = VarLabel::create( "NumLocalizedInCell",
                     CCVariable<int>::getTypeDescription() );

  numInCellLabel  = VarLabel::create( "NumInCell",
                     CCVariable<int>::getTypeDescription() );

  TotalVolumeDeformedLabel = VarLabel::create( "TotalVolumeDeformed",
                                 sum_vartype::getTypeDescription() );
  
  gradPAccNCLabel = VarLabel::create("gradPAccNC",
                        NCVariable<Vector>::getTypeDescription());

  dTdt_NCLabel     = VarLabel::create("dTdt_NC",
                        NCVariable<double>::getTypeDescription());

  massBurnFractionLabel  = VarLabel::create("massBurnFraction",
                        NCVariable<double>::getTypeDescription());
                     
  gSp_volLabel     =  VarLabel::create("g.sp_vol",
                        NCVariable<double>::getTypeDescription()); 
 
  gSp_vol_srcLabel =  VarLabel::create("g.sp_vol_src",
                        NCVariable<double>::getTypeDescription());

  // Interaction with Arches, Fluid Mechanics

  AccArchesNCLabel = VarLabel::create("AccArchesNC",
                        NCVariable<Vector>::getTypeDescription() );

  // Interaction with Arches, Heat Transfer

  heaTranSolid_NCLabel = VarLabel::create("heaTranSolid_NC",
                                         NCVariable<double>::getTypeDescription() );

  frictionalWorkLabel = VarLabel::create("frictionalWork",
                        NCVariable<double>::getTypeDescription());

  gNumNearParticlesLabel = VarLabel::create("NumNearParticles",
                        NCVariable<double>::getTypeDescription());

  gPosChargeLabel = VarLabel::create("g.PosCharge",
      NCVariable<double>::getTypeDescription());
  gPosChargeStarLabel = VarLabel::create("g.PosChargeStar",
        NCVariable<double>::getTypeDescription());
  gNegChargeLabel = VarLabel::create("g.NegCharge",
      NCVariable<double>::getTypeDescription());
  gNegChargeStarLabel = VarLabel::create("g.NegChargeStar",
          NCVariable<double>::getTypeDescription());
  gPosChargeNoBCLabel = VarLabel::create("g.PosChargeNoBC",
      NCVariable<double>::getTypeDescription());
  gNegChargeNoBCLabel = VarLabel::create("g.NegChargeNoBC",
      NCVariable<double>::getTypeDescription());

  gPosChargeRateLabel = VarLabel::create("g.PosChargeRate",
      NCVariable<double>::getTypeDescription());
  gNegChargeRateLabel = VarLabel::create("g.NegChargeRate",
        NCVariable<double>::getTypeDescription());

  // Reduction variables
  partCountLabel = VarLabel::create("particleCount",
                                   sumlong_vartype::getTypeDescription());

  StrainEnergyLabel = VarLabel::create( "StrainEnergy",
                        sum_vartype::getTypeDescription() );

  AccStrainEnergyLabel = VarLabel::create( "AccStrainEnergy",
                        max_vartype::getTypeDescription() );

  KineticEnergyLabel = VarLabel::create( "KineticEnergy",
                        sum_vartype::getTypeDescription() );

  ThermalEnergyLabel = VarLabel::create( "ThermalEnergy",
                        sum_vartype::getTypeDescription() );

  TotalMassLabel = VarLabel::create( "TotalMass",
                                 sum_vartype::getTypeDescription() );

  NeedAddMPMMaterialLabel = VarLabel::create("NeedAddMPMMaterial",
                                 sum_vartype::getTypeDescription());
  for(int iside=0;iside<6;iside++) {
      string label_name = Patch::getFaceName( (Patch::FaceType) iside ); // FIXME: assumes face indices
      
      BndyContactAreaLabel[iside] =
        VarLabel::create( std::string("BndyContactArea_"+label_name).c_str(),
                          sum_vartype::getTypeDescription() );
      BndyContactCellAreaLabel[iside] =
        VarLabel::create( std::string("BndyContactCellArea_"+label_name).c_str(),
                          sum_vartype::getTypeDescription() );
      BndyForceLabel[iside] =
        VarLabel::create( std::string("BndyForce_"+label_name).c_str(),
                          sumvec_vartype::getTypeDescription() );
      BndyTractionLabel[iside] =
        VarLabel::create( std::string("BndyTraction_"+label_name).c_str(),
                          sumvec_vartype::getTypeDescription() );
  }
  
  CenterOfMassPositionLabel = VarLabel::create( "CenterOfMassPosition",
                                 sumvec_vartype::getTypeDescription() );

  TotalMomentumLabel = VarLabel::create( "TotalMomentum",
                                 sumvec_vartype::getTypeDescription() );

  RigidReactionForceLabel = VarLabel::create( "RigidReactionForce",
                                 sumvec_vartype::getTypeDescription() );

  RigidReactionTorqueLabel = VarLabel::create( "RigidReactionTorque",
                                 sumvec_vartype::getTypeDescription() );

  TotalLocalizedParticleLabel = VarLabel::create("TotalLocalizedParticle",
                                   sumlong_vartype::getTypeDescription());


  // for assigning particle ids
  pCellNAPIDLabel =
    VarLabel::create("cellNAPID", CCVariable<int>::getTypeDescription());

  doMechLabel = VarLabel::create( "doMech", delt_vartype::getTypeDescription());

  // Implicit MPM labels

  gVelocityOldLabel = VarLabel::create("g.VelocityOld",
                                    NCVariable<Vector>::getTypeDescription()); 

  dispNewLabel = VarLabel::create("dispNew",
                                  NCVariable<Vector>::getTypeDescription()); 

  dispIncLabel = VarLabel::create("dispInc",
                                  NCVariable<Vector>::getTypeDescription()); 

  dispIncQNorm0 = VarLabel::create("dispIncQNorm0",
                                   sum_vartype::getTypeDescription());

  dispIncNormMax = VarLabel::create("dispIncNormMax",
                                    sum_vartype::getTypeDescription());

  dispIncQNorm = VarLabel::create("dispIncQNorm",
                                  sum_vartype::getTypeDescription());

  dispIncNorm = VarLabel::create("dispIncNorm",
                                 sum_vartype::getTypeDescription());
  
  pAccelerationLabel = VarLabel::create("p.acceleration",
                                   ParticleVariable<Vector>::getTypeDescription()); 

  pAccelerationLabel_preReloc = VarLabel::create("p.acceleration+",
                                   ParticleVariable<Vector>::getTypeDescription()); 

  // for Fracture ----------------------------
  pDispLabel = VarLabel::create("p.displacement",
                  ParticleVariable<Vector>::getTypeDescription());
  pDispLabel_preReloc = VarLabel::create( "p.displacement+",
                  ParticleVariable<Vector>::getTypeDescription());

  pDispGradsLabel = VarLabel::create("p.dispGrads",
                  ParticleVariable<Matrix3>::getTypeDescription());
  pDispGradsLabel_preReloc = VarLabel::create( "p.dispGrads+",
                  ParticleVariable<Matrix3>::getTypeDescription());

  pStrainEnergyDensityLabel = VarLabel::create("p.strainEnergyDensity",
                  ParticleVariable<double>::getTypeDescription());
  pStrainEnergyDensityLabel_preReloc = VarLabel::create( "p.strainEnergyDensity+",
                  ParticleVariable<double>::getTypeDescription());

  pgCodeLabel = VarLabel::create("p.gcode",
                  ParticleVariable<Short27>::getTypeDescription());

  pKineticEnergyDensityLabel = VarLabel::create("p.kineticEnergyDensity",
                  ParticleVariable<double>::getTypeDescription());

  gNumPatlsLabel = VarLabel::create("g.numPatls",
                  NCVariable<int>::getTypeDescription());

  GNumPatlsLabel = VarLabel::create("G.numPatls",
                  NCVariable<int>::getTypeDescription());

  gDisplacementLabel = VarLabel::create("g.displacement",
                  NCVariable<Vector>::getTypeDescription() );

  GDisplacementLabel = VarLabel::create("G.displacement",
                  NCVariable<Vector>::getTypeDescription() );

  gGridStressLabel = VarLabel::create("g.gridStress",
                  NCVariable<Matrix3>::getTypeDescription());
  GGridStressLabel = VarLabel::create("G.gridStress",
                  NCVariable<Matrix3>::getTypeDescription());

  gDispGradsLabel = VarLabel::create("g.dispGrads",
                  NCVariable<Matrix3>::getTypeDescription());
  GDispGradsLabel = VarLabel::create("G.dispGrads",
                  NCVariable<Matrix3>::getTypeDescription());

  gVelGradsLabel = VarLabel::create("g.velGrads",
                  NCVariable<Matrix3>::getTypeDescription());
  GVelGradsLabel = VarLabel::create("G.velGrads",
                  NCVariable<Matrix3>::getTypeDescription());

  gStrainEnergyDensityLabel = VarLabel::create("g.strainEnergyDensity",
                  NCVariable<double>::getTypeDescription());
  GStrainEnergyDensityLabel = VarLabel::create("G.strainEnergyDensity",
                  NCVariable<double>::getTypeDescription());

  gKineticEnergyDensityLabel = VarLabel::create("g.kineticEnergyDensity",
                  NCVariable<double>::getTypeDescription());
  GKineticEnergyDensityLabel = VarLabel::create("G.kineticEnergyDensity",
                  NCVariable<double>::getTypeDescription());

  GCrackNormLabel= VarLabel::create("G.cracknormal",
                  NCVariable<Vector>::getTypeDescription());

  GMassLabel = VarLabel::create("G.mass",
                  NCVariable<double>::getTypeDescription());

  GVolumeLabel = VarLabel::create("G.volume",
                  NCVariable<double>::getTypeDescription());

  GVelocityLabel = VarLabel::create("G.velocity",
                  NCVariable<Vector>::getTypeDescription());
  
  GTemperatureLabel = VarLabel::create("G.temperature",
                  NCVariable<double>::getTypeDescription());

  GTemperatureNoBCLabel = VarLabel::create("G.temperatureiNoBC",
                  NCVariable<double>::getTypeDescription());

  GExternalForceLabel = VarLabel::create("G.externalforce",
                  NCVariable<Vector>::getTypeDescription());

  GExternalHeatRateLabel = VarLabel::create("G.externalheatrate",
                  NCVariable<double>::getTypeDescription());

  GThermalContactTemperatureRateLabel = VarLabel::create(
                  "G.thermalContactTemperatureRate",
                  NCVariable<double>::getTypeDescription());

  GInternalForceLabel = VarLabel::create("G.internalforce",
                  NCVariable<Vector>::getTypeDescription());

  GdTdtLabel = VarLabel::create("G.dTdt",
                  NCVariable<double>::getTypeDescription());   

  GTemperatureRateLabel = VarLabel::create("G.temperatureRate",
                  NCVariable<double>::getTypeDescription());
     
  GTemperatureStarLabel = VarLabel::create("G.temperatureStar",
                  NCVariable<double>::getTypeDescription());
                   
  GVelocityStarLabel = VarLabel::create("G.velocityg_star",
                  NCVariable<Vector>::getTypeDescription() );

  GAccelerationLabel = VarLabel::create( "G.acceleration",
                  NCVariable<Vector>::getTypeDescription() );

  GSp_volLabel = VarLabel::create("G.sp_vol",
                  NCVariable<double>::getTypeDescription());

  GSp_vol_srcLabel =  VarLabel::create("G.sp_vol_src",
                  NCVariable<double>::getTypeDescription()); 
  // ------------------------------------------------------

  // Material point erosion algorithms
  pErosionLabel = VarLabel::create("p.erosion",
                            ParticleVariable<double>::getTypeDescription());
  pErosionLabel_preReloc = VarLabel::create("p.erosion+",
                            ParticleVariable<double>::getTypeDescription());
  
  // MPM Physical BC labels (permanent particle state)
  materialPointsPerLoadCurveLabel = VarLabel::create("pointsPerCurve", 
                            sumlong_vartype::getTypeDescription());
  pLoadCurveIDLabel = VarLabel::create("p.loadCurveID",
                            ParticleVariable<IntVector>::getTypeDescription());
  pLoadCurveIDLabel_preReloc = VarLabel::create("p.loadCurveID+",
                            ParticleVariable<IntVector>::getTypeDescription());

  // Labels for shell materials
  pThickTopLabel = VarLabel::create("p.thickTop",
                     ParticleVariable<double>::getTypeDescription());
  pInitialThickTopLabel = VarLabel::create("p.thickTop0",
                     ParticleVariable<double>::getTypeDescription());
  pThickBotLabel = VarLabel::create("p.thickBot",
                     ParticleVariable<double>::getTypeDescription());
  pInitialThickBotLabel = VarLabel::create("p.thickBot0",
                     ParticleVariable<double>::getTypeDescription());
  pNormalLabel = VarLabel::create("p.normal",
                     ParticleVariable<Vector>::getTypeDescription());
  pInitialNormalLabel = VarLabel::create("p.normal0",
                     ParticleVariable<Vector>::getTypeDescription());

  pThickTopLabel_preReloc = VarLabel::create("p.thickTop+",
                     ParticleVariable<double>::getTypeDescription());
  pInitialThickTopLabel_preReloc = VarLabel::create("p.thickTop0+",
                     ParticleVariable<double>::getTypeDescription());
  pThickBotLabel_preReloc = VarLabel::create("p.thickBot+",
                     ParticleVariable<double>::getTypeDescription());
  pInitialThickBotLabel_preReloc = VarLabel::create("p.thickBot0+",
                     ParticleVariable<double>::getTypeDescription());
  pNormalLabel_preReloc = VarLabel::create("p.normal+",
                     ParticleVariable<Vector>::getTypeDescription());
  pInitialNormalLabel_preReloc = VarLabel::create("p.normal0+",
                     ParticleVariable<Vector>::getTypeDescription());

  pTypeLabel = VarLabel::create("p.type",
                     ParticleVariable<int>::getTypeDescription());
  pTypeLabel_preReloc = VarLabel::create("p.type+",
                     ParticleVariable<int>::getTypeDescription());

  gNormalRotRateLabel = VarLabel::create( "g.normalRotRate",
                     NCVariable<Vector>::getTypeDescription());
  gNormalRotMomentLabel = VarLabel::create( "g.normalRotMoment",
                     NCVariable<Vector>::getTypeDescription());

  gNormalRotMassLabel = VarLabel::create( "g.normalRotMass",
                     NCVariable<double>::getTypeDescription());
  gNormalRotAccLabel = VarLabel::create( "g.normalRotAcc",
                     NCVariable<Vector>::getTypeDescription());

  // For Cohesive Zones
  czAreaLabel = VarLabel::create("cz.area",
                     ParticleVariable<double>::getTypeDescription());
  czAreaLabel_preReloc = VarLabel::create("cz.area+",
                     ParticleVariable<double>::getTypeDescription());

  czNormLabel = VarLabel::create("cz.norm",
                     ParticleVariable<Vector>::getTypeDescription());
  czNormLabel_preReloc = VarLabel::create("cz.norm+",
                     ParticleVariable<Vector>::getTypeDescription());

  czTangLabel = VarLabel::create("cz.tang",
                     ParticleVariable<Vector>::getTypeDescription());
  czTangLabel_preReloc = VarLabel::create("cz.tang+",
                     ParticleVariable<Vector>::getTypeDescription());

  czDispTopLabel = VarLabel::create("cz.disptop",
                     ParticleVariable<Vector>::getTypeDescription());
  czDispTopLabel_preReloc = VarLabel::create("cz.disptop+",
                     ParticleVariable<Vector>::getTypeDescription());

  czDispBottomLabel = VarLabel::create("cz.dispbottom",
                     ParticleVariable<Vector>::getTypeDescription());
  czDispBottomLabel_preReloc = VarLabel::create("cz.dispbottom+",
                     ParticleVariable<Vector>::getTypeDescription());

  czSeparationLabel = VarLabel::create("cz.separation",
                     ParticleVariable<Vector>::getTypeDescription());
  czSeparationLabel_preReloc = VarLabel::create("cz.separation+",
                     ParticleVariable<Vector>::getTypeDescription());

  czForceLabel = VarLabel::create("cz.force",
                     ParticleVariable<Vector>::getTypeDescription());
  czForceLabel_preReloc = VarLabel::create("cz.force+",
                     ParticleVariable<Vector>::getTypeDescription());

  czTopMatLabel = VarLabel::create("cz.topmat",
                     ParticleVariable<int>::getTypeDescription());
  czTopMatLabel_preReloc = VarLabel::create("cz.topmat+",
                     ParticleVariable<int>::getTypeDescription());

  czBotMatLabel = VarLabel::create("cz.botmat",
                     ParticleVariable<int>::getTypeDescription());
  czBotMatLabel_preReloc = VarLabel::create("cz.botmat+",
                     ParticleVariable<int>::getTypeDescription());

  czFailedLabel = VarLabel::create("cz.failed",
                     ParticleVariable<int>::getTypeDescription());
  czFailedLabel_preReloc = VarLabel::create("cz.failed+",
                     ParticleVariable<int>::getTypeDescription());

  czIDLabel = VarLabel::create("cz.CZID",
                        ParticleVariable<long64>::getTypeDescription() );

  czIDLabel_preReloc = VarLabel::create("cz.CZID+",
                        ParticleVariable<long64>::getTypeDescription() );

  // for assigning particle ids
  pCellNACZIDLabel =
    VarLabel::create("cellNACZID", CCVariable<short int>::getTypeDescription());

  // Hydro-mechanical coupling MPM
  ccPorosity =
      VarLabel::create("ccPorosity", CCVariable<double>::getTypeDescription());
  ccPorePressure = VarLabel::create("ccPorePressure",
      CCVariable<double>::getTypeDescription());
  ccPorePressureOld = VarLabel::create(
      "ccPorePressureOld", CCVariable<double>::getTypeDescription());
  ccRHS_FlowEquation = VarLabel::create(
      "ccRHS_FlowEquation", CCVariable<double>::getTypeDescription());
  ccTransmissivityMatrix = VarLabel::create(
      "ccTransmissivityMatrix", CCVariable<Stencil7>::getTypeDescription());

  pFluidMassLabel = VarLabel::create(
      "p.fluidmass", ParticleVariable<double>::getTypeDescription());
  pSolidMassLabel = VarLabel::create(
      "p.solidmass", ParticleVariable<double>::getTypeDescription());
  pPorosityLabel = VarLabel::create(
      "p.porosity", ParticleVariable<double>::getTypeDescription());
  pPorePressureLabel = VarLabel::create(
      "p.porepressure", ParticleVariable<double>::getTypeDescription());
  pPorePressureFilterLabel = VarLabel::create(
      "p.porepressureFilter", ParticleVariable<double>::getTypeDescription());
  pStressRateLabel = VarLabel::create(
      "p.stressrate", ParticleVariable<double>::getTypeDescription());
  pFluidVelocityLabel = VarLabel::create(
      "p.fluidvelocity", ParticleVariable<Vector>::getTypeDescription());
  pFluidAccelerationLabel = VarLabel::create(
      "p.fluidacceleration", ParticleVariable<Vector>::getTypeDescription());

  pPrescribedPorePressureLabel =
      VarLabel::create("p.prescribedporepressure",
          ParticleVariable<Vector>::getTypeDescription());

  gFluidMassLabel =
      VarLabel::create("g.fluidmass", NCVariable<double>::getTypeDescription());
  gFluidMassBarLabel = VarLabel::create(
      "g.fluidmassbar", NCVariable<double>::getTypeDescription());

  gFluidVelocityLabel = VarLabel::create(
      "g.fluidvelocity", NCVariable<Vector>::getTypeDescription());
  FluidVelInc =
      VarLabel::create("fluidvelinc", NCVariable<Vector>::getTypeDescription());

  gFluidVelocityStarLabel = VarLabel::create(
      "g.fluidvelocitystar", NCVariable<Vector>::getTypeDescription());
  gFluidAccelerationLabel = VarLabel::create(
      "g.fluidacceleration", NCVariable<Vector>::getTypeDescription());
  gExternalFluidForceLabel = VarLabel::create(
      "g.externalfluidforce", NCVariable<Vector>::getTypeDescription());
  gInternalFluidForceLabel = VarLabel::create(
      "g.internalfluidforce", NCVariable<Vector>::getTypeDescription());
  gInternalDragForceLabel = VarLabel::create(
      "g.internaldragforce", NCVariable<Vector>::getTypeDescription());
  gFlowInertiaForceLabel = VarLabel::create(
      "g.flowinertiaforce", NCVariable<Vector>::getTypeDescription());

  pFluidMassLabel_preReloc = VarLabel::create(
      "p.fluidmass+", ParticleVariable<double>::getTypeDescription());
  pSolidMassLabel_preReloc = VarLabel::create(
      "p.solidmass+", ParticleVariable<double>::getTypeDescription());
  pPorosityLabel_preReloc = VarLabel::create(
      "p.porosity+", ParticleVariable<double>::getTypeDescription());
  pPorePressureLabel_preReloc = VarLabel::create(
      "p.porepressure+", ParticleVariable<double>::getTypeDescription());
  pPorePressureFilterLabel_preReloc = VarLabel::create(
      "p.porepressureFilter+", ParticleVariable<double>::getTypeDescription());
  pStressRateLabel_preReloc = VarLabel::create(
      "p.stressrate+", ParticleVariable<double>::getTypeDescription());
  pFluidVelocityLabel_preReloc = VarLabel::create(
      "p.fluidvelocity+", ParticleVariable<Vector>::getTypeDescription());
  pFluidAccelerationLabel_preReloc = VarLabel::create(
      "p.fluidacceleration+", ParticleVariable<Vector>::getTypeDescription());

  gFluidMassLabel_preReloc = VarLabel::create(
      "g.fluidmass+", NCVariable<double>::getTypeDescription());
  gFluidMassBarLabel_preReloc = VarLabel::create(
      "g.fluidmassbar+", NCVariable<double>::getTypeDescription());
  gFluidVelocityLabel_preReloc = VarLabel::create(
      "g.fluidvelocity+", NCVariable<Vector>::getTypeDescription());
  gFluidVelocityStarLabel_preReloc = VarLabel::create(
      "g.fluidvelocitystar+", NCVariable<Vector>::getTypeDescription());
  gFluidAccelerationLabel_preReloc = VarLabel::create(
      "g.fluidacceleration+", NCVariable<Vector>::getTypeDescription());
  gPorePressureLabel = VarLabel::create(
      "g.PorePressure", NCVariable<double>::getTypeDescription());
  gPorePressureFilterLabel = VarLabel::create(
      "g.PorePressureFilter", NCVariable<double>::getTypeDescription());

  // MPM Hydrostatic BC label
  boundaryPointsPerCellLabel = VarLabel::create(
      "boundarypointsPerCell", CCVariable<int>::getTypeDescription());
} 

MPMLabel::~MPMLabel()
{
  if (diffusion) {
    delete diffusion;
  }

  VarLabel::destroy(timeStepLabel);
  VarLabel::destroy(simulationTimeLabel);
  VarLabel::destroy(delTLabel);

  VarLabel::destroy(fAllParticlesUpdated);

  VarLabel::destroy(heatRate_CCLabel);
  //non PermanentParticleState
  VarLabel::destroy(pVolumeDeformedLabel);
  VarLabel::destroy(pTempCurrentLabel); // for thermal stress
  VarLabel::destroy(pXXLabel);

  //PermanentParticleState
  VarLabel::destroy(pDeformationMeasureLabel);
  VarLabel::destroy(pDeformationMeasureLabel_preReloc);
  VarLabel::destroy(pVelGradLabel);
  VarLabel::destroy(pVelGradLabel_preReloc);
  VarLabel::destroy(pStressLabel);
  VarLabel::destroy(pStressLabel_preReloc);
  VarLabel::destroy(pStressVizualLabel);
  VarLabel::destroy(pStressVizualLabel_preReloc);
  VarLabel::destroy(pVolumeLabel);
  VarLabel::destroy(pVolumeLabel_preReloc);
  VarLabel::destroy(pMassLabel);
  VarLabel::destroy(pMassLabel_preReloc);
  VarLabel::destroy(pVelocityLabel);
  VarLabel::destroy(pVelocityLabel_preReloc);
  VarLabel::destroy(pVelocitySSPlusLabel);
  VarLabel::destroy(pExternalForceLabel);
  VarLabel::destroy(pExternalForceCorner1Label);
  VarLabel::destroy(pExternalForceCorner2Label);
  VarLabel::destroy(pExternalForceCorner3Label);
  VarLabel::destroy(pExternalForceCorner4Label);
  VarLabel::destroy(pExtForceLabel_preReloc);
  VarLabel::destroy(pXLabel);
  VarLabel::destroy(pXLabel_preReloc);
  VarLabel::destroy(pTemperatureLabel);
  VarLabel::destroy(pTemperatureLabel_preReloc);
  VarLabel::destroy(pTempPreviousLabel); // for thermal stress
  VarLabel::destroy(pTempPreviousLabel_preReloc); // for thermal stress
  VarLabel::destroy(pdTdtLabel);
  VarLabel::destroy(pdTdtLabel_preReloc);
  VarLabel::destroy(pExternalHeatRateLabel);
  VarLabel::destroy(pExternalHeatRateLabel_preReloc);
  VarLabel::destroy(pExternalHeatFluxLabel);
  VarLabel::destroy(pExternalHeatFluxLabel_preReloc);
  VarLabel::destroy(pSurfLabel);
  VarLabel::destroy(pSurfLabel_preReloc);
  VarLabel::destroy(pSurfGradLabel);
  VarLabel::destroy(pSurfGradLabel_preReloc);
  VarLabel::destroy(pLastLevelLabel);
  VarLabel::destroy(pLastLevelLabel_preReloc);
  VarLabel::destroy(pParticleIDLabel);
  VarLabel::destroy(pParticleIDLabel_preReloc);
  VarLabel::destroy(czIDLabel);
  VarLabel::destroy(czIDLabel_preReloc);
  VarLabel::destroy(pPressureLabel);
  VarLabel::destroy(pScratchVecLabel);
  VarLabel::destroy(pScratchLabel);
  VarLabel::destroy(pScaleFactorLabel);
  VarLabel::destroy(pScaleFactorLabel_preReloc);
  VarLabel::destroy(pLocalizedMPMLabel);
  VarLabel::destroy(pLocalizedMPMLabel_preReloc);
  VarLabel::destroy(pRefinedLabel);
  VarLabel::destroy(pRefinedLabel_preReloc);
  VarLabel::destroy(pSizeLabel);
  VarLabel::destroy(pSizeLabel_preReloc);
  VarLabel::destroy(pCurSizeLabel);
  VarLabel::destroy(pFiberDirLabel_preReloc);
  VarLabel::destroy(pFiberDirLabel);
  VarLabel::destroy(p_qLabel);
  VarLabel::destroy(p_qLabel_preReloc);
  VarLabel::destroy(pTemperatureGradientLabel);
  VarLabel::destroy(pTemperatureGradientLabel_preReloc);
  VarLabel::destroy(pESPotential);
  VarLabel::destroy(pESGradPotential);
  VarLabel::destroy(pPartitionUnityLabel);
  VarLabel::destroy(pPosChargeLabel);
  VarLabel::destroy(pPosChargeLabel_preReloc);
  VarLabel::destroy(pNegChargeLabel);
  VarLabel::destroy(pNegChargeLabel_preReloc);
  VarLabel::destroy(pPosChargeGradLabel);
  VarLabel::destroy(pPosChargeGradLabel_preReloc);
  VarLabel::destroy(pNegChargeGradLabel);
  VarLabel::destroy(pNegChargeGradLabel_preReloc);
  VarLabel::destroy(pPermittivityLabel);
  VarLabel::destroy(pPermittivityLabel_preReloc);
  VarLabel::destroy(pPosChargeFluxLabel);
  VarLabel::destroy(pPosChargeFluxLabel_preReloc);
  VarLabel::destroy(pNegChargeFluxLabel);
  VarLabel::destroy(pNegChargeFluxLabel_preReloc);

  VarLabel::destroy(gAccelerationLabel);
  VarLabel::destroy(gColorLabel);
  VarLabel::destroy(gMassLabel);
  VarLabel::destroy(gMassAllLabel);
  VarLabel::destroy(gMassF0Label);
  VarLabel::destroy(gMassF1Label);
  VarLabel::destroy(gVelocityF0Label);
  VarLabel::destroy(gVelocityF1Label);
  VarLabel::destroy(gVelocityStarF0Label);
  VarLabel::destroy(gVelocityStarF1Label);
  VarLabel::destroy(gInternalForceF0Label);
  VarLabel::destroy(gInternalForceF1Label);
  VarLabel::destroy(gExternalForceF0Label);
  VarLabel::destroy(gExternalForceF1Label);
  VarLabel::destroy(gAccelerationF0Label);
  VarLabel::destroy(gAccelerationF1Label);
  VarLabel::destroy(gVelocityLabel);
  VarLabel::destroy(gVelocityBCLabel);
  VarLabel::destroy(gVelSPSSPLabel);
  VarLabel::destroy(gMatlProminenceLabel);
  VarLabel::destroy(gAlphaMaterialLabel);
  VarLabel::destroy(gNormAlphaToBetaLabel);
  VarLabel::destroy(gPositionLabel);
  VarLabel::destroy(gPositionF0Label);
  VarLabel::destroy(gPositionF1Label);
  VarLabel::destroy(gExternalForceLabel);
  VarLabel::destroy(gInternalForceLabel);
  VarLabel::destroy(gContactLabel);
  VarLabel::destroy(gVelocityStarLabel);
  VarLabel::destroy(gNormTractionLabel);
  VarLabel::destroy(gNormTractionF0Label);
  VarLabel::destroy(gNormTractionF1Label);
  VarLabel::destroy(gStressLabel);
  VarLabel::destroy(gStressVizualLabel);
  VarLabel::destroy(gStressF0Label);
  VarLabel::destroy(gStressF1Label);
  VarLabel::destroy(gSurfNormLabel);
  VarLabel::destroy(gSurfNormF0Label);
  VarLabel::destroy(gSurfNormF1Label);
  VarLabel::destroy(gSurfLabel);
  VarLabel::destroy(gSurfGradLabel);
  VarLabel::destroy(gTemperatureLabel);
  VarLabel::destroy(gSp_volLabel); 
  VarLabel::destroy(gSp_vol_srcLabel); 
  VarLabel::destroy(gTemperatureNoBCLabel);
  VarLabel::destroy(gTemperatureStarLabel);
  VarLabel::destroy(gTemperatureRateLabel);
  VarLabel::destroy(gdTdtLabel);
  VarLabel::destroy(gHeatFluxLabel);
  VarLabel::destroy(gExternalHeatRateLabel);
  VarLabel::destroy(gExternalHeatFluxLabel);
  VarLabel::destroy(NC_CCweightLabel);
  VarLabel::destroy(gThermalContactTemperatureRateLabel);
  VarLabel::destroy(gStressForSavingLabel);
  VarLabel::destroy(gVolumeLabel);
  VarLabel::destroy(gVolumeF0Label);
  VarLabel::destroy(gVolumeF1Label);
  VarLabel::destroy(gZOILabel);
  VarLabel::destroy(MPMRefineCellLabel);
  VarLabel::destroy(cVolumeLabel);
  VarLabel::destroy(numLocInCellLabel);
  VarLabel::destroy(numInCellLabel);
  VarLabel::destroy(gradPAccNCLabel);
  VarLabel::destroy(dTdt_NCLabel);
  VarLabel::destroy(massBurnFractionLabel);
  VarLabel::destroy(AccArchesNCLabel);
  VarLabel::destroy(heaTranSolid_NCLabel);
  VarLabel::destroy(frictionalWorkLabel);
  VarLabel::destroy(gNumNearParticlesLabel);

  VarLabel::destroy(gPosChargeLabel);
  VarLabel::destroy(gPosChargeStarLabel);
  VarLabel::destroy(gNegChargeLabel);
  VarLabel::destroy(gNegChargeStarLabel);
  VarLabel::destroy(gPosChargeNoBCLabel);
  VarLabel::destroy(gNegChargeNoBCLabel);
  VarLabel::destroy(gPosChargeRateLabel);
  VarLabel::destroy(gNegChargeRateLabel);

  VarLabel::destroy(partCountLabel);
  VarLabel::destroy(doMechLabel);

  VarLabel::destroy(AccStrainEnergyLabel);
  VarLabel::destroy(StrainEnergyLabel);
  VarLabel::destroy(KineticEnergyLabel);
  VarLabel::destroy(ThermalEnergyLabel);
  VarLabel::destroy(TotalMassLabel);
  VarLabel::destroy(NeedAddMPMMaterialLabel);
  VarLabel::destroy(TotalVolumeDeformedLabel);
  for(int iside=0;iside<6;iside++) {
      VarLabel::destroy(BndyContactAreaLabel[iside]);
      VarLabel::destroy(BndyContactCellAreaLabel[iside]);
      VarLabel::destroy(BndyForceLabel[iside]);
      VarLabel::destroy(BndyTractionLabel[iside]);
  }
  VarLabel::destroy(CenterOfMassPositionLabel);
  VarLabel::destroy(TotalMomentumLabel);
  VarLabel::destroy(RigidReactionForceLabel);
  VarLabel::destroy(RigidReactionTorqueLabel);
  VarLabel::destroy(TotalLocalizedParticleLabel);
  VarLabel::destroy(pCellNAPIDLabel);
  VarLabel::destroy(pCellNACZIDLabel);

  VarLabel::destroy(gVelocityOldLabel);
  VarLabel::destroy(dispNewLabel);
  VarLabel::destroy(dispIncLabel);
  VarLabel::destroy(pAccelerationLabel);
  VarLabel::destroy(dispIncQNorm0);
  VarLabel::destroy(dispIncNormMax);
  VarLabel::destroy(dispIncQNorm);
  VarLabel::destroy(dispIncNorm);
  VarLabel::destroy(pAccelerationLabel_preReloc);

 // for Fracture --------------
  VarLabel::destroy(pDispLabel);
  VarLabel::destroy(pDispLabel_preReloc);
  VarLabel::destroy(pDispGradsLabel);
  VarLabel::destroy(pDispGradsLabel_preReloc);
  VarLabel::destroy(pStrainEnergyDensityLabel);
  VarLabel::destroy(pStrainEnergyDensityLabel_preReloc);
  VarLabel::destroy(pKineticEnergyDensityLabel);

  VarLabel::destroy(pgCodeLabel);
  VarLabel::destroy(gNumPatlsLabel);
  VarLabel::destroy(GNumPatlsLabel);
  VarLabel::destroy(gDisplacementLabel);
  VarLabel::destroy(GDisplacementLabel);
  VarLabel::destroy(gGridStressLabel);
  VarLabel::destroy(GGridStressLabel);
  VarLabel::destroy(gDispGradsLabel);
  VarLabel::destroy(GDispGradsLabel);
  VarLabel::destroy(gVelGradsLabel);
  VarLabel::destroy(GVelGradsLabel);
  VarLabel::destroy(gStrainEnergyDensityLabel);
  VarLabel::destroy(GStrainEnergyDensityLabel);
  VarLabel::destroy(gKineticEnergyDensityLabel);
  VarLabel::destroy(GKineticEnergyDensityLabel);

  VarLabel::destroy(GCrackNormLabel);
  VarLabel::destroy(GMassLabel);
  VarLabel::destroy(GVolumeLabel);
  VarLabel::destroy(GVelocityLabel);
  VarLabel::destroy(GTemperatureLabel);
  VarLabel::destroy(GTemperatureNoBCLabel);
  VarLabel::destroy(GExternalForceLabel);
  VarLabel::destroy(GExternalHeatRateLabel);
  VarLabel::destroy(GThermalContactTemperatureRateLabel);
  VarLabel::destroy(GInternalForceLabel);
  VarLabel::destroy(GdTdtLabel);
  VarLabel::destroy(GTemperatureRateLabel);
  VarLabel::destroy(GTemperatureStarLabel);
  VarLabel::destroy(GVelocityStarLabel);
  VarLabel::destroy(GAccelerationLabel);
  VarLabel::destroy(GSp_volLabel);
  VarLabel::destroy(GSp_vol_srcLabel);
  // --------------------------------

  // Destroy Material point erosion labels
  VarLabel::destroy(pErosionLabel);
  VarLabel::destroy(pErosionLabel_preReloc);

  // Destroy the MPM Physical BC pointer labels
  VarLabel::destroy(materialPointsPerLoadCurveLabel);
  VarLabel::destroy(pLoadCurveIDLabel);
  VarLabel::destroy(pLoadCurveIDLabel_preReloc);

  // Destroy Labels for shell materials
  VarLabel::destroy(pThickTopLabel);
  VarLabel::destroy(pInitialThickTopLabel);
  VarLabel::destroy(pThickBotLabel);
  VarLabel::destroy(pInitialThickBotLabel);
  VarLabel::destroy(pNormalLabel);
  VarLabel::destroy(pInitialNormalLabel);

  VarLabel::destroy(pThickTopLabel_preReloc);
  VarLabel::destroy(pInitialThickTopLabel_preReloc);
  VarLabel::destroy(pThickBotLabel_preReloc);
  VarLabel::destroy(pInitialThickBotLabel_preReloc);
  VarLabel::destroy(pNormalLabel_preReloc);
  VarLabel::destroy(pInitialNormalLabel_preReloc);

  VarLabel::destroy(pTypeLabel);
  VarLabel::destroy(pTypeLabel_preReloc);

  VarLabel::destroy(gNormalRotRateLabel); 
  VarLabel::destroy(gNormalRotMomentLabel); 
  VarLabel::destroy(gNormalRotMassLabel); 
  VarLabel::destroy(gNormalRotAccLabel); 
  
   // Debugging labels
  VarLabel::destroy(pColorLabel);
  VarLabel::destroy(pColorLabel_preReloc);

  // For Cohesive Zones
  VarLabel::destroy(czAreaLabel);
  VarLabel::destroy(czAreaLabel_preReloc);
  VarLabel::destroy(czNormLabel);
  VarLabel::destroy(czNormLabel_preReloc);
  VarLabel::destroy(czTangLabel);
  VarLabel::destroy(czTangLabel_preReloc);
  VarLabel::destroy(czDispTopLabel);
  VarLabel::destroy(czDispTopLabel_preReloc);
  VarLabel::destroy(czDispBottomLabel);
  VarLabel::destroy(czDispBottomLabel_preReloc);
  VarLabel::destroy(czSeparationLabel);
  VarLabel::destroy(czSeparationLabel_preReloc);
  VarLabel::destroy(czForceLabel);
  VarLabel::destroy(czForceLabel_preReloc);
  VarLabel::destroy(czTopMatLabel);
  VarLabel::destroy(czTopMatLabel_preReloc);
  VarLabel::destroy(czBotMatLabel);
  VarLabel::destroy(czBotMatLabel_preReloc);
  VarLabel::destroy(czFailedLabel);
  VarLabel::destroy(czFailedLabel_preReloc);

  // Hydro-mechanical coupling MPM
  VarLabel::destroy(ccPorePressure);
  VarLabel::destroy(ccPorosity);
  VarLabel::destroy(ccPorePressureOld);
  VarLabel::destroy(ccRHS_FlowEquation);
  VarLabel::destroy(ccTransmissivityMatrix);
  VarLabel::destroy(pFluidVelocityLabel);
  VarLabel::destroy(pFluidAccelerationLabel);
  VarLabel::destroy(pSolidMassLabel);
  VarLabel::destroy(pFluidMassLabel);
  VarLabel::destroy(pPorosityLabel);
  VarLabel::destroy(pStressRateLabel);
  VarLabel::destroy(pPorosityLabel_preReloc);
  VarLabel::destroy(pStressRateLabel_preReloc);

  VarLabel::destroy(gFluidMassBarLabel);
  VarLabel::destroy(gFluidMassLabel);
  VarLabel::destroy(pPrescribedPorePressureLabel);
  VarLabel::destroy(pPorePressureLabel);
  VarLabel::destroy(pPorePressureLabel_preReloc);
  VarLabel::destroy(pPorePressureFilterLabel);
  VarLabel::destroy(pPorePressureFilterLabel_preReloc);
  VarLabel::destroy(gFluidVelocityLabel);
  VarLabel::destroy(FluidVelInc);
  VarLabel::destroy(gFluidAccelerationLabel);
  VarLabel::destroy(gInternalDragForceLabel);
  VarLabel::destroy(gExternalFluidForceLabel);
  VarLabel::destroy(gInternalFluidForceLabel);
  VarLabel::destroy(gFlowInertiaForceLabel);
  VarLabel::destroy(gPorePressureLabel);
  VarLabel::destroy(gPorePressureFilterLabel);

  VarLabel::destroy(pFluidMassLabel_preReloc);
  VarLabel::destroy(pFluidVelocityLabel_preReloc);
  VarLabel::destroy(pFluidAccelerationLabel_preReloc);
  VarLabel::destroy(pSolidMassLabel_preReloc);
  VarLabel::destroy(gFluidMassBarLabel_preReloc);
  VarLabel::destroy(gFluidMassLabel_preReloc);
  VarLabel::destroy(gFluidVelocityLabel_preReloc);
  VarLabel::destroy(gFluidAccelerationLabel_preReloc);

  VarLabel::destroy(gFluidVelocityStarLabel);
  VarLabel::destroy(gFluidVelocityStarLabel_preReloc);

  // MPM Hydrostatic BC label
  VarLabel::destroy(boundaryPointsPerCellLabel);
}
