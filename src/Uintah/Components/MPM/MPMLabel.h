#ifndef UINTAH_HOMEBREW_MPMLABEL_H
#define UINTAH_HOMEBREW_MPMLABEL_H


#include <Uintah/Grid/VarLabel.h>

namespace Uintah {
  namespace MPM {
    class MPMLabel {
    public:

      MPMLabel();
      ~MPMLabel();

      static const MPMLabel* getLabels();

      const VarLabel* delTLabel;
      
      const VarLabel* pDeformationMeasureLabel;
      const VarLabel* pStressLabel;
      const VarLabel* pVolumeLabel;
      const VarLabel* pVolumeDeformedLabel;
      const VarLabel* pMassLabel;
      const VarLabel* pVelocityLabel;
      const VarLabel* pExternalForceLabel;
      const VarLabel* pXLabel;
      const VarLabel* pSurfLabel;
      const VarLabel* pSurfaceNormalLabel; //for fracture
      const VarLabel* pAverageMicrocrackLength; //for fracture
      const VarLabel* pTemperatureLabel; //for heat conduction
      const VarLabel* pTemperatureGradientLabel; //for heat conduction
      const VarLabel* pTemperatureRateLabel; //for heat conduction
      const VarLabel* pParticleIDLabel;
      const VarLabel* pIsIgnitedLabel; //for burn models
      const VarLabel* pMassRateLabel; //for burn models
      
      const VarLabel* gMassLabel;
      const VarLabel* gAccelerationLabel;
      const VarLabel* gMomExedAccelerationLabel;
      const VarLabel* gVelocityLabel;
      const VarLabel* gMomExedVelocityLabel;
      const VarLabel* gVelocityStarLabel;
      const VarLabel* gMomExedVelocityStarLabel;
      const VarLabel* gExternalForceLabel;
      const VarLabel* gInternalForceLabel;
      const VarLabel* gSelfContactLabel; //for fracture
      const VarLabel* gTemperatureRateLabel; //for heat conduction
      const VarLabel* gTemperatureLabel; //for heat conduction
      const VarLabel* gInternalHeatRateLabel;
      const VarLabel* gExternalHeatRateLabel;
      
      const VarLabel* cSelfContactLabel; //for fracture, CCVariable
      const VarLabel* cSurfaceNormalLabel; //for fracture, CCVariable
      
      const VarLabel* StrainEnergyLabel;
      const VarLabel* KineticEnergyLabel;

      const VarLabel* ppNAPIDLabel;
      
    };
  } // end namepsace MPM
} // end namespace Uintah


// $Log$
// Revision 1.8  2000/06/08 16:56:51  guilkey
// Added tasks and VarLabels for HE burn model stuff.
//
// Revision 1.7  2000/06/06 03:17:42  tan
// Added particle variable lable pAverageMicrocrackLength for fracture simulation.
//
// Revision 1.6  2000/06/02 23:16:32  guilkey
// Added ParticleID labels.
//
// Revision 1.5  2000/05/31 22:15:38  guilkey
// Added VarLabels for some integrated quantities.
//
// Revision 1.4  2000/05/31 16:11:11  tan
// gTemperatureLabel included
//
// Revision 1.3  2000/05/30 17:07:34  dav
// Removed commented out labels.  Other MPI fixes.  Changed delt to delT so I would stop thinking of it as just delta.
//
// Revision 1.2  2000/05/30 04:27:33  tan
// Added gTemperatureRateLabel for heat conduction computations.
//
// Revision 1.1  2000/05/26 21:37:30  jas
// Labels are now created and accessed using Singleton class MPMLabel.
//
#endif
