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

#include <CCA/Components/MPM/Materials/ConstitutiveModel/SoilModels/MatsuokaNakai.h>
#include <CCA/Components/MPM/Materials/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Math/MinMax.h>

#include <sci_defs/uintah_defs.h>

#include <iostream>
#include <string>

#include <unistd.h>
#include <Core/Exceptions/InvalidValue.h>

using namespace std; using namespace Uintah;

MatsuokaNakai::MatsuokaNakai(ProblemSpecP& ps,MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_NBASICINPUTS=5;
  d_NMGDC=0;

// Total number of properties
  d_NDMMPROP=d_NBASICINPUTS+d_NMGDC;

  // pre-initialize all of the user inputs to zero.
  for(int i = 0; i<d_NDMMPROP; i++){
     UI[i] = 0.;
  }
  // Read model parameters from the input file
  getInputParameters(ps);

 // Check that model parameters are valid and allow model to change if needed
 CheckModel(UI);

 //Create VarLabels for GeoModel internal state variables (ISVs)
 int nx;
 nx = d_NBASICINPUTS;

 for (int i = 0; i < nx; i++)
 {
     rinit[i] = UI[i];
 }
 d_NINSV = nx;

  d_NINSV=nx;
  //  cout << "d_NINSV = " << d_NINSV << endl;

  initializeLocalMPMLabels();
}

#if 0
MatsuokaNakai::MatsuokaNakai(const MatsuokaNakai* cm) : ConstitutiveModel(cm)
{
  for(int i=0;i<d_NDMMPROP;i++){
    UI[i] = cm->UI[i];
  }

  //Create VarLabels for MatsuokaNakai internal state variables (ISVs)
  initializeLocalMPMLabels();
}
#endif

MatsuokaNakai::~MatsuokaNakai()
{
   for (unsigned int i = 0; i< ISVLabels.size();i++){
     VarLabel::destroy(ISVLabels[i]);
   }
}

void MatsuokaNakai::CheckModel(double UI[])
{
    if (UI[0] <= 0.0) cerr << "Shear Modulus in the Matsuoka Nakai Model equal to: " << UI[0] << " This will cause the code to malfuncion. Any results obtained are invalid." << endl;
    if (UI[1] < 0.0) cerr << "Bulk Modulus in the Matsuoka Nakai Model is set to negative value: " << UI[1] << " This will cause the code to malfuncion. Any results obtained are invalid." << endl;
    if (UI[2] < 0.0) cerr << "Friction angle in the Matsuoka Nakai Model is set to: " << UI[2] << " This will cause the code to malfuncion. Any results obtained are invalid." << endl;
    if (UI[4] < 0.0) cerr << "lamda_c in the Matsuoka Nakai Model is set to: " << UI[4] << " This will cause the code to malfuncion. Any results obtained are invalid." << endl;
}

void MatsuokaNakai::Calculate_Stress(int& ninsv, double& dt,
    double UI[], double stress[], double D[], double svarg[], double& USM)

    /*
    copied from DIAMM, giving the input required

    int &ninsv, double &dt,double UI[], double stress[], double D[],
    double svarg[], double &USM
    C
    C     input arguments
    C     ===============
    C      NINSV      int                   Number of internal state vars
    C      DTARG      dp                    Current time increment
    C      UI       dp,ar(nprop)            User inputs
    C      D          dp,ar(6)              Strain rate
    C
    C     input output arguments
    C     ======================
    C      STRESS   dp,ar(6)                stress
    C      SVARG    dp,ar(ninsv)            state variables
    C
    C     output arguments
    C     ================
    C      USM      dp                      uniaxial strain modulus
    C***********************************************************************
    C
    C      stresss and strains, plastic strain tensors
    C          11, 22, 33, 12, 23, 13
    C
    C***********************************************************************
    */
{
    double deps[6]; // Strain increment (it needs to change the sign as the MN code read compression as possitive)
    for (int i = 0; i < 6; i++) {
        deps[i] = -D[i] * dt;
    }

    double K = UI[1];                                       // Bulk modulus
    double G = UI[0];                                       // Shear modulus
    USM = 1 * (G + 0.3 * K) / 3.0;                          // uniaxial strain modulus
    double Phi = svarg[2];                                  // Friction angle
    double M = 2 * sqrt(2) * tan(Phi * 3.1415 / 180);       // Friction state variable
    double N_ini = svarg[3];                                // Dilatancy
    double N = N_ini;
    double lamda_c = UI[4];                                 // Dilatancy parameter

    //double Plastic_Multiplier = 0;                          // Plastic multiplier
    double deps_plastic[6];                                 // plastic strain
    double deps_elastic[6];                                 // elastic strain

    for (int i = 0; i < 6; i++) {
        deps_plastic[i] = 0;
    }

    double K43G = K + 4.0 * G / 3.0;
    double K23G = K - 2.0 * G / 3.0;

    double ds_trial[6]; // Stress increment 
    ds_trial[0] = K43G * deps[0] + K23G * (deps[1] + deps[2]);
    ds_trial[1] = K43G * deps[1] + K23G * (deps[0] + deps[2]);
    ds_trial[2] = K43G * deps[2] + K23G * (deps[0] + deps[1]);
    ds_trial[3] = G * deps[3];
    ds_trial[4] = G * deps[4];
    ds_trial[5] = G * deps[5];

    // Trial stress (it needs to change the sign of previous stress as the MN code read compression as possitive)
    double stress_trial[6];
    //[0] xx // [1] yy // [2] zz // [3] xy // [4] yz // [5] xz
    for (int i = 0; i < 6; i++) {
        stress_trial[i] = -stress[i] + ds_trial[i];
    }

    // Definition
    //inline Matrix3(double v00, double v01, double v02,
    //  double v10, double v11, double v12,
    //  double v20, double v21, double v22);

    Matrix3 Stress_Tensor(stress_trial[0], stress_trial[3], stress_trial[5],
        stress_trial[3], stress_trial[1], stress_trial[4],
        stress_trial[5], stress_trial[4], stress_trial[2]);

    Matrix3 Stress_Tensor_power_2 = Stress_Tensor * Stress_Tensor;

   double I1 = Stress_Tensor.Trace();                                       // First stress invariants
   double I2 = 0.5 * (I1 * I1 - Stress_Tensor_power_2.Trace());             // Second stress invariants     
   double I3 = Stress_Tensor.Determinant();                                 // Third stress invariants
   
   double X = sqrt(I1 * I2 / I3 - 9);

   double f_trial = X - (M + N);      // Yield function

   // Cut off condition
   if (I3 <= 0.000 || I1 * I2 / I3 <= 9.000)
   {
       f_trial = 0;
       for (int i = 0; i < 6; i++) {
           stress_trial[i] = 0;
       }

       //cerr << "Cut off condition for I3 " << I3 << " and I1 * I2 / I3 <= 9 " << I1 * I2 / I3 << endl;
   }

   if (f_trial > 0)  // Elasto-plastic
   {
       // interate
       int count = 0;      
       double stress_inter[6];
       double f = 0;
       do {
           ++count;        
          
           // Stress Tensor Inverse
           Matrix3 Stress_Tensor_Inverse = Stress_Tensor.Inverse();
           Matrix3 Identity; Identity.Identity();
           double Stress_Tensor_Inverse_Trace = Stress_Tensor_Inverse.Trace();
           Matrix3 Stress_Tensor_Inverse_Dev = Stress_Tensor_Inverse - Stress_Tensor_Inverse_Trace / 3 * Identity;
           double Stress_Tensor_Inverse_Dev_Vector[6];
           Stress_Tensor_Inverse_Dev_Vector[0] = Stress_Tensor_Inverse_Dev(0, 0);
           Stress_Tensor_Inverse_Dev_Vector[1] = Stress_Tensor_Inverse_Dev(1, 1);
           Stress_Tensor_Inverse_Dev_Vector[2] = Stress_Tensor_Inverse_Dev(2, 2);
           Stress_Tensor_Inverse_Dev_Vector[3] = Stress_Tensor_Inverse_Dev(0, 1);
           Stress_Tensor_Inverse_Dev_Vector[4] = Stress_Tensor_Inverse_Dev(1, 2);
           Stress_Tensor_Inverse_Dev_Vector[5] = Stress_Tensor_Inverse_Dev(2, 0);

           // Flow rule
           double dQdSigma[6];
           for (int i = 0; i < 3; i++) { // For xx yy zz
               dQdSigma[i] = -N / 3 - I1 / X * Stress_Tensor_Inverse_Dev_Vector[i];
           }
           for (int i = 3; i < 6; i++) { // For xy yz xz
               dQdSigma[i] = -I1 / X * Stress_Tensor_Inverse_Dev_Vector[i];
           }

           double A = lamda_c * (M + N) * N;

           // Derivative
           double dYdI1 = I2 / (I3 * 2 * X);
           double dYdI2 = I1 / (I3 * 2 * X);
           double dYdI3 = -I1 * I2 / (I3 * I3 * 2 * X);

           Matrix3 dI1dSigma = Identity;

           Matrix3 S = Stress_Tensor; // Reduce Notation
           Matrix3 dI2dSigma(S(1, 1) + S(2, 2)      , -2 * S(0, 1)              , -2 * S(0, 2),
               -2 * S(0, 1)                         , S(0, 0) + S(2, 2)         , -2 * S(1, 2),
               -2 * S(0, 2)                         , -2 * S(1, 2)              , S(1, 1) + S(0, 0));

           Matrix3 dI3dSigma (S(1, 1) * S(2, 2) - S(1, 2) * S(1, 2),      2 * (S(1, 2) * S(0, 2) - S(2, 2) * S(0, 1))         , 2 * (S(1, 2) * S(0, 1) - S(1, 1) * S(0, 2)),
               2 * (S(1, 2) * S(0, 2) - S(2, 2) * S(0, 1)),                 S(0, 0) * S(2, 2) - S(0, 2) * S(0, 2)               , 2 * (S(0, 2) * S(0, 1) - S(0, 0) * S(1, 2)),
               2 * (S(1, 2) * S(0, 1) - S(1, 1) * S(0, 2)),                 2 * (S(0, 2) * S(0, 1) - S(0, 0) * S(1, 2))         , S(1, 1) * S(0, 0) - S(0, 2) * S(0, 2));

           if (I1 > 5000) {
               //cerr << count << endl;

               //cerr << " dI3dSigma(0,0) " << dI3dSigma(0, 0) << endl;

               //cerr << " dI3dSigma " << dI3dSigma << endl;

               //cerr << " S(1, 1) * S(2, 2) - S(1, 2) * S(1, 2) " << dI3dSigma_xx << endl;
               //cerr << " S(0, 0) * S(2, 2) - S(0, 2) * S(0, 2)  " << dI3dSigma_yy << endl;
               //cerr << " S(1, 1) * S(0, 0) - S(0, 2) * S(0, 2) " << dI3dSigma_zz << endl;
               //cerr << " 2 * (S(1, 2) * S(0, 2) - S(2, 2) * S(0, 1)) " << dI3dSigma_xy << endl;
              // cerr << " 2 * (S(0, 2) * S(0, 1) - S(0, 0) * S(1, 2)) " << dI3dSigma_yz << endl;
               //cerr << " 2 * (S(1, 2) * S(0, 1) - S(1, 1) * S(0, 2)) " << dI3dSigma_xz << endl;
           }

           Matrix3 dYdSigma = dYdI1 * dI1dSigma + dYdI2 * dI2dSigma + dYdI3 * dI3dSigma;

           if (I1 > 5000) {
               cerr << count << endl;

               cerr << " dYdI1 * dI1dSigma  " << dYdI1 * dI1dSigma << endl;
               cerr << "  dYdI2 * dI2dSigma " << dYdI2 * dI2dSigma << endl;
               cerr << " dYdI3 * dI3dSigma " << dYdI3 * dI3dSigma << endl;
               cerr << " dYdSigma " << dYdSigma << endl;
           }

           double dYdSigma_Vector[6];
           dYdSigma_Vector[0] = dYdSigma(0, 0);
           dYdSigma_Vector[1] = dYdSigma(1, 1);
           dYdSigma_Vector[2] = dYdSigma(2, 2);
           dYdSigma_Vector[3] = dYdSigma(0, 1);
           dYdSigma_Vector[4] = dYdSigma(1, 2);
           dYdSigma_Vector[5] = dYdSigma(2, 0);

           // Update plastic multiplier 
           double dYdSigma_x_ElasticMatrix[6];
           dYdSigma_x_ElasticMatrix[0] = K43G * dYdSigma_Vector[0] + K23G * (dYdSigma_Vector[1] + dYdSigma_Vector[2]);
           dYdSigma_x_ElasticMatrix[1] = K43G * dYdSigma_Vector[1] + K23G * (dYdSigma_Vector[0] + dYdSigma_Vector[2]);
           dYdSigma_x_ElasticMatrix[2] = K43G * dYdSigma_Vector[2] + K23G * (dYdSigma_Vector[0] + dYdSigma_Vector[1]);
           dYdSigma_x_ElasticMatrix[3] = G * dYdSigma_Vector[3];
           dYdSigma_x_ElasticMatrix[4] = G * dYdSigma_Vector[4];
           dYdSigma_x_ElasticMatrix[5] = G * dYdSigma_Vector[5];

           double B = 0; // dYdSigma_x_ElasticMatrix_x_dQdSigma 
           for (int i = 0; i < 6; i++) {
               B = B + dYdSigma_x_ElasticMatrix[i] * dQdSigma[i];
           }

           double Plastic_Multiplier_increment = f / (A + B);

           if (I1 > 5000) {

               for (int i = 0; i < 6; i++) {
                   cerr << " dQdSigma[i] " << dQdSigma[i] << endl;
                   cerr << " deps[i] " << deps[i] << endl;
                   cerr << " stress[i] " << -stress[i] << endl;
               }

               cerr << " X " << X << endl;
               cerr << " I1 " << I1 << endl;
               cerr << " N " << N << endl;
               cerr << " Plastic_Multiplier_increment " << Plastic_Multiplier_increment << endl;
               cerr << " B " << B << endl;

               throw InvalidValue("Stop to check", __FILE__, __LINE__);
           }

           // Update plastic strain
           for (int i = 0; i < 6; i++) {
               deps_plastic[i] = deps_plastic[i] + Plastic_Multiplier_increment * dQdSigma[i];
               deps_elastic[i] = deps[i] - deps_plastic[i];
           }

           // Update stress
           double ds[6]; // Stress increment 
           ds[0] = K43G * deps_elastic[0] + K23G * (deps_elastic[1] + deps_elastic[2]);
           ds[1] = K43G * deps_elastic[1] + K23G * (deps_elastic[0] + deps_elastic[2]);
           ds[2] = K43G * deps_elastic[2] + K23G * (deps_elastic[0] + deps_elastic[1]);
           ds[3] = G * deps_elastic[3];
           ds[4] = G * deps_elastic[4];
           ds[5] = G * deps_elastic[5];
          
           //[0] xx // [1] yy // [2] zz // [3] xy // [4] yz // [5] xz
           for (int i = 0; i < 6; i++) {
               stress_inter[i] = -stress[i] + ds[i];
           }

           // Update state variables
           double deps_volume = deps_plastic[0] + deps_plastic[1] + deps_plastic[2];
           double N = N_ini + lamda_c * (M + N) * deps_volume;

           // Recompute yield function
           Matrix3 Stress_Tensor(   stress_inter[0], stress_inter[3], stress_inter[5],
                                    stress_inter[3], stress_inter[1], stress_inter[4],
                                    stress_inter[5], stress_inter[4], stress_inter[2]);

           Stress_Tensor_power_2 = Stress_Tensor * Stress_Tensor;

           I1 = Stress_Tensor.Trace();                                       // First stress invariants
           I2 = 0.5 * (I1 * I1 - Stress_Tensor_power_2.Trace());             // Second stress invariants     
           I3 = Stress_Tensor.Determinant();                                 // Third stress invariants

           X = sqrt(I1 * I2 / I3 - 9);

           f = X - (M + N);      // Yield function


           if (count > 100) {
               //cerr << " Plastic_Multiplier_increment " << Plastic_Multiplier_increment << endl;
               for (int i = 0; i < 4; i++) {

                   cerr << " deps_plastic[i] " << deps_plastic[i] << endl;
                   cerr << " deps_elastic[i] " << deps_elastic[i] << endl;
                   cerr << " stress_inter[i] " << stress_inter[i] << endl;
                   cerr << " stress_trial[i] " << stress_trial[i] << endl;
                   
               }
               cerr << " f_trial " << f_trial << endl;

               cerr << " f " << f << endl;
               cerr << " X " << X << endl;
               cerr << " N " << N << endl;
               cerr << " Plastic_Multiplier_increment " << Plastic_Multiplier_increment << endl;
               cerr << " B " << B << endl;

               for (int i = 0; i < 6; i++) {
                   cerr << " dYdSigma_x_ElasticMatrix[i] " << dYdSigma_x_ElasticMatrix[i] << endl;
               }
           }

           if (count > 110) {             
               throw InvalidValue("More than 100 interations", __FILE__, __LINE__);
           }

           // Cut off condition
           if (I3 <= 0.000 || I1 * I2 / I3 <= 9.000)
           {
               f = 0;
               for (int i = 0; i < 6; i++) {
                   stress_inter[i] = 0;
               }

               //cerr << "Cut off condition for I3 " << I3 << " and I1 * I2 / I3 <= 9 " << I1 * I2 / I3 << endl;
           }

           // Liquefraction cut-ogg
           if (I1 <= 1000.000)
           {
               f = 0;
               for (int i = 0; i < 6; i++) {
                   stress_inter[i] = 0;
               }

               //cerr << "Cut off condition for liquefraction " << I1 << endl;
           }

       } while (abs(f)>0.001);

       // Update stress and state variables
       svarg[2] = Phi;  // Friction
       svarg[3] = N;    // Dilatancy

       // Update stress (it needs to change the sign of previous stress as the MN code read compression as possitive)
       for (int i = 0; i < 6; i++) {
           stress[i] = -stress_inter[i];
       }
   }
   else {
   // Update stress (it needs to change the sign of previous stress as the MN code read compression as possitive)
   for (int i = 0; i < 6; i++) {
       stress[i] = -stress_trial[i];
   }
    }


}


void MatsuokaNakai::outputProblemSpec(ProblemSpecP& ps,bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type","MatsuokaNakai");
  }

  cm_ps->appendElement("G",UI[0]);          // initial shear modulus 
  cm_ps->appendElement("K",UI[1]);          // initial bulk modulus 
  cm_ps->appendElement("Phi",UI[2]);        // initial friction angle 

  cm_ps->appendElement("N",UI[3]);          // initial dilatancy state variable 
  cm_ps->appendElement("lamda_c",UI[4]);    // parameter for dilatancy
}

MatsuokaNakai* MatsuokaNakai::clone()
{
  return scinew MatsuokaNakai(*this);
}

void MatsuokaNakai::initializeCMData(const Patch* patch,
                               const MPMMaterial* matl,
                               DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  std::vector<ParticleVariable<double> > ISVs(d_NINSV+1);

  cout << "In initializeCMData" << endl;
  for(int i=0;i<d_NINSV;i++){
    new_dw->allocateAndPut(ISVs[i],ISVLabels[i], pset);
    ParticleSubset::iterator iter = pset->begin();
    for(;iter != pset->end(); iter++){
      ISVs[i][*iter] = rinit[i];
    }
  }

  computeStableTimeStep(patch, matl, new_dw);
}

void MatsuokaNakai::addParticleState(std::vector<const VarLabel*>& from,
                               std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  for(int i=0;i<d_NINSV;i++){
    from.push_back(ISVLabels[i]);
    to.push_back(ISVLabels_preReloc[i]);
  }
}

void MatsuokaNakai::computeStableTimeStep(const Patch* patch,
                                    const MPMMaterial* matl,
                                    DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pmass, pvolume;
  constParticleVariable<Vector> pvelocity;

  new_dw->get(pmass,     lb->pMassLabel,     pset);
  new_dw->get(pvolume,   lb->pVolumeLabel,   pset);
  new_dw->get(pvelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12,1.e-12,1.e-12);

  double bulk = UI[0];
  double G = UI[3];
  for(ParticleSubset::iterator iter = pset->begin();iter != pset->end();iter++){
     particleIndex idx = *iter;

     // Compute wave speed at each particle, store the maximum
     c_dil = sqrt((bulk + 4.*G/3.)*pvolume[idx]/pmass[idx]);
     WaveSpeed=Vector(Max(c_dil+fabs(pvelocity[idx].x()),WaveSpeed.x()),
                      Max(c_dil+fabs(pvelocity[idx].y()),WaveSpeed.y()),
                      Max(c_dil+fabs(pvelocity[idx].z()),WaveSpeed.z()));
  }
  //UI[14]=matl->getInitialDensity();
  //UI[15]=matl->getRoomTemperature();
  //UI[14]=bulk/matl->getInitialDensity();  ??tim
  //UI[19]=matl->getInitialCv();
  WaveSpeed = dx/WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void MatsuokaNakai::computeStressTensor(const PatchSubset* patches,
                                  const MPMMaterial* matl,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  double rho_orig = matl->getInitialDensity();
  for(int p=0;p<patches->size();p++){
    double se = 0.0;
    const Patch* patch = patches->get(p);

    Matrix3 Identity; Identity.Identity();
    double c_dil=0.0;
    Vector WaveSpeed(1.e-12,1.e-12,1.e-12);
    Vector dx = patch->dCell();

    int dwi = matl->getDWIndex();
    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    constParticleVariable<Matrix3> deformationGradient, pstress;
    ParticleVariable<Matrix3> pstress_new;
    constParticleVariable<Matrix3> deformationGradient_new, velGrad;
    constParticleVariable<double> pmass, pvolume, ptemperature;
    constParticleVariable<double> pvolume_new;
    constParticleVariable<Vector> pvelocity;
    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    old_dw->get(pstress,             lb->pStressLabel,             pset);
    old_dw->get(pmass,               lb->pMassLabel,               pset);
    old_dw->get(pvolume,             lb->pVolumeLabel,             pset);
    old_dw->get(pvelocity,           lb->pVelocityLabel,           pset);
    old_dw->get(ptemperature,        lb->pTemperatureLabel,        pset);
    old_dw->get(deformationGradient, lb->pDeformationMeasureLabel, pset);

    std::vector<constParticleVariable<double> > ISVs(d_NINSV+1);
    for(int i=0;i<d_NINSV;i++){
      old_dw->get(ISVs[i],           ISVLabels[i],                 pset);
    }

    ParticleVariable<double> pdTdt,p_q;

    new_dw->allocateAndPut(pstress_new,     lb->pStressLabel_preReloc,   pset);
    new_dw->allocateAndPut(pdTdt,           lb->pdTdtLabel,              pset);
    new_dw->allocateAndPut(p_q,             lb->p_qLabel_preReloc,       pset);
    new_dw->get(deformationGradient_new,
                                 lb->pDeformationMeasureLabel_preReloc,  pset);
    new_dw->get(pvolume_new,     lb->pVolumeLabel_preReloc,              pset);
    new_dw->get(velGrad,         lb->pVelGradLabel_preReloc,             pset);

    std::vector<ParticleVariable<double> > ISVs_new(d_NINSV+1);
    for(int i=0;i<d_NINSV;i++){
      new_dw->allocateAndPut(ISVs_new[i],ISVLabels_preReloc[i], pset);
    }

    for(ParticleSubset::iterator iter = pset->begin();
                                        iter != pset->end(); iter++){
      particleIndex idx = *iter;

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Calculate rate of deformation D, and deviatoric rate DPrime,
      Matrix3 D = (velGrad[idx] + velGrad[idx].Transpose())*.5;

      // get the volumetric part of the deformation
      double J = deformationGradient_new[idx].Determinant();
      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        cerr << getpid() ;
        constParticleVariable<long64> pParticleID;
        old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
        cerr << "**ERROR** Negative Jacobian of deformation gradient"
             << " in particle " << pParticleID[idx] << endl;
        cerr << "l = " << velGrad[idx] << endl;
        cerr << "F_old = " << deformationGradient[idx] << endl;
        cerr << "F_new = " << deformationGradient_new[idx] << endl;
        cerr << "J = " << J << endl;
        throw InternalError("Negative Jacobian",__FILE__,__LINE__);
      }

      // Compute the local sound speed
      double rho_cur = rho_orig/J;

      // NEED TO FIND R
      Matrix3 tensorR, tensorU;

      // Look into using Rebecca's PD algorithm
      deformationGradient_new[idx].polarDecompositionRMB(tensorU, tensorR);

      // This is the previous timestep Cauchy stress
      // unrotated tensorSig=R^T*pstress*R
      Matrix3 tensorSig = (tensorR.Transpose())*(pstress[idx]*tensorR);

      // Load into 1-D array for the fortran code
      double sigarg[6];
      sigarg[0]=tensorSig(0,0);
      sigarg[1]=tensorSig(1,1);
      sigarg[2]=tensorSig(2,2);
      sigarg[3]=tensorSig(0,1);
      sigarg[4]=tensorSig(1,2);
      sigarg[5]=tensorSig(2,0);

      // UNROTATE D: S=R^T*D*R
      D=(tensorR.Transpose())*(D*tensorR);

      // Load into 1-D array for the fortran code
      double Darray[6];
      Darray[0]=D(0,0);
      Darray[1]=D(1,1);
      Darray[2]=D(2,2);
      Darray[3]=D(0,1);
      Darray[4]=D(1,2);
      Darray[5]=D(2,0);
      double svarg[d_NINSV];
      double USM=9e99;
      double dt = delT;

      // Load ISVs into a 1D array for fortran code
      for(int i=0;i<d_NINSV;i++){
        svarg[i]=ISVs[i][idx];
      }

      Calculate_Stress(d_NINSV, dt, UI, sigarg, Darray, svarg, USM);

      // Unload ISVs from 1D array into ISVs_new
      for(int i=0;i<d_NINSV;i++){
        ISVs_new[i][idx]=svarg[i];
      }

      // This is the Cauchy stress, still unrotated
      tensorSig(0,0) = sigarg[0];
      tensorSig(1,1) = sigarg[1];
      tensorSig(2,2) = sigarg[2];
      tensorSig(0,1) = sigarg[3];
      tensorSig(1,0) = sigarg[3];
      tensorSig(2,1) = sigarg[4];
      tensorSig(1,2) = sigarg[4];
      tensorSig(2,0) = sigarg[5];
      tensorSig(0,2) = sigarg[5];

      // ROTATE pstress_new: S=R*tensorSig*R^T
      pstress_new[idx] = (tensorR*tensorSig)*(tensorR.Transpose());

      c_dil = sqrt(USM/rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 AvgStress = (pstress_new[idx] + pstress[idx])*.5;

      double e = (D(0,0)*AvgStress(0,0) +
                  D(1,1)*AvgStress(1,1) +
                  D(2,2)*AvgStress(2,2) +
              2.*(D(0,1)*AvgStress(0,1) +
                  D(0,2)*AvgStress(0,2) +
                  D(1,2)*AvgStress(1,2))) * pvolume_new[idx]*delT;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pvelocity_idx = pvelocity[idx];
      WaveSpeed=Vector(Max(c_dil+fabs(pvelocity_idx.x()),WaveSpeed.x()),
                       Max(c_dil+fabs(pvelocity_idx.y()),WaveSpeed.y()),
                       Max(c_dil+fabs(pvelocity_idx.z()),WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z())/3.0;
        double c_bulk = sqrt(UI[0]/rho_cur);
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    }  // end loop over particles

    WaveSpeed = dx/WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se),     lb->StrainEnergyLabel);
    }
  }
}

void MatsuokaNakai::carryForward(const PatchSubset* patches,
                           const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    int dwi = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    std::vector<constParticleVariable<double> > ISVs(d_NINSV+1);
    std::vector<ParticleVariable<double> > ISVs_new(d_NINSV+1);

    for(int i=0;i<d_NINSV;i++){
      old_dw->get(ISVs[i],ISVLabels[i], pset);
      new_dw->allocateAndPut(ISVs_new[i],ISVLabels_preReloc[i], pset);
      ISVs_new[i].copyData(ISVs[i]);
  }

    // Don't affect the strain energy or timestep size
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());
    
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.),     lb->StrainEnergyLabel);
    }
  }

}

void MatsuokaNakai::addInitialComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet* ) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  cout << "In add InitialComputesAnd" << endl;

  // Other constitutive model and input dependent computes and requires
  for(int i=0;i<d_NINSV;i++){
    task->computes(ISVLabels[i], matlset);
  }
}

void MatsuokaNakai::addComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  // Computes and requires for internal state data
  for(int i=0;i<d_NINSV;i++){
    task->requires(Task::OldDW, ISVLabels[i],          matlset, Ghost::None);
    task->computes(             ISVLabels_preReloc[i], matlset);
  }
}

void MatsuokaNakai::addComputesAndRequires(Task*,
                                     const MPMMaterial*,
                                     const PatchSet*,
                                     const bool ) const
{
}

double MatsuokaNakai::computeRhoMicroCM(double pressure,
                                  const double p_ref,
                                  const MPMMaterial* matl, 
                                  double temperature,
                                  double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double p_gauge = pressure - p_ref;
  double rho_cur;
  double bulk = UI[0];

  rho_cur = rho_orig/(1-p_gauge/bulk);

  return rho_cur;

#if 1
  cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR MatsuokaNakai" << endl;
#endif
}

void MatsuokaNakai::computePressEOSCM(double rho_cur, double& pressure,
                                double p_ref,
                                double& dp_drho,      double& tmp,
                                const MPMMaterial* matl, 
                                double temperature)
{

  double bulk = UI[0];
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk*(1.0 - rho_orig/rho_cur);
  pressure = p_ref + p_g;
  dp_drho  = bulk*rho_orig/(rho_cur*rho_cur);
  tmp = bulk/rho_cur;  // speed of sound squared

#if 1
  cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR MatsuokaNakai" << endl;
#endif
}

double MatsuokaNakai::getCompressibility()
{
  return 1.0/UI[0];
}

void
MatsuokaNakai::getInputParameters(ProblemSpecP& ps)
{
  ps->getWithDefault("G",UI[0],0.0);              // initial shear modulus 
  ps->getWithDefault("K",UI[1],0.0);              // initial bulk modulus
  ps->getWithDefault("Phi",UI[2],0.0);            // initial friction angle modulus 
  ps->getWithDefault("N",UI[3],0.0);              // initial dilatancy state variable 
  ps->getWithDefault("lamda_c",UI[4],0.0);        // parameter for dilatancy
}

void
MatsuokaNakai::initializeLocalMPMLabels()
{
  vector<string> ISVNames;

  ISVNames.push_back("G");
  ISVNames.push_back("K");
  ISVNames.push_back("Phi");
  ISVNames.push_back("N");
  ISVNames.push_back("lamda_c");

  for(int i=0;i<d_NINSV;i++){
    ISVLabels.push_back(VarLabel::create(ISVNames[i],
                          ParticleVariable<double>::getTypeDescription()));
    ISVLabels_preReloc.push_back(VarLabel::create(ISVNames[i]+"+",
                          ParticleVariable<double>::getTypeDescription()));
  }
}
