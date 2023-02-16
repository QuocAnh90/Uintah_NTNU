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

//  MatsuokaNakai.h
//  class ConstitutiveModel ConstitutiveModel data type -- 3D -
//  holds ConstitutiveModel
//    This is for calling the Diamm model
//    Features:
//      Usage:



#ifndef __MATSUOKA_NAKAI_H__
#define __MATSUOKA_NAKAI_H__

#include <cmath>
#include "CCA/Components/MPM/Materials/ConstitutiveModel/ConstitutiveModel.h"
#include <Core/Math/Matrix3.h>
#include <vector>
#include <Core/Grid/Variables/VarLabel.h>

namespace Uintah {
  class MatsuokaNakai : public ConstitutiveModel {
  public:

    int d_NDMMPROP,d_NMGDC;
    int d_NBASICINPUTS;
    double rinit[100];
    double UI[190];

    std::vector<const VarLabel*> ISVLabels;
    std::vector<const VarLabel*> ISVLabels_preReloc;
    int d_NINSV;

  private:
    // Prevent copying of this class
    // copy constructor
      MatsuokaNakai& operator=(const MatsuokaNakai&cm);

    void getInputParameters(ProblemSpecP& ps);

    void initializeLocalMPMLabels();

    void CheckModel(double UI[]);

    void Calculate_Stress(int& nblk, int& ninsv, double& dt,
        double UI[], double stress[], double D[], double svarg[], double& USM);


  public:
    // constructors
      MatsuokaNakai(ProblemSpecP& ps, MPMFlags* flag);
    //    Diamm(const Diamm* cm);

    // destructor
    virtual ~MatsuokaNakai();

    virtual void outputProblemSpec(ProblemSpecP& ps,bool output_cm_tag = true);

    // clone
    MatsuokaNakai* clone();

    // compute stable timestep for this patch
    virtual void computeStableTimeStep(const Patch* patch,
                                       const MPMMaterial* matl,
                                       DataWarehouse* new_dw);

    // compute stress at each particle in the patch
    virtual void computeStressTensor(const PatchSubset* patches,
                                     const MPMMaterial* matl,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw);


    // carry forward CM data for RigidMPM
    virtual void carryForward(const PatchSubset* patches,
                              const MPMMaterial* matl,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

    // initialize  each particle's constitutive model data
    virtual void initializeCMData(const Patch* patch,
                                  const MPMMaterial* matl,
                                  DataWarehouse* new_dw);

    virtual void addInitialComputesAndRequires(Task* task,
                                               const MPMMaterial* matl,
                                               const PatchSet* patches) const;

    virtual void addComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet* patches) const;

    virtual void addComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet* patches,
                                        const bool recursion) const;

    virtual double computeRhoMicroCM(double pressure,
                                     const double p_ref,
                                     const MPMMaterial* matl, 
                                     double temperature,
                                     double rho_guess);


    virtual void computePressEOSCM(double rho_m, double& press_eos,
                                   double p_ref,
                                   double& dp_drho, double& ss_new,
                                   const MPMMaterial* matl, 
                                   double temperature);

    virtual double getCompressibility();


    virtual void addParticleState(std::vector<const VarLabel*>& from,
                                  std::vector<const VarLabel*>& to);


  };

} // End namespace Uintah

#endif  // __MATSUOKA_NAKAI_H__
