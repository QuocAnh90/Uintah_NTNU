//  HyperElasticDamage.cc 
//  
//  class ConstitutiveModel ConstitutiveModel data type -- 3D - 
//  holds ConstitutiveModel
//  information for the FLIP technique:
//    This constitutive model is based on theorem of Continuum Damage 
//    Mechanics (CDM) with following assumptions:
//	i.   Free energy equation is uncoupled as volumetric and deviatoric
//	     parts and the equation is the same as Neo-Hookean model (refer:
//	     CompNeoHook.h and CompNeoHook.cc).
//	ii.  The damage mechanism is associated with maximum distortional
//	     energy and is independent of hydrostatic pressure.
//	iii. The damage criterion is only function of equivalent strain. 
//	iv.  The damage process character function has exponential form.
//
//    Material property constants:
//	Young's modulus, Poisson's ratio;
//	Damage parameters: alpha[0, infinite), beta[0, 1].
//      Maximum equivalent strain: strainmax -- there will be no damage
//                                 when strain is less than this value.
//
//  Reference: "ON A FULLY THREE-DIMENSIONAL FINITE-STRAIN VISCOELASTIC
//		DAMAGE MODEL: FORMULATION AND COMPUTATIONAL ASPECTS", by
//		J.C.Simo, Computer Methods in Applied Mechanics and
//		Engineering 60 (1987) 153-173.
//
//		"COMPUTATIONAL INELASTICITY" by J.C.Simo & T.J.R.Hughes,
//		Springer, 1997.
//     
//
//    Features:
//     
//      
//      Usage:


#include <SCICore/Exceptions/InternalError.h>
#include "ConstitutiveModelFactory.h"
#include "HyperElasticDamage.h"
#include <SCICore/Malloc/Allocator.h>
#include <fstream>
#include <iostream>
using namespace std;
using namespace Uintah::MPM;
using SCICore::Exceptions::InternalError;

HyperElasticDamage::HyperElasticDamage(ProblemSpecP& ps)
{
  // Constructor
  // Initialize deformationGradient
  ps->require("bulk_modulus",d_Bulk);
  ps->require("shear_modulus",d_Shear);
  ps->require("alpha",d_Alpha);
  ps->require("beta",d_Beta);
  ps->require("max_equiv_strain",maxEquivStrain);

  deformationGradient.Identity();
  bElBar.Identity();
  damageG = 1.0;		// No damage at beginning
  E_bar.set(0.0);
 
}

HyperElasticDamage::HyperElasticDamage(double bulk, double shear,
			       double alpha, double beta, double strainmax): 
  d_Bulk(bulk),d_Shear(shear),d_Alpha(alpha),d_Beta(beta),maxEquivStrain(strainmax)
{
  // Main constructor
  // Initialize deformationGradient

  deformationGradient.Identity();
  bElBar.Identity();
  damageG = 1.0;		// No damage at beginning
  E_bar.set(0.0);

}

HyperElasticDamage::HyperElasticDamage(const HyperElasticDamage &cm):
  deformationGradient(cm.deformationGradient),
  bElBar(cm.bElBar),
  E_bar(cm.E_bar),
  current_E_bar(cm.current_E_bar),
  stressTensor(cm.stressTensor),
  d_Bulk(cm.d_Bulk),
  d_Shear(cm.d_Shear),
  d_Alpha(cm.d_Alpha),
  d_Beta(cm.d_Beta),
  damageG(cm.damageG),
  maxEquivStrain(cm.maxEquivStrain)
 
{
  // Copy constructor
 
}

HyperElasticDamage::~HyperElasticDamage()
{
  // Destructor
 
}

void HyperElasticDamage::setBulk(double bulk)
{
  // Assign HyperElasticDamage Bulk Modulus

  d_Bulk = bulk;
}

void HyperElasticDamage::setShear(double shear)
{
  // Assign HyperElasticDamage Shear Modulus

  d_Shear = shear;
}

void HyperElasticDamage::setDamageParameters(double alpha, double beta)
{
  // Assign HyperElasticDamage Damage Parameters

  d_Alpha = alpha;
  d_Beta = beta;
}

void HyperElasticDamage::setMaxEquivStrain(double strainmax)
{
  // Assign Maximun equivalent strain

  maxEquivStrain = strainmax;
}

void HyperElasticDamage::setStressTensor(Matrix3 st)
{
  // Assign the stress tensor (3 x 3 Matrix)

  stressTensor = st;

}

void HyperElasticDamage::setDeformationMeasure(Matrix3 dg) 
{
  // Assign the deformation gradient tensor (3 x 3 Matrix)

  deformationGradient = dg;

}

Matrix3 HyperElasticDamage::getStressTensor() const
{
  // Return the stress tensor (3 x 3 Matrix)

  return stressTensor;

}

Matrix3 HyperElasticDamage::getDeformationMeasure() const
{
  // Return the strain tensor (3 x 3 Matrix)

  return deformationGradient;

}

std::vector<double> HyperElasticDamage::getMechProps() const
{
  // Return bulk and shear modulus

  std::vector<double> props(6,0.0);

  props[0] = d_Bulk;
  props[1] = d_Shear;
  props[2] = d_Alpha;
  props[3] = d_Beta;
  props[4] = maxEquivStrain;
  props[5] = damageG;

  return props;

}


void HyperElasticDamage::addParticleState(std::vector<const VarLabel*>& from,
					std::vector<const VarLabel*>& to)
{
   throw InternalError("HyperElasticDamage won't work");
}

void HyperElasticDamage::computeStressTensor(const Patch* patch,
					     const MPMMaterial* matl,
					     DataWarehouseP& old_dw,
					     DataWarehouseP& new_dw)
{

#ifdef WONT_COMPILE_YET
  Matrix3 bElBarTrial, shearTrial,fbar,F_bar,C_bar;
  double J,p;
  double equiv_strain;
  Matrix3 damage_normal,deformationGradientInc;
  Matrix3 Ebar_increament;
  double onethird = (1.0/3.0);
  double damage_flag = 0.0;
  Matrix3 Identity;

  Identity.Identity();

  // Calculate the stress Tensor (symmetric 3 x 3 Matrix) given the
  // time step and the velocity gradient and the material constants

  
  // Compute the deformation gradient increment using the time_step
  // velocity gradient
  // F_n^np1 = dudx * dt + Identity
  deformationGradientInc = velocityGradient * time_step + Identity;

  // Update the deformation gradient tensor to its time n+1 value.
  deformationGradient = deformationGradientInc * deformationGradient;

  fbar = deformationGradientInc*pow(deformationGradientInc.Determinant(),
					-onethird);

  F_bar = deformationGradient*pow(deformationGradient.Determinant(),-onethird);

  C_bar = F_bar * F_bar.Transpose();

  // Deviatoric-Elastic part of the Lagrangian Strain Tensor
  current_E_bar = (C_bar - Identity)/2.0;

  // Calculate equivalent strain
  equiv_strain = sqrt(d_Shear*(C_bar.Trace() - 3.0));
  maxEquivStrain = (maxEquivStrain > equiv_strain) ? 
				maxEquivStrain : equiv_strain;

  // Normal of the damage surface
  damage_normal = Identity * d_Shear / (2.0*equiv_strain);

  // Check damage criterion
  if ( equiv_strain == maxEquivStrain && maxEquivStrain != 0.0 )   
  // on damage surface
  {
    Ebar_increament = current_E_bar - E_bar;

    for( int damage_i=1; damage_i<=3; damage_i++) 
       for( int damage_j=1; damage_j<=3; damage_j++)
     	  damage_flag += damage_normal(damage_i, damage_j) * 
			Ebar_increament(damage_i, damage_j);

  if (damage_flag > 0.0)	// loading: further damage 
     damageG = d_Beta + (1-d_Beta) * (1-exp(-maxEquivStrain/d_Alpha))
    			/ (maxEquivStrain/d_Alpha);

  } // end if
  E_bar = current_E_bar;

  bElBarTrial = fbar*bElBar*fbar.Transpose();

  // shearTrial is equal to the shear modulus times dev(bElBar)

  shearTrial = (bElBarTrial - Identity*onethird*bElBarTrial.Trace())
     		* d_Shear * damageG;

  // get the volumetric part of the deformation
  J = deformationGradient.Determinant();

  // get the hydrostatic part of the stress
  p = 0.5*d_Bulk*(J - 1.0/J);

  // compute the total stress (volumetric + deviatoric)

  stressTensor = Identity*J*p + shearTrial;

  bElBar = bElBarTrial;
#endif
}

//for fracture
void HyperElasticDamage::computeCrackSurfaceContactForce(const Patch* patch,
                                           const MPMMaterial* matl,
                                           DataWarehouseP& old_dw,
                                           DataWarehouseP& new_dw)
{
}

//for fracture
void HyperElasticDamage::addComputesAndRequiresForCrackSurfaceContact(
	                                     Task* task,
					     const MPMMaterial* matl,
					     const Patch* patch,
					     DataWarehouseP& old_dw,
					     DataWarehouseP& new_dw) const
{
}

double HyperElasticDamage::computeStrainEnergy(const Patch* patch,
					       const MPMMaterial* matl,
					       DataWarehouseP& new_dw)
{
#ifdef WONT_COMPILE_YET
  double strainenergy = 1;

  return strainenergy;
#endif
}

void HyperElasticDamage::initializeCMData(const Patch* patch,
					  const MPMMaterial* matl,
					  DataWarehouseP& new_dw)
{

}


void HyperElasticDamage::addComputesAndRequires(Task* task,
						const MPMMaterial* matl,
						const Patch* patch,
						DataWarehouseP& old_dw,
						DataWarehouseP& new_dw) const
{
   cerr << "HyperElasticDamage::addComputesAndRequires needs to be filled in\n";
}

double HyperElasticDamage::getLambda() const
{
  // Return the Lame constant lambda

  double lambda;

  lambda = d_Bulk - .6666666667*d_Shear;

  return lambda;

}

double HyperElasticDamage::getMu() const
{
  // Return the Lame constant mu

  return d_Shear;

}

void HyperElasticDamage::readParameters(ProblemSpecP ps, double *p_array)
{

  ps->require("bulk_modulus",p_array[0]);
  ps->require("shear_modulus",p_array[1]);
  ps->require("alpha",p_array[2]);
  ps->require("beta",p_array[3]);
  ps->require("max_equiv_strain",p_array[4]);

}

void HyperElasticDamage::writeParameters(ofstream& out, double *p_array)
{
  out << p_array[0] << " " << p_array[1] << " ";
  out << p_array[2] << " " << p_array[3] << " ";
  out << p_array[4] << " ";
}

ConstitutiveModel* HyperElasticDamage::readParametersAndCreate(ProblemSpecP ps)
{
  double p_array[5];
  readParameters(ps, p_array);
  return(create(p_array));
}

void HyperElasticDamage::writeRestartParameters(ofstream& out) const
{
  out << getType() << " ";
  out << d_Bulk << " " << d_Shear << " ";
  out << d_Alpha << " " << d_Beta << " ";
  out << maxEquivStrain << " ";
  out << (getDeformationMeasure())(1,1) << " "
      << (getDeformationMeasure())(1,2) << " "
      << (getDeformationMeasure())(1,3) << " "
      << (getDeformationMeasure())(2,1) << " "
      << (getDeformationMeasure())(2,2) << " "
      << (getDeformationMeasure())(2,3) << " "
      << (getDeformationMeasure())(3,1) << " "
      << (getDeformationMeasure())(3,2) << " "
      << (getDeformationMeasure())(3,3) << endl;
}

ConstitutiveModel* HyperElasticDamage::readRestartParametersAndCreate(ProblemSpecP ps)
{
#if 0
  Matrix3 dg(0.0);
  ConstitutiveModel *cm = readParametersAndCreate(in);
  
  in >> dg(1,1) >> dg(1,2) >> dg(1,3)
     >> dg(2,1) >> dg(2,2) >> dg(2,3)
     >> dg(3,1) >> dg(3,2) >> dg(3,3);
  cm->setDeformationMeasure(dg);
  
  return(cm);
#endif
}

ConstitutiveModel* HyperElasticDamage::create(double *p_array)
{
  return(scinew HyperElasticDamage(p_array[0],p_array[1],p_array[2],
				p_array[3],p_array[4]));
}

int HyperElasticDamage::getType() const
{
  //  return(ConstitutiveModelFactory::CM_HYPER_ELASTIC_DAMAGE);
}

std::string HyperElasticDamage::getName() const
{
  return("Hyper-Elastic-Damage");
}

int HyperElasticDamage::getNumParameters() const
{
  return(5);
}

void HyperElasticDamage::printParameterNames(ofstream& out) const
{
  out << "bulk" << endl
      << "shear" << endl
      << "d_Alpha" << endl
      << "d_Beta" << endl
      << "maxEquivStrain" << endl;
}

ConstitutiveModel* HyperElasticDamage::copy() const
{
  return( scinew HyperElasticDamage(*this) );
}

int HyperElasticDamage::getSize() const
{
  int s = 0;
  s += sizeof(double) * 5;  // properties
  s += sizeof(double) * 9;  // deformation gradient elements
  s += sizeof(double) * 6;  // stress tensor elements
  s += sizeof(double) * 9;  // bElBar elements
  s += sizeof(double) * 9;  // E_bar elements
  s += sizeof(double) * 9;  // current_E_bar elements
  s += sizeof(double) * 1;  // damG
  s += sizeof(int) * 1;     // type
  return(s);
}



//
// $Log$
// Revision 1.13  2000/09/12 16:52:10  tan
// Reorganized crack surface contact force algorithm.
//
// Revision 1.12  2000/06/16 05:03:05  sparker
// Moved timestep multiplier to simulation controller
// Fixed timestep min/max clamping so that it really works now
// Implemented "override" for reduction variables that will
//   allow the value of a reduction variable to be overridden
//
// Revision 1.11  2000/06/15 21:57:05  sparker
// Added multi-patch support (bugzilla #107)
// Changed interface to datawarehouse for particle data
// Particles now move from patch to patch
//
// Revision 1.10  2000/06/09 21:07:33  jas
// Added code to get the fudge factor directly into the constitutive model
// inititialization.
//
// Revision 1.9  2000/05/30 20:19:04  sparker
// Changed new to scinew to help track down memory leaks
// Changed region to patch
//
// Revision 1.8  2000/05/11 20:10:15  dav
// adding MPI stuff.  The biggest change is that old_dws cannot be const and so a large number of declarations had to change.
//
// Revision 1.7  2000/05/07 06:02:04  sparker
// Added beginnings of multiple patch support and real dependencies
//  for the scheduler
//
// Revision 1.6  2000/04/26 06:48:17  sparker
// Streamlined namespaces
//
// Revision 1.5  2000/04/25 18:42:34  jas
// Revised the factory method and constructor to take a ProblemSpec argument
// to create a new constitutive model.
//
// Revision 1.4  2000/04/19 21:15:56  jas
// Changed BoundedArray to vector<double>.  More stuff to compile.  Critical
// functions that need access to data warehouse still have WONT_COMPILE_YET
// around the methods.
//
// Revision 1.3  2000/04/14 17:34:42  jas
// Added ProblemSpecP capabilities.
//
// Revision 1.2  2000/03/20 17:17:09  sparker
// Made it compile.  There are now several #idef WONT_COMPILE_YET statements.
//
// Revision 1.1  2000/03/14 22:11:49  jas
// Initial creation of the constitutive model directory with the legacy
// constitutive model classes.
//
// Revision 1.1  2000/02/24 06:11:57  sparker
// Imported homebrew code
//
// Revision 1.1  2000/01/24 22:48:51  sparker
// Stuff may actually work someday...
//
// Revision 1.5  1999/12/17 22:05:24  guilkey
// Changed all constitutive models to take in velocityGradient and dt as
// arguments.  This allowed getting rid of velocityGradient as stored data
// in the constitutive model.  Also, in all hyperelastic models,
// deformationGradientInc was also removed from the private data.
//
// Revision 1.4  1999/11/17 22:26:36  guilkey
// Added guts to computeStrainEnergy functions for CompNeoHook CompNeoHookPlas
// and CompMooneyRivlin.  Also, made the computeStrainEnergy function non consted
// for all models.
//
// Revision 1.3  1999/11/17 20:08:47  guilkey
// Added a computeStrainEnergy function to each constitutive model
// so that we can have a valid strain energy calculation for functions
// other than the Elastic Model.  This is called from printParticleData.
// Currently, only the ElasticConstitutiveModel version gives the right
// answer, but that was true before as well.  The others will be filled in.
//
// Revision 1.2  1999/10/07 05:20:59  guilkey
// Fixed copy constructor, pack and unpack stream, and added logging.
//
//
