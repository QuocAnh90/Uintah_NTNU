#include "ParticlesNeighbor.h"

#include "CellsNeighbor.h"
#include "Lattice.h"
#include "Cell.h"
#include "LeastSquare.h"
#include <SCICore/Exceptions/InternalError.h>

#include <Uintah/Components/MPM/Util/Matrix3.h>

#include <Uintah/Grid/Patch.h>

#include <iostream>

namespace Uintah {
namespace MPM {

using SCICore::Exceptions::InternalError;

ParticlesNeighbor::ParticlesNeighbor()
: std::vector<particleIndex>()
{
}

void ParticlesNeighbor::buildIn(const IntVector& cellIndex,const Lattice& lattice)
{
  CellsNeighbor cellsNeighbor;
  cellsNeighbor.buildIncluding(cellIndex,lattice);
  
  for(CellsNeighbor::const_iterator iter_cell = cellsNeighbor.begin();
    iter_cell != cellsNeighbor.end();
    ++iter_cell )
  {
    std::vector<particleIndex>& parts = lattice[*iter_cell].particles;
    for( std::vector<particleIndex>::const_iterator iter_p = parts.begin();
         iter_p != parts.end();
         ++iter_p )
    {
      push_back(*iter_p);
    }
  }
}

/*
void  ParticlesNeighbor::interpolateVector(LeastSquare& ls,
                          const particleIndex& pIdx,
                          const ParticleVariable<Vector>& pVector,
                          Vector& data,
                          Matrix3& gradient) const
{
  Vector v;
  for(int i=0;i<3;++i) {
    ls.clean();
    for(const_iterator pIter = begin(); pIter != end(); pIter++) {
      ls.input( (*d_pX)[*pIter]-(*d_pX)[pIdx], pVector[*pIter](i) );
    }
    ls.output( data(i),v );
    for(int j=0;j<3;++j) {
      gradient(i,j) = v(j);
    }
  }
}

void  ParticlesNeighbor::interpolatedouble(LeastSquare& ls,
                          const particleIndex& pIdx,
                          const ParticleVariable<double>& pdouble,
                          double& data,
                          Vector& gradient) const
{
  ls.clean();
  for(const_iterator pIter = begin(); pIter != end(); pIter++) {
    ls.input( (*d_pX)[*pIter]-(*d_pX)[pIdx], pdouble[*pIter] );
  }
  ls.output( data,gradient );
}

void  ParticlesNeighbor::interpolateInternalForce(LeastSquare& ls,
                          const particleIndex& pIdx,
                          const ParticleVariable<Matrix3>& pStress,
                          Vector& pInternalForce) const
{
  Vector v;
  double data;
  for(int i=0;i<3;++i) {
    pInternalForce(i) = 0;
    for(int j=0;j<3;++j) {
      ls.clean();
      for(const_iterator pIter = begin(); pIter != end(); pIter++) {
        ls.input( (*d_pX)[*pIter]-(*d_pX)[pIdx], pStress[*pIter](i,j) );
      }
    }
    ls.output( data,v );
    pInternalForce(i) -= v(i);
  }
}
*/

bool ParticlesNeighbor::visible(particleIndex idx,
                                const Point& B,
				const ParticleVariable<Point>& pX,
				const ParticleVariable<int>& pIsBroken,
				const ParticleVariable<Vector>& pCrackSurfaceNormal,
				const ParticleVariable<double>& pVolume) const
{
  const Point& A = pX[idx];
  if(pIsBroken[idx]) {
    if( Dot( B - A, pCrackSurfaceNormal[idx] ) > 0 ) return false;
  }
  
  for(int i=0; i<(int)size(); i++) {
      int index = (*this)[i];
      //cout<<"neighbor distance"<<(pX[index] - A).length()<<endl;
      
      if( index != idx && pIsBroken[index] ) {
        const Vector& N = pCrackSurfaceNormal[index];
        double size2 = pow(pVolume[index],0.666666667) /4 * 2;
        const Point& O = pX[index];

	double A_N = Dot(A,N);
	
        double a = A_N - Dot(O,N);
        double b = A_N - Dot(B,N);
	
	if(b != 0) {
	  double lambda = a/b;
	  if( lambda>=0 && lambda<=1 ) {
	    Point p( A.x() * (1-lambda) + B.x() * lambda,
                     A.y() * (1-lambda) + B.y() * lambda,
		     A.z() * (1-lambda) + B.z() * lambda );
 	    if( (p - O).length2() < size2 ) return false;
	  }
	}
      }
  }
  return true;
}

} //namespace MPM
} //namespace Uintah

// $Log$
// Revision 1.14  2000/09/25 20:23:21  sparker
// Quiet g++ warnings
//
// Revision 1.13  2000/09/22 07:18:57  tan
// MPM code works with fracture in three point bending.
//
// Revision 1.12  2000/09/16 04:18:04  tan
// Modifications to make fracture works well.
//
// Revision 1.11  2000/09/12 16:52:11  tan
// Reorganized crack surface contact force algorithm.
//
// Revision 1.10  2000/09/11 00:15:00  tan
// Added calculations on random distributed microcracks in broken particles.
//
// Revision 1.9  2000/09/08 18:25:44  tan
// Added visibility calculation to fracture broken cell shape function
// interpolation.
//
// Revision 1.8  2000/07/06 16:59:34  tan
// Least square interpolation added for particle velocities and stresses
// updating.
//
// Revision 1.7  2000/07/06 06:23:23  tan
// Added Least Square interpolation of double (such as temperatures),
// vector (such as velocities) and stresses for particles in the
// self-contact cells.
//
// Revision 1.6  2000/06/27 23:11:05  jas
// Added in grid bcs.
//
// Revision 1.5  2000/06/23 21:56:40  tan
// Use vector instead of list for cells-neighbor and particles-neighbor.
//
// Revision 1.4  2000/06/15 21:57:10  sparker
// Added multi-patch support (bugzilla #107)
// Changed interface to datawarehouse for particle data
// Particles now move from patch to patch
//
// Revision 1.3  2000/06/06 01:58:25  tan
// Finished functions build particles neighbor for a given particle
// index.
//
// Revision 1.2  2000/06/05 22:30:11  tan
// Added interpolateVector and interpolatedouble for least-square approximation.
//
// Revision 1.1  2000/06/05 21:15:36  tan
// Added class ParticlesNeighbor to handle neighbor particles searching.
//
