#ifndef Uintah_Component_Arches_Convection_CQMOM_h
#define Uintah_Component_Arches_Convection_CQMOM_h

#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Exceptions/InvalidValue.h>
#include <CCA/Components/Arches/Directives.h>
#include <Core/Parallel/Parallel.h>

//#define cqmom_transport_dbg
//==========================================================================

/**
 * @class Convection_CQMOM
 * @author Alex Abboud
 * @date June, 2014
 *
 * @brief A class containing the functions need for calculating convective fluxes specific to CQMOM eqns, when
 *        one or more of the internal coordinates is specified as particle velocity.
 *        These equations require convection to be based on the abscissas rather than the moments directly,
 *        and will be different enough from regular scalar transport to require its own class.
 *        The operator splitting requires doing the convection in each direction separtely, rather than calculating
 *        a single convective flux for each cell once.
 *        The various helper functions have been reused or adapted from Discretization_new.
 *
 *        The convection term is given by \f$ F_conv = G_{i+1/2} - G_{i-1/2} \f$, with G as
 *        \f$ G_{i+1/2} \eqiv G( M_{i+1/2,l} , M_{i-1/2,r} )
 *        G( M_l, M_r ) = H^+ (M_l) + H^- (M_r) with
 *        H^+ = \sum_i w_i max(u_i,0) u_i^k with k = moment order (expanded for multivariate)
 *        H^- = \sum_i w_i min(u_i,0) u_i^k \f$
 *
 */

//NOTE: placing this here lets both cqmomeqn & this functino use it, possibly move it later?
struct constCCVarWrapper {
  Uintah::constCCVariable<double> data;
};

namespace Uintah{
  class Convection_CQMOM {
    
    // ** These functions should only be used for the CQMOM transport eqns
    // and only when the internal coordinates contain the particle velocity
    // To use diffusion in the CQMOM terms, the regular discretization class is used
    
  public:
    
    Convection_CQMOM();
    ~Convection_CQMOM();
    
    //---------------------------------------------------------------------------
    // Convection functions for CQMOM transport equations
    // --------------------------------------------------------------------------
    
    /** @brief Computes the x-convection term. */
    template<class fT > void
    doConvX( const Patch* p, fT& Fconv, constCCVariable<Vector>& areaFraction, const std::string d_convScheme,
             std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
             const int M, const int nNodes, const int uVelIndex, const std::vector<int> momentIndex );
    
    /** @brief Computes the y-convection term. */
    template<class fT > void
    doConvY( const Patch* p, fT& Fconv, constCCVariable<Vector>& areaFraction, const std::string d_convScheme,
             std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
             const int M, const int nNodes, const int yVelIndex, const std::vector<int> momentIndex );
    
    /** @brief Computes the z-convection term. */
    template<class fT > void
    doConvZ( const Patch* p, fT& Fconv, constCCVariable<Vector>& areaFraction, const std::string d_convScheme,
             std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
             const int M, const int nNodes, const int wVelIndex, const std::vector<int> momentIndex );
    
    
    //---------------------------------------------------------------------------
    // Custom Data Managers
    // --------------------------------------------------------------------------
    
    //  A boolean for marking faces of the cell as touching a boundary
    struct cqFaceBoundaryBool
    {
      bool minus;
      bool plus;
    };
    
    // Stores values in a one-D line on each face for a given cell
    struct cqFaceData1D {
      //plus and minus face values for vars like area frac which are the same
      double minus;
      double plus;
      
      double minus_right; // minus face(interior)
      double minus_left;  // minus face(exterior)
      
      double plus_right;  // plus face(exterior)
      double plus_left;   // plus face(interior)
    };
    
    //---------------------------------------------------------------------------
    // Helper/Utility type functions
    // --------------------------------------------------------------------------
    
    /** @brief Returns an iterator for all cells not touching the domain boundaries */
    inline CellIterator getInteriorCellIterator( const Patch* p ) const
    {
      
      bool xminus = p->getBCType(Patch::xminus) == Patch::None;
      bool xplus =  p->getBCType(Patch::xplus)  == Patch::None;
      bool yminus = p->getBCType(Patch::yminus) == Patch::None;
      bool yplus =  p->getBCType(Patch::yplus)  == Patch::None;
      bool zminus = p->getBCType(Patch::zminus) == Patch::None;
      bool zplus =  p->getBCType(Patch::zplus)  == Patch::None;
      
      IntVector clow  = p->getCellLowIndex();
      IntVector chigh = p->getCellHighIndex();
      IntVector clow_mod = clow;
      IntVector chigh_mod = chigh;
      
      if (xminus)
        clow_mod = clow_mod + IntVector(1,0,0);
      if (xplus)
        chigh_mod = chigh_mod - IntVector(1,0,0);
      if (yminus)
        clow_mod = clow_mod + IntVector(0,1,0);
      if (yplus)
        chigh_mod = chigh_mod - IntVector(0,1,0);
      if (zminus)
        clow_mod = clow_mod + IntVector(0,0,1);
      if (zplus)
        chigh_mod = chigh_mod - IntVector(0,0,1);
      
      CellIterator the_iterator = CellIterator(clow_mod, chigh_mod);
      
      the_iterator.reset();
      return the_iterator;
    }
    
    /** @brief Returns an iterator for all cell touching a domain boundary */
    inline CellIterator getInteriorBoundaryCellIterator( const Patch* p,
                                                        const std::vector<Patch::FaceType>::const_iterator bf_iter ) const
    {
      
      Patch::FaceType face = *bf_iter;
      IntVector l,h;
      p->getFaceCells( face, 0, l, h );
      
      if ( face == Patch::xminus || face == Patch::xplus ){
        
        bool yminus = p->getBCType(Patch::yminus) == Patch::None;
        bool yplus =  p->getBCType(Patch::yplus)  == Patch::None;
        bool zminus = p->getBCType(Patch::zminus) == Patch::None;
        bool zplus =  p->getBCType(Patch::zplus)  == Patch::None;
        
        if (yminus)
          l[1] += 1;
        if (yplus)
          h[1] -= 1;
        if (zminus)
          l[2] += 1;
        if (zplus)
          h[2] -= 1;
        
      } else if ( face == Patch::yminus || face == Patch::yplus ){
        
        bool xminus = p->getBCType(Patch::xminus) == Patch::None;
        bool xplus =  p->getBCType(Patch::xplus)  == Patch::None;
        bool zminus = p->getBCType(Patch::zminus) == Patch::None;
        bool zplus =  p->getBCType(Patch::zplus)  == Patch::None;
        
        if (xminus)
          l[0] += 2;
        if (xplus)
          h[0] -= 2;
        if (zminus)
          l[2] += 1;
        if (zplus)
          h[2] -= 1;
      } else if ( face == Patch::zminus || face == Patch::zplus ){
        
        bool yminus = p->getBCType(Patch::yminus) == Patch::None;
        bool yplus =  p->getBCType(Patch::yplus)  == Patch::None;
        bool xminus = p->getBCType(Patch::xminus) == Patch::None;
        bool xplus =  p->getBCType(Patch::xplus)  == Patch::None;
        
        if (yminus)
          l[1] += 2;
        if (yplus)
          h[1] -= 2;
        if (xminus)
          l[0] += 2;
        if (xplus)
          h[0] -= 2;
      }
      
      CellIterator the_iterator = CellIterator( l, h );
      return the_iterator;
      
    }

    /** @brief Checks if cell has a face on a boundary in a given absolute normal direction (coord) */
    inline cqFaceBoundaryBool checkFacesForBoundaries( const Patch* p, const IntVector c, const IntVector coord )
    {
      
      cqFaceBoundaryBool b;
      b.minus = false;
      b.plus  = false;
      
      IntVector l = p->getCellLowIndex();
      IntVector h = p->getCellHighIndex();
      
      bool fplus  = false;
      bool fminus = false;
      
      int I = -1;
      if ( coord[0] == 1 ) {
        fminus = p->getBCType(Patch::xminus) != Patch::Neighbor;
        fplus  = p->getBCType(Patch::xplus ) != Patch::Neighbor;
        I = 0;
      } else if ( coord[1] == 1 ){
        fminus = p->getBCType(Patch::yminus) != Patch::Neighbor;
        fplus  = p->getBCType(Patch::yplus ) != Patch::Neighbor;
        I = 1;
      } else if ( coord[2] == 1 ){
        fminus = p->getBCType(Patch::zminus) != Patch::Neighbor;
        fplus  = p->getBCType(Patch::zplus ) != Patch::Neighbor;
        I = 2;
      }
      
      if ( fminus && c[I] == l[I] ) {
        b.minus = true;
      }
      if ( fplus  && c[I] == h[I]-1 ) {
        b.plus = true;
      }
      
      return b;
    }

    //---------------------------------------------------------------------------
    // Assemblers
    // These functions assemble other terms
    // --------------------------------------------------------------------------
    //shouln't need flux versino with densoty for particle transport
    /** @brief Computes the flux term, \f$ int_A div{u \phi} \cdot dA \f$, where u is the velocity
     *          in the normal (coord) direction.  Note version does not have density. */
    inline double getFlux( const double area, cqFaceData1D GPhi, constCCVariable<Vector>& areaFraction, IntVector coord, IntVector c )
    {
      double F;
      cqFaceData1D areaFrac;
      IntVector cp = c + coord;
      Vector curr_areaFrac = areaFraction[c];
      Vector plus_areaFrac = areaFraction[cp];
      // may want to just pass dim in for efficiency sake
      int dim = 0;
      if (coord[0] == 1)
        dim =0;
      else if (coord[1] == 1)
        dim = 1;
      else
        dim = 2;
      areaFrac.plus  = plus_areaFrac[dim];
      areaFrac.minus = curr_areaFrac[dim];
      
      return F = area * (  areaFrac.plus * GPhi.plus
                         - areaFrac.minus * GPhi.minus );
    }
    
    inline cqFaceData1D sumNodes( std::vector<cqFaceData1D>& w, std::vector<cqFaceData1D>& a, const int& nNodes, const int& M,
                                  const int& velIndex, const std::vector<int>& momentIndex )
    {
      cqFaceData1D gPhi;
      double hPlus = 0.0, hMinus = 0.0;
      
      //plus cell face
      for (int i = 0; i < nNodes; i++) {
        if (a[i + nNodes*velIndex].plus_left > 0.0) {
          double nodeVal = 1.0;
          for (int m = 0 ; m < M ; m++) {
            nodeVal *= pow( (a[i + nNodes*m].plus_left),momentIndex[m] );
          }
          hPlus += w[i].plus_left * a[i + nNodes*velIndex].plus_left * nodeVal; //add hplus
        }
        if ( a[i + nNodes*velIndex].plus_right < 0.0) {
          double nodeVal = 1.0;
          for (int m = 0 ; m < M ; m++) {
            nodeVal *= pow( (a[i + nNodes*m].plus_right),momentIndex[m] );
          }
          hMinus += w[i].plus_right * a[i + nNodes*velIndex].plus_right * nodeVal; //add hminus
        }
      }
      
#ifdef cqmom_transport_dbg
      std::cout << "Pface h+: " << hPlus << "  h-: " << hMinus << std::endl;

      std::cout.precision(15);
      for (int i = 0; i < nNodes; i++) {
        std::cout << "vel +r (" << i+ nNodes*velIndex << ") " << a[i+ nNodes*velIndex].plus_right << std::endl;
        std::cout << "w +r " << w[i].plus_right << std::endl;
        for (int m = 0 ; m < M ; m++) {
          std::cout << "a[" << i << "][" << m << "]+r " << a[i+nNodes*m].plus_right << " " << momentIndex[m] << std::endl;
        }
      }
      for (int i = 0; i < nNodes; i++) {
        std::cout << "vel +l " << a[i+ nNodes*velIndex].plus_left << std::endl;
        std::cout << "w +l " << w[i].plus_left << std::endl;
        for (int m = 0 ; m < M ; m++) {
          std::cout << "a[" << i << "][" << m << "]+l " << a[i+ nNodes*m].plus_left << " " << momentIndex[m] << std::endl;
        }
      }
      std::cout.precision(5);
#endif
      
      gPhi.plus = hPlus + hMinus;
      
      //minus cell face
      hPlus = 0.0; hMinus = 0.0;
      for (int i = 0; i < nNodes; i++) {
        if (a[i + nNodes*velIndex].minus_left > 0.0) {
          double nodeVal = 1.0;
          for (int m = 0 ; m < M ; m++) {
            nodeVal *= pow( (a[i + nNodes*m].minus_left),momentIndex[m] );
          }
          hPlus += w[i].minus_left * a[i + nNodes*velIndex].minus_left * nodeVal; //add hplus
        }
        if ( a[i + nNodes*velIndex].minus_right < 0.0) {
          double nodeVal = 1.0;
          for (int m = 0 ; m < M ; m++) {
            nodeVal *= pow( (a[i + nNodes*m].minus_right),momentIndex[m] );
          }
          hMinus += w[i].minus_right * a[i + nNodes*velIndex].minus_right * nodeVal; //add hminus
        }
      }
      gPhi.minus = hPlus + hMinus;
#ifdef cqmom_transport_dbg
      std::cout << "Mface h+: " << hPlus << "  h-: " << hMinus << std::endl;
      std::cout << "G+: " << gPhi.plus << " G-: " << gPhi.minus << " Tot: " << gPhi.plus - gPhi.minus << std::endl;
#endif
      return gPhi;
    }
    
    
    //---------------------------------------------------------------------------
    // Interpolation
    // These functions interpolate
    // --------------------------------------------------------------------------
    
    
    /** @brief first-order CQMOM -- abscissas are equal across cell */
    template< class phiT >
    inline cqFaceData1D firstOrder( const IntVector c, const IntVector coord, phiT& phi )
    {
      IntVector cxp = c + coord;
      IntVector cxm = c - coord;
      
      cqFaceData1D face_values;
      
      face_values.minus_right =  phi[c];
      face_values.plus_left  =  phi[c];
      
      face_values.minus_left = phi[cxm];
      face_values.plus_right = phi[cxp];
      
      return face_values;
    }
    
    //---------------------------------------------------------------------------
    // Interpolation Class
    //
    /** @brief Calls the specific interpolant */
    template <typename operT, typename phiT>
    class ConvHelper1 {
    public:
      ConvHelper1<operT, phiT>( operT* opr, phiT& phi ) : _opr(opr), _phi(phi){};
      ~ConvHelper1(){};

      
      template<class fT>
      void do_convection( const Patch* p, fT& Fconv, std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
                         const int& M, const int& nNodes, const int& velIndex, const std::vector<int>& momentIndex, const int& dim,
                         constCCVariable<Vector>& area_fraction, Convection_CQMOM* D)
      {
        //set conv direction
        IntVector coord = IntVector(0,0,0);
        double area;
        Vector Dx = p->dCell();
        
        if (dim == 0) {
          coord[0] = 1;
          area = Dx.y()*Dx.z();
        } else if ( dim == 1) {
          coord[1] = 1;
          area = Dx.x()*Dx.z();
        } else {
          coord[2] = 1;
          area = Dx.x()*Dx.y();
        }
        
#ifdef cqmom_transport_dbg
        std::cout << "Interior Cell Iterator" << std::endl;
#endif
        
        CellIterator iIter  = D->getInteriorCellIterator( p );
        //-------------------- Interior
        for (iIter.begin(); !iIter.done(); iIter++){
          
          IntVector c   = *iIter;
          
          Convection_CQMOM::cqFaceData1D gPhi;
          double aSize = M*nNodes;
          std::vector<Convection_CQMOM::cqFaceData1D> faceAbscissas (aSize);
          std::vector<Convection_CQMOM::cqFaceData1D> faceWeights (nNodes);
          
          int ii = 0;
          for (std::vector<constCCVarWrapper>::iterator iter = weights.begin(); iter != weights.end(); ++iter) {
            faceWeights[ii] = _opr->no_bc( c, coord, (iter->data) );
            ii++;
          }
          
          ii = 0;
          for (std::vector<constCCVarWrapper>::iterator iter = abscissas.begin(); iter != abscissas.end(); ++iter) {
            faceAbscissas[ii] = _opr->no_bc( c, coord, (iter->data) );
            ii++;
          }
#ifdef cqmom_transport_dbg
          std::cout << "Cell location: " << c << " in dimension " << dim << std::endl;
          std::cout << "____________________________" << std::endl;
#endif
          gPhi     = D->sumNodes( faceWeights, faceAbscissas, nNodes, M, velIndex, momentIndex );
          Fconv[c] = D->getFlux( area, gPhi, area_fraction, coord, c );
        }
 
#ifdef cqmom_transport_dbg
        std::cout << "Boundary Cell Iterator" << std::endl;
#endif
        
        //--------------- Boundaries
        std::vector<Patch::FaceType> bf;
        std::vector<Patch::FaceType>::const_iterator bf_iter;
        p->getBoundaryFaces(bf);
        
        // Loop over all boundary faces on this patch
        for (bf_iter = bf.begin(); bf_iter != bf.end(); bf_iter++){
          
          Patch::FaceType face = *bf_iter;
          IntVector inside = p->faceDirection(face);
          CellIterator c_iter = D->getInteriorBoundaryCellIterator( p, bf_iter );
          Convection_CQMOM::cqFaceBoundaryBool faceIsBoundary;
          
          for (c_iter.begin(); !c_iter.done(); c_iter++){
            
            IntVector c = *c_iter - inside;
            Convection_CQMOM::cqFaceData1D gPhi;
            
            faceIsBoundary = D->checkFacesForBoundaries( p, c, coord );
            double aSize = M*nNodes;
            std::vector<Convection_CQMOM::cqFaceData1D> faceAbscissas (aSize);
            std::vector<Convection_CQMOM::cqFaceData1D> faceWeights (nNodes);
            
            int ii = 0;
            for (std::vector<constCCVarWrapper>::iterator iter = weights.begin(); iter != weights.end(); ++iter) {
              faceWeights[ii] = _opr->with_bc( c, coord, (iter->data), area_fraction, faceIsBoundary );
              ii++;
            }
            
            ii = 0;
            for (std::vector<constCCVarWrapper>::iterator iter = abscissas.begin(); iter != abscissas.end(); ++iter) {
              faceAbscissas[ii] = _opr->with_bc( c, coord, (iter->data), area_fraction, faceIsBoundary );
              ii++;
            }
#ifdef cqmom_transport_dbg
            std::cout << "Cell location: " << c << " in dimension " << dim << std::endl;
            std::cout << "____________________________" << std::endl;
#endif
            gPhi     = D->sumNodes( faceWeights, faceAbscissas, nNodes, M, velIndex, momentIndex );
            Fconv[c] = D->getFlux( area,  gPhi, area_fraction, coord, c );
          }
        }
      }

      
    private:
      operT* _opr;
      phiT& _phi;
      
    };
    
    // ---------------------------------------------------------------------------
    // First Order CQMOM Interpolation
    //
    /** @brief first order interolation - face values = CC values
        This is pseudo upwind, and bases the flux direction on each quadrature node */
    template <typename phiT>
    class FirstOrderInterpolation {
    public:
      FirstOrderInterpolation(){};
      ~FirstOrderInterpolation(){};
      
      cqFaceData1D no_bc( const IntVector c, const IntVector coord, constCCVariable<double>& phi)
      {
        Convection_CQMOM::cqFaceData1D face_values;
        IntVector cxp = c + coord;
        IntVector cxm = c - coord;
        
        face_values.minus_right =  phi[c];
        face_values.plus_left  =  phi[c];
        
        face_values.minus_left = phi[cxm];
        face_values.plus_right = phi[cxp];
        
        return face_values;
      };
      
      cqFaceData1D inline with_bc ( const IntVector c, const IntVector coord, constCCVariable<double>& phi,
                                  constCCVariable<Vector>& areaFraction, cqFaceBoundaryBool isBoundary)
      {
        Convection_CQMOM::cqFaceData1D face_values;
        IntVector cxp = c + coord;
        IntVector cxm = c - coord;
        /*NOTE: since the inlet boundary condition is equal to the face, don't find the average
        here as the regular scalars do when theres a dirichlet(or any) condition. Leaving the if statement and the
         average as placeholder in case the boundary condtions are handled differently later it will be easy to switch*/
        
        // - face
        if (isBoundary.minus) {
  //        face_values.minus_right = 0.5*(phi[c] + phi[cxm]);
  //        face_values.minus_left = 0.5*(phi[c] + phi[cxm]);
          face_values.minus_right =  phi[c];
          face_values.minus_left = phi[cxm];
        } else {
          face_values.minus_right =  phi[c];
          face_values.minus_left = phi[cxm];
        }
        
        // + face
        if (isBoundary.plus) {
//          face_values.plus_left = 0.5*(phi[c] + phi[cxp]);
//          face_values.plus_right = 0.5*(phi[c] + phi[cxp]);
          face_values.plus_right =  phi[cxp];
          face_values.plus_left = phi[c];
        } else {
          face_values.plus_right  =  phi[cxp];
          face_values.plus_left = phi[c];
        }
        
        return face_values;
      };
    };

    
    // ---------------------------------------------------------------------------
    // Second Order CQMOM Interpolation
    //
    /** @brief second order interolation with the pseudo secodn order scheme */
    template <typename phiT>
    class SecondOrderInterpolation {
    public:
      cqFaceData1D no_bc( const IntVector c, const IntVector coord, phiT& phi) {
        //placeholder
      }
      
      cqFaceData1D inline with_bc( const IntVector c, const IntVector coord, phiT& phi,
                                   constCCVariable<Vector>& areaFraction, cqFaceBoundaryBool isBoundary)
      {
        //placeholder
      }
    };
    
  }; // class Convection_CQMOM
  
  
  struct cqFaceData1D {
    //plus and minus face values for vars like area frac which are the same
    double minus;
    double plus;
    
    double minus_right; // minus face(interior)
    double minus_left;  // minus face(exterior)
    
    double plus_right;  // plus face(exterior)
    double plus_left;   // plus face(interior)
  };
  
  //========================= Convection ======================================
  
  
  //---------------------------------------------------------------------------
  // Method: Compute the convection term for each physical direction
  //---------------------------------------------------------------------------

  //x direction convection
  template<class fT > void
  Convection_CQMOM::doConvX( const Patch* p, fT& Fconv, constCCVariable<Vector>& areaFraction, const std::string d_convScheme,
                             std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
                             const int M, const int nNodes, const int uVelIndex, const std::vector<int> momentIndex )
  {
    int dim = 0;
    
    //CQMOM with IC as velocity needs its own set of convection schemes, for now labled as simply first/second
    if (d_convScheme == "first") {
      
      FirstOrderInterpolation<fT>* the_interpolant = scinew FirstOrderInterpolation<fT>();
      ConvHelper1<FirstOrderInterpolation<fT>, fT>* convection_helper =
      scinew ConvHelper1<FirstOrderInterpolation<fT>, fT>(the_interpolant, Fconv);
      
      convection_helper->do_convection( p, Fconv, weights, abscissas, M, nNodes, uVelIndex, momentIndex, dim, areaFraction, this );
      
      delete convection_helper;
      delete the_interpolant;
      
    } else if (d_convScheme == "second") {
      //insert second order later
      throw InvalidValue("Error: Second order scheme not yet implemented", __FILE__, __LINE__);
      
    } else {
      
      throw InvalidValue("Error: Convection scheme not recognized. Check UPS file and try again.", __FILE__, __LINE__);
      
    }
  }
  
  //y direction convection
  template<class fT> void
  Convection_CQMOM::doConvY( const Patch* p, fT& Fconv, constCCVariable<Vector>& areaFraction, const std::string d_convScheme,
                            std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
                            const int M, const int nNodes, const int vVelIndex, const std::vector<int> momentIndex )
  {
    int dim = 1;
    
    if (d_convScheme == "first") {
      FirstOrderInterpolation<fT>* the_interpolant = scinew FirstOrderInterpolation<fT>();
      ConvHelper1<FirstOrderInterpolation<fT>, fT>* convection_helper =
      scinew ConvHelper1<FirstOrderInterpolation<fT>, fT>(the_interpolant, Fconv);
      
      convection_helper->do_convection( p, Fconv, weights, abscissas, M, nNodes, vVelIndex, momentIndex, dim, areaFraction, this );
      
      delete convection_helper;
      delete the_interpolant;
      
    } else if (d_convScheme == "second") {
      throw InvalidValue("Error: Second order scheme not yet implemented", __FILE__, __LINE__);
      
    } else {
      
      throw InvalidValue("Error: Convection scheme not recognized. Check UPS file and try again.", __FILE__, __LINE__);
      
    }
  }
  
  //z directino convection
  template<class fT> void
  Convection_CQMOM::doConvZ( const Patch* p, fT& Fconv, constCCVariable<Vector>& areaFraction, const std::string d_convScheme,
                            std::vector<constCCVarWrapper> weights, std::vector<constCCVarWrapper> abscissas,
                            const int M, const int nNodes, const int wVelIndex, const std::vector<int> momentIndex )
  {
    int dim = 2;
    
    if (d_convScheme == "first") {
      FirstOrderInterpolation<fT>* the_interpolant = scinew FirstOrderInterpolation<fT>();
      ConvHelper1<FirstOrderInterpolation<fT>, fT>* convection_helper =
      scinew ConvHelper1<FirstOrderInterpolation<fT>, fT>(the_interpolant, Fconv);
      
      convection_helper->do_convection( p, Fconv, weights, abscissas, M, nNodes, wVelIndex, momentIndex, dim, areaFraction, this );
      
      delete convection_helper;
      delete the_interpolant;
      
    } else if (d_convScheme == "second") {
      throw InvalidValue("Error: Second order scheme not yet implemented", __FILE__, __LINE__);
      
    } else {
      
      throw InvalidValue("Error: Convection scheme not recognized. Check UPS file and try again.", __FILE__, __LINE__);
      
    }
  }
  
  //========================= Diffusion ======================================
  
  //Note: no diffusion for particle transport
  //If diffusion is to be used for testing use the normal diffusive flux in Discretization_new
  
}//namespace Uintah
#endif