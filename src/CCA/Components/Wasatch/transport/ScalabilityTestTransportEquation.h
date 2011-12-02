#ifndef Wasatch_ScalabilityTestTransportEquation_h
#define Wasatch_ScalabilityTestTransportEquation_h

//-- ExprLib includes --//
#include <CCA/Components/Wasatch/transport/TransportEquation.h>

//-- Wasatch includes --//
#include <CCA/Components/Wasatch/FieldTypes.h>
#include <CCA/Components/Wasatch/Expressions/ScalarRHS.h>
#include <CCA/Components/Wasatch/BCHelperTools.h>

//-- Uintah includes --//
#include <Core/ProblemSpec/ProblemSpecP.h>


namespace Wasatch{

  /**
   *  \ingroup WasatchCore
   *  \class ScalabilityTestTransportEquation
   *  \date April, 2011
   *  \author Tony Saad
   *
   *  \brief Special transport equation used for scalability testing.
   *
   *  Sets up solution for a transport equation of the form:
   *
   *  \f[
   *    \frac{\partial \rho \phi}{\partial t} =
   *    - \frac{\partial \rho \phi u_x }{\partial x}
   *    - \frac{\partial \rho \phi u_y }{\partial y}
   *    - \frac{\partial \rho \phi u_z }{\partial z}
   *    - \frac{\partial J_{\phi,x}}{\partial x}
   *    - \frac{\partial J_{\phi,y}}{\partial y}
   *    - \frac{\partial J_{\phi,z}}{\partial z}
   *    + s_\phi
   *  \f]
   *
   *  Any or all of the terms in the RHS above may be activated
   *  through the input file.
   *
   *  \par Notes & Restrictions
   *
   *  - Currently, only basic forms for the scalar diffusive velocity are
   *    supported.  Specifically, either an expression for the
   *    diffusion coefficient, \f$\Gamma_\phi\f$ is required or the
   *    diffusion coefficient must be a constant value.  See
   *    DiffusiveVelocity and DiffusiveVelocity2 classes.
   *
   *  - Source terms can only be added if the expression to evaluate
   *    them has been constructed elsewhere.
   *
   *  \todo Need to hook in parser support for boundary and initial conditions.
   */
  template<typename FieldT>
  class ScalabilityTestTransportEquation : public Wasatch::TransportEquation
  {
  public:

    // these typedefs are provided for convenience.
    typedef typename FaceTypes<FieldT>::XFace  XFaceT; ///< The type of field on the x-faces of the volume.
    typedef typename FaceTypes<FieldT>::YFace  YFaceT; ///< The type of field on the y-faces of the volume.
    typedef typename FaceTypes<FieldT>::ZFace  ZFaceT; ///< The type of field on the z-faces of the volume.

    /**
     *  \brief Construct a ScalabilityTestTransportEquation
     *  \param basePhiName This equation will be created n-times where n is a user
     *         specified number in the input file. The basePhiName refers to the
     *         base name of the solution variable. The n-equations that are created
     *         will correspond to basePhiName0, basePhiName1, etc...
     *  \param thisPhiName The name of the solution variable for this ScalarTransportEquation
     *  \param id The Expr::ExpressionID for the RHS expression for this ScalarTransportEquation
     *
     *  Note that the static member method get_rhs_expr_id can be useful to
     *  obtain the appropriate input arguments here.
     */
    ScalabilityTestTransportEquation( const std::string basePhiName,
                                      const std::string thisPhiName,
                                      const Expr::ExpressionID id );

    ~ScalabilityTestTransportEquation();

    /**
     *  \brief apply the boundary conditions on the initial condition
     *         associated with this transport equation
     */
    void setup_initial_boundary_conditions( const GraphHelper& graphHelper,
                                           const Uintah::PatchSet* const localPatches,
                                           const PatchInfoMap& patchInfoMap,
                                           const Uintah::MaterialSubset* const materials);


    /**
     *  \brief setup the boundary conditions associated with this transport equation
     */
    void setup_boundary_conditions( const GraphHelper& graphHelper,
                                   const Uintah::PatchSet* const localPatches,
                                   const PatchInfoMap& patchInfoMap,
                                   const Uintah::MaterialSubset* const materials);

    /**
     *  \brief setup the initial conditions for this transport equation.
     */
    Expr::ExpressionID initial_condition( Expr::ExpressionFactory& icFactory );

    /**
     * \brief Parse the input file to determine the rhs expression id.
     *        Also registers convective flux, diffusive velocity, and
     *        source term expressions.
     *
     *  \param thisPhiName the name of the solution variable for this
     *         equation.
     *
     *  \param factory The Expr::ExpressionFactory object that
     *         terms associated with the RHS of this transport
     *         equation should be registered on.
     *
     *  \param params The Uintah::ProblemSpec XML description for this
     *         equation.  Scope should be within the ScalabilityTest tag.
     */
    static Expr::ExpressionID get_rhs_expr_id( const std::string thisPhiName,
                                               Expr::ExpressionFactory& factory,
                                               Uintah::ProblemSpecP params );
  private:

  };


  template< typename FieldT >
  void setup_diffusive_velocity_expression( const std::string dir,
                                            Expr::ExpressionFactory& factory,
                                            typename ScalarRHS<FieldT>::FieldTagInfo& info );

  template< typename FieldT >
  void setup_convective_flux_expression( const std::string dir,
                                         const std::string thisPhiName,
                                         const Expr::Tag advVelocityTag,
                                         Expr::ExpressionFactory& factory,
                                         typename ScalarRHS<FieldT>::FieldTagInfo& info );

} // namespace Wasatch
#endif // Wasatch_ScalabilityTestTransportEquation_h
