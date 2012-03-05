/*----------------------------------------------------------------------*/
/*!
\file fluid3_boundary_impl.H

\brief Internal implementation of fluid3 boundary elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_BOUNDARY_IMPL_H
#define FLUID3_BOUNDARY_IMPL_H

#include "fluid3.H"
#include "../drt_f3_impl/fluid3_impl_parameter.H"

namespace DRT
{
namespace ELEMENTS
{

/// Interface base class for Fluid3BoundaryImpl
/*!
  This class exists to provide a common interface for all template
  versions of Fluid3BoundaryImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of Fluid3BoundaryImpl.
 */
class Fluid3BoundaryImplInterface
{

public:
  /// Empty constructor
  Fluid3BoundaryImplInterface() {}
  /// Empty destructor
  virtual ~Fluid3BoundaryImplInterface() {}
  /// Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in Fluid3BoundaryImpl.
   */
  virtual int Evaluate(DRT::ELEMENTS::Fluid3Boundary* ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  virtual void EvaluateInternalFacesUsingNeighborData(DRT::ELEMENTS::Fluid3Boundary* bele,
                       DRT::ELEMENTS::Fluid3*               ele1,
                       DRT::ELEMENTS::Fluid3*               ele2,
                       std::vector<int>&                    nds1,
                       std::vector<int>&                    nds2,
                       ParameterList&                       params       ,
                       DRT::Discretization&                 discretization,
                       RCP<LINALG::SparseOperator>          systemmatrix1,
                       RCP<LINALG::SparseOperator>          systemmatrix2,
                       RCP<Epetra_Vector>                   systemvector1,
                       RCP<Epetra_Vector>                   systemvector2,
                       RCP<Epetra_Vector>                   systemvector3) = 0;

  virtual int EvaluateInternalFaces(DRT::ELEMENTS::Fluid3Boundary* bele,
                                    DRT::ELEMENTS::Fluid3*                      ele1,
                                    DRT::ELEMENTS::Fluid3*                      ele2,
                                    std::vector<int>&                           nds1,
                                    std::vector<int>&                           nds2,
                                    ParameterList&                              params,
                                    DRT::Discretization&                        discretization,
                                    vector<int>&                                lm,
                                    vector<int>&                                lm1,
                                    vector<int>&                                lm2,
                                    Epetra_SerialDenseMatrix&                   elemat1,
                                    Epetra_SerialDenseVector&                   elevec1) = 0;




  /// Evaluate a Neumann boundary condition
  /*!
    This class does not provide a definition for this function, it
    must be defined in Fluid3BoundaryImpl.
   */
  virtual int EvaluateNeumann(DRT::ELEMENTS::Fluid3Boundary*      ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1_epetra,
                              Epetra_SerialDenseMatrix*  elemat1
                              ) = 0;

  /// Internal implementation class for Fluid3Boundary elements
  static Fluid3BoundaryImplInterface* Impl(const DRT::Element* ele);

};

  /// Internal Fluid3Boundary element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the Fluid3Boundary element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid3Boundary)
    from the mathematical contents (this class). Of course there are
    different implementations of the Fluid3Boundary element, this is just one
    such implementation.

    The Fluid3Boundary element will allocate exactly one object of this class
    for all Fluid3Boundary elements with the same number of nodes in the mesh.
    This allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>History</h3>

    The implementation here is the standard convection-diffusion element
    capable of dealing with systems of transported scalars.

    Right now we do not read any stabilization parameters from the
    input file but have a fixed version.

    \author gjb
    \date 08/08
    \author ehrl
    \date 03/10
  */
template<DRT::Element::DiscretizationType distype>
class Fluid3BoundaryImpl: public Fluid3BoundaryImplInterface
{
  friend class Fluid3ImplParameter;

public:

  /// Singleton access method
  static Fluid3BoundaryImpl<distype> * Instance();

  /// Constructor
  Fluid3BoundaryImpl();

  //! number of element nodes
  static const int bdrynen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions of the Fluid3Boundary element
  static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of space dimensions of the parent element
  static const int nsd_ = bdrynsd_+1;

  //! number of dof's per node
  static const int numdofpernode_ = nsd_ + 1;

  //! Evaluate
  virtual int Evaluate(DRT::ELEMENTS::Fluid3Boundary* ele,
                       ParameterList&             params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       );

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(DRT::ELEMENTS::Fluid3Boundary*    ele,
                              ParameterList&             params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1_epetra,
                              Epetra_SerialDenseMatrix*  elemat1
                              );

  virtual void EvaluateInternalFacesUsingNeighborData(DRT::ELEMENTS::Fluid3Boundary* bele,
                       DRT::ELEMENTS::Fluid3*               ele1,
                       DRT::ELEMENTS::Fluid3*               ele2,
                       std::vector<int>&                    nds1,
                       std::vector<int>&                    nds2,
                       ParameterList&                       params       ,
                       DRT::Discretization&                 discretization,
                       RCP<LINALG::SparseOperator>          systemmatrix1,
                       RCP<LINALG::SparseOperator>          systemmatrix2,
                       RCP<Epetra_Vector>                   systemvector1,
                       RCP<Epetra_Vector>                   systemvector2,
                       RCP<Epetra_Vector>                   systemvector3);

  virtual int EvaluateInternalFaces(DRT::ELEMENTS::Fluid3Boundary* bele,
                       DRT::ELEMENTS::Fluid3*                      ele1,
                       DRT::ELEMENTS::Fluid3*                      ele2,
                       std::vector<int>&                           nds1,
                       std::vector<int>&                           nds2,
                       ParameterList&                              params,
                       DRT::Discretization&                        discretization,
                       vector<int>&                                lm,
                       vector<int>&                                lm1,
                       vector<int>&                                lm2,
                       Epetra_SerialDenseMatrix&                   elemat1,
                       Epetra_SerialDenseVector&                   elevec1
                       );

private:

    void CreateLinesTri(const int& nline,
                        const int& nnode);

    void CreateLinesQuad(const int& nline,
                         const int& nnode);

    /*!
    \brief integrate shapefunction over surface element

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,
    */

    virtual void IntegrateShapeFunction(
                                    DRT::ELEMENTS::Fluid3Boundary*    ele,
                                    ParameterList&                    params,
                                    DRT::Discretization&              discretization,
                                    vector<int>&                      lm,
                                    Epetra_SerialDenseVector&         elevec1,
                                    const std::vector<double>&        edispnp);

    /*!
      \brief Calculate node normals

      Acc. to Wall (7.13), surface integral.

      \param elevec1 (out)      : Weighted element normal for assembly
      \param edispnp (in)       : Displacement-vector
    */
    void ElementNodeNormal(
                         DRT::ELEMENTS::Fluid3Boundary*   ele,
                         ParameterList&                 params,
                         DRT::Discretization&           discretization,
                         vector<int>&                   lm,
                         Epetra_SerialDenseVector&      elevec1,
                         const std::vector<double>&     edispnp);
    /*!
      \brief Calculate mean curvature H. Interpolate the results to achieve better
      results in the surface tension algorithm (c0 field -> c1 field).

      \param elevec1  (out)     : Nodal values of mean curvature
      \param edispnp  (in)      : Displacement-vector
      \param enormals (in)      : Node normals
    */
    void ElementMeanCurvature(
                              DRT::ELEMENTS::Fluid3Boundary*    ele,
                              ParameterList&                    params,
                              DRT::Discretization&              discretization,
                              vector<int>&                      lm,
                              Epetra_SerialDenseVector&         elevec1,
                              const std::vector<double>&        edispnp,
                              std::vector<double>&              enormals);


    /*!
      \brief Integrate surface tension

      \param elevec1   (out)    : RHS contribution of surface tension effect
      \param edispnp   (in)     : Displacement-vector
      \param enormals  (in)     : Node normals
      \param curvature (in)     : Interpolated curvature
    */
    void ElementSurfaceTension(
                               DRT::ELEMENTS::Fluid3Boundary*   ele,
                               ParameterList&                   params,
                               DRT::Discretization&             discretization,
                               vector<int>&                     lm,
                               Epetra_SerialDenseVector&        elevec1,
                               const std::vector<double>&       edispnp,
                               std::vector<double>&             enormals,
                               std::vector<double>&             ecurvature);

    /*!
    brief integrate elemental areas over a surface

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    */
    void AreaCaculation(DRT::ELEMENTS::Fluid3Boundary*    ele,
                        ParameterList&                    params,
                        DRT::Discretization&              discretization,
                        vector<int>&                      lm);


   /*!
    brief integrate elemental pressure over a surface

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    */
    void IntegratedPressureParameterCalculation(
      DRT::ELEMENTS::Fluid3Boundary*    ele,
      ParameterList&                    params,
      DRT::Discretization&              discretization,
      vector<int>&                      lm);

   /*!
    brief center of mass of a surface

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    */
    void CenterOfMassCalculation(
      DRT::ELEMENTS::Fluid3Boundary*    ele,
      ParameterList&                    params,
      DRT::Discretization&              discretization,
      vector<int>&                      lm);


    /*!
    \brief calculate the velocity component of the traction

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,
    */
    void CalcTractionVelocityComponent(
      DRT::ELEMENTS::Fluid3Boundary*    ele,
      ParameterList&                    params,
      DRT::Discretization&              discretization,
      vector<int>&                      lm,
      Epetra_SerialDenseVector&         elevec1);


    /*!
    \brief calculate the integral of the neumann inflow component

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,
    */
     void ComputeNeumannUvIntegral(
                     DRT::ELEMENTS::Fluid3Boundary*    ele,
                     ParameterList&                    params,
                     DRT::Discretization&              discretization,
                     vector<int>&                      lm,
                     Epetra_SerialDenseVector&         elevec1);


    /*!
     brief integrate elemental flow rates over a surface

     \param params (in/out)    : ParameterList for communication between control routine
                                 and elements
     \param discretization (in): A reference to the underlying discretization
     \param lm (in)            : location vector of this element
     */
     void ComputeFlowRate(
                     DRT::ELEMENTS::Fluid3Boundary*    ele,
                     ParameterList&                    params,
                     DRT::Discretization&              discretization,
                     vector<int>&                      lm,
                     Epetra_SerialDenseVector&         elevec1);

    /*!
    \brief determine elemental flow rate and its derivatives wrt dofs

    \param params (in/out)    : parameterList for communication between control routine
                                and elements
    \param discretization (in): reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : 2nd derivative of flowrates wrt
                                velocities and displacements
    \param elemat2 (out)      : 2nd derivative of flowrates wrt displacements
    \param elevec1 (out)      : derivative of flowrates wrt velocities
    \param elevec2 (out)      : derivative of flowrate wrt displacements
    \param elevec3 (out)      : flowrate

    */
    void FlowRateDeriv(
                       DRT::ELEMENTS::Fluid3Boundary*   ele,
                       ParameterList&                   params,
                       DRT::Discretization&             discretization,
                       vector<int>&                     lm,
                       Epetra_SerialDenseMatrix&        elemat1,
                       Epetra_SerialDenseMatrix&        elemat2,
                       Epetra_SerialDenseVector&        elevec1,
                       Epetra_SerialDenseVector&        elevec2,
                       Epetra_SerialDenseVector&        elevec3);


    /*!
    \brief apply impedance boundary condition (outlet pressures)

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,
    */
    void ImpedanceIntegration(
            DRT::ELEMENTS::Fluid3Boundary*    ele,
            ParameterList&                    params,
            DRT::Discretization&              discretization,
            vector<int>&                      lm,
            Epetra_SerialDenseVector&         elevec1);

    /*!
    \brief apply outflow boundary condition which is necessary for the conservative
    element formulation (since the convective term was partially integrated)

    \param params (in)        : ParameterList for communication between control routine
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : matrix to be filled by element.
    \param elevec1 (out)      : vector to be filled by element.

    */
    void ConservativeOutflowConsistency(
      DRT::ELEMENTS::Fluid3Boundary* ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseVector&  elevec1);

    /*!
    \brief compute potential Neumann inflow

    \param params (in)        : ParameterList for communication between control routine
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                                the controling method does not epxect the element to fill
                                this matrix.                          and elements
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,

    */
    void NeumannInflow(
      DRT::ELEMENTS::Fluid3Boundary*  ele,
      ParameterList&             params,
      DRT::Discretization&       discretization,
      vector<int>&               lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseVector&  elevec1);

    /*!
    \brief Evaluate Boundary ele and parent for mixed hybrid Dirichlet conditions

    \param params (in)        : ParameterList for communication between control routine
    \param discretization (in): A reference to the underlying discretization
    \param lm (in)            : location vector of this element
    \param elemat1 (out)      : matrix to be filled by element.
    \param elevec1 (out)      : vector to be filled by element.

    */
  template <DRT::Element::DiscretizationType bndydistype,
            DRT::Element::DiscretizationType pdistype>
  void MixHybDirichlet(
      DRT::ELEMENTS::Fluid3Boundary*  ele,
      ParameterList&                  params,
      DRT::Discretization&            discretization,
      vector<int>&                    lm,
      Epetra_SerialDenseMatrix&       elemat,
      Epetra_SerialDenseVector&       elevec);


    enum ActionType
    {
      none,
      integrate_Shapefunction,
      areacalc,
      flowratecalc,
      calc_flowrate,
      flowratederiv,
      Outletimpedance,
      calc_node_normal,
      calc_node_curvature,
      calc_surface_tension,
      enforce_weak_dbc,
      mixed_hybrid_dbc,
      conservative_outflow_bc,
      calc_Neumann_inflow,
      integ_pressure_calc,
      center_of_mass_calc,
      traction_velocity_component,
      traction_Uv_integral_component
    };

  //! evaluate shape functions/derivatives, metric tensor and
  //! optionally, the infinitesimal area element and the normalized normal vector the at integration point
  void EvalShapeFuncAtBouIntPoint(
      const DRT::UTILS::IntPointsAndWeights<bdrynsd_>& intpoints,  ///< integration points (in)
      const int                                        iquad,      ///< actual Gauss point (in)
      const std::vector<Epetra_SerialDenseVector>*     myknots,    ///< nurbs specific knotes (in)
      const Epetra_SerialDenseVector*                  weights
  );

  //! evaluate Nurbs specific knots and weights
  bool GetKnotVectorAndWeightsForNurbs(
      DRT::ELEMENTS::Fluid3Boundary*              ele,              ///< pointer to the boundary element
      DRT::Discretization&                        discretization,   ///< A reference to the underlying discretization
      std::vector<Epetra_SerialDenseVector>&      mypknots,         ///< nurbs specific knots
      std::vector<Epetra_SerialDenseVector>&      myknots,          ///< nurbs specific ...
      Epetra_SerialDenseVector&                   weights,          ///< nurbs specific weights
      double&                                     normalfac       ///< nurbs specific factor for normal
  );

  //-----------------------------------------------------------------
  //  evaluate the residual of Spaldings law of the wall
  //                                           (private) gammi 11/09
  //-----------------------------------------------------------------
  double SpaldingResidual(
    const double y     ,
    const double visc  ,
    const double utau  ,
    const double normu
    )
    {
      // get dimensionless velocity
      const double up = normu/utau;
      
      const double chi=0.4;
      const double B  =5.5;

      //      +
      // get y , a dimensionless boundary layer thickness
      const double yp = y*utau/visc;
      
      /*
      // Evaluate Spaldings law of the wall
      //                                /                                                  \
      //                               |                           /     +\ 2    /     +\ 3 |
      //                               |       +                  | chi*u  |    | chi*u  |  |
      //  +     / +\     +    -chi*B   |  chi*u               +    \      /      \      /   |
      // y = f | u  | = u  + e       * | e       - 1.0 - chi*u  - ----------- - ----------- |
      //        \  /                   |                              2.0           6.0     |
      //                               |                                                    |
      //                                \                                                  /
      */
      return(yp-(up+exp(-chi*B)*(exp(chi*up)-1.0-chi*up*(1+chi*up/2.0*(1.0+chi*up/3.0)))));
    }

  //-----------------------------------------------------------------
  //     evaluate the residual of Spaldings law of the wall
  //                                           (private) gammi 11/09
  //-----------------------------------------------------------------
  double JacobianSpaldingResidual_utau(
    const double y     ,
    const double visc  ,
    const double utau  ,
    const double normu )
    {
      // get dimensionless velocity
      const double up = normu/utau;
      
      const double chi=0.4;
      const double B  =5.5;

      //                    +
      // get derivative of y , a dimensionless boundary layer thickness
      const double dyplus_dutau=y/visc;

      //                                   +
      // get derivative of function f wrt u
      const double df_duplus=1+exp(-chi*B)*(chi*exp(chi*up)-chi-chi*chi*up-chi*chi*chi*up*up/2.0);

      //                    +
      // get derivative of u  wrt u_tau
      const double duplus_dutau=-normu/(utau*utau);
      
      return(dyplus_dutau-df_duplus*duplus_dutau);
    }

  //-----------------------------------------------------------------
  //     evaluate the residual of Spaldings law of the wall
  //                                           (private) gammi 11/09
  //-----------------------------------------------------------------
  double JacobianSpaldingResidual_u(
    const double y     ,
    const double visc  ,
    const double utau  ,
    const double normu )
    {
      // get dimensionless velocity
      const double up = normu/utau;

      const double chi=0.4;
      const double B  =5.5;

      //                    +
      // get derivative of y , a dimensionless boundary layer thickness
      const double dyplus_du=0.0;

      //                                   +
      // get derivative of function f wrt u
      const double df_duplus=1+exp(-chi*B)*(chi*exp(chi*up)-chi-chi*chi*up-chi*chi*chi*up*up/2.0);

      //                    +
      // get derivative of u  wrt u
      const double duplus_du=1.0/utau;
      
      return(dyplus_du-df_duplus*duplus_du);
    }



  //-----------------------------------------------------------------
  //     evaluate the residual of Spaldings law of the wall
  //                                           (private) gammi 11/09
  //-----------------------------------------------------------------
  double JacobianSpaldingResidual_uplus(
    const double y     ,
    const double visc  ,
    const double utau  ,
    const double normu )
    {
      // get dimensionless velocity
      const double up = normu/utau;

      const double chi=0.4;
      const double B  =5.5;

      //                    +
      // get derivative of y , a dimensionless boundary layer thickness
      const double dyplus_duplus=0.0;
      /*
      // Spaldings law of the wall
      //                     /                                                  \
      //                    |                           /     +\ 2    /     +\ 3 |
      //                    |       +                  | chi*u  |    | chi*u  |  |
      //  +   +    -chi*B   |  chi*u               +    \      /      \      /   |
      // y - u  - e       * | e       - 1.0 - chi*u  - ----------- - ----------- | = 0
      //                    |                              2.0           6.0     |
      //                    |                                                    |
      //                     \                                                  /
      */

      //                                   +
      // get derivative of function f wrt u
      const double df_duplus=1+exp(-chi*B)*(chi*exp(chi*up)-chi-chi*chi*up-chi*chi*chi*up*up/2.0);
      
      return(dyplus_duplus-df_duplus);
    }


  //! get material parameters
  void GetMaterialParams(Teuchos::RCP<const MAT::Material> material,  ///< reference pointer to material
                         const LINALG::Matrix<bdrynen_,1>& escaaf,    ///< scalar at time n+alpha_f / n+1
                         const double                      thermpressaf    ///< thermodynamic pressure at time n+alpha_f / n+1
                         );

  //! pointer to parameter list
  DRT::ELEMENTS::Fluid3ImplParameter* f3Parameter_;

  //! node coordinates for boundary element
  LINALG::Matrix<nsd_,bdrynen_> xyze_;
  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi_;
  //! array for shape functions for boundary element
  LINALG::Matrix<bdrynen_,1> funct_;
  //! array for shape function derivatives for boundary element
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv_;
  //! velocity vector at integration point
  LINALG::Matrix<nsd_,1> unitnormal_;
  //! velocity vector at integration point
  LINALG::Matrix<nsd_,1> velint_;
  //! infinitesimal area element drs
  double drs_;
  //! integration factor
  double fac_;
  //! physical viscosity
  double visc_;
  //! density at t_(n+alpha_F) or t_(n+1)
  double densaf_;
  //! density factor for Neumann boundary conditions (set according to problem type)
  double densfac_;

  };  // end class fluid3boundaryImpl

}
}

#endif

#endif // CCADISCRET
#endif // D_FLUID3

