/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_std.H

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_H
#define SCATRA_ELE_BOUNDARY_CALC_H

#include "scatra_ele_boundary_interface.H"
#include "scatra_ele_calc_utils.H"
#include "scatra_ele_action.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_elch.H"


  /// Internal Scalar transport element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the transport element. Additionally the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Transport)
    from the mathematical contents (this class). Of course there are
    different implementations of the Transport element, this is just one
    such implementation.

    The Transport element will allocate exactly one object of this class
    for all transport elements with the same number of nodes in the mesh.
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
  */

namespace DRT
{

namespace ELEMENTS
{
class ScaTraEleParameter;
class ScaTraEleParameterTimInt;

template<DRT::Element::DiscretizationType distype>
class ScaTraBoundaryImpl: public ScaTraBoundaryInterface
{
public:

  //! Constructor
  ScaTraBoundaryImpl(const int numdofpernode, const int numscal);

  //! number of element nodes (nomenclature: T. Hughes, The finite element method)
  static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of boundary(!) space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  // schott
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  //! compute largest element diameter for reinitialization pseudo time step size
  template<DRT::Element::DiscretizationType pdistype, class M1>
  double getEleDiameter(const M1& xyze)
  {
    double elediam = 0.0;

    // number of nodes of this element
  //  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::;
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

    // check all possible connections between nodes of an element
    // node 1 to 2
    //    :
    // node 1 to 8 = numnode
    // node 2 to 3
    //    :
    // node 2 to 8
    //    :
    //    :
    //    :
    // node 7 to 8
    for(size_t i_start=0; i_start< numnode-2; ++i_start)
    {
      for(size_t i_end= i_start+1; i_end < numnode-1; ++i_end)
      {
        LINALG::Matrix<3,1> direction;
        direction.Clear();
        direction(0) = xyze(0, i_start) - xyze(0, i_end);
        direction(1) = xyze(1, i_start) - xyze(1, i_end);
        direction(2) = xyze(2, i_start) - xyze(2, i_end);

        // update elediam
        if (direction.Norm2() > elediam) elediam=direction.Norm2();
      }
    }

    return elediam;
  }

  //! Evaluate Boundary ActionSetup some basic variables needed by all methods
  int SetupCalc(DRT::ELEMENTS::TransportBoundary* ele,
                 Teuchos::ParameterList&           params,
                 DRT::Discretization&              discretization);

  //! Evaluate Boundary Action
  int EvaluateAction( DRT::ELEMENTS::TransportBoundary*  ele1,
                      Teuchos::ParameterList&         params,
                      DRT::Discretization&            discretization,
                      SCATRA::BoundaryAction          action,
                      std::vector<int>&               lm,
                      Epetra_SerialDenseMatrix&       elemat1,
                      Epetra_SerialDenseMatrix&       elemat2,
                      Epetra_SerialDenseVector&       elevec1,
                      Epetra_SerialDenseVector&       elevec2,
                      Epetra_SerialDenseVector&       elevec3);

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(DRT::ELEMENTS::TransportBoundary* ele,
                              Teuchos::ParameterList&           params,
                              DRT::Discretization&              discretization,
                              DRT::Condition&                   condition,
                              std::vector<int>&                 lm,
                              Epetra_SerialDenseVector&         elevec1_epetra
                              );


  //!computes the factor for the wall shear stress dependent interface flux
  std::vector<double> WSSinfluence(
      std::vector<double>  ewss,          ///<vector of euclidean norm of shear stress tensor of each node of the element
      const bool           wssonoff,      ///<flag if WSS has influence on interface permeability
      const std::vector<double>* coeffs); ///<coefficients of the log law to determine the influence of WSS on concentration flux

  protected:

  //! calculate normal vectors
  void CalcNormalVectors(
      Teuchos::ParameterList&            params,
      DRT::ELEMENTS::TransportBoundary*  ele
      );

  //! evaluate loma thermal press
  void CalcLomaThermPress(
      DRT::ELEMENTS::TransportBoundary* ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm
      );

  //! evaluate shape functions and derivatives at int. point
  double EvalShapeFuncAndIntFac(
      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
      const int                                    iquad,      ///< id of current Gauss point
      const int                                    eleid,      ///< the element id
      LINALG::Matrix<1 + nsd_,1>*  normalvec = NULL            ///< normal vector at Gauss point(optional)
  );

  //! evaluate scatra-scatra interface coupling condition
  void EvaluateS2ICoupling(
      const DRT::Element*         ele,              ///< current boundary element
      Teuchos::ParameterList&     params,           ///< parameter list
      DRT::Discretization&        discretization,   ///< discretization
      std::vector<int>&           lm,               ///< location vector
      Epetra_SerialDenseMatrix&   emat,             ///< element matrix
      Epetra_SerialDenseMatrix&   eauxmat,          ///< auxiliary element matrix (for coupling terms)
      Epetra_SerialDenseVector&   erhs              ///< element right-hand side
  );

    //! evaluate Kedem-Katchalsky interface condition
    void EvaluateKedemKatchalsky(
        const DRT::Element*        ele,         ///< the actual boundary element
        const std::vector<double>& ephinp,      ///< scalar values at element nodes
        const std::vector<double>& epressure,   ///< pressure values at element nodes
        const std::vector<double>& ephibar,     ///< mean concentration values at element nodes
        const std::vector<double>& f_wss,       ///< factor for WSS at element nodes
        Epetra_SerialDenseMatrix&  emat,        ///< element-matrix
        Epetra_SerialDenseVector&  erhs,        ///< element-rhs
        const double               timefac,     ///< time factor
        const std::vector<int>*    onoff,       ///<flag for dofs to be considered by membrane equations of Kedem and Katchalsky
        const double               perm,        ///< surface/interface permeability coefficient
        const double               conductivity,///< hydraulic conductivity at interface
        const double               sigma ///< Staverman filtration coefficient at interface
    );

     //! integral of shape functions over boundary surface
     void IntegrateShapeFunctions(
         const DRT::Element*       ele,    ///< the actual boundary element
         Teuchos::ParameterList&   params, ///< the parameter list
         Epetra_SerialDenseVector& elevec1,///< result vector (to be assembled)
         const bool                addarea ///< flag for area calculation
     );

  //! evaluate surface/interface permeability
  void EvaluateSurfacePermeability(
      const DRT::Element*        ele,    ///< current boundary element
      const std::vector<double>& ephinp, ///< scalar values at element nodes
      const std::vector<double>& f_wss,  ///< factor for WSS at element nodes
      Epetra_SerialDenseMatrix&  emat,   ///< element-matrix
      Epetra_SerialDenseVector&  erhs,   ///< element-rhs
      const double               timefac,///< time factor
      const std::vector<int>*    onoff,  ///<flag for dofs to be considered by membrane equations of Kedem and Katchalsky
      const double               perm    ///< surface/interface permeability coefficient
  );

  //! integral of normal diffusive flux and velocity over boundary surface
  void NormDiffFluxAndVelIntegral(
      const DRT::Element*             ele,
      Teuchos::ParameterList&         params,
      const std::vector<double>&      enormdiffflux,
      const std::vector<double>&      enormvel);

  //! compute integral of convective mass/heat flux over boundary surface
  std::vector<double> CalcConvectiveFlux(
      const DRT::Element*                 ele,
      const std::vector<double>&          ephinp,
      const LINALG::Matrix<nsd_+1,nen_>&  evelnp,
      Epetra_SerialDenseVector&           erhs);

  //! Compute a constant normal vector for a boundary element
  void GetConstNormal(
      LINALG::Matrix<nsd_+1,1>&            normal, ///< the constant normal vector
      const LINALG::Matrix<nsd_+1,nen_>&   xyze    ///< element node coordinates
      );

  //! calculate potential Neumann inflow terms
  void NeumannInflow(
      const DRT::Element*                 ele,
      Teuchos::RCP<const MAT::Material>   material,
      const std::vector<double>&          ephinp,
      const LINALG::Matrix<nsd_+1,nen_>&  evelnp,
      Epetra_SerialDenseMatrix&           emat,
      Epetra_SerialDenseVector&           erhs,
      const double                        timefac);

  //! calculate boundary condition due to convective heat transfer
  void ConvectiveHeatTransfer(
      const DRT::Element*                 ele,
      Teuchos::RCP<const MAT::Material>   material,
      const std::vector<double>&           ephinp,
      Epetra_SerialDenseMatrix&           emat,
      Epetra_SerialDenseVector&           erhs,
      const double                        heatranscoeff,
      const double                        surtemp,
      const double                        timefac);

  /*!
  \brief Evaluate weak Dirichlet boundary conditions

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param material (in)      : material of this element
  \param elemat1 (out)      : matrix to be filled by element.
  \param elevec1 (out)      : vector to be filled by element.

  */
  template <DRT::Element::DiscretizationType bdistype,DRT::Element::DiscretizationType pdistype>
  void WeakDirichlet(DRT::ELEMENTS::TransportBoundary*  ele,
                     Teuchos::ParameterList&            params,
                     DRT::Discretization&               discretization,
                     Teuchos::RCP<const MAT::Material>  material,
                     Epetra_SerialDenseMatrix&          elemat_epetra,
                     Epetra_SerialDenseVector&          elevec_epetra);


  //! calculate boundary conditions for impl. Characteristic Galerkin time integration, just for the reinitialization equation
  template <DRT::Element::DiscretizationType bdistype,DRT::Element::DiscretizationType pdistype>
  void ReinitCharacteristicGalerkinBoundary(
      DRT::ELEMENTS::TransportBoundary*  ele,                  //!< transport element
      Teuchos::ParameterList&            params,               //!< parameter list
      DRT::Discretization&               discretization,       //!< discretization
      Teuchos::RCP<const MAT::Material>  material,             //!< material
      Epetra_SerialDenseMatrix&          elemat_epetra,        //!< ele sysmat
      Epetra_SerialDenseVector&          elevec_epetra         //!< ele rhs
      );

  //! pointer to parameter list for time integration
  DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparamstimint_;
  //! pointer to parameter list
  DRT::ELEMENTS::ScaTraEleParameter* scatraparams_;

   //! number of dof per node
  const int numdofpernode_;
  //! number of transported scalars (numscal_ <= numdofpernode_)
  const int numscal_;

  //! flag for using ALE formulation
  bool isale_;
  //! flag for stationary algorithm
  bool is_stationary_;
  //! flag for generalized-alpha time integration
  bool is_genalpha_;
  //! flag for incremental solver
  bool is_incremental_;

  //! node coordinates
  LINALG::Matrix<nsd_+1,nen_> xyze_;
  //! weights for nurbs elements
  LINALG::Matrix<nen_,1> weights_;
  //! knot vector for nurbs elements
  std::vector<Epetra_SerialDenseVector> myknots_;
  //! knot vector of corresponding parent element
  std::vector<Epetra_SerialDenseVector> mypknots_;
  //! for nurbs elements the normal vector is multiplied with normalfac_!
  double normalfac_;
  //! nodal displacement values for ALE
  LINALG::Matrix<nsd_+1,nen_> edispnp_;
  //! diffusivity / diffusivities (in case of systems) / thermal conductivity
  std::vector<double> diffus_;
  //! charge numbers of all ion species (ELCH specific)
  // vector<double> valence_;
  //! specific heat capacity at constant pressure (in case of temperature eq.)
  double shcacp_;
  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<nsd_,1> xsi_;
  //! array for shape functions
  LINALG::Matrix<nen_,1> funct_;
  //! array for shape function derivatives w.r.t r,s,t
  LINALG::Matrix<nsd_,nen_> deriv_;
  //! global derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<nsd_,nen_> derxy_;
  //! unit normal vector at integration point
  LINALG::Matrix<nsd_+1,1> normal_;
  //! velocity vector in gausspoint
  LINALG::Matrix<nsd_+1,1> velint_;
  //! metric tensor at integration point
  LINALG::Matrix<nsd_,nsd_>  metrictensor_;

  //! thermodynamic pressure (LOMA specific)
  double thermpress_;
};
}
}

#endif
