/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc.H

\brief evaluation of scatra boundary terms at integration points

\level 1

<pre>
\maintainer Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_H
#define SCATRA_ELE_BOUNDARY_CALC_H

#include "scatra_ele_action.H"
#include "scatra_ele_boundary_interface.H"
#include "scatra_ele_calc_utils.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_inpar/inpar_elch.H"

#include "../drt_fluid_ele/fluid_ele.H"

namespace FLD
{
  template<DRT::Element::DiscretizationType distype,int numdofpernode, DRT::ELEMENTS::Fluid::EnrichmentType enrtype> class RotationallySymmetricPeriodicBC;
}

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
class ScaTraEleParameterStd;
class ScaTraEleParameterTimInt;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleBoundaryCalc: public ScaTraBoundaryInterface
{
public:

  //! Constructor
  ScaTraEleBoundaryCalc(const int numdofpernode, const int numscal, const std::string& disname);

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

  //! setup element evaluation
  int SetupCalc( DRT::FaceElement*                 ele,
                 Teuchos::ParameterList&           params,
                 DRT::Discretization&              discretization);

  //! extract element based or nodal values
  //  return extracted values of phinp
  virtual void ExtractElementAndNodeValues(
    DRT::FaceElement*             ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
  );

  //! Evaluate the element (using location array)
  virtual int Evaluate(
      DRT::FaceElement*             ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      Epetra_SerialDenseMatrix&     elemat1_epetra,
      Epetra_SerialDenseMatrix&     elemat2_epetra,
      Epetra_SerialDenseVector&     elevec1_epetra,
      Epetra_SerialDenseVector&     elevec2_epetra,
      Epetra_SerialDenseVector&     elevec3_epetra
      );

  //! evaluate action
  virtual int EvaluateAction(
      DRT::FaceElement*                   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      SCATRA::BoundaryAction              action,
      DRT::Element::LocationArray&        la,
      Epetra_SerialDenseMatrix&           elemat1_epetra,
      Epetra_SerialDenseMatrix&           elemat2_epetra,
      Epetra_SerialDenseVector&           elevec1_epetra,
      Epetra_SerialDenseVector&           elevec2_epetra,
      Epetra_SerialDenseVector&           elevec3_epetra
      );

  //! evaluate Neumann boundary condition
  virtual int EvaluateNeumann(
      DRT::FaceElement*                   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      DRT::Condition&                     condition,
      DRT::Element::LocationArray&        la,
      Epetra_SerialDenseVector&           elevec1,
      const double                        scalar
      );

  //! evaluate scatra-scatra interface coupling condition at integration point
  template<DRT::Element::DiscretizationType distype_master>
  static void EvaluateS2ICouplingAtIntegrationPoint(
      DRT::Condition&                               s2icondition,
      const std::vector<LINALG::Matrix<nen_,1> >&   eslavephinp,
      const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,1> >&   emasterphinp,
      const LINALG::Matrix<nen_,1>&                 funct_slave,
      const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,1>&   funct_master,
      const LINALG::Matrix<nen_,1>&                 test_slave,
      const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,1>&   test_master,
      const int                                     numscal,
      const double                                  timefacfac,
      const double                                  timefacrhsfac,
      Epetra_SerialDenseMatrix*                     k_ss,     ///< element matrix for slave side
      Epetra_SerialDenseMatrix*                     k_sm,    ///< element matrix for master side
      Epetra_SerialDenseMatrix*                     k_ms,     ///< element matrix for slave side
      Epetra_SerialDenseMatrix*                     k_mm,    ///< element matrix for master side
      Epetra_SerialDenseVector*                     r_s,    ///< element residual for slave side
      Epetra_SerialDenseVector*                     r_m    ///< element residual for slave side
      );

  //! computes mechanoresponsive scalar transport
  void CalcMechanotransduction(
      DRT::FaceElement*                 ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm,
      Epetra_SerialDenseVector&         elevec1_epetra
      );

  private:

  // add nodal displacements to point coordinates
  void UpdateNodeCoordinates()
  {
    xyze_ += edispnp_;
    return;
  };

  protected:

  //! calculate normal vectors
  void CalcNormalVectors(
      Teuchos::ParameterList&            params,
      DRT::FaceElement*                  ele
      );

  //! evaluate shape functions and derivatives at int. point
  double EvalShapeFuncAndIntFac(
      const DRT::UTILS::IntPointsAndWeights<nsd_>&   intpoints,         ///< integration points
      const int                                      iquad,             ///< id of current Gauss point
      LINALG::Matrix<1 + nsd_,1>*                    normalvec = NULL   ///< normal vector at Gauss point(optional)
  );

  //! evaluate scatra-scatra interface coupling condition
  virtual void EvaluateS2ICoupling(
      const DRT::FaceElement*        ele,              ///< current boundary element
      Teuchos::ParameterList&        params,           ///< parameter list
      DRT::Discretization&           discretization,   ///< discretization
      DRT::Element::LocationArray&   la,               ///< location array
      Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
      Epetra_SerialDenseMatrix&      emastermatrix,    ///< element matrix for master side
      Epetra_SerialDenseVector&      eslaveresidual    ///< element residual for slave side
  );

  //! calculate boundary integral, i.e., surface area of boundary element
  void CalcBoundaryIntegral(
      const DRT::FaceElement*     ele,     //!< the element we are dealing with
      Epetra_SerialDenseVector&   scalar   //!< result vector for scalar integral to be computed
      );


  //! evaluate Kedem-Katchalsky interface condition
  void EvaluateKedemKatchalsky(
      const DRT::FaceElement*                       ele,            ///< the actual boundary element
      const std::vector<LINALG::Matrix<nen_,1> >&   ephinp,         ///< scalar values at element nodes
      const std::vector<LINALG::Matrix<nen_,1> >&   ephibar,        ///< mean concentration values at element nodes
      const LINALG::Matrix<nen_,1>&                 epressure,      ///< pressure values at element nodes
      const LINALG::Matrix<nen_,1>&                 f_ewss,          ///< factor for WSS at element nodes
      Epetra_SerialDenseMatrix&                     emat,           ///< element-matrix
      Epetra_SerialDenseVector&                     erhs,           ///< element-rhs
      const std::vector<int>*                       onoff,          ///<flag for dofs to be considered by membrane equations of Kedem and Katchalsky
      const double                                  perm,           ///< surface/interface permeability coefficient
      const double                                  conductivity,   ///< hydraulic conductivity at interface
      const double                                  sigma           ///< Staverman filtration coefficient at interface
  );

  //! integral of shape functions over boundary surface
  void IntegrateShapeFunctions(
      const DRT::FaceElement*   ele,    ///< the actual boundary element
      Teuchos::ParameterList&   params, ///< the parameter list
      Epetra_SerialDenseVector& elevec1,///< result vector (to be assembled)
      const bool                addarea ///< flag for area calculation
  );

  //! evaluate scatra-scatra interface coupling condition
  void EvaluateSurfacePermeability(
      const DRT::FaceElement*        ele,              ///< current boundary element
      Teuchos::ParameterList&        params,           ///< parameter list
      DRT::Discretization&           discretization,   ///< discretization
      DRT::Element::LocationArray&   la,               ///< location array
      Epetra_SerialDenseMatrix&      elemat1,     ///< element matrix for slave side
      Epetra_SerialDenseVector&      elevec1    ///< element residual for slave side
  );

  //!computes the factor for the wall shear stress dependent interface flux
  void WSSinfluence(
      LINALG::Matrix<nen_,1>&            f_ewss,
      const LINALG::Matrix<nsd_+1,nen_>& ewss,          ///<vector of euclidean norm of shear stress tensor of each node of the element
      const bool                         wssonoff,      ///<flag if WSS has influence on interface permeability
      const std::vector<double>*         coeffs);       ///<coefficients of the log law to determine the influence of WSS on concentration flux

  //! Factor needed for the calculation of reference concentrations
  virtual double FacForRefConc(
      const int                      iquad,             ///< current boundary integration point
      const DRT::FaceElement*        bele,              ///< current boundary element
      Teuchos::ParameterList&        params,            ///< parameter list
      DRT::Discretization&           discretization,    ///< discretization
      DRT::Element::LocationArray&   la,                ///< location array
      Epetra_SerialDenseMatrix&      elemat1,           ///< element matrix for slave side
      Epetra_SerialDenseVector&      elevec1            ///< element residual for slave side
  ){return 1.0;};

  //! compute integral of convective mass/heat flux over boundary surface
  virtual std::vector<double> CalcConvectiveFlux(
      const DRT::FaceElement*                       ele,
      const std::vector<LINALG::Matrix<nen_,1> >&   ephinp,
      const LINALG::Matrix<nsd_+1,nen_>&            evelnp,
      Epetra_SerialDenseVector&                     erhs
      );

  //! Compute a constant normal vector for a boundary element
  void GetConstNormal(
      LINALG::Matrix<nsd_+1,1>&            normal, ///< the constant normal vector
      const LINALG::Matrix<nsd_+1,nen_>&   xyze    ///< element node coordinates
      );

  //! calculate potential Neumann inflow terms
  virtual void NeumannInflow(
      const DRT::FaceElement*                   ele,
      Teuchos::ParameterList&                   params,
      DRT::Discretization&                      discretization,
      DRT::Element::LocationArray&              la,
      Epetra_SerialDenseMatrix&                 emat,
      Epetra_SerialDenseVector&                 erhs
      );

  //! get density at integration point
  virtual double GetDensity(
      Teuchos::RCP<const MAT::Material>             material,
      const std::vector<LINALG::Matrix<nen_,1> >&   ephinp,
      const int                                     k
      );

  //! calculate boundary condition due to convective heat transfer
  void ConvectiveHeatTransfer(
      const DRT::FaceElement*                       ele,
      Teuchos::RCP<const MAT::Material>             material,
      const std::vector<LINALG::Matrix<nen_,1> >&   ephinp,
      Epetra_SerialDenseMatrix&                     emat,
      Epetra_SerialDenseVector&                     erhs,
      const double                                  heatranscoeff,
      const double                                  surtemp
      );

  /*!
  \brief Evaluate weak Dirichlet boundary conditions

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param material (in)      : material of this element
  \param elemat1 (out)      : matrix to be filled by element.
  \param elevec1 (out)      : vector to be filled by element.

  */
  template <DRT::Element::DiscretizationType bdistype,DRT::Element::DiscretizationType pdistype>
  void WeakDirichlet(DRT::FaceElement*                  ele,
                     Teuchos::ParameterList&            params,
                     DRT::Discretization&               discretization,
                     Teuchos::RCP<const MAT::Material>  material,
                     Epetra_SerialDenseMatrix&          elemat_epetra,
                     Epetra_SerialDenseVector&          elevec_epetra);


  //! calculate boundary conditions for impl. Characteristic Galerkin time integration, just for the reinitialization equation
  template <DRT::Element::DiscretizationType bdistype,DRT::Element::DiscretizationType pdistype>
  void ReinitCharacteristicGalerkinBoundary(
      DRT::FaceElement*                  ele,                  //!< transport element
      Teuchos::ParameterList&            params,               //!< parameter list
      DRT::Discretization&               discretization,       //!< discretization
      Teuchos::RCP<const MAT::Material>  material,             //!< material
      Epetra_SerialDenseMatrix&          elemat_epetra,        //!< ele sysmat
      Epetra_SerialDenseVector&          elevec_epetra         //!< ele rhs
      );

  //! evaluate Robin boundary condition
  virtual void CalcRobinBoundary(
      DRT::FaceElement*                 ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm,
      Epetra_SerialDenseMatrix&         elemat1_epetra,
      Epetra_SerialDenseVector&         elevec1_epetra,
      const double                      scalar
  );

  //! pointer to parameter list for time integration
  DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparamstimint_;
  //! pointer to parameter list
  DRT::ELEMENTS::ScaTraEleParameterStd* scatraparams_;

   //! number of dof per node
  const int numdofpernode_;
  //! number of transported scalars (numscal_ <= numdofpernode_)
  const int numscal_;

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
  //! for the handling of rotationally symmetric periodic boundary conditions
  Teuchos::RCP<FLD::RotationallySymmetricPeriodicBC<distype,nsd_+2,DRT::ELEMENTS::Fluid::none> > rotsymmpbc_;
};
}
}

#endif
