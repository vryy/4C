/*----------------------------------------------------------------------*/
/*!
\file meshfree_scatra_boundary_impl.H

\brief Internal implementation of meshfree scalar transport boundary cells

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MESHFREE_SCATRA_BOUNDARY_IMPL_H
#define MESHFREE_SCATRA_BOUNDARY_IMPL_H


#include <vector>
#include "../linalg/linalg_serialdensematrix.H"  // class variables
#include "../linalg/linalg_serialdensevector.H"  // class variables

#include "../drt_lib/drt_element.H"              // Element::DiscretizationType
#include "../drt_inpar/inpar_scatra.H"           // enums ScaTraType and FluxType
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // DRT::UTILS::DisTypeToDim

namespace DRT
{
namespace ELEMENTS
{

// forward declarations
class MeshfreeTransportBoundary;

/// Interface base class for ScaTraImpl
/*!
  This class exists to provide a common interface for all template
  versions of ScaTraImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of ScaTraImpl.
 */
class MeshfreeScaTraBoundaryImplInterface
{

public:
  /// Empty constructor
  MeshfreeScaTraBoundaryImplInterface() {}
  /// Empty destructor
  virtual ~MeshfreeScaTraBoundaryImplInterface() {}
  /// Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in ScaTraBoundaryImpl.
   */
  virtual int Evaluate(DRT::ELEMENTS::MeshfreeTransportBoundary* ele,
                       Teuchos::ParameterList&    params,
                       DRT::Discretization&       discretization,
                       std::vector<int>&          lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  /// Evaluate a Neumann boundary condition
  /*!
    This class does not provide a definition for this function, it
    must be defined in ScaTraBoundaryImpl.
   */
  virtual int EvaluateNeumann(DRT::ELEMENTS::MeshfreeTransportBoundary* ele,
                              Teuchos::ParameterList&           params,
                              DRT::Discretization&              discretization,
                              DRT::Condition&                   condition,
                              std::vector<int>&                 lm,
                              Epetra_SerialDenseVector&         elevec1_epetra
                              ) = 0;

  /// Internal implementation class for scalar transport elements
  static MeshfreeScaTraBoundaryImplInterface* Impl(
      const DRT::Element* ele,
      const enum INPAR::SCATRA::ScaTraType scatratype);

};

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
template<DRT::Element::DiscretizationType distype>
class MeshfreeScaTraBoundaryImpl: public MeshfreeScaTraBoundaryImplInterface
{
public:
  /// Constructor
  MeshfreeScaTraBoundaryImpl(int numdofpernode, int numscal);

  /// Singleton access method
  static MeshfreeScaTraBoundaryImpl<distype> * Instance(const int numdofpernode,const int numscal=1);

  //! number of element knots (nomenclature: T. Hughes, The finite element method)
  static const int nek_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of boundary(!) space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! Evaluate
  virtual int Evaluate(DRT::ELEMENTS::MeshfreeTransportBoundary* ele,
                       Teuchos::ParameterList&    params,
                       DRT::Discretization&       discretization,
                       std::vector<int>&          lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       );

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(DRT::ELEMENTS::MeshfreeTransportBoundary* ele,
                              Teuchos::ParameterList&           params,
                              DRT::Discretization&              discretization,
                              DRT::Condition&                   condition,
                              std::vector<int>&                 lm,
                              Epetra_SerialDenseVector&         elevec1_epetra
                              );

  private:

  /*========================================================================*/
  //! @name dofs and nodes
  /*========================================================================*/

   //! number of dof per node
  const int numdofpernode_;
  //! number of transported scalars (numscal_ <= numdofpernode_)
  const int numscal_;
  //! numer of element nodes
  int nen_;

  /*========================================================================*/
  //! @name velocity, pressure, and related
  /*========================================================================*/

  //! velocity vector in gausspoint
  LINALG::SerialDenseVector velint_;
  //! nodal displacement values for ALE: size(nsd_ x nen_)
  LINALG::SerialDenseMatrix edispnp_;

  /*========================================================================*/
  //! @name Galerkin approximation and related
  /*========================================================================*/

  //! array for shape functions: (nen_)
  LINALG::SerialDenseVector funct_;
  //! global derivatives of shape functions w.r.t x,y,z: (nsd_ x nen_)
  LINALG::SerialDenseMatrix deriv_;


  /*========================================================================*/
  //! @name coefficients and related
  /*========================================================================*/
  //! diffusivity / diffusivities (in case of systems) / thermal conductivity
  std::vector<double> diffus_;

  /*========================================================================*/
  //! @name flags
  /*========================================================================*/

  //! flag for using ALE formulation
  bool isale_;
  //! flag for stationary algorithm
  bool is_stationary_;
  //! flag for generalized-alpha time integration
  bool is_genalpha_;
  //! flag for incremental solver
  bool is_incremental_;

  };

}
}

#endif
