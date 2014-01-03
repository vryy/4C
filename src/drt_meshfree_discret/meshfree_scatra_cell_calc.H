/*!-------------------------------------------------------------------------*\
 * \file meshfree_scatra_cell_calc.H
 *
 * \brief Internal implementation of meshfree scalar transport problems
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#ifndef MESHFREE_SCATRA_IMPL_H
#define MESHFREE_SCATRA_IMPL_H


#include "meshfree_scatra_cell_interface.H"

#include <vector>
#include "drt_meshfree_cell_utils.H"

#include "../linalg/linalg_serialdensematrix.H"  // class variables
#include "../linalg/linalg_serialdensevector.H"  // class variables

#include "../drt_lib/drt_element.H"              // Element::DiscretizationType
#include "../drt_lib/drt_singletondestruction.H" // inheritance
#include "../drt_inpar/inpar_scatra.H"           // enums ScaTraType and FluxType
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // DRT::UTILS::DisTypeToDim

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;

namespace MESHFREE
{
  class MeshfreeDiscretization;
  class Cell;
}

namespace ELEMENTS
{
  class ScaTraEleParameter;
  class ScaTraEleParameterTimInt;

/*--------------------------------------------------------------------------*/
/*!
 * \brief Internal meshfree scalar transport cell implementation for
 *        backgroundmesh integration
 *
 * This internal class keeps all the working arrays needed to
 * calculate the transport element. Additionally the method Sysmat()
 * provides a clean and fast element implementation.
 *
 * <h3>Purpose</h3>
 *
 * The idea is to separate the element maintenance (class Transport)
 * from the mathematical contents (this class). Of course there are
 * different implementations of the Transport element, this is just one
 * such implementation.
 *
 * The Transport element will allocate exactly one object of this class
 * for all transport elements with the same number of nodes in the mesh.
 * This allows us to use exactly matching working arrays (and keep them
 * around.)
 *
 * The code is meant to be as clean as possible. This is the only way
 * to keep it fast. The number of working arrays has to be reduced to
 * a minimum so that the element fits into the cache. (There might be
 * room for improvements.)
 *
 * <h3>History</h3>
 *
 * The implementation here is the standard convection-diffusion element
 * capable of dealing with systems of transported scalars.
 *
 * Right now we do not read any stabilization parameters from the
 * input file but have a fixed version.
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/

template<DRT::Element::DiscretizationType distype>
class MeshfreeScaTraCellCalc: public MeshfreeScaTraCellInterface
{
  friend class ScaTraEleParameter;
  friend class ScaTraEleParameterTimInt;

public:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Destructor
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeScaTraCellCalc() {};

  /*========================================================================*/
  //! @name Query method
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual int Evaluate(
    DRT::ELEMENTS::MeshfreeTransport* ele,            //!< the element to which sysmat belongs
    Teuchos::ParameterList&           params,         //!< ?
    DRT::Discretization&              discretization, //!< ?
    const std::vector<int>&           lm,             //!< ?
    Epetra_SerialDenseMatrix&         elemat1_epetra, //!< ?
    Epetra_SerialDenseMatrix&         elemat2_epetra, //!< ?
    Epetra_SerialDenseVector&         elevec1_epetra, //!< ?
    Epetra_SerialDenseVector&         elevec2_epetra, //!< ?
    Epetra_SerialDenseVector&         elevec3_epetra  //!< ?
    );

  virtual int EvaluateService(
    DRT::ELEMENTS::MeshfreeTransport* ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization &             discretization,
    const std::vector<int> &          lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra){return 0;};

protected:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor, private since we are Singleton
   *///                                                  (private) nis Jan14
  /*------------------------------------------------------------------------*/
  MeshfreeScaTraCellCalc(const int numdofpernode, const int numscal);

  /*========================================================================*/
  //! @name Auxiliary method - private
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Calculate matrix and rhs. Here the whole thing is hidden.
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void Sysmat(
    const DRT::MESHFREE::Cell*            cell,    //!< the element to which sysmat belongs
    Epetra_SerialDenseMatrix&             sys_mat, //!< element matrix to calculate
    Epetra_SerialDenseVector&             residual
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief get the body force
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void BodyForce(
    const DRT::MESHFREE::Cell* cell, //!< the element we are dealing with
    const double               time //!< current simulation time
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief get the material parameters
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void GetMaterialParams(
    const DRT::MESHFREE::Cell*            ele       //!< the element we are dealing with
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate the Laplacian (weak form) - inline
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  void GetLaplacianWeakForm(
      double& val,
      const LINALG::SerialDenseMatrix& deriv, // (nsd_ x nen_)
      const int vi,
      const int ui
    )
  {
    val = 0.0;
    for (int j = 0; j<nsd_; j++)
      val += deriv(j, vi)*deriv(j, ui);
    return;
  };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate the RHS of Laplacian (weak form) - inline
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  void GetLaplacianWeakFormRHS(
      double& val,
      const LINALG::SerialDenseMatrix& deriv, // (nsd_ x nen_)
      const LINALG::SerialDenseVector& gradphi,
      const int vi
    )
  {
    val = 0.0;
    for (int j = 0; j<nsd_; j++)
      val += deriv(j,vi)*gradphi(j);
    return;
  };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate divergence of vector field (e.g., velocity) - inline
   *
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  void GetDivergence(double&                          vdiv, //!< (o) divergence
                     const LINALG::SerialDenseMatrix& evel, //!< (i) element velocities
                     const LINALG::SerialDenseMatrix& deriv //!< (i) derivatives of basis functions
    )
  {
    LINALG::SerialDenseMatrix vderiv(nsd_,nsd_);
    vderiv.Multiply('N','T',1.0,evel,deriv,0.0);

    vdiv = 0.0;
    // compute vel x,x  + vel y,y +  vel z,z at integration pointer
    for (int j = 0; j<nsd_; ++j)
      vdiv += vderiv(j,j);
    return;
  };

  /*========================================================================*/
  //! @name Static member variables
  /*========================================================================*/

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
  //! number of element knots (geometric description of integration cell)
  static const int nek_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  //! number of Gauss points
  static const int ngp_ = DRT::MESHFREE::DisTypeToNumGaussPoints<distype>::ngp;

  /*========================================================================*/
  //! @name parameter lists
  /*========================================================================*/

  //! pointer to parameter lists
  DRT::ELEMENTS::ScaTraEleParameter* scatrapara_;
  //! pointer to parameter list for time integration
  DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparatimint_;

  /*========================================================================*/
  //! @name Dofs and nodes
  /*========================================================================*/

  //! number of dof per node
  const int numdofpernode_;
  //! number of transported scalars (numscal_ <= numdofpernode_)
  const int numscal_;
  //! numer of element nodes
  int nen_;
  //! reference to discretisation
  DRT::MESHFREE::MeshfreeDiscretization* discret_;
  //! knot coordinates for boundary cell
  LINALG::SerialDenseMatrix kxyz_;
  //! Gauss point coordinates for boundary cell
  LINALG::SerialDenseMatrix gxyz_;
  //! Gauss point weights for boundary cell
  LINALG::SerialDenseVector gw_;

  /*========================================================================*/
  //! @name Scalar degrees of freedom and related
  /*========================================================================*/

  //! scalar at t_(n-1): size(nen_)
  std::vector<LINALG::SerialDenseVector> ephin_;
  //! scalar at t_(n+1) or t_(n+alpha_F): size(nen_)
  std::vector<LINALG::SerialDenseVector> ephinp_;
  //! scalar at t_(n+alpha_M): size(nen_)
  std::vector<LINALG::SerialDenseVector> ephiam_;
  //! history vector of transported scalars: size(nen_)
  std::vector<LINALG::SerialDenseVector> ehist_;

  /*========================================================================*/
  //! @name Velocity, pressure, and related
  /*========================================================================*/

  //! nodal velocity values at t_(n+1) or t_(n+alpha_F): size(nsd_ x nen_)
  LINALG::SerialDenseMatrix evelnp_;
  //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F): size(nsd_ x nen_)
  LINALG::SerialDenseMatrix econvelnp_;
  //! nodal acceleration values at t_(n+1) or t_(n+alpha_F): size(nsd_ x nen_)
  LINALG::SerialDenseMatrix eaccnp_;
  //! nodal pressure values at t_(n+1) or t_(n+alpha_F): size(nen_)
  LINALG::SerialDenseVector eprenp_;
  //! nodal displacement values for ALE: size(nsd_ x nen_)
  LINALG::SerialDenseMatrix edispnp_;


  /*========================================================================*/
  //! @name Galerkin approximation and related
  /*========================================================================*/

  //! array for shape functions: (nen_)
  LINALG::SerialDenseVector funct_;
  //! global derivatives of shape functions w.r.t x,y,z: (nsd_ x nen_)
  LINALG::SerialDenseMatrix deriv_;
  //! bodyforce of node dofs: (numdofpernode_)
  std::vector<LINALG::SerialDenseVector> bodyforce_;

  /*========================================================================*/
  //! @name Coefficients and related
  /*========================================================================*/

  //! density at t_(n) - (vector over scalars)
  std::vector<double> densn_;
  //! density at t_(n+1) or t_(n+alpha_F) - (vector over scalars)
  std::vector<double> densnp_;
  //! density at t_(n+alpha_M) - (vector over scalars)
  std::vector<double> densam_;
  //! diffusivity / diffusivities (in case of systems) / thermal conductivity - (vector over scalars)
  std::vector<double> diffus_;
  //! flag for reactivity
  bool is_reactive_;
  //! reaction coefficient - (vector over scalars)
  std::vector<double> reacoeff_;

}; // end class MeshfreeScaTraCellCalc

} // end namespace MESHFREE
} // ens namespace DRT

#endif // end MESHFREE_SCATRA_IMPL_H
