/*!-------------------------------------------------------------------------*\
 * \file meshfree_scatra_impl.H
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

/*--------------------------------------------------------------------------*/
/*!
 * \brief Interface base class for ScaTraImpl
 *
 * This class exists to provide a common interface for all template
 * versions of ScaTraImpl. The only function
 * this class actually defines is Impl, which returns a pointer to
 * the appropriate version of ScaTraImpl.
 * Deriving from DRT::SingletonDestruction provides proper registration
 * of singleton at DRT::Problem instance and - as a result -
 * correct memory handling at the end of a BACI run.
 * (correctly delete allocated ScaTraImplInterface instances)
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeScaTraImplInterface : public DRT::SingletonDestruction
{

public:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/
  /*------------------------------------------------------------------------*/
  /*!
   * \brief Empty constructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeScaTraImplInterface() {}

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Empty destructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeScaTraImplInterface() {}


  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate the cell - virtual function
   *
   * This class does not provide a definition for this function, it must be
   * defined in ScaTraImpl.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int Evaluate(DRT::Element*              ele,
                       Teuchos::ParameterList&    params,
                       DRT::Discretization&       discretization,
                       std::vector<int>&          lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Internal implementation for meshfree scatra cells
   *
   * This function instantiates a meshfree scatra implementation, which then
   * evaluates the entire meshfree cell.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  static MeshfreeScaTraImplInterface* Impl(
      const DRT::Element* ele,
      const enum INPAR::SCATRA::ScaTraType& scatratype);

};

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
class MeshfreeScaTraImpl: public MeshfreeScaTraImplInterface
{
public:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  MeshfreeScaTraImpl(const int numdofpernode, const int numscal = 1);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Singleton access method
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  static MeshfreeScaTraImpl<distype> * Instance(const int numdofpernode, const int numscal = 1, bool create=true );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Called upon singleton destruction
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void Done();

  /*========================================================================*/
  //! @name Query method
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual int Evaluate(DRT::Element*              ele,            //!< the element to which sysmat belongs
                       Teuchos::ParameterList&    params,         //!< ?
                       DRT::Discretization&       discretization, //!< ?
                       std::vector<int>&          lm,             //!< ?
                       Epetra_SerialDenseMatrix&  elemat1_epetra, //!< ?
                       Epetra_SerialDenseMatrix&  elemat2_epetra, //!< ?
                       Epetra_SerialDenseVector&  elevec1_epetra, //!< ?
                       Epetra_SerialDenseVector&  elevec2_epetra, //!< ?
                       Epetra_SerialDenseVector&  elevec3_epetra  //!< ?
    );

private:

  /*========================================================================*/
  //! @name Auxiliary method - private
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Calculate matrix and rhs. Here the whole thing is hidden.
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void Sysmat(
    const DRT::MESHFREE::Cell*                  ele,       //!< the element to which sysmat belongs
    Epetra_SerialDenseMatrix&             sys_mat,   //!< element matrix to calculate
    Epetra_SerialDenseVector&             residual,  //!< element rhs to calculate
    const double &                        time,      //!< current simulation time
    const double &                        dt,        //!< current time-step length
    const double &                        timefac,   //!< time discretization factor
    const double &                        alphaF,    //!< factor for generalized-alpha time integration
    const enum INPAR::SCATRA::ScaTraType& scatratype //!< type of scalar transport problem
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief get the body force
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void BodyForce(
    const DRT::MESHFREE::Cell*   ele, //!< the element we are dealing with
    const double&                 time //!< current simulation time
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief get the material parameters
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void GetMaterialParams(
    const DRT::MESHFREE::Cell*            ele,       //!< the element we are dealing with
    const enum INPAR::SCATRA::ScaTraType& scatratype //!< type of scalar transport problem
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate element matrix and rhs vector
   *        (One-Step-Theta, BDF2, Gen.-alpha)
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void CalMatAndRHS(
    Epetra_SerialDenseMatrix&             emat,     //!< element stiffness matrix
    Epetra_SerialDenseVector&             erhs,     //!< element right-hand-side vector
    const double &                        fac,      //!< integration factor (Gauss weight including Jacobian)
    const double &                        timefac,  //!< integrator specific time factor
    const double &                        dt,       //!< time step length
    const double &                        alphaF,   //!< time step length
    const int &                           dofindex  //!< index of scalar (for multiple scalar transport)
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate mass matrix and rhs for determining initial time
   *        derivative
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void InitialTimeDerivative(
    const DRT::MESHFREE::Cell*            ele,       //!< the element
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    const enum INPAR::SCATRA::ScaTraType& scatratype //!< type of scalar transport problem
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate mass matrix and rhs for determining time derivative
   *        after reinitialisation
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void TimeDerivativeReinit(
      const DRT::MESHFREE::Cell*            ele,       //!< the element
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
      const double&                         dt,        //!< current time-step length
      const double&                         timefac,   //!< time discretization factor
      const enum INPAR::SCATRA::ScaTraType& scatratype //!< type of scalar transport problem
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate flux vector for actual solution
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void CalculateFlux(
      LINALG::SerialDenseMatrix&      flux,
      const DRT::MESHFREE::Cell*      ele,
      const INPAR::SCATRA::FluxType&  fluxtype,
      const int&                      dofindex,
      const enum INPAR::SCATRA::ScaTraType& scatratype //!< type of scalar transport problem
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate scalars and domain integral
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void CalculateScalars(
    const DRT::MESHFREE::Cell* ele,      //!< pointer to the element
    const std::vector<double>& ephinp,   //!< current scalar field
    Epetra_SerialDenseVector&  scalars,  //!< vector for computed integrals
    const bool                 inverting //!< flag for invertuing
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief calculate domain and bodyforce integral
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void CalculateDomainAndBodyforce(
    Epetra_SerialDenseVector&  scalars, //!< vector for the computed integrals
    const DRT::MESHFREE::Cell* ele,     //!< pointer to the element
    const double&              time     //!< current simulation time
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief integral of shape functions over the element
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  void IntegrateShapeFunctions(
    const DRT::MESHFREE::Cell*         ele,    //!< the current element
    Epetra_SerialDenseVector&          elevec1,//!< result vector (to be assembled)
    const Epetra_IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
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

  /*------------------------------------------------------------------------*/
  /*!
   * \brief  Do a finite difference check for a given element id.
   *         Meant for debugging only!
   *///                                                  (private) nis Mar12
  /*------------------------------------------------------------------------*/
  void FDcheck(
    const DRT::MESHFREE::Cell*            ele,       //!< the element those matrix is calculated
    Epetra_SerialDenseMatrix&             sys_mat,   //!< element matrix to calculate
    Epetra_SerialDenseVector&             residual,  //!< element rhs to calculate
    const double&                         time,      //!< current simulation time
    const double&                         dt,        //!< current time-step length
    const double&                         timefac,   //!< time discretization factor
    const double&                         alphaF,    //!< factor for generalized-alpha time integration
    const enum INPAR::SCATRA::ScaTraType& scatratype //!< type of scalar transport problem
  );

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

  //! transported scalar in gausspoint
  double phi_;
  //! vector containing all values from previous timelevel n in gausspoint
  std::vector<double> hist_;
  //! (global) gradient of transported scalar at integration point
  LINALG::SerialDenseVector gradphi_;
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

  //! velocity vector in gausspoint
  LINALG::SerialDenseVector velint_;
  //! convective velocity vector in gausspoint
  LINALG::SerialDenseVector convelint_;
  //! velocity divergence in gausspoint
  double vdiv_;

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
  //! rhs force in Gauss point
  std::vector<double> rhs_;
  //! rhs of temperature equation in case of reactive equation system
  std::vector<double> reatemprhs_;
  //! residual of scalar transport equation
  std::vector<double> scatrares_;
  //! convective term (part of residual)
  std::vector<double> conv_phi_;
  //! diffusive term (part of residual)
  std::vector<double> diff_phi_;
  //! reactive term (part of residual)
  std::vector<double> rea_phi_;
  //! (global) Laplacian of shape functions at integration point: (nen_)
  LINALG::SerialDenseVector laplace_;
  //! reconstructed flux values via L2 projection: (nsd_ x nen_)
  std::vector<LINALG::SerialDenseMatrix> efluxreconstr_;

  /*========================================================================*/
  //! @name Coefficients and related
  /*========================================================================*/

  //! density at t_(n) - (vector over scalars)
  std::vector<double> densn_;
  //! density at t_(n+1) or t_(n+alpha_F) - (vector over scalars)
  std::vector<double> densnp_;
  //! density at t_(n+alpha_M) - (vector over scalars)
  std::vector<double> densam_;
  //! factor for density derivative (at t_(n+1) or t_(n+alpha_F)) - (vector over scalars)
  std::vector<double> densgradfac_;
  //! diffusivity / diffusivities (in case of systems) / thermal conductivity - (vector over scalars)
  std::vector<double> diffus_;
  //! reaction coefficient - (vector over scalars)
  std::vector<double> reacoeff_;
  //! reaction coefficent derivative (needed for matrix contribution)
  std::vector<double> reacoeffderiv_;
  //! specific heat capacity (either at constant pressure or at constant volume)
  double shc_;
  //! fluid viscosity
  double visc_; // not initialized by constructor
  //! convective term: (nen_)
  LINALG::SerialDenseVector conv_;
  //! thermodynamic pressure at t_(n+1) or t_(n+alpha_F) (LOMA specific)
  double thermpressnp_;
  //! thermodynamic pressure at t_(n+alpha_M) (LOMA specific)
  double thermpressam_;
  //! time derivative of thermodynamic pressure (LOMA specific)
  double thermpressdt_;

  /*========================================================================*/
  //! @name Flags
  /*========================================================================*/

  //! flag for using ALE formulation
  bool is_ale_;
  //! flag for using ELCH formulation
  bool is_elch_;
  //! flag for inclusion of reactive terms
  bool is_reactive_;  // set to false by constructor - modified later
  //! flag for stationary algorithm
  bool is_stationary_;
  //! flag for generalized-alpha time integration
  bool is_genalpha_;
  //! flag for incremental solver
  bool is_incremental_;
  //! flag for conservative form of convective term
  bool is_conservative_;
  //! flag for material evaluation at Gaussian integration points
  bool mat_gp_;

}; // end class MeshfreeScaTraImpl

} // end namespace MESHFREE
} // ens namespace DRT

#endif // end MESHFREE_SCATRA_IMPL_H
