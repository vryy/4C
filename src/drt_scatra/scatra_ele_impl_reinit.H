/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_impl_reinit.H

\brief Internal implementation of pde-based
       reinitialization of level-set methods

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_IMPL_REINIT_H
#define SCATRA_ELE_IMPL_REINIT_H

#include "scatra_ele_impl_interface.H"
#include "scatra_ele_impl_utils.H"
#include "scatra_ele_impl_reinit_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
namespace ELEMENTS
{

//! Internal Scalar transport element implementation
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
class ReInitImpl: public ScaTraImplInterface
{
public:
  //! Constructor
  ReInitImpl(const int numdofpernode, const int numscal);

  //! Singleton access method
  static ReInitImpl<distype> * Instance(const int numdofpernode, const int numscal, bool create=true );

  //! called upon destruction
  virtual void Done();

  //! Evaluate
  virtual int Evaluate(
    DRT::Element*              ele,            //!< the element we are dealing with
    Teuchos::ParameterList&    params,         //!< ?
    DRT::Discretization&       discretization, //!< the discretization the element is part of
    std::vector<int>&          lm,             //!< ?
    Epetra_SerialDenseMatrix&  elemat1_epetra, //!< ?
    Epetra_SerialDenseMatrix&  elemat2_epetra, //!< ?
    Epetra_SerialDenseVector&  elevec1_epetra, //!< ?
    Epetra_SerialDenseVector&  elevec2_epetra, //!< ?
    Epetra_SerialDenseVector&  elevec3_epetra  //!< ?
  );

  /*========================================================================*/
  //! @name static member variables
  /*========================================================================*/

  //! number of element nodes (nomenclature: T. Hughes, The finite element method)
  static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! element-type specific flag if second derivatives are needed
  static const bool use2ndderiv_ = SCATRA::Use2ndDerivs<distype>::use;

  //! element-type specific flag if second derivatives are needed
  static const bool use2ndderivReinitialization_ = SCATRA::Use2ndDerivsReinitialization<distype>::use;

  //! number of components necessary to store second derivatives
  // 1 component  for nsd=1:  (N,xx)
  // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
  // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
  static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;


private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element*                   ele,       //!< the element we are dealing with
    const enum INPAR::SCATRA::ScaTraType  scatratype //!< type of scalar transport problem
  );

  //! evaluate shape functions and their derivatives at current integration point
  virtual double EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, //!< integration points
    const int                                    iquad,     //!< id of current Gauss point
    const int                                    eleid      //!< the element id
  );

  //! calculate the Laplacian for all shape functions(strong form)
  inline void GetLaplacianStrongForm(
    LINALG::Matrix<nen_,1>&                diff,  //!< ?
    const LINALG::Matrix<numderiv2_,nen_>& deriv2 //!< ?
  );

  //! calculate the Laplacian (weak form)
  inline void GetLaplacianWeakForm(
    double&                          val,   //!< ?
    const LINALG::Matrix<nsd_,nen_>& derxy, //!< ?
    const int                        vi,    //!< ?
    const int                        ui     //!< ?
  );

  //! calculate the Laplacian (weak form)
  inline void GetLaplacianWeakFormRHS(
    double& val,                              //!< ?
    const LINALG::Matrix<nsd_,nen_>& derxy,   //!< shape function derivative in xyz
    const LINALG::Matrix<nsd_,1>&    gradphi, //!< ?
    const int vi                              //!< ?
  );

  //! calculate divergence of vector field (e.g., velocity)
  inline void GetDivergence(
    double&                          vdiv,
    const LINALG::Matrix<nsd_,nen_>& evel,
    const LINALG::Matrix<nsd_,nen_>& derxy
  );

  /*========================================================================*/
  //! @name stabilization
  /*========================================================================*/

  //! calculate stabilization parameter
  virtual void CalTau(
    DRT::Element*  ele,           //!< the element we are dealing with
    double         diffus,        //!< ?
    const double   dt,            //!< current time-step length
    const double   timefac,       //!< time discretization factor
    const double   vol,           //!< element volume
    const int      k,             //!< index of current scalar
    const double   frt,           //!< factor F/RT needed for ELCH calculations
    const bool     migrationintau //!< ?
  );

  //! calculate characteristic element length
  virtual double CalcCharEleLength(
    const double  vol,     //!< element volume
    const double  vel_norm //!?
  );



  /*========================================================================*/
  //! @name expl. Taylor-Galerkin and impl. Characteristic Galerkin time integration
  /*========================================================================*/

  //! Calculate matrix and rhs. Here the whole thing is hidden.
  virtual void Sysmat_Transport_TG(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    const double                          dt,        //!< current time-step length
    const enum INPAR::SCATRA::TimeIntegrationScheme timealgo //!< time algorithm
  );

  //! calculate element matrix and rhs vector (implicit characteristic Galerkin method 2nd order)
  virtual void CalMatAndRHS_Transport_ICG2(
    Epetra_SerialDenseMatrix&  emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector&  erhs,     //!< element rhs to calculate
    const double               fac,      //!< integration factor
    const int                  k,        //!< index of current scalar
    DRT::Element*              ele,      //!< the element we are dealing with
    const double               dt        //!< current time-step length
  );

  //! calculate element matrix and rhs vector (TaylorGalerkin3)
  virtual void CalMatAndRHS_Transport_TG3(
    Epetra_SerialDenseMatrix&  emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector&  erhs,     //!< element rhs to calculate
    const double               fac,      //!< integration factor
    const int                  k,        //!< index of current scalar
    DRT::Element*              ele,      //!< the element we are dealing with
    const double               dt        //!< current time-step length
  );


  /*========================================================================*/
  //! @name Reinitialization (integrate element matrices and rhs)
  /*========================================================================*/

  //! evaluate shape functions and their derivatives at current integration point with second order derivatives
  virtual double EvalShapeFuncAndDerivsAtIntPointREINIT(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, //!< integration points
    const int                                    iquad,     //!< id of current Gauss point
    const int                                    eleid      //!< the element id
  );

  //! calculate error for steady state check
  virtual void CalErrorsReinitialization(
    const DRT::Element*          ele,   //!< the element
    Teuchos::ParameterList&      params //!< parameter list
  );

  //! calculate element matrix and rhs vector (One-Step-Theta)
  virtual void CalMatAndRHS_REINIT_ICG2(
    Epetra_SerialDenseMatrix&             emat,                     //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                     //!< element rhs to calculate
    const double                          fac,                      //!< integration factor
    const int                             k,                        //!< index of current scalar
    DRT::Element*                         ele,                      //!< the element we are dealing with
    double                                pseudo_timestep_size,     //!< pseudo time step size
    double                                mesh_size,                //!< meshsize
    const INPAR::SCATRA::PenaltyMethod    penalty_method,           //!< penalty method to keep the interface the interface position
    const double                          penalty_interface_reinit, //!< penalty parameter to keep the interface position
    const double                          epsilon_bandwidth,        //!< epsilon bandwith for the smoothed sign function
    const INPAR::SCATRA::SmoothedSignType smoothedSignType          //!< type of smoothing the sign function
  );

  //! calculate element matrix and rhs vector for the reinitialization equation (One-Step-Theta)
  virtual void CalMatAndRHS_REINIT_OST(
    Epetra_SerialDenseMatrix&             emat,                        //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                        //!< element rhs to calculate
    const double                          fac,                         //!< integration factor
    const int                             k,                           //!< index of current scalar
    DRT::Element*                         ele,                         //!< the element we are dealing with
    double                                pseudo_timestep_size_factor, //!< pseudo time step factor
    double                                meshsize,                    //!< meshsize
    const INPAR::SCATRA::PenaltyMethod    penalty_method,              //!< penalty method to keep the interface the interface position
    const double                          penalty_interface_reinit,    //!< penalty parameter to keep the interface position
    const double                          epsilon_bandwidth,           //!< epsilon bandwith for the smoothed sign function
    const INPAR::SCATRA::SmoothedSignType smoothedSignType,            //!< type of smoothing the sign function
    const bool                            shock_capturing,             //!< shock capturing switched on/off
    const double                          shock_capturing_diffusivity, //!< shock capturing diffusivity
    double                                timefac                      //!< time discretization factor
  );

  //! calculate element matrix and rhs vector for the intersection penalty method
  virtual void CalMatAndRHS_REINIT_Penalty(
    Epetra_SerialDenseMatrix&             emat,                     //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                     //!< element rhs to calculate
    const int                             k,                        //!< index of current scalar
    DRT::Element*                         ele,                      //!< the element we are dealing with
    const double                          penalty_interface_reinit  //!< penalty parameter to keep the interface position
  );

  //! calculate system matrix and rhs for the first order reinitialization equation based on an implicit characteristic Galerkin scheme
  virtual void Sysmat_Reinit_TG(
      DRT::Element*                         ele,                               //!< the element those matrix is calculated
      Epetra_SerialDenseMatrix&             emat,                              //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs,                              //!< element rhs to calculate
      const bool                            reinitswitch,                      //!< flag if scatra object is a reinitializer
      const double                          reinit_pseudo_timestepsize_factor, //!< pseudo timestep factor
      const INPAR::SCATRA::SmoothedSignType smoothedSignType,                  //!< type of smoothed sign function
      const double                          epsilon_bandwidth,                 //!< epsilon bandwith
      const INPAR::SCATRA::PenaltyMethod    penalty_method,                    //!< penalty method to keep the interface position
      const double                          penalty_interface_reinit           //!< penalty parameter to keep the interface position
  );

  //! Calculate matrix and rhs for the reinitializaton equation based on standard OST time discretization
  virtual void Sysmat_Reinit_OST(
    DRT::Element*                         ele,           //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,          //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,          //!< element rhs to calculate
    const double                          dt,            //!< current time-step length
    const double                          timefac,       //!< time discretization factor
    const double                          meshsize,      //!< meshsize
    const bool                            reinitswitch,  //!< flag if scatra object is a reinitialization object
    const double                          reinit_pseudo_timestepsize_factor, //!< pseudo time step factor
    const INPAR::SCATRA::SmoothedSignType smoothedSignType,                  //!< type of smoothed sign function
    const INPAR::SCATRA::PenaltyMethod    penalty_method,                    //!< penalty method to keep the interface position
    const double                          penalty_interface_reinit,          //!< penalty parameter to keep the interface position
    const double                          epsilon_bandwidth,                 //!< epsilon bandwith
    const bool                            shock_capturing,                   //!< shockcapturing operator turned on/off
    const double                          shock_capturing_diffusivity,       //!< shockcapturing diffusivity
    const enum INPAR::SCATRA::ScaTraType  scatratype                         //!< type of scalar transport problem
  );

  /*========================================================================*/
  //! @name dofs and nodes
  /*========================================================================*/

  //! number of dof per node
  const int numdofpernode_;
  //! number of transported scalars (numscal_ <= numdofpernode_)
  const int numscal_;

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! flag for using ALE formulation
  bool is_ale_;
  //! flag for inclusion of reactive terms
  bool is_reactive_;
  //! (sign) factor for diffusive and reactive stabilization terms
  double diffreastafac_;

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
  //! flag for stabilization parameter evaluation at Gaussian integration points
  bool tau_gp_;
  //! flag for improvement of consistency of stabilized method
  bool betterconsistency_;

  //! stabilization parameter definition
  enum INPAR::SCATRA::TauType whichtau_;

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  // ? scalar
  // double phi_

  //! (global) gradient of transported scalar at integration point
  LINALG::Matrix<nsd_,1> gradphi_;

  //! scalar at t_(n)
  std::vector<LINALG::Matrix<nen_,1> > ephin_;
  //! scalar at t_(n+1) or t_(n+alpha_F)
  std::vector<LINALG::Matrix<nen_,1> > ephinp_;
  //! scalar at t_(n+alpha_M)
  std::vector<LINALG::Matrix<nen_,1> > ephiam_;

  //! vector containing all values from previous timelevel n in gausspoint
  std::vector<double> hist_;
  //! history vector of transported scalars
  std::vector<LINALG::Matrix<nen_,1> > ehist_;

  //! scalar at t_(0) used for reinitialization as reference
  std::vector<LINALG::Matrix<nen_,1> > ephi0_Reinit_Reference_;
  //! scalar at t_(0) used for reinitialization as reference
  std::vector<LINALG::Matrix<nen_,1> > ephi0_penalty_;

  //! fine-scale solution ?
  std::vector<LINALG::Matrix<nen_,1> > fsphinp_;

  //! gradient of phi at global time for reinitialization(initialization for pseudo time loop)
  LINALG::Matrix<nsd_,1> grad_phi_0;

  /*========================================================================*/
  //! @name velocity, pressure, and related
  /*========================================================================*/

  //! nodal velocity values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<nsd_,nen_> evelnp_;
  //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<nsd_,nen_> econvelnp_;
  //! nodal fine-scale velocity values at t_(n+1) or t_(n+alpha_F)
  //! required for fine-scale subgrid diffusivity of type smagorinsky_small and multifractal subgrid scales
  LINALG::Matrix<nsd_,nen_> efsvel_;
  //! nodal acceleration values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<nsd_,nen_> eaccnp_;
  //! nodal displacement values for ALE
  LINALG::Matrix<nsd_,nen_> edispnp_;

  //! velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> velint_;
  //! convective velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> convelint_;
  //! subgrid-scale velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> sgvelint_;
  //! fine-scale velocity values in gausspoint
  //! required for fine-scale subgrid diffusivity of type smagorinsky_small and multifractal subgrid scales
  LINALG::Matrix<nsd_,1> fsvelint_;
  //! multifractal subgrid-scale velocity vector in gausspoint
  LINALG::Matrix<nsd_,1> mfsgvelint_;
  //! migration velocity vector (divided by D_k*z_k) in gausspoint (ELCH)
  LINALG::Matrix<nsd_,1> migvelint_;

  //! convective term ?
  LINALG::Matrix<nen_,1> conv_;
  //! subgrid-scale convective term ?
  LINALG::Matrix<nen_,1> sgconv_;

  //! velocity divergence
  double vdiv_;

  //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<nen_,1> eprenp_;

  /*========================================================================*/
  //! @name coefficients and related
  /*========================================================================*/

  //! density at t_(n)
  std::vector<double> densn_;
  //! density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp_;
  //! density at t_(n+alpha_M)
  std::vector<double> densam_;
  //! factor for density derivative (at t_(n+1) or t_(n+alpha_F))
  std::vector<double> densgradfac_;

  //! diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  std::vector<double> diffus_;
  //! all-scale or fine-scale subgrid diffusivities (one per transported scalar)
  std::vector<double> sgdiff_;
  //! reaction coefficient
  std::vector<double> reacoeff_;
  //! reaction coefficent derivative (needed for matrix contribution)
  std::vector<double> reacoeffderiv_;

  //! specific heat capacity (either at constant pressure or at constant volume)
  double shc_;
  //! fluid viscosity
  double visc_;

  //! diffusive part used in stabilization terms
  LINALG::Matrix<nen_,1> diff_;
  //! convective part of migration term (divided by D_k*z_k) used in stabilization terms (ELCH)
  LINALG::Matrix<nen_,1> migconv_;
  // reactive part of migration term used in stabilization terms
  LINALG::Matrix<nen_,1> migrea_;

  /*========================================================================*/
  //! @name Galerkin approximation and related
  /*========================================================================*/

  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<nsd_,1> xsi_;
  //! node coordinates
  LINALG::Matrix<nsd_,nen_> xyze_;

  //! array for shape functions
  LINALG::Matrix<nen_,1> funct_;
  //! array for shape function derivatives w.r.t r,s,t
  LINALG::Matrix<nsd_,nen_> deriv_;
  //! array for second derivatives of shape function w.r.t r,s,t
  LINALG::Matrix<numderiv2_,nen_> deriv2_;
  //! global derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<nsd_,nen_> derxy_;
  //! global second derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<numderiv2_,nen_> derxy2_;

  //! transposed jacobian "dx/ds"
  LINALG::Matrix<nsd_,nsd_> xjm_;
  //! inverse of transposed jacobian "ds/dx"
  LINALG::Matrix<nsd_,nsd_> xij_;
  //! 2nd derivatives of coord.-functions w.r.t r,s,t
  LINALG::Matrix<numderiv2_,nsd_> xder2_;

  //! (global) Laplacian of shape functions at integration point
  LINALG::Matrix<nen_,1> laplace_;
  //! rhs force in Gauss point
  std::vector<double> rhs_;
  //! rhs of temperature equation in case of reactive equation system
  std::vector<double> reatemprhs_;
  //! bodyforce in element nodesx
  std::vector<LINALG::Matrix<nen_,1> > bodyforce_;
  //! residual of scalar transport equation
  std::vector<double> scatrares_;
  //! convective term (part of residual) ?
  std::vector<double> conv_phi_;
  //! diffusive term (part of residual) ?
  std::vector<double> diff_phi_;
  //! reactive term (part of residual) ?
  std::vector<double> rea_phi_;

  //! the stabilisation parameters (one per transported scalar)
  std::vector<double> tau_;
  //! derivative of stabilization parameter with respect to el. potential (ELCH)
  std::vector<LINALG::Matrix<nen_,1> > tauderpot_;

  //! reconstructed flux values via L2 projection
  std::vector<LINALG::Matrix<nsd_,nen_> > efluxreconstr_;

  /*========================================================================*/
  //! @name Miscellaneous - to be done!!!
  /*========================================================================*/



}; // end class DRT::ELEMENTS::ReInitImpl

} // end namespace ELEMENTS
} // end namespace DRT

#endif // SCATRA_ELE_IMPL_H
