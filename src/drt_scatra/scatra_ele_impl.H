/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_impl.H

\brief Internal implementation of scalar transport elements

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_IMPL_H
#define SCATRA_ELE_IMPL_H

#include "scatra_ele_impl_interface.H"
#include "scatra_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"
#include <Epetra_IntSerialDenseVector.h>
#include "../drt_lib/drt_singletondestruction.H"

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
class ScaTraImpl: public ScaTraImplInterface
{
public:
  //! Constructor
  ScaTraImpl(const int numdofpernode, const int numscal);

  //! Singleton access method
  static ScaTraImpl<distype> * Instance(const int numdofpernode, const int numscal, bool create=true );

  //! called upon destruction
  virtual void Done();

  //! Evaluate
  virtual int Evaluate(
    DRT::Element*              ele,            //!< the element we are dealing with
    ParameterList&             params,         //!< the parameter list
    DRT::Discretization&       discretization, //!< the discretization the element is part of
    vector<int>&               lm,             //!< location matrix
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

  //! number of components necessary to store second derivatives
  // 1 component  for nsd=1:  (N,xx)
  // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
  // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
  static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;


private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! Calculate matrix and rhs. Here the whole thing is hidden.
  virtual void Sysmat(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff, //!< subgrid-diff.-scaling vector
    const double                          time,      //!< current simulation time
    const double                          dt,        //!< current time-step length
    const double                          timefac,   //!< time discretization factor
    const double                          alphaF,    //!< factor for generalized-alpha time integration
    const enum INPAR::SCATRA::AssgdType   whichassgd,//!< all-scale subgrid-diffusivity definition
    const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd,//!< fine-scale subgrid-diffusivity definition
    const bool                            assgd,     //!< all-scale subgrid-diff. flag
    const bool                            fssgd,     //!< fine-scale subgrid-diff. flag
    const double                          Cs,        //!< Smagorinsky constant
    const double                          tpn,       //!< turbulent Prandtl number
    const double                          Csgs_sgvel,//!< parameter of multifractal subgrid-scales (velocity)
    const double                          alpha,     //!< grid-filter to test-filter ratio
    const bool                            calc_N,    //!< flag to activate calculation of N
    const double                          N_vel,     //!< value for N if not calculated
    const enum INPAR::FLUID::RefVelocity  refvel,    //!< reference velocity
    const enum INPAR::FLUID::RefLength    reflength, //!< reference length
    const double                          c_nu,      //!< scaling for Re
    const bool                            nwl,       //!< flag to activate near-wall limit
    const double                          Csgs_sgphi,//!< parameter of multifractal subgrid-scales (phi)
    const double                          c_diff,    //!< scaling for Re*Pr
    const bool                            BD_gp,     //!< evaluation of model coefficient at gp
    const double                          frt,       //!< factor F/RT needed for ELCH calculations
    const enum INPAR::SCATRA::ScaTraType  scatratype //!< type of scalar transport problem
  );

  //! get the body force
  virtual void BodyForce(
    const DRT::Element*  ele, //!< the element we are dealing with
    const double         time //!< current simulation time
  );

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element*                   ele,        //!< the element we are dealing with
    const enum INPAR::SCATRA::ScaTraType  scatratype, //!< type of scalar transport problem
    const double                          dt          //!< current time-step length
  );

  //! calculate element matrix and rhs vector (One-Step-Theta, BDF2, Gen.-alpha)
  virtual void CalMatAndRHS(
    Epetra_SerialDenseMatrix&  emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector&  erhs,     //!< element rhs to calculate
    const double               fac,      //!< integration factor
    const bool                 fssgd,    //!< fine-scale subgrid-diff. flag
    const double               timefac,  //!< time discretization factor
    const double               dt,       //!< current time-step length
    const double               alphaF,   //!< factor for generalized-alpha time integration
    const int                  k         //!< index of current scalar
  );

  //! integral of shape functions over the element
  void IntegrateShapeFunctions(
    const DRT::Element*                ele,    //!< the current element
    Epetra_SerialDenseVector&          elevec1,//!< result vector (to be assembled)
    const Epetra_IntSerialDenseVector& dofids  //!< for which dof we need to integrate?
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

  //! calculate mass matrix and rhs for determining initial time derivative
  virtual void CalcInitialTimeDerivative(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    const enum INPAR::SCATRA::ScaTraType  scatratype, //!< type of scalar transport problem
    ParameterList&                        params,         //!< the parameter list
    DRT::Discretization&                  discretization, //!< the discretization the element is part of
    vector<int>&                          lm              //!< location matrix
  );

  //! calculate error compared to analytical solutions
  virtual void CalErrorComparedToAnalytSolution(
    const DRT::Element*          ele,                 //!< the element
    const enum INPAR::SCATRA::ScaTraType  scatratype, //!< type of scalar transport problem
    ParameterList&               params,              //!< parameter list
    Epetra_SerialDenseVector&    errors               //!< computed error contributions
  );

  //! calculate flux vector for actual solution
  virtual void CalculateFlux(
    LINALG::Matrix<3,nen_>&         flux,     //!< flux to calculate
    const DRT::Element*             ele,      //!< the element we are dealing with
    const double                    frt,      //!< factor F/RT needed for ELCH calculations
    const INPAR::SCATRA::FluxType   fluxtype, //!< type of flux to calculate
    const int                       k,        //!< index of current scalar
    const enum INPAR::SCATRA::ScaTraType  scatratype, //!< type of scalar transport problem
    const double                    dt        //!< current time-step length
  );

  //! calculate scalars and domain integral
  virtual void CalculateScalars(
    const DRT::Element*       ele,      //!< the element we are dealing with
    const vector<double>&     ephinp,   //!< current scalar field
    Epetra_SerialDenseVector& scalars,  //!< vector for computed integrals
    const bool                inverting //!< flag for invertuing
  );

   //! calculate domain and bodyforce integral
  virtual void CalculateDomainAndBodyforce(
    Epetra_SerialDenseVector&  scalars, //!< vector for the computed integrals
    const DRT::Element*        ele,     //!< pointer to the element
    const double               time     //!< current simulation time
  );

  //! Finite difference check for debugging purposes
  void FDcheck(
    DRT::Element*                         ele,          //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,         //!< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff,    //!< subgrid-diff.-scaling vector
    const double                          time,         //!< current simulation time
    const double                          dt,           //!< current time-step length
    const double                          timefac,      //!< time discretization factor
    const double                          alphaF,       //!< factor for generalized-alpha time integration
    const enum INPAR::SCATRA::AssgdType   whichassgd,   //!< all-scale subgrid-diffusivity definition
    const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd,   //!< fine-scale subgrid-diffusivity definition
    const bool                            assgd,        //!< all-scale subgrid-diff. flag
    const bool                            fssgd,        //!< fine-scale subgrid-diff. flag
    const enum INPAR::FLUID::TurbModelAction turbmodel, //!< turbulence model action
    const double                          Cs,           //!< Smagorinsky constant
    const double                          tpn,          //!< turbulent Prandtl number
    const double                          frt,          //!< factor F/RT needed for ELCH calculations
    const enum INPAR::SCATRA::ScaTraType  scatratype    //!< type of scalar transport problem
  );


  /*========================================================================*/
  //! @name stabilization
  /*========================================================================*/

  //! calculate normalized subgrid-diffusivity matrix
  virtual void CalcSubgrDiffMatrix(
    const DRT::Element*        ele,    //!< the element we are dealing with
    Epetra_SerialDenseMatrix&  emat,   //!< element matrix to calculate
    const double               timefac //!< time discretization factor
  );

  //! update material parameters based on inclusion of subgrid-scale part of
  //! scalar (only for special materials)
  virtual void UpdateMaterialParams(
    const DRT::Element*  ele,   //!< the element we are dealing with
    const double         sgphi, //!< approximation of subgrid scalar
    const int            k      //!< index of current scalar
  );

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

  //! calculate subgrid-scale velocity
  virtual void CalcSubgrVelocity(
    DRT::Element*                         ele,       //!< the element we are dealing with
    const double                          time,      //!< current simulation time
    const double                          dt,        //!< current time-step length
    const double                          timefac,   //!< time discretization factor
    const int                             k,         //!< index of current scalar
    const enum INPAR::SCATRA::ScaTraType  scatratype //!< type of scalar transport problem
  );


  //! calculate residual of scalar transport equation and subgrid-scale part of scalar
  virtual void CalcResidualAndSubgrScalar(
    const double  dt,       //!< current time-step length
    const double  timefac,  //!< time discretization factor
    const int     k         //!< index of current scalar
  );

  //! calculate rate of strain of (fine-scale) velocity
  inline double GetStrainRate(
    const LINALG::Matrix<nsd_,nen_>& evel
  )
  {
    // evel is tranferred here since the evaluation of the strain rate can be performed
    // for various velocities such as velint_, fsvel_, ...

    double rateofstrain=0;

    // get velocity derivatives at integration point
    //
    //              +-----  dN (x)
    //   dvel (x)    \        k
    //   -------- =   +     ------ * vel
    //      dx       /        dx        k
    //        j     +-----      j
    //              node k
    //
    // j : direction of derivative x/y/z
    //
    LINALG::Matrix<nsd_,nsd_> velderxy;
    velderxy.MultiplyNT(evel,derxy_);

    // compute (resolved) rate of strain
    //
    //          +-                                 -+ 1
    //          |          /   \           /   \    | -
    //          | 2 * eps | vel |   * eps | vel |   | 2
    //          |          \   / ij        \   / ij |
    //          +-                                 -+
    //
    LINALG::Matrix<nsd_,nsd_> two_epsilon;
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
      }
    }

    for(int rr=0;rr<nsd_;rr++)
    {
      for(int mm=0;mm<nsd_;mm++)
      {
        rateofstrain += two_epsilon(rr,mm)*two_epsilon(mm,rr);
      }
    }

    return(sqrt(rateofstrain/2.0));
  };


  /*========================================================================*/
  //! @name ELCH ?
  /*========================================================================*/

  //! calculate element matrix and rhs vector (Electrochemistry problems)
  virtual void CalMatElch(
    Epetra_SerialDenseMatrix&    emat,   //!< element matrix to calculate
    Epetra_SerialDenseVector&    erhs,   //!< element rhs to calculate
    const double                 frt,    //!< factor F/RT needed for ELCH calculations
    const double                 timefac,//!< time discretization factor
    const double                 alphaF, //!< factor for gen. alpha time integration
    const double                 fac,    //!< integration factor
    const enum INPAR::SCATRA::ScaTraType  scatratype //!< type of scalar transport problem
  );

  //! calculate conductivity of electrolyte solution
  void CalculateConductivity(
    const DRT::Element*                  ele,        //!< the current element
    const double                         frt,        //!< factor F/RT needed for ELCH calculations
    const enum INPAR::SCATRA::ScaTraType scatratype, //!< type of scalar transport problem
    Epetra_SerialDenseVector&            sigma       //!< values for different conductivity expressions
  );

  //! calculate el. potential field based on ion concentrations
  void CalculateElectricPotentialField(
    const DRT::Element*         ele,  //!<the current element
    const double                frt,  //!< factor F/RT needed for ELCH calculations
    const enum INPAR::SCATRA::ScaTraType  scatratype, //!< type of scalar transport problem
    Epetra_SerialDenseMatrix&   emat, //!< matrix (to be assembled)
    Epetra_SerialDenseVector&   erhs  //!< rhs-vector (to be assembled)
  );

  /*========================================================================*/
  //! @name Poro-Scatra ?
  /*========================================================================*/

  //! modify residual of scalar transport equation for poroelastic problem
  virtual void CalcResidual_PoroScatraMod(
      const double        dt,
      const double        timefac,
      const int           k,
      const double        porosity,
      const double            dporodt,
      LINALG::Matrix<3,1>&   gradporosity
  );

  //! modify element matrix and rhs vector for Poroelastic+Scatra problem (One-Step-Theta, BDF2, Gen.-alpha)
  virtual void CalMatAndRHS_PoroScatraMod(
      Epetra_SerialDenseMatrix&             emat,     //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs,     //!< element rhs to calculate
      const double                          fac,      //!< integration factor
      const double                          timefac,  //!< time discretization factor
      const int                             k,        //!< index of current scalar
      const int                             iquad,    //!< id of current Gauss point
      const int                             eleid     //!< the element id
  );

  /*========================================================================*/
  //! @name Turbulence -> implemented in scatra_ele_imlp_turbulence_service.cpp
  /*========================================================================*/

  //! calculate filtered quantities for dynamic Smagorinsky model
  void scatra_apply_box_filter(
          const double               thermpress,
          double&                    dens_hat,
          double&                    temp_hat,
          double&                    dens_temp_hat,
          RCP<vector<double> >       vel_hat,
          RCP<vector<double> >       densvel_hat,
          RCP<vector<double> >       densveltemp_hat,
          RCP<vector<double> >       densstraintemp_hat,
          double&                    volume,
          const DRT::Element*        ele
  );

  //! calculate turbulent Prandtl number for dynamic Smagorinsky model
  void scatra_calc_smag_const_LkMk_and_MkMk(
          RCP<Epetra_MultiVector>&  col_filtered_vel,
          RCP<Epetra_MultiVector>&  col_filtered_dens_vel,
          RCP<Epetra_MultiVector>&  col_filtered_dens_vel_temp,
          RCP<Epetra_MultiVector>&  col_filtered_dens_rateofstrain_temp,
          RCP<Epetra_Vector>&       col_filtered_temp,
          RCP<Epetra_Vector>&       col_filtered_dens,
          RCP<Epetra_Vector>&       col_filtered_dens_temp,
          double&                   LkMk,
          double&                   MkMk,
          double&                   xcenter,
          double&                   ycenter,
          double&                   zcenter,
          const DRT::Element*       ele
  );

  //! calculate mean turbulent Prandtl number
  void GetMeanPrtOfHomogenousDirection(
          ParameterList&             turbmodelparams,
          double&                    inv_Prt,
          int&                       nlayer
          );

  //! calculate all-scale subgrid diffusivity
  virtual void CalcSubgrDiff(  // ? subgrid = subgr or subgrid?
    const double                          dt,        //!< current time-step length
    const double                          timefac,   //!< time discretization factor
    const enum INPAR::SCATRA::AssgdType   whichassgd,//!< all-scale subgrid-diffusivity definition
    const bool                            assgd,     //!< all-scale subgrid-diff. flag
    const double                          Cs,        //!< Smagorinsky constant
    const double                          tpn,       //!< turbulent Prandtl number
    const double                          vol,       //!< element volume
    const int                             k          //!< index of current scalar
  );

  //! calculate fine-scale subgrid diffusivity
  virtual void CalcFineScaleSubgrDiff(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseVector&             subgrdiff, //!< subgrid-diff.-scaling vector ?
    const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd,//!< all-scale subgrid-diffusivity definition
    const double                          Cs,        //!< Smagorinsky constant
    const double                          tpn,       //!< turbulent Prandtl number
    const double                          vol,       //!< element volume
    const int                             k          //!< index of current scalar
  );

  //! calculate coefficients of mutifractal subgrid-scale modeling
  virtual void CalcBAndDForMultifracSubgridScales(
    LINALG::Matrix<nsd_,1>&               B_mfs,     //!< coefficient for fine-scale velocity (will be filled)
    double &                              D_mfs,     //!< coefficient for fine-scale scalar (will be filled)
    const double                          Csgs_sgvel,//!< parameter of multifractal subgrid-scales (velocity)
    const double                          alpha,     //!< grid-filter to test-filter ratio
    const bool                            calc_N,    //!< flag to activate calculation of N
    const double                          N_vel,     //!< value for N if not calculated
    const enum INPAR::FLUID::RefVelocity  refvel,    //!< reference velocity
    const enum INPAR::FLUID::RefLength    reflength, //!< reference length
    const double                          c_nu,      //!< scaling for Re
    const bool                            nwl,       //!< flag to activate near-wall limit
    const double                          Csgs_sgphi,//!< parameter of multifractal subgrid-scales (phi)
    const double                          c_diff,    //!< scaling for Re*Pr
    const double                          vol,       //!< volume of element
    const int                             k          //!< index of current scalar
  );

  //! calculate reference length for element Reynolds number
  double CalcRefLength(const enum INPAR::FLUID::RefLength reflength, //!< reference length
                       const double vol                              //!< volume of element
  );

  //! store turbulent Prandtl number and effective diffusivity for output
  virtual void StoreModelParametersForOutput(
    const bool                            isowned,        //!< owner of element
    ParameterList&                        turbulencelist, //!< turbulence parameter list
    const int                             nlayer,         //!< element layer
    const double                          tpn             //! turbulent Prandlt number
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

  //! flag for electrochemistry problems (ELCH)
  const bool is_elch_;
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

  //! flag for inclusion of subgrid-scale velocity
  bool sgvel_;

  //! flag for material evaluation at Gaussian integration points
  bool mat_gp_;
  //! flag for stabilization parameter evaluation at Gaussian integration points
  bool tau_gp_;
  //! flag for improvement of consistency of stabilized method
  bool betterconsistency_;

  //! flag for migration velocity in formula for stabilization parameter (ELCH)
  bool migrationintau_;
  //! flag for migration velocity in the stabilization operator (ELCH)
  bool migrationstab_;
  //! flag for migration velocity in residual of strong form (ELCH)
  bool migrationinresidual_;

  //! update material with subgrid-scale temperature (loma specific)
  bool update_mat_;

  //! stabilization parameter definition
  enum INPAR::SCATRA::TauType whichtau_;
  //! turbulence model
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! scalar at integration point
  vector<double> phi_;
  //! subgrid-scale part of scalar
  vector<double> sgphi_;
  //! multifractal subgrid-scale part of scalar
  vector<double> mfssgphi_;

  //! (global) gradient of transported scalar at integration point
  LINALG::Matrix<nsd_,1> gradphi_;
  //! (global) gradient of fine-scale part of scalar at integration point
  LINALG::Matrix<nsd_,1> fsgradphi_;
  //! (global) gradient of multifractal subgrid-scale scalar at integration point
  LINALG::Matrix<nsd_,1> mfsggradphi_;

  //! scalar at t_(n)
  vector<LINALG::Matrix<nen_,1> > ephin_;
  //! scalar at t_(n+1) or t_(n+alpha_F)
  vector<LINALG::Matrix<nen_,1> > ephinp_;
  //! scalar at t_(n+alpha_M)
  vector<LINALG::Matrix<nen_,1> > ephiam_;

  //! vector containing all values from previous timelevel n in gausspoint
  vector<double> hist_;
  //! history vector of transported scalars
  vector<LINALG::Matrix<nen_,1> > ehist_;

  //! scalar at t_(0) used for reinitialization as reference
  vector<LINALG::Matrix<nen_,1> > ephi0_Reinit_Reference_;
  //! scalar at t_(0) used for reinitialization as reference
  vector<LINALG::Matrix<nen_,1> > ephi0_penalty_;

  //! fine-scale solution ?
  vector<LINALG::Matrix<nen_,1> > fsphinp_;

  //! gradient of phi at global time for reinitialization(initialization for pseudo time loop)
  LINALG::Matrix<nsd_,1> grad_phi_0;

  /*========================================================================*/
  //! @name potential (ELCH)
  /*========================================================================*/

  //! values of all transported scalars in gausspoint (ELCH)
  vector<double> conint_;

  //! electric potential
  LINALG::Matrix<nen_,1> epotnp_;
  //! magnetic field at t_(n+1)
  LINALG::Matrix<nsd_,nen_> emagnetnp_;

  //! (global) gradient of el. potential at integration point
  LINALG::Matrix<nsd_,1> gradpot_;

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
  //! multifractal subgrid-scale velocity divergence
  double mfsvdiv_;

  //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<nen_,1> eprenp_;
  //! thermodynamic pressure at t_(n+1) or t_(n+alpha_F) (LOMA specific)
  double thermpressnp_;
  //! thermodynamic pressure at t_(n+alpha_M) (LOMA specific)
  double thermpressam_;
  //! time derivative of thermodynamic pressure (LOMA specific)
  double thermpressdt_;

  /*========================================================================*/
  //! @name coefficients and related
  /*========================================================================*/

  //! density at t_(n)
  vector<double> densn_;
  //! density at t_(n+1) or t_(n+alpha_F)
  vector<double> densnp_;
  //! density at t_(n+alpha_M)
  vector<double> densam_;
  //! factor for density derivative (at t_(n+1) or t_(n+alpha_F))
  vector<double> densgradfac_;

  //! diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  vector<double> diffus_;
  //! all-scale or fine-scale subgrid diffusivities (one per transported scalar)
  vector<double> sgdiff_;
  //! general reaction term R(phi) (without density)
  vector<double> reacterm_;
  //! reaction coefficient
  vector<double> reacoeff_;
  //! reaction coefficent derivative (needed for matrix contribution)
  vector<double> reacoeffderiv_;

  //! charge numbers of all ion species (ELCH specific)
  vector<double> valence_;
  //! diffusivity times valence for each ion species (ELCH specific)
  vector<double> diffusvalence_;

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
  vector<double> rhs_;
  //! rhs of temperature equation in case of reactive equation system
  vector<double> reatemprhs_;
  //! bodyforce in element nodesx
  vector<LINALG::Matrix<nen_,1> > bodyforce_;
  //! residual of scalar transport equation
  vector<double> scatrares_;
  //! convective term (part of residual) ?
  vector<double> conv_phi_;
  //! diffusive term (part of residual) ?
  vector<double> diff_phi_;
  //! reactive term (part of residual) ?
  vector<double> rea_phi_;

  //! the stabilisation parameters (one per transported scalar)
  vector<double> tau_;
  //! derivative of stabilization parameter with respect to el. potential (ELCH)
  vector<LINALG::Matrix<nen_,1> > tauderpot_;

  //! reconstructed flux values via L2 projection
  vector<LINALG::Matrix<nsd_,nen_> > efluxreconstr_;

  //! weights for nurbs elements
  LINALG::Matrix<nen_,1> weights_;
  //! knot vector for nurbs elements
  std::vector<Epetra_SerialDenseVector> myknots_;

  /*========================================================================*/
  //! @name Miscellaneous - to be done!!!
  /*========================================================================*/



}; // end class DRT::ELEMENTS::ScaTraImpl

} // end namespace ELEMENTS
} // end namespace DRT

#endif // SCATRA_ELE_IMPL_H
