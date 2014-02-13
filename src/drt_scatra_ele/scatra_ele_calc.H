/*!----------------------------------------------------------------------
\file scatra_ele_calc.H

\brief main file containing routines for calculation of scatra element

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_H
#define SCATRA_ELE_CALC_H

#include "scatra_ele_interface.H"

#include "scatra_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include <Epetra_IntSerialDenseVector.h>

// forward declarations

namespace DRT
{
namespace ELEMENTS
{
  class ScaTraEleParameter;
  class ScaTraEleParameterTimInt;

  class ScaTraEleDiffManager;
  class ScaTraEleReaManager;

  /// Scatra element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the scatra element. Additionally, the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid) from the
    mathematical contents (this class). There are different
    implementations of the scatra element, this is just one such
    implementation.

    The scatra element will allocate exactly one object of this class for all
    scatra elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>Usability</h3>

    The calculations are done by the Evaluate() method. There are two
    version. The virtual method that is inherited from ScatraEleInterface
    (and called from Scatra) and the non-virtual one that does the actual
    work. The non-virtual Evaluate() method must be callable without an actual
    Scatra object.
  */

  template<DRT::Element::DiscretizationType distype>
  class ScaTraEleCalc : public ScaTraEleInterface
  {
  protected:

    /// (private) protected constructor, since we are a Singleton.
    /// this constructor is called from a derived class
    /// -> therefore, it has to be protected instead of private
    ScaTraEleCalc(const int numdofpernode,const int numscal);

  public:

    virtual ~ScaTraEleCalc() {}

    /// called upon destruction
    virtual void Done() = 0;

    /// In this class we do not define a static ScaTraEle...* Instance
    /// since only derived child classes are free to be allocated!!

    /// Evaluate the element
     /*!
       Generic virtual interface function. Called via base pointer.
      */
     virtual int Evaluate(DRT::ELEMENTS::Transport*     ele,
                          Teuchos::ParameterList&       params,
                          DRT::Discretization &         discretization,
                          const std::vector<int> &      lm,
                          Epetra_SerialDenseMatrix&     elemat1_epetra,
                          Epetra_SerialDenseMatrix&     elemat2_epetra,
                          Epetra_SerialDenseVector&     elevec1_epetra,
                          Epetra_SerialDenseVector&     elevec2_epetra,
                          Epetra_SerialDenseVector&     elevec3_epetra);

     virtual int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                          Teuchos::ParameterList&       params,
                          DRT::Discretization &         discretization,
                          const std::vector<int> &      lm,
                          Epetra_SerialDenseMatrix&     elemat1_epetra,
                          Epetra_SerialDenseMatrix&     elemat2_epetra,
                          Epetra_SerialDenseVector&     elevec1_epetra,
                          Epetra_SerialDenseVector&     elevec2_epetra,
                          Epetra_SerialDenseVector&     elevec3_epetra);

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

  protected:

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! extract element based or nodal values
    //  return extracted values of phinp
    std::vector<double> ExtractElementAndNodeValues(
      DRT::ELEMENTS::Transport*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      const std::vector<int>&    lm
    );

    //! extract turbulence approach
    void ExtractTurbulenceApproach(
      DRT::ELEMENTS::Transport*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      const std::vector<int>&    lm,
      int&                       nlayer
    );

    //! calculate matrix and rhs. Here the whole thing is hidden.
    virtual void Sysmat(
      DRT::Element*                         ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
      Epetra_SerialDenseVector&             subgrdiff  //!< subgrid-diff.-scaling vector
    );

    //! get the body force
    virtual void BodyForce(
      const DRT::Element*  ele //!< the element we are dealing with
    );

    //! further node-based source terms not given via Neumann volume condition
    void OtherNodeBasedSourceTerms(
      const std::vector<int>&    lm,                //!< location vector
      DRT::Discretization&       discretization,    //!< discretization
      Teuchos::ParameterList&    params             //!< parameterlist
    );

    //! evaluate shape functions and their derivatives at element center
    const double EvalShapeFuncAndDerivsAtEleCenter();

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsAtIntPoint(
      const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints, //!< integration points
      const int                                    iquad      //!< id of current Gauss point
    );

    /*========================================================================*/
    //! @name routines for additional element evaluations (called from EvaluateService)
    /*========================================================================*/
    
    //! calculate mass matrix and rhs for determining initial time derivative
    virtual void CalcInitialTimeDerivative(
      DRT::ELEMENTS::Transport*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
      Teuchos::ParameterList&               params,         //!< the parameter list
      DRT::Discretization&                  discretization, //!< the discretization the element is part of
      const std::vector<int>&               lm              //!< location matrix
    );

    //!  integrate shape functions over domain
    void IntegrateShapeFunctions(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseVector&       elevec1,   //!< rhs vector
      const Epetra_IntSerialDenseVector& dofids  //!< index of current scalar
    );

    //!  calculate weighted mass flux (no reactive flux so far)
    virtual void CalculateFlux(
      LINALG::Matrix<3,nen_>&         flux,      //!< flux to be computed
      const DRT::Element*             ele,       //!< the element we are dealing with
      const INPAR::SCATRA::FluxType   fluxtype,  //!< type fo flux
      const int                       k          //!< index of current scalar
    );

    //! calculate scalar(s) and domain integral
    void CalculateScalars(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseVector&       scalars,   //!< scalar to be computed
      const bool                      inverting  //!< bool indicating inversion
    );

    /*========================================================================*/
    //! @name material and related and related functions
    /*========================================================================*/

    //! get the material parameters
    virtual void GetMaterialParams(
      const DRT::Element* ele,       //!< the element we are dealing with
      double&             densn,     //!< density at t_(n)
      double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
      double&             densam,    //!< density at t_(n+alpha_M)
      Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
      Teuchos::RCP<ScaTraEleReaManager>  reamanager,   //!< reaction manager
      double&             visc,       //!< fluid viscosity
      const int           iquad=-1    //!< id of current gauss point (default = -1)
    );

    //! evaluate material
    virtual void Materials(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int                               k,        //!< id of current scalar
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
      Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
      double&                                 visc,      //!< fluid viscosity
      const int                               iquad=-1   //!< id of current gauss point (default = -1)
    );

    //! material ScaTra
    virtual void MatScaTra(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int                               k,        //!< id of current scalar
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
      Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
      double&                                 visc,      //!< fluid viscosity
      const int                               iquad=-1   //!< id of current gauss point (default = -1)
    );

    //! material MyoCard
    void MatMyoCard(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      double&                                 diffus,   //!< diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
      double&                                 visc      //!< fluid viscosity
    );

    /*========================================================================*/
    //! @name stabilization and related functions
    /*========================================================================*/

    //! calculate stabilization parameter
    virtual void CalcTau(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol,       //!< element volume
      const int                    k          //!< index of current scalar
    );

    //! calculate characteristic element length
    virtual double CalcCharEleLength(
      const double                  vol,        //!< element volume
      const double                  vel_norm,   //!< norm of velocity
      const LINALG::Matrix<nsd_,1>& convelint   //!< convective velocity at integration point
    );

    //! calculation of tau according to Taylor, Hughes and Zarins
    virtual void CalcTauTaylorHughesZarins(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint  //!< convective velocity at integration point
    );

    //! calculation of tau according to Franca and Valentin
    virtual void CalcTauFrancaValentin(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of tau according to Franca, Shakib and Codina
    virtual void CalcTauFrancaShakibCodina(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of tau according to Codina
    virtual void CalcTauCodina(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of tau according to Franca, Madureira and Valentin
    virtual void CalcTauFrancaMadureiraValentin(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const double                 vol        //!< element volume
    );

    //! exact calculation of tau for 1D
    virtual void CalcTau1DExact(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of residual and residual-based subgrid-scale scalar for stabilization
    void CalcResidualAndSubgrScalar(
      const int      k,          //!< index of current scalar
      double&        scatrares,  //!< residual of convection-diffusion-reaction eq
      double&        sgphi,      //!< residual-based subgrid-scale scalar
      const double   densam,     //!< density at t_(n+am)
      const double   densnp,     //!< density at t_(n+1)
      const double   phinp,      //!< scalar at t_(n+1)
      const double   hist,       //!< history of time integartion
      const double   conv_phi,   //!< convective contribution
      const double   diff_phi,   //!< diffusive contribution
      const double   rea_phi,    //!< reactive contribution
      const double   rhsint,     //!< rhs at integration point
      const double   tau         //!< the stabilisation parameter
    );

    //! calculate artificial diffusivity
    void CalcArtificialDiff(
      const double                  vol,        //!< element volume
      const int                     k,          //!< id of current scalar
      Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager
      const double                  densnp,     //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint,  //!< convective velocity at integration point
      const LINALG::Matrix<nsd_,1>& gradphi,    //!< scalar gradient
      const double                  conv_phi,   //!< convective contribution
      const double                  scatrares,  //!< residual of convection-diffusion-reaction eq
      const double                  tau         //!< the stabilisation parameter
    );

    //! calculate subgrid-scale velocity
    void CalcSubgrVelocity(
      const DRT::Element*           ele,        //!< the element we are dealing with
      LINALG::Matrix<nsd_,1>&       sgvelint,   //!< subgrid velocity at integration point
      const double                  densam,     //!< density at t_(n+am)
      const double                  densnp,     //!< density at t_(n+1)
      const double                  visc,       //!< fluid viscosity
      const LINALG::Matrix<nsd_,1>& convelint,  //!< convective velocity at integration point
      const double                  tau         //!< the stabilisation parameter
    );

    /*========================================================================*/
    //! @name turbulence and related functions
    /*========================================================================*/

    //! calculate mean turbulent Prandtl number
    void GetMeanPrtOfHomogenousDirection(
      Teuchos::ParameterList&    turbmodelparams,  //!< turbulence parameter list
      int&                       nlayer            //!< layer of homogeneous plane
    );

    //! output of model parameters
    void StoreModelParametersForOutput(
      const DRT::Element*                   ele,             //!< the element we are dealing with
      const bool                            isowned,         //!< owner
      Teuchos::ParameterList&               turbulencelist,  //!< turbulence parameter list
      const int                             nlayer           //!< layer of homogeneous plane
    );

    //! calculate all-scale art. subgrid diffusivity
    void CalcSubgrDiff(
      Teuchos::RCP<ScaTraEleDiffManager>    diffmanager,  //!< diffusion manager
      double&                               visc,      //!< fluid viscosity
      const double                          vol,       //!< element volume
      const int                             k,         //!< index of current scalar
      const double                          densnp     //!< density at t_(n+1)
   );

    //!  calculate fine-scale art. subgrid diffusivity
    void CalcFineScaleSubgrDiff(
      double&                               sgdiff,    //!< subgrid-scale diffusion
      Epetra_SerialDenseVector&             subgrdiff, //!< subgrid-scale diffusion vector
      DRT::Element*                         ele,       //!< the element we are dealing with
      const double                          vol,       //!< element volume
      const int                             k,         //!< index of current scalar
      const double                          densnp,    //!< density at t_(n+1)
      const double                          diffus,    //!< diffusion
      const LINALG::Matrix<nsd_,1>          convelint  //!< convective velocity
    );

    //! calculation of coefficients B and D for multifractal subgrid-scales
    void CalcBAndDForMultifracSubgridScales(
      LINALG::Matrix<nsd_,1>&                     B_mfs,      //!< coefficient for fine-scale velocity (will be filled)
      double &                                    D_mfs,      //!< coefficient for fine-scale scalar (will be filled)
      const double                                vol,        //!< volume of element
      const int                                   k,          //!< index of current scalar
      const double                                densnp,     //!< density at t_(n+1)
      const double                                diffus,     //!< diffusivity
      const double                                visc,       //!< fluid vicosity
      const LINALG::Matrix<nsd_,1>                convelint,  //!< convective velocity
      const LINALG::Matrix<nsd_,1>                fsvelint    //!< fine-scale velocity
    );

    // calculate reference length for multifractal subgrid-scales
    double CalcRefLength(
      const double                       vol,        //!< volume of element
      const LINALG::Matrix<nsd_,1>       convelint  //!< convective velocity
    );

    //! calculate filtered quantities for dynamic Smagorinsky model
    void scatra_apply_box_filter(
      const double               thermpress,
      double&                    dens_hat,
      double&                    temp_hat,
      double&                    dens_temp_hat,
      double&                    phi2_hat,
      double&                    phiexpression_hat,
      Teuchos::RCP<std::vector<double> >       vel_hat,
      Teuchos::RCP<std::vector<double> >       densvel_hat,
      Teuchos::RCP<std::vector<double> >       densveltemp_hat,
      Teuchos::RCP<std::vector<double> >       densstraintemp_hat,
      Teuchos::RCP<std::vector<double> >       phi_hat,
      Teuchos::RCP<std::vector<std::vector<double> > > alphaijsc_hat,
      double&                    volume,
      const DRT::Element*        ele
    );

    //! calculate turbulent Prandtl number for dynamic Smagorinsky model
    void scatra_calc_smag_const_LkMk_and_MkMk(
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_vel,
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_vel,
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_vel_temp,
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_rateofstrain_temp,
      Teuchos::RCP<Epetra_Vector>&       col_filtered_temp,
      Teuchos::RCP<Epetra_Vector>&       col_filtered_dens,
      Teuchos::RCP<Epetra_Vector>&       col_filtered_dens_temp,
      double&                   LkMk,
      double&                   MkMk,
      double&                   xcenter,
      double&                   ycenter,
      double&                   zcenter,
      const DRT::Element*       ele
    );

    void scatra_calc_vreman_dt(
      Teuchos::RCP<Epetra_MultiVector>& col_filtered_phi,
      Teuchos::RCP<Epetra_Vector>& col_filtered_phi2,
      Teuchos::RCP<Epetra_Vector>&   col_filtered_phiexpression,
      Teuchos::RCP<Epetra_MultiVector>& col_filtered_alphaijsc,
      double& dt_numerator,
      double& dt_denominator,
      const DRT::Element* ele
    );

    //! calculate normalized subgrid-diffusivity matrix
    virtual void CalcSubgrDiffMatrix(
      const DRT::Element*        ele,    //!< the element we are dealing with
      Epetra_SerialDenseMatrix&  emat   //!< element matrix to calculate
    );

    //! calculate dissipation introduced by stabilization and turbulence models
    void CalcDissipation(
      Teuchos::ParameterList&               params,          //!< parameter list
      DRT::ELEMENTS::Transport*             ele,             //!< pointer to element
      DRT::Discretization&                  discretization,  //!< scatra discretization
      const std::vector<int>&               lm               //!< location matrix
    );

    /*========================================================================*/
    //! @name methods for evaluation of individual terms
    /*========================================================================*/

    //! calculate the Laplacian for all shape functions(strong form)
    void GetLaplacianStrongForm(
      LINALG::Matrix<nen_,1>&                diff  //!< laplace term to be computed
    );

    //! calculate divergence of vector field (e.g., velocity)
    void GetDivergence(
      double&                          vdiv,
      const LINALG::Matrix<nsd_,nen_>& evel
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

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakForm(
      double&                          val,   //!< ?
      const int                        vi,    //!< ?
      const int                        ui     //!< ?
    )
    {
      val = 0.0;
      for (int j = 0; j<nsd_; j++)
      {
        val += derxy_(j, vi)*derxy_(j, ui);
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakForm(
      double&                          val,   //!< ?
      const LINALG::Matrix<nsd_,nsd_>& diffus3, //!< ?
      const int                        vi,    //!< ?
      const int                        ui     //!< ?
    )
    {
      val = 0.0;
      for (int j = 0; j<nsd_; j++)
      {
        for (int i = 0; i<nsd_; i++)
        {
          val += derxy_(j, vi)*diffus3(j,i)*derxy_(i, ui);
        }
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakFormRHS(
      double& val,                              //!< ?
      const LINALG::Matrix<nsd_,1>&    gradphi, //!< ?
      const int vi                              //!< ?
    )
    {
      val = 0.0;
      for (int j = 0; j<nsd_; j++)
      {
        val += derxy_(j,vi)*gradphi(j);
      }
      return;
    };

    //! compute rhs containing bodyforce
    virtual void GetRhsInt(
      double&      rhsint,  //!< rhs containing bodyforce at integration point
      const double densnp,  //!< density at t_(n+1)
      const int    k        //!< index of current scalar
    );

    //! calculation of convective element matrix in convective form
    void CalcMatConv(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timefacfac,  //!< domain-integration factor times time-integration factor
      const double                  densnp,      //!< density  at time_(n+1)
      const LINALG::Matrix<nen_,1>& conv,        //!< convective operator
      const LINALG::Matrix<nen_,1>& sgconv       //!< subgrid-scale convective operator
    );

    //! calculation of convective element matrix: add conservative contributions
    virtual void CalcMatConvAddCons(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timefacfac,  //!< domain-integration factor times time-integration factor
      const LINALG::Matrix<nsd_,1>& convelint,   //!< convective velocity at Gauss point
      const LINALG::Matrix<nsd_,1>& gradphi,     //!< scalar gradient at Gauss point
      const double                  vdiv,        //!< velocity divergence
      const double                  densnp,      //!< density at time_(n+1)
      const double                  visc        //!< viscosity
    );

    //! calculation of diffusive element matrix
    virtual void CalcMatDiff(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timefacfac,  //!< domain-integration factor times time-integration factor
      Teuchos::RCP<ScaTraEleDiffManager>  diffmanager  //!< diffusion manager
    );

    //!calculation of stabilization element matrix
    void CalcMatTransConvDiffStab(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timetaufac,  //!< domain-integration factor times time-integration factor times tau
      const double                  densnp,      //!< density  at time_(n+1)
      const LINALG::Matrix<nen_,1>& conv,        //!< convective operator
      const LINALG::Matrix<nen_,1>& sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>& diff         //!< laplace term
    );

    //! calculation of mass element matrix
    void CalcMatMass(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  fac,         //!< domain-integration factor
      const double                  densam       //!< density at time_(n+am)
    );

    //! calculation of stabilization mass element matrix
    void CalcMatMassStab(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  taufac,      //!< domain-integration factor times tau
      const double                  densam,      //!< density at time_(n+am)
      const double                  densnp,      //!< density at time_(n+1)
      const LINALG::Matrix<nen_,1>& conv,        //!< convective operator
      const LINALG::Matrix<nen_,1>& sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>& diff         //!< laplace term
    );

    //! calculation of reactive element matrix
    virtual void CalcMatReact(
      Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
      const int                          k,           //!< index of current scalar
      const double                       timefacfac,  //!< domain-integration factor times time-integration factor
      const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
      const double                       taufac,      //!< domain-integration factor times tau
      const double                       densnp,      //!< density at time_(n+1)
      const double                       phinp,       //!< scalar at time_(n+1)
      Teuchos::RCP<ScaTraEleReaManager>  reamanager,  //!< reaction manager
      const LINALG::Matrix<nen_,1>&      conv,        //!< convective operator
      const LINALG::Matrix<nen_,1>&      sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>&      diff         //!< laplace term
    );

    //! calculation of linearized mass rhs vector
    void CalcRHSLinMass(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  fac,         //!< domain-integration factor
      const double                  densam,      //!< density at time_(n+am)
      const double                  densnp,      //!< density at time_(n+1)
      const double                  phinp,       //!< scalar at time_(n+1)
      const double                  hist         //!< history of time integartion
    );

    //! adaption of rhs with respect to time integration
    void ComputeRhsInt(
      double&                       rhsint,      //!< rhs at Gauss point
      const double                  densam,      //!< density at time_(n+am)
      const double                  densnp,      //!< density at time_(n+1)
      const double                  hist         //!< history of time integartion
    );

    //! adaption of residual with respect to time integration
    void RecomputeScatraResForRhs(
      double&                       scatrares,  //!< residual of convection-diffusion-reaction eq
      const int                     k,          //!< index of current scalar
      const LINALG::Matrix<nsd_,1>&  convelint, //!< convective velocity at Gauss point
      LINALG::Matrix<nsd_,1>&        gradphi,   //!< scalar gradient at Gauss point
      const LINALG::Matrix<nen_,1>&  diff,      //!< laplace term
      const double                  densn,      //!< density at time_(n)
      const double                  densnp,     //!< density at time_(n+1)
      double&                       conv_phi,   //!< convective contribution
      double&                       rea_phi,    //!< reactive contribution
      const double                  phin,       //!< scalar at time_(n)
      Teuchos::RCP<ScaTraEleReaManager> reamanager,  //!< reaction manager
      const double                  rhsint      //!< rhs at Gauss point
    );

    //! adaption of convective term for rhs
    virtual void RecomputeConvPhiForRhs(
      double&                       conv_phi,   //!< convective contribution
      const int                     k,          //!< index of current scalar
      const LINALG::Matrix<nsd_,1>&  sgvelint,  //!< subgrid-scale velocity at Gauss point
      const LINALG::Matrix<nsd_,1>&  gradphi,   //!< scalar gradient at Gauss point
      const double                  densnp,     //!< density at time_(n+1)
      const double                  densn,      //!< density at time_(n)
      const double                  phinp,      //!< scalar at time_(n+1)
      const double                  phin,       //!< scalar at time_(n)
      const double                  vdiv       //!< velocity divergence
    );

    //! standard Galerkin transient, old part of rhs and source term
    void CalcRHSHistAndSource(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  fac,         //!< domain-integration factor
      const double                  rhsint       //!< rhs at Gauss point
    );

    //! standard Galerkin convective term on right hand side
    void CalcRHSConv(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  conv_phi     //!< convective contribution
    );

    //! standard Galerkin diffusive term on right hand side
    void CalcRHSDiff(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      Teuchos::RCP<ScaTraEleDiffManager>  diffmanager,  //!< diffusion manager
      const LINALG::Matrix<nsd_,1>& gradphi      //!< scalar gradient at Gauss point
    );

    //! transient, convective and diffusive stabilization terms on right hand side
    void CalcRHSTransConvDiffStab(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhstaufac,   //!< time-integration factor for rhs times domain-integration factor times tau
      const double                  densnp,     //!< density at time_(n+1)
      const double                  scatrares,   //!< residual of convection-diffusion-reaction eq
      const LINALG::Matrix<nen_,1>& conv,        //!< convective operator
      const LINALG::Matrix<nen_,1>& sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>& diff         //!< laplace term
    );

    //! reactive terms (standard Galerkin and stabilization) on rhs
    void CalcRHSReact(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  rhstaufac,   //!< time-integration factor for rhs times domain-integration factor times tau
      const double                  rea_phi,     //!< reactive term
      const double                  densnp,      //!< density at time_(n+1)
      const double                  scatrares,   //!< residual of convection-diffusion-reaction eq
      Teuchos::RCP<ScaTraEleReaManager>  reamanager //!< reaction manager
    );

    //! fine-scale subgrid-diffusivity term on right hand side
    void CalcRHSFSSGD(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  sgdiff,      //!< subgrid-scale diffusivity
      const LINALG::Matrix<nsd_,1>  fsgradphi    //!< gardient of fine-scale velocity
    );

    //! multifractal subgrid-scale modeling on right hand side only rasthofer 11/13  |
    void CalcRHSMFS(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  densnp,      //!< density at time_(n+1)
      const LINALG::Matrix<nsd_,1>  convelint,   //!< convective velocity
      const LINALG::Matrix<nsd_,1>  mfsggradphi, //!< gradient of multifractal fine-scale scalar
      const LINALG::Matrix<nsd_,1>  mfsgvelint,  //!< multifractal fine-scale velocity
      const LINALG::Matrix<nsd_,1>  gradphi,     //!< scalar gradient
      const double                  phinp,       //!< scalar at time_(n+1)
      const double                  mfssgphi,    //!< multifractal fine-scale scalar
      const double                  mfsvdiv      //!< divergence of fine-scale velocity
    );

    /*========================================================================*/
    //! @name dofs and nodes
    /*========================================================================*/

    //! number of dof per node
    const int numdofpernode_;
    //! number of transported scalars (numscal_ <= numdofpernode_)
    const int numscal_;

    /*========================================================================*/
    //! @name parameter lists
    /*========================================================================*/

    //! pointer to parameter lists
    DRT::ELEMENTS::ScaTraEleParameter* scatrapara_;
    //! pointer to parameter list for time integration
    DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparatimint_;

    /*========================================================================*/
    //! @name manager classes for efficient application to various problems
    /*========================================================================*/

    //! manager for diffusion
    Teuchos::RCP<ScaTraEleDiffManager> diffmanager_;
    //! manager for reaction
    Teuchos::RCP<ScaTraEleReaManager> reamanager_;

    /*========================================================================*/
    //! @name scalar degrees of freedom and related
    /*========================================================================*/

    //! scalar at t_(n)
    std::vector<LINALG::Matrix<nen_,1> > ephin_;
    //! scalar at t_(n+1) or t_(n+alpha_F)
    std::vector<LINALG::Matrix<nen_,1> > ephinp_;
    //! history vector of transported scalars
    std::vector<LINALG::Matrix<nen_,1> > ehist_;
    //! fine-scale solution ?
    std::vector<LINALG::Matrix<nen_,1> > fsphinp_;

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
    //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nen_,1> eprenp_;


    /*========================================================================*/
    //! @name element coefficients and related extracted in evaluate
    /*========================================================================*/

    //! turbulent Prandtl number
    double tpn_;

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

    //! bodyforce in element nodesx
    std::vector<LINALG::Matrix<nen_,1> > bodyforce_;
//
    //! weights for nurbs elements
    LINALG::Matrix<nen_,1> weights_;
    //! knot vector for nurbs elements
    std::vector<Epetra_SerialDenseVector> myknots_;

    /*========================================================================*/
    //! @name can be very useful
    /*========================================================================*/

    double eid_;

  };


  /// Scatra diffusion manager
  /*!
      This is a basic class to handle diffusion. It exclusively contains
      the isotropic diffusion coefficient. For anisotropic diffusion or
      more advanced diffusion laws, e.g., nonlinear ones, a derived class
      has to be constructed in the problem-dependent subclass for element
      evaluation.
  */
  class ScaTraEleDiffManager
  {
  public:

    ScaTraEleDiffManager(int numscal):
      diff_(numscal,0.0),
      sgdiff_(numscal,0.0)
      {return;}

    virtual ~ScaTraEleDiffManager(){};

    //! Set the isotropic diffusion coefficient
    virtual void SetIsotropicDiff(const double diff, const int k)
    {
      if (diff < 0.0) dserror("negative (physical) diffusivity: %f",0,diff);

      diff_[k] = diff;
      return;
    }

    //! Set the isotropic subgrid diffusion coefficient
    virtual void SetIsotropicSubGridDiff(const double sgdiff, const int k)
    {
      diff_[k] += sgdiff;
      sgdiff_[k] += sgdiff;
      return;
    }

    //! Return the stored isotropic diffusion coefficient
    virtual double GetIsotropicDiff(const int k){return diff_[k];}

    //! Return the stored sub-grid diffusion coefficient
    virtual double GetSubGrDiff(const int k){return sgdiff_[k];}

  protected:

    //! scalar diffusion coefficient
    std::vector<double> diff_;

    //! subgrid diffusion coefficient
    std::vector<double> sgdiff_;
  };


  /// Scatra reaction manager
  /*!
      This is a basic class to handle reaction. It exclusively contains
      the a constant reaction coefficient. For more advanced reaction laws,
      e.g., nonlinear ones, a derived class has to be constructed in the
      problem-dependent subclass for element evaluation.
  */
  class ScaTraEleReaManager
    {
    public:

      ScaTraEleReaManager(int numscal):
        reacoeff_(numscal,0.0),
        reacoeffderivmatrix_(numscal,std::vector<double>(numscal,0.0 )), // size of matrix + initialized to zero
        reabodyforce_(numscal,0.0), // size of vector + initialized to zero
        reabodyforcederiv_(numscal,std::vector<double>(numscal,0.0 )),
        include_me_(false)
        {return;}

      virtual ~ScaTraEleReaManager(){};

      //! @name set routines

      //! Set the reaction coefficient
      virtual void SetReaCoeff(const double reacoeff, const int k)
      {
        if (reacoeff < 0.0) dserror("Reaction coefficient is not positive: %f",0,reacoeff);

        reacoeff_[k] = reacoeff;
        if( reacoeff!=0.0 )
          include_me_=true;

        return;
      }

      virtual void SetReaCoeffDerivMatrix(const double reacoeffderivmatrix, const int k, const int j)
      {
        (reacoeffderivmatrix_[k])[j] = reacoeffderivmatrix;
        return;
      }

      //! Set the derivative of reaction coefficient
      virtual void SetReaBodyForce(const double reabodyforce, const int k)
      {
        reabodyforce_[k] = reabodyforce;
        return;
      }

      virtual void SetReaBodyForceDerivMatrix(const double reabodyforcederiv, const int k, const int j)
      {
        (reabodyforcederiv_[k])[j] = reabodyforcederiv;
        return;
      }

      //! @name access routines

      //! Return the reaction coefficient
      virtual double GetReaCoeff(const int k){return reacoeff_[k];}

      //! Return the reaction coefficient
      virtual double GetReaCoeffDerivMatrix(const int k , const int j){return (reacoeffderivmatrix_[k])[j];}

      //! Return the reaction coefficient
      virtual double GetReaBodyForce(const int k){return reabodyforce_[k];}

      //! Return the reaction coefficient
      virtual double GetReaBodyForceDerivMatrix(const int k, const int j){return (reabodyforcederiv_[k])[j];}

      //! return flag: reaction activated
      virtual bool Active(){return include_me_;}

    private:

      //! scalar reaction coefficient
      std::vector<double> reacoeff_;

      //! reaction derivative for coupled problems
      std::vector<std::vector<double> > reacoeffderivmatrix_;

      //! scalar reaction coefficient
      std::vector<double> reabodyforce_;

      //! scalar reaction coefficient
      std::vector<std::vector<double> > reabodyforcederiv_;

      //! flag for reaction
      bool include_me_;
    };

  }
}


#endif
