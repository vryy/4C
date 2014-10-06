/*!----------------------------------------------------------------------
\file scatra_ele_calc.cpp

\brief main file containing routines for calculation of scatra element

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/


#include "scatra_ele_calc.H"

#include "scatra_ele.H"

#include "scatra_ele_parameter.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_globalproblem.H"  // for time curve in body force
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"


/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(const int numdofpernode, const int numscal)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    ephin_(numscal_),   // size of vector
    ephinp_(numscal_),  // size of vector
    ehist_(numscal_),   // size of vector
    fsphinp_(numscal_), // size of vector
    evelnp_(true),      // initialized to zero
    econvelnp_(true),   // initialized to zero
    efsvel_(true),      // initialized to zero
    eaccnp_(true),      // initialized to zero
    edispnp_(true),     // initialized to zero
    eprenp_(true),      // initialized to zero
    tpn_(0.0),
    xsi_(true),     // initialized to zero
    xyze_(true),    // initialized to zero
    funct_(true),   // initialized to zero
    deriv_(true),   // initialized to zero
    deriv2_(true),  // initialized to zero
    derxy_(true),   // initialized to zero
    derxy2_(true),  // initialized to zero
    xjm_(true),     // initialized to zero
    xij_(true),     // initialized to zero
    xder2_(true),   // initialized to zero
    bodyforce_(numdofpernode_), // size of vector
    weights_(true),      // initialized to zero
    myknots_(nsd_),      // size of vector
    eid_(0)
{
  // get parameter lists
  scatrapara_ = DRT::ELEMENTS::ScaTraEleParameterStd::Instance();
  scatraparatimint_ = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance();

  // get diffusion manager for diffusivity / diffusivities (in case of systems)
  //                           or (thermal conductivity/specific heat) in case of loma
  diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManager(numscal_));
  // get diffusion manager
  reamanager_ = Teuchos::rcp(new ScaTraEleReaManager(numscal_));

  return;
}


/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalc<distype>::Evaluate(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra
  )
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  //get element coordinates
  ReadElementCoordinatesAndProject(ele);

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,myknots_,weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractElementAndNodeValues(ele,params,discretization,lm);

  //--------------------------------------------------------------------------------
  // prepare turbulence models
  //--------------------------------------------------------------------------------

  int nlayer = 0;
  ExtractTurbulenceApproach(ele,params,discretization,lm,nlayer);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  Sysmat(
    ele,
    elemat1_epetra,
    elevec1_epetra,
    elevec2_epetra);

#if 0
  // for debugging of matrix entries
  if((ele->Id()==2) and (time < 1.3 and time > 1.1))
  {
    FDcheck(
      ele,
      elemat1_epetra,
      elevec1_epetra,
      elevec2_epetra,
      time,
      dt,
      timefac,
      alphaF,
      whichassgd,
      whichfssgd,
      assgd,
      fssgd,
      turbmodel_,
      Cs,
      tpn,
      frt,
      scatratype);
  }
#endif

  // ---------------------------------------------------------------------
  // output values of Prt, diffeff and Cs_delta_sq_Prt (channel flow only)
  // ---------------------------------------------------------------------

  if (scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky and scatrapara_->CsAv())
  {
    Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
    StoreModelParametersForOutput(ele,ele->Owner() == discretization.Comm().MyPID(),turbulencelist,nlayer);
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | extract element based or nodal values                     ehrl 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double>  DRT::ELEMENTS::ScaTraEleCalc<distype>::ExtractElementAndNodeValues(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm
)
{
  // set element id
  eid_ = ele->Id();

  // get convective (velocity - mesh displacement) velocity at nodes
  const Teuchos::RCP<Epetra_MultiVector> convelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field");
  DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,convelocity,nsd_);

  if (scatrapara_->IsAle())
  {
    // get velocity at nodes
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,velocity,nsd_);

    const Teuchos::RCP<Epetra_MultiVector> dispnp = params.get< Teuchos::RCP<Epetra_MultiVector> >("dispnp");
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_);
    // add nodal displacements to point coordinates
    xyze_ += edispnp_;
  }
  else
  {
    edispnp_.Clear();

    // velocity = convective velocity for the non-ale case
    evelnp_ = econvelnp_;
  }

  // get data required for subgrid-scale velocity: acceleration and pressure
  if (scatrapara_->RBSubGrVel())
  {
    const Teuchos::RCP<Epetra_MultiVector> accpre = params.get< Teuchos::RCP<Epetra_MultiVector> >("acceleration/pressure field");
    LINALG::Matrix<nsd_+1,nen_> eaccprenp;
    DRT::UTILS::ExtractMyNodeBasedValues(ele,eaccprenp,accpre,nsd_+1);

    // split acceleration and pressure values
    for (int i=0;i<nen_;++i)
    {
      for (int j=0;j<nsd_;++j)
      {
        eaccnp_(j,i) = eaccprenp(j,i);
      }
      eprenp_(i) = eaccprenp(nsd_,i);
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist==Teuchos::null || phinp==Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'");
  std::vector<double> myhist(lm.size());
  std::vector<double> myphinp(lm.size());
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

  // fill all element arrays
  for (int i=0;i<nen_;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each transported scalar, insert into element arrays
      ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
    }
    for (int k = 0; k< numscal_; ++k)
    {
      // the history vectors contains information of time step t_n
      ehist_[k](i,0) = myhist[k+(i*numdofpernode_)];
    }
  } // for i

  if (scatraparatimint_->IsGenAlpha() and not scatraparatimint_->IsIncremental())
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    if (phin==Teuchos::null) dserror("Cannot get state vector 'phin'");
    std::vector<double> myphin(lm.size());
    DRT::UTILS::ExtractMyValues(*phin,myphin,lm);

    // fill element array
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephin_[k](i,0) = myphin[k+(i*numdofpernode_)];
      }
    } // for i
  }

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  OtherNodeBasedSourceTerms(lm,discretization,params);

  // return extracted values of phinp
  return myphinp;
}


/*----------------------------------------------------------------------*
 | extract turbulence approach                          rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::ExtractTurbulenceApproach(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm,
  int&                       nlayer
)
{
  if (scatrapara_->TurbModel()!=INPAR::FLUID::no_model or (scatraparatimint_->IsIncremental() and scatrapara_->FSSGD()))
  {
    // do some checks first
    if (numscal_!=1 or numdofpernode_!=1)
      dserror("For the time being, turbulence approaches only support one scalar field!");
  }

  // set turbulent Prandt number to value given in parameterlist
  tpn_ = scatrapara_->TPN();

  // if we have a dynamic model,we overwrite this value by a local element-based one here
  if (scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
  {
    Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
    // remark: for dynamic estimation, this returns (Cs*h)^2 / Pr_t
    Teuchos::RCP<Epetra_Vector> ele_prt = turbulencelist.get<Teuchos::RCP<Epetra_Vector> >("col_ele_Prt");
    const int id = ele->LID();
    tpn_ = (*ele_prt)[id];

    // when no averaging was done, we just keep the calculated (clipped) value
    if (scatrapara_->CsAv())
      GetMeanPrtOfHomogenousDirection(params.sublist("TURBULENCE MODEL"),nlayer);
  }

  // get fine-scale values
  if ((scatraparatimint_->IsIncremental() and
      (scatrapara_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_all or scatrapara_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_small))
      or scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // get fine scale scalar field
    Teuchos::RCP<const Epetra_Vector> gfsphinp = discretization.GetState("fsphinp");
    if (gfsphinp==Teuchos::null) dserror("Cannot get state vector 'fsphinp'");

    std::vector<double> myfsphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*gfsphinp,myfsphinp,lm);

    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        fsphinp_[k](i,0) = myfsphinp[k+(i*numdofpernode_)];
      }
    }

    // get fine-scale velocity at nodes
    if (scatrapara_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_small or scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      const Teuchos::RCP<Epetra_MultiVector> fsvelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("fine-scale velocity field");
      DRT::UTILS::ExtractMyNodeBasedValues(ele,efsvel_,fsvelocity,nsd_);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 g.bau 08/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::Sysmat(
  DRT::Element*                         ele, ///< the element whose matrix is calculated
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs, ///< element rhs to calculate
  Epetra_SerialDenseVector&             subgrdiff ///< subgrid-diff.-scaling vector
  )
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  const double vol=EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center are also necessary
  // even if the stabilization parameter is evaluated at the element center
  if (not scatrapara_->MatGP())
    GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);

  //----------------------------------------------------------------------
  // calculation of subgrid diffusivity and stabilization parameter(s)
  // at element center
  //----------------------------------------------------------------------

  // the stabilization parameters (one per transported scalar)
  std::vector<double> tau(numscal_,0.0);
  // subgrid-scale diffusion coefficient
  double sgdiff(0.0);

  if (not scatrapara_->TauGP())
  {
    // get velocity at element center
    LINALG::Matrix<nsd_,1> convelint(true);
    convelint.Multiply(econvelnp_,funct_);

    for (int k = 0;k<numscal_;++k) // loop of each transported scalar
    {
      // calculation of all-scale subgrid diffusivity (by, e.g.,
      // Smagorinsky model) at element center
      if (scatrapara_->TurbModel() == INPAR::FLUID::smagorinsky
          or scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky
          or scatrapara_->TurbModel() == INPAR::FLUID::dynamic_vreman)
      {
        CalcSubgrDiff(diffmanager_,visc,vol,k,densnp);
      }

      // calculation of fine-scale artificial subgrid diffusivity at element center
      if (scatrapara_->FSSGD()) CalcFineScaleSubgrDiff(sgdiff,subgrdiff,ele,vol,k,densnp,diffmanager_->GetIsotropicDiff(k),convelint);

      // calculation of stabilization parameter at element center
      CalcTau(tau[k],diffmanager_->GetIsotropicDiff(k),reamanager_->GetReaCoeff(k),densnp,convelint,vol);
    }
  }

  // prepare multifractal subgrid-scale modeling
  // calculation of model coefficients B (velocity) and D (scalar)
  // at element center
  // coefficient B of fine-scale velocity
  LINALG::Matrix<nsd_,1> B_mfs(true);
  // coefficient D of fine-scale scalar
  double D_mfs = 0.0;
  if (scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    if (not scatrapara_->BD_Gp())
    {
      // make sure to get material parameters at element center
      // hence, determine them if not yet available
      if (scatrapara_->MatGP()) GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);
      // provide necessary velocities and gradients at element center
      // get velocity at element center
      LINALG::Matrix<nsd_,1> convelint(true);
      LINALG::Matrix<nsd_,1> fsvelint(true);
      convelint.Multiply(econvelnp_,funct_);
      fsvelint.Multiply(efsvel_,funct_);

      // calculate model coefficients
      for (int k = 0;k<numscal_;++k) // loop of each transported scalar
        CalcBAndDForMultifracSubgridScales(B_mfs,D_mfs,vol,k,densnp,diffmanager_->GetIsotropicDiff(k),visc,convelint,fsvelint);
    }
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP())
      GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc,iquad);

    // get velocity at integration point
    //LINALG::Matrix<nsd_,1> velint(true);
    LINALG::Matrix<nsd_,1> convelint(true);
    //velint.Multiply(evelnp_,funct_);
    convelint.Multiply(econvelnp_,funct_);

    // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
    LINALG::Matrix<nen_,1> conv(true);
    conv.MultiplyTN(derxy_,convelint);

    // velocity divergence required for conservative form
    double vdiv(0.0);
    if (scatrapara_->IsConservative()) GetDivergence(vdiv,evelnp_);

    // get fine-scale velocity and its derivatives at integration point
    LINALG::Matrix<nsd_,1> fsvelint(true);
    if (scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
      fsvelint.Multiply(efsvel_,funct_);

    // loop all scalars
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // scalar at integration point at time step n+1
      const double phinp = funct_.Dot(ephinp_[k]);
      // scalar at integration point at time step n
      const double phin = funct_.Dot(ephin_[k]);

      // gradient of current scalar value at integration point
      LINALG::Matrix<nsd_,1> gradphi(true);
      gradphi.Multiply(derxy_,ephinp_[k]);

      // convective term using current scalar value
      double conv_phi(0.0);
      conv_phi = convelint.Dot(gradphi);

      // diffusive part used in stabilization terms
      double diff_phi(0.0);
      LINALG::Matrix<nen_,1> diff(true);
      // diffusive term using current scalar value for higher-order elements
      if (use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff);
        diff.Scale(diffmanager_->GetIsotropicDiff(k));
        diff_phi = diff.Dot(ephinp_[k]);
      }

      // reactive part of the form: (reaction coefficient)*phi
      double rea_phi(0.0);
      rea_phi = densnp*phinp*reamanager_->GetReaCoeff(k);

      // compute gradient of fine-scale part of scalar value
      LINALG::Matrix<nsd_,1> fsgradphi(true);
      if (scatrapara_->FSSGD())
        fsgradphi.Multiply(derxy_,fsphinp_[k]);

      // get history data (or acceleration)
      double hist(0.0);
      hist = funct_.Dot(ehist_[k]);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      GetRhsInt(rhsint,densnp,k);

      //--------------------------------------------------------------------
      // calculation of (fine-scale) subgrid diffusivity, subgrid-scale
      // velocity and stabilization parameter(s) at integration point
      //--------------------------------------------------------------------

      // subgrid-scale convective term
      LINALG::Matrix<nen_,1> sgconv(true);
      // subgrid-scale velocity vector in gausspoint
      LINALG::Matrix<nsd_,1> sgvelint(true);

      if (scatrapara_->TauGP())
      {
        // artificial diffusion / shock capturing: adaption of diffusion coefficient
        if (scatrapara_->ASSGD())
        {
          // residual of convection-diffusion-reaction eq
          double scatrares(0.0);
          // residual-based subgrid-scale scalar (just a dummy here)
          double sgphi(0.0);

          // compute residual of scalar transport equation
          // (subgrid-scale part of scalar, which is also computed, not required)
          CalcResidualAndSubgrScalar(k,scatrares,sgphi,densam,densnp,phinp,hist,conv_phi,diff_phi,rea_phi,rhsint,tau[k]);

          // pre-calculation of stabilization parameter at integration point need for some forms of artificial diffusion
          CalcTau(tau[k],diffmanager_->GetIsotropicDiff(k),reamanager_->GetReaCoeff(k),densnp,convelint,vol);

          // compute artificial diffusion
          CalcArtificialDiff(vol,k,diffmanager_,densnp,convelint,gradphi,conv_phi,scatrares,tau[k]);

          // adapt diffusive term using current scalar value for higher-order elements,
          // since diffus -> diffus + sgdiff
          if (use2ndderiv_)
          {
            // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
            diff.Clear();
            GetLaplacianStrongForm(diff);
            diff.Scale(diffmanager_->GetIsotropicDiff(k));
            diff_phi = diff.Dot(ephinp_[k]);
          }
        }

      // calculation of all-scale subgrid diffusivity (by, e.g.,
      // Smagorinsky model) at element center
      if (scatrapara_->TurbModel() == INPAR::FLUID::smagorinsky
          or scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky
          or scatrapara_->TurbModel() == INPAR::FLUID::dynamic_vreman)
      {
        CalcSubgrDiff(diffmanager_,visc,vol,k,densnp);

        // adapt diffusive term using current scalar value for higher-order elements,
        // since diffus -> diffus + sgdiff
        if (use2ndderiv_)
        {
          // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
          diff.Clear();
          GetLaplacianStrongForm(diff);
          diff.Scale(diffmanager_->GetIsotropicDiff(k));
          diff_phi = diff.Dot(ephinp_[k]);
        }
      }

      // calculation of fine-scale artificial subgrid diffusivity at element center
      if (scatrapara_->FSSGD()) CalcFineScaleSubgrDiff(sgdiff,subgrdiff,ele,vol,k,densnp,diffmanager_->GetIsotropicDiff(k),convelint);

        // calculation of subgrid-scale velocity at integration point if required
        if (scatrapara_->RBSubGrVel())
        {
          // calculation of stabilization parameter related to fluid momentum
          // equation at integration point
          CalcTau(tau[k],visc,0.0,densnp,convelint,vol);
          // calculation of residual-based subgrid-scale velocity
          CalcSubgrVelocity(ele,sgvelint,densam,densnp,visc,convelint,tau[k]);

          // calculation of subgrid-scale convective part
          sgconv.MultiplyTN(derxy_,sgvelint);
        }

        // calculation of stabilization parameter at integration point
        CalcTau(tau[k],diffmanager_->GetIsotropicDiff(k),reamanager_->GetReaCoeff(k),densnp,convelint,vol);
      }

      // residual of convection-diffusion-reaction eq
      double scatrares(0.0);
      // residual-based subgrid-scale scalar (just a dummy here)
      double sgphi(0.0);

      // compute residual of scalar transport equation and
      // subgrid-scale part of scalar
      CalcResidualAndSubgrScalar(k,scatrares,sgphi,densam,densnp,phinp,hist,conv_phi,diff_phi,rea_phi,rhsint,tau[k]);

      // prepare multifractal subgrid-scale modeling
      // calculation of model coefficients B (velocity) and D (scalar)
      // at Gauss point as well as calculation
      // of multifractal subgrid-scale quantities
      LINALG::Matrix<nsd_,1> mfsgvelint(true);
      double mfsvdiv(0.0);
      double mfssgphi(0.0);
      LINALG::Matrix<nsd_,1> mfsggradphi(true);
      if (scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
      {
        if (scatrapara_->BD_Gp())
          // calculate model coefficients
          CalcBAndDForMultifracSubgridScales(B_mfs,D_mfs,vol,k,densnp,diffmanager_->GetIsotropicDiff(k),visc,convelint,fsvelint);

        // calculate fine-scale velocity, its derivative and divergence for multifractal subgrid-scale modeling
        for (int idim=0; idim<nsd_; idim++)
          mfsgvelint(idim,0) = fsvelint(idim,0) * B_mfs(idim,0);
        // required for conservative formulation in the context of passive scalar transport
        if (scatrapara_->MfsConservative() or scatrapara_->IsConservative())
        {
          // get divergence of subgrid-scale velocity
          LINALG::Matrix<nsd_,nsd_> mfsvderxy;
          mfsvderxy.MultiplyNT(efsvel_,derxy_);
          for (int idim = 0; idim<nsd_; idim++)
            mfsvdiv += mfsvderxy(idim,idim) * B_mfs(idim,0);
        }

        // calculate fine-scale scalar and its derivative for multifractal subgrid-scale modeling
        mfssgphi = D_mfs * funct_.Dot(fsphinp_[k]);
        mfsggradphi.Multiply(derxy_,fsphinp_[k]);
        mfsggradphi.Scale(D_mfs);
      }

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      // stabilization parameter and integration factors
      const double taufac     = tau[k]*fac;
      const double timefacfac = scatraparatimint_->TimeFac()*fac;
      const double timetaufac = scatraparatimint_->TimeFac()*taufac;

      //----------------------------------------------------------------
      // 1) element matrix: stationary terms
      //----------------------------------------------------------------

      // calculation of convective element matrix in convective form
      CalcMatConv(emat,k,timefacfac,densnp,conv,sgconv);

      // add conservative contributions
      if (scatrapara_->IsConservative())
        CalcMatConvAddCons(emat,k,timefacfac,vdiv,densnp);

      // calculation of diffusive element matrix
      CalcMatDiff(emat,k,timefacfac,diffmanager_);

      //----------------------------------------------------------------
      // convective stabilization term
      //----------------------------------------------------------------

      // convective stabilization of convective term (in convective form)
      // transient stabilization of convective term (in convective form)
      if(scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
        CalcMatTransConvDiffStab(emat,k,timetaufac,densnp,conv,sgconv,diff);

      //----------------------------------------------------------------
      // 2) element matrix: instationary terms
      //----------------------------------------------------------------

      if (not scatraparatimint_->IsStationary())
      {
        CalcMatMass(emat,k,fac,densam);

        if(scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
          CalcMatMassStab(emat,k,taufac,densam,densnp,conv,sgconv,diff);
      }

      //----------------------------------------------------------------
      // 3) element matrix: reactive term
      //----------------------------------------------------------------

      // including stabilization
      // if (reamanager_->Active())
        CalcMatReact(emat,k,timefacfac,timetaufac,taufac,densnp,phinp,reamanager_,conv,sgconv,diff);

      //----------------------------------------------------------------
      // 4) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac    = scatraparatimint_->TimeFacRhs() * fac;
      double rhstaufac = scatraparatimint_->TimeFacRhsTau() * taufac;

      if (scatraparatimint_->IsIncremental() and not scatraparatimint_->IsStationary())
        CalcRHSLinMass(erhs,k,rhsfac,fac,densam,densnp,phinp,hist);

      // the order of the following three functions is important
      // and must not be changed
      ComputeRhsInt(rhsint,densam,densnp,hist);

      RecomputeScatraResForRhs(scatrares,k,convelint,gradphi,diff,densn,densnp,conv_phi,rea_phi,phin,reamanager_,rhsint);

      RecomputeConvPhiForRhs(conv_phi,k,sgvelint,gradphi,densnp,densn,phinp,phin,vdiv);

      //----------------------------------------------------------------
      // standard Galerkin transient, old part of rhs and bodyforce term
      //----------------------------------------------------------------
      CalcRHSHistAndSource(erhs,k,fac,rhsint);

      //----------------------------------------------------------------
      // standard Galerkin terms on right hand side
      //----------------------------------------------------------------

      // convective term
      CalcRHSConv(erhs,k,rhsfac,conv_phi);

      // diffusive term
      CalcRHSDiff(erhs,k,rhsfac,diffmanager_,gradphi);

      //----------------------------------------------------------------
      // stabilization terms
      //----------------------------------------------------------------
      if (scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
        CalcRHSTransConvDiffStab(erhs,k,rhstaufac,densnp,scatrares,conv,sgconv,diff);

      //----------------------------------------------------------------
      // reactive terms (standard Galerkin and stabilization) on rhs
      //----------------------------------------------------------------

      if (reamanager_->Active())
        CalcRHSReact(erhs,k,rhsfac,rhstaufac,rea_phi,densnp,scatrares,reamanager_);

      //----------------------------------------------------------------
      // 5) advanced turbulence models
      //----------------------------------------------------------------

      //----------------------------------------------------------------
      // fine-scale subgrid-diffusivity term on right hand side
      //----------------------------------------------------------------
      if (scatraparatimint_->IsIncremental() and scatrapara_->FSSGD())
        CalcRHSFSSGD(erhs,k,rhsfac,sgdiff,fsgradphi);

      //---------------------------------------------------------------
      // multifractal subgrid-scale modeling on right hand side only
      //---------------------------------------------------------------
      if (scatraparatimint_->IsIncremental() and scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
        CalcRHSMFS(erhs,k,rhsfac,densnp,convelint,mfsggradphi,mfsgvelint,gradphi,phinp,mfssgphi,mfsvdiv);

    }// end loop all scalars

  }// end loop Gauss points

  return;
}


/*----------------------------------------------------------------------*
  |  get the body force  (private)                              gjb 06/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::BodyForce(
  const DRT::Element*    ele
  )
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique Neumann condition
  switch(nsd_)
  {
  case 3:
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
    break;
  case 2:
    DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
    break;
  case 1:
    DRT::UTILS::FindElementConditions(ele, "LineNeumann", myneumcond);
    break;
  default:
    dserror("Illegal number of spatial dimensions: %d",nsd_);
    break;
  }

  if (myneumcond.size()>1)
    dserror("More than one Neumann condition on one node!");

  if (myneumcond.size()==1)
  {
    // check for potential time curve
    const std::vector<int>* curve  = myneumcond[0]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];

    // initialization of time-curve factor
    double curvefac(0.0);

    // compute potential time curve or set time-curve factor to one
    if (curvenum >= 0)
    {
      // time factor (negative time indicating error)
      if (scatraparatimint_->Time() >= 0.0)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(scatraparatimint_->Time());
      else dserror("Negative time in bodyforce calculation: time = %f",scatraparatimint_->Time());
    }
    else curvefac = 1.0;

    // get values and switches from the condition
    const std::vector<int>*    onoff = myneumcond[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = myneumcond[0]->Get<std::vector<double> >("val"  );

    // set this condition to the bodyforce array
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      for (int jnode=0; jnode<nen_; jnode++)
      {
        (bodyforce_[idof])(jnode) = (*onoff)[idof]*(*val)[idof]*curvefac;
      }
    }
  }
  else
  {
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      // no bodyforce
      bodyforce_[idof].Clear();
    }
  }

  return;

} //ScaTraEleCalc::BodyForce


/*------------------------------------------------------------------------*
 | further node-based source terms not given via Neumann volume condition |
 |                                                        rasthofer 12/13 |
 *------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::OtherNodeBasedSourceTerms(
  const std::vector<int>&    lm,
  DRT::Discretization&       discretization,
  Teuchos::ParameterList&    params
)
{
  // set externally calculated source term instead of body force by volume
  // Neumann boundary condition of input file
  if (scatrapara_->ScalarForcing()==INPAR::FLUID::scalarforcing_isotropic)
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> source = discretization.GetState("forcing");
    std::vector<double> mysource(lm.size());
    DRT::UTILS::ExtractMyValues(*source,mysource,lm);

    // fill element array
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numdofpernode_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        bodyforce_[k](i,0) = mysource[k+(i*numdofpernode_)];
      }
    } // for i
  }
  // special forcing mean scalar gradient
  else if (scatrapara_->ScalarForcing()==INPAR::FLUID::scalarforcing_mean_scalar_gradient)
  {
    // get mean-scalar gradient
    const double grad_phi = params.sublist("TURBULENCE MODEL").get<double>("MEAN_SCALAR_GRADIENT");

    // fill element array
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numdofpernode_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        bodyforce_[k](i,0) = -grad_phi*evelnp_(2,i);
      }
    } // for i
  }

  return;
}

/*----------------------------------------------------------------------*
 | read element coordinates, assuming they are all 3D and then project to
 | the respective lower dimensional space   bertoglio 08/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::ReadElementCoordinatesAndProject(
    const DRT::ELEMENTS::Transport*  ele
    )
{
  if(nsd_ == 1){
    double lengthLine3D;
    LINALG::Matrix<3,nen_> xyze3D;
    GEO::fillInitialPositionArray<distype,3,LINALG::Matrix<3,nen_> >(ele,xyze3D);
    xyze_(0,0)=0.0;
    for(int i=1;i<nen_;i++){
      lengthLine3D = sqrt( (xyze3D(0,i-1)-xyze3D(0,i))*(xyze3D(0,i-1)-xyze3D(0,i)) +
                              (xyze3D(1,i-1)-xyze3D(1,i))*(xyze3D(1,i-1)-xyze3D(1,i)) +
                              (xyze3D(2,i-1)-xyze3D(2,i))*(xyze3D(2,i-1)-xyze3D(2,i)));
      xyze_(0,i)=xyze_(0,i-1)+lengthLine3D;
    }
  }
 // else if(nsd_ == 2)
 // {
    // TO DO
//  }
  else
  {
    // Directly copy the coordinates since in 3D the transformation is just the identity
    GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  }
  return;
} //ScaTraEleCalc::ReadElementCoordinatesAndProject

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at ele. center   ehrl 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::ScaTraEleCalc<distype>::EvalShapeFuncAndDerivsAtEleCenter()
{
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0);

  return vol;

} //ScaTraImpl::EvalShapeFuncAndDerivsAtEleCenter


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     gjb 08/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalc<distype>::EvalShapeFuncAndDerivsAtIntPoint(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
  const int                                    iquad       ///< id of current Gauss point
  )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  if (not DRT::NURBS::IsNurbs(distype))
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
    if (use2ndderiv_)
    {
      // get the second derivatives of standard element at current GP
      DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
    }
  }
  else // nurbs elements are always somewhat special...
  {
    if (use2ndderiv_)
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv_deriv2
        (funct_  ,
         deriv_  ,
         deriv2_ ,
         xsi_    ,
         myknots_,
         weights_,
         distype );
    }
    else
    {
      DRT::NURBS::UTILS::nurbs_get_funct_deriv
        (funct_  ,
         deriv_  ,
         xsi_    ,
         myknots_,
         weights_,
         distype );
    }
  } // IsNurbs()

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */

  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad]*det;

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  // compute second global derivatives (if needed)
  if (use2ndderiv_)
  {
    // get global second derivatives
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else
    derxy2_.Clear();

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

} //ScaTraImpl::EvalShapeFuncAndDerivsAtIntPoint


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                      gjb 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>  reamanager,   //!< reaction manager
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
// get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

// get diffusivity / diffusivities
  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
    }
  }
  else
    Materials(material,0,densn,densnp,densam,diffmanager,reamanager,visc,iquad);

  return;
} //ScaTraEleCalc::GetMaterialParams


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point

  )
{
  if (material->MaterialType() == INPAR::MAT::m_scatra)
    MatScaTra(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
  else dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad   //!< id of current gauss point (default = -1)
  )
{

  int leleid = -1;
  if(DRT::Problem::Instance()->ProblemType()==prb_acou) leleid = DRT::Problem::Instance()->GetDis("scatra")->ElementColMap()->LID(eid_);

  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  // get constant diffusivity
  diffmanager->SetIsotropicDiff(actmat->Diffusivity(leleid),k);

  // get reaction coefficient
  reamanager->SetReaCoeff(actmat->ReaCoeff(leleid),k);

  // in case of multifractal subgrid-scales, read Schmidt number
  if (scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales or scatrapara_->RBSubGrVel()
      or scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
  {
    //access fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
    fluiddis = DRT::Problem::Instance()->GetDis("fluid");
    //get corresponding fluid element (it has the same global ID as the scatra element)
    DRT::Element* fluidele = fluiddis->gElement(eid_);
    if (fluidele == NULL)
      dserror("Fluid element %i not on local processor", eid_);

    // get fluid material
    Teuchos::RCP<MAT::Material> fluidmat = fluidele->Material();
    if(fluidmat->MaterialType() != INPAR::MAT::m_fluid)
      dserror("Invalid fluid material for passive scalar transport in turbulent flow!");

    const Teuchos::RCP<const MAT::NewtonianFluid>& actfluidmat = Teuchos::rcp_dynamic_cast<const MAT::NewtonianFluid>(fluidmat);

    // get constant dynamic viscosity
    visc = actfluidmat->Viscosity();
    densn = actfluidmat->Density();
    densnp = actfluidmat->Density();
    densam = actfluidmat->Density();

    if (densam != 1.0 or densnp != 1.0 or densn != 1.0)
       dserror("Check your diffusivity! Dynamic diffusivity required!");
   }

  return;
} // ScaTraEleCalc<distype>::MatScaTra


/*---------------------------------------------------------------------------------------*
 |  calculate the Laplacian in strong form for all shape functions (private)   gjb 04/10 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::GetLaplacianStrongForm(
  LINALG::Matrix<nen_,1>& diff
  )
{
  diff.Clear();
  // compute N,xx  +  N,yy +  N,zz for each shape function at integration point
  for (int i=0; i<nen_; ++i)
  {
    for (int j = 0; j<nsd_; ++j)
    {
      diff(i) += derxy2_(j,i);
    }
  }
  return;
} // ScaTraEleCalc<distype>::GetLaplacianStrongForm


/*-----------------------------------------------------------------------------*
 |  calculate divergence of vector field (e.g., velocity)  (private) gjb 04/10 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::GetDivergence(
  double&                          vdiv,
  const LINALG::Matrix<nsd_,nen_>& evel)
{
  LINALG::Matrix<nsd_,nsd_> vderxy;
  vderxy.MultiplyNT(evel,derxy_);

  vdiv = 0.0;
  // compute vel x,x  + vel y,y +  vel z,z at integration point
  for (int j = 0; j<nsd_; ++j)
  {
    vdiv += vderxy(j,j);
  }
  return;
} // ScaTraEleCalc<distype>::GetDivergence


/*-----------------------------------------------------------------------------*
 | compute rhs containing bodyforce                                 ehrl 11/13 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::GetRhsInt(
  double&      rhsint,  //!< rhs containing bodyforce at Gauss point
  const double densnp,  //!< density at t_(n+1)
  const int    k        //!< index of current scalar
  )
{
  // compute rhs containing bodyforce (divided by specific heat capacity) and,
  // for temperature equation, the time derivative of thermodynamic pressure,
  // if not constant, and for temperature equation of a reactive
  // equation system, the reaction-rate term
  rhsint = bodyforce_[k].Dot(funct_);

  return;
} // GetRhsInt


/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix in convective form     ehrl 11/13 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatConv(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  const double                  densnp,
  const LINALG::Matrix<nen_,1>& conv,
  const LINALG::Matrix<nen_,1>& sgconv
  )
{
  // convective term in convective form
  const double densfac = timefacfac*densnp;
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = densfac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*(conv(ui)+sgconv(ui));
    }
  }
  return;
} // ScaTraEleCalc<distype>::CalcMatConv


/*------------------------------------------------------------------------------------------*
 |  calculation of convective element matrix: add conservative contributions     ehrl 11/13 |
 *------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatConvAddCons(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  const double                  vdiv,
  const double                  densnp
  )
{
  const double consfac = timefacfac*densnp*vdiv;
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = consfac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*funct_(ui);
    }
  }

  return;
}


/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatDiff(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  Teuchos::RCP<ScaTraEleDiffManager>  diffmanager
  )
{
  // diffusive term
  const double fac_diffus = timefacfac * diffmanager->GetIsotropicDiff(k);
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf,ui,vi);
      emat(fvi,fui) += fac_diffus*laplawf;
    }
  }
  return;
}


/*------------------------------------------------------------------- *
 |  calculation of stabilization element matrix            ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatTransConvDiffStab(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timetaufac,
  const double                  densnp,
  const LINALG::Matrix<nen_,1>& conv,
  const LINALG::Matrix<nen_,1>& sgconv,
  const LINALG::Matrix<nen_,1>& diff
  )
{
  const double dens2taufac = timetaufac*densnp*densnp;
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = dens2taufac*(conv(vi)+sgconv(vi)+scatrapara_->USFEMGLSFac()*1.0/scatraparatimint_->TimeFac()*funct_(vi));
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*conv(ui);
    }
  }

  //----------------------------------------------------------------
  // stabilization terms for higher-order elements
  //----------------------------------------------------------------
  if (use2ndderiv_)
  {
    const double denstaufac = timetaufac*densnp;
    // convective stabilization of diffusive term (in convective form)
    // transient stabilization of diffusive term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = denstaufac*(conv(vi)+sgconv(vi)+scatrapara_->USFEMGLSFac()*1.0/scatraparatimint_->TimeFac()*funct_(vi));
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) -= v*diff(ui);
      }
    }

    const double densdifftaufac = scatrapara_->USFEMGLSFac()*denstaufac;
    // diffusive stabilization of convective term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densdifftaufac*diff(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) -= v*conv(ui);
      }
    }

    const double difftaufac = scatrapara_->USFEMGLSFac()*timetaufac;
    // diffusive stabilization of diffusive term
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = difftaufac*diff(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*diff(ui);
      }
    }
  }
  return;
}


/*------------------------------------------------------------------- *
 |  calculation of mass element matrix                    ehrl 11/13  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatMass(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  fac,
  const double                  densam
  )
{
  const double densamfac = fac*densam;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = densamfac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*funct_(ui);
    }
  }
  return;
}


/*------------------------------------------------------------------- *
 |  calculation of stabilization mass element matrix      ehrl 11/13  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatMassStab(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  taufac,
  const double                  densam,
  const double                  densnp,
  const LINALG::Matrix<nen_,1>& conv,
  const LINALG::Matrix<nen_,1>& sgconv,
  const LINALG::Matrix<nen_,1>& diff
  )
{
  const double densamnptaufac = taufac*densam*densnp;
  //----------------------------------------------------------------
  // stabilization of transient term
  //----------------------------------------------------------------
  // convective stabilization of transient term (in convective form)
  // transient stabilization of transient term
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = densamnptaufac*(conv(vi)+sgconv(vi)+scatrapara_->USFEMGLSFac()*1.0/scatraparatimint_->TimeFac()*funct_(vi));
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*funct_(ui);
    }
  }

  if (use2ndderiv_)
  {
    const double densamreataufac = scatrapara_->USFEMGLSFac()*taufac*densam;
    // diffusive stabilization of transient term
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densamreataufac*diff(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) -= v*funct_(ui);
      }
    }
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of reactive element matrix                ehrl 11/13  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcMatReact(
  Epetra_SerialDenseMatrix&          emat,
  const int                          k,
  const double                       timefacfac,
  const double                       timetaufac,
  const double                       taufac,
  const double                       densnp,
  const double                       phinp,
  Teuchos::RCP<ScaTraEleReaManager>  reamanager,
  const LINALG::Matrix<nen_,1>&      conv,
  const LINALG::Matrix<nen_,1>&      sgconv,
  const LINALG::Matrix<nen_,1>&      diff
  )
{
  const double fac_reac        = timefacfac*densnp*reamanager->GetReaCoeff(k);
  const double timetaufac_reac = timetaufac*densnp*reamanager->GetReaCoeff(k);

  //----------------------------------------------------------------
  // standard Galerkin reactive term
  //----------------------------------------------------------------
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = fac_reac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*funct_(ui);
    }
  }

  //----------------------------------------------------------------
  // stabilization of reactive term
  //----------------------------------------------------------------
  if(scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
  {
    double densreataufac = timetaufac_reac*densnp;
    // convective stabilization of reactive term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densreataufac*(conv(vi)+sgconv(vi)+scatrapara_->USFEMGLSFac()*1.0/scatraparatimint_->TimeFac()*funct_(vi));
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*funct_(ui);
      }
    }

    if (use2ndderiv_)
    {
      // diffusive stabilization of reactive term
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = scatrapara_->USFEMGLSFac()*timetaufac_reac*diff(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) -= v*funct_(ui);
        }
      }
    }

    //----------------------------------------------------------------
    // reactive stabilization
    //----------------------------------------------------------------
    densreataufac = scatrapara_->USFEMGLSFac()*timetaufac_reac*densnp;

    // reactive stabilization of convective (in convective form) and reactive term
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densreataufac*funct_(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*(conv(ui)+reamanager_->GetReaCoeff(k)*funct_(ui));
      }
    }

    if (use2ndderiv_)
    {
      // reactive stabilization of diffusive term
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = scatrapara_->USFEMGLSFac()*timetaufac_reac*funct_(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) -= v*diff(ui);
        }
      }
    }


    if (not scatraparatimint_->IsStationary())
    {
      // reactive stabilization of transient term
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = scatrapara_->USFEMGLSFac()*taufac*densnp*reamanager->GetReaCoeff(k)*densnp*funct_(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) += v*funct_(ui);
        }
      }

      if (use2ndderiv_ and reamanager->GetReaCoeff(k)!=0.0)
        dserror("Second order reactive stabilization is not fully implemented!! ");
    }
  }
 //}

  return;
}


/*------------------------------------------------------------------- *
 |  calculation of linearized mass rhs vector              ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSLinMass(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  const double                  fac,
  const double                  densam,
  const double                  densnp,
  const double                  phinp,
  const double                  hist
  )
{
  double vtrans = 0.0;

  if (scatraparatimint_->IsGenAlpha())
    vtrans = rhsfac*densam*hist;
  else
  {
    // compute scalar at integration point
    vtrans = fac*densnp*phinp;
  }

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    erhs[fvi] -= vtrans*funct_(vi);
  }

  return;
}


/*------------------------------------------------------------------- *
 | adaption of rhs with respect to time integration        ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::ComputeRhsInt(
  double&                       rhsint,
  const double                  densam,
  const double                  densnp,
  const double                  hist
  )
{
  if (scatraparatimint_->IsGenAlpha())
  {
    if (not scatraparatimint_->IsIncremental())
      rhsint   += densam*hist*(scatraparatimint_->AlphaF()/scatraparatimint_->TimeFac());

    rhsint   *= (scatraparatimint_->TimeFac()/scatraparatimint_->AlphaF());
  }
  else // OST, BDF2
  {
    if (not scatraparatimint_->IsStationary())
    {
      rhsint *= scatraparatimint_->TimeFac();
      rhsint += densnp*hist;
    }
  }

  return;
}


/*------------------------------------------------------------------- *
 | adaption of residual with respect to time integration   ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::RecomputeScatraResForRhs(
  double&                       scatrares,
  const int                     k,
  const LINALG::Matrix<nsd_,1>& convelint,
  LINALG::Matrix<nsd_,1>&       gradphi,
  const LINALG::Matrix<nen_,1>& diff,
  const double                  densn,
  const double                  densnp,
  double&                       conv_phi,
  double&                       rea_phi,
  const double                  phin,
  Teuchos::RCP<ScaTraEleReaManager> reamanager,
  const double                  rhsint
  )
{
  if (scatraparatimint_->IsGenAlpha())
  {
    if (not scatraparatimint_->IsIncremental())
    {
      // for this case, gradphi_ (i.e. the gradient
      // at time n+1) is overwritten by the gradient at time n
      // analogously, conv_phi_ at time n+1 is replace by its
      // value at time n
      // gradient of scalar value at n
      gradphi.Multiply(derxy_,ephin_[k]);

      // convective term using scalar value at n
      conv_phi = convelint.Dot(gradphi);

      // diffusive term using current scalar value for higher-order elements
      double diff_phin = 0.0;
      if (use2ndderiv_) diff_phin = diff.Dot(ephin_[k]);

      // reactive term using scalar value at n
      // if no reaction is chosen, GetReaCoeff(k) returns 0.0
      rea_phi = densnp*reamanager->GetReaCoeff(k)*phin;
      // reacterm_[k] must be evaluated at t^n to be used in the line above!

      scatrares = (1.0-scatraparatimint_->AlphaF()) * (densn*conv_phi
                    - diff_phin + rea_phi) - rhsint*scatraparatimint_->AlphaF()/scatraparatimint_->TimeFac();
    }
  }
  else if (scatraparatimint_->IsIncremental() and not scatraparatimint_->IsGenAlpha())
  {
    if (not scatraparatimint_->IsStationary())
      scatrares *= scatraparatimint_->Dt();
  }
  else
    scatrares = -rhsint;

  return;
}


/*------------------------------------------------------------------- *
 | adaption of convective term for rhs                     ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::RecomputeConvPhiForRhs(
  double&                       conv_phi,
  const int                     k,
  const LINALG::Matrix<nsd_,1>& sgvelint,
  const LINALG::Matrix<nsd_,1>& gradphi,
  const double                  densnp,
  const double                  densn,
  const double                  phinp,
  const double                  phin,
  const double                  vdiv
  )
{
  if (scatraparatimint_->IsIncremental())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (scatrapara_->IsConservative())
    {
      // convective term in conservative form
      conv_phi += phinp*vdiv;
    }

    // multiply convective term by density
    conv_phi *= densnp;
  }
  else if (not scatraparatimint_->IsIncremental() and scatraparatimint_->IsGenAlpha())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (scatrapara_->IsConservative())
    {
      // convective term in conservative form
      // caution: velocity divergence is for n+1 and not for n!
      // -> hopefully, this inconsistency is of small amount
      conv_phi += phin*vdiv;
    }

    // multiply convective term by density
    conv_phi *= densn;
  }

  return;
}


/*-------------------------------------------------------------------------------------- *
 |  standard Galerkin transient, old part of rhs and source term              ehrl 11/13 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSHistAndSource(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  fac,
  const double                  rhsint
  )
{
  double vrhs = fac*rhsint;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    erhs[fvi] += vrhs*funct_(vi);
  }

  return;
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin convective term on right hand side    ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSConv(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  const double                  conv_phi
  )
{
  double vrhs = rhsfac*conv_phi;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    erhs[fvi] -= vrhs*funct_(vi);
  }

  return;
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSDiff(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  Teuchos::RCP<ScaTraEleDiffManager>  diffmanager,
  const LINALG::Matrix<nsd_,1>& gradphi
  )
{
  double vrhs = rhsfac*diffmanager->GetIsotropicDiff(k);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf,gradphi,vi);
    erhs[fvi] -= vrhs*laplawf;
  }

  return;
}


/*--------------------------------------------------------------------------------------------*
 |  transient, convective and diffusive stabilization terms on right hand side     ehrl 11/13 |
 *--------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSTransConvDiffStab(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhstaufac,
  const double                  densnp,
  const double                  scatrares,
  const LINALG::Matrix<nen_,1>& conv,
  const LINALG::Matrix<nen_,1>& sgconv,
  const LINALG::Matrix<nen_,1>& diff
  )
{
  // convective rhs stabilization (in convective form)
  double vrhs = rhstaufac*scatrares*densnp;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    erhs[fvi] -= vrhs*(conv(vi)+sgconv(vi)+scatrapara_->USFEMGLSFac()*1.0/scatraparatimint_->TimeFac()*funct_(vi));
  }

// diffusive rhs stabilization
  if (use2ndderiv_)
  {
    vrhs = rhstaufac*scatrares;
    // diffusive stabilization of convective temporal rhs term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      erhs[fvi] += scatrapara_->USFEMGLSFac()*vrhs*diff(vi);
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | reactive terms (standard Galerkin and stabilization) on rhs   ehrl 11/13  |
 *--------------------------------------------------------------------       */
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSReact(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  const double                  rhstaufac,
  const double                  rea_phi,
  const double                  densnp,
  const double                  scatrares,
  Teuchos::RCP<ScaTraEleReaManager>  reamanager
  )
{

  // standard Galerkin term
  double vrhs = rhsfac*rea_phi;

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    erhs[fvi] -= vrhs*funct_(vi);
  }

  // reactive rhs stabilization
  if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
  {
    vrhs = scatrapara_->USFEMGLSFac()*rhstaufac*densnp*reamanager_->GetReaCoeff(k)*scatrares;
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      erhs[fvi] -= vrhs*funct_(vi);
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | fine-scale subgrid-diffusivity term on right hand side          vg 11/13  |
 *--------------------------------------------------------------------       */
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSFSSGD(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  const double                  sgdiff,
  const LINALG::Matrix<nsd_,1>  fsgradphi
  )
{
  const double vrhs = rhsfac*sgdiff;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf,fsgradphi,vi);
    erhs[fvi] -= (vrhs*laplawf);
  }

  return;
}


/*------------------------------------------------------------------------------*
 | multifractal subgrid-scale modeling on right hand side only rasthofer 11/13  |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcRHSMFS(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  const double                  densnp,
  const LINALG::Matrix<nsd_,1>  convelint,
  const LINALG::Matrix<nsd_,1>  mfsggradphi,
  const LINALG::Matrix<nsd_,1>  mfsgvelint,
  const LINALG::Matrix<nsd_,1>  gradphi,
  const double                  phinp,
  const double                  mfssgphi,
  const double                  mfsvdiv
  )
{
  if (nsd_<3) dserror("Turbulence is 3D!");
  // fixed-point iteration only (i.e. beta=0.0 assumed), cf
  // turbulence part in Evaluate()
  {
   double cross = convelint.Dot(mfsggradphi) + mfsgvelint.Dot(gradphi);
   double reynolds = mfsgvelint.Dot(mfsggradphi);

   // conservative formulation
   double conserv = 0.0;
   if (scatrapara_->MfsConservative() or scatrapara_->IsConservative())
   {
     double convdiv = 0.0;
     GetDivergence(convdiv,econvelnp_);

     conserv = mfssgphi * convdiv + phinp * mfsvdiv + mfssgphi * mfsvdiv;
   }

   const double vrhs = rhsfac*densnp*(cross+reynolds+conserv);

   for (int vi=0; vi<nen_; ++vi)
   {
     const int fvi = vi*numdofpernode_+k;
     //erhs[fvi] -= rhsfac*densnp_[k]*funct_(vi)*(cross+reynolds);
     //erhs[fvi] -= rhsfac*densnp_[k]*funct_(vi)*(cross+reynolds+conserv);
     erhs[fvi] -= vrhs*funct_(vi);
   }
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::nurbs27>;


