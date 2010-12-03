/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl.cpp

\brief Internal implementation of Fluid3 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_impl.H"
#include "fluid3_impl_parameter.H"

#include "../drt_lib/drt_element.H"

#include "../drt_f3/fluid3_stabilization.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"

#include <Epetra_SerialDenseSolver.h>

#ifdef DEBUG
#endif


//----------------------------------------------------------------------*
//
//----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3ImplInterface* DRT::ELEMENTS::Fluid3ImplInterface::Impl(DRT::Element::DiscretizationType distype)
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return Fluid3Impl<DRT::Element::hex8>::Instance();
  }
  case DRT::Element::hex20:
  {
    return Fluid3Impl<DRT::Element::hex20>::Instance();
  }
  case DRT::Element::hex27:
  {
    return Fluid3Impl<DRT::Element::hex27>::Instance();
  }
  case DRT::Element::tet4:
  {
    return Fluid3Impl<DRT::Element::tet4>::Instance();
  }
  case DRT::Element::tet10:
  {
    return Fluid3Impl<DRT::Element::tet10>::Instance();
  }
  case DRT::Element::wedge6:
  {
    return Fluid3Impl<DRT::Element::wedge6>::Instance();
  }
  /* wedge15 cannot be used since no mesh generator exists
  case DRT::Element::wedge15:
  {
    return Fluid3Impl<DRT::Element::wedge15>::Instance();
  }
  */
  case DRT::Element::pyramid5:
  {
    return Fluid3Impl<DRT::Element::pyramid5>::Instance();
  }
  case DRT::Element::quad4:
  {
    return Fluid3Impl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return Fluid3Impl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return Fluid3Impl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return Fluid3Impl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return Fluid3Impl<DRT::Element::tri6>::Instance();
  }
  // no 1D elements
  //nurbs are not available yet
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3Impl<distype> * DRT::ELEMENTS::Fluid3Impl<distype>::Instance()
{
  static Fluid3Impl<distype> * instance;
  if ( instance==NULL )
    instance = new Fluid3Impl<distype>();
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid3Impl<distype>::Fluid3Impl()
  : xyze_(true),
    funct_(true),
    deriv_(true),
    deriv2_(true),
    xjm_(true),
    xji_(true),
    vderxy_(true),
    fsvderxy_(true),
    derxy_(true),
    derxy2_(true),
    bodyforce_(true),
    histmom_(true),
    velino_(true),
    velint_(true),
    fsvelint_(true),
    sgvelint_(true),
    convvelint_(true),
    accint_(true),
    gradp_(true),
    tau_(true),
    viscs2_(true),
    conv_c_(true),
    sgconv_c_(true),  // initialize to zero
    vdiv_(0.0),
    rhsmom_(true),
    conv_old_(true),
    visc_old_(true),
    momres_old_(true),  // initialize to zero
    conres_old_(true),
    xder2_(true),
    vderiv_(true),
    det_(0.0),
    fac_(0.0),
    visc_(0.0),
    sgvisc_(0.0),
    visceff_(0.0),
    fssgvisc_(0.0),
    rhscon_(true),
    densaf_(1.0),         // initialized to 1.0 (filled in Fluid3::GetMaterialParams)
    densam_(1.0),         // initialized to 1.0 (filled in Fluid3::GetMaterialParams)
    densn_(1.0),          // initialized to 1.0 (filled in Fluid3::GetMaterialParams)
    scadtfac_(0.0),       // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    scaconvfacaf_(0.0),   // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    scaconvfacn_(0.0),    // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    thermpressadd_(0.0),  // initialized to 0.0 (filled in Fluid3::GetMaterialParams)
    velintn_(true),
    vderxyn_(true),
    grad_scaaf_(true),
    grad_scan_(true),
    conv_scaaf_(0.0),
    conv_scan_(0.0),
    rotsymmpbc_(NULL),
    //flags
    is_higher_order_ele_(false)
{
  rotsymmpbc_= new FLD::RotationallySymmetricPeriodicBC<distype>();

  // pointer to class Fluid3ImplParameter (access to the general parameter)
  f3Parameter_ = DRT::ELEMENTS::Fluid3ImplParameter::Instance();
}

/*----------------------------------------------------------------------*
 * Action type: Integrate shape function
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3Impl<distype>::IntegrateShapeFunction(
    DRT::ELEMENTS::Fluid3*    ele,
    DRT::Discretization&      discretization,
    vector<int>&              lm            ,
    Epetra_SerialDenseVector& elevec1       )
{
  // --------------------------------------------------
  // construct views
  LINALG::Matrix<numdofpernode_*nen_,    1> vector(elevec1.A(),true);

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

//------------------------------------------------------------------
//                       INTEGRATION LOOP
//------------------------------------------------------------------

  for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    for (int ui=0; ui<nen_; ++ui) // loop rows  (test functions)
    {
      // integrated shape function is written into the pressure dof
      int fuippp=numdofpernode_*ui+nsd_;
      vector(fuippp)+=fac_*funct_(ui);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid3*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elemat1_epetra,
                                                 Epetra_SerialDenseMatrix&  elemat2_epetra,
                                                 Epetra_SerialDenseVector&  elevec1_epetra,
                                                 Epetra_SerialDenseVector&  elevec2_epetra,
                                                 Epetra_SerialDenseVector&  elevec3_epetra )
{
  // rotationally symmetric periodic bc's: do setup for current element
  rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat1(elemat1_epetra,true);
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat2(elemat2_epetra,true);
  LINALG::Matrix<(nsd_+1)*nen_,            1> elevec1(elevec1_epetra,true);
  // elevec2 and elevec3 are currently not in use

  LINALG::Matrix<nsd_,nen_> edeadaf(true);
  BodyForce(ele, f3Parameter_, edeadaf);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  double * saccn = NULL;
  double * sveln = NULL;
  double * svelnp = NULL;

  // if not available, the arrays for the subscale quantities have to be
  // resized and initialised to zero
  if ( f3Parameter_->tds_==INPAR::FLUID::subscales_time_dependent )
  {
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);
    ele->ActivateTDS( intpoints.IP().nquad, nsd_, &saccn, &sveln, &svelnp );
  }

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, scalar,
  // acceleration/scalar time derivative and history
  // velocity/pressure and scalar values are at time n+alpha_F/n+alpha_M
  // for generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration/scalar time derivative values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------

  // fill the local element vector/matrix with the global values
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  LINALG::Matrix<nen_,1> epreaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, &epreaf,"velaf");

  LINALG::Matrix<nen_,1> escaaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &escaaf,"scaaf");

  LINALG::Matrix<nsd_,nen_> emhist(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &emhist, NULL,"hist");

  LINALG::Matrix<nsd_,nen_> eaccam(true);
  LINALG::Matrix<nen_,1> escadtam(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &eaccam, &escadtam,"accam");

  LINALG::Matrix<nsd_,nen_> eveln(true);
  LINALG::Matrix<nen_,1> escaam(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &eveln, &escaam,"scaam");

  if (f3Parameter_->is_genalpha_)
    eveln.Clear();
  else
    eaccam.Clear();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------

  LINALG::Matrix<nsd_, nen_> edispnp(true);
  LINALG::Matrix<nsd_, nen_> egridv(true);

  if (ele->IsAle())
  {
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &egridv, NULL,"gridv");
  }

  // ---------------------------------------------------------------------
  // get additional state vector for AVM3 case: fine-scale velocity
  // values are at time n+alpha_F for generalized-alpha scheme and at
  // time n+1 for all other schemes
  // ---------------------------------------------------------------------

  // fine-scale velocity at time n+alpha_F/n+1
  LINALG::Matrix<nsd_,nen_> fsevelaf(true);
  if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
  {
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &fsevelaf, NULL,"fsvelaf");
  }

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // Call the inner evaluate that does not know about the DRT element or the
  // discretization object.

  int result = Evaluate(
    ele->Id(),
    params,
    edeadaf,
    elemat1,
    elemat2,
    elevec1,
    evelaf,
    epreaf,
    escaaf,
    emhist,
    eaccam,
    escadtam,
    eveln,
    escaam,
    edispnp,
    egridv,
    fsevelaf,
    mat,
    ele->IsAle(),
    ele->Owner()==discretization.Comm().MyPID(),
    ele->CsDeltaSq(),
    saccn,
    sveln,
    svelnp);

  // rotate matrices and vectors if we have a rotationally symmetric problem
  rotsymmpbc_->RotateMatandVecIfNecessary(elemat1,elemat2,elevec1);

  return result;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid3Impl<distype>::Evaluate(
  int                                           eid,
  Teuchos::ParameterList&                       params,
  const LINALG::Matrix<nsd_,nen_> &             edeadaf,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> & elemat1,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> & elemat2,
  LINALG::Matrix<(nsd_+1)*nen_,            1> & elevec1,
  const LINALG::Matrix<nsd_,nen_> & evelaf,
  const LINALG::Matrix<nen_,1>    & epreaf,
  const LINALG::Matrix<nen_,1>    & escaaf,
  const LINALG::Matrix<nsd_,nen_> & emhist,
  const LINALG::Matrix<nsd_,nen_> & eaccam,
  const LINALG::Matrix<nen_,1>    & escadtam,
  const LINALG::Matrix<nsd_,nen_> & eveln,
  const LINALG::Matrix<nen_,1>    & escaam,
  const LINALG::Matrix<nsd_,nen_> & edispnp,
  const LINALG::Matrix<nsd_,nen_> & egridv,
  const LINALG::Matrix<nsd_,nen_> & fsevelaf,
  Teuchos::RCP<MAT::Material>               mat,
  bool                                      isale,
  bool                                      isowned,
  double CsDeltaSq,
  double * saccn,
  double * sveln,
  double * svelnp)
{
  // flag for higher order elements
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (f3Parameter_->is_inconsistent_ == true) is_higher_order_ele_ = false;

  // stationary formulation does not support Ale formulation
  if (isale and f3Parameter_->is_stationary_) dserror("No ALE support within stationary fluid solver.");

  // set thermodynamic pressure at n+1/n+alpha_F and n+alpha_M/n and
  // its time derivative at n+alpha_M/n+1
  const double thermpressaf   = params.get<double>("thermpress at n+alpha_F/n+1");
  const double thermpressam   = params.get<double>("thermpress at n+alpha_M/n");
  const double thermpressdtam = params.get<double>("thermpressderiv at n+alpha_M/n+1");

  // ---------------------------------------------------------------------
  // set parameters for classical turbulence models
  // ---------------------------------------------------------------------
  ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

  double Cs_delta_sq   = 0.0;
  visceff_  = 0.0;

  // remember the layer of averaging for the dynamic Smagorinsky model
  int  nlayer=0;

  GetTurbulenceParams(turbmodelparams,
                      Cs_delta_sq,
                      nlayer,
                      CsDeltaSq);

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------

  Sysmat(eid,
         edeadaf,
         evelaf,
         eveln,
         fsevelaf,
         epreaf,
         eaccam,
         escaaf,
         escaam,
         escadtam,
         emhist,
         edispnp,
         egridv,
         elemat1,
         elemat2,  // -> emesh
         elevec1,
         thermpressaf,
         thermpressam,
         thermpressdtam,
         mat,
         Cs_delta_sq,
         isale,
         saccn,
         sveln,
         svelnp);

  // ---------------------------------------------------------------------
  // output values of Cs, visceff and Cs_delta_sq
  // ---------------------------------------------------------------------
  //
  // do the fastest test first
  if (isowned)
  {
    if (f3Parameter_->turb_mod_action_ == INPAR::FLUID::dynamic_smagorinsky
        or
        f3Parameter_->turb_mod_action_ == INPAR::FLUID::smagorinsky_with_van_Driest_damping
      )
    {
      if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
      {
        if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
            ==
            "channel_flow_of_height_2")
        {
          // Cs was changed in Sysmat (Cs->sqrt(Cs_delta_sq)/pow((vol),(1.0/3.0)))
          // to compare it with the standard Smagorinsky Cs
          (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_sum")))         [nlayer]+=f3Parameter_->Cs_;
          (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_delta_sq_sum")))[nlayer]+=Cs_delta_sq;
          (*(turbmodelparams.get<RCP<vector<double> > >("local_visceff_sum")))    [nlayer]+=visceff_;
        }
      }
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::Sysmat(
  int                                           eid,
  const LINALG::Matrix<nsd_,nen_>&              edeadaf,
  const LINALG::Matrix<nsd_,nen_>&              evelaf,
  const LINALG::Matrix<nsd_,nen_>&              eveln,
  const LINALG::Matrix<nsd_,nen_>&              fsevelaf,
  const LINALG::Matrix<nen_,1>&                 epreaf,
  const LINALG::Matrix<nsd_,nen_>&              eaccam,
  const LINALG::Matrix<nen_,1>&                 escaaf,
  const LINALG::Matrix<nen_,1>&                 escaam,
  const LINALG::Matrix<nen_,1>&                 escadtam,
  const LINALG::Matrix<nsd_,nen_>&              emhist,
  const LINALG::Matrix<nsd_,nen_>&              edispnp,
  const LINALG::Matrix<nsd_,nen_>&              egridv,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
  LINALG::Matrix<(nsd_+1)*nen_,1>&              eforce,
  const double                                  thermpressaf,
  const double                                  thermpressam,
  const double                                  thermpressdtam,
  Teuchos::RCP<const MAT::Material>             material,
  double&                                       Cs_delta_sq,
  bool                                          isale,
  double * saccn,
  double * sveln,
  double * svelnp
  )
{
  LINALG::Matrix<nen_*nsd_,nen_*nsd_>     estif_u(true);
  LINALG::Matrix<nen_*nsd_,nen_>          estif_p_v(true);
  LINALG::Matrix<nen_, nen_*nsd_>         estif_q_u(true);
  LINALG::Matrix<nen_,nen_>               ppmat(true);

  LINALG::Matrix<nen_,1>                  preforce(true);
  LINALG::Matrix<nsd_,nen_>               velforce(true);

  //! linearisation of residual of momentum equation wrt velocities (precomputed)
  LINALG::Matrix<nsd_*nsd_,nen_>          lin_resM_Du(true);
  //! linearisation of residual of momentum equation wrt velocities (precomputed)
  LINALG::Matrix<nsd_,1>                  resM_Du(true);

  // add displacement when fluid nodes move in the ALE case
  if (isale) xyze_ += edispnp;

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(eid);

  // element aera or volume
  const double vol = fac_;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (f3Parameter_->vstab_ == INPAR::FLUID::viscous_stab_usfem or
      f3Parameter_->vstab_ == INPAR::FLUID::viscous_stab_usfem_only_rhs)   vstabfac =  1.0;
  else if(f3Parameter_->vstab_ == INPAR::FLUID::viscous_stab_gls or
      f3Parameter_->vstab_ == INPAR::FLUID::viscous_stab_gls_only_rhs) vstabfac = -1.0;

  //----------------------------------------------------------------------
  // get material parameters at element center
  //----------------------------------------------------------------------
  if (not f3Parameter_->mat_gp_ or not f3Parameter_->tau_gp_)
    GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);

  if (not f3Parameter_->tau_gp_)
  {
    // ---------------------------------------------------------------------
    // calculate all-scale or fine-scale subgrid viscosity at element center
    // ---------------------------------------------------------------------
    visceff_ = visc_;
    if (f3Parameter_->turb_mod_action_ != INPAR::FLUID::no_model)
    {
      CalcSubgrVisc(evelaf,vol,f3Parameter_->Cs_,Cs_delta_sq,f3Parameter_->l_tau_);

      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
      CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,f3Parameter_->Cs_);

    // get velocity at element center
    velint_.Multiply(evelaf,funct_);

    // ---------------------------------------------------------------------
    // calculate stabilization parameter at element center
    // ---------------------------------------------------------------------
    CalcStabParameter(f3Parameter_->timefac_,vol);
  }

  // Gaussian integration points
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  //------------------------------------------------
  //------------------------------------------------
  //    INTEGRATION LOOP
  //------------------------------------------------
  //------------------------------------------------

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {

  //----------------------------------------------------------------------
  // evaluate shape functions and derivatives at integration point
  //----------------------------------------------------------------------

    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,eid);

  //----------------------------------------------------------------------
  // get material parameters (evaluation at integration point)
  //----------------------------------------------------------------------

    if (f3Parameter_->mat_gp_)
      GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);

  // ---------------------------------------------------------------------
  // CALCULATE all-scale / fine-scale subgrid viscosity and
  // stabilization parameter at integration point
  // ---------------------------------------------------------------------

    if (f3Parameter_->tau_gp_)
    {
      visceff_ = visc_;
      if (f3Parameter_->turb_mod_action_ != INPAR::FLUID::no_model)
      {
        CalcSubgrVisc(evelaf,vol,f3Parameter_->Cs_,Cs_delta_sq,f3Parameter_->l_tau_);

        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
        CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,f3Parameter_->Cs_);

      // Stabilization parameter
      CalcStabParameter(f3Parameter_->timefac_,vol);
    }

  //--------------------------------------------------------------------
  // COMPUTE numerical representation of some single operators
  //--------------------------------------------------------------------

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf,funct_);

    // get momentum history data at integration point
    histmom_.Multiply(emhist,funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelaf,derxy_);

    // get fine-scale velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv) fsvderxy_.MultiplyNT(fsevelaf,derxy_);
    else                           fsvderxy_.Clear();

    // get convective velocity at integration point
    // We handle the ale case very implicitly here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_.Update(velint_);
    if (isale) convvelint_.Multiply(-1.0, egridv, funct_, 1.0);

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    gradp_.Multiply(derxy_,epreaf);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = funct_.Dot(epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(edeadaf,funct_);

    // get second derivative of the viscous term:
    // div(epsilon(u))
    if (is_higher_order_ele_)
    {
      CalcDivEps(evelaf);
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }

    // convective term from previous iteration
    conv_old_.Multiply(vderxy_,convvelint_);

    // compute convective operator
    conv_c_.MultiplyTN(derxy_,convvelint_);


    // velocity divergence from previous iteration
    vdiv_ = 0.0;
    for (int idim = 0; idim <nsd_; ++idim)
    {
      vdiv_ += vderxy_(idim, idim);
    }

  //--------------------------------------------------------------------
  // factors for stabilization, time integration
  // and fine-scale subgrid-viscosity
  //--------------------------------------------------------------------

    const double tau_M       = tau_(0)*fac_;
    const double tau_Mp      = tau_(1)*fac_;
    const double tau_C       = tau_(2)*fac_;

    const double timefacfac  = f3Parameter_->timefac_ * fac_;
    const double timetauM    = f3Parameter_->timefac_ * tau_M;

    double rhsfac            = fac_;

  //--------------------------------------------------------------------
  // The following computations are performed depending on
  // time-integration, that is, whether it is generalized-alpha or not,
  // since several terms differ with respect to the scheme:
  //
  // 1) calculation of rhs for momentum equation and momentum residual
  // -> different for generalized-alpha and other schemes
  //
  // 2) calculation of additional subgrid-scale velocity when cross-
  //    and Reynolds-stress are included:
  // - Cross- and Reynolds-stress are always included simultaneously.
  // - They are included in a complete form on left- and right-hand side.
  // - For this purpose, a subgrid-scale convective term is computed.
  // - Within a Newton linearization, the present formulation is not
  //   consistent for the reactive terms.
  // - To turn them off, both flags must be "no".
  //
  // 3) calculation of convective scalar term, rhs for continuity
  //    equation and residual of continuity equations
  // -> only required for low-Mach-number flow
  // -> for incompressible flow, residual is velocity divergence only
  //--------------------------------------------------------------------

  /*-------------------------------------------------------------------*
  *                                                                   *
  *                  get residual of momentum equation                *
  *                                                                   *
  *-------------------------------------------------------------------*/

    GetResidualMomentumEq(eaccam,
                          f3Parameter_->timefac_,
                          rhsfac);

  /*-------------------------------------------------------------------*
  *                                                                   *
  *                  update of SUBSCALE VELOCITY                      *
  *                                                                   *
  *-------------------------------------------------------------------*/
    double fac1   =0.0;
    double fac2   =0.0;
    double fac3   =0.0;
    double facMtau=0.0;

    UpdateSubscaleVelocity(fac1,
                           fac2,
                           fac3,
                           facMtau,
                           iquad,
                           saccn,
                           sveln,
                           svelnp);


  /*-------------------------------------------------------------------*
  *                                                                   *
  *                 get residual of continuity equation               *
  *                                                                   *
  *-------------------------------------------------------------------*/


    GetResidualContinuityEq(eveln,
                            escaaf,
                            escaam,
                            escadtam,
                            f3Parameter_->timefac_);



  //------------------------------------------------------------------------
  // perform integration for element matrix and right hand side
  //------------------------------------------------------------------------

    lin_resM_Du.Clear();
    resM_Du.Clear();

  //----------------------------------------------------------------------
  //   PROVIDE LINEARISATION OF GALERKIN MOMENTUM RESIDUAL WRT VELOCITIES
  //-------------------------------------------------------------------------

    LinGalMomResU(lin_resM_Du,
                  timefacfac);

  //----------------------------------------------------------------------
  //   RESCALE GALERKIN RESIDUAL
  //-------------------------------------------------------------------------

    if(f3Parameter_->tds_      ==INPAR::FLUID::subscales_time_dependent
       &&
       f3Parameter_->transient_==INPAR::FLUID::inertia_stab_keep)
    {
      LinGalMomResU_subscales(estif_p_v,
                              lin_resM_Du,
                              resM_Du,
                              timefacfac,
                              facMtau);
    }

  //----------------------------------------------------------------------
  //   INERTIA & CONVECTION (CONVECTIVE AND REACTIVE PART)
  //   including rhs contribution
  //----------------------------------------------------------------------

    InertiaAndConvectionGalPart(estif_u,
                                velforce,
                                lin_resM_Du,
                                resM_Du,
                                rhsfac);

  //----------------------------------------------------------------------
  //   VISCOUS GALERKIN PART
  //   including viscous stress computation
  //
  //   excluding viscous Galerkin Part of the LOMA formulation
  //-----------------------------------------------------------------------

    // viscous stresses
    LINALG::Matrix<nsd_,nsd_> viscstress(true);

    ViscousGalPart(estif_u,
                   viscstress,
                   timefacfac,
                   rhsfac);

  //----------------------------------------------------------------------
  //   STABILISATION OF CONTINUITY PART AND RHS OF VISCOUS TERM
  //   adding LOMA contribution to the viscous stress tensor
  //----------------------------------------------------------------------

    ContStab_and_ViscousTermRhs(estif_u,
                                velforce,
                                viscstress,
                                f3Parameter_->timefac_,
                                timefacfac,
                                tau_C,
                                rhsfac);


  //----------------------------------------------------------------------
  //   COMPUTATION OF PRESSURE TERM
  //   including rhs contribution
  //----------------------------------------------------------------------

    PressureGalPart(estif_p_v,
                    velforce,
                    timefacfac,
                    rhsfac,
                    press);

  //----------------------------------------------------------------------
  //   COMPUTATION OF CONTINUITY TERM
  //   including rhs contribution
  //----------------------------------------------------------------------

    ContinuityGalPart(estif_q_u,
                      preforce,
                      timefacfac,
                      rhsfac);

  //----------------------------------------------------------------------
  // COMPUTATION OF BODY FORCE TERM on right-hand side
  //----------------------------------------------------------------------

      BodyForceRhsTerm(velforce);

  //----------------------------------------------------------------------
  //   CONSERVATIVE FORMULATION (GENERAL & LOMA & VARYING DENSITY)
  //   including rhs contribution
  //----------------------------------------------------------------------

    if (f3Parameter_->is_conservative_)
    {
      ConservativeFormulation(estif_u,
                              velforce,
                              timefacfac,
                              rhsfac);
    }

  //----------------------------------------------------------------------
  //   COMPUTATION OF ADITIONAL TERMS (continuity equation) for LOMA flow
  //   including rhs contribution
  //----------------------------------------------------------------------

    if (f3Parameter_->physicaltype_ == INPAR::FLUID::loma)
    //if(loma_)
    {
      LomaGalPart(estif_q_u,
                  preforce,
                  timefacfac,
                  rhsfac);
    }

  //----------------------------------------------------------------------
  //  PROVIDE LINEARISATION OF MOMENTUM RESIDUAL WRT VELOCITIES FOR ALL
  //       STABILISATION PARTS (EXTEND GALERKIN RES IF NECESSARY)
  //----------------------------------------------------------------------

    StabLinGalMomResU(lin_resM_Du,
                      timefacfac);

  //----------------------------------------------------------------------
  //   PRESSURE STABILISATION PART
  //   including rhs contribution
  //----------------------------------------------------------------------

    if (f3Parameter_->pspg_ == INPAR::FLUID::pstab_use_pspg)
    {
      PSPG(estif_q_u,
           ppmat,
           preforce,
           lin_resM_Du,
           fac3,
           timefacfac,
           tau_Mp,
           rhsfac);
    }

  //----------------------------------------------------------------------
  //    SUPG STABILISATION PART
  //    including rhs socntribution
  //-----------------------------------------------------------------------

    if(f3Parameter_->supg_ == INPAR::FLUID::convective_stab_supg)
    {
      SUPG(estif_u,
           estif_p_v,
           velforce,
           lin_resM_Du,
           fac3,
           timefacfac,
           timetauM);
    }

  //----------------------------------------------------------------------
  //                       STABILISATION, VISCOUS PART
  //----------------------------------------------------------------------


    if (is_higher_order_ele_ && (f3Parameter_->vstab_ != INPAR::FLUID::viscous_stab_none))
    {
      ViscStab(estif_u,
               estif_p_v,
               velforce,
               lin_resM_Du,
               f3Parameter_->timefac_,
               tau_Mp,
               vstabfac,
               fac3);
    }


  //----------------------------------------------------------------------
  //    CROSS-STRESS STABILISATION PART
  //    including rhs contribution
  //----------------------------------------------------------------------

    if(f3Parameter_->cross_ != INPAR::FLUID::cross_stress_stab_none)
    {
      CrossStressStab(estif_u,
                      estif_p_v,
                      velforce,
                      lin_resM_Du,
                      f3Parameter_->timefac_,
                      timefacfac,
                      tau_Mp,
                      fac3);
    }

  //----------------------------------------------------------------------
  //    REYNOLDS-STRESS STABILISATION PART
  //    including rhs contribution
  //----------------------------------------------------------------------

    if(f3Parameter_->reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
    {
      ReynoldsStressStab(estif_u,
                         estif_p_v,
                         lin_resM_Du,
                         timefacfac,
                         fac3);
    }

  //----------------------------------------------------------------------
  //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)
  //----------------------------------------------------------------------

    if(f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
    {
      const double fssgviscfac = fssgvisc_*rhsfac;

      FineScaleSubGridViscosityTerm(velforce,
                                    fssgviscfac);
    }

  //--------------------------------------------------------------------
  // LINEARIZATION WITH RESPECT TO MESH MOTION
  //---------------------------------------------------------------------
    if (emesh.IsInitialized())
    {
      if (nsd_ == 3)
        LinMeshMotion_3D(emesh,
                        evelaf,
                        press,
                        f3Parameter_->timefac_,
                        timefacfac);
      else if(nsd_ == 2)
        LinMeshMotion_2D(emesh,
                         evelaf,
                         press,
                         f3Parameter_->timefac_,
                         timefacfac);
      else
        dserror("Linearization of the mesh motion is not available in 1D");
    }
  } // end loop gausspoints


  //--------------------------------------------------------------------
  //   FILL ELEMENT MATRIX
  //---------------------------------------------------------------------

  // add pressure part of residual to force vector
  for (int vi=0; vi<nen_; ++vi)
  {
    eforce(numdofpernode_*vi+nsd_)+=preforce(vi);
  }

  // add velocity part of residual to force vector
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      eforce(numdofpernode_*vi+idim)+=velforce(idim,vi);
    }
  }

  // add pressure/pressure part of matrix to the element stiffness
  for (int ui=0; ui<nen_; ++ui)
  {
    const int fuippp = numdofpernode_*ui+nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int numdof_vi_p_nsd = numdofpernode_*vi+nsd_;

      estif(numdof_vi_p_nsd,fuippp)+=ppmat(vi,ui);
    }
  }

  // sort estif_ into the global matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int numdof_vi = numdofpernode_*vi;
        const int nsd_vi = nsd_*vi;

        for (int idim=0; idim <nsd_; ++idim)
        {
          estif(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // sort estif_p_v_ into the global matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui_nsd = numdofpernode_*ui + nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int nsd_vi = nsd_*vi;
      const int numdof_vi = numdofpernode_*vi;

      for (int idim=0; idim <nsd_; ++idim)
      {
        estif(numdof_vi+idim, numdof_ui_nsd) += estif_p_v(nsd_vi+idim, ui);
      }
    }
  }

  // sort estif_q_u_ into the global matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
        estif(numdofpernode_*vi+nsd_, numdof_ui_jdim) += estif_q_u(vi, nsd_ui_jdim);
    }
  }

  return;
}


/*!
      \brief Do a finite difference check for a given element id ---
      this function is for debugging purposes only

      \param ele              (i) the element those matrix is calculated
                                  (pass-through)
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1 (pass-through)
      \param eveln            (i) nodal velocities at n (pass-through)
      \param fsevelaf         (i) fine-scale nodal velocities at n+alpha_F/n+1
                                  (pass-through)
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1 (pass-through)
      \param eaccam           (i) nodal accelerations at n+alpha_M (pass-through)
      \param escaaf           (i) nodal scalar at n+alpha_F/n+1 (pass-through)
      \param escaam           (i) nodal scalar at n+alpha_M/n (pass-through)
      \param escadtam         (i) nodal scalar derivatives at n+alpha_M/n+1
                                  (pass-through)
      \param emhist           (i) time rhs for momentum equation (pass-through)
      \param edispnp          (i) nodal displacements (on moving mesh)
                                  (pass-through)
      \param egridv           (i) grid velocity (on moving mesh) (pass-through)
      \param estif            (i) element matrix to calculate (pass-through)
      \param emesh            (i) linearization wrt mesh motion (pass-through)
      \param eforce           (i) element rhs to calculate (pass-through)
      \param material         (i) fluid material (pass-through)
      \param time             (i) current simulation time (pass-through)
      \param timefac          (i) time discretization factor (pass-through)
      \param newton           (i) boolean flag for linearisation (pass-through)
      \param loma             (i) boolean flag for potential low-Mach-number solver
                                  (pass-through)
      \param conservative     (i) boolean flag for conservative form (pass-through)
      \param is_genalpha      (i) boolean flag for generalized-alpha time
                                  integration (pass-through)
      \param higher_order_ele (i) keep or drop second derivatives (pass-through)
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity
                                  (pass-through)
      \param pspg             (i) boolean flag for stabilisation (pass-through)
      \param supg             (i) boolean flag for stabilisation (pass-through)
      \param vstab            (i) boolean flag for stabilisation (pass-through)
      \param cstab            (i) boolean flag for stabilisation (pass-through)
      \param cross            (i) boolean flag for stabilisation (pass-through)
      \param reynolds         (i) boolean flag for stabilisation (pass-through)
      \param turb_mod_action  (i) selecting turbulence model (none, Smagorisky,
                                  dynamic Smagorinsky, Smagorinsky with van Driest
                                  damping for channel flows) (pass-through)
      \param Cs               (i) Smagorinsky model parameter (pass-through)
      \param Cs_delta_sq      (i) Model parameter computed by dynamic Smagorinsky
                                  approach (Cs*h*h) (pass-through)
      \param l_tau            (i) viscous length scale, required for van driest
                                  damping function and defined on input (pass-through)
*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::FDcheck(
  int                                                   eid,
  const LINALG::Matrix<nsd_,nen_>&                      evelaf,
  const LINALG::Matrix<nsd_,nen_>&                      eveln,
  const LINALG::Matrix<nsd_,nen_>&                      fsevelaf,
  const LINALG::Matrix<nen_,1>&                         epreaf,
  const LINALG::Matrix<nsd_,nen_>&                      eaccam,
  const LINALG::Matrix<nen_,1>&                         escaaf,
  const LINALG::Matrix<nen_,1>&                         escaam,
  const LINALG::Matrix<nen_,1>&                         escadtam,
  const LINALG::Matrix<nsd_,nen_>&                      emhist,
  const LINALG::Matrix<nsd_,nen_>&                      edispnp,
  const LINALG::Matrix<nsd_,nen_>&                      egridv,
  const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    estif,
  const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    emesh,
  const LINALG::Matrix<(nsd_+1)*nen_,    1>&            eforce,
  const double                                          thermpressaf,
  const double                                          thermpressam,
  const double                                          thermpressdtam,
  const Teuchos::RCP<const MAT::Material>               material,
  const double                                          timefac,
  const double&                                         Cs,
  const double&                                         Cs_delta_sq,
  const double&                                         l_tau)
{
  // magnitude of dof perturbation
  const double epsilon=1e-14;

  // make a copy of all input parameters potentially modified by Sysmat
  // call --- they are not intended to be modified
  double copy_Cs         =Cs;
  double copy_Cs_delta_sq=Cs_delta_sq;
  double copy_l_tau      =l_tau;

  Teuchos::RCP<const MAT::Material> copy_material=material;

  // allocate arrays to compute element matrices and vectors at perturbed
  // positions
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> checkmat1(true);
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> checkmat2(true);
  LINALG::Matrix<(nsd_+1)*nen_,            1> checkvec1(true);

  // alloc the vectors that will contain the perturbed velocities or
  // pressures
  LINALG::Matrix<nsd_,nen_>                   checkevelaf(true);
  LINALG::Matrix<nsd_,nen_>                   checkeaccam(true);
  LINALG::Matrix<nen_,1>                      checkepreaf(true);

  // echo to screen
  printf("+-------------------------------------------+\n");
  printf("| FINITE DIFFERENCE CHECK FOR ELEMENT %5d |\n",eid);
  printf("+-------------------------------------------+\n");
  printf("\n");
  // loop columns of matrix by looping nodes and then dof per nodes

  // loop nodes
  for(int nn=0;nn<nen_;++nn)
  {
    printf("-------------------------------------\n");
    printf("-------------------------------------\n");
    printf("NODE of element local id %d\n",nn);
    // loop dofs
    for(int rr=0;rr<(nsd_+1);++rr)
    {
      // number of the matrix column to check
      int dof=nn*(nsd_+1)+rr;

      // clear element matrices and vectors to assemble
      checkmat1.Clear();
      checkmat2.Clear();
      checkvec1.Clear();

      // copy velocities and pressures to perturbed arrays
      for(int mm=0;mm<nen_;++mm)
      {
        for(int dim=0;dim<nsd_;++dim)
        {
          checkevelaf(dim,mm)=evelaf(dim,mm);

          checkeaccam(dim,mm)=eaccam(dim,mm);
        }

        checkepreaf(  mm)=epreaf(  mm);
      }

      // perturb the respective elemental quantities
      if(rr==nsd_)
      {
        printf("pressure dof (%d) %f\n",nn,epsilon);

        if (f3Parameter_->is_genalpha_)
        {
          checkepreaf(nn)+=f3Parameter_->alphaF_*epsilon;
        }
        else
        {
          checkepreaf(nn)+=epsilon;
        }
      }
      else
      {
        printf("velocity dof %d (%d)\n",rr,nn);

        if (f3Parameter_->is_genalpha_)
        {
          checkevelaf(rr,nn)+=f3Parameter_->alphaF_*epsilon;
          checkeaccam(rr,nn)+=f3Parameter_->alphaM_/(f3Parameter_->gamma_*f3Parameter_->dt_)*epsilon;
        }
        else
        {
          checkevelaf(rr,nn)+=epsilon;
        }
      }

      // calculate the right hand side for the perturbed vector
      Sysmat2D3D(checkevelaf,
                 eveln,
                 fsevelaf,
                 checkepreaf,
                 checkeaccam,
                 escaaf,
                 escaam,
                 escadtam,
                 emhist,
                 edispnp,
                 egridv,
                 checkmat1,
                 checkmat2,
                 checkvec1,
                 thermpressaf,
                 thermpressam,
                 thermpressdtam,
                 copy_material,
                 timefac,
                 copy_Cs,
                 copy_Cs_delta_sq,
                 copy_l_tau);

      // compare the difference between linaer approximation and
      // (nonlinear) right hand side evaluation

      // note that it makes more sense to compare these quantities
      // than to compare the matrix entry to the difference of the
      // the right hand sides --- the latter causes numerical problems
      // do to deletion

      for(int mm=0;mm<(nsd_+1)*nen_;++mm)
      {
        double val;
        double lin;
        double nonlin;

        // For af-generalized-alpha scheme, the residual vector for the
        // solution rhs is scaled on the time-integration level...
        if (f3Parameter_->is_genalpha_)
        {
          val   =-(eforce(mm)   /(epsilon))*(f3Parameter_->gamma_*f3Parameter_->dt_)/(f3Parameter_->alphaM_);
          lin   =-(eforce(mm)   /(epsilon))*(f3Parameter_->gamma_*f3Parameter_->dt_)/(f3Parameter_->alphaM_)+estif(mm,dof);
          nonlin=-(checkvec1(mm)/(epsilon))*(f3Parameter_->gamma_*f3Parameter_->dt_)/(f3Parameter_->alphaM_);
        }
        else
        {
          val   =-eforce(mm)/epsilon;
          lin   =-eforce(mm)/epsilon+estif(mm,dof);
          nonlin=-checkvec1(mm)/epsilon;
        }

        double norm=abs(lin);
        if(norm<1e-12)
        {
          norm=1e-12;
        }

        // output to screen
        printf("relerr         %+12.5e ",(lin-nonlin)/norm);
        printf("abserr         %+12.5e ",lin-nonlin);
        printf("orig. value    %+12.5e ",val);
        printf("lin. approx.   %+12.5e ",lin);
        printf("nonlin. funct. %+12.5e ",nonlin);
        printf("matrix entry   %+12.5e ",estif(mm,dof));
        printf("\n");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadaf only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::BodyForce(
           DRT::ELEMENTS::Fluid3*               ele,
           DRT::ELEMENTS::Fluid3ImplParameter*  f3Parameter,
           LINALG::Matrix<nsd_,nen_>&           edeadaf)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  if (nsd_==3)
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
  else if (nsd_==2)
    DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
  else
    dserror("Body force for a 1D problem is not yet implemented");

  if (myneumcond.size()>1)
    dserror("more than one VolumeNeumann cond on one node");

  if (myneumcond.size()==1)
  {
    // find out whether we will use a time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if (f3Parameter->time_ >= 0.0)
      {
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(f3Parameter->time_);
      }
      else
      {
        // do not compute an "alternative" curvefac here since a negative time value
        // indicates an error.
        dserror("Negative time value in body force calculation: time = %f", f3Parameter->time_);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );
    const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the edeadaf array
    for (int isd=0;isd<nsd_;isd++)
    {
      // get factor given by spatial function
      if (functions) functnum = (*functions)[isd];
      else functnum = -1;

      double num = (*onoff)[isd]*(*val)[isd]*curvefac;

      for ( int jnode=0; jnode<nen_; ++jnode )
      {
        if (functnum>0)
        {
          // evaluate function at the position of the current node
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,
                                                                             (ele->Nodes()[jnode])->X(),
                                                                             f3Parameter->time_,
                                                                             NULL);
        }
        else functionfac = 1.0;

        edeadaf(isd,jnode) = num*functionfac;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at element center  vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::EvalShapeFuncAndDerivsAtEleCenter(
  const int  eleid
)
{
  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_stab(DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);

  // coordinates of the current integration point
  const double* gpcoord = (intpoints_stab.IP().qxg)[0];
  for (int idim=0;idim<nsd_;idim++)
  {
	  xsi_(idim) = gpcoord[idim];
  }
  const double wquad = intpoints_stab.IP().qwgt[0];

  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);


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

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  det_ = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = wquad*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point   vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,		// integration points
    const int                              iquad,				// actual integration point
    const int                              eleid				// element ID
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {
	  xsi_(idim) = gpcoord[idim];
  }

  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // get Jacobian matrix and determinant
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
  det_ = xji_.Invert(xjm_);

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = intpoints.IP().qwgt[iquad]*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::shape_function_deriv2<distype>(xsi_, deriv2_);
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}


/*----------------------------------------------------------------------*
 |  compute material parameters                                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::GetMaterialParams(
  Teuchos::RCP<const MAT::Material>  material,
  const LINALG::Matrix<nsd_,nen_>&       evelaf,
  const LINALG::Matrix<nen_,1>&       escaaf,
  const LINALG::Matrix<nen_,1>&       escaam,
  const double                        thermpressaf,
  const double                        thermpressam,
  const double                        thermpressdtam
)
{
// initially set density values and values with respect to continuity rhs
// -> will only be changed for low-Mach-number flow "material" below
densam_        = 1.0;
densaf_        = 1.0;
densn_         = 1.0;
scadtfac_      = 0.0;
scaconvfacaf_  = 0.0;
scaconvfacn_   = 0.0;
thermpressadd_ = 0.0;

if (material->MaterialType() == INPAR::MAT::m_fluid)
{
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

  // get constant viscosity
  visc_ = actmat->Viscosity();

  // Varying Density
  if (f3Parameter_->physicaltype_ == INPAR::FLUID::varying_density)
  {
	  densaf_ = funct_.Dot(escaaf);
	  densam_ = densaf_;

	  densn_ = funct_.Dot(escaam);
  }
  // Boussinesq approximation: Calculation of delta rho
  else if (f3Parameter_->physicaltype_ == INPAR::FLUID::boussinesq)
  {
    const double density_0 = actmat->Density();

    if(escaaf(0) < EPS12)
      dserror("Boussinesq approximation: density in escaaf is zero");

	  deltadens_ =  (funct_.Dot(escaaf)- density_0)/ density_0;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_carreauyasuda)
{
  const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(material.get());

  double nu_0   = actmat->Nu0();    // parameter for zero-shear viscosity
  double nu_inf = actmat->NuInf();  // parameter for infinite-shear viscosity
  double lambda = actmat->Lambda();  // parameter for characteristic time
  double a      = actmat->AParam(); // constant parameter
  double b      = actmat->BParam(); // constant parameter

  // compute rate of strain at n+alpha_F or n+1
  double rateofstrain   = -1.0e30;
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
  // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
  const double tmp = pow(lambda*rateofstrain,b);
  visc_ = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
}
else if (material->MaterialType() == INPAR::MAT::m_modpowerlaw)
{
  const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(material.get());

  // get material parameters
  double m     = actmat->MCons();     // consistency constant
  double delta = actmat->Delta();      // safety factor
  double a     = actmat->AExp();      // exponent

  // compute rate of strain at n+alpha_F or n+1
  double rateofstrain   = -1.0e30;
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  // compute viscosity according to a modified power law model for shear-thinning fluids
  // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
  visc_ = m * pow((delta + rateofstrain), (-1)*a);
}
else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
{
  const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

  // compute mixture fraction at n+alpha_F or n+1
  const double mixfracaf = funct_.Dot(escaaf);

  // compute dynamic viscosity at n+alpha_F or n+1 based on mixture fraction
  visc_ = actmat->ComputeViscosity(mixfracaf);

  // compute density at n+alpha_F or n+1 based on mixture fraction
  densaf_ = actmat->ComputeDensity(mixfracaf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = actmat->EosFacA()*densaf_;

  if (f3Parameter_->is_genalpha_)
  {
    // compute density at n+alpha_M based on mixture fraction
    const double mixfracam = funct_.Dot(escaam);
    densam_ = actmat->ComputeDensity(mixfracam);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = actmat->EosFacA()*densam_;
  }
  else
  {
    // compute density at n based on mixture fraction
    const double mixfracn = funct_.Dot(escaam);
    densn_ = actmat->ComputeDensity(mixfracn);

    // factor for convective scalar term at n
    scaconvfacn_ = actmat->EosFacA()*densn_;

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // factor for scalar time derivative
    scadtfac_ = scaconvfacaf_;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_sutherland)
{
  const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

  // compute temperature at n+alpha_F or n+1
  const double tempaf = funct_.Dot(escaaf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on temperature
  // and thermodynamic pressure
  densaf_ = actmat->ComputeDensity(tempaf,thermpressaf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = 1.0/tempaf;

  if (f3Parameter_->is_genalpha_)
  {
    // compute temperature at n+alpha_M
    const double tempam = funct_.Dot(escaam);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = 1.0/tempam;

    // compute density at n+alpha_M based on temperature
    densam_ = actmat->ComputeDensity(tempam,thermpressam);

    // addition due to thermodynamic pressure at n+alpha_M
    thermpressadd_ = -thermpressdtam/thermpressam;
  }
  else
  {
    // compute temperature at n
    const double tempn = funct_.Dot(escaam);

    // compute density at n based on temperature at n and
    // (approximately) thermodynamic pressure at n+1
    densn_ = actmat->ComputeDensity(tempn,thermpressaf);

    // factor for convective scalar term at n
    scaconvfacn_ = 1.0/tempn;

    // factor for scalar time derivative
    scadtfac_ = scaconvfacaf_;

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // addition due to thermodynamic pressure
    thermpressadd_ = -(thermpressaf-thermpressam)/thermpressaf;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
{
  const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

  // get progress variable at n+alpha_F or n+1
  const double provaraf = funct_.Dot(escaaf);

  // compute temperature based on progress variable at n+alpha_F or n+1
  const double tempaf = actmat->ComputeTemperature(provaraf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on progress variable
  densaf_ = actmat->ComputeDensity(provaraf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = actmat->ComputeFactor(provaraf);

  if (f3Parameter_->is_genalpha_)
  {
    // compute density at n+alpha_M based on progress variable
    const double provaram = funct_.Dot(escaam);
    densam_ = actmat->ComputeDensity(provaram);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = actmat->ComputeFactor(provaram);
  }
  else
  {
    // compute density at n based on progress variable
    const double provarn = funct_.Dot(escaam);
    densn_ = actmat->ComputeDensity(provarn);

    // factor for convective scalar term at n
    scaconvfacn_ = actmat->ComputeFactor(provarn);

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // factor for scalar time derivative
    scadtfac_ = scaconvfacaf_;
  }
}
else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
{
  const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

  // get progress variable at n+alpha_F or n+1
  const double provaraf = funct_.Dot(escaaf);

  // compute temperature based on progress variable at n+alpha_F or n+1
  const double tempaf = actmat->ComputeTemperature(provaraf);

  // compute viscosity according to Sutherland law
  visc_ = actmat->ComputeViscosity(tempaf);

  // compute density at n+alpha_F or n+1 based on progress variable
  densaf_ = actmat->ComputeDensity(provaraf);

  // factor for convective scalar term at n+alpha_F or n+1
  scaconvfacaf_ = actmat->ComputeFactor(provaraf);

  if (f3Parameter_->is_genalpha_)
  {
    // compute density at n+alpha_M based on progress variable
    const double provaram = funct_.Dot(escaam);
    densam_ = actmat->ComputeDensity(provaram);

    // factor for scalar time derivative at n+alpha_M
    scadtfac_ = actmat->ComputeFactor(provaram);
  }
  else
  {
    // compute density at n based on progress variable
    const double provarn = funct_.Dot(escaam);
    densn_ = actmat->ComputeDensity(provarn);

    // factor for convective scalar term at n
    scaconvfacn_ = actmat->ComputeFactor(provarn);

    // set density at n+1 at location n+alpha_M as well
    densam_ = densaf_;

    // factor for scalar time derivative
    scadtfac_ = scaconvfacaf_;
  }
}
else dserror("Material type is not supported");

// check whether there is zero or negative (physical) viscosity
if (visc_ < EPS15) dserror("zero or negative (physical) diffusivity");

return;
} // Fluid3Impl::GetMaterialParams


/*----------------------------------------------------------------------*
 |  compute turbulence parameters                            ehrl 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::GetTurbulenceParams(
                               ParameterList&             turbmodelparams,
                               double&                    Cs_delta_sq,
                               int&                       nlayer,
                               double CsDeltaSq)
{
  if(f3Parameter_->turb_mod_action_ != INPAR::FLUID::no_model and nsd_ == 2)
    dserror("turbulence and 2D flow does not make any sense");

  // classical smagorinsky does only have constant parameter
  if (f3Parameter_->turb_mod_action_ == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
  {
    // this will be the y-coordinate of a point in the element interior
    // we will determine the element layer in which he is contained to
    // be able to do the output of visceff etc.
    double center = 0;

    for(int inode=0;inode<nen_;inode++)
    {
      center += xyze_( 1, inode );
    }
    center/=nen_;

    // node coordinates of plane to the element layer
    RefCountPtr<vector<double> > planecoords
      =
      turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_");

    bool found = false;
    for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
    {
      if(center<(*planecoords)[nlayer+1])
      {
        found = true;
        break;
      }
      nlayer++;
    }
    if (found ==false)
    {
      dserror("could not determine element layer");
    }
  }
  // --------------------------------------------------
  // Smagorinsky model with dynamic Computation of Cs
  //else if (physical_turbulence_model == "Dynamic_Smagorinsky")
  else if (f3Parameter_->turb_mod_action_ == INPAR::FLUID::dynamic_smagorinsky)
  {
    //turb_mod_action_ = Fluid3::dynamic_smagorinsky;

    // for turbulent channel flow, use averaged quantities
    if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
        ==
        "channel_flow_of_height_2")
    {
      RCP<vector<double> > averaged_LijMij
        =
        turbmodelparams.get<RCP<vector<double> > >("averaged_LijMij_");
      RCP<vector<double> > averaged_MijMij
        =
        turbmodelparams.get<RCP<vector<double> > >("averaged_MijMij_");

      //this will be the y-coordinate of a point in the element interior
      // here, the layer is determined in order to get the correct
      // averaged value from the vector of averaged (M/L)ijMij
      double center = 0;
      for(int inode=0;inode<nen_;inode++)
      {
        center += xyze_( 1, inode );
      }
      center/=nen_;

      RCP<vector<double> > planecoords
        =
        turbmodelparams.get<RCP<vector<double> > >("planecoords_");

      bool found = false;
      for (nlayer=0;nlayer < static_cast<int>((*planecoords).size()-1);)
      {
        if(center<(*planecoords)[nlayer+1])
        {
          found = true;
          break;
        }
        nlayer++;
      }
      if (found ==false)
      {
        dserror("could not determine element layer");
      }

      // Cs_delta_sq is set by the averaged quantities
      Cs_delta_sq = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;

      // clipping to get algorithm stable
      if (Cs_delta_sq<0)
      {
        Cs_delta_sq=0;
      }
    }
    else
    {
      // when no averaging was done, we just keep the calculated (clipped) value
      Cs_delta_sq = CsDeltaSq;
    }
  }
  return;
} // Fluid3Impl::GetTurbulenceParams


/*----------------------------------------------------------------------*
 |  calculation of (all-scale) subgrid viscosity               vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcSubgrVisc(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const double                            vol,
  double&                                 Cs,
  double&                                 Cs_delta_sq,
  double&                                 l_tau
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);
  //
  // SMAGORINSKY MODEL
  // -----------------
  //                                   +-                                 -+ 1
  //                               2   |          / h \           / h \    | -
  //    visc          = dens * lmix  * | 2 * eps | u   |   * eps | u   |   | 2
  //        turbulent           |      |          \   / ij        \   / ij |
  //                            |      +-                                 -+
  //                            |
  //                            |      |                                   |
  //                            |      +-----------------------------------+
  //                            |           'resolved' rate of strain
  //                    mixing length
  // -> either provided by dynamic modeling procedure and stored in Cs_delta_sq
  // -> or computed based on fixed Smagorinsky constant Cs:
  //             Cs = 0.17   (Lilly --- Determined from filter
  //                          analysis of Kolmogorov spectrum of
  //                          isotropic turbulence)
  //             0.1 < Cs < 0.24 (depending on the flow)
  //

  // compute (all-scale) rate of strain
  double rateofstrain = -1.0e30;
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  if (f3Parameter_->turb_mod_action_ == INPAR::FLUID::dynamic_smagorinsky)
  {
    // subgrid viscosity
    sgvisc_ = densaf_ * Cs_delta_sq * rateofstrain;

    // for evaluation of statistics: remember the 'real' Cs
    Cs = sqrt(Cs_delta_sq)/pow((vol),(1.0/3.0));
  }
  else
  {
    if (f3Parameter_->turb_mod_action_ == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
    {
      // since the Smagorinsky constant is only valid if hk is in the inertial
      // subrange of turbulent flows, the mixing length is damped in the
      // viscous near wall region using the van Driest damping function
      /*
                                     /         /   y+ \ \
                   lmix = Cs * hk * | 1 - exp | - ---- | |
                                     \         \   A+ / /
      */
      // A+ is a constant parameter, y+ the distance from the wall in wall
     // units
      const double A_plus = 26.0;
      double y_plus;

      // the integration point coordinate is defined by the isometric approach
      /*
                  +-----
                   \
              x =   +      N (x) * x
                   /        j       j
                  +-----
                  node j
      */

      LINALG::Matrix<nsd_,1> centernodecoord;
      centernodecoord.Multiply(xyze_,funct_);

      if (centernodecoord(1,0)>0) y_plus=(1.0-centernodecoord(1,0))/l_tau;
      else                        y_plus=(1.0+centernodecoord(1,0))/l_tau;

      //   lmix *= (1.0-exp(-y_plus/A_plus));
      // multiply with van Driest damping function
      Cs *= (1.0-exp(-y_plus/A_plus));
    }

    // get characteristic element length for Smagorinsky model for 2D and 3D
    // 3D: hk = V^1/3
    // 2D: hk = A^1/2
    const double hk = pow(vol,(1.0/dim));

    // mixing length set proportional to grid witdh: lmix = Cs * hk
    double lmix = Cs * hk;

    Cs_delta_sq = lmix * lmix;

    // subgrid viscosity
    sgvisc_ = densaf_ * Cs_delta_sq * rateofstrain;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of fine-scale subgrid viscosity                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcFineScaleSubgrVisc(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const LINALG::Matrix<nsd_,nen_>&        fsevelaf,
  const double                            vol,
  double&                                 Cs
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  //     // get characteristic element length for Smagorinsky model for 2D and 3D
  // 3D: hk = V^1/3
  // 2D: hk = A^1/2
  const double hk = pow(vol,(1.0/dim));

  if (f3Parameter_->fssgv_ == INPAR::FLUID::smagorinsky_all)
  {
    //
    // ALL-SCALE SMAGORINSKY MODEL
    // ---------------------------
    //                                      +-                                 -+ 1
    //                                  2   |          / h \           / h \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // compute (all-scale) rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

    fssgvisc_ = densaf_ * Cs * Cs * hk * hk * rateofstrain;
  }
  else if (f3Parameter_->fssgv_ == INPAR::FLUID::smagorinsky_small)
  {
    //
    // FINE-SCALE SMAGORINSKY MODEL
    // ----------------------------
    //                                      +-                                 -+ 1
    //                                  2   |          /    \          /   \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | fsu |   * eps | fsu |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // fine-scale rate of strain
    double fsrateofstrain = -1.0e30;
    fsrateofstrain = GetStrainRate(fsevelaf,derxy_,fsvderxy_);

    fssgvisc_ = densaf_ * Cs * Cs * hk * hk * fsrateofstrain;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter                     vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcStabParameter(
  const double               timefac,
  const double               vol
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  // get element-type constant for tau
  const double mk = DRT::ELEMENTS::MK<distype>();

  // get velocity norm
  const double vel_norm = velint_.Norm2();

  // ---------------------------------------------------------------
  // computation of stabilization parameter tau
  // ---------------------------------------------------------------
  //a)  Franca, Barrenechea, Valentin, Wall
  //b)  Franca, Barrenechea, Valentin, Wall without dt
  //c)  Bazilevs
  //e)  Codina
  //d)  Bazilevs without dt
  //f)  Oberai

  //----------------------------------------------------------------
  // definiton of element size 'strle' for tau_Mu
  //----------------------------------------------------------------

// a) streamlength (based on velocity vector at element centre) -> default
      // Tezeduyar 1992

      // normed velocity at element centre
      if (vel_norm>=1e-6) velino_.Update(1.0/vel_norm,velint_);
      else
      {
        velino_.Clear();
        velino_(0,0) = 1;
      }

      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_, velino_);
      const double val = tmp.Norm1();
      const double strle = 2.0/val;

//  b) volume-equival. diameter -> not default
      // warning: 3D formula
      // const double strle = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

//  c) cubic/square root of the element volume/area  -> not default

    //  strle = pow(vol,1/dim);
      //const double strle = pow(vol,1/dim);

   // else dserror("elment length calculation is not implemented for 1D flow");


  //----------------------------------------------------------------
  // definiton of element size 'hk' for tau_Mp and tau_C
  //----------------------------------------------------------------

    double hk = 0.0;
// a) streamlength (based on velocity vector at element centre) -> default

      // normed velocity at element centre
 /*     if (vel_norm>=1e-6) velino_.Update(1.0/vel_norm,velint_);
      else
      {
        velino_.Clear();
        velino_(0,0) = 1;
      }

      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_, velino_);
      const double val = tmp.Norm1();
      const double hk = 2.0/val;*/

//  b) volume-equival. diameter -> not default
      // warning: 3D formula
      if (nsd_==3)
        hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
        //const double hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

//  c) cubic/square root of the element volume/area  -> not default
      else if (nsd_==2)
        hk = pow(vol,1/dim);
        //const double hk = pow(vol,1/dim);

      else dserror("elment length calculation is not implemented for 1D flow");

  switch (f3Parameter_->whichtau_)
  {
  case INPAR::FLUID::tautype_franca_barrenechea_valentin_wall:
  {
    /*----------------------------------------------------- compute tau_Mu ---*/
    /* stability parameter definition according to

    Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
    element method for a generalized Stokes problem. Numerische
    Mathematik, Vol. 92, pp. 652-677, 2002.
    http://www.lncc.br/~valentin/publication.htm

    and:

    Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
    Finite Element Method for the Advective-Reactive-Diffusive
    Equation. Computer Methods in Applied Mechanics and Enginnering,
    Vol. 190, pp. 1785-1800, 2000.
    http://www.lncc.br/~valentin/publication.htm                   */

    // viscous : reactive forces
    const double re01 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(strle));

    // convective : viscous forces
    const double re02 = mk * densaf_ * vel_norm * strle / (2.0 * visceff_);

    const double xi01 = DMAX(re01,1.0);
    const double xi02 = DMAX(re02,1.0);

    tau_(0) = timefac*DSQR(strle)/(DSQR(strle)*densaf_*xi01+(4.0*timefac*visceff_/mk)*xi02);

    // compute tau_Mp
    //    stability parameter definition according to Franca and Valentin (2000)
    //                                       and Barrenechea and Valentin (2002)

    // viscous : reactive forces
    const double re11 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(hk));

    // convective : viscous forces
    const double re12 = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);

    const double xi11 = DMAX(re11,1.0);
    const double xi12 = DMAX(re12,1.0);

    /*
                  xi1,xi2 ^
                          |      /
                          |     /
                          |    /
                        1 +---+
                          |
                          |
                          |
                          +--------------> re1,re2
                              1
    */

    tau_(1) = timefac*DSQR(hk)/(DSQR(hk)*densaf_*xi11+(4.0*timefac*visceff_/mk)*xi12);

    /*------------------------------------------------------ compute tau_C ---*/
    // PhD thesis Wall (1999)
    /*
                      xi2 ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> Re2
                              1
    */

    // compute tau_c with hk instead of streamlength
    const double xi_tau_c = DMIN(re12,1.0);
    tau_(2) = densaf_ * vel_norm * hk * 0.5 * xi_tau_c;
  }
  break; // end franca_barrenechea_valentin_wall

  case INPAR::FLUID::tautype_franca_barrenechea_valentin_wall_wo_dt:
  {
    // compute tau_Mu
    /* convective : viscous forces */
    const double re_tau_mu = mk * densaf_ * vel_norm * strle / (2.0 * visceff_);
    const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
    tau_(0) = (DSQR(strle)*mk)/(4.0*visceff_*xi_tau_mu);

    // compute tau_Mp
    /* convective : viscous forces */
    const double re_tau_mp = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);
    const double xi_tau_mp = DMAX(re_tau_mp,1.0);
    tau_(1) = (DSQR(hk)*mk)/(4.0*visceff_*xi_tau_mp);

    // compute tau_C
    const double re_tau_c = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);
    const double xi_tau_c = DMIN(re_tau_c, 1.0);
    tau_(2) = 0.5 * densaf_ * vel_norm *hk * xi_tau_c;

    break;
  }

  case INPAR::FLUID::tautype_bazilevs:
  {
    /*

    tau_M: Bazilevs et al. (2007)
    e.g.: Variational multiscale residual-based turbulence modeling for
          large eddy simulation of incompressible flows
                                                                              1.0
                 +-                                                      -+ - ---
                 |        2                                               |   2.0
                 | 4.0*rho         n+1             n+1          2         |
          tau  = | -------  + rho*u     * G * rho*u     + C * mu  * G : G |
             M   |     2                  -                I        -   - |
                 |   dt                   -                         -   - |
                 +-                                                      -+

    */
    /*            +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
    */
    /*            +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
    */
    /*                                 +----
               n+1             n+1     \         n+1              n+1
          rho*u     * G * rho*u     =   +   rho*u    * G   * rho*u
                      -                /         i     -ij        j
                      -               +----        -
                                        i,j
    */

    double G;
    double normG = 0;
    double Gnormu = 0;

    const double dens_sqr = densaf_*densaf_;
    for (int nn=0;nn<nsd_;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*velint_(nn);
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1; mm<nsd_; ++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG+=G*G;
        Gnormu+=dens_sqr_velint_nn*G*velint_(rr);
      }
    }

    // definition of constant:
    // 12.0/m_k = 36.0 for linear elements and 144.0 for quadratic elements
    // (differently defined, e.g., in Akkerman et al. (2008))
    const double CI = 12.0/mk;

    tau_(0) = 1.0/(sqrt((4.0*dens_sqr)/(f3Parameter_->dt_*f3Parameter_->dt_)+Gnormu+CI*visceff_*visceff_*normG));
    tau_(1) = tau_(0);

    /*
      tau_C: Bazilevs et al. (2007), derived from fine-scale complement Shur
                                     operator of the pressure equation

                                  1.0
                    tau  = -----------------
                       C            /     \
                            tau  * | g * g |
                               M    \-   -/
    */
    /*           +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
    */
    /*           +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
    */
    double g;
    double normgsq = 0;
    for (int rr=0;rr<nsd_;++rr)
    {
      g = xji_(rr,0);
      for(int mm=1;mm<nsd_;++mm)
      {
        g += xji_(rr,mm);
      }
      normgsq += g*g;
    }

    tau_(2) = 1.0/(tau_(0)*normgsq);
  }
  break;  // end Bazilev

  case INPAR::FLUID::tautype_franca_barrenechea_valentin_codina:
  {
    /*----------------------------------------------------- compute tau_Mu ---*/
    /* stability parameter definition according to

    Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
    element method for a generalized Stokes problem. Numerische
    Mathematik, Vol. 92, pp. 652-677, 2002.
    http://www.lncc.br/~valentin/publication.htm

    and:

    Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
    Finite Element Method for the Advective-Reactive-Diffusive
    Equation. Computer Methods in Applied Mechanics and Enginnering,
    Vol. 190, pp. 1785-1800, 2000.
    http://www.lncc.br/~valentin/publication.htm                   */

    /* viscous : reactive forces */
    const double re01 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(strle));

    /* convective : viscous forces */
    const double re02 = mk * densaf_ * vel_norm * strle / (2.0 * visceff_);

    const double xi01 = DMAX(re01,1.0);
    const double xi02 = DMAX(re02,1.0);

    tau_(0) = timefac*DSQR(strle)/(DSQR(strle)*densaf_*xi01+(4.0*timefac*visceff_/mk)*xi02);

    if(nsd_==2)
      //TODO: in 2D hk and strlng was defined in the same way
      tau_(1) = tau_(0);

    else if(nsd_==3)
    {

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)

      /* viscous : reactive forces */
      const double re11 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(hk));

      /* convective : viscous forces */
      const double re12 = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);

      const double xi11 = DMAX(re11,1.0);
      const double xi12 = DMAX(re12,1.0);

      /*
                    xi1,xi2 ^
                            |      /
                            |     /
                            |    /
                          1 +---+
                            |dens
                            |
                            |
                            +--------------> re1,re2
                                1
      */
      tau_(1) = timefac*DSQR(hk)/(DSQR(hk)*densaf_*xi11+(4.0*timefac*visceff_/mk)*xi12);
    }
    else
      dserror("stabilization 'codina' is not implemented for 1D flow");

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    tau_(2) = sqrt(DSQR(visceff_)+DSQR(0.5*densaf_*vel_norm*hk));
  }
  break; // end Codina


  case INPAR::FLUID::tautype_bazilevs_wo_dt:
  {
    /*
    tau_M by Bazilevs et al. (2007) adapted for the time-dependent subgrid-scale
    approach by Codina. Note that the density dependency is included in the ODE
    for the unresolved scales

                                                        1.0
                 +-                                -+ - ---
                 |                                  |   2.0
                 |  n+1      n+1          2         |
          tau  = | u     * Gu     + C * nu  * G : G |
             M   |         -         I        -   - |
                 |         -                  -   - |
                 +-                                -+

    */
    /*            +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
    */
    /*            +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
    */
    /*                        +----
          n+af        n+af     \     n+af          n+af
         u     * G * u      =   +   u     * G   * u
                 -             /     i      -ij    j
                 -            +----         -
                               i,j
    */

    double G;
    double normG = 0;
    double Gnormu = 0;

    for (int nn=0;nn<nsd_;++nn)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1; mm<nsd_; ++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG+=G*G;
        Gnormu+=velint_(nn)*G*velint_(rr);
      }
    }

    const double dens_sqr = densaf_*densaf_;

    // definition of constant:
    // 12.0/m_k = 36.0 for linear elements and 144.0 for quadratic elements
    // (differently defined, e.g., in Akkerman et al. (2008))
    const double CI = 12.0/mk;

    tau_(0) = 1.0/(sqrt(Gnormu+CI*(visceff_*visceff_/dens_sqr)*normG));
    tau_(1) = tau_(0);

    /*
      tau_C: Bazilevs et al. (2007), derived from fine-scale complement Shur
                                     operator of the pressure equation

                                  1.0
                    tau  = -----------------
                       C            /     \
                            tau  * | g * g |
                               M    \-   -/
    */
    /*           +-     -+   +-     -+   +-     -+
                 |       |   |       |   |       |
                 |  dr   |   |  ds   |   |  dt   |
            g  = |  ---  | + |  ---  | + |  ---  |
             i   |  dx   |   |  dx   |   |  dx   |
                 |    i  |   |    i  |   |    i  |
                 +-     -+   +-     -+   +-     -+
    */
    /*           +----
                  \
         g * g =   +   g * g
         -   -    /     i   i
                 +----
                   i
    */
    double g;
    double normgsq = 0;
    for (int rr=0;rr<nsd_;++rr)
    {
      g = xji_(rr,0);
      for(int mm=1;mm<nsd_;++mm)
      {
        g += xji_(rr,mm);
      }
      normgsq += g*g;
    }

    tau_(2) = 1.0/(tau_(0)*normgsq);
  }
  break;  // end Bazilevs for tds

  case INPAR::FLUID::tautype_oberai:
  {
    /*----------------------------------------------------- compute tau_Mu ---*/
    /* stability parameter definition according to Oberai */

    // compute element Reynolds number
    // (caution: element Reynolds number defined here without division by 2)
    const double ere = densaf_ * vel_norm * hk / visceff_;

    // compute CFL number
    const double cfl = vel_norm * f3Parameter_->dt_ / hk;

    // compute k_bar
    const double pp = exp(ere/200.0);
    const double pm = exp(-ere/200.0);
    const double k_bar = ere*(1.0 + 8.0*((pp-pm)/(pp+pm))/(1.0+((pp-pm)/(pp+pm))));

    // compute factors gamma, kappa, alpha and mu
    const double kp = exp(k_bar);
    const double km = exp(-k_bar);
    const double fac_gamma = exp(k_bar*cfl*(1.0-(k_bar/ere)));
    const double fac_kappa = (2.0/(kp+km))-1.0;
    double fac_alpha = 1.0;
    if (k_bar < 700.0) fac_alpha = (kp-km)/(kp+km);
    const double fac_mu = (2.0*(2.0/(kp+km))+1.0)/3.0;

    // compute non-dimensional form of tau
    tau_(0) = (fac_mu*(fac_gamma-1.0)-(fac_gamma+1.0)*cfl*((fac_alpha/2.0)+(fac_kappa/ere)))/(fac_alpha*(fac_gamma-1.0)+cfl*fac_kappa*(fac_gamma+1.0));

    // compute dimensional (final) form of tau
    tau_(0) *= hk/vel_norm;
    tau_(1) = tau_(0);

    /*------------------------------------------------------ compute tau_C ---*/
    // PhD thesis Wall (1999)
    /*
                      xi2 ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> Re2
                              1
    */
    // convective : viscous forces
    const double re02 = mk * densaf_ * vel_norm * strle / (2.0 * visceff_);

    const double xi_tau_c = DMIN(re02,1.0);
    tau_(2) = densaf_ * vel_norm * hk * 0.5 * xi_tau_c;
  }
  break; // end oberai

  default: dserror("unknown definition of tau\n %i  ", f3Parameter_->whichtau_);
  }  // end switch

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcDivEps(
    const LINALG::Matrix<nsd_,nen_>&      evelaf)
{
  /*--- viscous term: div(epsilon(u)) --------------------------------*/
  /*   /                                                \
       |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
     1 |                                                |
     - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
     2 |                                                |
       |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
       \                                                /

       with N_x .. x-line of N
       N_y .. y-line of N                                             */

  /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
  /*   /                            \
       |  N_x,xx + N_y,yx + N_z,zx  |
     1 |                            |
  -  - |  N_x,xy + N_y,yy + N_z,zy  |
     3 |                            |
       |  N_x,xz + N_y,yz + N_z,zz  |
       \                            /

         with N_x .. x-line of N
         N_y .. y-line of N                                             */

  // set visc_old to zero
  visc_old_.Clear();

  double prefac;
  if(f3Parameter_->physicaltype_ == INPAR::FLUID::loma)
  //if(loma_)
  {
    prefac = 1.0/3.0;
    derxy2_.Scale(prefac);
 }
  else prefac = 1.0;

  if (nsd_==3)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      double sum = (derxy2_(0,inode)+derxy2_(1,inode)+derxy2_(2,inode))/prefac;
      viscs2_(0,inode) = 0.5 * (sum + derxy2_(0,inode));
      viscs2_(1,inode) = 0.5 *  derxy2_(3,inode);
      viscs2_(2,inode) = 0.5 *  derxy2_(4,inode);
      viscs2_(3,inode) = 0.5 *  derxy2_(3,inode);
      viscs2_(4,inode) = 0.5 * (sum + derxy2_(1,inode));
      viscs2_(5,inode) = 0.5 *  derxy2_(5,inode);
      viscs2_(6,inode) = 0.5 *  derxy2_(4,inode);
      viscs2_(7,inode) = 0.5 *  derxy2_(5,inode);
      viscs2_(8,inode) = 0.5 * (sum + derxy2_(2,inode));

      for (int idim=0; idim<nsd_; ++idim)
      {
        const int nsd_idim = idim*nsd_;
        for (int jdim=0; jdim<nsd_; ++jdim)
        {
          visc_old_(idim) += viscs2_(nsd_idim+jdim,inode)*evelaf(jdim,inode);
        }
      }
    }
  }
  else if (nsd_==2)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
    double sum = (derxy2_(0,inode)+derxy2_(1,inode))/prefac;
    viscs2_(0,inode) = 0.5 * (sum + derxy2_(0,inode));
    viscs2_(1,inode) = 0.5 * derxy2_(2,inode);
    viscs2_(2,inode) = 0.5 * derxy2_(2,inode);
    viscs2_(3,inode) = 0.5 * (sum + derxy2_(1,inode));

    for (int idim=0; idim<nsd_; ++idim)
    {
      const int nsd_idim = idim*nsd_;
      for (int jdim=0; jdim<nsd_; ++jdim)
      {
        visc_old_(idim) += viscs2_(nsd_idim+jdim,inode)*evelaf(jdim,inode);
      }
    }
    }
  }
  else dserror("Epsilon(N) is not implemented for the 1D case");

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::GetResidualMomentumEq(
    const LINALG::Matrix<nsd_,nen_>&              eaccam,
    const double &                                timefac,
    double &                                      rhsfac)
{
  if (f3Parameter_->is_genalpha_)
  {
    // rhs of momentum equation: density*bodyforce at n+alpha_F
    rhsmom_.Update(densaf_,bodyforce_,0.0);

    // get acceleration at time n+alpha_M at integration point
    accint_.Multiply(eaccam,funct_);

    // evaluate momentum residual once for all stabilization right hand sides
    for (int rr=0;rr<nsd_;++rr)
    {
      momres_old_(rr) =
        densam_*accint_(rr)
        +
        densaf_*conv_old_(rr)
        +
        gradp_(rr)
        -
        2*visceff_*visc_old_(rr)
        -
        densaf_*bodyforce_(rr);
    }
  }
  else
  {
    // rhs of momentum equation:
    // density*timefac*bodyforce at n+1 + density*histmom at n

    // in the case of a Boussinesq approximation:
    //
    //                    f = (rho - rho_0)/rho_0 *g
    //
    // else:
    //
    //                           f = rho * g

    if (f3Parameter_->physicaltype_ == INPAR::FLUID::boussinesq)
    //if(boussinesq_)
    {
      rhsmom_.Update(densn_,histmom_,deltadens_*timefac,bodyforce_);
    }
    else
    {
      rhsmom_.Update(densn_,histmom_,densaf_*timefac,bodyforce_);
    }

    // modify integration factor for Galerkin rhs
    rhsfac *= timefac;

    // evaluate momentum residual once for all stabilization right hand sides

    // stationary  : momres_old = theta ( ... ) -  bodyforce_
    // instationary: momres_old = u_(n+1) + theta ( ... ) - histmom_ - bodyforce_
    for (int rr=0;rr<nsd_;++rr)
    {
      momres_old_(rr) =
        timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))
        -
        rhsmom_(rr);
    }

    if (f3Parameter_->is_stationary_ == false)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) += densaf_*velint_(rr);
      }
    }
  }

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::UpdateSubscaleVelocity(
    double &        fac1,
    double &        fac2,
    double &        fac3,
    double &        facMtau,
    int    &        iquad,
    double *        saccn,
    double *        sveln,
    double *        svelnp
  )
{
  if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
  {
    /*-----------------------------------------------------------------*
     *                                                                 *
     *                    classical subgrid closure                    *
     *                                                                 *
     *-----------------------------------------------------------------*/
    if (f3Parameter_->cross_    != INPAR::FLUID::cross_stress_stab_none or
        f3Parameter_->reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
    {
      if (f3Parameter_->is_genalpha_)
      {
        // compute subgrid-scale velocity
        sgvelint_.Update(-tau_(1),momres_old_,0.0);
      }
      else
      {
        // compute subgrid-scale velocity
        sgvelint_.Update(-(tau_(1)/f3Parameter_->dt_),momres_old_,0.0);
      }
    }
  } // end classical subgrid closure
  else
  {
    /*-----------------------------------------------------------------*
     *                                                                 *
     *              time dependent subgrid scale closure               *
     *                                                                 *
     *-----------------------------------------------------------------*/
    if(f3Parameter_->is_stationary_)
    {
      dserror("there is no time dependent subgrid scale closure for stationary problems\n");
    }
    if ( saccn==NULL or sveln==NULL or svelnp==NULL )
    {
      dserror( "no subscale array provided" );
    }

    double alphaF = f3Parameter_->alphaF_;
    double alphaM = f3Parameter_->alphaM_;
    double gamma  = f3Parameter_->gamma_;
    double dt     = f3Parameter_->dt_;

    /*
                                            1.0
       facMtau =  -------------------------------------------------------
                     n+aM                      n+aF
                  rho     * alphaM * tauM + rho     * alphaF * gamma * dt
    */
    facMtau = 1.0/(densam_*alphaM*tau_(1)+densaf_*f3Parameter_->afgdt_);

    /*
       factor for old subgrid velocities:

                 n+aM                      n+aF
       fac1 = rho     * alphaM * tauM + rho     * gamma * dt * (alphaF-1)
    */
    fac1=(densam_*alphaM*tau_(1)+densaf_*gamma*dt*(alphaF-1.0))*facMtau;
    /*
      factor for old subgrid accelerations

                 n+aM
       fac2 = rho     * tauM * dt * (alphaM-gamma)
    */
    fac2=(densam_*dt*tau_(1)*(alphaM-gamma))*facMtau;
    /*
      factor for residual in current subgrid velocities:

       fac3 = gamma * dt * tauM
    */
    fac3=(gamma*dt*tau_(1))*facMtau;

    // if no generalised alpha time integration is used, the momres_old_
    // contains a scaled residual
    if (!f3Parameter_->is_genalpha_)
    {
      dserror("the time-dependent subgrid closure requires a genalpha time integration\n");
    }

    /*         +-                                       -+
        ~n+1   |        ~n           ~ n            n+1  |
        u    = | fac1 * u  + fac2 * acc  -fac3 * res     |
         (i)   |                                    (i)  |
               +-                                       -+
    */

    /* compute the intermediate value of subscale velocity

            ~n+af            ~n+1                   ~n
            u     = alphaF * u     + (1.0-alphaF) * u
             (i)              (i)

    */

    for (int rr=0;rr<nsd_;++rr)
    {
//       ele->UpdateSvelnpInOneDirection(
//         fac1           ,
//         fac2           ,
//         fac3           ,
//         momres_old_(rr),
//         f3Parameter_->alphaF_        ,
//         rr             ,
//         iquad          ,
//         sgvelint_(rr)  );

      int pos = rr + nsd_*iquad;

      /*
       *  ~n+1           ~n           ~ n            n+1
       *  u    =  fac1 * u  + fac2 * acc  -fac3 * res
       *   (i)
       *
       */

      svelnp[pos] =
        fac1*sveln[pos]
        +
        fac2*saccn[pos]
        -
        fac3*momres_old_(rr);

      /* compute the intermediate value of subscale velocity
       *
       *          ~n+af            ~n+1                   ~n
       *          u     = alphaF * u     + (1.0-alphaF) * u
       *           (i)              (i)
       *
       */
      sgvelint_(rr) =
        alphaF      *svelnp[pos]
        +
        (1.0-alphaF)*sveln [pos];
    }
  } // end time dependent subgrid scale closure

  /*-------------------------------------------------------------------*
   *                                                                   *
   *       include computed subgrid velocity in convective term        *
   *                                                                   *
   *-------------------------------------------------------------------*/

  if (f3Parameter_->cross_    != INPAR::FLUID::cross_stress_stab_none or
      f3Parameter_->reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
  {
    // compute subgrid-scale convective operator
    sgconv_c_.MultiplyTN(derxy_,sgvelint_);
  }
  else
  {
    sgconv_c_.Clear();
  }
}




template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::GetResidualContinuityEq(
    const LINALG::Matrix<nsd_,nen_>&          eveln,
    const LINALG::Matrix<nen_,1>&             escaaf,
    const LINALG::Matrix<nen_,1>&             escaam,
    const LINALG::Matrix<nen_,1>&             escadtam,
    const double &                            timefac)
{
  if (f3Parameter_->is_genalpha_)
  {
    // "incompressible" part of continuity residual: velocity divergence
    conres_old_ = vdiv_;

    if(f3Parameter_->physicaltype_ == INPAR::FLUID::loma)
    //if(loma_)
    {
      // time derivative of scalar at n+alpha_M
      const double tder_sca = funct_.Dot(escadtam);

      // gradient of scalar value at n+alpha_F
      grad_scaaf_.Multiply(derxy_,escaaf);

      // convective scalar term at n+alpha_F
      conv_scaaf_ = velint_.Dot(grad_scaaf_);

      // add subgrid-scale velocity part also to convective scalar term
      // -> currently not taken into account
      /*if (cross    != Fluid3::cross_stress_stab_none or
          reynolds != Fluid3::reynolds_stress_stab_none)
        conv_scaaf_ += sgvelint_.Dot(grad_scaaf_);*/

      /*

               /                                                dp   \
              |         1     / dT     /         \   \     1      th  |
              |    q , --- * | ---- + | u o nabla | T | - --- * ----  |
              |         T     \ dt     \         /   /    p      dt   |
               \                                           th        /
                      +---------------------------------------------+
                                        rhscon_
      */

      // rhs of continuity equation (only relevant for low-Mach-number flow)
      rhscon_ = scadtfac_*tder_sca + scaconvfacaf_*conv_scaaf_ + thermpressadd_;

      // residual of continuity equation
      conres_old_ -= rhscon_;
    }
  }
  else
  {
    // "incompressible" part of continuity residual: velocity divergence
    conres_old_ = timefac*vdiv_;

    if (f3Parameter_->physicaltype_ == INPAR::FLUID::loma or f3Parameter_->physicaltype_ == INPAR::FLUID::varying_density)
    //if(loma_ and varyingdensity_)
    {
      // get velocity at n
      velintn_.Multiply(eveln,funct_);

      // get velocity derivatives at n
      vderxyn_.MultiplyNT(eveln,derxy_);

      // velocity divergence at n
      //const double vdivn = vderxyn_(0, 0) + vderxyn_(1, 1) + vderxyn_(2, 2);

      double vdivn = 0.0;
      for (int idim = 0; idim<nsd_; ++idim)
      {
        vdivn += vderxyn_(idim,idim);
      }

      // scalar value at n+1
      const double scaaf = funct_.Dot(escaaf);

      // gradient of scalar value at n+1
      grad_scaaf_.Multiply(derxy_,escaaf);

      // convective scalar term at n+1
      conv_scaaf_ = velint_.Dot(grad_scaaf_);

      // scalar value at n
      const double scan = funct_.Dot(escaam);

      // gradient of scalar value at n
      grad_scan_.Multiply(derxy_,escaam);

      // convective scalar term at n
      conv_scan_ = velintn_.Dot(grad_scan_);

      // add subgrid-scale velocity part also to convective scalar term
      // (subgrid-scale velocity at n+1 also approximately used at n)
      // -> currently not taken into account
      /*if (cross    != Fluid3::cross_stress_stab_none or
          reynolds != Fluid3::reynolds_stress_stab_none)
      {
        conv_scaaf_ += sgvelint_.Dot(grad_scaaf_);
        conv_scan_  += sgvelint_.Dot(grad_scan_);
      }*/

      // rhs of continuity equation (only relevant for low-Mach-number flow)
      rhscon_ = scadtfac_*(scaaf-scan) +
                timefac*scaconvfacaf_*conv_scaaf_ +
                f3Parameter_->omtheta_*f3Parameter_->dt_*(scaconvfacn_*conv_scan_-vdivn) +
                thermpressadd_;

      // residual of continuity equation
      conres_old_ -= rhscon_;
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::LinGalMomResU(
                     LINALG::Matrix<nsd_*nsd_,nen_> &    lin_resM_Du,
                     const double &                      timefacfac)
{
  /*
      stationary                            cross-stress, part 1
       +-----+                             +-------------------+
       |     |                             |                   |

                 /       n+1       \        /      ~n+1       \
       rho*Du + |   rho*u   o nabla | Du + |   rho*u   o nabla | Du +
                 \      (i)        /        \      (i)        /

                 /                \  n+1
              + |   rho*Du o nabla | u
                 \                /   (i)
                |                        |
                +------------------------+
                        Newton
  */

  int idim_nsd_p_idim[nsd_];

  for (int idim = 0; idim <nsd_; ++idim)
  {
    idim_nsd_p_idim[idim]=idim*nsd_+idim;
  }

  if (f3Parameter_->is_stationary_ == false)
  {
    const double fac_densam=fac_*densam_;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double v=fac_densam*funct_(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }
  }

  const double timefacfac_densaf=timefacfac*densaf_;

  for (int ui=0; ui<nen_; ++ui)
  {
    const double v=timefacfac_densaf*conv_c_(ui);

    for (int idim = 0; idim <nsd_; ++idim)
    {
      lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
    }
  }

  if(f3Parameter_->is_newton_)
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const double temp=timefacfac_densaf*funct_(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        const int idim_nsd=idim*nsd_;

        for(int jdim=0;jdim<nsd_;++jdim)
        {
          lin_resM_Du(idim_nsd+jdim,ui)+=temp*vderxy_(idim,jdim);
        }
      }
    }
  }

  if(f3Parameter_->cross_==INPAR::FLUID::cross_stress_stab)
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const double v=timefacfac_densaf*sgconv_c_(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }
  }

  return;
}




template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::LinGalMomResU_subscales(
            LINALG::Matrix<nen_*nsd_,nen_>      estif_p_v,
            LINALG::Matrix<nsd_*nsd_,nen_> &    lin_resM_Du,
            LINALG::Matrix<nsd_,1> &            resM_Du,
            const double &                      timefacfac,
            const double &                      facMtau)
{
  // rescale Galerkin residual of all terms which have not been
  // integrated by parts

  const double C_saccGAL=densaf_*f3Parameter_->afgdt_*facMtau;

  for (int ui=0; ui<nen_; ++ui)
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int idim_nsd=idim*nsd_;

      for(int jdim=0;jdim<nsd_;++jdim)
      {
        lin_resM_Du(idim_nsd+jdim,ui)*=C_saccGAL;
      }
    }
  }

  // include all contributions which have been integrated by parts
  // and thus can not be rescaled

  /* viscous term (intermediate) */
  /*  factor:
                                rhoaM*alphaM*tauM                 gamma*dt
          2*nu*alphaF*---------------------------------------,  * --------
                      rhoaM*alphaM*tauM+rhoaf*alphaF*gamma*dt      alphaM


             /                         \
            |               /    \      |
            |  nabla o eps | Dacc | , v |
            |               \    /      |
             \                         /

  */

  if (is_higher_order_ele_)
  {
    const double v = 2.0*visceff_*timefacfac*(1.0-C_saccGAL);
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int nsd_idim=nsd_*idim;

      for(int jdim=0;jdim<nsd_;++jdim)
      {
        const int nsd_idim_p_jdim=nsd_idim+jdim;

        for (int ui=0; ui<nen_; ++ui)
        {
          lin_resM_Du(nsd_idim_p_jdim,ui)+=v*viscs2_(nsd_idim_p_jdim, ui);
        }
      }
    }
  }

  /*  factor:
                              rhoaM*alphaM*tauM                gamma*dt
          alphaF * ---------------------------------------,  * --------
                   rhoaM*alphaM*tauM+rhoaF*alphaF*gamma*dt      alphaM

                       /               \
                      |                 |
                      |  nabla Dp ,  v  |
                      |                 |
                       \               /
  */
  for (int ui=0; ui<nen_; ++ui)
  {
    const double v=(1.0-C_saccGAL)*timefacfac;
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = nsd_*vi;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_p_v(fvi + idim,ui)-= v*derxy_(idim,ui)*funct_(vi);
      }
    }
  }

  /*  factor: +1

           /                       \
          |     n+am    ~ n+am      |
          |  rho     * acc     , v  |
          |               (i)       |
           \                       /


         using
                                  n+af             /
           n+am    ~ n+am      rho        ~n+af   |    n+am      n+am
        rho     * acc     = - --------- * u     - | rho     * acc     +
                     (i)           n+af    (i)    |              (i)
                               tau_M               \

                                  n+af    / n+af        \   n+af            n+1
                             + rho     * | c     o nabla | u     + nabla o p    -
                                          \ (i)         /   (i)             (i)

                                                        / n+af \
                             - 2 * mu * grad o epsilon | u      | -
                                                        \ (i)  /
                                               \
                                  n+af    n+af  |
                             - rho     * f      |
                                                |
                                               /
  */
  for(int idim = 0; idim <nsd_; ++idim)
  {
    resM_Du(idim)=fac_*(-densaf_*sgvelint_(idim)/tau_(1)-momres_old_(idim));
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::InertiaAndConvectionGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &        lin_resM_Du,
    LINALG::Matrix<nsd_,1> &                resM_Du,
    const double &                          rhsfac)
{
  /* inertia (contribution to mass matrix) if not is_stationary */
  /*
            /              \
           |                |
           |    rho*Du , v  |
           |                |
            \              /
  */
  /* convection, convective part (convective form) */
  /*
            /                             \
           |  /       n+1       \          |
           | |   rho*u   o nabla | Du , v  |
           |  \      (i)        /          |
            \                             /
  */
  /*  convection, reactive part (convective form)
            /                               \
           |  /                \   n+1       |
           | |  rho*Du o nabla  | u     , v  |
           |  \                /   (i)       |
            \                               /
  */
  if(f3Parameter_->is_newton_ || (is_higher_order_ele_ && f3Parameter_->tds_==INPAR::FLUID::subscales_time_dependent))
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        const int idim_nsd=idim*nsd_;

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = nsd_*vi;

          const int fvi_p_idim = fvi+idim;

          for (int jdim= 0; jdim<nsd_;++jdim)
          {
            estif_u(fvi_p_idim,fui+jdim) += funct_(vi)*lin_resM_Du(idim_nsd+jdim,ui);
          } // end for (jdim)
        } // end for (idim)
      } //vi
    } // ui
  }
  else
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi;

        for (int idim = 0; idim <nsd_; ++idim)
        {
          estif_u(fvi+idim,fui+idim) += funct_(vi)*lin_resM_Du(idim*nsd_+idim,ui);
        } // end for (idim)
      } //vi
    } // ui
  }

  // inertia terms on the right hand side for instationary fluids
  if (f3Parameter_->is_stationary_ == false)
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      if (f3Parameter_->is_genalpha_)
      {
        resM_Du(idim)+=fac_*densam_*accint_(idim);
      }
      else
      {
        resM_Du(idim)+=fac_*densaf_*velint_(idim);
      }
    }
  }  // end if(stationary)

  for (int idim = 0; idim <nsd_; ++idim)
  {
    resM_Du(idim)+=rhsfac*densaf_*conv_old_(idim);
  }  // end for(idim)

  for (int vi=0; vi<nen_; ++vi)
  {
    for(int idim = 0; idim <nsd_; ++idim)
    {
      velforce(idim,vi)-=resM_Du(idim)*funct_(vi);
    }
  }
  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::ViscousGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nsd_> &             viscstress,
    const double &                          timefacfac,
    const double &                          rhsfac)
{
  const double visceff_timefacfac = visceff_*timefacfac;

  /* viscosity term */
  /*
                   /                        \
                  |       /  \         / \   |
            2 mu  |  eps | Du | , eps | v |  |
                  |       \  /         \ /   |
                   \                        /
  */

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi   = nsd_*vi;

    for (int jdim= 0; jdim<nsd_;++jdim)
    {
      const double temp=visceff_timefacfac*derxy_(jdim,vi);

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = nsd_*ui;

        for (int idim = 0; idim <nsd_; ++idim)
        {
          const int fvi_p_idim = fvi+idim;

          estif_u(fvi_p_idim,fui+jdim) += temp*derxy_(idim, ui);
        } // end for (jdim)
      } // end for (idim)
    } // ui
  } //vi

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi   = nsd_*vi;

    for (int jdim= 0; jdim<nsd_;++jdim)
    {
      const double temp=visceff_timefacfac*derxy_(jdim,vi);

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = nsd_*ui;

        for (int idim = 0; idim <nsd_; ++idim)
        {
          const int fvi_p_idim = fvi+idim;

          estif_u(fvi_p_idim,fui+idim) += temp*derxy_(jdim, ui);
        } // end for (jdim)
      } // end for (idim)
    } // ui
  } //vi

  const double v = visceff_*rhsfac;

  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      viscstress(idim,jdim)=v*(vderxy_(jdim,idim)+vderxy_(idim,jdim));
    }
  }

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::ContStab_and_ViscousTermRhs(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nsd_,nsd_> &               viscstress,
    const double &                            timefac,
    const double &                            timefacfac,
    const double &                            tau_C,
    const double &                            rhsfac)
{
  // In the case no continuity stabilization and no LOMA:
  // the factors 'conti_stab_and_vol_visc_fac' and 'conti_stab_and_vol_visc_rhs' are zero
  // therefore there is no contribution to the element stiffness matrix and
  // the viscous stress tensor is NOT altered!!
  //
  // ONLY
  // the rhs contribution of the viscous term is added!!

  double conti_stab_and_vol_visc_fac=0.0;
  double conti_stab_and_vol_visc_rhs=0.0;

  if(f3Parameter_->cstab_ == INPAR::FLUID::continuity_stab_yes)
  {
    conti_stab_and_vol_visc_fac+=timefac*tau_C;
    conti_stab_and_vol_visc_rhs-=tau_C*conres_old_;
  }
  if (f3Parameter_->physicaltype_ == INPAR::FLUID::loma)
  //if(loma_)
  {
    conti_stab_and_vol_visc_fac-=(2.0/3.0)*visceff_*timefacfac;
    conti_stab_and_vol_visc_rhs+=(2.0/3.0)*visceff_*rhsfac*vdiv_;
  }

  /* continuity stabilisation on left hand side */
  /*
              /                        \
             |                          |
        tauC | nabla o Du  , nabla o v  |
             |                          |
              \                        /
  */
  /* viscosity term - subtraction for low-Mach-number flow */
  /*
             /                             \             /                        \
            |  1                      / \   |     2 mu  |                          |
     - 2 mu |  - (nabla o u) I , eps | v |  | = - ----- | nabla o Du  , nabla o v  |
            |  3                      \ /   |       3   |                          |
             \                             /             \                        /
  */
  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui;

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int fui_p_idim = fui+idim;
      const double v0 = conti_stab_and_vol_visc_fac*derxy_(idim,ui);
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi;

        for(int jdim=0;jdim<nsd_;++jdim)
        {
          estif_u(fvi+jdim,fui_p_idim) += v0*derxy_(jdim, vi) ;
        }
      }
    } // end for(idim)
  }

  for(int idim=0;idim<nsd_;++idim)
  {
    viscstress(idim,idim)-=conti_stab_and_vol_visc_rhs;
  }

  //computation of right-hand-side viscosity term
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        /* viscosity term on right-hand side */
        velforce(idim,vi)-= viscstress(idim,jdim)*derxy_(jdim,vi);
      }
    }
  }

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::PressureGalPart(
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfac,
    const double &                            rhsfac,
    const double &                            press)
{
  for (int ui=0; ui<nen_; ++ui)
  {
    const double v = -timefacfac*funct_(ui);
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = nsd_*vi;
      /* pressure term */
      /*
           /                \
          |                  |
          |  Dp , nabla o v  |
          |                  |
           \                /
      */
      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_p_v(fvi + idim,ui) += v*derxy_(idim, vi);
      }
    }
  }

  const double pressfac = press*rhsfac;

  for (int vi=0; vi<nen_; ++vi)
  {
    /* pressure term on right-hand side */
    for (int idim = 0; idim <nsd_; ++idim)
    {
      velforce(idim,vi)+= pressfac*derxy_(idim, vi) ;
    }
  }  //end for(idim)

  return;
}




template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::ContinuityGalPart(
    LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,
    LINALG::Matrix<nen_,1> &                  preforce,
    const double &                            timefacfac,
    const double &                            rhsfac)
{
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = timefacfac*funct_(vi);
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        /* continuity term */
        /*
             /                \
            |                  |
            | nabla o Du  , q  |
            |                  |
             \                /
        */
        estif_q_u(vi,fui+idim) += v*derxy_(idim,ui);
      }
    }
  }  // end for(idim)

  const double rhsfac_vdiv = -rhsfac * vdiv_;
  for (int vi=0; vi<nen_; ++vi)
  {
    // continuity term on right-hand side
    preforce(vi) += rhsfac_vdiv*funct_(vi);
  }

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::BodyForceRhsTerm(
    LINALG::Matrix<nsd_,nen_> &               velforce)
{
  for (int idim = 0; idim <nsd_; ++idim)
  {
    const double scaled_rhsmom=fac_*rhsmom_(idim);

    for (int vi=0; vi<nen_; ++vi)
    {
      velforce(idim,vi)+=scaled_rhsmom*funct_(vi);
    }
  }  // end for(idim)

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::ConservativeFormulation(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfac,
    const double &                            rhsfac)
{
  //----------------------------------------------------------------------
  // computation of additions to convection term (convective and
  // reactive part) for conservative form of convection term including
  // right-hand-side contribution
  //----------------------------------------------------------------------

  // TODO: cleaning

  for (int idim = 0; idim <nsd_; ++idim)
   {
     for (int ui=0; ui<nen_; ++ui)
     {
       const int fui   = nsd_*ui + idim;
       //const int fui   = 4*ui;
       //const int fuip  = fui+1;
       //const int fuipp = fui+2;
       double v = timefacfac*densaf_*funct_(ui)*vdiv_;
       if (f3Parameter_->physicaltype_ == INPAR::FLUID::loma) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
       //if (loma_) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
       // only with linear density-concentration correlation
       else if(f3Parameter_->physicaltype_ == INPAR::FLUID::varying_density) v += timefacfac*conv_scaaf_;
       //else if(varyingdensity_) v += timefacfac*conv_scaaf_;
       for (int vi=0; vi<nen_; ++vi)
       {
         const int fvi   = nsd_*vi + idim;
         //const int fvi   = 4*vi;
         //const int fvip  = fvi+1;
         //const int fvipp = fvi+2;
         /* convection, convective part (conservative addition) */
         /*
           /                                                \
           |      /              n+1    n+1           \      |
           |  Du | rho*nabla o u    +  u   *nabla rho | , v  |
           |      \             (i)     (i)          /       |
           \                                                 /
         */
         double v2 = v*funct_(vi) ;
         estif_u(fvi  , fui  ) += v2;
         //estif(fvi  , fui  ) += v2;
         //estif(fvip , fuip ) += v2;
         //estif(fvipp, fuipp) += v2;
       }
     }

     if (f3Parameter_->is_newton_)
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         const int fvi   = nsd_*vi + idim;
         //const int fvi   = 4*vi;
         //const int fvip  = fvi+1;
         //const int fvipp = fvi+2;
         const double v_idim = timefacfac*densaf_*velint_(idim)*funct_(vi);
         //const double v0 = timefacfac*densaf_*velint_(0)*funct_(vi);
         //const double v1 = timefacfac*densaf_*velint_(1)*funct_(vi);
         //const double v2 = timefacfac*densaf_*velint_(2)*funct_(vi);
         for (int ui=0; ui<nen_; ++ui)
         {
           const int fui   = nsd_*ui;
           //const int fui   = 4*ui;
           //const int fuip  = fui+1;
           //const int fuipp = fui+2;
           /*  convection, reactive part (conservative addition) */
           /*
             /                              \
             |  n+1  /               \      |
             | u    | rho*nabla o Du | , v  |
             |  (i)  \              /       |
             \                             /
           */
           for(int jdim=0; jdim<nsd_;++jdim)
           estif_u(fvi,  fui+jdim  ) += v_idim*derxy_(jdim, ui) ;

           //estif(fvi,  fui  ) += v0*derxy_(0, ui) ;
           //estif(fvi,  fuip ) += v0*derxy_(1, ui) ;
           //estif(fvi,  fuipp) += v0*derxy_(2, ui) ;
           //estif(fvip, fui  ) += v1*derxy_(0, ui) ;
           //estif(fvip, fuip ) += v1*derxy_(1, ui) ;
           //estif(fvip, fuipp) += v1*derxy_(2, ui) ;
           //estif(fvipp,fui  ) += v2*derxy_(0, ui) ;
           //estif(fvipp,fuip ) += v2*derxy_(1, ui) ;
           //estif(fvipp,fuipp) += v2*derxy_(2, ui) ;
         }
       }

       if (f3Parameter_->physicaltype_ == INPAR::FLUID::loma)
       //if(loma_)
       {
         for (int vi=0; vi<nen_; ++vi)
         {
           const int fvi   = nsd_*vi + idim;
           //const int fvi   = 4*vi;
           //const int fvip  = fvi+1;
           //const int fvipp = fvi+2;
           const double v_idim = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(idim)*velint_(idim)*funct_(vi);
           //const double v0 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(0)*velint_(0)*funct_(vi);
           //const double v1 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(1)*velint_(1)*funct_(vi);
           //const double v2 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(2)*velint_(2)*funct_(vi);
           for (int ui=0; ui<nen_; ++ui)
           {
             const int fui   = nsd_*ui;
             //const int fui   = 4*ui;
             //const int fuip  = fui+1;
             //const int fuipp = fui+2;
             /*  convection, reactive part (conservative addition) */
             /*
               /                           \
               |  n+1  /             \      |
               | u    | Du*nabla rho | , v  |
               |  (i)  \            /       |
               \                           /
             */
             for(int jdim=0;jdim<nsd_;++jdim)
               estif_u(fvi,  fui +jdim  ) += v_idim*funct_(ui) ;

             //estif(fvi,  fui  ) += v0*funct_(ui) ;
             //estif(fvi,  fuip ) += v0*funct_(ui) ;
             //estif(fvi,  fuipp) += v0*funct_(ui) ;
             //estif(fvip, fui  ) += v1*funct_(ui) ;
             //estif(fvip, fuip ) += v1*funct_(ui) ;
             //estif(fvip, fuipp) += v1*funct_(ui) ;
             //estif(fvipp,fui  ) += v2*funct_(ui) ;
             //estif(fvipp,fuip ) += v2*funct_(ui) ;
             //estif(fvipp,fuipp) += v2*funct_(ui) ;
           }
         }
       }
       if (f3Parameter_->physicaltype_ == INPAR::FLUID::varying_density)
       //if(varyingdensity_)
       {
         for (int vi=0; vi<nen_; ++vi)
         {
           const int fvi   = nsd_*vi + idim;
           //const int fvi   = 4*vi;
           //const int fvip  = fvi+1;
           //const int fvipp = fvi+2;
           const double v_idim = +timefacfac*grad_scaaf_(idim)*velint_(idim)*funct_(vi);
           //const double v0 = +timefacfac*grad_scaaf_(0)*velint_(0)*funct_(vi);
           //const double v1 = +timefacfac*grad_scaaf_(1)*velint_(1)*funct_(vi);
           //const double v2 = +timefacfac*grad_scaaf_(2)*velint_(2)*funct_(vi);
           for (int ui=0; ui<nen_; ++ui)
           {
             const int fui   = nsd_*ui;
             //const int fui   = 4*ui;
             //const int fuip  = fui+1;
             //const int fuipp = fui+2;
             /*  convection, reactive part (conservative addition) */
             /*
               /                           \
               |  n+1  /             \      |
               | u    | Du*nabla rho | , v  |
               |  (i)  \            /       |
               \                           /
             */
             for(int jdim=0;jdim<nsd_;++jdim)
               estif_u(fvi,  fui+jdim  ) += v_idim*funct_(ui) ;

           }
         }
       }
     }

     for (int vi=0; vi<nen_; ++vi)
     {
       //const int fvi   = 4*vi;
       /* convection (conservative addition) on right-hand side */
       double v = -rhsfac*densaf_*funct_(vi)*vdiv_;
       velforce(idim, vi    ) += v*velint_(idim) ;
       //eforce(fvi    ) += v*velint_(0) ;
       //eforce(fvi + 1) += v*velint_(1) ;
       //eforce(fvi + 2) += v*velint_(2) ;
     }

     if (f3Parameter_->physicaltype_ == INPAR::FLUID::loma)
     //if(loma_)
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         //const int fvi   = 4*vi;
         /* convection (conservative addition) on rhs for low-Mach-number flow */
         double v = rhsfac*densaf_*scaconvfacaf_*conv_scaaf_*funct_(vi);
         velforce(idim, vi    ) += v*velint_(idim) ;
         //eforce(fvi    ) += v*velint_(0) ;
         //eforce(fvi + 1) += v*velint_(1) ;
         //eforce(fvi + 2) += v*velint_(2) ;
       }
     }
     if (f3Parameter_->physicaltype_ == INPAR::FLUID::varying_density)
     //if(varyingdensity_)
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         //const int fvi   = 4*vi;
         /* convection (conservative addition) on rhs for low-Mach-number flow */
         double v = -rhsfac*conv_scaaf_*funct_(vi);
         velforce(idim, vi    ) += v*velint_(idim) ;
       }
     }
   }  // end for(idim)

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::LomaGalPart(
    LINALG::Matrix<nen_, nen_*nsd_> &       estif_q_u,
    LINALG::Matrix<nen_,1> &                preforce,
    const double &                          timefacfac,
    const double &                          rhsfac)
{
  //----------------------------------------------------------------------
  // computation of additional terms for low-Mach-number flow:
  // 2) additional rhs term of continuity equation
  //----------------------------------------------------------------------

  if (f3Parameter_->is_newton_)
  {
    const double timefacfac_scaconvfacaf=timefacfac*scaconvfacaf_;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui=nsd_*ui;

      const double timefacfac_scaconvfacaf_funct_ui=timefacfac_scaconvfacaf*funct_(ui);

      for(int jdim=0;jdim<nsd_;++jdim)
      {
        const double temp=timefacfac_scaconvfacaf_funct_ui*grad_scaaf_(jdim);

        for (int vi=0; vi<nen_; ++vi)
        {
          //const int fvippp= numdofpernode_*vi+nsd_;


          /*
              factor afgtd/am

                      /                    \
                1    |       /         \    |
               --- * |  q , | Du o grad | T |
                T    |       \         /    |
                      \                    /
          */
          estif_q_u(vi,fui+jdim) -= temp*funct_(vi);
        }
      }
    }
  } // end if (is_newton_)

  const double fac_rhscon = fac_*rhscon_;
  for (int vi=0; vi<nen_; ++vi)
  {
    /* additional rhs term of continuity equation */
    preforce(vi) += fac_rhscon*funct_(vi) ;
  }

  /*
  if(ele->Id()==100 )
  {
    printf("rhscon_       %17.12e\n",rhscon_);
  }
  */

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::StabLinGalMomResU(
    LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
    const double &                            timefacfac)
{

  /*
                 /       n+1       \        /                \  n+1
       rho*Du + |   rho*u   o nabla | Du + |   rho*Du o nabla | u   +
                 \      (i)        /        \                /   (i)

                             /  \
              + nabla o eps | Du |
                             \  /
  */
  if(f3Parameter_->tds_==INPAR::FLUID::subscales_time_dependent
     ||
     f3Parameter_->cross_==INPAR::FLUID::cross_stress_stab)
  {
    //----------------------------------------------------------------------
    /* GALERKIN residual was rescaled and cannot be reused; so rebuild it */

    lin_resM_Du.Clear();

    int idim_nsd_p_idim[nsd_];

    for (int idim = 0; idim <nsd_; ++idim)
    {
      idim_nsd_p_idim[idim]=idim*nsd_+idim;
    }

    if (f3Parameter_->is_stationary_ == false)
    {
      const double fac_densam=fac_*densam_;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double v=fac_densam*funct_(ui);

        for (int idim = 0; idim <nsd_; ++idim)
        {
          lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
        }
      }
    }

    const double timefacfac_densaf=timefacfac*densaf_;

    for (int ui=0; ui<nen_; ++ui)
    {
      // deleted +sgconv_c_(ui)
      const double v=timefacfac_densaf*conv_c_(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }

    if(f3Parameter_->is_newton_)
    {
//
//
// dr_j   d    /    du_j \          du_j         dN_B
// ----= ---- | u_i*----  | = N_B * ---- + u_i * ---- * d_jk
// du_k  du_k  \    dx_i /          dx_k         dx_i

      for (int ui=0; ui<nen_; ++ui)
      {
        const double temp=timefacfac_densaf*funct_(ui);

        for (int idim = 0; idim <nsd_; ++idim)
        {
          const int idim_nsd=idim*nsd_;

          for(int jdim=0;jdim<nsd_;++jdim)
          {
            lin_resM_Du(idim_nsd+jdim,ui)+=temp*vderxy_(idim,jdim);
          }
        }
      }
    }
  }

  if (is_higher_order_ele_)
  {
    const double v = -2.0*visceff_*timefacfac;
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int nsd_idim=nsd_*idim;

      for(int jdim=0;jdim<nsd_;++jdim)
      {
        const int nsd_idim_p_jdim=nsd_idim+jdim;

        for (int ui=0; ui<nen_; ++ui)
        {
          lin_resM_Du(nsd_idim_p_jdim,ui)+=v*viscs2_(nsd_idim_p_jdim, ui);
        }
      }
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::PSPG(
    LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,
    LINALG::Matrix<nen_,nen_> &               ppmat,
    LINALG::Matrix<nen_,1> &                  preforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
    const double &                            fac3,
    const double &                            timefacfac,
    const double &                            tau_Mp,
    const double &                            rhsfac)
{
  // conservative, stabilization terms are neglected (Hughes)

  /* pressure stabilisation:                                            */
  /*
              /                 \
             |  ~n+af            |
           - |  u     , nabla q  |
             |                   |
              \                 /
  */

    double scal_grad_q=0.0;

    if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
    {
      scal_grad_q=tau_(1);
    }
    else
    {
      scal_grad_q=f3Parameter_->alphaF_*fac3;
    }

    /* pressure stabilisation: inertia if not stationary*/
    /*
              /                  \
             |                    |
             |  rho*Du , nabla q  |
             |                    |
              \                  /
    */
    /* pressure stabilisation: convection, convective part */
    /*
              /                                   \
             |  /       n+1       \                |
             | |   rho*u   o nabla | Du , nabla q  |
             |  \      (i)        /                |
              \                                   /
    */
    /* pressure stabilisation: convection, reactive part if Newton */
    /*
              /                                   \
             |  /                \   n+1           |
             | |   rho*Du o nabla | u     , grad q |
             |  \                /   (i)           |
              \                                   /
    */
    /* pressure stabilisation: viscosity (-L_visc_u) */
    /*
              /                              \
             |               /  \             |
             |  nabla o eps | Du | , nabla q  |
             |               \  /             |
              \                              /
    */

    if (is_higher_order_ele_ || f3Parameter_->is_newton_)
    {
      for(int jdim=0;jdim<nsd_;++jdim)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui_p_jdim   = nsd_*ui + jdim;

          for(int idim=0;idim<nsd_;++idim)
          {
            const int nsd_idim=nsd_*idim;

            for (int vi=0; vi<nen_; ++vi)
            {
              const double temp_vi_idim=derxy_(idim,vi)*scal_grad_q;

              estif_q_u(vi,fui_p_jdim) += lin_resM_Du(nsd_idim+jdim,ui)*temp_vi_idim;
            } // jdim
          } // vi
        } // ui
      } //idim
    } // end if (is_higher_order_ele_) or (newton_)
    else
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for(int idim=0;idim<nsd_;++idim)
        {
          const int nsd_idim=nsd_*idim;

          const double temp_vi_idim=derxy_(idim, vi)*scal_grad_q;

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui_p_idim   = nsd_*ui + idim;

            estif_q_u(vi,fui_p_idim) += lin_resM_Du(nsd_idim+idim,ui)*temp_vi_idim;
          } // vi
        } // ui
      } //idim
    } // end if not (is_higher_order_ele_) nor (newton_)


    for (int ui=0; ui<nen_; ++ui)
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        const double v=timefacfac*derxy_(idim,ui)*scal_grad_q;

        for (int vi=0; vi<nen_; ++vi)
        {
          /* pressure stabilisation: pressure( L_pres_p) */
          /*
               /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
               \                    /
          */
          ppmat(vi,ui)+=v*derxy_(idim,vi);
        } // vi
      } // end for(idim)
    }  // ui

    if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        const double temp= tau_Mp*momres_old_(idim);

        for (int vi=0; vi<nen_; ++vi)
        {
          // pressure stabilisation
          preforce(vi) -= temp*derxy_(idim, vi);
        }
      } // end for(idim)
    }
    else
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        const double temp= rhsfac*sgvelint_(idim);

        for (int vi=0; vi<nen_; ++vi)
        {
          // pressure stabilisation
          preforce(vi) += temp*derxy_(idim, vi);
        }
      } // end for(idim)
    }

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::SUPG(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
    const double &                            fac3,
    const double &                            timefacfac,
    const double &                            timetauM)
{
  /*
                    /                                \
                   |  ~n+af    /     n+af       \     |
                 - |  u     , | rho*u    o nabla | v  |
                   |           \     (i)        /     |
                    \                                /
   */

     LINALG::Matrix<nsd_,1> temp;
     double supgfac;

     if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
     {
       supgfac=densaf_*tau_(0);
     }
     else
     {
       supgfac=densaf_*f3Parameter_->alphaF_*fac3;
     }

     LINALG::Matrix<nen_,1> supg_test;

     for (int vi=0; vi<nen_; ++vi)
     {
       supg_test(vi)=supgfac*conv_c_(vi);
     }

     if(f3Parameter_->reynolds_ == INPAR::FLUID::reynolds_stress_stab)
     {
       double reyfac;

       if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
       {
         reyfac=densaf_*tau_(1);
       }
       else
       {
         reyfac=densaf_*f3Parameter_->alphaF_*fac3;
       }

       for (int vi=0; vi<nen_; ++vi)
       {
         supg_test(vi)+=reyfac*sgconv_c_(vi);
       }
     }

     /* supg stabilisation: inertia if not stationary */
     /*
            /                                \
           |            /     n+1       \     |
           |  rho*Du , | rho*u   o nabla | v  |
           |            \     (i)       /     |
            \                                /
     */
     /* supg stabilisation: convective part ( L_conv_u) , convective term */
     /*
            /                                                     \
           |    /       n+1        \        /      n+1       \     |
           |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
           |    \       (i)        /        \      (i)       /     |
            \                                                     /
     */
     /* supg stabilisation: viscous part  (-L_visc_u) if is_higher_order_ele_ */
     /*
            /                                              \
           |               /  \    /       n+1        \     |
           |  nabla o eps | Du |, |   rho*u    o nabla | v  |
           |               \  /    \       (i)        /     |
            \                                              /
     */
     /* supg stabilisation: convective part ( L_conv_u) , reactive term if Newton */
     /*
            /                                                     \
           |    /       n+1        \        /     n+1        \     |
           |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
           |    \       (i)        /        \     (i)        /     |
            \                                                     /
     */
     if (is_higher_order_ele_ || f3Parameter_->is_newton_)
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         for(int idim=0;idim<nsd_;++idim)
         {
           const int nsd_idim=nsd_*idim;

           const int fvi_p_idim = nsd_*vi+idim;

           for(int jdim=0;jdim<nsd_;++jdim)
           {
             const int nsd_idim_p_jdim=nsd_idim+jdim;
             for (int ui=0; ui<nen_; ++ui)
             {
               const int fui_p_jdim   = nsd_*ui + jdim;

               estif_u(fvi_p_idim,fui_p_jdim) += lin_resM_Du(nsd_idim_p_jdim,ui)*supg_test(vi);
             } // jdim
           } // vi
         } // ui
       } //idim
     } // end if (is_higher_order_ele_) or (newton_)
     else
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         for(int idim=0;idim<nsd_;++idim)
         {
           const int fvi_p_idim = nsd_*vi+idim;

           const int nsd_idim=nsd_*idim;

           for (int ui=0; ui<nen_; ++ui)
           {
             const int fui_p_idim   = nsd_*ui + idim;

             estif_u(fvi_p_idim,fui_p_idim) += lin_resM_Du(nsd_idim+idim,ui)*supg_test(vi);
           } // ui
         } //idim
       } // vi
     } // end if not (is_higher_order_ele_) nor (newton_)

     /* supg stabilisation: pressure part  ( L_pres_p) */
     /*
              /                                    \
             |              /       n+1       \     |
             |  nabla Dp , |   rho*u   o nabla | v  |
             |              \       (i)       /     |
              \                                    /
     */
     for (int vi=0; vi<nen_; ++vi)
     {
       const double v = timefacfac*supg_test(vi);

       for (int idim = 0; idim <nsd_; ++idim)
       {
         const int fvi   = nsd_*vi + idim;

         for (int ui=0; ui<nen_; ++ui)
         {
           estif_p_v(fvi,ui) += v*derxy_(idim, ui);
         }
       }
     }  // end for(idim)

     /* supg stabilisation: inertia, linearisation of testfunction  */
     /*
                 /                                       \
                |         n+1       /              \      |
                |    rho*u      ,  | rho*Du o nabla | v   |
                |         (i)       \              /      |
                 \                                       /
     */
     /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
     /*
                 /                                                       \
                |    /       n+1        \   n+1     /              \      |
                |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
                |    \       (i)        /   (i)     \              /      |
                 \                                                       /
     */
     /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
     /*
                /                                     \
               |         n+1    /                \     |
               |  nabla p    , |   rho*Du o nabla | v  |
               |         (i)    \                /     |
                \                                     /
     */
     /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
     /*
                /                                               \
               |               / n+1 \    /               \      |
               |  nabla o eps | u     |, |  rho*Du o nabla | v   |
               |               \ (i) /    \               /      |
                \                                               /
     */
     /* supg stabilisation: bodyforce part, linearisation of test function */
     /*
                /                                      \
               |                  /               \     |
               |  rho*rhsint   , |  rho*Du o nabla | v  |
               |                  \               /     |
                \                                      /
     */
     if (f3Parameter_->is_newton_)
     {
       if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
       {
         for(int jdim=0;jdim<nsd_;++jdim)
         {
           temp(jdim)=densaf_*timetauM*momres_old_(jdim);
         }
       }
       else
       {
         for(int jdim=0;jdim<nsd_;++jdim)
         {
           temp(jdim)=-timefacfac*sgvelint_(jdim)*densaf_;
         }
       }
       /*
       if(ele->Id()==100)
       {
         printf("timetauM       %17.12e\n",timetauM);
         printf("densaf_        %17.12e\n",densaf_);
         printf("momres_old_(0) %17.12e\n",momres_old_(0));
         printf("momres_old_(1) %17.12e\n",momres_old_(1));
         printf("momres_old_(2) %17.12e\n",momres_old_(2));
         printf("temp(0)        %17.12e\n",temp(0));
         printf("temp(1)        %17.12e\n",temp(1));
         printf("temp(2)        %17.12e\n",temp(2));
       }
       */
       for(int jdim=0;jdim<nsd_;++jdim)
       {
         for (int vi=0; vi<nen_; ++vi)
         {
           const int fvi_p_jdim = nsd_*vi+jdim;

           for(int idim=0;idim<nsd_;++idim)
           {
             const double v=temp(jdim)*derxy_(idim,vi);

             for (int ui=0; ui<nen_; ++ui)
             {
               const int fui_p_idim   = nsd_*ui + idim;

               estif_u(fvi_p_jdim,fui_p_idim) += v*funct_(ui);
             } // jdim
           } // vi
         } // ui
       } //idim
     }

     // NOTE: Here we have a difference to the previous version of this
     // element! Before we did not care for the mesh velocity in this
     // term. This seems unreasonable and wrong.

     if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
     {
       for(int jdim=0;jdim<nsd_;++jdim)
       {
         temp(jdim)=fac_*momres_old_(jdim);
       }
     }
     else
     {
       for(int jdim=0;jdim<nsd_;++jdim)
       {
         temp(jdim)=-1.0/supgfac*fac_*densaf_*sgvelint_(jdim);
       }
     }

     for (int idim = 0; idim <nsd_; ++idim)
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         // supg stabilisation
         velforce(idim,vi) -= temp(idim)*supg_test(vi);
       }
     }  // end for(idim)

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::ViscStab(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
    const double &                            timefac,
    const double &                            tau_Mp,
    const double &                            vstabfac,
    const double &                            fac3)
{
  // viscous stabilization either on left hand side or on right hand side
   if (f3Parameter_->vstab_ == INPAR::FLUID::viscous_stab_gls || f3Parameter_->vstab_ == INPAR::FLUID::viscous_stab_usfem)
   {
     double two_visc_tau;

     if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
     {
       two_visc_tau      = vstabfac*2.0*visc_*tau_(1);
     }
     else
     {
       two_visc_tau      = vstabfac*2.0*visc_*f3Parameter_->alphaF_*fac3;
     }


     /* viscous stabilisation, inertia part if not stationary */
     /*
                /                        \
               |                          |
           +/- |    rho*Du , div eps (v)  |
               |                          |
                \                        /
     */
     /* viscous stabilisation, convective part, convective type */
     /*
              /                                      \
             |  /       n+1       \                   |
         +/- | |   rho*u   o nabla | Du , div eps (v) |
             |  \       (i)       /                   |
              \                                      /
     */
     /* viscous stabilisation, reactive part of convection */
     /*
              /                                       \
             |  /                \   n+1               |
         +/- | |   rho*Du o nabla | u    , div eps (v) |
             |  \                /   (i)               |
              \                                       /
     */
     /* viscous stabilisation, viscous part (-L_visc_u) */
     /*
              /                                 \
             |               /  \                |
        -/+  |  nabla o eps | Du | , div eps (v) |
             |               \  /                |
              \                                 /
     */
     for(int jdim=0;jdim<nsd_;++jdim)
     {
       for (int ui=0; ui<nen_; ++ui)
       {
         const int fui_p_jdim   = nsd_*ui + jdim;

         for(int idim=0;idim<nsd_;++idim)
         {
           for(int kdim=0;kdim<nsd_;++kdim)
           {
             for (int vi=0; vi<nen_; ++vi)
             {
               const int fvi_p_idim = nsd_*vi+idim;

               estif_u(fvi_p_idim,fui_p_jdim) += two_visc_tau*lin_resM_Du(nsd_*kdim+jdim,ui)*viscs2_(nsd_*idim+kdim,vi);
             } // vi
           } // kdim
         } // idim
       } // ui
     } //jdim


     /* viscous stabilisation, pressure part ( L_pres_p) */
     /*
              /                        \
             |                          |
        +/-  |  nabla Dp , div eps (v)  |
             |                          |
              \                        /
     */
     for (int idim=0;idim<nsd_; ++idim)
     {
       for (int ui=0; ui<nen_; ++ui)
       {
         //const int fui = ui*numdofpernode_ + nsd_;
         for (int vi=0; vi<nen_; ++vi)
         {
           //const int fvi = vi*numdofpernode_ + idim;

           for(int jdim=0;jdim<nsd_;++jdim)
           {
             estif_p_v(vi*nsd_ + idim,ui)
               //ppmat(vi, ui)
               += two_visc_tau*timefac*fac_*(derxy_(jdim, ui))*viscs2_(jdim+(idim*nsd_), vi);
           }
         }
       }
     } // end for(idim)

     if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
     {
       const double two_visc_tauMp = vstabfac*2.0*visc_*tau_Mp;

       for (int idim =0;idim<nsd_;++idim)
       {
         for (int vi=0; vi<nen_; ++vi)
         {
           /* viscous stabilisation */
           for (int jdim=0;jdim<nsd_;++jdim)
           {
             velforce(idim,vi)-= two_visc_tauMp*momres_old_(jdim)*viscs2_(jdim+(idim*nsd_),vi);
           }
         }
       } // end for(idim)
     }
     else
     {
       for (int idim =0;idim<nsd_;++idim)
       {
         for (int vi=0; vi<nen_; ++vi)
         {
           /* viscous stabilisation */
           for (int jdim=0;jdim<nsd_;++jdim)
           {
             velforce(idim,vi)+= fac_*vstabfac*2.0*visc_*sgvelint_(jdim)*viscs2_(jdim+(idim*nsd_),vi);
           }
         }
       } // end for(idim)

     }

   } // end if viscous stabilization on left hand side

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CrossStressStab(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
    const double &                            timefac,
    const double &                            timefacfac,
    const double &                            tau_Mp,
    const double &                            fac3)
{
  /*
                               this part is linearised in
                              combination with the standard
                                  Galerkin term above
                                          +----+    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
                                          |    |
                    /                                \
                   |   /    ~n+af       \   n+af      |
                 + |  | rho*u    o nabla | u     , v  |
                   |   \     (i)        /   (i)       |
                    \                                /
                        |       |
                        +-------+
                     linearisation of
                  this part is performed
                     in the following

   */

     double crossfac;

     if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
     {
       crossfac=densaf_*tau_(1);
     }
     else
     {
       crossfac=densaf_*f3Parameter_->alphaF_*fac3;
     }

     // Stabilization of lhs and the rhs
     if(f3Parameter_->cross_ == INPAR::FLUID::cross_stress_stab)
     {
       /*
              /                         \
             |  /          \   n+af      |
             | | Du o nabla | u     , v  |
             |  \          /             |
              \                         /
       */
       /*
              /                                              \
             |  / / /          \   n+af \         \   n+af    |
             | | | | Du o nabla | u      | o nabla | u   , v  |
             |  \ \ \          /        /         /           |
              \                                              /
       */
       /*
              /                                               \
             |  / / / n+af        \     \         \   n+af     |
             | | | | u     o nabla | Du  | o nabla | u    , v  |
             |  \ \ \             /     /         /            |
              \                                               /
       */
       /*
              /                                             \
             |  / /             /  \ \         \   n+af      |
             | | | nabla o eps | Du | | o nabla | u     , v  |
             |  \ \             \  / /         /             |
              \                                             /
       */
       for(int jdim=0;jdim<nsd_;++jdim)
       {
         for (int ui=0; ui<nen_; ++ui)
         {
           const int fui_p_jdim   = nsd_*ui + jdim;

           for(int idim=0;idim<nsd_;++idim)
           {
             for (int vi=0; vi<nen_; ++vi)
             {
               const int fvi_p_idim = nsd_*vi+idim;

               for(int kdim=0;kdim<nsd_;++kdim)
               {
                 estif_u(fvi_p_idim,fui_p_jdim) -= crossfac*lin_resM_Du(nsd_*kdim+jdim,ui)*vderxy_(idim,kdim)*funct_(vi);
               }
             } // jdim
           } // vi
         } // ui
       } //idim

       /*
                       /                               \
                      |  /                \   n+af      |
                      | | nabla Dp o nabla | u     , v  |
                      |  \                /             |
                       \                               /
       */
       for (int vi=0; vi<nen_; ++vi)
       {
         for (int idim = 0; idim <nsd_; ++idim)
         {
           const int fvi   = nsd_*vi + idim;

           for (int ui=0; ui<nen_; ++ui)
           {
             for(int kdim=0;kdim<nsd_;++kdim)
             {
               estif_p_v(fvi,ui) -= crossfac*timefacfac*vderxy_(idim,kdim)*derxy_(kdim,ui)*funct_(vi);
             }
           }
         }  // end for(idim)
       } // vi
     } // end if(cross_ == INPAR::FLUID::cross_stress_stab)

     // Stabilization only of the rhs
     LINALG::Matrix<nsd_,1> temp;

     temp.Clear();

     for(int jdim=0;jdim<nsd_;++jdim)
     {
       for(int kdim=0;kdim<nsd_;++kdim)
       {
         temp(jdim)+=densaf_*fac_*sgvelint_(kdim)*vderxy_(jdim,kdim);
       }
     }

     for (int idim = 0; idim <nsd_; ++idim)
     {
       for (int vi=0; vi<nen_; ++vi)
       {
         velforce(idim,vi) -= temp(idim)*funct_(vi);
       }
     }  // end for(idim)


  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::ReynoldsStressStab(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_*nsd_,nen_> &          lin_resM_Du,
    const double &                            timefacfac,
    const double &                            fac3)
{
  /*
                            linearisation of
                         this part is performed
                            in the following
                                +--------+
                                |        |
                   /                                 \
                  |  ~n+af     /    ~n+af       \     |
                - |  u     ,  | rho*u    o nabla | v  |
                  |   (i)      \     (i)        /     |
                   \                                 /
                     |   |
                     +---+
            this part is linearised
          in combination with the SUPG
                  term above

  */

    double reyfac;

    if(f3Parameter_->tds_==INPAR::FLUID::subscales_quasistatic)
    {
      reyfac=densaf_*tau_(1);
    }
    else
    {
      reyfac=densaf_*f3Parameter_->alphaF_*fac3;
    }

    if(f3Parameter_->reynolds_ == INPAR::FLUID::reynolds_stress_stab)
    {
      /*
          /                          \
         |  ~n+af                     |
         |  u     , ( Du o nabla ) v  |
         |                            |
          \                          /
      */
      /*
          /                                                 \
         |  ~n+af    / / / n+af        \     \         \     |
         |  u     , | | | u     o nabla | Du  | o nabla | v  |
         |           \ \ \             /     /         /     |
          \                                                 /
      */
      /*
          /                                                 \
         |  ~n+af    / / /          \   n+af \         \     |
         |  u     , | | | Du o nabla | u      | o nabla | v  |
         |           \ \ \          /        /         /     |
          \                                                 /
      */
      /*
          /                                               \
         |  ~n+af    / /             /  \  \         \     |
         |  u     , | | nabla o eps | Du |  | o nabla | v  |
         |           \ \             \  /  /         /     |
          \                                               /
      */
      for(int jdim=0;jdim<nsd_;++jdim)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui_p_jdim   = nsd_*ui + jdim;

          for(int idim=0;idim<nsd_;++idim)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi_p_idim = nsd_*vi+idim;

              for(int kdim=0;kdim<nsd_;++kdim)
              {
                estif_u(fvi_p_idim,fui_p_jdim) += reyfac*lin_resM_Du(nsd_*kdim+jdim,ui)*sgvelint_(idim)*derxy_(idim,vi);
              }
            } // jdim
          } // vi
        } // ui
      } //idim

      /*
          /                                \
         |  ~n+af    /                \     |
         |  u     , | nabla Dp o nabla | v  |
         |           \                /     |
          \                                /
      */
      for (int vi=0; vi<nen_; ++vi)
      {
        for (int idim = 0; idim <nsd_; ++idim)
        {
          const int fvi   = nsd_*vi + idim;

          for (int ui=0; ui<nen_; ++ui)
          {
            for(int kdim=0;kdim<nsd_;++kdim)
            {
              estif_p_v(fvi,ui) += reyfac*timefacfac*sgvelint_(idim)*derxy_(kdim,ui)*derxy_(kdim,vi);
            }
          }
        }  // end for(idim)
      } // vi
    } // end if(reynolds_ == INPAR::FLUID::reynolds_stress_stab)


#if 0
    //!!!!!!!!!!!!!!!!!!!!
    // the rhs contribution is already contained in the supg term!!!
    //!!!!!!!!!!!!!!!!!!!!



    LINALG::Matrix<nsd_,nsd_> temp;

    temp.Clear();

    for(int idim=0;idim<nsd_;++idim)
    {
      for(int kdim=0;kdim<nsd_;++kdim)
      {
        temp(idim,kdim)+=densaf_*fac_*sgvelint_(idim)*sgvelint_(kdim);
      }
    }

    for (int idim = 0; idim <nsd_; ++idim)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for(int kdim=0;kdim<nsd_;++kdim)
        {
          // reynolds stabilisation
          velforce(idim,vi) += temp(idim,kdim)*derxy_(kdim,vi);
        }
      }
    }  // end for(idim)
#endif

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::FineScaleSubGridViscosityTerm(
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          fssgviscfac)
{
  if (nsd_ == 2)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      /* fine-scale subgrid-viscosity term on right hand side */
      /*
                          /                          \
                         |       /    \         / \   |
         - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                         |       \    /         \ /   |
                          \                          /
      */
      velforce(0, vi) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                     +    derxy_(1, vi)*fsvderxy_(0, 1)
                                     +    derxy_(1, vi)*fsvderxy_(1, 0)) ;
      velforce(1, vi) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                     +    derxy_(0, vi)*fsvderxy_(1, 0)
                                     +2.0*derxy_(1, vi)*fsvderxy_(1, 1)) ;
    }
  }
  else if(nsd_ == 3)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      /* fine-scale subgrid-viscosity term on right hand side */
      /*
                            /                          \
                           |       /    \         / \   |
           - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                           |       \    /         \ /   |
                            \                          /
      */
      velforce(0, vi) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                     +    derxy_(1, vi)*fsvderxy_(0, 1)
                                     +    derxy_(1, vi)*fsvderxy_(1, 0)
                                     +    derxy_(2, vi)*fsvderxy_(0, 2)
                                     +    derxy_(2, vi)*fsvderxy_(2, 0)) ;
      velforce(1, vi) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                     +    derxy_(0, vi)*fsvderxy_(1, 0)
                                     +2.0*derxy_(1, vi)*fsvderxy_(1, 1)
                                     +    derxy_(2, vi)*fsvderxy_(1, 2)
                                     +    derxy_(2, vi)*fsvderxy_(2, 1)) ;
      velforce(2, vi) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 2)
                                     +    derxy_(0, vi)*fsvderxy_(2, 0)
                                     +    derxy_(1, vi)*fsvderxy_(1, 2)
                                     +    derxy_(1, vi)*fsvderxy_(2, 1)
                                     +2.0*derxy_(2, vi)*fsvderxy_(2, 2)) ;
    }
  }
  else
    dserror("FineScaleSubGridViscosity is not implemented for 1D");

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::LinMeshMotion_2D(
    LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
    const LINALG::Matrix<nsd_,nen_>&              evelaf,
    const double &                                press,
    const double &                                timefac,
    const double &                                timefacfac)
{
  // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
  // xGderiv_ == xjm_

  // mass + rhs
  for (int vi=0; vi<nen_; ++vi)
  {
    const int tvi   = 3*vi;
    const int tvip  = tvi + 1;

    const double v = fac_*funct_(vi);
    for (int ui=0; ui<nen_; ++ui)
    {
      const int tui   = 3*ui;
      const int tuip  = tui + 1;

      emesh(tvi,   tui ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
      emesh(tvi,   tuip) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);

      emesh(tvip,  tui ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
      emesh(tvip,  tuip) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);
    }
  }

  vderiv_.MultiplyNT(evelaf, deriv_);

//#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

//#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
//#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
//#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
//#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))

  for (int vi=0; vi<nen_; ++vi)
  {
    const int tvi  = 3*vi;
    const int tvip = tvi+1;
    const double v = timefacfac/det_*funct_(vi);
    for (int ui=0; ui<nen_; ++ui)
    {
      const int tui  = 3*ui;
      const int tuip = tui+1;

      emesh(tvi , tui ) += v*(
      + convvelint_(1)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
      );

      emesh(tvi , tuip) += v*(
      + convvelint_(0)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
      );

      emesh(tvip, tui ) += v*(
      + convvelint_(1)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
      );

      emesh(tvip, tuip) += v*(
      + convvelint_(0)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
      );
    }
  }

  // pressure
  for (int vi=0; vi<nen_; ++vi)
  {
    const int tvi  = 3*vi;
    const int tvip = tvi+1;
    const double v = press*timefacfac/det_;
    for (int ui=0; ui<nen_; ++ui)
    {
      const int tui = 3*ui;
      emesh(tvi,  tui + 1) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
      emesh(tvip, tui    ) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
    }
  }

  // div u
  for (int vi=0; vi<nen_; ++vi)
  {
    const int tvipp = 3*vi + 2;
    const double v = timefacfac/det_*funct_(vi);
    for (int ui=0; ui<nen_; ++ui)
    {
      const int tui = 3*ui;
      emesh(tvipp, tui) += v*(
      deriv_(0,ui)*vderiv_(1,1) - deriv_(1,ui)*vderiv_(1,0)
      ) ;

      emesh(tvipp, tui + 1) += v*(
      deriv_(0,ui)*vderiv_(0,1) - deriv_(1,ui)*vderiv_(0,0)
      ) ;
    }
  }


  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::LinMeshMotion_3D(
    LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
    const LINALG::Matrix<nsd_,nen_>&              evelaf,
    const double &                                press,
    const double &                                timefac,
    const double &                                timefacfac)
{
  // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
  // xGderiv_ == xjm_

  // mass + rhs
  for (int vi=0; vi<nen_; ++vi)
  {
    double v = fac_*funct_(vi,0);
    for (int ui=0; ui<nen_; ++ui)
    {
      emesh(vi*4    , ui*4    ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
      emesh(vi*4    , ui*4 + 1) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);
      emesh(vi*4    , ui*4 + 2) += v*(velint_(0)-rhsmom_(0))*derxy_(2, ui);

      emesh(vi*4 + 1, ui*4    ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
      emesh(vi*4 + 1, ui*4 + 1) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);
      emesh(vi*4 + 1, ui*4 + 2) += v*(velint_(1)-rhsmom_(1))*derxy_(2, ui);

      emesh(vi*4 + 2, ui*4    ) += v*(velint_(2)-rhsmom_(2))*derxy_(0, ui);
      emesh(vi*4 + 2, ui*4 + 1) += v*(velint_(2)-rhsmom_(2))*derxy_(1, ui);
      emesh(vi*4 + 2, ui*4 + 2) += v*(velint_(2)-rhsmom_(2))*derxy_(2, ui);
    }
  }

  //vderiv_  = sum(evelaf(i,k) * deriv_(j,k), k);
  vderiv_.MultiplyNT(evelaf,deriv_);

#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
#define derxjm_002(ui) (deriv_(1, ui)*xjm_(2, 1) - deriv_(2, ui)*xjm_(1, 1))

#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
#define derxjm_102(ui) (deriv_(2, ui)*xjm_(1, 0) - deriv_(1, ui)*xjm_(2, 0))

#define derxjm_200(ui) (deriv_(2, ui)*xjm_(1, 1) - deriv_(1, ui)*xjm_(2, 1))
#define derxjm_201(ui) (deriv_(1, ui)*xjm_(2, 0) - deriv_(2, ui)*xjm_(1, 0))

#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
#define derxjm_012(ui) (deriv_(2, ui)*xjm_(0, 1) - deriv_(0, ui)*xjm_(2, 1))

#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))
#define derxjm_112(ui) (deriv_(0, ui)*xjm_(2, 0) - deriv_(2, ui)*xjm_(0, 0))

#define derxjm_210(ui) (deriv_(0, ui)*xjm_(2, 1) - deriv_(2, ui)*xjm_(0, 1))
#define derxjm_211(ui) (deriv_(2, ui)*xjm_(0, 0) - deriv_(0, ui)*xjm_(2, 0))

#define derxjm_021(ui) (deriv_(1, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(1, 2))
#define derxjm_022(ui) (deriv_(0, ui)*xjm_(1, 1) - deriv_(1, ui)*xjm_(0, 1))

#define derxjm_120(ui) (deriv_(0, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(0, 2))
#define derxjm_122(ui) (deriv_(1, ui)*xjm_(0, 0) - deriv_(0, ui)*xjm_(1, 0))

#define derxjm_220(ui) (deriv_(1, ui)*xjm_(0, 1) - deriv_(0, ui)*xjm_(1, 1))
#define derxjm_221(ui) (deriv_(0, ui)*xjm_(1, 0) - deriv_(1, ui)*xjm_(0, 0))

  for (int ui=0; ui<nen_; ++ui)
  {
    double v00 = + convvelint_(1)*(vderiv_(0, 0)*derxjm_(0,0,1,ui) + vderiv_(0, 1)*derxjm_(0,1,1,ui) + vderiv_(0, 2)*derxjm_(0,2,1,ui))
                 + convvelint_(2)*(vderiv_(0, 0)*derxjm_(0,0,2,ui) + vderiv_(0, 1)*derxjm_(0,1,2,ui) + vderiv_(0, 2)*derxjm_(0,2,2,ui));
    double v01 = + convvelint_(0)*(vderiv_(0, 0)*derxjm_(1,0,0,ui) + vderiv_(0, 1)*derxjm_(1,1,0,ui) + vderiv_(0, 2)*derxjm_(1,2,0,ui))
                 + convvelint_(2)*(vderiv_(0, 0)*derxjm_(1,0,2,ui) + vderiv_(0, 1)*derxjm_(1,1,2,ui) + vderiv_(0, 2)*derxjm_(1,2,2,ui));
    double v02 = + convvelint_(0)*(vderiv_(0, 0)*derxjm_(2,0,0,ui) + vderiv_(0, 1)*derxjm_(2,1,0,ui) + vderiv_(0, 2)*derxjm_(2,2,0,ui))
                 + convvelint_(1)*(vderiv_(0, 0)*derxjm_(2,0,1,ui) + vderiv_(0, 1)*derxjm_(2,1,1,ui) + vderiv_(0, 2)*derxjm_(2,2,1,ui));
    double v10 = + convvelint_(1)*(vderiv_(1, 0)*derxjm_(0,0,1,ui) + vderiv_(1, 1)*derxjm_(0,1,1,ui) + vderiv_(1, 2)*derxjm_(0,2,1,ui))
                 + convvelint_(2)*(vderiv_(1, 0)*derxjm_(0,0,2,ui) + vderiv_(1, 1)*derxjm_(0,1,2,ui) + vderiv_(1, 2)*derxjm_(0,2,2,ui));
    double v11 = + convvelint_(0)*(vderiv_(1, 0)*derxjm_(1,0,0,ui) + vderiv_(1, 1)*derxjm_(1,1,0,ui) + vderiv_(1, 2)*derxjm_(1,2,0,ui))
                 + convvelint_(2)*(vderiv_(1, 0)*derxjm_(1,0,2,ui) + vderiv_(1, 1)*derxjm_(1,1,2,ui) + vderiv_(1, 2)*derxjm_(1,2,2,ui));
    double v12 = + convvelint_(0)*(vderiv_(1, 0)*derxjm_(2,0,0,ui) + vderiv_(1, 1)*derxjm_(2,1,0,ui) + vderiv_(1, 2)*derxjm_(2,2,0,ui))
                 + convvelint_(1)*(vderiv_(1, 0)*derxjm_(2,0,1,ui) + vderiv_(1, 1)*derxjm_(2,1,1,ui) + vderiv_(1, 2)*derxjm_(2,2,1,ui));
    double v20 = + convvelint_(1)*(vderiv_(2, 0)*derxjm_(0,0,1,ui) + vderiv_(2, 1)*derxjm_(0,1,1,ui) + vderiv_(2, 2)*derxjm_(0,2,1,ui))
                 + convvelint_(2)*(vderiv_(2, 0)*derxjm_(0,0,2,ui) + vderiv_(2, 1)*derxjm_(0,1,2,ui) + vderiv_(2, 2)*derxjm_(0,2,2,ui));
    double v21 = + convvelint_(0)*(vderiv_(2, 0)*derxjm_(1,0,0,ui) + vderiv_(2, 1)*derxjm_(1,1,0,ui) + vderiv_(2, 2)*derxjm_(1,2,0,ui))
                 + convvelint_(2)*(vderiv_(2, 0)*derxjm_(1,0,2,ui) + vderiv_(2, 1)*derxjm_(1,1,2,ui) + vderiv_(2, 2)*derxjm_(1,2,2,ui));
    double v22 = + convvelint_(0)*(vderiv_(2, 0)*derxjm_(2,0,0,ui) + vderiv_(2, 1)*derxjm_(2,1,0,ui) + vderiv_(2, 2)*derxjm_(2,2,0,ui))
                 + convvelint_(1)*(vderiv_(2, 0)*derxjm_(2,0,1,ui) + vderiv_(2, 1)*derxjm_(2,1,1,ui) + vderiv_(2, 2)*derxjm_(2,2,1,ui));

    for (int vi=0; vi<nen_; ++vi)
    {
      double v = timefacfac/det_*funct_(vi);

      emesh(vi*4 + 0, ui*4 + 0) += v*v00;
      emesh(vi*4 + 0, ui*4 + 1) += v*v01;
      emesh(vi*4 + 0, ui*4 + 2) += v*v02;

      emesh(vi*4 + 1, ui*4 + 0) += v*v10;
      emesh(vi*4 + 1, ui*4 + 1) += v*v11;
      emesh(vi*4 + 1, ui*4 + 2) += v*v12;

      emesh(vi*4 + 2, ui*4 + 0) += v*v20;
      emesh(vi*4 + 2, ui*4 + 1) += v*v21;
      emesh(vi*4 + 2, ui*4 + 2) += v*v22;
    }
  }

  // viscosity

#define xji_00 xji_(0,0)
#define xji_01 xji_(0,1)
#define xji_02 xji_(0,2)
#define xji_10 xji_(1,0)
#define xji_11 xji_(1,1)
#define xji_12 xji_(1,2)
#define xji_20 xji_(2,0)
#define xji_21 xji_(2,1)
#define xji_22 xji_(2,2)

#define xjm(i,j) xjm_(i,j)

  // part 1: derivative of 1/det

  double v = visceff_*timefac*fac_;
  for (int ui=0; ui<nen_; ++ui)
  {
    double derinvJ0 = -v*(deriv_(0,ui)*xji_00 + deriv_(1,ui)*xji_01 + deriv_(2,ui)*xji_02);
    double derinvJ1 = -v*(deriv_(0,ui)*xji_10 + deriv_(1,ui)*xji_11 + deriv_(2,ui)*xji_12);
    double derinvJ2 = -v*(deriv_(0,ui)*xji_20 + deriv_(1,ui)*xji_21 + deriv_(2,ui)*xji_22);
    for (int vi=0; vi<nen_; ++vi)
    {
      double visres0 =   2.0*derxy_(0, vi)* vderxy_(0, 0)
                         +     derxy_(1, vi)*(vderxy_(0, 1) + vderxy_(1, 0))
                         +     derxy_(2, vi)*(vderxy_(0, 2) + vderxy_(2, 0)) ;
      double visres1 =         derxy_(0, vi)*(vderxy_(0, 1) + vderxy_(1, 0))
                         + 2.0*derxy_(1, vi)* vderxy_(1, 1)
                         +     derxy_(2, vi)*(vderxy_(1, 2) + vderxy_(2, 1)) ;
      double visres2 =         derxy_(0, vi)*(vderxy_(0, 2) + vderxy_(2, 0))
                         +     derxy_(1, vi)*(vderxy_(1, 2) + vderxy_(2, 1))
                         + 2.0*derxy_(2, vi)* vderxy_(2, 2) ;
      emesh(vi*4 + 0, ui*4 + 0) += derinvJ0*visres0;
      emesh(vi*4 + 1, ui*4 + 0) += derinvJ0*visres1;
      emesh(vi*4 + 2, ui*4 + 0) += derinvJ0*visres2;

      emesh(vi*4 + 0, ui*4 + 1) += derinvJ1*visres0;
      emesh(vi*4 + 1, ui*4 + 1) += derinvJ1*visres1;
      emesh(vi*4 + 2, ui*4 + 1) += derinvJ1*visres2;

      emesh(vi*4 + 0, ui*4 + 2) += derinvJ2*visres0;
      emesh(vi*4 + 1, ui*4 + 2) += derinvJ2*visres1;
      emesh(vi*4 + 2, ui*4 + 2) += derinvJ2*visres2;
    }
  }

  // part 2: derivative of viscosity residual

  v = timefacfac*visceff_/det_;
  for (int ui=0; ui<nen_; ++ui)
  {
    double v0 = - vderiv_(0,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
                - vderiv_(0,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
                - vderiv_(0,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
                - vderiv_(1,0)*(derxjm_100(ui)*xji_00)
                - vderiv_(1,1)*(derxjm_100(ui)*xji_01)
                - vderiv_(1,2)*(derxjm_100(ui)*xji_02)
                - vderiv_(2,0)*(derxjm_200(ui)*xji_00)
                - vderiv_(2,1)*(derxjm_200(ui)*xji_01)
                - vderiv_(2,2)*(derxjm_200(ui)*xji_02);
    double v1 = - vderiv_(0,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
                - vderiv_(0,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
                - vderiv_(0,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
                - vderiv_(1,0)*(derxjm_110(ui)*xji_00)
                - vderiv_(1,1)*(derxjm_110(ui)*xji_01)
                - vderiv_(1,2)*(derxjm_110(ui)*xji_02)
                - vderiv_(2,0)*(derxjm_210(ui)*xji_00)
                - vderiv_(2,1)*(derxjm_210(ui)*xji_01)
                - vderiv_(2,2)*(derxjm_210(ui)*xji_02);
    double v2 = - vderiv_(0,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
                - vderiv_(0,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
                - vderiv_(0,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
                - vderiv_(1,0)*(derxjm_120(ui)*xji_00)
                - vderiv_(1,1)*(derxjm_120(ui)*xji_01)
                - vderiv_(1,2)*(derxjm_120(ui)*xji_02)
                - vderiv_(2,0)*(derxjm_220(ui)*xji_00)
                - vderiv_(2,1)*(derxjm_220(ui)*xji_01)
                - vderiv_(2,2)*(derxjm_220(ui)*xji_02);

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 0, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(2*derxjm_001(ui)*xji_00 + 2*derxjm_001(ui)*xji_00 + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
         - vderiv_(0,1)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
         - vderiv_(0,2)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
         - vderiv_(1,0)*(derxjm_001(ui)*xji_10)
         - vderiv_(1,1)*(derxjm_011(ui)*xji_10)
         - vderiv_(1,2)*(derxjm_021(ui)*xji_10)
         - vderiv_(2,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
         - vderiv_(2,1)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
         - vderiv_(2,2)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20);
    v1 = - vderiv_(0,0)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
         - vderiv_(0,1)*(2*derxjm_011(ui)*xji_01 + 2*derxjm_011(ui)*xji_01 + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
         - vderiv_(0,2)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
         - vderiv_(1,0)*(derxjm_001(ui)*xji_11)
         - vderiv_(1,1)*(derxjm_011(ui)*xji_11)
         - vderiv_(1,2)*(derxjm_021(ui)*xji_11)
         - vderiv_(2,0)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
         - vderiv_(2,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
         - vderiv_(2,2)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21);
    v2 = - vderiv_(0,0)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
         - vderiv_(0,1)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
         - vderiv_(0,2)*(2*derxjm_021(ui)*xji_02 + 2*derxjm_021(ui)*xji_02 + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
         - vderiv_(1,0)*(derxjm_001(ui)*xji_12)
         - vderiv_(1,1)*(derxjm_011(ui)*xji_12)
         - vderiv_(1,2)*(derxjm_021(ui)*xji_12)
         - vderiv_(2,0)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
         - vderiv_(2,1)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
         - vderiv_(2,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22);

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 0, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(2*derxjm_002(ui)*xji_00 + 2*derxjm_002(ui)*xji_00 + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
         - vderiv_(0,1)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
         - vderiv_(0,2)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
         - vderiv_(1,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
         - vderiv_(1,1)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
         - vderiv_(1,2)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
         - vderiv_(2,0)*(derxjm_002(ui)*xji_20)
         - vderiv_(2,1)*(derxjm_012(ui)*xji_20)
         - vderiv_(2,2)*(derxjm_022(ui)*xji_20);
    v1 = - vderiv_(0,0)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
         - vderiv_(0,1)*(2*derxjm_012(ui)*xji_01 + 2*derxjm_012(ui)*xji_01 + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
         - vderiv_(0,2)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
         - vderiv_(1,0)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
         - vderiv_(1,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
         - vderiv_(1,2)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
         - vderiv_(2,0)*(derxjm_002(ui)*xji_21)
         - vderiv_(2,1)*(derxjm_012(ui)*xji_21)
         - vderiv_(2,2)*(derxjm_022(ui)*xji_21);
    v2 = - vderiv_(0,0)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
         - vderiv_(0,1)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
         - vderiv_(0,2)*(2*derxjm_022(ui)*xji_02 + 2*derxjm_022(ui)*xji_02 + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui))
         - vderiv_(1,0)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
         - vderiv_(1,1)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
         - vderiv_(1,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
         - vderiv_(2,0)*(derxjm_002(ui)*xji_22)
         - vderiv_(2,1)*(derxjm_012(ui)*xji_22)
         - vderiv_(2,2)*(derxjm_022(ui)*xji_22);

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 0, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(derxjm_100(ui)*xji_00)
         - vderiv_(0,1)*(derxjm_110(ui)*xji_00)
         - vderiv_(0,2)*(derxjm_120(ui)*xji_00)
         - vderiv_(1,0)*(2*xji_10*derxjm_100(ui) + 2*xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
         - vderiv_(1,1)*(2*xji_11*derxjm_100(ui) + 2*xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
         - vderiv_(1,2)*(2*xji_12*derxjm_100(ui) + 2*xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
         - vderiv_(2,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
         - vderiv_(2,1)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
         - vderiv_(2,2)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20);
    v1 = - vderiv_(0,0)*(derxjm_100(ui)*xji_01)
         - vderiv_(0,1)*(derxjm_110(ui)*xji_01)
         - vderiv_(0,2)*(derxjm_120(ui)*xji_01)
         - vderiv_(1,0)*(2*xji_10*derxjm_110(ui) + 2*xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
         - vderiv_(1,1)*(2*xji_11*derxjm_110(ui) + 2*xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
         - vderiv_(1,2)*(2*xji_12*derxjm_110(ui) + 2*xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
         - vderiv_(2,0)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
         - vderiv_(2,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
         - vderiv_(2,2)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21);
    v2 = - vderiv_(0,0)*(derxjm_100(ui)*xji_02)
         - vderiv_(0,1)*(derxjm_110(ui)*xji_02)
         - vderiv_(0,2)*(derxjm_120(ui)*xji_02)
         - vderiv_(1,0)*(2*xji_10*derxjm_120(ui) + 2*xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
         - vderiv_(1,1)*(2*xji_11*derxjm_120(ui) + 2*xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
         - vderiv_(1,2)*(2*xji_12*derxjm_120(ui) + 2*xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
         - vderiv_(2,0)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
         - vderiv_(2,1)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
         - vderiv_(2,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22);

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 1, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(derxjm_001(ui)*xji_10)
         - vderiv_(0,1)*(derxjm_001(ui)*xji_11)
         - vderiv_(0,2)*(derxjm_001(ui)*xji_12)
         - vderiv_(1,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
         - vderiv_(1,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
         - vderiv_(1,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
         - vderiv_(2,0)*(derxjm_201(ui)*xji_10)
         - vderiv_(2,1)*(derxjm_201(ui)*xji_11)
         - vderiv_(2,2)*(derxjm_201(ui)*xji_12);
    v1 = - vderiv_(0,0)*(derxjm_011(ui)*xji_10)
         - vderiv_(0,1)*(derxjm_011(ui)*xji_11)
         - vderiv_(0,2)*(derxjm_011(ui)*xji_12)
         - vderiv_(1,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + xji_20*derxjm_211(ui) + xji_21*derxjm_201(ui))
         - vderiv_(1,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
         - vderiv_(1,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + xji_22*derxjm_211(ui) + xji_21*derxjm_221(ui))
         - vderiv_(2,0)*(derxjm_211(ui)*xji_10)
         - vderiv_(2,1)*(derxjm_211(ui)*xji_11)
         - vderiv_(2,2)*(derxjm_211(ui)*xji_12);
    v2 = - vderiv_(0,0)*(derxjm_021(ui)*xji_10)
         - vderiv_(0,1)*(derxjm_021(ui)*xji_11)
         - vderiv_(0,2)*(derxjm_021(ui)*xji_12)
         - vderiv_(1,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + xji_20*derxjm_221(ui) + xji_22*derxjm_201(ui))
         - vderiv_(1,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
         - vderiv_(1,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
         - vderiv_(2,0)*(derxjm_221(ui)*xji_10)
         - vderiv_(2,1)*(derxjm_221(ui)*xji_11)
         - vderiv_(2,2)*(derxjm_221(ui)*xji_12);

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 1, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
         - vderiv_(0,1)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
         - vderiv_(0,2)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
         - vderiv_(1,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + 2*xji_10*derxjm_102(ui) + 2*xji_10*derxjm_102(ui))
         - vderiv_(1,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + 2*xji_11*derxjm_102(ui) + 2*xji_10*derxjm_112(ui))
         - vderiv_(1,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + 2*xji_12*derxjm_102(ui) + 2*xji_10*derxjm_122(ui))
         - vderiv_(2,0)*(derxjm_102(ui)*xji_20)
         - vderiv_(2,1)*(derxjm_112(ui)*xji_20)
         - vderiv_(2,2)*(derxjm_122(ui)*xji_20);
    v1 = - vderiv_(0,0)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
         - vderiv_(0,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
         - vderiv_(0,2)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
         - vderiv_(1,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + 2*xji_10*derxjm_112(ui) + 2*xji_11*derxjm_102(ui))
         - vderiv_(1,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + 2*xji_11*derxjm_112(ui) + 2*xji_11*derxjm_112(ui))
         - vderiv_(1,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + 2*xji_12*derxjm_112(ui) + 2*xji_11*derxjm_122(ui))
         - vderiv_(2,0)*(derxjm_102(ui)*xji_21)
         - vderiv_(2,1)*(derxjm_112(ui)*xji_21)
         - vderiv_(2,2)*(derxjm_122(ui)*xji_21);
    v2 = - vderiv_(0,0)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
         - vderiv_(0,1)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
         - vderiv_(0,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
         - vderiv_(1,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + 2*xji_10*derxjm_122(ui) + 2*xji_12*derxjm_102(ui))
         - vderiv_(1,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + 2*xji_11*derxjm_122(ui) + 2*xji_12*derxjm_112(ui))
         - vderiv_(1,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + 2*xji_12*derxjm_122(ui) + 2*xji_12*derxjm_122(ui))
         - vderiv_(2,0)*(derxjm_102(ui)*xji_22)
         - vderiv_(2,1)*(derxjm_112(ui)*xji_22)
         - vderiv_(2,2)*(derxjm_122(ui)*xji_22);

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 1, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(derxjm_200(ui)*xji_00)
         - vderiv_(0,1)*(derxjm_210(ui)*xji_00)
         - vderiv_(0,2)*(derxjm_220(ui)*xji_00)
         - vderiv_(1,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
         - vderiv_(1,1)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
         - vderiv_(1,2)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
         - vderiv_(2,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + 2*xji_20*derxjm_200(ui) + 2*xji_20*derxjm_200(ui))
         - vderiv_(2,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + 2*xji_21*derxjm_200(ui) + 2*xji_20*derxjm_210(ui))
         - vderiv_(2,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + 2*xji_22*derxjm_200(ui) + 2*xji_20*derxjm_220(ui));
    v1 = - vderiv_(0,0)*(derxjm_200(ui)*xji_01)
         - vderiv_(0,1)*(derxjm_210(ui)*xji_01)
         - vderiv_(0,2)*(derxjm_220(ui)*xji_01)
         - vderiv_(1,0)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
         - vderiv_(1,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
         - vderiv_(1,2)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
         - vderiv_(2,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + 2*xji_20*derxjm_210(ui) + 2*xji_21*derxjm_200(ui))
         - vderiv_(2,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + 2*xji_21*derxjm_210(ui) + 2*xji_21*derxjm_210(ui))
         - vderiv_(2,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + 2*xji_22*derxjm_210(ui) + 2*xji_21*derxjm_220(ui));
    v2 = - vderiv_(0,0)*(derxjm_200(ui)*xji_02)
         - vderiv_(0,1)*(derxjm_210(ui)*xji_02)
         - vderiv_(0,2)*(derxjm_220(ui)*xji_02)
         - vderiv_(1,0)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20)
         - vderiv_(1,1)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21)
         - vderiv_(1,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22)
         - vderiv_(2,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + 2*xji_20*derxjm_220(ui) + 2*xji_22*derxjm_200(ui))
         - vderiv_(2,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + 2*xji_21*derxjm_220(ui) + 2*xji_22*derxjm_210(ui))
         - vderiv_(2,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + 2*xji_22*derxjm_220(ui) + 2*xji_22*derxjm_220(ui));

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 2, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
         - vderiv_(0,1)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
         - vderiv_(0,2)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
         - vderiv_(1,0)*(derxjm_201(ui)*xji_10)
         - vderiv_(1,1)*(derxjm_211(ui)*xji_10)
         - vderiv_(1,2)*(derxjm_221(ui)*xji_10)
         - vderiv_(2,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + 2*xji_20*derxjm_201(ui) + 2*xji_20*derxjm_201(ui))
         - vderiv_(2,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + 2*xji_21*derxjm_201(ui) + 2*xji_20*derxjm_211(ui))
         - vderiv_(2,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + 2*xji_22*derxjm_201(ui) + 2*xji_20*derxjm_221(ui));
    v1 = - vderiv_(0,0)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
         - vderiv_(0,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
         - vderiv_(0,2)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
         - vderiv_(1,0)*(derxjm_201(ui)*xji_11)
         - vderiv_(1,1)*(derxjm_211(ui)*xji_11)
         - vderiv_(1,2)*(derxjm_221(ui)*xji_11)
         - vderiv_(2,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + 2*xji_20*derxjm_211(ui) + 2*xji_21*derxjm_201(ui))
         - vderiv_(2,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + 2*xji_21*derxjm_211(ui) + 2*xji_21*derxjm_211(ui))
         - vderiv_(2,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + 2*xji_22*derxjm_211(ui) + 2*xji_21*derxjm_221(ui));
    v2 = - vderiv_(0,0)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20)
         - vderiv_(0,1)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21)
         - vderiv_(0,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22)
         - vderiv_(1,0)*(derxjm_201(ui)*xji_12)
         - vderiv_(1,1)*(derxjm_211(ui)*xji_12)
         - vderiv_(1,2)*(derxjm_221(ui)*xji_12)
         - vderiv_(2,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + 2*xji_20*derxjm_221(ui) + 2*xji_22*derxjm_201(ui))
         - vderiv_(2,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + 2*xji_21*derxjm_221(ui) + 2*xji_22*derxjm_211(ui))
         - vderiv_(2,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + 2*xji_22*derxjm_221(ui) + 2*xji_22*derxjm_221(ui));

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 2, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }

    ////////////////////////////////////////////////////////////////

    v0 = - vderiv_(0,0)*(derxjm_002(ui)*xji_20)
         - vderiv_(0,1)*(derxjm_002(ui)*xji_21)
         - vderiv_(0,2)*(derxjm_002(ui)*xji_22)
         - vderiv_(1,0)*(derxjm_102(ui)*xji_20)
         - vderiv_(1,1)*(derxjm_102(ui)*xji_21)
         - vderiv_(1,2)*(derxjm_102(ui)*xji_22)
         - vderiv_(2,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
         - vderiv_(2,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
         - vderiv_(2,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui));
    v1 = - vderiv_(0,0)*(derxjm_012(ui)*xji_20)
         - vderiv_(0,1)*(derxjm_012(ui)*xji_21)
         - vderiv_(0,2)*(derxjm_012(ui)*xji_22)
         - vderiv_(1,0)*(derxjm_112(ui)*xji_20)
         - vderiv_(1,1)*(derxjm_112(ui)*xji_21)
         - vderiv_(1,2)*(derxjm_112(ui)*xji_22)
         - vderiv_(2,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + xji_10*derxjm_112(ui) + xji_11*derxjm_102(ui))
         - vderiv_(2,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
         - vderiv_(2,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + xji_12*derxjm_112(ui) + xji_11*derxjm_122(ui));
    v2 = - vderiv_(0,0)*(derxjm_022(ui)*xji_20)
         - vderiv_(0,1)*(derxjm_022(ui)*xji_21)
         - vderiv_(0,2)*(derxjm_022(ui)*xji_22)
         - vderiv_(1,0)*(derxjm_122(ui)*xji_20)
         - vderiv_(1,1)*(derxjm_122(ui)*xji_21)
         - vderiv_(1,2)*(derxjm_122(ui)*xji_22)
         - vderiv_(2,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + xji_10*derxjm_122(ui) + xji_12*derxjm_102(ui))
         - vderiv_(2,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
         - vderiv_(2,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui));

    for (int vi=0; vi<nen_; ++vi)
    {
      emesh(vi*4 + 2, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
    }
  }


  // pressure
  for (int vi=0; vi<nen_; ++vi)
  {
    double v = press*timefacfac/det_;
    for (int ui=0; ui<nen_; ++ui)
    {
      emesh(vi*4    , ui*4 + 1) += v*(deriv_(0, vi)*derxjm_(0,0,1,ui) + deriv_(1, vi)*derxjm_(0,1,1,ui) + deriv_(2, vi)*derxjm_(0,2,1,ui)) ;
      emesh(vi*4    , ui*4 + 2) += v*(deriv_(0, vi)*derxjm_(0,0,2,ui) + deriv_(1, vi)*derxjm_(0,1,2,ui) + deriv_(2, vi)*derxjm_(0,2,2,ui)) ;

      emesh(vi*4 + 1, ui*4 + 0) += v*(deriv_(0, vi)*derxjm_(1,0,0,ui) + deriv_(1, vi)*derxjm_(1,1,0,ui) + deriv_(2, vi)*derxjm_(1,2,0,ui)) ;
      emesh(vi*4 + 1, ui*4 + 2) += v*(deriv_(0, vi)*derxjm_(1,0,2,ui) + deriv_(1, vi)*derxjm_(1,1,2,ui) + deriv_(2, vi)*derxjm_(1,2,2,ui)) ;

      emesh(vi*4 + 2, ui*4 + 0) += v*(deriv_(0, vi)*derxjm_(2,0,0,ui) + deriv_(1, vi)*derxjm_(2,1,0,ui) + deriv_(2, vi)*derxjm_(2,2,0,ui)) ;
      emesh(vi*4 + 2, ui*4 + 1) += v*(deriv_(0, vi)*derxjm_(2,0,1,ui) + deriv_(1, vi)*derxjm_(2,1,1,ui) + deriv_(2, vi)*derxjm_(2,2,1,ui)) ;
    }
  }

  // div u
  for (int vi=0; vi<nen_; ++vi)
  {
    double v = timefacfac/det_*funct_(vi,0);
    for (int ui=0; ui<nen_; ++ui)
    {
      emesh(vi*4 + 3, ui*4 + 0) += v*(
        + vderiv_(1, 0)*derxjm_(0,0,1,ui) + vderiv_(1, 1)*derxjm_(0,1,1,ui) + vderiv_(1, 2)*derxjm_(0,2,1,ui)
        + vderiv_(2, 0)*derxjm_(0,0,2,ui) + vderiv_(2, 1)*derxjm_(0,1,2,ui) + vderiv_(2, 2)*derxjm_(0,2,2,ui)
        ) ;

      emesh(vi*4 + 3, ui*4 + 1) += v*(
        + vderiv_(0, 0)*derxjm_(1,0,0,ui) + vderiv_(0, 1)*derxjm_(1,1,0,ui) + vderiv_(0, 2)*derxjm_(1,2,0,ui)
        + vderiv_(2, 0)*derxjm_(1,0,2,ui) + vderiv_(2, 1)*derxjm_(1,1,2,ui) + vderiv_(2, 2)*derxjm_(1,2,2,ui)
        ) ;

      emesh(vi*4 + 3, ui*4 + 2) += v*(
        + vderiv_(0, 0)*derxjm_(2,0,0,ui) + vderiv_(0, 1)*derxjm_(2,1,0,ui) + vderiv_(0, 2)*derxjm_(2,2,0,ui)
        + vderiv_(1, 0)*derxjm_(2,0,1,ui) + vderiv_(1, 1)*derxjm_(2,1,1,ui) + vderiv_(1, 2)*derxjm_(2,2,1,ui)
        ) ;
    }
  }

  return;
}

#endif
#endif

