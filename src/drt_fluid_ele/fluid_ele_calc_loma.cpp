/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_loma.cpp

\brief Low-Mach-number flow routines for calculation of fluid element

\level 2

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_loma.H"
#include "fluid_ele.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "../drt_lib/drt_element_integration_select.H"

#include "../drt_geometry/position_array.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcLoma<distype> * DRT::ELEMENTS::FluidEleCalcLoma<distype>::Instance( bool create )
{
  static FluidEleCalcLoma<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcLoma<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcLoma<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcLoma<distype>::FluidEleCalcLoma()
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
{
  // we use the standard parameter list here, since there are not any additional
  // loma-specific parameters required in this derived class
  my::fldpara_ =  DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcLoma<distype>::Evaluate(DRT::ELEMENTS::Fluid*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elemat1_epetra,
                                                 Epetra_SerialDenseMatrix&  elemat2_epetra,
                                                 Epetra_SerialDenseVector&  elevec1_epetra,
                                                 Epetra_SerialDenseVector&  elevec2_epetra,
                                                 Epetra_SerialDenseVector&  elevec3_epetra,
                                                 bool                       offdiag)
{
  if (not offdiag)
    return my::Evaluate( ele, discretization, lm, params, mat,
                     elemat1_epetra, elemat2_epetra,
                     elevec1_epetra, elevec2_epetra, elevec3_epetra,
                     my::intpoints_);
  else
    return EvaluateOD( ele, discretization, lm, params, mat,
                         elemat1_epetra, elemat2_epetra,
                         elevec1_epetra, elevec2_epetra, elevec3_epetra,
                         my::intpoints_);
}


/*----------------------------------------------------------------------*
 * evaluation of off-diagonal matrix block for monolithic loma solver (2)
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcLoma<distype>::EvaluateOD(
                                                     DRT::ELEMENTS::Fluid*       ele,
                                                     DRT::Discretization & discretization,
                                                     const std::vector<int> &     lm,
                                                     Teuchos::ParameterList&       params,
                                                     Teuchos::RCP<MAT::Material> & mat,
                                                     Epetra_SerialDenseMatrix&  elemat1_epetra,
                                                     Epetra_SerialDenseMatrix&  elemat2_epetra,
                                                     Epetra_SerialDenseVector&  elevec1_epetra,
                                                     Epetra_SerialDenseVector&  elevec2_epetra,
                                                     Epetra_SerialDenseVector&  elevec3_epetra,
                                                     const DRT::UTILS::GaussIntegration & intpoints )
{
  // rotationally symmetric periodic bc's: do setup for current element
  my::rotsymmpbc_->Setup(ele);

  // construct view
  LINALG::Matrix<(my::nsd_+1)*my::nen_,my::nen_> elemat1(elemat1_epetra,true);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes,
  // with pressure gradient prescribed as body force included for turbulent
  // channel flow and with scatra body force included for variable-density flow
  // (evaluation at time n+alpha_F for generalized-alpha scheme,
  //  and at time n+1 otherwise)
  // ---------------------------------------------------------------------
  LINALG::Matrix<my::nsd_,my::nen_> ebofoaf(true);
  LINALG::Matrix<my::nsd_,my::nen_> eprescpgaf(true);
  LINALG::Matrix<my::nen_,1>    escabofoaf(true);
  my::BodyForce(ele,ebofoaf,eprescpgaf,escabofoaf);

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, scalar,
  // acceleration/scalar time derivative and history
  // velocity/pressure and scalar values are at time n+alpha_F/n+alpha_M
  // for generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration/scalar time derivative values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F and n+alpha_M
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1>    epreaf(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evelaf, &epreaf,"velaf");

  LINALG::Matrix<my::nsd_,my::nen_> evelam(true);
  LINALG::Matrix<my::nen_,1>    epream(true);
  if (my::fldpara_->PhysicalType() == INPAR::FLUID::weakly_compressible && my::fldparatimint_->IsGenalpha())
  {
    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evelam, &epream,"velam");
  }
  if (my::fldpara_->PhysicalType() == INPAR::FLUID::weakly_compressible_stokes && my::fldparatimint_->IsGenalpha())
  {
    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evelam, &epream,"velam");
  }

  LINALG::Matrix<my::nen_,1> escaaf(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, NULL, &escaaf,"scaaf");

  LINALG::Matrix<my::nsd_,my::nen_> emhist(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &emhist, NULL,"hist");

  LINALG::Matrix<my::nsd_,my::nen_> eaccam(true);
  LINALG::Matrix<my::nen_,1>    escadtam(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &eaccam, &escadtam,"accam");

  LINALG::Matrix<my::nsd_,my::nen_> eveln(true);
  LINALG::Matrix<my::nen_,1>    escaam(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &eveln, &escaam,"scaam");

  if (not my::fldparatimint_->IsGenalpha()) eaccam.Clear();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);

  if (ele->IsAle())
  {
    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &edispnp, NULL,"dispnp");
    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &egridv, NULL,"gridv");
  }

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if(my::isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,my::myknots_,my::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
    return(0);
  } // Nurbs specific stuff

  //----------------------------------------------------------------
  // prepare some dynamic Smagorinsky related stuff
  //----------------------------------------------------------------
  double CsDeltaSq = 0.0;
  double CiDeltaSq = 0.0;
  if (my::fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
  {
    Teuchos::RCP<Epetra_Vector> ele_CsDeltaSq = params.sublist("TURBULENCE MODEL").get<Teuchos::RCP<Epetra_Vector> >("col_Cs_delta_sq");
    Teuchos::RCP<Epetra_Vector> ele_CiDeltaSq = params.sublist("TURBULENCE MODEL").get<Teuchos::RCP<Epetra_Vector> >("col_Ci_delta_sq");
    const int id = ele->LID();
    CsDeltaSq = (*ele_CsDeltaSq)[id];
    CiDeltaSq = (*ele_CiDeltaSq)[id];
  }

  // set element id
  my::eid_ = ele->Id();
  // call inner evaluate (does not know about DRT element or discretization object)
  int result = EvaluateOD(
    params,
    ebofoaf,
    eprescpgaf,
    elemat1,
    evelaf,
    epreaf,
    epream,
    escaaf,
    emhist,
    eaccam,
    escadtam,
    escabofoaf,
    eveln,
    escaam,
    edispnp,
    egridv,
    mat,
    ele->IsAle(),
    CsDeltaSq,
    CiDeltaSq,
    intpoints);

  return result;
}


/*----------------------------------------------------------------------*
 * evaluation of off-diagonal matrix block for monolithic loma solver (3)
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcLoma<distype>::EvaluateOD(
  Teuchos::ParameterList&              params,
  const LINALG::Matrix<my::nsd_,my::nen_> &    ebofoaf,
  const LINALG::Matrix<my::nsd_,my::nen_> &    eprescpgaf,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,my::nen_> & elemat1,
  const LINALG::Matrix<my::nsd_,my::nen_> &    evelaf,
  const LINALG::Matrix<my::nen_,1>    &    epreaf,
  const LINALG::Matrix<my::nen_,1>    &    epream,
  const LINALG::Matrix<my::nen_,1>    &    escaaf,
  const LINALG::Matrix<my::nsd_,my::nen_> &    emhist,
  const LINALG::Matrix<my::nsd_,my::nen_> &    eaccam,
  const LINALG::Matrix<my::nen_,1>    &    escadtam,
  const LINALG::Matrix<my::nen_,1>    &    escabofoaf,
  const LINALG::Matrix<my::nsd_,my::nen_> &    eveln,
  const LINALG::Matrix<my::nen_,1>    &    escaam,
  const LINALG::Matrix<my::nsd_,my::nen_> &    edispnp,
  const LINALG::Matrix<my::nsd_,my::nen_> &    egridv,
  Teuchos::RCP<MAT::Material>          mat,
  bool                                 isale,
  double                               CsDeltaSq,
  double                               CiDeltaSq,
  const DRT::UTILS::GaussIntegration & intpoints)
{
  // flag for higher order elements
  my::is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (my::fldpara_->IsInconsistent() == true) my::is_higher_order_ele_ = false;

  // stationary formulation does not support ALE formulation
  if (isale and my::fldparatimint_->IsStationary())
    dserror("No ALE support within stationary fluid solver.");

  // set thermodynamic pressure at n+1/n+alpha_F and n+alpha_M/n and
  // its time derivative at n+alpha_M/n+1
  const double thermpressaf   = params.get<double>("thermpress at n+alpha_F/n+1");
  const double thermpressam   = params.get<double>("thermpress at n+alpha_M/n");
  const double thermpressdtaf = params.get<double>("thermpressderiv at n+alpha_F/n+1");
  const double thermpressdtam = params.get<double>("thermpressderiv at n+alpha_M/n+1");

  // ---------------------------------------------------------------------
  // set parameters for classical turbulence models
  // ---------------------------------------------------------------------
  Teuchos::ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");

  double Cs_delta_sq   = 0.0;
  double Ci_delta_sq   = 0.0;
  my::visceff_  = 0.0;

  // remember the layer of averaging for the dynamic Smagorinsky model
  int  nlayer=0;

  my::GetTurbulenceParams(turbmodelparams,
                      Cs_delta_sq,
                      Ci_delta_sq,
                      nlayer,
                      CsDeltaSq,
                      CiDeltaSq);

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix
  // ---------------------------------------------------------------------
  SysmatOD(
                        ebofoaf,
                        eprescpgaf,
                        evelaf,
                        eveln,
                        epreaf,
                        epream,
                        eaccam,
                        escaaf,
                        escaam,
                        escadtam,
                        escabofoaf,
                        emhist,
                        edispnp,
                        egridv,
                        elemat1,
                        thermpressaf,
                        thermpressam,
                        thermpressdtaf,
                        thermpressdtam,
                        mat,
                        Cs_delta_sq,
                        Ci_delta_sq,
                        isale,
                        intpoints);

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix for off-diagonal matrix block              |
 |  for monolithic low-Mach-number solver                      vg 10/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcLoma<distype>::SysmatOD(
  const LINALG::Matrix<my::nsd_,my::nen_>&     ebofoaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&     eprescpgaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&     evelaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&     eveln,
  const LINALG::Matrix<my::nen_,1>&        epreaf,
  const LINALG::Matrix<my::nen_,1>&        epream,
  const LINALG::Matrix<my::nsd_,my::nen_>&     eaccam,
  const LINALG::Matrix<my::nen_,1>&        escaaf,
  const LINALG::Matrix<my::nen_,1>&        escaam,
  const LINALG::Matrix<my::nen_,1>&        escadtam,
  const LINALG::Matrix<my::nen_,1>&        escabofoaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&     emhist,
  const LINALG::Matrix<my::nsd_,my::nen_>&     edispnp,
  const LINALG::Matrix<my::nsd_,my::nen_>&     egridv,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,my::nen_>&  estif,
  const double                         thermpressaf,
  const double                         thermpressam,
  const double                         thermpressdtaf,
  const double                         thermpressdtam,
  Teuchos::RCP<const MAT::Material>    material,
  double&                              Cs_delta_sq,
  double&                              Ci_delta_sq,
  bool                                 isale,
  const DRT::UTILS::GaussIntegration & intpoints
  )
{
  // definition of temperature-based residual vector for continuity
  // and energy-conservation equation
  LINALG::Matrix<my::nen_,1> lin_resC_DT(true);
  LINALG::Matrix<my::nen_,1> lin_resE_DT(true);

  // add displacement when fluid nodes move in the ALE case
  if (isale) my::xyze_ += edispnp;

  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter();

  // set element area or volume
  const double vol = my::fac_;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // get material parameters at element center
  if (not my::fldpara_->MatGp() or not my::fldpara_->TauGp())
  {
    my::GetMaterialParams(material,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

    // calculate all-scale subgrid viscosity at element center
    my::visceff_ = my::visc_;
    if (my::fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
        my::fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
        my::fldpara_->TurbModAction() == INPAR::FLUID::vreman)
    {
      my::CalcSubgrVisc(evelaf,vol,Cs_delta_sq,Ci_delta_sq);
      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      my::visceff_ += my::sgvisc_;
    }
  }

  // calculate stabilization parameter at element center
  if (not my::fldpara_->TauGp())
  {
    // get convective velocity at element center for evaluation of
    // stabilization parameter
    my::velint_.Multiply(evelaf,my::funct_);
    my::convvelint_.Update(my::velint_);
    if (isale) my::convvelint_.Multiply(-1.0,egridv,my::funct_,1.0);

    // calculate stabilization parameters at element center
    my::CalcStabParameter(vol);
  }

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // get convective velocity at integration point
    // (including grid velocity in ALE case,
    // values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::convvelint_.Multiply(evelaf,my::funct_);
    if (isale)
    {
      my::gridvelint_.Multiply(egridv,my::funct_);
      my::convvelint_.Update(-1.0,my::gridvelint_,1.0);
    }

    //----------------------------------------------------------------------
    // potential evaluation of material parameters, subgrid viscosity
    // and/or stabilization parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    if (my::fldpara_->MatGp())
    {
      my::GetMaterialParams(material,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);
      // calculate all-scale or fine-scale subgrid viscosity at integration point
      my::visceff_ = my::visc_;
      if (my::fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
          my::fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
          my::fldpara_->TurbModAction() == INPAR::FLUID::vreman)
      {
        my::CalcSubgrVisc(evelaf,vol,Cs_delta_sq,Ci_delta_sq);
        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        my::visceff_ += my::sgvisc_;
      }
    }

    // calculate stabilization parameter at integration point
    if (my::fldpara_->TauGp())
      my::CalcStabParameter(vol);

    // evaluation of convective operator
    my::conv_c_.MultiplyTN(my::derxy_,my::convvelint_);

    // compute additional Galerkin terms on right-hand side of continuity equation
    // -> different for generalized-alpha and other time-integration schemes
    my::ComputeGalRHSContEq(eveln,escaaf,escaam,escadtam,isale);

    // compute subgrid-scale part of scalar
    // -> different for generalized-alpha and other time-integration schemes
    my::ComputeSubgridScaleScalar(escaaf,escaam);

    // update material parameters including subgrid-scale part of scalar
    if (my::fldpara_->UpdateMat())
    {
      if (my::fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or my::fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
          my::fldpara_->TurbModAction() == INPAR::FLUID::vreman)
        dserror("No material update in combination with smagorinsky model!");
      my::UpdateMaterialParams(material,evelaf,epreaf,epream,escaaf,escaam,thermpressaf,thermpressam,my::sgscaint_);
      my::visceff_ = my::visc_;
      if (my::fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or my::fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
          my::fldpara_->TurbModAction() == INPAR::FLUID::vreman)
        my::visceff_ += my::sgvisc_;
    }
    //----------------------------------------------------------------------
    //  evaluate temperature-based residual vector for continuity equation
    //----------------------------------------------------------------------
    // set residual vector for continuity equation to zero
    lin_resC_DT.Clear();

    // transient term
    if (not my::fldparatimint_->IsStationary())
    {
      const double scadtfacfac = my::scadtfac_*my::fac_;
      for (int ui=0; ui<my::nen_; ++ui)
      {
        lin_resC_DT(ui) += scadtfacfac*my::funct_(ui);
      }
    }

    // convective term
    const double timefac_scaconvfacaf = my::fldparatimint_->TimeFac()*my::fac_*my::scaconvfacaf_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      lin_resC_DT(ui) += timefac_scaconvfacaf*my::conv_c_(ui);
    }

    //----------------------------------------------------------------------
    // subgrid-scale-velocity term (governed by cross-stress flag here)
    //----------------------------------------------------------------------
    if (my::fldpara_->ContiCross()    == INPAR::FLUID::cross_stress_stab or
        my::fldpara_->ContiReynolds() == INPAR::FLUID::reynolds_stress_stab)
    {
      //----------------------------------------------------------------------
      //  evaluation of various values at integration point:
      //  1) velocity derivatives
      //  2) pressure (including derivatives)
      //  3) body-force vector
      //  4) "history" vector for momentum equation
      //----------------------------------------------------------------------
      // get velocity derivatives at integration point
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      my::vderxy_.MultiplyNT(evelaf,my::derxy_);

      // get pressure gradient at integration point
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      my::gradp_.Multiply(my::derxy_,epreaf);

      // get bodyforce at integration point
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      my::bodyforce_.Multiply(ebofoaf,my::funct_);
      // get prescribed pressure gradient acting as body force
      // (required for turbulent channel flow)
      my::generalbodyforce_.Multiply(eprescpgaf,my::funct_);

      // get momentum history data at integration point
      // (only required for one-step-theta and BDF2 time-integration schemes)
      my::histmom_.Multiply(emhist,my::funct_);

      // convective term from previous iteration
      my::conv_old_.Multiply(my::vderxy_,my::convvelint_);

      // compute viscous term from previous iteration
      if (my::is_higher_order_ele_) my::CalcDivEps(evelaf,eveln);
      else my::visc_old_.Clear();

      // compute residual of momentum equation and subgrid-scale velocity
      // -> residual of momentum equation different for generalized-alpha
      //    and other time-integration schemes
      // -> no time-dependent subgrid scales considered here
      double fac1    = 0.0;
      double fac2    = 0.0;
      double fac3    = 0.0;
      double facMtau = 0.0;
      double * saccn = NULL;
      double * sveln = NULL;
      double * svelnp = NULL;
      my::ComputeSubgridScaleVelocity(eaccam,fac1,fac2,fac3,facMtau,*iquad,saccn,sveln,svelnp);

      if (my::fldpara_->ContiCross()==INPAR::FLUID::cross_stress_stab)
      {
        // evaluate subgrid-scale-velocity term
        for (int ui=0; ui<my::nen_; ++ui)
        {
          lin_resC_DT(ui) += timefac_scaconvfacaf*my::sgconv_c_(ui);
        }
      }
    }

    //----------------------------------------------------------------------
    // computation of standard Galerkin contributions to element matrix:
    // transient and convective term (potentially incl. cross-stress term)
    //----------------------------------------------------------------------
    /*
            /                                        \
           |         1     / dT     /         \   \   |
       -   |    q , --- * | ---- + | u o nabla | T |  |
           |         T     \ dt     \         /   /   |
            \                                        /
    */
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int numdof_vi_p_nsd = my::numdofpernode_*vi+my::nsd_;
      for (int ui=0; ui<my::nen_; ++ui)
      {
          estif(numdof_vi_p_nsd,ui) -= my::funct_(vi)*lin_resC_DT(ui);
      }
    }

    //----------------------------------------------------------------------
    // computation of SUPG and contributions to element matrix
    // (potentially including Reynolds-stress term)
    //----------------------------------------------------------------------
    if (my::fldpara_->SUPG())
    {
      // weighting functions for SUPG term
      LINALG::Matrix<my::nen_,1> supg_rey_weight;
      const double prefac = my::scaconvfacaf_*my::tau_(0);
      for (int vi=0; vi<my::nen_; ++vi)
      {
        supg_rey_weight(vi) = prefac*my::conv_c_(vi);
      }

      // weighting functions for Reynolds-stress term
      if (my::fldpara_->Reynolds() == INPAR::FLUID::reynolds_stress_stab)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          supg_rey_weight(vi) += prefac*my::sgconv_c_(vi);
        }
      }

      //----------------------------------------------------------------------
      //  evaluate residual vector for energy-conservation equation
      //----------------------------------------------------------------------
      // set residual vector for energy-conservation equation to zero
      lin_resE_DT.Clear();

      // transient term
      if (not my::fldparatimint_->IsStationary())
      {
        const double densamfac = my::fac_*my::densam_;
        for (int ui=0; ui<my::nen_; ++ui)
        {
          lin_resE_DT(ui) += densamfac*my::funct_(ui);
        }
      }

      // convective term
      const double denstimefac = my::fldparatimint_->TimeFac()*my::fac_*my::densaf_;
      for (int ui=0; ui<my::nen_; ++ui)
      {
        lin_resE_DT(ui) += denstimefac*my::conv_c_(ui);
      }

      // diffusive term
      if (my::is_higher_order_ele_)
      {
        // compute second derivatives of shape functions
        LINALG::Matrix<my::nen_,1> diff;
        diff.Clear();
        // compute N,xx + N,yy + N,zz for each shape function
        for (int i=0; i<my::nen_; ++i)
        {
          for (int j = 0; j<my::nsd_; ++j)
          {
            diff(i) += my::derxy2_(j,i);
          }
        }

        const double difftimefac = my::fldparatimint_->TimeFac()*my::fac_*my::diffus_;
        for (int ui=0; ui<my::nen_; ++ui)
        {
          lin_resE_DT(ui) -= difftimefac*diff(ui);
        }
      }

      /*    SUPG/Reynolds-stress term
          /                                                                      \
         |   /         \            dDT          /          \                     |
     -   |  | u o nabla | q , rho * ---- + rho * | u o nabla | DT - diff * lap DT |
         |   \         /             dt          \           /                    |
          \                                                                      /
      */
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int numdof_vi_p_nsd = my::numdofpernode_*vi+my::nsd_;
        for (int ui=0; ui<my::nen_; ++ui)
        {
          estif(numdof_vi_p_nsd,ui) -= supg_rey_weight(vi)*lin_resE_DT(ui);
        }
      }
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  return;
}


// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcLoma<DRT::Element::nurbs27>;



