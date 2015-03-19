/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_calc_poro_p2.cpp

 \brief Internal implementation of poro Fluid element (p2 poro fluid)

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "fluid_ele_calc_poro_p2.H"

#include "fluid_ele.H"
#include "fluid_ele_parameter_poro.H"
#include "../drt_lib/drt_element_integration_select.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"
#include "fluid_ele_action.H"

#include "../drt_geometry/position_array.H"


/*----------------------------------------------------------------------*
 *  create/delete instance (public)                            vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP2<distype>* DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::Instance( bool create )
{
  static FluidEleCalcPoroP2<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcPoroP2<distype>();
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
 *  called upon destruction (public)                            vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}


/*----------------------------------------------------------------------*
 *  constructor (protected)                                   vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::FluidEleCalcPoroP2()
  : DRT::ELEMENTS::FluidEleCalcPoro<distype>::FluidEleCalcPoro(),
    W_(0.0),
    dW_dJ_(0.0)
{

}

/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (2)       vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::Evaluate(
    DRT::ELEMENTS::Fluid*                       ele,
    DRT::Discretization &                       discretization,
    const std::vector<int> &                    lm,
    Teuchos::ParameterList&                     params,
    Teuchos::RCP<MAT::Material> &               mat,
    Epetra_SerialDenseMatrix&                   elemat1_epetra,
    Epetra_SerialDenseMatrix&                   elemat2_epetra,
    Epetra_SerialDenseVector&                   elevec1_epetra,
    Epetra_SerialDenseVector&                   elevec2_epetra,
    Epetra_SerialDenseVector&                   elevec3_epetra,
    const DRT::UTILS::GaussIntegration &        intpoints)
{
  // set element id
  my::eid_ = ele->Id();
  //get structure material
  my::GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  my::rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat1(elemat1_epetra,true);
  //LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat2(elemat2_epetra,true);
  LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> elevec1(elevec1_epetra, true);
  // elevec2 and elevec3 are currently not in use

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
  this->BodyForce(ele,ebofoaf,eprescpgaf,escabofoaf);

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, acceleration
  // and history
  // velocity/pressure values are at time n+alpha_F/n+alpha_M for
  // generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<my::nsd_, my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_, 1> epreaf(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelaf,
      &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<my::nsd_, my::nen_> evelnp(true);
  LINALG::Matrix<my::nen_, 1> eprenp(true);
  if (FluidEleCalc<distype>::fldparatimint_->IsGenalphaNP())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelnp,
        &eprenp, "velnp");

  LINALG::Matrix<my::nsd_, my::nen_> emhist(true);
  LINALG::Matrix<my::nen_, 1> echist(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &emhist,
      &echist, "hist");

  LINALG::Matrix<my::nsd_, my::nen_> eaccam(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eaccam,
      NULL, "accam");

  LINALG::Matrix<my::nen_, 1> epren(true);
  LINALG::Matrix<my::nsd_, my::nen_> eveln(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln,
      &epren, "veln");

  LINALG::Matrix<my::nen_, 1> epressnp_timederiv(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
      &epressnp_timederiv, "accnp");

  LINALG::Matrix<my::nen_,1> escaaf(true);

  LINALG::Matrix<my::nen_,1> eporosity(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, NULL, &eporosity,"scaaf");

  if (not my::fldparatimint_->IsGenalpha())
    eaccam.Clear();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridvn(true);
  LINALG::Matrix<my::nsd_, my::nen_> edispn(true);

 // LINALG::Matrix<my::nen_, 1> eporositynp(true);
 // LINALG::Matrix<my::nen_, 1> eporositydot(true);
 // LINALG::Matrix<my::nen_, 1> eporositydotn(true);

  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispnp,
      NULL, "dispnp");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridv,
      NULL, "gridv");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridvn,
      NULL, "gridvn");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispn,
      NULL, "dispn");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> >(
      ele, my::xyze_);

  my::PreEvaluate(params,ele,discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  return my::Evaluate(
                  params,
                  ebofoaf,
                  elemat1,
                  elevec1,
                  evelaf,
                  epreaf,
                  evelnp,
                  eveln,
                  eprenp,
                  epren,
                  emhist,
                  echist,
                  epressnp_timederiv,
                  eaccam,
                  edispnp,
                  edispn,
                  egridv,
                  egridvn,
                  escaaf,
                  &eporosity,
                  NULL,
                  NULL,
                  mat,
                  ele->IsAle(),
                  intpoints);

}

/*---------------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow off diagonal (2)  vuong 06/11|
 *-------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::EvaluateOD(
    DRT::ELEMENTS::Fluid*                 ele,
    DRT::Discretization &                 discretization,
    const std::vector<int> &              lm,
    Teuchos::ParameterList&               params,
    Teuchos::RCP<MAT::Material> &         mat,
    Epetra_SerialDenseMatrix&             elemat1_epetra,
    Epetra_SerialDenseMatrix&             elemat2_epetra,
    Epetra_SerialDenseVector&             elevec1_epetra,
    Epetra_SerialDenseVector&             elevec2_epetra,
    Epetra_SerialDenseVector&             elevec3_epetra,
    const DRT::UTILS::GaussIntegration &  intpoints)
{
  // set element id
  my::eid_ = ele->Id();

  //get structure material
  my::GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  my::rotsymmpbc_->Setup(ele);

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
  this->BodyForce(ele,ebofoaf,eprescpgaf,escabofoaf);

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure, acceleration
  // and history
  // velocity/pressure values are at time n+alpha_F/n+alpha_M for
  // generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<my::nsd_, my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_, 1> epreaf(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelaf,
      &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<my::nsd_, my::nen_> evelnp(true);
  LINALG::Matrix<my::nen_, 1> eprenp(true);
  if (my::fldparatimint_->IsGenalphaNP())
    this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelnp,
        &eprenp, "velnp");

  LINALG::Matrix<my::nsd_, my::nen_> eveln(true);
  LINALG::Matrix<my::nen_, 1> epren(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln,
      &epren, "veln");

  LINALG::Matrix<my::nen_, 1> epressnp_timederiv(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
      &epressnp_timederiv, "accnp");

  LINALG::Matrix<my::nen_,1> escaaf(true);

  LINALG::Matrix<my::nen_,1> eporosity(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, NULL, &eporosity,"scaaf");

  LINALG::Matrix<my::nsd_, my::nen_> emhist(true);
  LINALG::Matrix<my::nen_, 1> echist(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &emhist,
      &echist, "hist");

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);
  LINALG::Matrix<my::nsd_, my::nen_> edispn(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridvn(true);

  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispnp,
      NULL, "dispnp");
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridv,
      NULL, "gridv");
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispn,
      NULL, "dispn");
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridvn,
      NULL, "gridvn");

  //ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &initporosity_, "initporosity");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> >(
      ele, my::xyze_);

  my::PreEvaluate(params,ele,discretization);

  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  switch(act)
  {
  case FLD::calc_porousflow_fluid_coupling:
  {
    // construct views
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nsd_ * my::nen_> elemat1(elemat1_epetra, true);
    //  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat2(elemat2_epetra,true);
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> elevec1(elevec1_epetra, true);
    // elevec2 and elevec3 are currently not in use

    // call inner evaluate (does not know about DRT element or discretization object)
    return my::EvaluateOD(params,
        ebofoaf,
        elemat1,
        elevec1,
        evelaf,
        epreaf,
        evelnp,
        eveln,
        eprenp,
        epren,
        epressnp_timederiv,
        edispnp,
        edispn,
        egridv,
        egridvn,
        escaaf,
        emhist,
        echist,
        &eporosity,
        mat,
        ele->IsAle(),
        intpoints);
  }
  break;
  case FLD::calc_poroscatra_mono_odblock:
  {
    // construct views
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_> elemat1(elemat1_epetra, true);
    //  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat2(elemat2_epetra,true);
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> elevec1(elevec1_epetra, true);
    // elevec2 and elevec3 are currently not in use

    // call inner evaluate (does not know about DRT element or discretization object)
    return EvaluateODPoroScatra(params,
        elemat1,
        elevec1,
        evelaf,
        epreaf,
        evelnp,
        eprenp,
        epressnp_timederiv,
        edispnp,
        egridv,
        escaaf,
        &eporosity,
        mat,
        ele->IsAle(),
        intpoints);
  }
  break;
  default:
    dserror("Unknown type of action for poro Fluid");
    break;
  }
  return -1;
}

/*----------------------------------------------------------------------*
 *  Compute Porosity                                          vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::ComputePorosity(
                                          Teuchos::ParameterList&             params,
                                          const double&                       press,
                                          const double&                       J,
                                          const int&                          gp,
                                          const LINALG::Matrix<my::nen_,1>&   shapfct,
                                          const LINALG::Matrix<my::nen_,1>*   myporosity,
                                          double&                             porosity,
                                          double*                             dphi_dp,
                                          double*                             dphi_dJ,
                                          double*                             dphi_dJdp,
                                          double*                             dphi_dJJ,
                                          double*                             dphi_dpp,
                                          bool                                save)
{
  if(myporosity == NULL)
    dserror("no porosity values given!!");
  else
    porosity = shapfct.Dot(*myporosity);

  return;
}

/*----------------------------------------------------------------------*
 * Compute Porosity gradient                                  vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::ComputePorosityGradient(
                        const double&                                      dphidp,
                        const double&                                      dphidJ,
                        const LINALG::Matrix<my::nsd_,1>&                  gradJ,
                        const LINALG::Matrix<my::nen_,1>*                  eporositynp)
{
  if(eporositynp == NULL)
    dserror("no porosity values given for calculation of porosity gradient!!");

  //if( (my::fldpara_->PoroContiPartInt() == false) or my::visceff_)
  {
    //--------------------- current porosity gradient
    my::grad_porosity_.Multiply(my::derxy_,*eporositynp);
  }
  my::refgrad_porosity_.Multiply(my::xjm_,my::grad_porosity_);
}


/*----------------------------------------------------------------------*
 *  PSPG contributions                                      vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::PSPG(
    LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &         estif_q_u,
    LINALG::Matrix<my::nen_,my::nen_> &                   ppmat,
    LINALG::Matrix<my::nen_,1> &                          preforce,
    const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &    lin_resM_Du,
    const LINALG::Matrix<my::nsd_,my::nen_> &             lin_resM_Dp,
    const double &                                        fac3,
    const double &                                        dphi_dp,
    const double &                                        timefacfac,
    const double &                                        timefacfacpre,
    const double &                                        rhsfac)
{
  dserror("PSPG for poro p2 not working");
//  my::PSPG( estif_q_u,
//            ppmat,
//            preforce,
//            lin_resM_Du,
//            fac3,
//            timefacfac,
//            timefacfacpre,
//            rhsfac);

  double scal_grad_q=0.0;

  if(my::fldpara_->Tds()==INPAR::FLUID::subscales_quasistatic)
  {
    scal_grad_q=0.0;//my::tau_(1);
  }
  else
  {
    scal_grad_q=my::fldparatimint_->AlphaF()*fac3;
  }

  for(int jdim=0;jdim<my::nsd_;++jdim)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui_p_jdim   = my::nsd_*ui + jdim;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double temp_vi_idim=my::derxy_(jdim,vi)*scal_grad_q;

        estif_q_u(vi,fui_p_jdim) += timefacfacpre*my::conres_old_*my::funct_(ui)*temp_vi_idim;
      } // vi
    } // ui
  } //jdim

  const double temp = rhsfac*scal_grad_q*my::conres_old_;
  for (int idim = 0; idim <my::nsd_; ++idim)
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      // pressure stabilisation
      preforce(vi) -= temp*my::derxy_(idim, vi)*my::velint_(idim);
    }
  } // end for(idim)
}

/*----------------------------------------------------------------------*
 * reactive stabilization                                      vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::ReacStab(
    LINALG::Matrix<my::nen_*my::nsd_,my::nen_*my::nsd_> &     estif_u,
    LINALG::Matrix<my::nen_*my::nsd_,my::nen_> &              estif_p_v,
    LINALG::Matrix<my::nsd_,my::nen_> &                       velforce,
    LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &              lin_resM_Du,
    const LINALG::Matrix<my::nsd_,my::nen_> &                 lin_resM_Dp,
    const double &                                            dphi_dp,
    const double &                                            timefacfac,
    const double &                                            timefacfacpre,
    const double &                                            rhsfac,
    const double &                                            fac3)
{
  my::ReacStab(estif_u,
           estif_p_v,
           velforce,
           lin_resM_Du,
           lin_resM_Dp,
           dphi_dp,
           timefacfac,
           timefacfacpre,
           rhsfac,
           fac3);

  double reac_tau;
  if (my::fldpara_->Tds()==INPAR::FLUID::subscales_quasistatic)
  {
    dserror("PSPG for poro p2 not checked");
    reac_tau = 0.0;//my::fldpara_->ViscReaStabFac()*my::reacoeff_*my::tau_(1);
  }
  else
  {
    dserror("Is this factor correct? Check for bugs!");
    reac_tau = my::fldpara_->ViscReaStabFac()*my::reacoeff_*my::fldparatimint_->AlphaF()*fac3;
  }

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const double v = reac_tau*my::funct_(vi)*my::conres_old_;

    for(int idim=0;idim<my::nsd_;++idim)
    {
      const int fvi_p_idim = my::nsd_*vi+idim;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui_p_idim   = my::nsd_*ui + idim;

        estif_u(fvi_p_idim,fui_p_idim) += v*my::funct_(idim);
      } // ui
    } //idim
  } // vi

  const double reac_fac = reac_tau*rhsfac;
  const double v = reac_fac*my::conres_old_;
  for (int idim =0;idim<my::nsd_;++idim)
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
        velforce(idim,vi) -= v*my::funct_(vi)*my::velint_(idim);
    }
  } // end for(idim)
}

/*----------------------------------------------------------------------*
 * Evaluate Pressure Equation                                 vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::EvaluatePressureEquation(
    Teuchos::ParameterList&                       params,
    const double&                                 timefacfacpre,
    const double&                                 rhsfac,
    const double&                                 dphi_dp,
    const double&                                 dphi_dJ,
    const double&                                 dphi_dJdp,
    const double&                                 dphi_dpp,
    const LINALG::Matrix<my::nen_,1>*             eporositydot,
    const LINALG::Matrix<my::nen_,1>*             eporositydotn,
    const LINALG::Matrix<my::nen_,1>&             echist,
    const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
    LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
    LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
    LINALG::Matrix<my::nen_,1>&                   preforce
    )
{
  double    dW_dp   = 0.0;
  my::structmat_->ConstitutiveDerivatives(params,
                                    my::press_,
                                    my::J_,
                                    my::porosity_,
                                    &dW_dp,
                                    NULL, // not needed
                                    &dW_dJ_,
                                    &W_
                                    );
  //--------------------------------------------------------

  // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

  for (int k=0; k<my::nen_; k++)
  {
    const double v = timefacfacpre* my::funct_(k);
    const double w = rhsfac* my::funct_(k);

    for (int i=0; i<my::nen_; i++)
      ppmat(k,i) += v * dW_dp * my::funct_(i);

    preforce(k) -= w*W_;
  }

  my::rhscon_ = W_;
  my::histcon_ = 0.0;
}

/*----------------------------------------------------------------------*
 * Evaluate Pressure Equation off diagonal                   vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::FillMatrixContiOD(
    const double&                                               timefacfacpre,
    const double &                                              dphi_dp,
    const double &                                              dphi_dJ,
    const double&                                               dphi_dJJ,
    const double&                                               dphi_dJdp,
    const double &                                              refporositydot,
    const LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&           dgradphi_dus,
    const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dphi_dus,
    const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dJ_dus,
    const LINALG::Matrix<my::nsd_, my::nen_>&                   egridv,
    const LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>&        lin_resM_Dus,
    LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&              ecoupl_p
    )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const double w = timefacfacpre*my::funct_(vi);
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = my::nsd_*ui;
      for(int idim = 0; idim <my::nsd_; ++idim)
      {
        ecoupl_p(vi,fui+idim)+=  w * dW_dJ_ * dJ_dus(fui+idim);
      }
    }
  } // end for(idim)

  //shape derivatives
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    const double v = timefacfacpre * my::funct_(vi) * W_;
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      const int ui_nsd  = ui * my::nsd_;
      for (int idim = 0; idim < my::nsd_; ++idim)
        ecoupl_p(vi, ui_nsd + idim) += v * my::derxy_(idim, ui);
    }
  }
}

/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow with scatra(3)    vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::EvaluateODPoroScatra(
    Teuchos::ParameterList&                                           params,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_> &             elemat1,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
    const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
    const LINALG::Matrix<my::nen_, 1> &                               epreaf,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
    const LINALG::Matrix<my::nen_, 1> &                               eprenp,
    const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
    const LINALG::Matrix<my::nen_,1>&                                 escaaf,
    const LINALG::Matrix<my::nen_,1>*                                 eporositynp,
    Teuchos::RCP<MAT::Material>                                       mat,
    bool                                                              isale,
    const DRT::UTILS::GaussIntegration &                              intpoints)
{
  // flag for higher order elements
  my::is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (my::fldpara_->IsInconsistent() == true)
    my::is_higher_order_ele_ = false;

  // stationary formulation does not support ALE formulation
  //if (isale and my::fldparatimint_->IsStationary())
  //  dserror("No ALE support within stationary fluid solver.");

    // ---------------------------------------------------------------------
    // call routine for calculating element matrix and right hand side
    // ---------------------------------------------------------------------
    SysmatPoroScatraOD(
        params,
        evelaf,
        evelnp,
        epreaf,
        eprenp,
        epressnp_timederiv,
        edispnp,
        egridv,
        escaaf,
        eporositynp,
        elemat1,
        elevec1,
        mat,
        isale,
        intpoints);

    return 0;
}

/*----------------------------------------------------------------------*
 |  calculate coupling matrix flow                          vuong 06/11 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP2<distype>::SysmatPoroScatraOD(
    Teuchos::ParameterList&                                         params,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
    const LINALG::Matrix<my::nen_, 1>&                              epreaf,
    const LINALG::Matrix<my::nen_, 1>&                              eprenp,
    const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
    const LINALG::Matrix<my::nen_,1>&                               escaaf,
    const LINALG::Matrix<my::nen_,1>*                               eporositynp,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_,my::nen_>&             ecoupl,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
    Teuchos::RCP<const MAT::Material>                               material,
    bool                                                            isale,
    const DRT::UTILS::GaussIntegration &                            intpoints)
{
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // evaluate shape function derivatives w.r.t. to material coordinates at integration point
    //const double det0 = my::SetupMaterialDerivatives();

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including my::derivatives and grid velocity)
    //  2) pressure (including my::derivatives)
    //  3) body-force vector
    //  4) "history" vector for momentum equation
    //----------------------------------------------------------------------
    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::velint_.Multiply(evelaf,my::funct_);

    // get velocity my::derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::vderxy_.MultiplyNT(evelaf,my::derxy_);

    // get convective velocity at integration point
    // (ALE case handled implicitly here using the (potential
    //  mesh-movement-dependent) convective velocity, avoiding
    //  various ALE terms used to be calculated before)
    // convvelint_.Update(my::velint_);
    my::convvelint_.Multiply(-1.0, egridv, my::funct_, 0.0);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    if(my::fldparatimint_->IsGenalphaNP())
      my::press_ = my::funct_.Dot(eprenp);
    else
      my::press_ = my::funct_.Dot(epreaf);

    // get pressure time my::derivative at integration point
    // (value at n+1 )
    my::pressdot_ = my::funct_.Dot(epressnp_timederiv);

    // get pressure gradient at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    if(my::fldparatimint_->IsGenalphaNP())
      my::gradp_.Multiply(my::derxy_,eprenp);
    else
      my::gradp_.Multiply(my::derxy_,epreaf);

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::gridvelint_.Multiply(egridv,my::funct_);

    // get displacement my::derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    LINALG::Matrix<my::nsd_,my::nsd_> gridvelderxy;
    gridvelderxy.MultiplyNT(egridv,my::derxy_);

    //----------------------------------------------------------------------
    //  evaluation of various partial operators at integration point
    //  1) convective term from previous iteration (mandatorily set to zero)
    //  2) viscous term from previous iteration and viscous operator
    //  3) divergence of velocity from previous iteration
    //----------------------------------------------------------------------
    // set convective term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    //  conv_old_.Clear();

    //set old convective term to ALE-Term only
    my::conv_old_.Multiply(my::vderxy_,my::convvelint_);
    my::conv_c_.MultiplyTN(my::derxy_,my::convvelint_);

    // set viscous term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    my::visc_old_.Clear();

    // compute divergence of velocity from previous iteration
    my::vdiv_ = 0.0;
    // double dispdiv=0.0;
    my::gridvdiv_ = 0.0;
    if (not my::fldparatimint_->IsGenalphaNP())
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        my::vdiv_ += my::vderxy_(idim, idim);

        my::gridvdiv_ += gridvelderxy(idim,idim);
      }
    }
    else
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        //get vdiv at time n+1 for np_genalpha,
        LINALG::Matrix<my::nsd_,my::nsd_> vderxy;
        vderxy.MultiplyNT(evelnp,my::derxy_);
        my::vdiv_ += vderxy(idim, idim);

        my::gridvdiv_ += gridvelderxy(idim,idim);
      }
    }

//    // -------------------------(material) deformation gradient F = d my::xyze_ / d XYZE = my::xyze_ * N_XYZ_^T
//    LINALG::Matrix<my::nsd_,my::nsd_> defgrd(false);
//    defgrd.MultiplyNT(my::xyze_,my::N_XYZ_);
//
//    // inverse deformation gradient F^-1
//    LINALG::Matrix<my::nsd_,my::nsd_> defgrd_inv(false);
//    defgrd_inv.Invert(defgrd);
//
//    //------------------------------------ build F^-T as vector 9x1
//    LINALG::Matrix<my::nsd_*my::nsd_,1> defgrd_IT_vec;
//    for(int i=0; i<my::nsd_; i++)
//      for(int j=0; j<my::nsd_; j++)
//        defgrd_IT_vec(i*my::nsd_+j) = defgrd_inv(j,i);
//
//    // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
//    my::J_ = my::det_/det0;
//
//    //************************************************auxilary variables for computing the porosity
//
//    double dphi_dp=0.0;
//    double dphi_dJ=0.0;
//    double dphi_dJdp=0.0;
//    double dphi_dJJ=0.0;
//    my::porosity_=0.0;
//
//    // compute scalar at n+alpha_F or n+1
//    const double scalaraf = my::funct_.Dot(escaaf);
//    params.set<double>("scalar",scalaraf);
//    ComputePorosity(  params,
//                      my::press_,
//                      my::J_,
//                      *(iquad),
//                      my::funct_,
//                      eporositynp,
//                      my::porosity_,
//                      &dphi_dp,
//                      &dphi_dJ,
//                      &dphi_dJdp,
//                      &dphi_dJJ,
//                      NULL, //dphi_dpp not needed
//                      false);
//
//    double refporositydot = my::so_interface_->RefPorosityTimeDeriv();

//    //---------------------------  dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx at gausspoint
//    LINALG::Matrix<my::nsd_,1> gradJ(true);
//    // spatial porosity gradient
//    LINALG::Matrix<my::nsd_,1>             grad_porosity(true);
//    //--------------------- linearization of porosity w.r.t. structure displacements
//    LINALG::Matrix<1,my::nsd_*my::nen_> dphi_dus(true);
//
//    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J * N_x
//    LINALG::Matrix<1,my::nsd_*my::nen_> dJ_dus(true);
//    //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) + d\phi/(dp*du_s) * dp/dx
//    LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> dgradphi_dus(true);
//
//    {
//      // dF/dx
//      LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_> F_x(true);
//
//      // dF/dX
//      LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_> F_X(true);
//
//      my::ComputeFDerivative( edispnp,
//                          defgrd_inv,
//                          F_x,
//                          F_X);
//
//      //compute gradients if needed
//      my::ComputeGradients(
//                       dphi_dp,
//                       dphi_dJ,
//                       defgrd_IT_vec,
//                       F_x,
//                       eporositynp,
//                       gradJ);
//
//      my::ComputeLinearizationOD(
//                              dphi_dJ,
//                              dphi_dJJ,
//                              dphi_dJdp,
//                              defgrd_inv,
//                              defgrd_IT_vec,
//                              F_x,
//                              F_X,
//                              gradJ,
//                              dJ_dus,
//                              dphi_dus,
//                              dgradphi_dus);
//    }
    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point

    double reacoeff = 0.0;

    if (my::fldpara_->MatGp())
    {
      Teuchos::RCP<const MAT::FluidPoro> actmat = Teuchos::rcp_static_cast<const MAT::FluidPoro>(material);

      // set density at n+alpha_F/n+1 and n+alpha_M/n+1
      //my::densaf_ = actmat->Density();
      //my::densam_ = my::densaf_;

      // calculate reaction coefficient
      reacoeff = actmat->ComputeReactionCoeff();
    }
    else dserror("Fluid material parameters have to be evaluated at gauss point for porous flow!");

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = my::fldparatimint_->TimeFac() * my::fac_;
    const double timefacfacpre = my::fldparatimint_->TimeFacPre() * my::fac_;


    for (int ui=0; ui<my::nen_; ++ui)
    {
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = (my::nsd_+1)*vi;
        const double tmp = my::funct_(vi)* reacoeff;
        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          ecoupl(fvi+idim,ui) += timefacfac * tmp * my::velint_(idim) * my::funct_(ui);
        } // end for (idim)ecoupl
      } //vi
    } // ui


    if (not my::fldparatimint_->IsStationary())
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = my::nsd_*vi;
          const double tmp = my::funct_(vi)* reacoeff;
          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            ecoupl(fvi+idim,ui) += timefacfac * tmp * (-my::gridvelint_(idim)) * my::funct_(ui);
          } // end for (idim)
        } //vi
      } // ui
    }

    for (int ui=0; ui<my::nen_; ++ui)
      for (int vi=0; vi<my::nen_; ++vi)
        ecoupl(vi,ui) +=   my::fac_ * my::funct_(vi) * my::funct_(ui);

    LINALG::Matrix<my::nen_,1>    derxy_convel(true);

    for (int i =0; i< my::nen_; i++)
      for (int j =0; j< my::nsd_; j++)
        derxy_convel(i) += my::derxy_(j,i) * my::velint_(j);

    if (not my::fldparatimint_->IsStationary())
    {
      for (int i =0; i< my::nen_; i++)
        for (int j =0; j< my::nsd_; j++)
          derxy_convel(i) += my::derxy_(j,i) * (-my::gridvelint_(j));
    }

    if( static_cast<DRT::ELEMENTS::FluidEleParameterPoro*>(my::fldpara_)->PoroContiPartInt() == false )
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          ecoupl(vi,ui) +=
                               + timefacfacpre * my::vdiv_ * my::funct_(vi) * my::funct_(ui)
                               + timefacfacpre * my::funct_(vi) * derxy_convel(ui)
                                 ;
        }
      }
    }
    else //my::fldpara_->PoroContiPartInt() == true
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          ecoupl(vi,ui) += -1.0 * timefacfacpre * derxy_convel(vi) * my::funct_(ui)
                                 ;
        }
      }

      if (not my::fldparatimint_->IsStationary())
      {
        for (int ui=0; ui<my::nen_; ++ui)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            ecoupl(vi,ui) += timefacfacpre * my::funct_(vi) * my::gridvdiv_ * my::funct_(ui)
                                   ;
          }
        }
      }
    }
  }//loop over gausspoints
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcPoroP2<DRT::Element::nurbs27>;



