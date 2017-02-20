/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_calc_poro_p1.cpp

 \brief Internal implementation of poro Fluid element (p1 poro fluid)

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
 *----------------------------------------------------------------------*/

#include "fluid_ele_calc_poro_p1.H"

#include "fluid_ele.H"
#include "fluid_ele_parameter_poro.H"
#include "../drt_lib/drt_element_integration_select.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_geometry/position_array.H"


/*----------------------------------------------------------------------*
 *  create/delete instance (public)                            vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP1<distype>* DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Instance( bool create )
{
  static FluidEleCalcPoroP1<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcPoroP1<distype>();
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
 * called upon destruction (public)                            vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}


/*----------------------------------------------------------------------*
 *   constructor (protected)                                   vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::FluidEleCalcPoroP1()
  : DRT::ELEMENTS::FluidEleCalcPoro<distype>::FluidEleCalcPoro()
{

}

/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (2)
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Evaluate(
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
  static LINALG::Matrix<my::nsd_,my::nen_> ebofoaf(true);
  ebofoaf.Clear();
  static LINALG::Matrix<my::nsd_,my::nen_> eprescpgaf(true);
  eprescpgaf.Clear();
  static LINALG::Matrix<my::nen_,1>    escabofoaf(true);
  escabofoaf.Clear();
  my::BodyForce(ele,ebofoaf,eprescpgaf,escabofoaf);

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
  static LINALG::Matrix<my::nsd_, my::nen_> evelaf(true);
  evelaf.Clear();
  static LINALG::Matrix<my::nen_, 1> epreaf(true);
  epreaf.Clear();
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelaf,
      &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  static LINALG::Matrix<my::nsd_, my::nen_> evelnp(true);
  evelnp.Clear();
  static LINALG::Matrix<my::nen_, 1> eprenp(true);
  eprenp.Clear();
  if (FluidEleCalc<distype>::fldparatimint_->IsGenalphaNP())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelnp,
        &eprenp, "velnp");

  static LINALG::Matrix<my::nsd_, my::nen_> emhist(true);
  emhist.Clear();
  static LINALG::Matrix<my::nen_, 1> echist(true);
  echist.Clear();
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &emhist,
      &echist, "hist");

  static LINALG::Matrix<my::nsd_, my::nen_> eaccam(true);
  static LINALG::Matrix<my::nen_,1> epressam_timederiv(true);
  eaccam.Clear();
  epressam_timederiv.Clear();

  if (my::fldparatimint_->IsGenalpha())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eaccam,
        &epressam_timederiv, "accam");

  static LINALG::Matrix<my::nen_,1> epressn_timederiv(true);
  epressn_timederiv.Clear();
  if (my::fldparatimint_->IsGenalpha())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
        &epressn_timederiv, "accn");

  static LINALG::Matrix<my::nen_, 1> epren(true);
  epren.Clear();
  static LINALG::Matrix<my::nsd_, my::nen_> eveln(true);
  eveln.Clear();
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln,
      &epren, "veln");

  static LINALG::Matrix<my::nen_, 1> epressnp_timederiv(true);
  epressnp_timederiv.Clear();
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
      &epressnp_timederiv, "accnp");

  static LINALG::Matrix<my::nen_,1> escaaf(true);
  escaaf.Clear();
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, NULL, &escaaf,"scaaf");

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  static LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  edispnp.Clear();
  static LINALG::Matrix<my::nsd_, my::nen_> egridv(true);
  egridv.Clear();
  static LINALG::Matrix<my::nsd_, my::nen_> egridvn(true);
  egridvn.Clear();
  static LINALG::Matrix<my::nsd_, my::nen_> edispn(true);
  edispn.Clear();

  static LINALG::Matrix<my::nen_, 1> eporositynp(true);
  eporositynp.Clear();
  static LINALG::Matrix<my::nen_, 1> eporositydot(true);
  eporositydot.Clear();
  static LINALG::Matrix<my::nen_, 1> eporositydotn(true);
  eporositydotn.Clear();

  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispnp,
      &eporositynp, "dispnp");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridv,
      &eporositydot, "gridv");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridvn,
      &eporositydotn, "gridvn");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispn,
      NULL, "dispn");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> >(
      ele, my::xyze_);

  // construct views
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat1(elemat1_epetra,true);
  //LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat2(elemat2_epetra,true);
  LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> elevec1(elevec1_epetra, true);
  // elevec2 and elevec3 are currently not in use

  my::PreEvaluate(params,ele,discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = my::Evaluate(
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
                  epressam_timederiv,
                  epressn_timederiv,
                  eaccam,
                  edispnp,
                  edispn,
                  egridv,
                  egridvn,
                  escaaf,
                  &eporositynp,
                  &eporositydot,
                  &eporositydotn,
                  mat,
                  ele->IsAle(),
                  intpoints);

  return result;
}


/*----------------------------------------------------------------------*
 *  Compute Porosity                                          vuong 06/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputePorosity(
                                          Teuchos::ParameterList& params,
                                          const double& press,
                                          const double& J,
                                          const int& gp,
                                          const LINALG::Matrix<my::nen_,1>&       shapfct,
                                          const LINALG::Matrix<my::nen_,1>*       myporosity,
                                          double& porosity,
                                          double* dphi_dp,
                                          double* dphi_dJ,
                                          double* dphi_dJdp,
                                          double* dphi_dJJ,
                                          double* dphi_dpp,
                                          bool save)
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
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputePorosityGradient(
    const double&                                      dphidp,
    const double&                                      dphidJ,
    const LINALG::Matrix<my::nsd_,1>&                  gradJ,
    const LINALG::Matrix<my::nsd_,1>&                  gradp,
    const LINALG::Matrix<my::nen_,1>*                  eporositynp,
    LINALG::Matrix<my::nsd_,1>&                        grad_porosity,
    LINALG::Matrix<my::nsd_,1>&                        refgrad_porosity)
{
  if(eporositynp == NULL)
    dserror("no porosity values given for calculation of porosity gradient!!");

  //if( (my::fldpara_->PoroContiPartInt() == false) or my::visceff_)
  {
    //--------------------- current porosity gradient
    grad_porosity.Multiply(my::derxy_,*eporositynp);
  }

  refgrad_porosity.Multiply(my::xjm_,grad_porosity);
}

/*----------------------------------------------------------------------*
 *  Evaluate Pressure Equation                                 vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::EvaluatePressureEquation(
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

  // first evaluate terms without porosity time derivative
  my::EvaluatePressureEquationNonTransient(params,
                                    timefacfacpre,
                                    rhsfac,
                                    dphi_dp,
                                    dphi_dJ,
                                    dphi_dJdp,
                                    dphi_dpp,
                                    dgradphi_dp,
                                    estif_q_u,
                                    ppmat,
                                    preforce);

  // now the porosity time derivative (different for standard poro and poro_p1 elements)
  if (my::porofldpara_->IsStationaryConti() == false)
  {
    // inertia terms on the right hand side for instationary fluids

    if(eporositydot)
    {
      double porositydot =  my::funct_.Dot(*eporositydot);
      //double porositydot =  my::funct_.Dot(*eporositydotn);

      for (int vi=0; vi<my::nen_; ++vi)
      {//TODO : check genalpha case
        preforce(vi)-=  rhsfac * porositydot * my::funct_(vi) ;
      }

      //no need for adding RHS form previous time step, as it is already included in 'porositydot'
      //(for the one-step-theta case at least)
      //ComputeContiTimeRHS(params,*eporositydotn,preforce,rhsfac,1.0);

      //just update internal variables, no contribution to rhs
      const double porositydotn = my::funct_.Dot(*eporositydotn);

      my::histcon_ = my::fldparatimint_->OmTheta() * my::fldparatimint_->Dt() * porositydotn;

      //rhs from last time step
      my::rhscon_ = 1.0/my::fldparatimint_->Dt()/my::fldparatimint_->Theta() * my::histcon_;

      //transient part of continuity equation residual
      my::conres_old_ += porositydot - my::rhscon_;
    }
    else
      dserror("no porosity time derivative given for poro_p1 element!");
  }

  return;
}

/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (2)          vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::EvaluateOD(
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

  // set element id
  my::eid_ = ele->Id();

  //get structure material
  my::GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  my::rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1)* my::nen_> elemat1(elemat1_epetra, true);
  //  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat2(elemat2_epetra,true);
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
  if (my::fldparatimint_->IsGenalphaNP())
    this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelnp,
        &eprenp, "velnp");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<my::nsd_, my::nen_> eveln(true);
  LINALG::Matrix<my::nen_, 1> epren(true);
  if (my::fldparatimint_->IsGenalphaNP())
    this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln,
        &epren, "veln");

  static LINALG::Matrix<my::nsd_, my::nen_> eaccam(true);
  static LINALG::Matrix<my::nen_,1> epressam_timederiv(true);
  eaccam.Clear();
  epressam_timederiv.Clear();

  if (my::fldparatimint_->IsGenalpha())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eaccam,
        &epressam_timederiv, "accam");

  static LINALG::Matrix<my::nen_,1> epressn_timederiv(true);
  epressn_timederiv.Clear();
  if (my::fldparatimint_->IsGenalpha())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
        &epressn_timederiv, "accn");

  LINALG::Matrix<my::nen_, 1> epressnp_timederiv(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
      &epressnp_timederiv, "accnp");

  LINALG::Matrix<my::nen_,1> escaaf(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, NULL, &escaaf,"scaaf");

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

  LINALG::Matrix<my::nen_, 1> eporositynp(true);

  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispnp,
      &eporositynp, "dispnp");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridv,
      NULL, "gridv");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispn,
      NULL, "dispn");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridvn,
      NULL, "gridvn");

  //ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &initporosity_, "initporosity");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> >(
      ele, my::xyze_);

  my::PreEvaluate(params,ele,discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = EvaluateOD(params,
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
      epressam_timederiv,
      epressn_timederiv,
      eaccam,
      edispnp,
      edispn,
      egridv,
      egridvn,
      escaaf,
      &eporositynp,
      mat,
      ele->IsAle(),
      intpoints);

  return result;
}

/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (3)          vuong 06/11|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::EvaluateOD(
    Teuchos::ParameterList&                                           params,
    const LINALG::Matrix<my::nsd_,my::nen_> &                         ebofoaf,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_> &  elemat1,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
    const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
    const LINALG::Matrix<my::nen_, 1> &                               epreaf,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        eveln,
    const LINALG::Matrix<my::nen_, 1> &                               eprenp,
    const LINALG::Matrix<my::nen_, 1> &                               epren,
    const LINALG::Matrix<my::nsd_,my::nen_> &                         emhist,
    const LINALG::Matrix<my::nen_,1>&                                 echist,
    const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
    const LINALG::Matrix<my::nen_, 1> &                               epressam_timederiv,
    const LINALG::Matrix<my::nen_, 1> &                               epressn_timederiv,
    const LINALG::Matrix<my::nsd_,my::nen_> &                         eaccam,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        edispn,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        egridvn,
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
  SysmatOD(
      params,
      ebofoaf,
      evelaf,
      evelnp,
      eveln,
      epreaf,
      eprenp,
      epren,
      emhist,
      echist,
      epressnp_timederiv,
      epressam_timederiv,
      epressn_timederiv,
      eaccam,
      edispnp,
      edispn,
      egridv,
      egridvn,
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
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::SysmatOD(
    Teuchos::ParameterList&                                         params,
    const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
    const LINALG::Matrix<my::nen_, 1>&                              epreaf,
    const LINALG::Matrix<my::nen_, 1>&                              eprenp,
    const LINALG::Matrix<my::nen_, 1>&                              epren,
    const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
    const LINALG::Matrix<my::nen_,1>&                               echist,
    const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
    const LINALG::Matrix<my::nen_, 1> &                             epressam_timederiv,
    const LINALG::Matrix<my::nen_, 1> &                             epressn_timederiv,
    const LINALG::Matrix<my::nsd_,my::nen_> &                       eaccam,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
    const LINALG::Matrix<my::nen_,1>&                               escaaf,
    const LINALG::Matrix<my::nen_,1>*                               eporositynp,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_,(my::nsd_ + 1) * my::nen_>&  ecoupl,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
    Teuchos::RCP<const MAT::Material>                               material,
    bool                                                            isale,
    const DRT::UTILS::GaussIntegration &                            intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  static LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_> ecoupl_u(true); // coupling matrix for momentum equation
  static LINALG::Matrix<my::nen_, my::nen_ * my::nsd_> ecoupl_p(true); // coupling matrix for continuity equation
  //LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_ * my::nsd_> emesh(true); // linearisation of mesh motion

  static LINALG::Matrix<my::nen_ * my::nsd_, my::nen_> ecouplp1_u(true); // coupling matrix for momentum equation
  static LINALG::Matrix<my::nen_, my::nen_> ecouplp1_p(true); // coupling matrix for continuity equation

  ecoupl_u.Clear();
  ecoupl_p.Clear();
  ecouplp1_u.Clear();
  ecouplp1_p.Clear();

  //material coordinates xyze0
  my::xyze0_ = my::xyze_;

  // add displacement when fluid nodes move in the ALE case (in poroelasticity this is always the case)
  my::xyze_ += edispnp;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter();

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  GaussPointLoopP1OD(  params,
                       ebofoaf,
                       evelaf,
                       evelnp,
                       eveln,
                       epreaf,
                       eprenp,
                       epren,
                       emhist,
                       echist,
                       epressnp_timederiv,
                       epressam_timederiv,
                       epressn_timederiv,
                       eaccam,
                       edispnp,
                       edispn,
                       egridv,
                       egridvn,
                       escaaf,
                       eporositynp,
                       eforce,
                       ecoupl_u,
                       ecoupl_p,
                       ecouplp1_u,
                       ecouplp1_p,
                       material,
                       intpoints);
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix
  //------------------------------------------------------------------------

  // add fluid velocity-structure displacement part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int nsd_ui = my::nsd_ *ui;
    const int nsdp1_ui = (my::nsd_ + 1)*ui;

    for (int jdim=0; jdim < my::nsd_;++jdim)
    {
      const int nsd_ui_jdim = nsd_ui+jdim;
      const int nsdp1_ui_jdim = nsdp1_ui+jdim;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int numdof_vi = my::numdofpernode_*vi;
        const int nsd_vi = my::nsd_*vi;

        for (int idim=0; idim <my::nsd_; ++idim)
        {
          ecoupl(numdof_vi+idim, nsdp1_ui_jdim) += ecoupl_u(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add fluid pressure-structure displacement part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int nsd_ui = my::nsd_ *ui;
    const int nsdp1_ui = (my::nsd_ + 1)*ui;

    for (int jdim=0; jdim < my::nsd_;++jdim)
    {
      const int nsd_ui_jdim = nsd_ui+jdim;
      const int nsdp1_ui_jdim = nsdp1_ui+jdim;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        ecoupl(my::numdofpernode_*vi+my::nsd_, nsdp1_ui_jdim) += ecoupl_p(vi, nsd_ui_jdim);
      }
    }
  }

  // add fluid velocity-structure porosity part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int nsdp1_ui = (my::nsd_ + 1)*ui;

    for (int idim=0; idim < my::nsd_;++idim)
    {
      const int nsdp1_ui_nsd = nsdp1_ui+my::nsd_;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int numdof_vi = my::numdofpernode_*vi;
        const int nsd_vi = my::nsd_*vi;

        ecoupl(numdof_vi+idim, nsdp1_ui_nsd) += ecouplp1_u(nsd_vi+idim, ui);
      }
    }
  }

  // add fluid pressure-structure porosity part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int nsdp1_ui_nsd = (my::nsd_ + 1)*ui+my::nsd_;

    for (int vi=0; vi<my::nen_; ++vi)
      ecoupl(my::numdofpernode_*vi+my::nsd_, nsdp1_ui_nsd) += ecouplp1_p(vi, ui);
  }

  return;
}    //SysmatOD

/*----------------------------------------------------------------------*
 |  calculate coupling matrix flow                          vuong 06/11 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::GaussPointLoopP1OD(
                        Teuchos::ParameterList&                                         params,
                        const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
                        const LINALG::Matrix<my::nen_, 1>&                              epreaf,
                        const LINALG::Matrix<my::nen_, 1>&                              eprenp,
                        const LINALG::Matrix<my::nen_, 1>&                              epren,
                        const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
                        const LINALG::Matrix<my::nen_,1>&                               echist,
                        const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
                        const LINALG::Matrix<my::nen_, 1> &                             epressam_timederiv,
                        const LINALG::Matrix<my::nen_, 1> &                             epressn_timederiv,
                        const LINALG::Matrix<my::nsd_,my::nen_>&                        eaccam,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
                        const LINALG::Matrix<my::nen_,1>&                               escaaf,
                        const LINALG::Matrix<my::nen_,1>*                               eporositynp,
                        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
                        LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&       ecoupl_u,
                        LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&                  ecoupl_p,
                        LINALG::Matrix<my::nen_ * my::nsd_, my::nen_>&                  ecouplp1_u,
                        LINALG::Matrix<my::nen_, my::nen_>&                             ecouplp1_p,
                        Teuchos::RCP<const MAT::Material>                               material,
                        const DRT::UTILS::GaussIntegration &                            intpoints)
{
  // definition of velocity-based momentum residual vectors
  static LINALG::Matrix< my::nsd_, my::nen_ * my::nsd_>  lin_resM_Dus(true);
  static LINALG::Matrix< my::nsd_, my::nen_ * my::nsd_>  lin_resM_Dus_gridvel(true);
  static LINALG::Matrix< my::nsd_, my::nen_>  lin_resM_Dphi(true);

  // set element area or volume
  const double vol = my::fac_;

  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    lin_resM_Dus.Clear();
    lin_resM_Dus_gridvel.Clear();
    lin_resM_Dphi.Clear();

    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // evaluate shape function derivatives w.r.t. to material coordinates at integration point
    my::SetupMaterialDerivatives();

    // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ * N_XYZ_^T
    static LINALG::Matrix<my::nsd_,my::nsd_>          defgrd(false);
    my::ComputeDefGradient(defgrd,my::N_XYZ_,my::xyze_);

    // inverse deformation gradient F^-1
    static LINALG::Matrix<my::nsd_,my::nsd_>          defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J and the volume change
    my::ComputeJacobianDeterminantVolumeChange(
        my::J_,
        volchange,
        defgrd,
        my::N_XYZ_,
        edispnp);

    my::EvaluateVariablesAtGaussPointOD(
          params,
          ebofoaf,
          evelaf,
          evelnp,
          eveln,
          epreaf,
          eprenp,
          epren,
          epressnp_timederiv,
          epressam_timederiv,
          epressn_timederiv,
          eaccam,
          edispnp,
          edispn,
          egridv,
          egridvn,
          escaaf,
          emhist,
          echist,
          eporositynp);

    //************************************************auxilary variables for computing the porosity

    double dphi_dp=0.0;
    double dphi_dJ=0.0;
    double dphi_dJdp=0.0;
    double dphi_dJJ=0.0;
    my::porosity_=0.0;

    // compute scalar at n+alpha_F or n+1
    const double scalaraf = my::funct_.Dot(escaaf);
    params.set<double>("scalar",scalaraf);
    ComputePorosity(  params,
                      my::press_,
                      volchange,
                      *(iquad),
                      my::funct_,
                      eporositynp,
                      my::porosity_,
                      &dphi_dp,
                      &dphi_dJ,
                      &dphi_dJdp,
                      &dphi_dJJ,
                      NULL, //dphi_dpp not needed
                      false);

    double refporositydot = my::structmat_->RefPorosityTimeDeriv();

    //---------------------------  dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx at gausspoint
    static LINALG::Matrix<my::nsd_,1> gradJ(false);
    // spatial porosity gradient
    static LINALG::Matrix<my::nsd_,1>             grad_porosity(false);
    //--------------------- linearization of porosity w.r.t. structure displacements
    static LINALG::Matrix<1,my::nsd_*my::nen_> dphi_dus(false);

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J * N_x
    static LINALG::Matrix<1,my::nsd_*my::nen_> dJ_dus(false);
    //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) + d\phi/(dp*du_s) * dp/dx
    static LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> dgradphi_dus(false);

    //------------------------------------ build F^-T as vector 9x1
    static LINALG::Matrix<my::nsd_*my::nsd_,1> defgrd_IT_vec(false);
    for (int i=0; i<my::nsd_; i++)
      for (int j=0; j<my::nsd_; j++)
        defgrd_IT_vec(i*my::nsd_+j) = defgrd_inv(j,i);

    // dF/dx
    static LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_> F_x(false);

    // dF/dX
    static LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_> F_X(false);

    my::ComputeFDerivative( edispnp,
                        defgrd_inv,
                        F_x,
                        F_X);

    //compute gradients if needed
    my::ComputeGradients(
                          my::J_,
                          dphi_dp,
                          dphi_dJ,
                          defgrd_IT_vec,
                          F_x,
                          my::gradp_,
                          eporositynp,
                          gradJ,
                          my::grad_porosity_,
                          my::refgrad_porosity_);

    ComputeLinearizationOD(
                            dphi_dJ,
                            dphi_dJJ,
                            dphi_dJdp,
                            defgrd_inv,
                            defgrd_IT_vec,
                            F_x,
                            F_X,
                            gradJ,
                            dJ_dus,
                            dphi_dus,
                            dgradphi_dus);

    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    my::GetMaterialParamters(material);

    // set viscous term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    my::visc_old_.Clear();
    my::viscs2_.Clear();
    // compute viscous term from previous iteration and viscous operator
    if (my::is_higher_order_ele_) my::CalcDivEps(evelaf);

    my::ComputeSpatialReactionTerms(material,defgrd_inv);

    //compute linearization of spatial reaction tensor w.r.t. structural displacements
    my::ComputeLinSpatialReactionTerms(material,defgrd_inv,&dJ_dus,NULL);

    // get stabilization parameters at integration point
    my::ComputeStabilizationParameters(vol);

    // compute old RHS of momentum equation and subgrid scale velocity
    my::ComputeOldRHSAndSubgridScaleVelocity();

    // compute old RHS of continuity equation
    my::ComputeOldRHSConti(dphi_dp);

    // compute strong residual of mixture (structural) equation
    if (  my::porofldpara_->StabBiot() and
          (not my::porofldpara_->IsStationaryConti()) and
          my::structmat_->PoroLawType() != INPAR::MAT::m_poro_law_constant
       )
      my::ComputeMixtureStrongResidual(params,defgrd,edispnp,edispn,F_X,true);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = my::fldparatimint_->TimeFac() * my::fac_;
    const double timefacfacpre = my::fldparatimint_->TimeFacPre() * my::fac_;

    //***********************************************************************************************
    // 1) coupling terms in momentum balance

    my::FillMatrixMomentumOD(
                          timefacfac,
                          evelaf,
                          egridv,
                          epreaf,
                          dgradphi_dus,
                          dphi_dp,
                          dphi_dJ,
                          dphi_dus,
                          refporositydot,
                          lin_resM_Dus,
                          lin_resM_Dus_gridvel,
                          ecoupl_u);

    //*************************************************************************************************************
    // 2) coupling terms in continuity equation

    my::FillMatrixContiOD(  timefacfacpre,
                        dphi_dp,
                        dphi_dJ,
                        dphi_dJJ,
                        dphi_dJdp,
                        refporositydot,
                        dgradphi_dus,
                        dphi_dus,
                        dJ_dus,
                        egridv,
                        lin_resM_Dus,
                        lin_resM_Dus_gridvel,
                        ecoupl_p);

    //*************************************************************************************************************
    // 3) additionale terms due to p1 approach (derivatives w.r.t. porosity)
    // 3.1) Momentum equation

    /*  reaction */
    /*
      /                           \
     |                             |
  -  |    sigma * v_f D(phi), v    |
     |                             |
      \                           /
     */
    {
      const double porosity_inv = 1.0/my::porosity_;
      for (int ui=0; ui<my::nen_; ++ui)
      {

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          lin_resM_Dphi(idim,ui) += timefacfac * porosity_inv * my::reavel_(idim) * my::funct_(ui);
        } // end for (idim)
      } // ui
    }


    if (not my::porofldpara_->IsStationaryMomentum())
    //transient terms
    /*  reaction  */
    /*
      /                           \
     |                             |
  -  |    sigma * v_s D(phi), v    |
     |                             |
      \                           /
     */
    {
      const double porosity_inv = 1.0/my::porosity_;
      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          lin_resM_Dphi(idim,ui) += timefacfac * porosity_inv * (-my::reagridvel_(idim)) * my::funct_(ui);
        } // end for (idim)
      } // ui
    }

    //viscous terms (brinkman terms)
    if(my::visceff_)
    {
      static LINALG::Matrix<my::nsd_,my::nsd_> viscstress(false);
      for (int jdim = 0; jdim < my::nsd_; ++jdim)
      {
        for (int idim = 0; idim < my::nsd_; ++idim)
        {
          viscstress(idim,jdim)=my::visceff_*(my::vderxy_(jdim,idim)+my::vderxy_(idim,jdim));
        }
      }

      static LINALG::Matrix<my::nsd_,1> viscstress_gradphi(false);
      viscstress_gradphi.Multiply(viscstress,my::grad_porosity_);

      static LINALG::Matrix<my::nsd_,my::nen_> viscstress_derxy(false);
      viscstress_derxy.Multiply(viscstress,my::derxy_);

      const double porosity_inv = 1.0/my::porosity_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          lin_resM_Dphi(idim,ui) += timefacfac * porosity_inv * ( porosity_inv*viscstress_gradphi(idim)*
                                                                  my::funct_(ui)
                                                                - viscstress_derxy(idim,ui)
                                                              )
                                      ;
        }
      }
    }

    for (int ui=0; ui<my::nen_; ++ui)
    {
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = my::nsd_*vi;
        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          ecouplp1_u(fvi+idim,ui) += my::funct_(vi) * lin_resM_Dphi(idim,ui);
        } // end for (idim)
      } //vi
    } // ui

    //*************************************************************************************************************
    //3.2) Continuity equation

    //transient terms
    /* time derivative*/
    /*
      /                           \
     |                             |
  -  |    D(phi), v                |
     |                             |
      \                           /
     */
    {
      for (int ui=0; ui<my::nen_; ++ui)
        for (int vi=0; vi<my::nen_; ++vi)
          ecouplp1_p(vi,ui) +=   my::fac_ * my::funct_(vi) * my::funct_(ui);
    }

    static LINALG::Matrix<my::nen_,1>    derxy_convel(false);
    derxy_convel.Clear();

    for (int i =0; i< my::nen_; i++)
      for (int j =0; j< my::nsd_; j++)
        derxy_convel(i) += my::derxy_(j,i) * my::velint_(j);

    if (not my::porofldpara_->IsStationaryConti())
    {
      for (int i =0; i< my::nen_; i++)
        for (int j =0; j< my::nsd_; j++)
          derxy_convel(i) += my::derxy_(j,i) * (-my::gridvelint_(j));
    }

    if( static_cast<DRT::ELEMENTS::FluidEleParameterPoro*>(my::fldpara_)->PoroContiPartInt() == false )
    {
      /*
        /                           \     /                             \
       |                             |    |                              |
       |    \nabla v_f D(phi), v     | +  |  (v_f-v_s) \nabla  D(phi), v |
       |                             |    |                              |
        \                           /     \                             /
       */
      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          ecouplp1_p(vi,ui) +=
                               + timefacfacpre * my::vdiv_ * my::funct_(vi) * my::funct_(ui)
                               + timefacfacpre * my::funct_(vi) * derxy_convel(ui)
                                 ;
        }
      }
    }
    else //my::fldpara_->PoroContiPartInt() == true
    {
      /*
          /                             \
          |                              |
       -  |  (v_f-v_s) \nabla  D(phi), v |
          |                              |
          \                             /
       */
      for (int ui=0; ui<my::nen_; ++ui)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          ecouplp1_p(vi,ui) += -1.0 * timefacfacpre * derxy_convel(vi) * my::funct_(ui)
                                 ;
        }
      }
      /*
          /                             \
          |                              |
          |  \nabla v_s D(phi), v        |
          |                              |
          \                             /
       */
      if (not my::porofldpara_->IsStationaryConti())
      {
        for (int ui=0; ui<my::nen_; ++ui)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            ecouplp1_p(vi,ui) += timefacfacpre * my::funct_(vi) * my::gridvdiv_ * my::funct_(ui)
                                   ;
          }
        }
      }
    }

    //*************************************************************************************************************
    // PSPG
    if(my::fldpara_->PSPG())
    {
      double scal_grad_q=0.0;

      if(my::fldpara_->Tds()==INPAR::FLUID::subscales_quasistatic)
      {
        scal_grad_q=my::tau_(1);
      }
      else
      {
        scal_grad_q= 0.0; //my::fldpara_->AlphaF()*fac3;
      }

      {
        const double v1 = -timefacfacpre * my::dtaudphi_(1)/scal_grad_q;
        for (int ui=0; ui<my::nen_; ++ui)
        {
          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            const double v= v1 * my::sgvelint_(idim) * my::funct_(ui);

            for (int vi=0; vi<my::nen_; ++vi)
            {
              ecouplp1_p(vi,ui) += v * my::derxy_(idim, vi);
            } // vi
          } // end for(idim)
        }  // ui
      }

      //linearization of residual in stabilization term w.r.t. porosity
      if (my::is_higher_order_ele_ || my::fldpara_->IsNewton())
      {
        static LINALG::Matrix<my::nen_ , my::nen_ >   temp(false);
        temp.Clear();

        for (int vi=0; vi<my::nen_; ++vi)
          for (int ui=0; ui<my::nen_; ++ui)
            for (int idim=0;idim<my::nsd_;++idim)
                temp(vi,ui) += my::derxy_(idim,vi)*lin_resM_Dphi(idim,ui);

        for (int ui=0; ui<my::nen_; ++ui)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            ecouplp1_p(vi,ui) +=  scal_grad_q * temp(vi,ui);
          } // vi
        } // ui
      } // end if (is_higher_order_ele_) or (newton_)
    }

    //*************************************************************************************************************
    if(my::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
    {
      double reac_tau;
      if (my::fldpara_->Tds()==INPAR::FLUID::subscales_quasistatic)
        reac_tau = my::fldpara_->ViscReaStabFac()*my::reacoeff_*my::tau_(1);
      else
      {
        dserror("Is this factor correct? Check for bugs!");
        reac_tau=0.0;
        //reac_tau = my::fldpara_->ViscReaStabFac()*my::reacoeff_*my::fldpara_->AlphaF()*fac3;
      }

      if (my::is_higher_order_ele_ or my::fldpara_->IsNewton())
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = reac_tau*my::funct_(vi);

          for (int idim=0;idim<my::nsd_;++idim)
          {
            const int fvi_p_idim = my::nsd_*vi+idim;

            for (int ui=0; ui<my::nen_; ++ui)
            {
              ecouplp1_u(fvi_p_idim,ui) += v*lin_resM_Dphi(idim,ui);
            } // vi
          } // ui
        } //idim
      } // end if (is_higher_order_ele_) or (newton_)

      {//linearization of stabilization parameter w.r.t. porosity
        const double v = timefacfac * my::fldpara_->ViscReaStabFac() * ( my::reacoeff_*my::dtaudphi_(1)/my::tau_(1)
        + my::reacoeff_/my::porosity_
        );
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double w = -1.0 * v * my::funct_(vi) ;

          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            const double w_sgvelint= w*my::sgvelint_(idim);
            const int fvi = my::nsd_*vi + idim;

            for (int ui=0; ui<my::nen_; ++ui)
            {
              ecouplp1_u(fvi,ui) += w_sgvelint*my::funct_(ui);
            }
          }
        }  // end for(idim)
      }
    }

  }//loop over gausspoints
}


/*----------------------------------------------------------------------------------------*
 *   Computation of Linearization of porosity gradient w.r.t. pressure         vuong 06/11 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputeLinearization(
                                  const double&                                      dphi_dp,
                                  const double&                                      dphi_dpp,
                                  const double&                                      dphi_dJdp,
                                  const LINALG::Matrix<my::nsd_,1>&                  gradJ,
                                  LINALG::Matrix<my::nsd_,my::nen_>&                 dgradphi_dp)
{
  //porosity is a primary variable -> d(grad(phi)/d(pressure)) is zero!
    dgradphi_dp.Clear();
}

/*----------------------------------------------------------------------------------------*
 *  Computation of Linearization of porosity gradient w.r.t. displacements    vuong 06/11 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputeLinearizationOD(
                            const double&                                      dphi_dJ,
                            const double&                                      dphi_dJJ,
                            const double&                                      dphi_dJp,
                            const LINALG::Matrix<my::nsd_,my::nsd_>&           defgrd_inv,
                            const LINALG::Matrix<my::nsd_*my::nsd_,1>&         defgrd_IT_vec,
                            const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&  F_x,
                            const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&  F_X,
                            const LINALG::Matrix<my::nsd_,1>&                  gradJ,
                            LINALG::Matrix<1,my::nsd_*my::nen_>&               dJ_dus,
                            LINALG::Matrix<1,my::nsd_*my::nen_>&               dphi_dus,
                            LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&        dgradphi_dus)
{
  //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J * N_x
  for (int i=0; i<my::nen_; i++)
    for (int j=0; j<my::nsd_; j++)
      dJ_dus(j+i*my::nsd_)=my::J_*my::derxy_(j,i);

  //--------------------- linearization of porosity w.r.t. structure displacements
  //porosity is a primary variable -> d(grad(phi)/d(displacement)) is zero!
  dphi_dus.Clear();

  //------------------ d( grad(\phi) ) / du_s is also zero
  dgradphi_dus.Clear();

  if( (static_cast<DRT::ELEMENTS::FluidEleParameterPoro*>(my::fldpara_)->PoroContiPartInt() == false) or my::visceff_)
  {
    //---------------------d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x

    //dF^-T/dus : dF/dx = - (F^-1. dN/dx . u_s)^T  : dF/dx
    static LINALG::Matrix<my::nsd_,my::nsd_*my::nen_> dFinvdus_dFdx(false);
    dFinvdus_dFdx.Clear();
    for (int i=0; i<my::nsd_; i++)
      for (int n =0; n<my::nen_; n++)
        for (int j=0; j<my::nsd_; j++)
        {
          const int gid = my::nsd_ * n +j;
          const double defgrd_inv_ij = defgrd_inv(i,j);
          for (int k=0; k<my::nsd_; k++)
          {
            const double derxy_kn = my::derxy_(k,n);
            for (int p=0; p<my::nsd_; p++)
              dFinvdus_dFdx(p, gid) += -defgrd_inv_ij * derxy_kn * F_x(k*my::nsd_+i,p);
          }
        }

    //F^-T : d(dF/dx)/dus =  F^-T : (N,XX * F^ -1 + dF/dX * F^-1 * N,x)
    static LINALG::Matrix<my::nsd_,my::nsd_*my::nen_>        FinvT_dFx_dus(false);
    FinvT_dFx_dus.Clear();

    for (int n =0; n<my::nen_; n++)
      for (int j=0; j<my::nsd_; j++)
      {
        const int gid = my::nsd_ * n +j;
        for (int p=0; p<my::nsd_; p++)
        {
          double val = 0.0;
          const double derxy_p_n = my::derxy_(p,n);
          for (int k=0; k<my::nsd_; k++)
          {
            const double defgrd_inv_kj = defgrd_inv(k,j);
            const double defgrd_inv_kp = defgrd_inv(k,p);
            for (int i=0; i<my::nsd_; i++)
            {
              val +=   defgrd_inv(i,j) * my::N_XYZ2full_(i*my::nsd_+k,n) * defgrd_inv_kp ;
              for (int l=0; l<my::nsd_; l++)
                val += - defgrd_inv(i,l) * F_X(i*my::nsd_+l,k) * defgrd_inv_kj * derxy_p_n ;
            }
          }
          FinvT_dFx_dus(p, gid) += val;
        }
      }

    //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
    static LINALG::Matrix<1,my::nsd_> temp;
    temp.MultiplyTN( defgrd_IT_vec, F_x);

    //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
    static LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> dgradJ_dus;

    dgradJ_dus.MultiplyTN(temp,dJ_dus);

    dgradJ_dus.Update(my::J_,dFinvdus_dFdx,1.0);

    dgradJ_dus.Update(my::J_,FinvT_dFx_dus,1.0);

  }

  return;
}

/*----------------------------------------------------------------------*
 *   PSPG contributions                                     vuong 06/11 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::PSPG(
    LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &         estif_q_u,
    LINALG::Matrix<my::nen_,my::nen_> &                   ppmat,
    LINALG::Matrix<my::nen_,1> &                          preforce,
    const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &    lin_resM_Du,
    const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&     lin_resMRea_Du,
    const LINALG::Matrix<my::nsd_,my::nen_> &             lin_resM_Dp,
    const double &                                        dphi_dp,
    const double &                                        fac3,
    const double &                                        timefacfac,
    const double &                                        timefacfacpre,
    const double &                                        rhsfac)
{
  my::PSPG( estif_q_u,
            ppmat,
            preforce,
            lin_resM_Du,
            lin_resMRea_Du,
            lin_resM_Dp,
            dphi_dp,
            fac3,
            timefacfac,
            timefacfacpre,
            rhsfac);
}

/*----------------------------------------------------------------------*
 *     contributions of reactive stabilization              vuong 06/11 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ReacStab(
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
}

/*--------------------------------------------------------------------------*
 *  compute fluid volume                                         vuong 06/11 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputeVolume(
    Teuchos::ParameterList&         params,
    DRT::ELEMENTS::Fluid*           ele,
    DRT::Discretization&            discretization,
    std::vector<int>&               lm,
    Epetra_SerialDenseVector&       elevec1)
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);
  // set element id
  my::eid_ = ele->Id();

  LINALG::Matrix<my::nsd_,my::nen_> edispnp(true);
  LINALG::Matrix<my::nen_, 1> eporositynp(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &edispnp, &eporositynp,"dispnp");

  LINALG::Matrix<my::nsd_,my::nen_> egridvnp(true);
  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &egridvnp, NULL,"gridv");

  // get new node positions of ALE mesh
  my::xyze_ += edispnp;

  // integration loop
  for ( DRT::UTILS::GaussIntegration::iterator iquad=my::intpoints_.begin(); iquad!=my::intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    //-----------------------------------computing the porosity
    my::porosity_=my::funct_.Dot(eporositynp);

    // get structure velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    LINALG::Matrix<my::nsd_,my::nsd_>              gridvelderxy;
    gridvelderxy.MultiplyNT(egridvnp,my::derxy_);

    my::gridvdiv_ = 0.0;
    for (int idim = 0; idim <my::nsd_; ++idim)
      my::gridvdiv_ += gridvelderxy(idim,idim);

    elevec1(0) += my::porosity_* my::fac_;
  } // end of integration loop

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1<DRT::Element::nurbs27>;
