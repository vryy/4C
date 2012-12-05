/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_poro.cpp

\brief Internal implementation of Fluid element

<pre>
Maintainer: Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_poro.H"
#include "fluid_ele.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_utils.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_fem_general/drt_utils_gder2.H"

#include "../drt_geometry/position_array.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/standardtypes_cpp.H"

//#include "Sacado.hpp"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoro<distype> * DRT::ELEMENTS::FluidEleCalcPoro<distype>::Instance( bool create )
{
  static FluidEleCalcPoro<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcPoro<distype>();
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
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoro<distype>::FluidEleCalcPoro()
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc(),
    visceff_(0.0)
{

}


/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::Evaluate(DRT::ELEMENTS::Fluid*    ele,
                                                 DRT::Discretization &          discretization,
                                                 const std::vector<int> &       lm,
                                                 Teuchos::ParameterList&        params,
                                                 Teuchos::RCP<MAT::Material> &  mat,
                                                 Epetra_SerialDenseMatrix&      elemat1_epetra,
                                                 Epetra_SerialDenseMatrix&      elemat2_epetra,
                                                 Epetra_SerialDenseVector&      elevec1_epetra,
                                                 Epetra_SerialDenseVector&      elevec2_epetra,
                                                 Epetra_SerialDenseVector&      elevec3_epetra,
                                                 bool                           offdiag)
{
  if (not offdiag)
    return Evaluate(  ele,
                      discretization,
                      lm,
                      params,
                      mat,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra,
                      my::intpoints_);
  else
    return EvaluateOD(  ele,
                        discretization,
                        lm,
                        params,
                        mat,
                        elemat1_epetra,
                        elemat2_epetra,
                        elevec1_epetra,
                        elevec2_epetra,
                        elevec3_epetra,
                        my::intpoints_);
}


/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (2)
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::Evaluate(
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
  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  my::rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat1(elemat1_epetra,true);
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> elemat2(elemat2_epetra,true);
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
  this->BodyForce(ele,my::fldpara_,ebofoaf,eprescpgaf,escabofoaf);

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
  if (my::fldpara_->IsGenalphaNP())
    my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelnp,
        NULL, "velnp");

  LINALG::Matrix<my::nsd_, my::nen_> emhist(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &emhist,
      NULL, "hist");

  LINALG::Matrix<my::nsd_, my::nen_> eaccam(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eaccam,
      NULL, "accam");

  LINALG::Matrix<my::nsd_, my::nen_> eveln(true);
  LINALG::Matrix<my::nen_, 1> epren(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln,
      &epren, "veln");

  LINALG::Matrix<my::nen_,1> epressn_timederiv(true);
  my::ExtractValuesFromGlobalVector(discretization,lm,*my::rotsymmpbc_,NULL,&epressn_timederiv,"accn");

  LINALG::Matrix<my::nen_, 1> epressnp_timederiv(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
      &epressnp_timederiv, "accnp");

  if (not my::fldpara_->IsGenalpha())
    eaccam.Clear();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);
  LINALG::Matrix<my::nsd_, my::nen_> edispn(true);

  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispnp,
      NULL, "dispnp");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridv,
      NULL, "gridv");
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispn,
      NULL, "dispn");

  //ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, NULL, &initporosity_, "initporosity");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> >(
      ele, my::xyze_);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = Evaluate(ele->Id(), params, ebofoaf, elemat1, elemat2,
      elevec1, evelaf, epreaf, evelnp, emhist, epren, epressn_timederiv, epressnp_timederiv,
      eaccam, edispnp, edispn, egridv, mat, ele->IsAle(), intpoints);

  return result;
}


/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (2)
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateOD(
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
  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  my::rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nsd_ * my::nen_> elemat1(elemat1_epetra, true);
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
  my::BodyForce(ele,my::fldpara_,ebofoaf,eprescpgaf,escabofoaf);

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
  if (my::fldpara_->IsGenalphaNP())
    this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evelnp,
        NULL, "velnp");

  LINALG::Matrix<my::nsd_, my::nen_> emhist(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &emhist,
      NULL, "hist");

  LINALG::Matrix<my::nsd_, my::nen_> eaccam(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eaccam,
      NULL, "accam");

  LINALG::Matrix<my::nsd_, my::nen_> eveln(true);
  LINALG::Matrix<my::nen_, 1> epren(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln,
      &epren, "veln");

  LINALG::Matrix<my::nen_, 1> epressnp_timederiv(true);
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, NULL,
      &epressnp_timederiv, "accnp");

  if (not my::fldpara_->IsGenalpha())
    eaccam.Clear();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);
  LINALG::Matrix<my::nsd_, my::nen_> edispn(true);

  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispnp,
      NULL, "dispnp");
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &egridv,
      NULL, "gridv");
  this->ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &edispn,
      NULL, "dispn");

  //ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &initporosity_, "initporosity");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> >(
      ele, my::xyze_);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = EvaluateOD(ele->Id(),
      params,
      ebofoaf,
      elemat1,
      elevec1,
      evelaf,
      epreaf,
      evelnp,
      emhist,
      epren,
      epressnp_timederiv,
      eaccam,
      edispnp,
      edispn,
      egridv,
      mat,
      ele->IsAle(),
      intpoints);

  return result;
}


/*----------------------------------------------------------------------*
 * evaluation of system matrix and residual for porous flow (3)
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::Evaluate(
  int                                                             eid,
  Teuchos::ParameterList&                                         params,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       ebofoaf,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat1,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat2,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,1> &                       elevec1,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelaf,
  const LINALG::Matrix<my::nen_,1>    &                           epreaf,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelnp,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
  const LINALG::Matrix<my::nen_,1>    &                           epren,
  const LINALG::Matrix<my::nen_,1>    &                           epressn_timederiv,
  const LINALG::Matrix<my::nen_,1>    &                           epressnp_timederiv,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       eaccam,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispnp,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispn,
  const LINALG::Matrix<my::nsd_,my::nen_> &                       egridv,
  Teuchos::RCP<MAT::Material>                                     mat,
  bool                                                            isale,
  const DRT::UTILS::GaussIntegration &                            intpoints)
{
  // flag for higher order elements
  my::is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (my::fldpara_->IsInconsistent() == true) my::is_higher_order_ele_ = false;

  // stationary formulation does not support ALE formulation
  if (isale and my::fldpara_->IsStationary())
    dserror("No ALE support within stationary fluid solver.");

  //initporosity_   = params.get<double>("initporosity");

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(eid,
       ebofoaf,
       evelaf,
       evelnp,
       epreaf,
       eaccam,
       emhist,
       epren,
       epressn_timederiv,
       epressnp_timederiv,
       edispnp,
       edispn,
       egridv,
       elemat1,
       elemat2,  // -> emesh
       elevec1,
       mat,
       isale,
       intpoints);

  return 0;
}




/*----------------------------------------------------------------------*
 * evaluation of coupling terms for porous flow (3)
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoro<distype>::EvaluateOD(
    int eid,
    Teuchos::ParameterList&                                           params,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        ebofoaf,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nsd_ * my::nen_> &  elemat1,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
    const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
    const LINALG::Matrix<my::nen_, 1> &                               epreaf,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        emhist,
    const LINALG::Matrix<my::nen_, 1> &                               epren,
    const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        eaccam,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        edispn,
    const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
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
  if (isale and my::fldpara_->IsStationary())
    dserror("No ALE support within stationary fluid solver.");

    // ---------------------------------------------------------------------
    // call routine for calculating element matrix and right hand side
    // ---------------------------------------------------------------------
    SysmatOD(eid,
        ebofoaf,
        evelaf,
        evelnp,
        epreaf,
        eaccam,
        emhist,
        epren,
        epressnp_timederiv,
        edispnp,
        edispn,
        egridv,
        elemat1,
        elevec1,
        mat,
        isale,
        intpoints);

    return 0;
  }


/*----------------------------------------------------------------------*
 |  calculate element matrix and rhs for porous flow         vuong 06/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::Sysmat(
  int                                                           eid,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      ebofoaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      evelaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      evelnp,
  const LINALG::Matrix<my::nen_,1>&                             epreaf,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      eaccam,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      emhist,
  const LINALG::Matrix<my::nen_,1>    &                         epren,
  const LINALG::Matrix<my::nen_,1>    &                         epressn_timederiv,
  const LINALG::Matrix<my::nen_,1>    &                         epressnp_timederiv,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      edispnp,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      edispn,
  const LINALG::Matrix<my::nsd_,my::nen_>&                      egridv,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  emesh,
  LINALG::Matrix<(my::nsd_+1)*my::nen_,1>&                      eforce,
  Teuchos::RCP<const MAT::Material>                             material,
  bool                                                          isale,
  const DRT::UTILS::GaussIntegration &                          intpoints
  )
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  LINALG::Matrix<my::nen_*my::nsd_,my::nen_*my::nsd_>  estif_u(true);
  LINALG::Matrix<my::nen_*my::nsd_,my::nen_>          estif_p_v(true);
  LINALG::Matrix<my::nen_, my::nen_*my::nsd_>         estif_q_u(true);
  LINALG::Matrix<my::nen_,my::nen_>                   ppmat(true);

  // definition of vectors
  LINALG::Matrix<my::nen_,1>     preforce(true);
  LINALG::Matrix<my::nsd_,my::nen_>  velforce(true);

  // definition of velocity-based momentum residual vectors
  LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>  lin_resM_Du(true);
  LINALG::Matrix<my::nsd_,1>          resM_Du(true);

  //material coordinates xyze0
  LINALG::Matrix<my::nsd_,my::nen_> xyze0= my::xyze_;

  // add displacement when fluid nodes move in the ALE case
  //if (isale)
  my::xyze_ += edispnp;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter(eid);

  // set element area or volume
  const double vol = my::fac_;

  //access structure discretization
  RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding structure element (it has the same global ID as the fluid element)
  DRT::Element* structele = structdis->gElement(eid);
  if(structele == NULL)
    dserror("Fluid element %i not on local processor", eid);
  //get fluid material
  MAT::StructPoro* structmat = static_cast<MAT::StructPoro*>((structele->Material()).get());
  if(structmat->MaterialType() != INPAR::MAT::m_structporo)
    dserror("invalid structure material for poroelasticity");

  // get material parameters at element center
  if (not my::fldpara_->MatGp() or not my::fldpara_->TauGp())
  {
    const MAT::FluidPoro* actmat = static_cast<const MAT::FluidPoro*>(material.get());
    if(actmat->MaterialType() != INPAR::MAT::m_fluidporo)
     dserror("invalid fluid material for poroelasticity");

    // set density at n+alpha_F/n+1 and n+alpha_M/n+1
    my::densaf_ = actmat->Density();
    my::densam_ = my::densaf_;

    // calculate reaction coefficient
    my::reacoeff_ = actmat->ComputeReactionCoeff();

    visceff_       = actmat->EffectiveViscosity();
  }

  // calculate stabilization parameters at element center
  if (not my::fldpara_->TauGp())
  {
    // check stabilization parameter definition for porous flow
    if (not (my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
             my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt))
      dserror("incorrect definition of stabilization parameter for porous flow");

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient
    double sigma_tot = my::reacoeff_;
    if (my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina)
      sigma_tot += 1.0/my::fldpara_->TimeFac();

    // calculate characteristic element length
    double strle  = 0.0;
    double hk     = 0.0;
    my::CalcCharEleLength(vol,0.0,strle,hk);

    // constants c_u and c_p as suggested in Badia and Codina (2010), method A
    const double c_u = 4.0;
    const double c_p = 4.0;

    // tau_Mu not required for porous flow
    my::tau_(0) = 0.0;
    my::tau_(1) = 1.0/(c_u*my::densaf_*sigma_tot);
    my::tau_(2) = c_p*DSQR(hk)*my::reacoeff_;
  }

  std::vector<double> dporodt_gp(intpoints.NumPoints(),0.0); // urrecha
  std::vector<double> rhscon(intpoints.NumPoints(),0.0);
  std::vector<double> J_gp(intpoints.NumPoints(),0.0);

  std::vector<LINALG::Matrix<my::nsd_,1> > gradporosity_gp (intpoints.NumPoints());
  for(int i = 0 ; i < intpoints.NumPoints() ; i++)
  {
    for(int j=0; j<my::nsd_; j++)
      gradporosity_gp.at(i)(j) = 0.0;
  }

  std::vector<LINALG::Matrix<1,my::nsd_> > gradJ_gp (intpoints.NumPoints());
  for(int i = 0 ; i < intpoints.NumPoints() ; i++)
  {
    for(int j=0; j<my::nsd_; j++)
      gradJ_gp.at(i)(j) = 0.0;
  }

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    if(my::is_higher_order_ele_)
    {
      // get the second derivatives of standard element at current GP
      DRT::UTILS::shape_function_deriv2<distype>(my::xsi_,my::deriv2_);

      // get the second derivatives of standard element at current GP w.r.t. xyz
      DRT::UTILS::gder2<distype>(my::xjm_,my::derxy_,my::deriv2_,my::xyze_,my::derxy2_);
    }
    else
    {
      my::deriv2_.Clear();
      my::derxy2_.Clear();
    }

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including derivatives and grid velocity)
    //  2) pressure (including derivatives)
    //  3) body-force vector
    //  4) "history" vector for momentum equation
    //----------------------------------------------------------------------
    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::velint_.Multiply(evelaf,my::funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::vderxy_.MultiplyNT(evelaf,my::derxy_);

    // get convective velocity at integration point
    // (ALE case handled implicitly here using the (potential
    //  mesh-movement-dependent) convective velocity, avoiding
    //  various ALE terms used to be calculated before)
    //convmy::velint_.Update(my::velint_);
    my::convvelint_.Multiply(-1.0, egridv, my::funct_, 0.0);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = my::funct_.Dot(epreaf);

    // get pressure time derivative at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press_dot = my::funct_.Dot(epressnp_timederiv);

    // get pressure time derivative at integration point
    // (value at n )
    //double pressn_dot = my::funct_.Dot(epressn_timederiv);

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::gradp_.Multiply(my::derxy_,epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::bodyforce_.Multiply(ebofoaf,my::funct_);

    // get momentum history data at integration point
    // (only required for one-step-theta and BDF2 time-integration schemes)
    my::histmom_.Multiply(emhist,my::funct_);

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    LINALG::Matrix<my::nsd_,1>             gridvelint;
    gridvelint.Multiply(egridv,my::funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    LINALG::Matrix<my::nsd_,my::nsd_>              gridvelderxy;
    gridvelderxy.MultiplyNT(egridv,my::derxy_);

    //----------------------------------------------------------------------
    //  evaluation of various partial operators at integration point
    //  1) convective term from previous iteration (mandatorily set to zero)
    //  2) viscous term from previous iteration and viscous operator
    //  3) divergence of velocity from previous iteration
    //----------------------------------------------------------------------
    // set convective term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    //conv_old_.Clear();

    //set old convective term to ALE-Term only
    my::conv_old_.Multiply(my::vderxy_,my::convvelint_);
    my::conv_c_.MultiplyTN(my::derxy_,my::convvelint_);

    // set viscous term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    my::visc_old_.Clear();

    // compute divergence of velocity from previous iteration
    my::vdiv_ = 0.0;

    double gridvdiv=0.0;

    if (not my::fldpara_->IsGenalphaNP())
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        my::vdiv_ += my::vderxy_(idim, idim);
        gridvdiv += gridvelderxy(idim,idim);
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

        gridvdiv += gridvelderxy(idim,idim);
      }
    }

    //------------------------get determinant of Jacobian dX / ds
    // transposed jacobian "dX/ds"
    LINALG::Matrix<my::nsd_,my::nsd_> xjm0(true);
    xjm0.MultiplyNT(my::deriv_,xyze0);

    // inverse of transposed jacobian "ds/dX"
    LINALG::Matrix<my::nsd_,my::nsd_> xji0(true);
    const double  det0= xji0.Invert(xjm0);

    // ----------------------compute derivatives N_XYZ at gp w.r.t. material coordinates
    LINALG::Matrix<my::nsd_,my::nen_>          N_XYZ(false);
    N_XYZ.Multiply(xji0,my::deriv_);

    // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ * N_XYZ^T
    LINALG::Matrix<my::nsd_,my::nsd_>          defgrd(false);
    defgrd.MultiplyNT(my::xyze_,N_XYZ);

    // inverse deformation gradient F^-1
    LINALG::Matrix<my::nsd_,my::nsd_>          defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //------------------------------------ build F^-1 as vector 9x1
    LINALG::Matrix<my::nsd_*my::nsd_,1> defgrd_inv_vec;
    for(int i=0; i<my::nsd_; i++)
      for(int j=0; j<my::nsd_; j++)
        defgrd_inv_vec(i*my::nsd_+j) = defgrd_inv(i,j);

    //------------------------------------ build F^-T as vector 9x1
    LINALG::Matrix<my::nsd_*my::nsd_,1> defgrd_IT_vec;
    for(int i=0; i<my::nsd_; i++)
      for(int j=0; j<my::nsd_; j++)
        defgrd_IT_vec(i*my::nsd_+j) = defgrd_inv(j,i);

    // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    double  J = my::det_/det0;

    //-----------------------------------auxilary variables for computing the porosity
    double dphi_dp=0.0;
    double dphi_dJ=0.0;
    double dphi_dJdp=0.0;
    double dphi_dJJ=0.0;
    double dphi_dpp=0.0;
    double porosity=0.0;

    structmat->ComputePorosity(press, J, *(iquad),porosity,dphi_dp,dphi_dJ,dphi_dJdp,dphi_dJJ,dphi_dpp);;

    //--------------------------- compute dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx at gausspoint
    //LINALG::Matrix<1,my::nsd_> gradJ_2(true);
    //gradJ_2.MultiplyTN(J, defgrd_IT_vec, F_x);

    //ask structure material for the spatial gradient of the jacobian
    LINALG::Matrix<1,my::nsd_> gradJ;
    structmat->GetGradJAtGP(gradJ,*(iquad));

    gradJ_gp[*(iquad)] = gradJ;
    J_gp[*(iquad)] = J;

    //--linearization of porosity gradient w.r.t. pressure at gausspoint
    //d(grad(phi))/dp = dphi/(dJdp)* dJ/dx + d^2phi/(dp)^2 * dp/dx + dphi/dp* N,x
    LINALG::Matrix<my::nsd_,my::nen_>             dgradphi_dp(true);

    // calculate spatial porosity gradient
    LINALG::Matrix<my::nsd_,1>             grad_porosity(true);

    if( (my::fldpara_->PoroContiPartInt() == false) or visceff_)
    {
      LINALG::Matrix<my::nsd_,1>             mat_grad_porosity;
      //-----get material porosity gradient from structure material
      structmat->GetGradPorosityAtGP(mat_grad_porosity,*(iquad));

      grad_porosity.MultiplyTN(defgrd_inv,mat_grad_porosity);

      gradporosity_gp[*(iquad)] = grad_porosity;    // trial urrecha.

      dgradphi_dp.MultiplyTT(dphi_dJdp,gradJ,my::funct_ );
      dgradphi_dp.MultiplyNT(dphi_dpp, my::gradp_,my::funct_,1.0);
      dgradphi_dp.Update(dphi_dp, my::derxy_,1.0);
    }

    dporodt_gp[*iquad]= dphi_dJ*J*gridvdiv  +  dphi_dp*press_dot;    // trial urrecha.

    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    if (my::fldpara_->MatGp())
    {
      const MAT::FluidPoro* actmat = static_cast<const MAT::FluidPoro*>(material.get());

      // set density at n+alpha_F/n+1 and n+alpha_M/n+1
      my::densaf_ = actmat->Density();
      my::densam_ = my::densaf_;

      // calculate reaction coefficient
      my::reacoeff_ = actmat->ComputeReactionCoeff()*porosity;

      visceff_       = actmat->EffectiveViscosity();
    }
    else dserror("Fluid material parameters have to be evaluated at gauss point for porous flow!");

    // calculate stabilization parameters at integration point
    if (my::fldpara_->TauGp())
    {
      // check stabilization parameter definition for porous flow
      if (not (my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
               my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt))
        dserror("incorrect definition of stabilization parameter for porous flow");

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient
      double sigma_tot = my::reacoeff_;
      if (my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina)
        sigma_tot += 1.0/my::fldpara_->TimeFac();

      // calculate characteristic element length
      double strle  = 0.0;
      double hk     = 0.0;
      my::CalcCharEleLength(vol,0.0,strle,hk);

      // constants c_u and c_p as suggested in Badia and Codina (2010), method A
      const double c_u = 4.0;
      const double c_p = 4.0;

      // tau_Mu not required for porous flow
      my::tau_(0) = 0.0;
      my::tau_(1) = 1.0/(c_u*my::densaf_*sigma_tot);
      my::tau_(2) = c_p*DSQR(hk)*my::reacoeff_;
    }
    else dserror("Fluid stabilization parameters have to be evaluated at gauss point for porous flow!");


    //******************* FAD ************************
    /*
    // sacado data type replaces "double"
    typedef Sacado::Fad::DFad<double> FAD;  // for first derivs
    // sacado data type replaces "double" (for first+second derivs)
    typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD;

    LINALG::TMatrix<FAD,my::nen_,1> fad_my::funct_;
    LINALG::TMatrix<FAD,my::nen_,1> fad_epreaf;
    LINALG::TMatrix<FAD,my::nen_,1> fad_epressnp_timederiv;
    LINALG::TMatrix<FAD,my::nsd_,my::nen_> fad_derxy_;
    LINALG::TMatrix<FAD,my::nsd_,1> fad_my::gradp_;
    for (int i=0; i<my::nen_; i++)
    {
      fad_my::funct_(i)=my::funct_(i);
      fad_epreaf(i)=epreaf(i);
      //fad_epressnp_timederiv(i)=epressnp_timederiv(i);
      fad_epreaf(i).diff(i,my::nen_);
      for(int j=0; j<my::nsd_; j++)
        fad_derxy_(j,i)=derxy_(j,i);
    }

    FAD fad_press = fad_my::funct_.Dot(fad_epreaf);
    fad_my::gradp_.Multiply(fad_derxy_,fad_epreaf);

    FAD fad_a     = ( bulkmodulus_/(1-initporosity_) + fad_press - penalty_/initporosity_ ) * J;
    FAD fad_b     = -fad_a + bulkmodulus_ + penalty_;
    FAD fad_c   = (fad_b/fad_a) * (fad_b/fad_a) + 4*penalty_/fad_a;
    FAD fad_d     = sign*sqrt(fad_c)*fad_a;

    FAD fad_porosity = 1/(2*fad_a)*(-fad_b+fad_d);

    FAD fad_d_p   =  J * (-fad_b+2*penalty_)/fad_d;
    FAD fad_d_J   =  fad_a/J * ( -fad_b + 2*penalty_ ) / fad_d;

    LINALG::TMatrix<FAD,1,my::nsd_>             fad_grad_porosity;
    FAD fad_dphi_dp=  - J * fad_porosity/fad_a + (J+fad_d_p)/(2*fad_a);
    FAD fad_dphi_dJ=  -fad_porosity/J+ 1/(2*J) + fad_d_J / (2*fad_a);

    for (int idim=0; idim<my::nsd_; ++idim)
    {
      fad_grad_porosity(idim)=fad_dphi_dp*fad_my::gradp_(idim)+fad_dphi_dJ*gradJ(idim);
    }

    for (int i=0; i<my::nen_; i++)
     for (int j=0; j<my::nsd_; j++)
     {
       if( (dgradphi_dp(j,i)-fad_grad_porosity(j).dx(i)) > 1e-8)
       {
         cout<<"dgradphi_dp("<<j<<","<<i<<"): "<<dgradphi_dp(j,i)<<endl;
         cout<<"fad_grad_porosity.dx("<<i<<"): "<<fad_grad_porosity(j).dx(i)<<endl;
         cout<<"dgradphi_dp:"<<endl<<dgradphi_dp<<endl;
         cout<<"fad_grad_porosity:"<<endl<<fad_grad_porosity<<endl;
         dserror("check dgradphi_dp failed!");
       }
     }
   cout<<"dgradphi_dp check done and ok"<<endl;
   */
    //******************* FAD ************************

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac    = my::fldpara_->TimeFac()    * my::fac_;
    const double timefacfacpre = my::fldpara_->TimeFacPre() * my::fac_;
    const double rhsfac        = my::fldpara_->TimeFacRhs() * my::fac_;

    //----------------------------------------------------------------------
    // computation of various residuals and residual-based values such as
    // the subgrid-scale velocity
    //----------------------------------------------------------------------
    // compute rhs for momentum equation and momentum residual
    // -> different for generalized-alpha and other time-integration schemes
    //GetResidualMomentumEq(eaccam,my::fldpara_->TimeFac());
    if (my::fldpara_->IsGenalpha())
    {
      if (my::fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
        dserror("The combination of generalized-alpha time integration and a Boussinesq approximation has not been implemented yet!");

      // rhs of momentum equation: density*bodyforce at n+alpha_F
      my::rhsmom_.Update(my::densaf_,my::bodyforce_,0.0);

      // get acceleration at time n+alpha_M at integration point
      my::accint_.Multiply(eaccam,my::funct_);

      // evaluate momentum residual once for all stabilization right hand sides
      for (int rr=0;rr<my::nsd_;++rr)
      {
        my::momres_old_(rr) = my::densam_*my::accint_(rr)+my::densaf_*my::conv_old_(rr)+my::gradp_(rr)
                         -2*visceff_*my::visc_old_(rr)+my::reacoeff_*(my::velint_(rr) + my::convvelint_(rr))-my::densaf_*my::bodyforce_(rr);
      }

  //    // modify integration factors for right-hand side such that they
  //    // are identical in case of generalized-alpha time integration:
  //    rhsfac   /= my::fldpara_->AlphaF();
  //    rhsresfac = rhsfac;
    }
    else
    {
      if (not my::fldpara_->IsStationary())
      {
        // rhs of instationary momentum equation:
        // density*theta*bodyforce at n+1 + density*(histmom/dt)
        // in the case of a Boussinesq approximation: f = (rho - rho_0)/rho_0 *g
        // else:                                      f = rho * g
        if (my::fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
          //my::rhsmom_.Update((densn_/my::fldpara_->Dt()),histmom_,deltadens_*my::fldpara_->Theta(),my::bodyforce_);
          my::rhsmom_.Update((my::densn_/my::fldpara_->Dt()/my::fldpara_->Theta()),my::histmom_,my::deltadens_,my::bodyforce_);
        else
          //my::rhsmom_.Update((my::densn_/my::fldpara_->Dt()),my::histmom_,my::densaf_*my::fldpara_->Theta(),my::bodyforce_);
          my::rhsmom_.Update((my::densn_/my::fldpara_->Dt()/my::fldpara_->Theta()),my::histmom_,my::densaf_,my::bodyforce_);

        // compute instationary momentum residual:
        // momres_old = u_(n+1)/dt + theta ( ... ) - my::histmom_/dt - theta*my::bodyforce_
        for (int rr=0;rr<my::nsd_;++rr)
        {
          /*my::momres_old_(rr) = my::densaf_*my::velint_(rr)/my::fldpara_->Dt()
                             +my::fldpara_->Theta()*(my::densaf_*conv_old_(rr)+my::gradp_(rr)
                             -2*visceff_*visc_old_(rr)+my::reacoeff_*my::velint_(rr))-my::rhsmom_(rr);*/
          my::momres_old_(rr) = ((my::densaf_*my::velint_(rr)/my::fldpara_->Dt()
                           +my::fldpara_->Theta()*(my::densaf_*my::conv_old_(rr)+my::gradp_(rr)
                           -2*visceff_*my::visc_old_(rr)+my::reacoeff_*(my::velint_(rr)+my::convvelint_(rr))))/my::fldpara_->Theta())-my::rhsmom_(rr);
       }

  //      // modify residual integration factor for right-hand side in instat. case:
  //      rhsresfac *= my::fldpara_->Dt();
      }
      else
      {
        // rhs of stationary momentum equation: density*bodyforce
        // in the case of a Boussinesq approximation: f = (rho - rho_0)/rho_0 *g
        // else:                                      f = rho * g
        if (my::fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
             my::rhsmom_.Update(my::deltadens_,my::bodyforce_, 0.0);
        else my::rhsmom_.Update(my::densaf_,my::bodyforce_,0.0);

        // compute stationary momentum residual:
        for (int rr=0;rr<my::nsd_;++rr)
        {
          my::momres_old_(rr) = my::densaf_*my::conv_old_(rr)+my::gradp_(rr)-2*visceff_*my::visc_old_(rr)
                           +my::reacoeff_*(my::velint_(rr)+my::convvelint_(rr))-my::rhsmom_(rr);
        }
      }
    }
    //-------------------------------------------------------

    // compute subgrid-scale velocity
    my::sgvelint_.Update(-my::tau_(1),my::momres_old_,0.0);

    // set velocity-based momentum residual vectors to zero
    lin_resM_Du.Clear();
    resM_Du.Clear();

    double grad_porosity_relvelint=0.0;
    for(int i=0; i<my::nsd_; i++)
      grad_porosity_relvelint += grad_porosity(i) * (my::velint_(i)-gridvelint(i));

    my::rhscon_ =0.0;

    // compute residual of continuity equation
    my::conres_old_ = (dphi_dp  * press_dot+ dphi_dJ  * J  * gridvdiv
          + porosity*my::vdiv_+grad_porosity_relvelint)/my::fldpara_->Theta()-my::rhscon_;

    rhscon[*iquad] = dporodt_gp[*iquad] + grad_porosity_relvelint + porosity * my::vdiv_ ;

    // compute first version of velocity-based momentum residual containing
    // inertia and reaction term
    int idim_nsd_p_idim[my::nsd_];
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      idim_nsd_p_idim[idim]=idim*my::nsd_+idim;
    }

    if (my::fldpara_->IsStationary() == false)
    {
      const double fac_densam=my::fac_*my::densam_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const double v=fac_densam*my::funct_(ui);

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
        }
      }
    }

   //reactive part
    const double fac_reac=timefacfac*my::reacoeff_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const double v=fac_reac*my::funct_(ui);

      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }

    //convective ALE-part

    const double timefacfac_densaf=timefacfac*my::densaf_;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const double v=timefacfac_densaf*my::conv_c_(ui);

      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }

    //----------------------------------------------------------------------
    // computation of standard Galerkin and stabilization contributions to
    // element matrix and right-hand-side vector
    //----------------------------------------------------------------------
    // 1) standard Galerkin inertia and reaction terms
    /* inertia (contribution to mass matrix) if not is_stationary */
    /*
            /              \
           |                |
           |    rho*Du , v  |
           |                |
            \              /
    */
    /*  reaction */
    /*
            /                \
           |                  |
           |    sigma*Du , v  |
           |                  |
            \                /
    */
    /* convection, convective ALE part  */
    /*
              /                             \
             |  /        n+1       \          |
             | |   rho*us   o nabla | Du , v  |
             |  \       (i)        /          |
              \                             /
    */

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui   = my::nsd_*ui;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi   = my::nsd_*vi;

        for (int idim = 0; idim <my::nsd_; ++idim)
          estif_u(fvi+idim,fui+idim) += my::funct_(vi)*lin_resM_Du(idim*my::nsd_+idim,ui);
      } //vi
    } // ui

    // inertia terms on the right hand side for instationary fluids
    if (not my::fldpara_->IsStationary())
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        if (my::fldpara_->IsGenalpha()) resM_Du(idim)+=rhsfac*my::densam_*my::accint_(idim);
        else                            resM_Du(idim)+=my::fac_*my::densaf_*my::velint_(idim);
      }

      //coupling part RHS
      // reacoeff * phi * v_s
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for(int idim = 0; idim <my::nsd_; ++idim)
          velforce(idim,vi) -= -rhsfac* my::funct_(vi) * my::reacoeff_ * gridvelint(idim) ;
      }
    }  // end if (not stationary)

    // convective ALE-part
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      resM_Du(idim)+=rhsfac*my::densaf_*my::conv_old_(idim);
    }  // end for(idim)

    // reactive part
    double rhsfac_rea =rhsfac*my::reacoeff_;
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      resM_Du(idim) += rhsfac_rea*my::velint_(idim);
    }

    for (int vi=0; vi<my::nen_; ++vi)
    {
      for(int idim = 0; idim <my::nsd_; ++idim)
      {
        velforce(idim,vi)-=resM_Du(idim)*my::funct_(vi);
      }
    }


/************************************************************************/
    /* Brinkman term: viscosity term */
    /*
                     /                        \
                    |       /  \         / \   |
              2 mu  |  eps | Du | , eps | v |  |
                    |       \  /         \ /   |
                     \                        /
    */

    if(visceff_)
    {
      LINALG::Matrix<my::nsd_,my::nsd_> viscstress(true);
      const double visceff_timefacfac = visceff_*timefacfac;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi   = my::nsd_*vi;
        const double temp2=visceff_timefacfac*my::funct_(vi)/porosity;

        for (int jdim= 0; jdim<my::nsd_;++jdim)
        {
          const double temp=visceff_timefacfac*my::derxy_(jdim,vi);

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui   = my::nsd_*ui;

            for (int idim = 0; idim <my::nsd_; ++idim)
            {
              const int fvi_p_idim = fvi+idim;

              estif_u(fvi_p_idim,fui+jdim) += temp*my::derxy_(idim, ui)
                                              + temp2*(my::derxy_(idim, ui)*grad_porosity(jdim));
            } // end for (jdim)
          } // end for (idim)
        } // ui
      } //vi

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi   = my::nsd_*vi;
        const double temp2=visceff_timefacfac*my::funct_(vi)/porosity;

        for (int jdim= 0; jdim<my::nsd_;++jdim)
        {
          const double temp=visceff_timefacfac*my::derxy_(jdim,vi);

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui   = my::nsd_*ui;

            for (int idim = 0; idim <my::nsd_; ++idim)
            {
              const int fvi_p_idim = fvi+idim;

              estif_u(fvi_p_idim,fui+idim) += temp*my::derxy_(jdim, ui)
                                              + temp2*(my::derxy_(jdim, ui)*grad_porosity(jdim));
            } // end for (jdim)
          } // end for (idim)
        } // ui
      } //vi

      //const double v = visceff_*rhsfac;

      for (int jdim = 0; jdim < my::nsd_; ++jdim)
      {
        for (int idim = 0; idim < my::nsd_; ++idim)
        {
          viscstress(idim,jdim)=visceff_*(my::vderxy_(jdim,idim)+my::vderxy_(idim,jdim));
        }
      }

      LINALG::Matrix<my::nsd_,1> viscstress_gradphi(true);
      viscstress_gradphi.Multiply(viscstress,grad_porosity);

      // computation of right-hand-side viscosity term
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for (int idim = 0; idim < my::nsd_; ++idim)
        {
          for (int jdim = 0; jdim < my::nsd_; ++jdim)
          {
            /* viscosity term on right-hand side */
            velforce(idim,vi)-= rhsfac*viscstress(idim,jdim)*my::derxy_(jdim,vi)+viscstress_gradphi(idim)*my::funct_(vi);
          }
        }
      }

      LINALG::Matrix<my::nsd_,my::nen_> viscstress_dgradphidp(true);
      viscstress_dgradphidp.Multiply(viscstress,dgradphi_dp);
      for (int ui=0; ui<my::nen_; ++ui)
      {
        const double v = timefacfacpre*my::funct_(ui);
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = my::nsd_*vi;
          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            estif_p_v(fvi + idim,ui) += v * 1/(porosity) * ( 1/porosity*viscstress_gradphi(idim)*
                                                            dphi_dp* my::funct_(vi)
                                                            - viscstress_dgradphidp(idim,ui)
                                                           )
                                        ;
          }
        }
      }
    }

/************************************************************************/
    // 2) stabilization of continuity equation
    if (my::fldpara_->CStab() == INPAR::FLUID::continuity_stab_yes)
    {
      dserror("continuity stabilization not implemented for poroelasticity");
      /*
      LINALG::Matrix<my::nsd_,my::nsd_> contstab(true);
      const double conti_stab_fac = timefacfacpre*my::tau_(2);
      const double conti_stab_rhs = rhsfac*my::tau_(2)*my::conres_old_;

      * continuity stabilisation on left hand side *
      *
                 /                        \
                |                          |
           tauC | nabla o Du  , nabla o v  |
                |                          |
                 \                        /


      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = my::nsd_*ui;

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          const int fui_p_idim = fui+idim;
          const double v0 = conti_stab_fac* (porosity*my::derxy_(idim,ui)
                            +  grad_porosity(idim) * my::funct_(ui));
          for (int vi=0; vi<my::nen_; ++vi)
          {
            const int fvi = my::nsd_*vi;

            for(int jdim=0;jdim<my::nsd_;++jdim)
            {
              estif_u(fvi+jdim,fui_p_idim) += v0*my::derxy_(jdim, vi) ;
            }
          }
        } // end for(idim)
      }

      //auxiliary variables
      double vel_grad_porosity = 0.0;
      LINALG::Matrix<1,my::nen_> dgradphi_dp_gridvel ;
      LINALG::Matrix<1,my::nen_>  dgradphi_dp_velint;
      dgradphi_dp_gridvel.MultiplyTN(gridvelint,dgradphi_dp);
      dgradphi_dp_velint.MultiplyTN(my::velint_,dgradphi_dp);

      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        vel_grad_porosity += grad_porosity(idim)*my::velint_(idim);
      }

      // pressure terms on left-hand side
      * poroelasticity term *
      *
           /                            \
          |                   n+1        |
          | d(grad(phi))/dp* u    Dp, q  |
          |                   (i)        |
           \                            /

           /                            \
          |                  n+1        |
       +  | d(phi)/dp * div u    Dp, q  |
          |                  (i)        |
           \                            /

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = my::nsd_*ui;
        const double v0 = conti_stab_fac* ( dphi_dp*my::vdiv_*my::funct_(ui)
                                          +  dgradphi_dp_velint(ui)
                                            );
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = my::nsd_*vi;

          for(int jdim=0;jdim<my::nsd_;++jdim)
          {
            estif_p_v(fvi+jdim,fui) += v0*my::derxy_(jdim, vi) ;
          }
        }
      }

      for(int idim=0;idim<my::nsd_;++idim)
      {
        contstab(idim,idim)-=conti_stab_rhs;
      }

      // computation of right-hand-side term
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for (int idim = 0; idim < my::nsd_; ++idim)
        {
          for (int jdim = 0; jdim < my::nsd_; ++jdim)
            velforce(idim,vi)-= contstab(idim,jdim)*my::derxy_(jdim,vi);
        }
      }
      */
    }

/************************************************************************/
    // 3) standard Galerkin pressure term + poroelasticity terms
    /* pressure term */
    /*
         /                \
        |                  |
        |  Dp , nabla o v  |
        |                  |
         \                /
    */
    /* poroelasticity pressure term */
    /*
         /                           \     /                           \
        |         n+1             |     |         n+1             |
        |  sigma*u  * dphi/dp*Dp , v  |  -  |  sigma*vs  * dphi/dp*Dp , v |
        |         (i)                 |     |         (i)                 |
         \                           /       \                           /
    */
    for (int ui=0; ui<my::nen_; ++ui)
     {
       const double v = timefacfacpre*my::funct_(ui);
       const double w = my::funct_(ui)*timefacfacpre*my::reacoeff_/porosity*dphi_dp;
       for (int vi=0; vi<my::nen_; ++vi)
       {
         const int fvi = my::nsd_*vi;
         for (int idim = 0; idim <my::nsd_; ++idim)
         {
           estif_p_v(fvi + idim,ui) += v * ( -my::derxy_(idim, vi) )
                     + w *
                     (
                         my::velint_(idim)
                      - gridvelint(idim)
                     )*my::funct_(vi)
                         ;
         }
       }
     }

     const double pressfac = press*rhsfac;

     for (int vi=0; vi<my::nen_; ++vi)
     {
       /* pressure term on right-hand side */
       for (int idim = 0; idim <my::nsd_; ++idim)
       {
         velforce(idim,vi)+= pressfac *  ( my::derxy_(idim, vi) );
       }
     }  //end for(idim)

/************************************************************************/
    // 4) standard Galerkin continuity term + poroelasticity terms

     if( my::fldpara_->PoroContiPartInt() == false )
     {
       for (int vi=0; vi<my::nen_; ++vi)
       {
         const double v = timefacfacpre*my::funct_(vi);
         for (int ui=0; ui<my::nen_; ++ui)
         {
           const int fui = my::nsd_*ui;

           for (int idim = 0; idim <my::nsd_; ++idim)
           {
             /* continuity term */
             /*
                  /                      \
                 |                        |
                 | phi * nabla o Du  , q  |
                 |                        |
                  \                      /
             */
               /* porosity gradient term */
               /*
                    /                   \
                   |                     |
                   | grad(phi)* Du  , q  |
                   |                     |
                    \                   /
               */
             estif_q_u(vi,fui+idim) += v * ( porosity * my::derxy_(idim,ui)
                         +  grad_porosity(idim) * my::funct_(ui)
                         );
           }
         }
       }  // end for(idim)

       //auxiliary variables
       double vel_grad_porosity = 0.0;
       LINALG::Matrix<1,my::nen_> dgradphi_dp_gridvel ;
       LINALG::Matrix<1,my::nen_>  dgradphi_dp_velint;
       dgradphi_dp_gridvel.MultiplyTN(gridvelint,dgradphi_dp);
       dgradphi_dp_velint.MultiplyTN(my::velint_,dgradphi_dp);

       for (int idim = 0; idim <my::nsd_; ++idim)
       {
         vel_grad_porosity += grad_porosity(idim)*my::velint_(idim);
       }

       // pressure terms on left-hand side
       /* poroelasticity term */
       /*
            /                            \
           |                   n+1        |
           | d(grad(phi))/dp* u    Dp, q  |
           |                   (i)        |
            \                            /

            /                            \
           |                  n+1        |
        +  | d(phi)/dp * div u    Dp, q  |
           |                  (i)        |
            \                            /
       */

       for (int vi=0; vi<my::nen_; ++vi)
       {
         const double v=timefacfacpre*my::funct_(vi);

         for (int ui=0; ui<my::nen_; ++ui)
         {
           ppmat(vi,ui)+= v * ( dphi_dp*my::vdiv_*my::funct_(ui)
                 +  dgradphi_dp_velint(ui)
                   );
         } // ui
       }  // vi

       //right-hand side
       const double rhsfac_vdiv = rhsfac * my::vdiv_;
       for (int vi=0; vi<my::nen_; ++vi)
       {
         // velocity term on right-hand side
         preforce(vi) -= rhsfac_vdiv * porosity * my::funct_(vi)
                     + rhsfac * vel_grad_porosity * my::funct_(vi)
                     ;
       }

       //transient porosity terms
       /*
            /                             \   /                                          \
           |                   n+1         |    |                      /   n+1  \             |
         - | d(grad(phi))/dp* vs    Dp, q  |  + | d(phi)/(dJdp) * J *div| vs       |  * Dp , q  |
           |                   (i)         |    |                      \  (i)   /             |
            \                             /    \                                             /

            /                    \       /                                \
           |                      |   |                    n+1           |
         + | d(phi)/dp *  Dp , q  | + | d^2(phi)/(dp)^2 * p   *  Dp , q  |
           |                      |   |                    (i)           |
            \                    /       \                                /

       */

       if (my::fldpara_->IsStationary() == false)
       {
         for (int vi=0; vi<my::nen_; ++vi)
         {
           const double v = timefacfacpre*my::funct_(vi);
           const double w = my::fac_ * my::funct_(vi);
           for (int ui=0; ui<my::nen_; ++ui)
           {
               ppmat(vi,ui) += - v * dgradphi_dp_gridvel(ui)
                   + v * ( dphi_dJdp * J * gridvdiv )* my::funct_(ui)
                   + w * my::funct_(ui) *  dphi_dp
                   + v * dphi_dpp * my::funct_(ui) * press_dot
                               ;
           }
         }  // end for(idim)

       // inertia terms on the right hand side for instationary fluids

         for (int vi=0; vi<my::nen_; ++vi)
         {
           preforce(vi)-= rhsfac * ( press_dot *
                       dphi_dp
                     ) * my::funct_(vi) ;
         }

         double    grad_porosity_gridvelint=0.0;
         for (int j =0; j< my::nsd_; j++)
         {
           grad_porosity_gridvelint += grad_porosity(j) * gridvelint(j);
         }
        //coupling term on right hand side
         for (int vi=0; vi<my::nen_; ++vi)
         {
           preforce(vi) -= rhsfac *my::funct_(vi) * (- grad_porosity_gridvelint );
           preforce(vi) -= rhsfac  *my::funct_(vi) *  dphi_dJ  * J * gridvdiv;
         }
       }  // end if (not stationary)
     }
     else
     {
       for (int vi=0; vi<my::nen_; ++vi)
       {
         for (int ui=0; ui<my::nen_; ++ui)
         {
           const int fui = my::nsd_*ui;

           for (int idim = 0; idim <my::nsd_; ++idim)
           {
               /* porosity convective term */
               /*
                    /                   \
                   |                     |
                   | phi * Du       , q  |
                   |                     |
                    \                   /
               */
             estif_q_u(vi,fui+idim) += timefacfacpre * my::derxy_(idim,vi) *
                                       (
                                           -1.0 * porosity * my::funct_(ui)
                                       );
           }
         }
       }  // end for(idim)

       LINALG::Matrix<1,my::nen_> deriv_vel ;
       deriv_vel.MultiplyTN(my::velint_,my::derxy_);
       //stationary right-hand side
       for (int vi=0; vi<my::nen_; ++vi)
       {
         // velocity term on right-hand side
         preforce(vi) -= - rhsfac * porosity * deriv_vel(vi)
                     ;
       }


       // pressure terms on left-hand side
       /*
            /                                   \
           |                   n+1               |
           | -d(phi)/dp      * u    Dp, grad(q)  |
           |                   (i)               |
            \                                   /
       */

       for (int vi=0; vi<my::nen_; ++vi)
       {
         for (int ui=0; ui<my::nen_; ++ui)
         {
           ppmat(vi,ui)+= - timefacfacpre * dphi_dp * deriv_vel(vi) * my::funct_(ui)
                   ;
         } // ui
       }  // vi

       if (my::fldpara_->IsStationary() == false)
       {
         //transient porosity terms
         /*
             /                                             \
            |                        /   n+1  \             |
            | d(phi)/(dJdp) * J *div| vs       |  * Dp , q  |
            |                        \  (i)   /             |
             \                                             /

              /                    \       /                                \
             |                      |   |                    n+1           |
           + | d(phi)/dp *  Dp , q  | + | d^2(phi)/(dp)^2 * p   *  Dp , q  |
             |                      |   |                    (i)           |
              \                    /       \                                /

              /                            \
             |                  n+1        |
          +  | d(phi)/dp * div vs   Dp, q  |
             |                  (i)        |
              \                            /

         */

         LINALG::Matrix<1,my::nen_> deriv_gridvel ;
         deriv_gridvel.MultiplyTN(gridvelint,my::derxy_);

         for (int vi=0; vi<my::nen_; ++vi)
         {
           const double v = timefacfacpre*my::funct_(vi);
           const double w = my::fac_ * my::funct_(vi);
           for (int ui=0; ui<my::nen_; ++ui)
           {
               ppmat(vi,ui) +=
                   timefacfacpre * dphi_dp * deriv_gridvel(vi) * my::funct_(ui)
                   + v * ( (dphi_dJdp * J + dphi_dp) * gridvdiv )* my::funct_(ui)
                   + w * my::funct_(ui) *  dphi_dp
                   + v * dphi_dpp * my::funct_(ui) * press_dot
                               ;
           }
         }  // end for(idim)

       // inertia terms on the right hand side for instationary fluids

         for (int vi=0; vi<my::nen_; ++vi)
         {
           preforce(vi)-= rhsfac * ( press_dot *
                       dphi_dp
                     ) * my::funct_(vi) ;
         }

        //coupling term on right hand side
         for (int vi=0; vi<my::nen_; ++vi)
         {
           preforce(vi) -= rhsfac * porosity * deriv_gridvel(vi);
           preforce(vi) -= rhsfac  *my::funct_(vi) *  (dphi_dJ  * J + porosity) * gridvdiv;
         }

       }  // end if (not stationary)
     }
/***********************************************************************************************************/

    // 5) standard Galerkin bodyforce term on right-hand side
    this->BodyForceRhsTerm(velforce,
                     rhsfac);

    // 6) PSPG term
    if (my::fldpara_->PSPG() == INPAR::FLUID::pstab_use_pspg)
    {
      PSPG(estif_q_u,
           ppmat,
           preforce,
           lin_resM_Du,
           0.0,
           timefacfac,
           timefacfacpre,
           rhsfac);
    }

    // 7) reactive stabilization term
    if (my::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
    {
      this->ReacStab(estif_u,
               estif_p_v,
               velforce,
               lin_resM_Du,
               timefacfac,
               timefacfacpre,
               rhsfac,
               0.0);
    }

  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  for (int vi=0; vi<my::nen_; ++vi)
  {
    eforce(my::numdofpernode_*vi+my::nsd_)+=preforce(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int idim=0; idim<my::nsd_; ++idim)
    {
      eforce(my::numdofpernode_*vi+idim)+=velforce(idim,vi);
    }
  }

  // add pressure-pressure part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int fuipp = my::numdofpernode_*ui+my::nsd_;

    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int numdof_vi_p_nsd = my::numdofpernode_*vi+my::nsd_;

      estif(numdof_vi_p_nsd,fuipp)+=ppmat(vi,ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int numdof_ui = my::numdofpernode_*ui;
    const int nsd_ui = my::nsd_*ui;

    for (int jdim=0; jdim < my::nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int numdof_vi = my::numdofpernode_*vi;
        const int nsd_vi = my::nsd_*vi;

        for (int idim=0; idim <my::nsd_; ++idim)
        {
          estif(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int numdof_ui_nsd = my::numdofpernode_*ui + my::nsd_;

    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int nsd_vi = my::nsd_*vi;
      const int numdof_vi = my::numdofpernode_*vi;

      for (int idim=0; idim <my::nsd_; ++idim)
      {
        estif(numdof_vi+idim, numdof_ui_nsd) += estif_p_v(nsd_vi+idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    const int numdof_ui = my::numdofpernode_*ui;
    const int nsd_ui = my::nsd_*ui;

    for (int jdim=0; jdim < my::nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<my::nen_; ++vi)
        estif(my::numdofpernode_*vi+my::nsd_, numdof_ui_jdim) += estif_q_u(vi, nsd_ui_jdim);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculate coupling matrix flow                          vuong 06/11 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::SysmatOD(
    int                                                             eid,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       ebofoaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
    const LINALG::Matrix<my::nen_, 1>&                              epreaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       eaccam,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       emhist,
    const LINALG::Matrix<my::nen_, 1> &                             epren,
    const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
    const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_,my::nsd_ * my::nen_>&  ecoupl,
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
    Teuchos::RCP<const MAT::Material>                               material,
    bool                                                            isale,
    const DRT::UTILS::GaussIntegration &                            intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_> ecoupl_u(true); // coupling matrix for momentum equation
  LINALG::Matrix<my::nen_, my::nen_ * my::nsd_> ecoupl_p(true); // coupling matrix for continuity equation
  LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nen_ * my::nsd_> emesh(true); // linearisation of mesh motion

  // definition of vectors
  LINALG::Matrix<my::nen_, 1> preforce(true);
  LINALG::Matrix<my::nsd_, my::nen_> velforce(true);

  //material coordinates xyze0
  LINALG::Matrix<my::nsd_, my::nen_> xyze0 = my::xyze_;


  // add displacement when fluid nodes move in the ALE case (in poroelasticity this is always the case)
  my::xyze_ += edispnp;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter(eid);

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  //for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)

  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {

    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    //------------------------get determinant of Jacobian dX / ds
    // transposed jacobian "dX/ds"
    LINALG::Matrix<my::nsd_,my::nsd_> xjm0;
    xjm0.MultiplyNT(my::deriv_,xyze0);

    // inverse of transposed jacobian "ds/dX"
    LINALG::Matrix<my::nsd_,my::nsd_> xji0(true);
    const double det0= xji0.Invert(xjm0);

    // ----------------------compute my::derivatives N_XYZ at gp w.r.t. material coordinates
    LINALG::Matrix<my::nsd_,my::nen_> N_XYZ(false);
    N_XYZ.Multiply(xji0,my::deriv_);
    LINALG::Matrix<my::numderiv2_,my::nen_> N_XYZ2;

    if(my::is_higher_order_ele_)
    {
      // get the second derivatives of standard element at current GP w.r.t. rst
      DRT::UTILS::shape_function_deriv2<distype>(my::xsi_,my::deriv2_);

      // get the second my::derivatives of standard element at current GP w.r.t. xyz
      DRT::UTILS::gder2<distype>(my::xjm_,my::derxy_,my::deriv2_,my::xyze_,my::derxy2_);

      // get the second derivatives of standard element at current GP w.r.t. XYZ
      DRT::UTILS::gder2<distype>(xjm0,N_XYZ,my::deriv2_,xyze0,N_XYZ2);
    }
    else
    {
      my::deriv2_.Clear();
      my::derxy2_.Clear();
    }

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
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = my::funct_.Dot(epreaf);

    // get pressure time my::derivative at integration point
    // (value at n+1 )
    double press_dot = my::funct_.Dot(epressnp_timederiv);

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::gradp_.Multiply(my::derxy_,epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::bodyforce_.Multiply(ebofoaf,my::funct_);

    // get momentum history data at integration point
    // (only required for one-step-theta and BDF2 time-integration schemes)
    my::histmom_.Multiply(emhist,my::funct_);

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    LINALG::Matrix<my::nsd_,1> gridvelint;
    gridvelint.Multiply(egridv,my::funct_);

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
    double gridvdiv = 0.0;
    if (not my::fldpara_->IsGenalphaNP())
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        my::vdiv_ += my::vderxy_(idim, idim);

        gridvdiv += gridvelderxy(idim,idim);
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

        gridvdiv += gridvelderxy(idim,idim);
      }
    }

    // -------------------------(material) deformation gradient F = d my::xyze_ / d XYZE = my::xyze_ * N_XYZ^T
    LINALG::Matrix<my::nsd_,my::nsd_> defgrd(false);
    defgrd.MultiplyNT(my::xyze_,N_XYZ);

    // inverse deformation gradient F^-1
    LINALG::Matrix<my::nsd_,my::nsd_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //------------------------------------ build F^-T as vector 9x1
    LINALG::Matrix<my::nsd_*my::nsd_,1> defgrd_IT_vec;
    for(int i=0; i<my::nsd_; i++)
      for(int j=0; j<my::nsd_; j++)
        defgrd_IT_vec(i*my::nsd_+j) = defgrd_inv(j,i);

    // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    double J = my::det_/det0;

    //************************************************auxilary variables for computing the porosity

    double dphi_dp=0.0;
    double dphi_dJ=0.0;
    double dphi_dJdp=0.0;
    double dphi_dJJ=0.0;
    double dphi_dpp=0.0;
    double porosity=0.0;

    //access structure discretization
    RCP<DRT::Discretization> structdis = Teuchos::null;
    structdis = DRT::Problem::Instance()->GetDis("structure");
    //get corresponding structure element (it has the same global ID as the fluid element)
    DRT::Element* structele = structdis->gElement(eid);
    if(structele == NULL)
      dserror("Structure element %i not on local processor", eid);
    //get structure material
    const MAT::StructPoro* structmat = static_cast<const MAT::StructPoro*>((structele->Material()).get());

    structmat->ComputePorosity(press, J, *(iquad),porosity,dphi_dp,dphi_dJ,dphi_dJdp,dphi_dJJ,dphi_dpp);
    //porosity = structmat->GetPorosityAtGP(*(iquad));

    //--------------------------- build d^2 N/(dX dx) at gausspoint (wrt xyt)
    //! second derivatives are orderd as followed: (N,Xx ; N,Yy ; N,Zz ; N,Xy ; N,Xz ; N,Yx ; N,Yz;  N,Zx ; N,Zy)
    LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> N_X_x(true);

    LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_> F_x(true);

    // F^-T : N_X_x
    LINALG::Matrix<my::nsd_,my::nsd_*my::nen_> Finv_N_X_x(true);

    CalcAuxiliaryDerivatives( edispnp,
                              defgrd,
                              defgrd_inv,
                              N_X_x,
                              F_x,
                              Finv_N_X_x);

    //-------------------- compute gradJ = dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx  at gausspoint
    //LINALG::Matrix<1,my::nsd_> gradJ;
    //gradJ2.MultiplyTN(J, defgrd_IT_vec, F_x);
    //gradJ.MultiplyTN(J, defgrd_IT_vec, F_X);
    //gradJ.Multiply(gradJ,defgrd_inv);

    //ask structure material for the spatial gradient of the jacobian
    LINALG::Matrix<1,my::nsd_> gradJ;
    structmat->GetGradJAtGP(gradJ,*(iquad));

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J * N_x
    LINALG::Matrix<1,my::nsd_*my::nen_> dJ_dus;

    for (int i=0; i<my::nen_; i++)
      for (int j=0; j<my::nsd_; j++)
        dJ_dus(j+i*my::nsd_)=J*my::derxy_(j,i);

    //---------------------d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x

    //dF^-T/dus : dF/dx = - (F^-1. dN/dx . u_s)^T  : dF/dx
    LINALG::Matrix<my::nsd_,my::nsd_*my::nen_> dFinvdus_dFdx(true);
    for (int i=0; i<my::nsd_; i++)
      for (int n =0; n<my::nen_; n++)
        for(int j=0; j<my::nsd_; j++)
        {
          const int gid = my::nsd_ * n +j;
          for (int k=0; k<my::nsd_; k++)
            for(int p=0; p<my::nsd_; p++)
              dFinvdus_dFdx(p, gid) += -defgrd_inv(i,j) * my::derxy_(k,n) * F_x(k*my::nsd_+i,p);
        }

    //--------------------- linearization of porosity w.r.t. structure displacements
    LINALG::Matrix<1,my::nsd_*my::nen_> dphi_dus;
    dphi_dus.Update( dphi_dJ , dJ_dus );


    //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) + d\phi/(dp*du_s) * dp/dx
    LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> dgradphi_dus(true);

    // spatial porosity gradient
    LINALG::Matrix<my::nsd_,1>             grad_porosity(true);

    if( (my::fldpara_->PoroContiPartInt() == false) or visceff_)
    {
      //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
      LINALG::Matrix<1,my::nsd_> temp2;
      temp2.MultiplyTN( defgrd_IT_vec, F_x);

      //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
      LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> dgradJ_dus;

      dgradJ_dus.MultiplyTN(temp2,dJ_dus);

      dgradJ_dus.Update(J,dFinvdus_dFdx,1.0);

      dgradJ_dus.Update(J,Finv_N_X_x,1.0);

      //--------------------- current porosity gradient
      //LINALG::Matrix<my::nsd_,1> grad_porosity;
      //for (int idim=0; idim<my::nsd_; ++idim)
      //  grad_porosity(idim)=dphi_dp*my::gradp_(idim)+dphi_dJ*gradJ(idim);

      LINALG::Matrix<my::nsd_,1>             mat_grad_porosity;
      //-----get material porosity gradient from structure material
      structmat->GetGradPorosityAtGP(mat_grad_porosity,*(iquad));

      grad_porosity.MultiplyTN(defgrd_inv,mat_grad_porosity);

      //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) + d\phi/(dp*du_s) * dp/dx
      LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> dgradphi_dus;
      dgradphi_dus.MultiplyTN(dphi_dJJ, gradJ ,dJ_dus);
      dgradphi_dus.Update(dphi_dJ, dgradJ_dus, 1.0);
      dgradphi_dus.Multiply(dphi_dJdp, my::gradp_, dJ_dus, 1.0);
    }

    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point

    if (my::fldpara_->MatGp())
    {
      const MAT::FluidPoro* actmat = static_cast<const MAT::FluidPoro*>(material.get());

      // set density at n+alpha_F/n+1 and n+alpha_M/n+1
      //my::densaf_ = actmat->Density();
      //my::densam_ = my::densaf_;

      // calculate reaction coefficient
      my::reacoeff_ = actmat->ComputeReactionCoeff() * porosity;
    }
    else dserror("Fluid material parameters have to be evaluated at gauss point for porous flow!");

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = my::fldpara_->TimeFac() * my::fac_;
    const double timefacfacpre = my::fldpara_->TimeFacPre() * my::fac_;


    //******************* FAD ************************
/*
     // sacado data type replaces "double"
     typedef Sacado::Fad::DFad<double> FAD;  // for first my::derivs
     // sacado data type replaces "double" (for first+second my::derivs)
     typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FADFAD;

     LINALG::TMatrix<FAD,my::nen_,1> fad_funct_;
     LINALG::TMatrix<FAD,my::nen_,1> fad_epreaf;
     LINALG::TMatrix<FAD,my::nen_,1> fad_epressnp_timederiv;
     LINALG::TMatrix<FAD,3,my::nen_> fad_derxy_;
     LINALG::TMatrix<FAD,3,my::nen_> fad_deriv_;
     //   LINALG::TMatrix<FAD,my::nsd_,my::nen_> fad_derxy2_;
     LINALG::TMatrix<FAD,my::nsd_,1> fad_gradp_;
     LINALG::TMatrix<FAD,3,my::nen_> fad_xyze_;
     LINALG::TMatrix<FAD,3,my::nen_> fad_xcurr_;
     LINALG::TMatrix<FAD,3,my::nen_> fad_xyze0;
     LINALG::TMatrix<FAD,3,my::nen_> fad_edispnp;
     LINALG::TMatrix<FAD,3,my::nen_> fad_N_XYZ;
     LINALG::TMatrix<FAD,9,my::nen_> fad_N_X_x;
     LINALG::TMatrix<FAD,3,3> fad_xji_;
     LINALG::TMatrix<FAD,3,3> fad_xjm_;

     for (int i=0; i<my::nen_; i++)
     {
     fad_funct_(i)=my::funct_(i);
     fad_epreaf(i)=epreaf(i);
     //fad_epressnp_timemy::deriv(i)=epressnp_timemy::deriv(i);
     //fad_epreaf(i).diff(i,my::nen_);
     for(int j=0; j<my::nsd_; j++)
     {
     fad_derxy_(j,i)=my::derxy_(j,i);
     fad_deriv_(j,i)=my::deriv_(j,i);
     fad_xyze0(j,i)= xyze0(j,i);
     fad_edispnp(j,i)= edispnp(j,i);
     fad_edispnp(j,i).diff(i*my::nsd_+j, my::nsd_*my::nen_);
     fad_N_XYZ(j,i) = N_XYZ(j,i);
     }
     for(int j=0; j<9; j++)
     fad_N_X_x(j,i)=N_X_x(j,i);
     }

     fad_xyze_.Update(1.0,fad_xyze0,1.0,fad_edispnp);

     FAD fad_press = fad_funct_.Dot(fad_epreaf);
     LINALG::TMatrix<FAD,3,1> fad_gradp;
     fad_gradp.Multiply(fad_derxy_,fad_epreaf);


     // compute F
     LINALG::TMatrix<FAD,3,3> fad_defgrd(false);
     fad_defgrd.MultiplyNT(fad_xyze_,fad_N_XYZ);
     FAD fad_J = Determinant3x3(fad_defgrd);


     LINALG::TMatrix<FAD,3,3>    fad_defgrd_inv(false);
     fad_defgrd_inv = fad_defgrd;
     Inverse3x3(fad_defgrd_inv);

     LINALG::TMatrix<FAD,3,3> fad_cauchygreen;
     fad_cauchygreen.MultiplyTN(fad_defgrd,fad_defgrd);

     LINALG::TMatrix<FAD,3,3>    fad_C_inv(false);
     fad_C_inv = fad_cauchygreen;
     Inverse3x3(fad_C_inv);

     LINALG::TMatrix<FAD,6,1> fad_C_inv_vec(false);
     fad_C_inv_vec(0) = fad_C_inv(0,0);
     fad_C_inv_vec(1) = fad_C_inv(1,1);
     fad_C_inv_vec(2) = fad_C_inv(2,2);
     fad_C_inv_vec(3) = fad_C_inv(0,1);
     fad_C_inv_vec(4) = fad_C_inv(1,2);
     fad_C_inv_vec(5) = fad_C_inv(2,0);


     LINALG::TMatrix<FAD,9,1> fad_defgrd_IT_vec;
     fad_defgrd_IT_vec(0)=fad_defgrd_inv(0,0);
     fad_defgrd_IT_vec(1)=fad_defgrd_inv(1,0);
     fad_defgrd_IT_vec(2)=fad_defgrd_inv(2,0);
     fad_defgrd_IT_vec(3)=fad_defgrd_inv(0,1);
     fad_defgrd_IT_vec(4)=fad_defgrd_inv(1,1);
     fad_defgrd_IT_vec(5)=fad_defgrd_inv(2,1);
     fad_defgrd_IT_vec(6)=fad_defgrd_inv(0,2);
     fad_defgrd_IT_vec(7)=fad_defgrd_inv(1,2);
     fad_defgrd_IT_vec(8)=fad_defgrd_inv(2,2);

     LINALG::TMatrix<FAD,9,my::nsd_> fad_F_x(true);
     for(int i=0; i<my::nsd_; i++)
     {
     for(int n=0; n<my::nen_; n++)
     {
     fad_F_x(i*my::nsd_+0, 0) +=   fad_N_X_x(0,n)*fad_edispnp(i,n);
     fad_F_x(i*my::nsd_+1, 0) +=   fad_N_X_x(5,n)*fad_edispnp(i,n);
     fad_F_x(i*my::nsd_+2, 0) +=   fad_N_X_x(7,n)*fad_edispnp(i,n);

     fad_F_x(i*my::nsd_+0, 1) +=   fad_N_X_x(3,n)*fad_edispnp(i,n) ;
     fad_F_x(i*my::nsd_+1, 1) +=   fad_N_X_x(1,n)*fad_edispnp(i,n) ;
     fad_F_x(i*my::nsd_+2, 1) +=   fad_N_X_x(8,n)*fad_edispnp(i,n) ;

     fad_F_x(i*my::nsd_+0, 2) +=   fad_N_X_x(4,n)*fad_edispnp(i,n);
     fad_F_x(i*my::nsd_+1, 2) +=   fad_N_X_x(6,n)*fad_edispnp(i,n);
     fad_F_x(i*my::nsd_+2, 2) +=   fad_N_X_x(2,n)*fad_edispnp(i,n);
     }
     }

     LINALG::TMatrix<FAD,1,my::nsd_> fad_gradJ;
     fad_gradJ.MultiplyTN(fad_J, fad_defgrd_IT_vec, fad_F_x);

     double initporosity_  = structmat->Initporosity();
     double bulkmodulus_ = structmat->Bulkmodulus();
     double penalty_ = structmat->Penaltyparameter();

     FAD fad_a     = ( bulkmodulus_/(1-initporosity_) + fad_press - penalty_/initporosity_ ) * fad_J;
     FAD fad_b     = -fad_a + bulkmodulus_ + penalty_;
     FAD fad_c    = (fad_b/fad_a) * (fad_b/fad_a) + 4*penalty_/fad_a;
     FAD fad_d     = sqrt(fad_c)*fad_a;

     FAD fad_sign = 1.0;

     FAD fad_test = 1 / (2 * fad_a) * (-fad_b + fad_d);
     if (fad_test >= 1.0 or fad_test < 0.0)
     {
       fad_sign = -1.0;
       fad_d = fad_sign * fad_d;
     }

     FAD fad_porosity = 1/(2*fad_a)*(-fad_b+fad_d);

     FAD fad_d_p   =  fad_J * (-fad_b+2*penalty_)/fad_d;
     FAD fad_d_J   =  fad_a/fad_J * ( -fad_b + 2*penalty_ ) / fad_d;

     FAD fad_dphi_dp=  - fad_J * fad_porosity/fad_a + (fad_J+fad_d_p)/(2*fad_a);
     FAD fad_dphi_dJ=  -fad_porosity/fad_J+ 1/(2*fad_J) + fad_d_J / (2*fad_a);

     LINALG::TMatrix<FAD,1,my::nsd_>             fad_grad_porosity;
     for (int idim=0; idim<my::nsd_; ++idim)
     {
     fad_grad_porosity(idim)=fad_dphi_dp*fad_gradp(idim)+fad_dphi_dJ*fad_gradJ(idim);
     }

     for (int i=0; i<my::nsd_*my::nen_; i++)
     {
     if( (dJ_dus(i)-fad_J.dx(i)) > 1e-8)
     {
     cout<<"dJdus("<<i<<"): "<<dJ_dus(i)<<endl;
     cout<<"fad_J.dx("<<i<<"): "<<fad_J.dx(i)<<endl;
     dserror("check dJdus failed!");
     }
     }
     cout<<"dJdus check done and ok"<<endl;


     for (int i=0; i<my::nsd_*my::nen_; i++)
     for (int j=0; j<my::nsd_; j++)
     {
     if( (dgradJ_dus(j,i)-fad_gradJ(j).dx(i)) > 1e-8)
     {
     cout<<"dgradJ_dus("<<i<<"): "<<dgradJ_dus(j,i)<<endl;
     cout<<"fad_gradJ.dx("<<i<<"): "<<fad_gradJ(j).dx(i)<<endl;
     cout<<"gradJ:"<<endl<<gradJ<<endl;
     cout<<"fad_gradJ:"<<endl<<fad_gradJ<<endl;
     dserror("check dgradJ_dus failed!");
     }
     }
     cout<<"dgradJ_dus check done and ok"<<endl;

     for (int i=0; i<my::nsd_*my::nen_; i++)
     if( (dphi_dus(i)-fad_porosity.dx(i)) > 1e-8)
     {
     cout<<"dphi_dus("<<i<<"): "<<dphi_dus(i)<<endl;
     cout<<"fad_porosity.dx("<<i<<"): "<<fad_porosity.dx(i)<<endl;
     cout<<"dphi_dus:"<<endl<<dphi_dus<<endl;
     cout<<"fad_porosity:"<<endl<<fad_porosity<<endl;
     dserror("check dgradJ_dus failed!");
     }
     cout<<"dphi_dus check done and ok"<<endl;

     for (int i=0; i<my::nsd_*my::nen_; i++)
     for (int j=0; j<my::nsd_; j++)
     {
     if( (dgradphi_dus(j,i)-fad_grad_porosity(j).dx(i)) > 1e-8)
     {
     cout<<"dgradphi_dus("<<i<<"): "<<dgradphi_dus(j,i)<<endl;
     cout<<"fad_grad_porosity.dx("<<i<<"): "<<fad_grad_porosity(j).dx(i)<<endl;
     cout<<"dgradphi_dus:"<<endl<<dgradphi_dus<<endl;
     cout<<"fad_grad_porosity:"<<endl<<fad_grad_porosity<<endl;
     dserror("check dgradphi_dus failed!");
     }
     }
     cout<<"dgradphi_dus check done and ok"<<endl;
*/
    //******************* FAD ************************


    //***********************************************************************************************
    // 1) coupling terms in momentum balance

    //stationary
    /*  reaction */
    /*
     /                                      \
     |                    n+1                |
     |    sigma * dphi/dus * u    * Dus , v  |
     |                        (i)            |
     \                                     /
     */
    /*  reactive ALE term */
    /*
     /                                  \
     |                  n+1             |
     |    - rho * grad u     * Dus , v  |
     |                  (i)             |
     \                                 /
     */

    const double fac_reac= timefacfac*my::reacoeff_/porosity;
    const double fac_densaf=my::fac_*my::densaf_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = my::nsd_*ui;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = my::nsd_*vi;

        for (int idim = 0; idim <my::nsd_; ++idim)
        {
          for (int jdim = 0; jdim <my::nsd_; ++jdim)
          {
            ecoupl_u(fvi+idim,fui+jdim) += my::funct_(vi)*fac_reac* my::velint_(idim)*dphi_dus(fui+jdim)
            - my::funct_(vi) * fac_densaf * my::vderxy_(idim, jdim) * my::funct_(ui)
            ;
          }
        } // end for (idim)
      } //vi
    } // ui

    //transient terms
    /*  reaction */
    /*
     /                            \        /                                           \
     |                             |      |                            n+1              |
  -  |    sigma * phi * D(v_s) , v |  -   |    sigma * d(phi)/d(us) * vs *  D(u_s) , v  |
     |                             |      |                             (i)             |
     \                           /         \                                           /
     */

    if (not my::fldpara_->IsStationary())
    {
      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = my::nsd_*ui;

        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = my::nsd_*vi;
          const double tmp = - my::funct_(vi)* my::reacoeff_/porosity;

          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            ecoupl_u(fvi+idim,fui+idim) += my::fac_ * tmp * porosity * my::funct_(ui);
            for (int jdim =0; jdim<my::nsd_; ++jdim)
              ecoupl_u(fvi+idim,fui+jdim) += timefacfac * tmp * gridvelint(idim) * dphi_dus(fui+jdim)
              ;
          } // end for (idim)
        } //vi
      } // ui
    }

    //viscous terms (brinkman terms)
    if(visceff_)
    {
      LINALG::Matrix<my::nsd_,my::nsd_> viscstress(true);
      //const double v = visceff_;

      for (int jdim = 0; jdim < my::nsd_; ++jdim)
      {
        for (int idim = 0; idim < my::nsd_; ++idim)
        {
          viscstress(idim,jdim)=visceff_*(my::vderxy_(jdim,idim)+my::vderxy_(idim,jdim));
        }
      }

      LINALG::Matrix<my::nsd_,1> viscstress_gradphi(true);
      viscstress_gradphi.Multiply(viscstress,grad_porosity);

      LINALG::Matrix<my::nsd_,my::nen_*my::nsd_> viscstress_dgradphidus(true);
      viscstress_dgradphidus.Multiply(viscstress,dgradphi_dus);

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = my::nsd_*ui;
        const double v = timefacfacpre*my::funct_(ui);
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = my::nsd_*vi;
          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            for (int jdim =0; jdim<my::nsd_; ++jdim)
              ecoupl_u(fvi + idim,fui+jdim) += v * 1/(porosity) * ( 1/porosity*viscstress_gradphi(idim)*
                                                                  dphi_dus(fui+jdim)
                                                              - viscstress_dgradphidus(idim,fui+jdim)
                                                             )
                                          ;
          }
        }
      }
    }

    //*************************************************************************************************************
    // 2) coupling terms in continuity equation

    if( my::fldpara_->PoroContiPartInt() == false )
    {
      //auxiliary variables
      LINALG::Matrix<1,my::nen_*my::nsd_> grad_porosity_us_velint;
      grad_porosity_us_velint.MultiplyTN(my::velint_,dgradphi_dus);

      // structure coupling terms on left-hand side
      /*  stationary */
      /*
        /                                 \      /                                    \
       |                 n+1              |     |                        n+1           |
       |   dphi/dus * div u    * Dus , v  |  +  |   d(grad(phi))/dus * u    * Dus , v  |
       |                  (i)             |     |                       (i)            |
        \                                /       \                                    /
       */
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double v=timefacfacpre*my::funct_(vi);

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = my::nsd_*ui;
          for(int idim = 0; idim <my::nsd_; ++idim)
          {
            ecoupl_p(vi,fui+idim)+= v * dphi_dus(fui+idim) * my::vdiv_
            + v * grad_porosity_us_velint(fui+idim)
            ;
          }
        } // ui
      } // vi

      //transient coupling terms
      if (my::fldpara_->IsStationary() == false)
      {
        LINALG::Matrix<1,my::nen_*my::nsd_> grad_porosity_us_gridvelint;
        grad_porosity_us_gridvelint.MultiplyTN(gridvelint,dgradphi_dus);

        /*
          /                            \       /                                                      \
         |                              |     |                                    n+1                 |
         |   dphi/dJ * J * div Dus , v  |   + |   d^2(phi)/(dJ)^2 * dJ/dus  * J * div vs    * Dus , v  |
         |                              |     |                                        (i)             |
          \                            /      \                                                       /

         /                                            \        /                                        \
         |                           n+1               |      |                           n+1            |
      +  |   dphi/dJ * dJ/dus * div vs    * Dus, v     |    - |   d(grad(phi))/d(us) *  vs    * Dus , v  |
         |                           (i)               |      |                           (i)            |
         \                                            /        \                                        /

            /                       \
           |                         |
         - |    grad phi * Dus , v   |
           |                         |
            \                       /

            /                                          \
           |                             n+1            |
         + |    dphi/(dpdJ) * dJ/dus  * p    * Dus , v  |
           |                             (i)            |
           \                                           /
         */

        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = my::fac_*my::funct_(vi);
          const double w = timefacfacpre*my::funct_(vi);
          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = my::nsd_*ui;
            for(int idim = 0; idim <my::nsd_; ++idim)
            {
              ecoupl_p(vi,fui+idim)+= v * ( dphi_dJ * J * my::derxy_(idim,ui) )
              + w * ( + gridvdiv * ( dphi_dJJ * J + dphi_dJ ) * dJ_dus(fui+idim)
                      - grad_porosity_us_gridvelint(fui+idim)
                    )
              - v * grad_porosity(idim) * my::funct_(ui)
              + w * dphi_dJdp * press_dot * dJ_dus(fui+idim)
              ;
            }
          }
        } // end for(idim)

      } // end if (not stationary)
    }
    else
    {
      LINALG::Matrix<1,my::nen_> deriv_vel ;
      deriv_vel.MultiplyTN(my::velint_,my::derxy_);

      // structure coupling terms on left-hand side
      /*  stationary */
      /*
        /                                    \
       |                 n+1                  |
       |   -dphi/dus *  u    * Dus , grad(v)  |
       |                  (i)                 |
        \                                    /
       */
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = my::nsd_*ui;
          for(int idim = 0; idim <my::nsd_; ++idim)
          {
            ecoupl_p(vi,fui+idim)+= timefacfacpre * (-1.0) * dphi_dus(fui+idim) * deriv_vel(vi)
            ;
          }
        } // ui
      } // vi

      //transient coupling terms
      if (my::fldpara_->IsStationary() == false)
      {

        /*
          /                                    \       /                                                      \
         |                                      |     |                                    n+1                 |
         |   (dphi/dJ * J + phi )* div Dus , v  |   + |   d^2(phi)/(dJ)^2 * dJ/dus  * J * div vs    * Dus , v  |
         |                                      |     |                                        (i)             |
          \                                    /      \                                                       /

         /                                            \
         |                           n+1               |
      +  |   dphi/dJ * dJ/dus * div vs    * Dus, v     |
         |                           (i)               |
         \                                            /

         /                                   \    /                        \
         |                    n+1            |    |                        |
      +  |   dphi/dus * div vs  * Dus, v     |  + |   phi *  Dus, grad(v)  |
         |                    (i)            |    |                        |
         \                                  /     \                       /

            /                                          \
           |                             n+1            |
         + |    dphi/(dpdJ) * dJ/dus  * p    * Dus , v  |
           |                             (i)            |
           \                                           /
         */

        LINALG::Matrix<1,my::nen_> deriv_gridvel ;
        deriv_gridvel.MultiplyTN(gridvelint,my::derxy_);

        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = my::fac_*my::funct_(vi);
          const double w = timefacfacpre*my::funct_(vi);
          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = my::nsd_*ui;
            for(int idim = 0; idim <my::nsd_; ++idim)
            {
              ecoupl_p(vi,fui+idim)+=   v * ( (dphi_dJ * J + porosity)* my::derxy_(idim,ui) )
                                      + my::fac_ * my::derxy_(idim,vi) * ( porosity * my::funct_(ui) )
                                      + w * ( + gridvdiv * (
                                                              ( dphi_dJJ * J + dphi_dJ ) * dJ_dus(fui+idim)
                                                              + dphi_dus(fui+idim)
                                                            )
                                            )
                                      + timefacfacpre * deriv_gridvel(vi) * dphi_dus(fui+idim)
                                      + w * dphi_dJdp * press_dot * dJ_dus(fui+idim)
              ;
            }
          }
        } // end for(idim)

      } // end if (not stationary)
    }// end if (partial integration)

    //*************************************************************************************************************
    // 3) shape derivatives

    if (my::nsd_ == 3)
      LinMeshMotion_3D_OD(
          emesh,
          evelaf,
          egridv,
          press,
          press_dot,
          porosity,
          dphi_dp,
          dphi_dJ,
          J,
          gradJ,
          my::fldpara_->TimeFac(),
          timefacfac);
    else if(my::nsd_ == 2)
      LinMeshMotion_2D_OD(
          emesh,
          evelaf,
          egridv,
          press,
          press_dot,
          porosity,
          dphi_dp,
          dphi_dJ,
          J,
          gradJ,
          my::fldpara_->TimeFac(),
          timefacfac);
    else
      dserror("Linearization of the mesh motion is only available in 3D");

  }//loop over gausspoints

  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  /*

   */

  // add fluid velocity-structure displacement part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    //   const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = my::nsd_*ui;

    for (int jdim=0; jdim < my::nsd_;++jdim)
    {
      //   const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int numdof_vi = my::numdofpernode_*vi;
        const int nsd_vi = my::nsd_*vi;

        for (int idim=0; idim <my::nsd_; ++idim)
        {
          ecoupl(numdof_vi+idim, nsd_ui_jdim) += ecoupl_u(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add fluid pressure-structure displacement part to matrix
  for (int ui=0; ui<my::nen_; ++ui)
  {
    //    const int numdof_ui = my::numdofpernode_*ui;
    const int nsd_ui = my::nsd_*ui;

    for (int jdim=0; jdim < my::nsd_;++jdim)
    {
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<my::nen_; ++vi)
      {
        ecoupl(my::numdofpernode_*vi+my::nsd_, nsd_ui_jdim) += ecoupl_p(vi, nsd_ui_jdim);
      }
    }
  }

  //add linearisation of mesh motion
  //ecoupl.Update(1.0,emesh,1.0);

  return;
}    //SysmatOD

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::LinMeshMotion_3D_OD(
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_) * my::nen_>& emesh,
    const LINALG::Matrix<my::nsd_, my::nen_>&                         evelaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                         egridv,
    const double &                                                    press,
    const double &                                                    press_dot,
    const double &                                                    porosity,
    const double &                                                    dphi_dp,
    const double &                                                    dphi_dJ,
    const double &                                                    J,
    LINALG::Matrix<1, my::nsd_>&                                      gradJ,
    const double & timefac, const double &                            timefacfac)
{

  //*************************** linearisation of mesh motion in momentum balance**********************************
  // mass + rhs

  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = my::fac_ * my::funct_(vi, 0);
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 4    , ui * 3    ) += v * (my::velint_(0) - my::rhsmom_(0) * my::fldpara_->Dt()
                                                     * my::fldpara_->Theta()) * my::derxy_(0, ui);
      emesh(vi * 4    , ui * 3 + 1) += v * (my::velint_(0) - my::rhsmom_(0)
                                                     * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(1, ui);
      emesh(vi * 4    , ui * 3 + 2) += v * (my::velint_(0) - my::rhsmom_(0)
                                                     * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(2, ui);

      emesh(vi * 4 + 1, ui * 3    ) += v * (my::velint_(1) - my::rhsmom_(1)
                                           * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(0, ui);
      emesh(vi * 4 + 1, ui * 3 + 1) += v * (my::velint_(1) - my::rhsmom_(1)
                                           * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(1, ui);
      emesh(vi * 4 + 1, ui * 3 + 2) += v * (my::velint_(1) - my::rhsmom_(1)
                                           * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(2, ui);

      emesh(vi * 4 + 2, ui * 3    ) += v * (my::velint_(2) - my::rhsmom_(2)
                                            * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(0, ui);
      emesh(vi * 4 + 2, ui * 3 + 1) += v * (my::velint_(2) - my::rhsmom_(2)
                                            * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(1, ui);
      emesh(vi * 4 + 2, ui * 3 + 2) += v * (my::velint_(2) - my::rhsmom_(2)
                                            * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(2, ui);
    }
  }

  LINALG::Matrix<my::nsd_, 1> gridvelint;
  gridvelint.Multiply(egridv, my::funct_);

  //---------reaction term (darcy term)
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = timefacfac * my::funct_(vi, 0) * my::reacoeff_;
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 4    , ui * 3    ) += v * (my::velint_(0) - gridvelint(0)) * my::derxy_(0, ui);
      emesh(vi * 4    , ui * 3 + 1) += v * (my::velint_(0) - gridvelint(0)) * my::derxy_(1, ui);
      emesh(vi * 4    , ui * 3 + 2) += v * (my::velint_(0) - gridvelint(0)) * my::derxy_(2, ui);

      emesh(vi * 4 + 1, ui * 3    ) += v * (my::velint_(1) - gridvelint(1)) * my::derxy_(0,ui);
      emesh(vi * 4 + 1, ui * 3 + 1) += v * (my::velint_(1) - gridvelint(1)) * my::derxy_(1, ui);
      emesh(vi * 4 + 1, ui * 3 + 2) += v * (my::velint_(1) - gridvelint(1)) * my::derxy_(2, ui);

      emesh(vi * 4 + 2, ui * 3    ) += v * (my::velint_(2) - gridvelint(2)) * my::derxy_(0,ui);
      emesh(vi * 4 + 2, ui * 3 + 1) += v * (my::velint_(2) - gridvelint(2)) * my::derxy_(1, ui);
      emesh(vi * 4 + 2, ui * 3 + 2) += v * (my::velint_(2) - gridvelint(2)) * my::derxy_(2, ui);
    }
  }

  //---------------convective term

  //vmy::deriv_  = sum(evelaf(i,k) * my::deriv_(j,k), k);
  my::vderiv_.MultiplyNT(evelaf, my::deriv_);
  my::convvelint_.Multiply(-1.0, egridv, my::funct_, 0.0);

#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

#define derxjm_001(ui) (my::deriv_(2, ui)*my::xjm_(1, 2) - my::deriv_(1, ui)*my::xjm_(2, 2))
#define derxjm_002(ui) (my::deriv_(1, ui)*my::xjm_(2, 1) - my::deriv_(2, ui)*my::xjm_(1, 1))

#define derxjm_100(ui) (my::deriv_(1, ui)*my::xjm_(2, 2) - my::deriv_(2, ui)*my::xjm_(1, 2))
#define derxjm_102(ui) (my::deriv_(2, ui)*my::xjm_(1, 0) - my::deriv_(1, ui)*my::xjm_(2, 0))

#define derxjm_200(ui) (my::deriv_(2, ui)*my::xjm_(1, 1) - my::deriv_(1, ui)*my::xjm_(2, 1))
#define derxjm_201(ui) (my::deriv_(1, ui)*my::xjm_(2, 0) - my::deriv_(2, ui)*my::xjm_(1, 0))

#define derxjm_011(ui) (my::deriv_(0, ui)*my::xjm_(2, 2) - my::deriv_(2, ui)*my::xjm_(0, 2))
#define derxjm_012(ui) (my::deriv_(2, ui)*my::xjm_(0, 1) - my::deriv_(0, ui)*my::xjm_(2, 1))

#define derxjm_110(ui) (my::deriv_(2, ui)*my::xjm_(0, 2) - my::deriv_(0, ui)*my::xjm_(2, 2))
#define derxjm_112(ui) (my::deriv_(0, ui)*my::xjm_(2, 0) - my::deriv_(2, ui)*my::xjm_(0, 0))

#define derxjm_210(ui) (my::deriv_(0, ui)*my::xjm_(2, 1) - my::deriv_(2, ui)*my::xjm_(0, 1))
#define derxjm_211(ui) (my::deriv_(2, ui)*my::xjm_(0, 0) - my::deriv_(0, ui)*my::xjm_(2, 0))

#define derxjm_021(ui) (my::deriv_(1, ui)*my::xjm_(0, 2) - my::deriv_(0, ui)*my::xjm_(1, 2))
#define derxjm_022(ui) (my::deriv_(0, ui)*my::xjm_(1, 1) - my::deriv_(1, ui)*my::xjm_(0, 1))

#define derxjm_120(ui) (my::deriv_(0, ui)*my::xjm_(1, 2) - my::deriv_(1, ui)*my::xjm_(0, 2))
#define derxjm_122(ui) (my::deriv_(1, ui)*my::xjm_(0, 0) - my::deriv_(0, ui)*my::xjm_(1, 0))

#define derxjm_220(ui) (my::deriv_(1, ui)*my::xjm_(0, 1) - my::deriv_(0, ui)*my::xjm_(1, 1))
#define derxjm_221(ui) (my::deriv_(0, ui)*my::xjm_(1, 0) - my::deriv_(1, ui)*my::xjm_(0, 0))

  for (int ui = 0; ui < my::nen_; ++ui)
  {
    double v00 = +my::convvelint_(1)* (my::vderiv_(0, 0) * derxjm_(0,0,1,ui)
                      + my::vderiv_(0, 1)* derxjm_(0,1,1,ui) + my::vderiv_(0, 2) * derxjm_(0,2,1,ui))
                 + my::convvelint_(2) * (my::vderiv_(0, 0) * derxjm_(0,0,2,ui)
                      + my::vderiv_(0, 1)* derxjm_(0,1,2,ui) + my::vderiv_(0, 2) * derxjm_(0,2,2,ui));
    double v01 = +my::convvelint_(0) * (my::vderiv_(0, 0) * derxjm_(1,0,0,ui)
                      + my::vderiv_(0, 1) * derxjm_(1,1,0,ui) + my::vderiv_(0, 2) * derxjm_(1,2,0,ui))
                 + my::convvelint_(2) * (my::vderiv_(0, 0) * derxjm_(1,0,2,ui)
                      + my::vderiv_(0, 1) * derxjm_(1,1,2,ui) + my::vderiv_(0, 2) * derxjm_(1,2,2,ui));
    double v02 = +my::convvelint_(0) * (my::vderiv_(0, 0) * derxjm_(2,0,0,ui)
                      + my::vderiv_(0, 1) * derxjm_(2,1,0,ui) + my::vderiv_(0, 2) * derxjm_(2,2,0,ui))
                 + my::convvelint_(1) * (my::vderiv_(0, 0) * derxjm_(2,0,1,ui)
                      + my::vderiv_(0, 1) * derxjm_(2,1,1,ui) + my::vderiv_(0, 2) * derxjm_(2,2,1,ui));
    double v10 = +my::convvelint_(1) * (my::vderiv_(1, 0) * derxjm_(0,0,1,ui)
                      + my::vderiv_(1, 1) * derxjm_(0,1,1,ui) + my::vderiv_(1, 2) * derxjm_(0,2,1,ui))
                 + my::convvelint_(2) * (my::vderiv_(1, 0) * derxjm_(0,0,2,ui)
                      + my::vderiv_(1, 1) * derxjm_(0,1,2,ui) + my::vderiv_(1, 2) * derxjm_(0,2,2,ui));
    double v11 = +my::convvelint_(0) * (my::vderiv_(1, 0) * derxjm_(1,0,0,ui)
                      + my::vderiv_(1, 1) * derxjm_(1,1,0,ui) + my::vderiv_(1, 2) * derxjm_(1,2,0,ui))
                 + my::convvelint_(2) * (my::vderiv_(1, 0) * derxjm_(1,0,2,ui)
                      + my::vderiv_(1, 1) * derxjm_(1,1,2,ui) + my::vderiv_(1, 2) * derxjm_(1,2,2,ui));
    double v12 = +my::convvelint_(0) * (my::vderiv_(1, 0) * derxjm_(2,0,0,ui)
                      + my::vderiv_(1, 1) * derxjm_(2,1,0,ui) + my::vderiv_(1, 2) * derxjm_(2,2,0,ui))
                 + my::convvelint_(1) * (my::vderiv_(1, 0) * derxjm_(2,0,1,ui)
                      + my::vderiv_(1, 1) * derxjm_(2,1,1,ui) + my::vderiv_(1, 2) * derxjm_(2,2,1,ui));
    double v20 = +my::convvelint_(1) * (my::vderiv_(2, 0) * derxjm_(0,0,1,ui)
                      + my::vderiv_(2, 1) * derxjm_(0,1,1,ui) + my::vderiv_(2, 2) * derxjm_(0,2,1,ui))
                 + my::convvelint_(2) * (my::vderiv_(2, 0) * derxjm_(0,0,2,ui)
                      + my::vderiv_(2, 1) * derxjm_(0,1,2,ui) + my::vderiv_(2, 2) * derxjm_(0,2,2,ui));
    double v21 = +my::convvelint_(0) * (my::vderiv_(2, 0) * derxjm_(1,0,0,ui)
                      + my::vderiv_(2, 1) * derxjm_(1,1,0,ui) + my::vderiv_(2, 2) * derxjm_(1,2,0,ui))
                 + my::convvelint_(2) * (my::vderiv_(2, 0) * derxjm_(1,0,2,ui)
                      + my::vderiv_(2, 1) * derxjm_(1,1,2,ui) + my::vderiv_(2, 2) * derxjm_(1,2,2,ui));
    double v22 = +my::convvelint_(0) * (my::vderiv_(2, 0) * derxjm_(2,0,0,ui)
                      + my::vderiv_(2, 1) * derxjm_(2,1,0,ui) + my::vderiv_(2, 2) * derxjm_(2,2,0,ui))
                 + my::convvelint_(1) * (my::vderiv_(2, 0) * derxjm_(2,0,1,ui)
                      + my::vderiv_(2, 1) * derxjm_(2,1,1,ui) + my::vderiv_(2, 2) * derxjm_(2,2,1,ui));

    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = my::densaf_*timefacfac / my::det_ * my::funct_(vi);

      emesh(vi * 4 + 0, ui * 3 + 0) += v * v00;
      emesh(vi * 4 + 0, ui * 3 + 1) += v * v01;
      emesh(vi * 4 + 0, ui * 3 + 2) += v * v02;

      emesh(vi * 4 + 1, ui * 3 + 0) += v * v10;
      emesh(vi * 4 + 1, ui * 3 + 1) += v * v11;
      emesh(vi * 4 + 1, ui * 3 + 2) += v * v12;

      emesh(vi * 4 + 2, ui * 3 + 0) += v * v20;
      emesh(vi * 4 + 2, ui * 3 + 1) += v * v21;
      emesh(vi * 4 + 2, ui * 3 + 2) += v * v22;
    }
  }

  // pressure;
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = press * timefacfac / my::det_;
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 4, ui * 3 + 1) += v * (my::deriv_(0, vi) * derxjm_(0,0,1,ui)
          + my::deriv_(1, vi) * derxjm_(0,1,1,ui) + my::deriv_(2, vi)
          *derxjm_(0,2,1,ui));
      emesh(vi * 4, ui * 3 + 2) += v * (my::deriv_(0, vi) * derxjm_(0,0,2,ui)
          + my::deriv_(1, vi) * derxjm_(0,1,2,ui) + my::deriv_(2, vi)
          *derxjm_(0,2,2,ui));

      emesh(vi * 4 + 1, ui * 3 + 0) += v * (my::deriv_(0, vi) * derxjm_(1,0,0,ui)
          + my::deriv_(1, vi) * derxjm_(1,1,0,ui) + my::deriv_(2, vi)
          * derxjm_(1,2,0,ui));
      emesh(vi * 4 + 1, ui * 3 + 2) += v * (my::deriv_(0, vi) * derxjm_(1,0,2,ui)
          + my::deriv_(1, vi) * derxjm_(1,1,2,ui) + my::deriv_(2, vi)
          * derxjm_(1,2,2,ui));

      emesh(vi * 4 + 2, ui * 3 + 0) += v * (my::deriv_(0, vi) * derxjm_(2,0,0,ui)
          + my::deriv_(1, vi) * derxjm_(2,1,0,ui) + my::deriv_(2, vi)
          * derxjm_(2,2,0,ui));
      emesh(vi * 4 + 2, ui * 3 + 1) += v * (my::deriv_(0, vi) * derxjm_(2,0,1,ui)
          + my::deriv_(1, vi) * derxjm_(2,1,1,ui) + my::deriv_(2, vi)
          * derxjm_(2,2,1,ui));
    }
  }

  // //---------viscous term (brinkman term)
#define xji_00 my::xji_(0,0)
#define xji_01 my::xji_(0,1)
#define xji_02 my::xji_(0,2)
#define xji_10 my::xji_(1,0)
#define xji_11 my::xji_(1,1)
#define xji_12 my::xji_(1,2)
#define xji_20 my::xji_(2,0)
#define xji_21 my::xji_(2,1)
#define xji_22 my::xji_(2,2)

#define xjm(i,j) my::xjm_(i,j)

  if(visceff_)
  {
    //-----------grad(phi) = dphi/dp gradp + dphi/dJ gradJ
    LINALG::Matrix<my::nsd_, 1> gradphi;
    for(int ui=0 ; ui<my::nsd_;++ui)
    {
      gradphi(ui) = dphi_dp*my::gradp_(ui)+dphi_dJ*gradJ(ui);
    }

    // part 1: derivative of 1/det

    double v = visceff_*timefac*my::fac_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double derinvJ0 = -v*(my::deriv_(0,ui)*xji_00 + my::deriv_(1,ui)*xji_01 + my::deriv_(2,ui)*xji_02);
      double derinvJ1 = -v*(my::deriv_(0,ui)*xji_10 + my::deriv_(1,ui)*xji_11 + my::deriv_(2,ui)*xji_12);
      double derinvJ2 = -v*(my::deriv_(0,ui)*xji_20 + my::deriv_(1,ui)*xji_21 + my::deriv_(2,ui)*xji_22);
      for (int vi=0; vi<my::nen_; ++vi)
      {
        double visres0 =   2.0*my::derxy_(0, vi)* my::vderxy_(0, 0)
                           +     my::derxy_(1, vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                           +     my::derxy_(2, vi)*(my::vderxy_(0, 2) + my::vderxy_(2, 0)) ;
        double visres1 =         my::derxy_(0, vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                           + 2.0*my::derxy_(1, vi)* my::vderxy_(1, 1)
                           +     my::derxy_(2, vi)*(my::vderxy_(1, 2) + my::vderxy_(2, 1)) ;
        double visres2 =         my::derxy_(0, vi)*(my::vderxy_(0, 2) + my::vderxy_(2, 0))
                           +     my::derxy_(1, vi)*(my::vderxy_(1, 2) + my::vderxy_(2, 1))
                           + 2.0*my::derxy_(2, vi)* my::vderxy_(2, 2) ;
        emesh(vi*4 + 0, ui*3 + 0) += derinvJ0*visres0;
        emesh(vi*4 + 1, ui*3 + 0) += derinvJ0*visres1;
        emesh(vi*4 + 2, ui*3 + 0) += derinvJ0*visres2;

        emesh(vi*4 + 0, ui*3 + 1) += derinvJ1*visres0;
        emesh(vi*4 + 1, ui*3 + 1) += derinvJ1*visres1;
        emesh(vi*4 + 2, ui*3 + 1) += derinvJ1*visres2;

        emesh(vi*4 + 0, ui*3 + 2) += derinvJ2*visres0;
        emesh(vi*4 + 1, ui*3 + 2) += derinvJ2*visres1;
        emesh(vi*4 + 2, ui*3 + 2) += derinvJ2*visres2;

        double visres0_poro =     2.0*gradphi(0)*my::funct_(vi)* my::vderxy_(0, 0)
                                +     gradphi(1)*my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                +     gradphi(2)*my::funct_(vi)*(my::vderxy_(0, 2) + my::vderxy_(2, 0)) ;
        double visres1_poro =         gradphi(0)*my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                + 2.0*gradphi(1)*my::funct_(vi)* my::vderxy_(1, 1)
                                +     gradphi(2)*my::funct_(vi)*(my::vderxy_(1, 2) + my::vderxy_(2, 1)) ;
        double visres2_poro =         gradphi(0)*my::funct_(vi)*(my::vderxy_(0, 2) + my::vderxy_(2, 0))
                                +     gradphi(1)*my::funct_(vi)*(my::vderxy_(1, 2) + my::vderxy_(2, 1))
                                + 2.0*gradphi(2)*my::funct_(vi)* my::vderxy_(2, 2) ;

        emesh(vi*4 + 0, ui*3 + 0) += derinvJ0/porosity*visres0_poro;
        emesh(vi*4 + 1, ui*3 + 0) += derinvJ0/porosity*visres1_poro;
        emesh(vi*4 + 2, ui*3 + 0) += derinvJ0/porosity*visres2_poro;

        emesh(vi*4 + 0, ui*3 + 1) += derinvJ1/porosity*visres0_poro;
        emesh(vi*4 + 1, ui*3 + 1) += derinvJ1/porosity*visres1_poro;
        emesh(vi*4 + 2, ui*3 + 1) += derinvJ1/porosity*visres2_poro;

        emesh(vi*4 + 0, ui*3 + 2) += derinvJ2/porosity*visres0_poro;
        emesh(vi*4 + 1, ui*3 + 2) += derinvJ2/porosity*visres1_poro;
        emesh(vi*4 + 2, ui*3 + 2) += derinvJ2/porosity*visres2_poro;

        double v0_poro =    // 2.0*my::funct_(vi)* my::vderxy_(0, 0)
                           +     my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                                                               * (   gradphi(0) * derxjm_(0,0,1,ui)
                                                                                   + gradphi(1) * derxjm_(0,1,1,ui)
                                                                                   + gradphi(2) * derxjm_(0,2,1,ui)
                                                                                 )
                           +     my::funct_(vi)*(my::vderxy_(0, 2) + my::vderxy_(2, 0))
                                                                               * (   gradphi(0) * derxjm_(0,0,2,ui)
                                                                                   + gradphi(1) * derxjm_(0,1,2,ui)
                                                                                   + gradphi(2) * derxjm_(0,2,2,ui));
        double v1_poro =         my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                                                               * (   gradphi(0) * derxjm_(1,0,0,ui)
                                                                                   + gradphi(1) * derxjm_(1,1,0,ui)
                                                                                   + gradphi(2) * derxjm_(1,2,0,ui))
                          // + 2.0*gradphi(1)*my::funct_(vi)* my::vderxy_(1, 1)
                           +     my::funct_(vi)*(my::vderxy_(1, 2) + my::vderxy_(2, 1))
                                                                               * (   gradphi(0) * derxjm_(1,0,2,ui)
                                                                                   + gradphi(1) * derxjm_(1,1,2,ui)
                                                                                   + gradphi(2) * derxjm_(1,2,2,ui));
        double v2_poro =         my::funct_(vi)*(my::vderxy_(0, 2) + my::vderxy_(2, 0))
                                                                               * (   gradphi(0) * derxjm_(2,0,0,ui)
                                                                                   + gradphi(1) * derxjm_(2,1,0,ui)
                                                                                   + gradphi(2) * derxjm_(2,2,0,ui))
                           +     my::funct_(vi)*(my::vderxy_(1, 2) + my::vderxy_(2, 1))
                                                                               * (   gradphi(0) * derxjm_(2,0,1,ui)
                                                                                   + gradphi(1) * derxjm_(2,1,1,ui)
                                                                                   + gradphi(2) * derxjm_(2,2,1,ui));
                          // + 2.0*gradphi(2)*my::funct_(vi)* my::vderxy_(2, 2) ;

        emesh(vi * 4 + 3, ui * 3 + 0) += v/porosity * v0_poro;
        emesh(vi * 4 + 3, ui * 3 + 1) += v/porosity * v1_poro;
        emesh(vi * 4 + 3, ui * 3 + 2) += v/porosity * v2_poro;
      }
    }

    // part 2: derivative of viscosity residual

     v = timefacfac*visceff_/my::det_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double v0 = - my::vderiv_(0,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
                  - my::vderiv_(0,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
                  - my::vderiv_(0,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
                  - my::vderiv_(1,0)*(derxjm_100(ui)*xji_00)
                  - my::vderiv_(1,1)*(derxjm_100(ui)*xji_01)
                  - my::vderiv_(1,2)*(derxjm_100(ui)*xji_02)
                  - my::vderiv_(2,0)*(derxjm_200(ui)*xji_00)
                  - my::vderiv_(2,1)*(derxjm_200(ui)*xji_01)
                  - my::vderiv_(2,2)*(derxjm_200(ui)*xji_02);
      double v1 = - my::vderiv_(0,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
                  - my::vderiv_(0,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
                  - my::vderiv_(0,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
                  - my::vderiv_(1,0)*(derxjm_110(ui)*xji_00)
                  - my::vderiv_(1,1)*(derxjm_110(ui)*xji_01)
                  - my::vderiv_(1,2)*(derxjm_110(ui)*xji_02)
                  - my::vderiv_(2,0)*(derxjm_210(ui)*xji_00)
                  - my::vderiv_(2,1)*(derxjm_210(ui)*xji_01)
                  - my::vderiv_(2,2)*(derxjm_210(ui)*xji_02);
      double v2 = - my::vderiv_(0,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
                  - my::vderiv_(0,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
                  - my::vderiv_(0,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
                  - my::vderiv_(1,0)*(derxjm_120(ui)*xji_00)
                  - my::vderiv_(1,1)*(derxjm_120(ui)*xji_01)
                  - my::vderiv_(1,2)*(derxjm_120(ui)*xji_02)
                  - my::vderiv_(2,0)*(derxjm_220(ui)*xji_00)
                  - my::vderiv_(2,1)*(derxjm_220(ui)*xji_01)
                  - my::vderiv_(2,2)*(derxjm_220(ui)*xji_02);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 0, ui*3 + 0) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(2*derxjm_001(ui)*xji_00 + 2*derxjm_001(ui)*xji_00 + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
           - my::vderiv_(0,1)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
           - my::vderiv_(0,2)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
           - my::vderiv_(1,0)*(derxjm_001(ui)*xji_10)
           - my::vderiv_(1,1)*(derxjm_011(ui)*xji_10)
           - my::vderiv_(1,2)*(derxjm_021(ui)*xji_10)
           - my::vderiv_(2,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
           - my::vderiv_(2,1)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
           - my::vderiv_(2,2)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20);
      v1 = - my::vderiv_(0,0)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
           - my::vderiv_(0,1)*(2*derxjm_011(ui)*xji_01 + 2*derxjm_011(ui)*xji_01 + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
           - my::vderiv_(0,2)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
           - my::vderiv_(1,0)*(derxjm_001(ui)*xji_11)
           - my::vderiv_(1,1)*(derxjm_011(ui)*xji_11)
           - my::vderiv_(1,2)*(derxjm_021(ui)*xji_11)
           - my::vderiv_(2,0)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
           - my::vderiv_(2,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
           - my::vderiv_(2,2)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21);
      v2 = - my::vderiv_(0,0)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
           - my::vderiv_(0,1)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
           - my::vderiv_(0,2)*(2*derxjm_021(ui)*xji_02 + 2*derxjm_021(ui)*xji_02 + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
           - my::vderiv_(1,0)*(derxjm_001(ui)*xji_12)
           - my::vderiv_(1,1)*(derxjm_011(ui)*xji_12)
           - my::vderiv_(1,2)*(derxjm_021(ui)*xji_12)
           - my::vderiv_(2,0)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
           - my::vderiv_(2,1)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
           - my::vderiv_(2,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 0, ui*3 + 1) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(2*derxjm_002(ui)*xji_00 + 2*derxjm_002(ui)*xji_00 + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
           - my::vderiv_(0,1)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
           - my::vderiv_(0,2)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
           - my::vderiv_(1,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
           - my::vderiv_(1,1)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
           - my::vderiv_(1,2)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
           - my::vderiv_(2,0)*(derxjm_002(ui)*xji_20)
           - my::vderiv_(2,1)*(derxjm_012(ui)*xji_20)
           - my::vderiv_(2,2)*(derxjm_022(ui)*xji_20);
      v1 = - my::vderiv_(0,0)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
           - my::vderiv_(0,1)*(2*derxjm_012(ui)*xji_01 + 2*derxjm_012(ui)*xji_01 + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
           - my::vderiv_(0,2)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
           - my::vderiv_(1,0)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
           - my::vderiv_(1,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
           - my::vderiv_(1,2)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
           - my::vderiv_(2,0)*(derxjm_002(ui)*xji_21)
           - my::vderiv_(2,1)*(derxjm_012(ui)*xji_21)
           - my::vderiv_(2,2)*(derxjm_022(ui)*xji_21);
      v2 = - my::vderiv_(0,0)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
           - my::vderiv_(0,1)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
           - my::vderiv_(0,2)*(2*derxjm_022(ui)*xji_02 + 2*derxjm_022(ui)*xji_02 + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui))
           - my::vderiv_(1,0)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
           - my::vderiv_(1,1)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
           - my::vderiv_(1,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
           - my::vderiv_(2,0)*(derxjm_002(ui)*xji_22)
           - my::vderiv_(2,1)*(derxjm_012(ui)*xji_22)
           - my::vderiv_(2,2)*(derxjm_022(ui)*xji_22);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 0, ui*3 + 2) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(derxjm_100(ui)*xji_00)
           - my::vderiv_(0,1)*(derxjm_110(ui)*xji_00)
           - my::vderiv_(0,2)*(derxjm_120(ui)*xji_00)
           - my::vderiv_(1,0)*(2*xji_10*derxjm_100(ui) + 2*xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
           - my::vderiv_(1,1)*(2*xji_11*derxjm_100(ui) + 2*xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
           - my::vderiv_(1,2)*(2*xji_12*derxjm_100(ui) + 2*xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
           - my::vderiv_(2,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
           - my::vderiv_(2,1)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
           - my::vderiv_(2,2)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20);
      v1 = - my::vderiv_(0,0)*(derxjm_100(ui)*xji_01)
           - my::vderiv_(0,1)*(derxjm_110(ui)*xji_01)
           - my::vderiv_(0,2)*(derxjm_120(ui)*xji_01)
           - my::vderiv_(1,0)*(2*xji_10*derxjm_110(ui) + 2*xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
           - my::vderiv_(1,1)*(2*xji_11*derxjm_110(ui) + 2*xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
           - my::vderiv_(1,2)*(2*xji_12*derxjm_110(ui) + 2*xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
           - my::vderiv_(2,0)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
           - my::vderiv_(2,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
           - my::vderiv_(2,2)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21);
      v2 = - my::vderiv_(0,0)*(derxjm_100(ui)*xji_02)
           - my::vderiv_(0,1)*(derxjm_110(ui)*xji_02)
           - my::vderiv_(0,2)*(derxjm_120(ui)*xji_02)
           - my::vderiv_(1,0)*(2*xji_10*derxjm_120(ui) + 2*xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
           - my::vderiv_(1,1)*(2*xji_11*derxjm_120(ui) + 2*xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
           - my::vderiv_(1,2)*(2*xji_12*derxjm_120(ui) + 2*xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
           - my::vderiv_(2,0)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
           - my::vderiv_(2,1)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
           - my::vderiv_(2,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 1, ui*3 + 0) += v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(derxjm_001(ui)*xji_10)
           - my::vderiv_(0,1)*(derxjm_001(ui)*xji_11)
           - my::vderiv_(0,2)*(derxjm_001(ui)*xji_12)
           - my::vderiv_(1,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
           - my::vderiv_(1,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
           - my::vderiv_(1,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
           - my::vderiv_(2,0)*(derxjm_201(ui)*xji_10)
           - my::vderiv_(2,1)*(derxjm_201(ui)*xji_11)
           - my::vderiv_(2,2)*(derxjm_201(ui)*xji_12);
      v1 = - my::vderiv_(0,0)*(derxjm_011(ui)*xji_10)
           - my::vderiv_(0,1)*(derxjm_011(ui)*xji_11)
           - my::vderiv_(0,2)*(derxjm_011(ui)*xji_12)
           - my::vderiv_(1,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + xji_20*derxjm_211(ui) + xji_21*derxjm_201(ui))
           - my::vderiv_(1,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
           - my::vderiv_(1,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + xji_22*derxjm_211(ui) + xji_21*derxjm_221(ui))
           - my::vderiv_(2,0)*(derxjm_211(ui)*xji_10)
           - my::vderiv_(2,1)*(derxjm_211(ui)*xji_11)
           - my::vderiv_(2,2)*(derxjm_211(ui)*xji_12);
      v2 = - my::vderiv_(0,0)*(derxjm_021(ui)*xji_10)
           - my::vderiv_(0,1)*(derxjm_021(ui)*xji_11)
           - my::vderiv_(0,2)*(derxjm_021(ui)*xji_12)
           - my::vderiv_(1,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + xji_20*derxjm_221(ui) + xji_22*derxjm_201(ui))
           - my::vderiv_(1,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
           - my::vderiv_(1,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
           - my::vderiv_(2,0)*(derxjm_221(ui)*xji_10)
           - my::vderiv_(2,1)*(derxjm_221(ui)*xji_11)
           - my::vderiv_(2,2)*(derxjm_221(ui)*xji_12);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 1, ui*3 + 1) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
           - my::vderiv_(0,1)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
           - my::vderiv_(0,2)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
           - my::vderiv_(1,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + 2*xji_10*derxjm_102(ui) + 2*xji_10*derxjm_102(ui))
           - my::vderiv_(1,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + 2*xji_11*derxjm_102(ui) + 2*xji_10*derxjm_112(ui))
           - my::vderiv_(1,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + 2*xji_12*derxjm_102(ui) + 2*xji_10*derxjm_122(ui))
           - my::vderiv_(2,0)*(derxjm_102(ui)*xji_20)
           - my::vderiv_(2,1)*(derxjm_112(ui)*xji_20)
           - my::vderiv_(2,2)*(derxjm_122(ui)*xji_20);
      v1 = - my::vderiv_(0,0)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
           - my::vderiv_(0,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
           - my::vderiv_(0,2)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
           - my::vderiv_(1,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + 2*xji_10*derxjm_112(ui) + 2*xji_11*derxjm_102(ui))
           - my::vderiv_(1,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + 2*xji_11*derxjm_112(ui) + 2*xji_11*derxjm_112(ui))
           - my::vderiv_(1,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + 2*xji_12*derxjm_112(ui) + 2*xji_11*derxjm_122(ui))
           - my::vderiv_(2,0)*(derxjm_102(ui)*xji_21)
           - my::vderiv_(2,1)*(derxjm_112(ui)*xji_21)
           - my::vderiv_(2,2)*(derxjm_122(ui)*xji_21);
      v2 = - my::vderiv_(0,0)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
           - my::vderiv_(0,1)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
           - my::vderiv_(0,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
           - my::vderiv_(1,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + 2*xji_10*derxjm_122(ui) + 2*xji_12*derxjm_102(ui))
           - my::vderiv_(1,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + 2*xji_11*derxjm_122(ui) + 2*xji_12*derxjm_112(ui))
           - my::vderiv_(1,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + 2*xji_12*derxjm_122(ui) + 2*xji_12*derxjm_122(ui))
           - my::vderiv_(2,0)*(derxjm_102(ui)*xji_22)
           - my::vderiv_(2,1)*(derxjm_112(ui)*xji_22)
           - my::vderiv_(2,2)*(derxjm_122(ui)*xji_22);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 1, ui*3 + 2) +=   v * (my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v * my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(derxjm_200(ui)*xji_00)
           - my::vderiv_(0,1)*(derxjm_210(ui)*xji_00)
           - my::vderiv_(0,2)*(derxjm_220(ui)*xji_00)
           - my::vderiv_(1,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
           - my::vderiv_(1,1)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
           - my::vderiv_(1,2)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
           - my::vderiv_(2,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + 2*xji_20*derxjm_200(ui) + 2*xji_20*derxjm_200(ui))
           - my::vderiv_(2,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + 2*xji_21*derxjm_200(ui) + 2*xji_20*derxjm_210(ui))
           - my::vderiv_(2,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + 2*xji_22*derxjm_200(ui) + 2*xji_20*derxjm_220(ui));
      v1 = - my::vderiv_(0,0)*(derxjm_200(ui)*xji_01)
           - my::vderiv_(0,1)*(derxjm_210(ui)*xji_01)
           - my::vderiv_(0,2)*(derxjm_220(ui)*xji_01)
           - my::vderiv_(1,0)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
           - my::vderiv_(1,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
           - my::vderiv_(1,2)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
           - my::vderiv_(2,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + 2*xji_20*derxjm_210(ui) + 2*xji_21*derxjm_200(ui))
           - my::vderiv_(2,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + 2*xji_21*derxjm_210(ui) + 2*xji_21*derxjm_210(ui))
           - my::vderiv_(2,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + 2*xji_22*derxjm_210(ui) + 2*xji_21*derxjm_220(ui));
      v2 = - my::vderiv_(0,0)*(derxjm_200(ui)*xji_02)
           - my::vderiv_(0,1)*(derxjm_210(ui)*xji_02)
           - my::vderiv_(0,2)*(derxjm_220(ui)*xji_02)
           - my::vderiv_(1,0)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20)
           - my::vderiv_(1,1)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21)
           - my::vderiv_(1,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22)
           - my::vderiv_(2,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + 2*xji_20*derxjm_220(ui) + 2*xji_22*derxjm_200(ui))
           - my::vderiv_(2,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + 2*xji_21*derxjm_220(ui) + 2*xji_22*derxjm_210(ui))
           - my::vderiv_(2,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + 2*xji_22*derxjm_220(ui) + 2*xji_22*derxjm_220(ui));

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 2, ui*3 + 0) += v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                         + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
           - my::vderiv_(0,1)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
           - my::vderiv_(0,2)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
           - my::vderiv_(1,0)*(derxjm_201(ui)*xji_10)
           - my::vderiv_(1,1)*(derxjm_211(ui)*xji_10)
           - my::vderiv_(1,2)*(derxjm_221(ui)*xji_10)
           - my::vderiv_(2,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + 2*xji_20*derxjm_201(ui) + 2*xji_20*derxjm_201(ui))
           - my::vderiv_(2,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + 2*xji_21*derxjm_201(ui) + 2*xji_20*derxjm_211(ui))
           - my::vderiv_(2,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + 2*xji_22*derxjm_201(ui) + 2*xji_20*derxjm_221(ui));
      v1 = - my::vderiv_(0,0)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
           - my::vderiv_(0,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
           - my::vderiv_(0,2)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
           - my::vderiv_(1,0)*(derxjm_201(ui)*xji_11)
           - my::vderiv_(1,1)*(derxjm_211(ui)*xji_11)
           - my::vderiv_(1,2)*(derxjm_221(ui)*xji_11)
           - my::vderiv_(2,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + 2*xji_20*derxjm_211(ui) + 2*xji_21*derxjm_201(ui))
           - my::vderiv_(2,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + 2*xji_21*derxjm_211(ui) + 2*xji_21*derxjm_211(ui))
           - my::vderiv_(2,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + 2*xji_22*derxjm_211(ui) + 2*xji_21*derxjm_221(ui));
      v2 = - my::vderiv_(0,0)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20)
           - my::vderiv_(0,1)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21)
           - my::vderiv_(0,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22)
           - my::vderiv_(1,0)*(derxjm_201(ui)*xji_12)
           - my::vderiv_(1,1)*(derxjm_211(ui)*xji_12)
           - my::vderiv_(1,2)*(derxjm_221(ui)*xji_12)
           - my::vderiv_(2,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + 2*xji_20*derxjm_221(ui) + 2*xji_22*derxjm_201(ui))
           - my::vderiv_(2,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + 2*xji_21*derxjm_221(ui) + 2*xji_22*derxjm_211(ui))
           - my::vderiv_(2,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + 2*xji_22*derxjm_221(ui) + 2*xji_22*derxjm_221(ui));

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 2, ui*3 + 1) += v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                         + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(derxjm_002(ui)*xji_20)
           - my::vderiv_(0,1)*(derxjm_002(ui)*xji_21)
           - my::vderiv_(0,2)*(derxjm_002(ui)*xji_22)
           - my::vderiv_(1,0)*(derxjm_102(ui)*xji_20)
           - my::vderiv_(1,1)*(derxjm_102(ui)*xji_21)
           - my::vderiv_(1,2)*(derxjm_102(ui)*xji_22)
           - my::vderiv_(2,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
           - my::vderiv_(2,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
           - my::vderiv_(2,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui));
      v1 = - my::vderiv_(0,0)*(derxjm_012(ui)*xji_20)
           - my::vderiv_(0,1)*(derxjm_012(ui)*xji_21)
           - my::vderiv_(0,2)*(derxjm_012(ui)*xji_22)
           - my::vderiv_(1,0)*(derxjm_112(ui)*xji_20)
           - my::vderiv_(1,1)*(derxjm_112(ui)*xji_21)
           - my::vderiv_(1,2)*(derxjm_112(ui)*xji_22)
           - my::vderiv_(2,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + xji_10*derxjm_112(ui) + xji_11*derxjm_102(ui))
           - my::vderiv_(2,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
           - my::vderiv_(2,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + xji_12*derxjm_112(ui) + xji_11*derxjm_122(ui));
      v2 = - my::vderiv_(0,0)*(derxjm_022(ui)*xji_20)
           - my::vderiv_(0,1)*(derxjm_022(ui)*xji_21)
           - my::vderiv_(0,2)*(derxjm_022(ui)*xji_22)
           - my::vderiv_(1,0)*(derxjm_122(ui)*xji_20)
           - my::vderiv_(1,1)*(derxjm_122(ui)*xji_21)
           - my::vderiv_(1,2)*(derxjm_122(ui)*xji_22)
           - my::vderiv_(2,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + xji_10*derxjm_122(ui) + xji_12*derxjm_102(ui))
           - my::vderiv_(2,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
           - my::vderiv_(2,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui));

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*4 + 2, ui*3 + 2) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 + my::deriv_(2,vi)*v2)
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 + gradphi(2)*v2);
      }
    }
  }//if(visceff_)
    //*************************** linearisation of mesh motion in continuity equation**********************************

  if( my::fldpara_->PoroContiPartInt() == false )
  {
    // (porosity)*div u
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = timefacfac / my::det_ * my::funct_(vi, 0) * porosity;
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        emesh(vi * 4 + 3, ui * 3    ) += v * (+ my::vderiv_(1, 0) * derxjm_(0,0,1,ui)
                                              + my::vderiv_(1, 1) * derxjm_(0,1,1,ui)
                                              + my::vderiv_(1, 2) * derxjm_(0,2,1,ui)
                                              + my::vderiv_(2, 0) * derxjm_(0,0,2,ui)
                                              + my::vderiv_(2, 1) * derxjm_(0,1,2,ui)
                                              + my::vderiv_(2, 2) * derxjm_(0,2,2,ui));

        emesh(vi * 4 + 3, ui * 3 + 1) += v * (+ my::vderiv_(0, 0) * derxjm_(1,0,0,ui)
                                              + my::vderiv_(0, 1) * derxjm_(1,1,0,ui)
                                              + my::vderiv_(0, 2) * derxjm_(1,2,0,ui)
                                              + my::vderiv_(2, 0) * derxjm_(1,0,2,ui)
                                              + my::vderiv_(2, 1) * derxjm_(1,1,2,ui)
                                              + my::vderiv_(2, 2) * derxjm_(1,2,2,ui));

        emesh(vi * 4 + 3, ui * 3 + 2) += v * (+ my::vderiv_(0, 0) * derxjm_(2,0,0,ui)
                                              + my::vderiv_(0, 1) * derxjm_(2,1,0,ui)
                                              + my::vderiv_(0, 2) * derxjm_(2,2,0,ui)
                                              + my::vderiv_(1, 0) * derxjm_(2,0,1,ui)
                                              + my::vderiv_(1, 1) * derxjm_(2,1,1,ui)
                                              + my::vderiv_(1, 2) * derxjm_(2,2,1,ui));
      }
    }

    LINALG::Matrix<my::nsd_, my::nsd_> gridvderiv;
    gridvderiv.MultiplyNT(egridv, my::deriv_);

    // (dphi_dJ*J)*div vs
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = timefacfac / my::det_ * my::funct_(vi, 0) * dphi_dJ * J;
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        emesh(vi * 4 + 3, ui * 3 + 0) += v * (+ gridvderiv(1, 0) * derxjm_(0,0,1,ui)
                                              + gridvderiv(1, 1) * derxjm_(0,1,1,ui)
                                              + gridvderiv(1, 2) * derxjm_(0,2,1,ui)
                                              + gridvderiv(2, 0) * derxjm_(0,0,2,ui)
                                              + gridvderiv(2, 1) * derxjm_(0,1,2,ui)
                                              + gridvderiv(2, 2) * derxjm_(0,2,2,ui));

        emesh(vi * 4 + 3, ui * 3 + 1) += v * (+ gridvderiv(0, 0) * derxjm_(1,0,0,ui)
                                              + gridvderiv(0, 1) * derxjm_(1,1,0,ui)
                                              + gridvderiv(0, 2) * derxjm_(1,2,0,ui)
                                              + gridvderiv(2, 0) * derxjm_(1,0,2,ui)
                                              + gridvderiv(2, 1) * derxjm_(1,1,2,ui)
                                              + gridvderiv(2, 2) * derxjm_(1,2,2,ui));

        emesh(vi * 4 + 3, ui * 3 + 2) += v * (+ gridvderiv(0, 0) * derxjm_(2,0,0,ui)
                                              + gridvderiv(0, 1) * derxjm_(2,1,0,ui)
                                              + gridvderiv(0, 2) * derxjm_(2,2,0,ui)
                                              + gridvderiv(1, 0) * derxjm_(2,0,1,ui)
                                              + gridvderiv(1, 1) * derxjm_(2,1,1,ui)
                                              + gridvderiv(1, 2) * derxjm_(2,2,1,ui));
      }
    }

    //-----------(u-vs)grad(phi) = (u-vs)dphi/dp gradp + (u-vs)dphi/dJ gradJ
    // (u-vs)dphi/dp gradp

    for (int ui = 0; ui < my::nen_; ++ui)
    {
      double v00 = + (my::velint_(1) - gridvelint(1)) * (   my::gradp_(0) * derxjm_(0,0,1,ui)
                                                          + my::gradp_(1) * derxjm_(0,1,1,ui)
                                                          + my::gradp_(2) * derxjm_(0,2,1,ui) )
                   + (my::velint_(2) - gridvelint(2)) * (   my::gradp_(0) * derxjm_(0,0,2,ui)
                                                          + my::gradp_(1) * derxjm_(0,1,2,ui)
                                                          + my::gradp_(2) * derxjm_(0,2,2,ui));
      double v01 = + (my::velint_(0) - gridvelint(0)) * (   my::gradp_(0) * derxjm_(1,0,0,ui)
                                                          + my::gradp_(1) * derxjm_(1,1,0,ui)
                                                          + my::gradp_(2) * derxjm_(1,2,0,ui))
                   + (my::velint_(2) - gridvelint(2)) * (   my::gradp_(0) * derxjm_(1,0,2,ui)
                                                          + my::gradp_(1) * derxjm_(1,1,2,ui)
                                                          + my::gradp_(2) * derxjm_(1,2,2,ui));
      double v02 = + (my::velint_(0) - gridvelint(0)) * (   my::gradp_(0) * derxjm_(2,0,0,ui)
                                                          + my::gradp_(1) * derxjm_(2,1,0,ui)
                                                          + my::gradp_(2) * derxjm_(2,2,0,ui))
                   + (my::velint_(1) - gridvelint(1)) * (   my::gradp_(0) * derxjm_(2,0,1,ui)
                                                          + my::gradp_(1) * derxjm_(2,1,1,ui)
                                                          + my::gradp_(2) * derxjm_(2,2,1,ui));

      for (int vi = 0; vi < my::nen_; ++vi)
      {
        double v = timefacfac / my::det_ * my::funct_(vi) * dphi_dp;

        emesh(vi * 4 + 3, ui * 3 + 0) += v * v00;
        emesh(vi * 4 + 3, ui * 3 + 1) += v * v01;
        emesh(vi * 4 + 3, ui * 3 + 2) += v * v02;
      }
    }

    // (u-vs)dphi/dJ gradJ
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      double v00 = + (my::velint_(1) - gridvelint(1)) * (   gradJ(0) * derxjm_(0,0,1,ui)
                                                          + gradJ(1) * derxjm_(0,1,1,ui)
                                                          + gradJ(2) * derxjm_(0,2,1,ui))
                   + (my::velint_(2) - gridvelint(2)) * (   gradJ(0) * derxjm_(0,0,2,ui)
                                                          + gradJ(1) * derxjm_(0,1,2,ui)
                                                          + gradJ(2) * derxjm_(0,2,2,ui));
      double v01 = + (my::velint_(0) - gridvelint(0)) * (   gradJ(0) * derxjm_(1,0,0,ui)
                                                          + gradJ(1) * derxjm_(1,1,0,ui)
                                                          + gradJ(2) * derxjm_(1,2,0,ui))
                   + (my::velint_(2) - gridvelint(2)) * (   gradJ(0) * derxjm_(1,0,2,ui)
                                                          + gradJ(1) * derxjm_(1,1,2,ui)
                                                          + gradJ(2) * derxjm_(1,2,2,ui));
      double v02 = + (my::velint_(0) - gridvelint(0)) * (   gradJ(0) * derxjm_(2,0,0,ui)
                                                          + gradJ(1) * derxjm_(2,1,0,ui)
                                                          + gradJ(2) * derxjm_(2,2,0,ui))
                   + (my::velint_(1) - gridvelint(1)) * (   gradJ(0) * derxjm_(2,0,1,ui)
                                                          + gradJ(1) * derxjm_(2,1,1,ui)
                                                          + gradJ(2) * derxjm_(2,2,1,ui));

      for (int vi = 0; vi < my::nen_; ++vi)
      {
        double v = timefacfac / my::det_ * my::funct_(vi) * dphi_dJ;

        emesh(vi * 4 + 3, ui * 3 + 0) += v * v00;
        emesh(vi * 4 + 3, ui * 3 + 1) += v * v01;
        emesh(vi * 4 + 3, ui * 3 + 2) += v * v02;
      }
    }
  }
  else
  {
    LINALG::Matrix<my::nsd_, my::nsd_> gridvderiv;
    gridvderiv.MultiplyNT(egridv, my::deriv_);

    // (dphi_dJ*J+phi)*div vs
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = timefacfac / my::det_ * my::funct_(vi, 0) * (dphi_dJ * J+porosity);
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        emesh(vi * 4 + 3, ui * 3 + 0) += v * (+ gridvderiv(1, 0) * derxjm_(0,0,1,ui)
                                              + gridvderiv(1, 1) * derxjm_(0,1,1,ui)
                                              + gridvderiv(1, 2) * derxjm_(0,2,1,ui)
                                              + gridvderiv(2, 0) * derxjm_(0,0,2,ui)
                                              + gridvderiv(2, 1) * derxjm_(0,1,2,ui)
                                              + gridvderiv(2, 2) * derxjm_(0,2,2,ui));

        emesh(vi * 4 + 3, ui * 3 + 1) += v * (+ gridvderiv(0, 0) * derxjm_(1,0,0,ui)
                                              + gridvderiv(0, 1) * derxjm_(1,1,0,ui)
                                              + gridvderiv(0, 2) * derxjm_(1,2,0,ui)
                                              + gridvderiv(2, 0) * derxjm_(1,0,2,ui)
                                              + gridvderiv(2, 1) * derxjm_(1,1,2,ui)
                                              + gridvderiv(2, 2) * derxjm_(1,2,2,ui));

        emesh(vi * 4 + 3, ui * 3 + 2) += v * (+ gridvderiv(0, 0) * derxjm_(2,0,0,ui)
                                              + gridvderiv(0, 1) * derxjm_(2,1,0,ui)
                                              + gridvderiv(0, 2) * derxjm_(2,2,0,ui)
                                              + gridvderiv(1, 0) * derxjm_(2,0,1,ui)
                                              + gridvderiv(1, 1) * derxjm_(2,1,1,ui)
                                              + gridvderiv(1, 2) * derxjm_(2,2,1,ui));
      }
    }

    //----------- phi * (u-vs)grad(vi)

    for (int ui = 0; ui < my::nen_; ++ui)
    {
      double v00 = + (my::velint_(1) - gridvelint(1)) * (   my::derxy_(0,ui) * derxjm_(0,0,1,ui)
                                                          + my::derxy_(1,ui) * derxjm_(0,1,1,ui)
                                                          + my::derxy_(2,ui) * derxjm_(0,2,1,ui) )
                   + (my::velint_(2) - gridvelint(2)) * (   my::derxy_(0,ui) * derxjm_(0,0,2,ui)
                                                          + my::derxy_(1,ui) * derxjm_(0,1,2,ui)
                                                          + my::derxy_(2,ui) * derxjm_(0,2,2,ui));
      double v01 = + (my::velint_(0) - gridvelint(0)) * (   my::derxy_(0,ui) * derxjm_(1,0,0,ui)
                                                          + my::derxy_(1,ui) * derxjm_(1,1,0,ui)
                                                          + my::derxy_(2,ui) * derxjm_(1,2,0,ui))
                   + (my::velint_(2) - gridvelint(2)) * (   my::derxy_(0,ui) * derxjm_(1,0,2,ui)
                                                          + my::derxy_(1,ui) * derxjm_(1,1,2,ui)
                                                          + my::derxy_(2,ui) * derxjm_(1,2,2,ui));
      double v02 = + (my::velint_(0) - gridvelint(0)) * (   my::derxy_(0,ui) * derxjm_(2,0,0,ui)
                                                          + my::derxy_(1,ui) * derxjm_(2,1,0,ui)
                                                          + my::derxy_(2,ui) * derxjm_(2,2,0,ui))
                   + (my::velint_(1) - gridvelint(1)) * (   my::derxy_(0,ui) * derxjm_(2,0,1,ui)
                                                          + my::derxy_(1,ui) * derxjm_(2,1,1,ui)
                                                          + my::derxy_(2,ui) * derxjm_(2,2,1,ui));

      for (int vi = 0; vi < my::nen_; ++vi)
      {
        double v = timefacfac / my::det_ * porosity;

        emesh(vi * 4 + 3, ui * 3 + 0) += v * v00;
        emesh(vi * 4 + 3, ui * 3 + 1) += v * v01;
        emesh(vi * 4 + 3, ui * 3 + 2) += v * v02;
      }
    }

  }//partial integration

  // dphi_dp*dp/dt
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = timefacfac * my::funct_(vi, 0) * dphi_dp * press_dot;
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 4 + 3, ui * 3)     += v * my::derxy_(0, ui);
      emesh(vi * 4 + 3, ui * 3 + 1) += v * my::derxy_(1, ui);
      emesh(vi * 4 + 3, ui * 3 + 2) += v * my::derxy_(2, ui);
    }
  }
  //-------------------

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::LinMeshMotion_2D_OD(
    LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_) * my::nen_>& emesh,
    const LINALG::Matrix<my::nsd_, my::nen_>&                         evelaf,
    const LINALG::Matrix<my::nsd_, my::nen_>&                         egridv,
    const double &                                                    press,
    const double &                                                    press_dot,
    const double &                                                    porosity,
    const double &                                                    dphi_dp,
    const double &                                                    dphi_dJ,
    const double &                                                    J,
    LINALG::Matrix<1, my::nsd_>&                                      gradJ,
    const double & timefac, const double &                            timefacfac)
{

  //*************************** linearisation of mesh motion in momentum balance**********************************
  // mass + rhs

  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = my::fac_ * my::funct_(vi, 0);
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 3    , ui * 2    ) += v * (my::velint_(0) - my::rhsmom_(0) * my::fldpara_->Dt()
                                                     * my::fldpara_->Theta()) * my::derxy_(0, ui);
      emesh(vi * 3    , ui * 2 + 1) += v * (my::velint_(0) - my::rhsmom_(0)
                                                     * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(1, ui);

      emesh(vi * 3 + 1, ui * 2    ) += v * (my::velint_(1) - my::rhsmom_(1)
                                           * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(0, ui);
      emesh(vi * 3 + 1, ui * 2 + 1) += v * (my::velint_(1) - my::rhsmom_(1)
                                           * my::fldpara_->Dt() * my::fldpara_->Theta()) * my::derxy_(1, ui);
    }
  }

  LINALG::Matrix<my::nsd_, 1> gridvelint;
  gridvelint.Multiply(egridv, my::funct_);

  //---------reaction term (darcy term)
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = timefacfac * my::funct_(vi, 0) * my::reacoeff_;
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 3    , ui * 2    ) += v * (my::velint_(0) - gridvelint(0)) * my::derxy_(0, ui);
      emesh(vi * 3    , ui * 2 + 1) += v * (my::velint_(0) - gridvelint(0)) * my::derxy_(1, ui);

      emesh(vi * 3 + 1, ui * 2    ) += v * (my::velint_(1) - gridvelint(1)) * my::derxy_(0,ui);
      emesh(vi * 3 + 1, ui * 2 + 1) += v * (my::velint_(1) - gridvelint(1)) * my::derxy_(1, ui);
    }
  }

  //---------------convective term

  //vmy::deriv_  = sum(evelaf(i,k) * my::deriv_(j,k), k);
  my::vderiv_.MultiplyNT(evelaf, my::deriv_);
  my::convvelint_.Multiply(-1.0, egridv, my::funct_, 0.0);

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int tvi  = 3*vi;
    const int tvip = tvi+1;
    const double v = my::densaf_*timefacfac/my::det_*my::funct_(vi);
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int tui  = 2*ui;
      const int tuip = tui+1;

      emesh(tvi , tui ) += v*(
      + my::convvelint_(1)*(-my::vderiv_(0, 0)*my::deriv_(1,ui) + my::vderiv_(0, 1)*my::deriv_(0,ui))
      );

      emesh(tvi , tuip) += v*(
      + my::convvelint_(0)*(-my::vderiv_(0, 0)*my::deriv_(1,ui) + my::vderiv_(0, 1)*my::deriv_(0,ui))
      );

      emesh(tvip, tui ) += v*(
      + my::convvelint_(1)*(-my::vderiv_(1, 0)*my::deriv_(1,ui) + my::vderiv_(1, 1)*my::deriv_(0,ui))
      );

      emesh(tvip, tuip) += v*(
      + my::convvelint_(0)*(-my::vderiv_(1, 0)*my::deriv_(1,ui) + my::vderiv_(1, 1)*my::deriv_(0,ui))
      );
    }
  }

  // pressure
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int tvi  = 3*vi;
    const int tvip = tvi+1;
    const double v = press*timefacfac/my::det_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int tui = 2*ui;
      emesh(tvi,  tui + 1) += v*(my::deriv_(0, vi)*my::deriv_(1, ui) - my::deriv_(0, ui)*my::deriv_(1, vi)) ;
      emesh(tvip, tui    ) += v*(my::deriv_(0, vi)*my::deriv_(1, ui) - my::deriv_(0, ui)*my::deriv_(1, vi)) ;
    }
  }

  // //---------viscous term (brinkman term)

  if(visceff_)
  {
    //-----------grad(phi) = dphi/dp gradp + dphi/dJ gradJ
    LINALG::Matrix<my::nsd_, 1> gradphi;
    for(int ui=0 ; ui<my::nsd_;++ui)
      gradphi(ui) = dphi_dp*my::gradp_(ui)+dphi_dJ*gradJ(ui);

    // part 1: derivative of 1/det

    double v = visceff_*timefac*my::fac_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double derinvJ0 = -v*(my::deriv_(0,ui)*my::xji_(0,0) + my::deriv_(1,ui)*my::xji_(0,1) );
      double derinvJ1 = -v*(my::deriv_(0,ui)*my::xji_(1,0) + my::deriv_(1,ui)*my::xji_(1,1) );
      for (int vi=0; vi<my::nen_; ++vi)
      {
        double visres0 =   2.0*my::derxy_(0, vi)* my::vderxy_(0, 0)
                           +     my::derxy_(1, vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0)) ;
        double visres1 =         my::derxy_(0, vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                           + 2.0*my::derxy_(1, vi)* my::vderxy_(1, 1) ;

        emesh(vi*3    , ui*2    ) += derinvJ0*visres0;
        emesh(vi*3 + 1, ui*2    ) += derinvJ0*visres1;

        emesh(vi*3    , ui*2 + 1) += derinvJ1*visres0;
        emesh(vi*3 + 1, ui*2 + 1) += derinvJ1*visres1;

        double visres0_poro =     2.0*gradphi(0)*my::funct_(vi)* my::vderxy_(0, 0)
                                +     gradphi(1)*my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0));
        double visres1_poro =         gradphi(0)*my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                + 2.0*gradphi(1)*my::funct_(vi)* my::vderxy_(1, 1) ;

        emesh(vi*3 + 0, ui*2 + 0) += derinvJ0/porosity*visres0_poro;
        emesh(vi*3 + 1, ui*2 + 0) += derinvJ0/porosity*visres1_poro;

        emesh(vi*3 + 0, ui*2 + 1) += derinvJ1/porosity*visres0_poro;
        emesh(vi*3 + 1, ui*2 + 1) += derinvJ1/porosity*visres1_poro;

        double v0_poro =   +     my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                                                               * (   gradphi(0) * my::deriv_(1,ui)
                                                                                   + gradphi(1) * my::deriv_(0,ui)
                                                                                 );
        double v1_poro =         my::funct_(vi)*(my::vderxy_(0, 1) + my::vderxy_(1, 0))
                                                                               * (   gradphi(0) * my::deriv_(1,ui)
                                                                                   + gradphi(1) * my::deriv_(0,ui)
                                                                                 );

        emesh(vi * 3 + 2, ui * 2 + 0) += v/porosity * v0_poro;
        emesh(vi * 3 + 2, ui * 2 + 1) += v/porosity * v1_poro;
      }
    }

    // part 2: derivative of viscosity residual

     v = timefacfac*visceff_/my::det_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double v0 = - my::vderiv_(0,0)*(my::xji_(1,0)*my::deriv_(1,ui) + my::xji_(1,0)*my::deriv_(1,ui) )
                  - my::vderiv_(0,1)*(my::xji_(1,1)*my::deriv_(1,ui) + my::xji_(1,0)*my::deriv_(0,ui) )
                  - my::vderiv_(1,0)*(my::deriv_(1,ui)*my::xji_(0,0))
                  - my::vderiv_(1,1)*(my::deriv_(1,ui)*my::xji_(0,1));
      double v1 = - my::vderiv_(0,0)*(my::xji_(1,0)*my::deriv_(0,ui) + my::xji_(1,1)*my::deriv_(1,ui) )
                  - my::vderiv_(0,1)*(my::xji_(1,1)*my::deriv_(0,ui) + my::xji_(1,1)*my::deriv_(0,ui) )
                  - my::vderiv_(1,0)*(my::deriv_(0,ui)*my::xji_(0,0))
                  - my::vderiv_(1,1)*(my::deriv_(0,ui)*my::xji_(0,1));

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*3 + 0, ui*2 + 0) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 )
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 );
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(2*my::deriv_(1,ui)*my::xji_(0,0) + 2*my::deriv_(1,ui)*my::xji_(0,0) )
           - my::vderiv_(0,1)*(2*my::deriv_(0,ui)*my::xji_(0,0) + 2*my::deriv_(1,ui)*my::xji_(0,1) )
           - my::vderiv_(1,0)*(my::deriv_(1,ui)*my::xji_(1,0))
           - my::vderiv_(1,1)*(my::deriv_(0,ui)*my::xji_(1,0));
      v1 = - my::vderiv_(0,0)*(2*my::deriv_(0,ui)*my::xji_(0,0) + 2*my::deriv_(1,ui)*my::xji_(0,1) )
           - my::vderiv_(0,1)*(2*my::deriv_(0,ui)*my::xji_(0,1) + 2*my::deriv_(0,ui)*my::xji_(0,1) )
           - my::vderiv_(1,0)*(my::deriv_(1,ui)*my::xji_(1,1))
           - my::vderiv_(1,1)*(my::deriv_(0,ui)*my::xji_(1,1));

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*3 + 0, ui*2 + 1) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 )
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 );
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(my::deriv_(1,ui)*my::xji_(0,0))
           - my::vderiv_(0,1)*(my::deriv_(0,ui)*my::xji_(0,0))
           - my::vderiv_(1,0)*(2*my::xji_(1,0)*my::deriv_(1,ui) + 2*my::xji_(1,0)*my::deriv_(1,ui) )
           - my::vderiv_(1,1)*(2*my::xji_(1,1)*my::deriv_(1,ui) + 2*my::xji_(1,0)*my::deriv_(0,ui) );
      v1 = - my::vderiv_(0,0)*(my::deriv_(1,ui)*xji_01)
           - my::vderiv_(0,1)*(my::deriv_(0,ui)*xji_01)
           - my::vderiv_(1,0)*(2*my::xji_(1,0)*my::deriv_(0,ui) + 2*my::xji_(1,1)*my::deriv_(1,ui) )
           - my::vderiv_(1,1)*(2*my::xji_(1,1)*my::deriv_(0,ui) + 2*my::xji_(1,1)*my::deriv_(0,ui) );

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*3 + 1, ui*2 + 0) += v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 )
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 );
      }

      ////////////////////////////////////////////////////////////////

      v0 = - my::vderiv_(0,0)*(my::deriv_(1,ui)*my::xji_(1,0))
           - my::vderiv_(0,1)*(my::deriv_(1,ui)*my::xji_(1,1))
           - my::vderiv_(1,0)*(my::xji_(0,0)*my::deriv_(1,ui) + my::xji_(0,0)*my::deriv_(1,ui) )
           - my::vderiv_(1,1)*(my::xji_(0,1)*my::deriv_(1,ui) + my::xji_(0,0)*my::deriv_(0,ui) );
      v1 = - my::vderiv_(0,0)*(my::deriv_(0,ui)*my::xji_(1,0))
           - my::vderiv_(0,1)*(my::deriv_(0,ui)*my::xji_(1,1))
           - my::vderiv_(1,0)*(my::xji_(0,0)*my::deriv_(0,ui) + my::xji_(0,1)*my::deriv_(1,ui) )
           - my::vderiv_(1,1)*(my::xji_(0,1)*my::deriv_(0,ui) + my::xji_(0,1)*my::deriv_(0,ui) );

      for (int vi=0; vi<my::nen_; ++vi)
      {
        emesh(vi*3 + 1, ui*2 + 1) +=   v*(my::deriv_(0,vi)*v0 + my::deriv_(1,vi)*v1 )
                                     + v*my::funct_(vi)/porosity*(gradphi(0)*v0 + gradphi(1)*v1 );
      }

    }
  }//if(visceff_)

    //*************************** linearisation of mesh motion in continuity equation**********************************

  // dphi_dp*dp/dt
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    double v = timefacfac * my::funct_(vi, 0) * dphi_dp * press_dot;
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      emesh(vi * 3 + 2, ui * 2    ) += v * my::derxy_(0, ui);
      emesh(vi * 3 + 2, ui * 2 + 1) += v * my::derxy_(1, ui);
    }
  }

  if( my::fldpara_->PoroContiPartInt() == false )
  {
    // (porosity)*div u
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = timefacfac / my::det_ * my::funct_(vi, 0) * porosity;
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        emesh(vi * 3 + 2, ui * 2    ) += v * (+ my::vderiv_(1, 0) * my::deriv_(1,ui)
                                              + my::vderiv_(1, 1) * my::deriv_(0,ui));

        emesh(vi * 3 + 2, ui * 2 + 1) += v * (+ my::vderiv_(0, 0) * my::deriv_(1,ui)
                                              + my::vderiv_(0, 1) * my::deriv_(0,ui));
      }
    }

    LINALG::Matrix<my::nsd_, my::nsd_> gridvderiv;
    gridvderiv.MultiplyNT(egridv, my::deriv_);

    // (dphi_dJ*J)*div vs
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = timefacfac / my::det_ * my::funct_(vi, 0) * dphi_dJ * J;
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        emesh(vi * 3 + 2, ui * 2    ) += v * (+ gridvderiv(1, 0) * my::deriv_(1,ui)
                                              + gridvderiv(1, 1) * my::deriv_(0,ui));

        emesh(vi * 3 + 2, ui * 2 + 1) += v * (+ gridvderiv(0, 0) * my::deriv_(1,ui)
                                              + gridvderiv(0, 1) * my::deriv_(0,ui));
      }
    }

    //-----------(u-vs)grad(phi) = (u-vs)dphi/dp gradp + (u-vs)dphi/dJ gradJ
    // (u-vs)dphi/dp gradp

    for (int ui = 0; ui < my::nen_; ++ui)
    {
      double v00 = + (my::velint_(1) - gridvelint(1)) * (   my::gradp_(0) * my::deriv_(1,ui)
                                                          + my::gradp_(1) * my::deriv_(0,ui) );
      double v01 = + (my::velint_(0) - gridvelint(0)) * (   my::gradp_(0) * my::deriv_(1,ui)
                                                          + my::gradp_(1) * my::deriv_(0,ui) );

      for (int vi = 0; vi < my::nen_; ++vi)
      {
        double v = timefacfac / my::det_ * my::funct_(vi) * dphi_dp;

        emesh(vi * 3 + 2, ui * 2    ) += v * v00;
        emesh(vi * 3 + 2, ui * 2 + 1) += v * v01;
      }
    }

    // (u-vs)dphi/dJ gradJ
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      double v00 = + (my::velint_(1) - gridvelint(1)) * (   gradJ(0) * my::deriv_(1,ui)
                                                          + gradJ(1) * my::deriv_(0,ui));
      double v01 = + (my::velint_(0) - gridvelint(0)) * (   gradJ(0) * my::deriv_(1,ui)
                                                          + gradJ(1) * my::deriv_(0,ui));

      for (int vi = 0; vi < my::nen_; ++vi)
      {
        double v = timefacfac / my::det_ * my::funct_(vi) * dphi_dJ;

        emesh(vi * 3 + 2, ui * 2    ) += v * v00;
        emesh(vi * 3 + 2, ui * 2 + 1) += v * v01;
      }
    }
  }
  else
  {
    LINALG::Matrix<my::nsd_, my::nsd_> gridvderiv;
    gridvderiv.MultiplyNT(egridv, my::deriv_);

    // (dphi_dJ*J+phi)*div vs
    for (int vi = 0; vi < my::nen_; ++vi)
    {
      double v = timefacfac / my::det_ * my::funct_(vi, 0) * (dphi_dJ * J+porosity);
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        emesh(vi * 3 + 2, ui * 2    ) += v * (+ gridvderiv(1, 0) * my::deriv_(1,ui)
                                              + gridvderiv(1, 1) * my::deriv_(0,ui));

        emesh(vi * 3 + 2, ui * 2 + 1) += v * (+ gridvderiv(0, 0) * my::deriv_(1,ui)
                                              + gridvderiv(0, 1) * my::deriv_(0,ui));
      }
    }

    //----------- phi * (u-vs)grad(vi)

    for (int ui = 0; ui < my::nen_; ++ui)
    {
      double v00 = + (my::velint_(1) - gridvelint(1)) * (   my::derxy_(0,ui) * my::deriv_(1,ui)
                                                          + my::derxy_(1,ui) * my::deriv_(0,ui));
      double v01 = + (my::velint_(0) - gridvelint(0)) * (   my::derxy_(0,ui) * my::deriv_(1,ui)
                                                          + my::derxy_(1,ui) * my::deriv_(0,ui));

      for (int vi = 0; vi < my::nen_; ++vi)
      {
        double v = timefacfac / my::det_ * porosity;

        emesh(vi * 3 + 2, ui * 2    ) += v * v00;
        emesh(vi * 3 + 2, ui * 2 + 1) += v * v01;
      }
    }

  }//partial integration
  //-------------------

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::PSPG(
    LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &         estif_q_u,
    LINALG::Matrix<my::nen_,my::nen_> &                   ppmat,
    LINALG::Matrix<my::nen_,1> &                          preforce,
    LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &          lin_resM_Du,
    const double &                                        fac3,
    const double &                                        timefacfac,
    const double &                                        timefacfacpre,
    const double &                                        rhsfac)
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

    if(my::fldpara_->Tds()==INPAR::FLUID::subscales_quasistatic)
    {
      scal_grad_q=my::tau_(1);
    }
    else
    {
      scal_grad_q=my::fldpara_->AlphaF()*fac3;
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
    /* pressure stabilisation: reaction if included */
    /*
              /                     \
             |                      |
             |  sigma*Du , nabla q  |
             |                      |
              \                    /
    */
    /* pressure stabilisation: viscosity (-L_visc_u) */
    /*
              /                              \
             |               /  \             |
         mu  |  nabla o eps | Du | , nabla q  |
             |               \  /             |
              \                              /
    */

    if (my::is_higher_order_ele_ || my::fldpara_->IsNewton())
    {
      for(int jdim=0;jdim<my::nsd_;++jdim)
      {
        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui_p_jdim   = my::nsd_*ui + jdim;

          for(int idim=0;idim<my::nsd_;++idim)
          {
            const int nsd_idim=my::nsd_*idim;

            for (int vi=0; vi<my::nen_; ++vi)
            {
              const double temp_vi_idim=my::derxy_(idim,vi)*scal_grad_q;

              estif_q_u(vi,fui_p_jdim) += lin_resM_Du(nsd_idim+jdim,ui)*temp_vi_idim;
            } // jdim
          } // vi
        } // ui
      } //idim
    } // end if (is_higher_order_ele_) or (newton_)
    else
    {
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for(int idim=0;idim<my::nsd_;++idim)
        {
          const int nsd_idim=my::nsd_*idim;

          const double temp_vi_idim=my::derxy_(idim, vi)*scal_grad_q;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui_p_idim   = my::nsd_*ui + idim;

            estif_q_u(vi,fui_p_idim) += lin_resM_Du(nsd_idim+idim,ui)*temp_vi_idim;
          } // vi
        } // ui
      } //idim
    } // end if not (is_higher_order_ele_) nor (newton_)


    for (int ui=0; ui<my::nen_; ++ui)
    {
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        const double v=timefacfacpre*my::derxy_(idim,ui)*scal_grad_q;

        for (int vi=0; vi<my::nen_; ++vi)
        {
          /* pressure stabilisation: pressure( L_pres_p) */
          /*
               /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
               \                    /
          */
          ppmat(vi,ui)+=v*my::derxy_(idim,vi);
        } // vi
      } // end for(idim)
    }  // ui

    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      const double temp = rhsfac*my::sgvelint_(idim);

      for (int vi=0; vi<my::nen_; ++vi)
      {
        // pressure stabilisation
        preforce(vi) += temp*my::derxy_(idim, vi);
      }
    } // end for(idim)

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoro<distype>::CalcAuxiliaryDerivatives(
                              const LINALG::Matrix<my::nsd_, my::nen_>& edispnp,
                              const LINALG::Matrix<my::nsd_,my::nsd_>& defgrd,
                              const LINALG::Matrix<my::nsd_,my::nsd_>& defgrd_inv,
                              LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>& N_X_x,
                              LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>& F_x,
                              LINALG::Matrix<my::nsd_,my::nsd_*my::nen_>& Finv_N_X_x)
{
  if(my::nsd_==3)
  {
    for (int i=0; i<my::nen_; ++i)
    {
      N_X_x(0,i) = my::derxy2_(0,i)*defgrd(0,0)+my::derxy2_(3,i)*defgrd(1,0)+my::derxy2_(4,i)*defgrd(2,0);
      N_X_x(1,i) = my::derxy2_(3,i)*defgrd(0,1)+my::derxy2_(1,i)*defgrd(1,1)+my::derxy2_(5,i)*defgrd(2,1);
      N_X_x(2,i) = my::derxy2_(4,i)*defgrd(0,2)+my::derxy2_(5,i)*defgrd(1,2)+my::derxy2_(2,i)*defgrd(2,2);

      N_X_x(3,i) = my::derxy2_(3,i)*defgrd(0,0)+my::derxy2_(1,i)*defgrd(1,0)+my::derxy2_(5,i)*defgrd(2,0);
      N_X_x(4,i) = my::derxy2_(4,i)*defgrd(0,0)+my::derxy2_(5,i)*defgrd(1,0)+my::derxy2_(2,i)*defgrd(2,0);

      N_X_x(5,i) = my::derxy2_(0,i)*defgrd(0,1)+my::derxy2_(3,i)*defgrd(1,1)+my::derxy2_(4,i)*defgrd(2,1);
      N_X_x(6,i) = my::derxy2_(4,i)*defgrd(0,1)+my::derxy2_(5,i)*defgrd(1,1)+my::derxy2_(2,i)*defgrd(2,1);

      N_X_x(7,i) = my::derxy2_(0,i)*defgrd(0,2)+my::derxy2_(3,i)*defgrd(1,2)+my::derxy2_(4,i)*defgrd(2,2);
      N_X_x(8,i) = my::derxy2_(3,i)*defgrd(0,2)+my::derxy2_(1,i)*defgrd(1,2)+my::derxy2_(5,i)*defgrd(2,2);
    }

    for(int i=0; i<my::nsd_; i++)
    {
      for(int n=0; n<my::nen_; n++)
      {
        //! second my::derivatives  are orderd as followed: (N,Xx ; N,Yy ; N,Zz ; N,Xy ; N,Xz ; N,Yx ; N,Yz;  N,Zx ; N,Zy)
        F_x(i*my::nsd_+0, 0) += N_X_x(0,n)*edispnp(i,n);
        F_x(i*my::nsd_+1, 0) += N_X_x(5,n)*edispnp(i,n);
        F_x(i*my::nsd_+2, 0) += N_X_x(7,n)*edispnp(i,n);

        F_x(i*my::nsd_+0, 1) += N_X_x(3,n)*edispnp(i,n);
        F_x(i*my::nsd_+1, 1) += N_X_x(1,n)*edispnp(i,n);
        F_x(i*my::nsd_+2, 1) += N_X_x(8,n)*edispnp(i,n);

        F_x(i*my::nsd_+0, 2) += N_X_x(4,n)*edispnp(i,n);
        F_x(i*my::nsd_+1, 2) += N_X_x(6,n)*edispnp(i,n);
        F_x(i*my::nsd_+2, 2) += N_X_x(2,n)*edispnp(i,n);
      }
    }

    //! second derivatives  are ordered as followed: (N,Xx ; N,Yy ; N,Zz ; N,Xy ; N,Xz ; N,Yx ; N,Yz;  N,Zx ; N,Zy)
    for (int n =0; n<my::nen_; n++)
    {
      int n_nsd=n*my::nsd_;
      Finv_N_X_x(0, n_nsd + 0) += defgrd_inv(0,0) * N_X_x(0,n) + defgrd_inv(1,0) * N_X_x(5,n)+ defgrd_inv(2,0) * N_X_x(7,n);
      Finv_N_X_x(0, n_nsd + 1) += defgrd_inv(0,1) * N_X_x(0,n) + defgrd_inv(1,1) * N_X_x(5,n)+ defgrd_inv(2,1) * N_X_x(7,n);
      Finv_N_X_x(0, n_nsd + 2) += defgrd_inv(0,2) * N_X_x(0,n) + defgrd_inv(1,2) * N_X_x(5,n)+ defgrd_inv(2,2) * N_X_x(7,n);

      Finv_N_X_x(1, n_nsd + 0) += defgrd_inv(0,0) * N_X_x(3,n) + defgrd_inv(1,0) * N_X_x(1,n)+ defgrd_inv(2,0) * N_X_x(8,n);
      Finv_N_X_x(1, n_nsd + 1) += defgrd_inv(0,1) * N_X_x(3,n) + defgrd_inv(1,1) * N_X_x(1,n)+ defgrd_inv(2,1) * N_X_x(8,n);
      Finv_N_X_x(1, n_nsd + 2) += defgrd_inv(0,2) * N_X_x(3,n) + defgrd_inv(1,2) * N_X_x(1,n)+ defgrd_inv(2,2) * N_X_x(8,n);

      Finv_N_X_x(2, n_nsd + 0) += defgrd_inv(0,0) * N_X_x(4,n) + defgrd_inv(1,0) * N_X_x(6,n)+ defgrd_inv(2,0) * N_X_x(2,n);
      Finv_N_X_x(2, n_nsd + 1) += defgrd_inv(0,1) * N_X_x(4,n) + defgrd_inv(1,1) * N_X_x(6,n)+ defgrd_inv(2,1) * N_X_x(2,n);
      Finv_N_X_x(2, n_nsd + 2) += defgrd_inv(0,2) * N_X_x(4,n) + defgrd_inv(1,2) * N_X_x(6,n)+ defgrd_inv(2,2) * N_X_x(2,n);
    }
  }
  else if(my::nsd_==2)
  {
    //! second my::derivatives  are orderd as followed: (N,Xx ; N,Yy ; N,Xy ; N,Yx )
    for (int i=0; i<my::nen_; ++i)
    {
      N_X_x(0,i) = my::derxy2_(0,i)*defgrd(0,0)+my::derxy2_(2,i)*defgrd(1,0);
      N_X_x(1,i) = my::derxy2_(2,i)*defgrd(0,1)+my::derxy2_(1,i)*defgrd(1,1);

      N_X_x(2,i) = my::derxy2_(2,i)*defgrd(0,0)+my::derxy2_(1,i)*defgrd(1,0);

      N_X_x(3,i) = my::derxy2_(0,i)*defgrd(0,1)+my::derxy2_(2,i)*defgrd(1,1);
    }

    for(int i=0; i<my::nsd_; i++)
    {
      for(int n=0; n<my::nen_; n++)
      {
        //! second my::derivatives  are orderd as followed: (N,Xx ; N,Yy ; N,Xy ; N,Yx )
        F_x(i*my::nsd_+0, 0) += N_X_x(0,n)*edispnp(i,n);
        F_x(i*my::nsd_+1, 0) += N_X_x(3,n)*edispnp(i,n);

        F_x(i*my::nsd_+0, 1) += N_X_x(2,n)*edispnp(i,n);
        F_x(i*my::nsd_+1, 1) += N_X_x(1,n)*edispnp(i,n);
      }
    }

    //! second derivatives  are ordered as followed: (N,Xx ; N,Yy ; N,Xy ; N,Yx )
    for (int n =0; n<my::nen_; n++)
    {
      int n_nsd=n*my::nsd_;
      Finv_N_X_x(0, n_nsd + 0) += defgrd_inv(0,0) * N_X_x(0,n) + defgrd_inv(1,0) * N_X_x(3,n);
      Finv_N_X_x(0, n_nsd + 1) += defgrd_inv(0,1) * N_X_x(0,n) + defgrd_inv(1,1) * N_X_x(3,n);

      Finv_N_X_x(1, n_nsd + 0) += defgrd_inv(0,0) * N_X_x(2,n) + defgrd_inv(1,0) * N_X_x(1,n);
      Finv_N_X_x(1, n_nsd + 1) += defgrd_inv(0,1) * N_X_x(2,n) + defgrd_inv(1,1) * N_X_x(1,n);
    }
  }
  else
    dserror("not implemented for 1D!");
  return;
}


// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcPoro<DRT::Element::nurbs27>;



