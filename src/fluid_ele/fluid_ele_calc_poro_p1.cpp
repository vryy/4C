/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of poro Fluid element (p1 poro fluid)


\level 2

*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_poro_p1.H"

#include "fluid_ele.H"
#include "fluid_ele_parameter_poro.H"
#include "nurbs_utils.H"

#include "fluidporo.H"
#include "structporo.H"

#include "fluid_rotsym_periodicbc.H"


template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP1<distype>* DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Instance(
    ::UTILS::SingletonAction action)
{
  static auto singleton_owner = ::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::FluidEleCalcPoroP1<distype>>(
            new DRT::ELEMENTS::FluidEleCalcPoroP1<distype>());
      });

  return singleton_owner.Instance(action);
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::FluidEleCalcPoroP1()
    : DRT::ELEMENTS::FluidEleCalcPoro<distype>::FluidEleCalcPoro()
{
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if (Base::isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size =
        DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, ele, Base::myknots_, Base::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return (0);
  }  // Nurbs specific stuff

  // set element id
  Base::eid_ = ele->Id();
  // get structure material
  Base::GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  Base::rotsymmpbc_->Setup(ele);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes,
  // with pressure gradient prescribed as body force included for turbulent
  // channel flow and with scatra body force included for variable-density flow
  // (evaluation at time n+alpha_F for generalized-alpha scheme,
  //  and at time n+1 otherwise)
  // ---------------------------------------------------------------------
  static LINALG::Matrix<Base::nsd_, Base::nen_> ebofoaf(true);
  ebofoaf.Clear();
  static LINALG::Matrix<Base::nsd_, Base::nen_> eprescpgaf(true);
  eprescpgaf.Clear();
  static LINALG::Matrix<Base::nen_, 1> escabofoaf(true);
  escabofoaf.Clear();
  Base::BodyForce(ele, ebofoaf, eprescpgaf, escabofoaf);

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
  static LINALG::Matrix<Base::nsd_, Base::nen_> evelaf(true);
  evelaf.Clear();
  static LINALG::Matrix<Base::nen_, 1> epreaf(true);
  epreaf.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  static LINALG::Matrix<Base::nsd_, Base::nen_> evelnp(true);
  evelnp.Clear();
  static LINALG::Matrix<Base::nen_, 1> eprenp(true);
  eprenp.Clear();
  if (FluidEleCalc<distype>::fldparatimint_->IsGenalphaNP())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &evelnp, &eprenp, "velnp");

  static LINALG::Matrix<Base::nsd_, Base::nen_> emhist(true);
  emhist.Clear();
  static LINALG::Matrix<Base::nen_, 1> echist(true);
  echist.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &emhist, &echist, "hist");

  static LINALG::Matrix<Base::nsd_, Base::nen_> eaccam(true);
  static LINALG::Matrix<Base::nen_, 1> epressam_timederiv(true);
  eaccam.Clear();
  epressam_timederiv.Clear();

  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &eaccam, &epressam_timederiv, "accam");

  static LINALG::Matrix<Base::nen_, 1> epressn_timederiv(true);
  epressn_timederiv.Clear();
  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, nullptr, &epressn_timederiv, "accn");

  static LINALG::Matrix<Base::nen_, 1> epren(true);
  epren.Clear();
  static LINALG::Matrix<Base::nsd_, Base::nen_> eveln(true);
  eveln.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &eveln, &epren, "veln");

  static LINALG::Matrix<Base::nen_, 1> epressnp_timederiv(true);
  epressnp_timederiv.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &epressnp_timederiv, "accnp");

  static LINALG::Matrix<Base::nen_, 1> escaaf(true);
  escaaf.Clear();
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &escaaf, "scaaf");

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  static LINALG::Matrix<Base::nsd_, Base::nen_> edispnp(true);
  edispnp.Clear();
  static LINALG::Matrix<Base::nsd_, Base::nen_> egridv(true);
  egridv.Clear();
  static LINALG::Matrix<Base::nsd_, Base::nen_> egridvn(true);
  egridvn.Clear();
  static LINALG::Matrix<Base::nsd_, Base::nen_> edispn(true);
  edispn.Clear();

  static LINALG::Matrix<Base::nen_, 1> eporositynp(true);
  eporositynp.Clear();
  static LINALG::Matrix<Base::nen_, 1> eporositydot(true);
  eporositydot.Clear();
  static LINALG::Matrix<Base::nen_, 1> eporositydotn(true);
  eporositydotn.Clear();

  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, &eporositynp, "dispnp");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridv, &eporositydot, "gridv");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridvn, &eporositydotn, "gridvn");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispn, nullptr, "dispn");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, Base::nsd_, LINALG::Matrix<Base::nsd_, Base::nen_>>(
      ele, Base::xyze_);

  // construct views
  LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, (Base::nsd_ + 1) * Base::nen_> elemat1(
      elemat1_epetra, true);
  LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, 1> elevec1(elevec1_epetra, true);
  // elemat2 and elevec2+3 are currently not in use

  Base::PreEvaluate(params, ele, discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = Base::Evaluate(params, ebofoaf, elemat1, elevec1, evelaf, epreaf, evelnp, eveln,
      eprenp, epren, emhist, echist, epressnp_timederiv, epressam_timederiv, epressn_timederiv,
      eaccam, edispnp, edispn, egridv, egridvn, escaaf, &eporositynp, &eporositydot, &eporositydotn,
      mat, ele->IsAle(), intpoints);

  return result;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputePorosity(Teuchos::ParameterList& params,
    const double& press, const double& J, const int& gp,
    const LINALG::Matrix<Base::nen_, 1>& shapfct, const LINALG::Matrix<Base::nen_, 1>* myporosity,
    double& porosity, double* dphi_dp, double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ,
    double* dphi_dpp, bool save)
{
  if (myporosity == nullptr)
    dserror("no porosity values given!!");
  else
    porosity = shapfct.Dot(*myporosity);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputePorosityGradient(const double& dphidp,
    const double& dphidJ, const LINALG::Matrix<Base::nsd_, 1>& gradJ,
    const LINALG::Matrix<Base::nsd_, 1>& gradp, const LINALG::Matrix<Base::nen_, 1>* eporositynp,
    LINALG::Matrix<Base::nsd_, 1>& grad_porosity, LINALG::Matrix<Base::nsd_, 1>& refgrad_porosity)
{
  if (eporositynp == nullptr)
    dserror("no porosity values given for calculation of porosity gradient!!");

  //--------------------- current porosity gradient
  grad_porosity.Multiply(Base::derxy_, *eporositynp);

  refgrad_porosity.Multiply(Base::xjm_, grad_porosity);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::EvaluatePressureEquation(
    Teuchos::ParameterList& params, const double& timefacfacpre, const double& rhsfac,
    const double& dphi_dp, const double& dphi_dJ, const double& dphi_dJdp, const double& dphi_dpp,
    const LINALG::Matrix<Base::nen_, 1>* eporositydot,
    const LINALG::Matrix<Base::nen_, 1>* eporositydotn, const LINALG::Matrix<Base::nen_, 1>& echist,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& dgradphi_dp,
    LINALG::Matrix<Base::nen_, Base::nen_ * Base::nsd_>& estif_q_u,
    LINALG::Matrix<Base::nen_, Base::nen_>& ppmat, LINALG::Matrix<Base::nen_, 1>& preforce)
{
  // first evaluate terms without porosity time derivative
  Base::EvaluatePressureEquationNonTransient(params, timefacfacpre, rhsfac, dphi_dp, dphi_dJ,
      dphi_dJdp, dphi_dpp, dgradphi_dp, estif_q_u, ppmat, preforce);

  // now the porosity time derivative (different for standard poro and poro_p1 elements)
  if (!Base::porofldpara_->IsStationaryConti())
  {
    // inertia terms on the right hand side for instationary fluids

    if (eporositydot)
    {
      double porositydot = Base::funct_.Dot(*eporositydot);
      // double porositydot =  Base::funct_.Dot(*eporositydotn);

      for (int vi = 0; vi < Base::nen_; ++vi)
      {
        // check genalpha case
        preforce(vi) -= rhsfac * porositydot * Base::funct_(vi);
      }

      // just update internal variables, no contribution to rhs
      const double porositydotn = Base::funct_.Dot(*eporositydotn);

      Base::hist_con_ = Base::fldparatimint_->OmTheta() * Base::fldparatimint_->Dt() * porositydotn;

      // rhs from last time step
      Base::rhscon_ =
          1.0 / Base::fldparatimint_->Dt() / Base::fldparatimint_->Theta() * Base::hist_con_;

      // transient part of continuity equation residual
      Base::conres_old_ += porositydot - Base::rhscon_;
    }
    else
      dserror("no porosity time derivative given for poro_p1 element!");
  }
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::EvaluateOD(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if (Base::isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size =
        DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, ele, Base::myknots_, Base::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return (0);
  }  // Nurbs specific stuff

  // set element id
  Base::eid_ = ele->Id();

  // get structure material
  Base::GetStructMaterial(ele);

  // rotationally symmetric periodic bc's: do setup for current element
  // (only required to be set up for routines "ExtractValuesFromGlobalVector")
  Base::rotsymmpbc_->Setup(ele);

  // construct views
  LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, (Base::nsd_ + 1) * Base::nen_> elemat1(
      elemat1_epetra, true);
  //  LINALG::Matrix<(Base::nsd_+1)*Base::nen_,(Base::nsd_+1)*Base::nen_>
  //  elemat2(elemat2_epetra,true);
  LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, 1> elevec1(elevec1_epetra, true);
  // elevec2 and elevec3 are currently not in use

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes,
  // with pressure gradient prescribed as body force included for turbulent
  // channel flow and with scatra body force included for variable-density flow
  // (evaluation at time n+alpha_F for generalized-alpha scheme,
  //  and at time n+1 otherwise)
  // ---------------------------------------------------------------------
  LINALG::Matrix<Base::nsd_, Base::nen_> ebofoaf(true);
  LINALG::Matrix<Base::nsd_, Base::nen_> eprescpgaf(true);
  LINALG::Matrix<Base::nen_, 1> escabofoaf(true);
  this->BodyForce(ele, ebofoaf, eprescpgaf, escabofoaf);

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
  LINALG::Matrix<Base::nsd_, Base::nen_> evelaf(true);
  LINALG::Matrix<Base::nen_, 1> epreaf(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<Base::nsd_, Base::nen_> evelnp(true);
  LINALG::Matrix<Base::nen_, 1> eprenp(true);
  if (Base::fldparatimint_->IsGenalphaNP())
    this->ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &evelnp, &eprenp, "velnp");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<Base::nsd_, Base::nen_> eveln(true);
  LINALG::Matrix<Base::nen_, 1> epren(true);
  if (Base::fldparatimint_->IsGenalphaNP())
    this->ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &eveln, &epren, "veln");

  static LINALG::Matrix<Base::nsd_, Base::nen_> eaccam(true);
  static LINALG::Matrix<Base::nen_, 1> epressam_timederiv(true);
  eaccam.Clear();
  epressam_timederiv.Clear();

  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, &eaccam, &epressam_timederiv, "accam");

  static LINALG::Matrix<Base::nen_, 1> epressn_timederiv(true);
  epressn_timederiv.Clear();
  if (Base::fldparatimint_->IsGenalpha())
    Base::ExtractValuesFromGlobalVector(
        discretization, lm, *Base::rotsymmpbc_, nullptr, &epressn_timederiv, "accn");

  LINALG::Matrix<Base::nen_, 1> epressnp_timederiv(true);
  this->ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &epressnp_timederiv, "accnp");

  LINALG::Matrix<Base::nen_, 1> escaaf(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, nullptr, &escaaf, "scaaf");

  LINALG::Matrix<Base::nsd_, Base::nen_> emhist(true);
  LINALG::Matrix<Base::nen_, 1> echist(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &emhist, &echist, "hist");

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<Base::nsd_, Base::nen_> edispnp(true);
  LINALG::Matrix<Base::nsd_, Base::nen_> egridv(true);
  LINALG::Matrix<Base::nsd_, Base::nen_> edispn(true);
  LINALG::Matrix<Base::nsd_, Base::nen_> egridvn(true);

  LINALG::Matrix<Base::nen_, 1> eporositynp(true);

  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, &eporositynp, "dispnp");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridv, nullptr, "gridv");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispn, nullptr, "dispn");
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridvn, nullptr, "gridvn");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, Base::nsd_, LINALG::Matrix<Base::nsd_, Base::nen_>>(
      ele, Base::xyze_);

  Base::PreEvaluate(params, ele, discretization);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = EvaluateOD(params, ebofoaf, elemat1, elevec1, evelaf, epreaf, evelnp, eveln, eprenp,
      epren, emhist, echist, epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam,
      edispnp, edispn, egridv, egridvn, escaaf, &eporositynp, mat, ele->IsAle(), intpoints);

  return result;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::EvaluateOD(Teuchos::ParameterList& params,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& ebofoaf,
    LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, (Base::nsd_ + 1) * Base::nen_>& elemat1,
    LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, 1>& elevec1,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& evelaf,
    const LINALG::Matrix<Base::nen_, 1>& epreaf,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& evelnp,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& eveln,
    const LINALG::Matrix<Base::nen_, 1>& eprenp, const LINALG::Matrix<Base::nen_, 1>& epren,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& emhist,
    const LINALG::Matrix<Base::nen_, 1>& echist,
    const LINALG::Matrix<Base::nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<Base::nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<Base::nen_, 1>& epressn_timederiv,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& eaccam,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& edispnp,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& edispn,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& egridv,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& egridvn,
    const LINALG::Matrix<Base::nen_, 1>& escaaf, const LINALG::Matrix<Base::nen_, 1>* eporositynp,
    Teuchos::RCP<MAT::Material> mat, bool isale, const DRT::UTILS::GaussIntegration& intpoints)
{
  // flag for higher order elements
  Base::is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (Base::fldpara_->IsInconsistent()) Base::is_higher_order_ele_ = false;

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  SysmatOD(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren, emhist, echist,
      epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn, egridv,
      egridvn, escaaf, eporositynp, elemat1, elevec1, mat, isale, intpoints);

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::SysmatOD(Teuchos::ParameterList& params,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& ebofoaf,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& evelaf,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& evelnp,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& eveln,
    const LINALG::Matrix<Base::nen_, 1>& epreaf, const LINALG::Matrix<Base::nen_, 1>& eprenp,
    const LINALG::Matrix<Base::nen_, 1>& epren,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& emhist,
    const LINALG::Matrix<Base::nen_, 1>& echist,
    const LINALG::Matrix<Base::nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<Base::nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<Base::nen_, 1>& epressn_timederiv,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& eaccam,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& edispnp,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& edispn,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& egridv,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& egridvn,
    const LINALG::Matrix<Base::nen_, 1>& escaaf, const LINALG::Matrix<Base::nen_, 1>* eporositynp,
    LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, (Base::nsd_ + 1) * Base::nen_>& ecoupl,
    LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, 1>& eforce,
    Teuchos::RCP<const MAT::Material> material, bool isale,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  static LINALG::Matrix<Base::nen_ * Base::nsd_, Base::nen_ * Base::nsd_> ecoupl_u(
      true);  // coupling matrix for momentum equation
  static LINALG::Matrix<Base::nen_, Base::nen_ * Base::nsd_> ecoupl_p(
      true);  // coupling matrix for continuity equation

  static LINALG::Matrix<Base::nen_ * Base::nsd_, Base::nen_> ecouplp1_u(
      true);  // coupling matrix for momentum equation
  static LINALG::Matrix<Base::nen_, Base::nen_> ecouplp1_p(
      true);  // coupling matrix for continuity equation

  ecoupl_u.Clear();
  ecoupl_p.Clear();
  ecouplp1_u.Clear();
  ecouplp1_p.Clear();

  // material coordinates xyze0
  Base::xyze0_ = Base::xyze_;

  // add displacement when fluid nodes move in the ALE case (in poroelasticity this is always the
  // case)
  Base::xyze_ += edispnp;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // evaluate shape functions and derivatives at element center
  Base::EvalShapeFuncAndDerivsAtEleCenter();

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  GaussPointLoopP1OD(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp, epren, emhist, echist,
      epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn, egridv,
      egridvn, escaaf, eporositynp, eforce, ecoupl_u, ecoupl_p, ecouplp1_u, ecouplp1_p, material,
      intpoints);
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix
  //------------------------------------------------------------------------

  // add fluid velocity-structure displacement part to matrix
  for (int ui = 0; ui < Base::nen_; ++ui)
  {
    const int nsd_ui = Base::nsd_ * ui;
    const int nsdp1_ui = (Base::nsd_ + 1) * ui;

    for (int jdim = 0; jdim < Base::nsd_; ++jdim)
    {
      const int nsd_ui_jdim = nsd_ui + jdim;
      const int nsdp1_ui_jdim = nsdp1_ui + jdim;

      for (int vi = 0; vi < Base::nen_; ++vi)
      {
        const int numdof_vi = Base::numdofpernode_ * vi;
        const int nsd_vi = Base::nsd_ * vi;

        for (int idim = 0; idim < Base::nsd_; ++idim)
        {
          ecoupl(numdof_vi + idim, nsdp1_ui_jdim) += ecoupl_u(nsd_vi + idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add fluid pressure-structure displacement part to matrix
  for (int ui = 0; ui < Base::nen_; ++ui)
  {
    const int nsd_ui = Base::nsd_ * ui;
    const int nsdp1_ui = (Base::nsd_ + 1) * ui;

    for (int jdim = 0; jdim < Base::nsd_; ++jdim)
    {
      const int nsd_ui_jdim = nsd_ui + jdim;
      const int nsdp1_ui_jdim = nsdp1_ui + jdim;

      for (int vi = 0; vi < Base::nen_; ++vi)
      {
        ecoupl(Base::numdofpernode_ * vi + Base::nsd_, nsdp1_ui_jdim) += ecoupl_p(vi, nsd_ui_jdim);
      }
    }
  }

  // add fluid velocity-structure porosity part to matrix
  for (int ui = 0; ui < Base::nen_; ++ui)
  {
    const int nsdp1_ui = (Base::nsd_ + 1) * ui;

    for (int idim = 0; idim < Base::nsd_; ++idim)
    {
      const int nsdp1_ui_nsd = nsdp1_ui + Base::nsd_;

      for (int vi = 0; vi < Base::nen_; ++vi)
      {
        const int numdof_vi = Base::numdofpernode_ * vi;
        const int nsd_vi = Base::nsd_ * vi;

        ecoupl(numdof_vi + idim, nsdp1_ui_nsd) += ecouplp1_u(nsd_vi + idim, ui);
      }
    }
  }

  // add fluid pressure-structure porosity part to matrix
  for (int ui = 0; ui < Base::nen_; ++ui)
  {
    const int nsdp1_ui_nsd = (Base::nsd_ + 1) * ui + Base::nsd_;

    for (int vi = 0; vi < Base::nen_; ++vi)
      ecoupl(Base::numdofpernode_ * vi + Base::nsd_, nsdp1_ui_nsd) += ecouplp1_p(vi, ui);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::GaussPointLoopP1OD(Teuchos::ParameterList& params,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& ebofoaf,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& evelaf,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& evelnp,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& eveln,
    const LINALG::Matrix<Base::nen_, 1>& epreaf, const LINALG::Matrix<Base::nen_, 1>& eprenp,
    const LINALG::Matrix<Base::nen_, 1>& epren,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& emhist,
    const LINALG::Matrix<Base::nen_, 1>& echist,
    const LINALG::Matrix<Base::nen_, 1>& epressnp_timederiv,
    const LINALG::Matrix<Base::nen_, 1>& epressam_timederiv,
    const LINALG::Matrix<Base::nen_, 1>& epressn_timederiv,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& eaccam,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& edispnp,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& edispn,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& egridv,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& egridvn,
    const LINALG::Matrix<Base::nen_, 1>& escaaf, const LINALG::Matrix<Base::nen_, 1>* eporositynp,
    LINALG::Matrix<(Base::nsd_ + 1) * Base::nen_, 1>& eforce,
    LINALG::Matrix<Base::nen_ * Base::nsd_, Base::nen_ * Base::nsd_>& ecoupl_u,
    LINALG::Matrix<Base::nen_, Base::nen_ * Base::nsd_>& ecoupl_p,
    LINALG::Matrix<Base::nen_ * Base::nsd_, Base::nen_>& ecouplp1_u,
    LINALG::Matrix<Base::nen_, Base::nen_>& ecouplp1_p, Teuchos::RCP<const MAT::Material> material,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  // definition of velocity-based momentum residual vectors
  static LINALG::Matrix<Base::nsd_, Base::nen_ * Base::nsd_> lin_resM_Dus(true);
  static LINALG::Matrix<Base::nsd_, Base::nen_ * Base::nsd_> lin_resM_Dus_gridvel(true);
  static LINALG::Matrix<Base::nsd_, Base::nen_> lin_resM_Dphi(true);

  // set element area or volume
  const double vol = Base::fac_;

  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    lin_resM_Dus.Clear();
    lin_resM_Dus_gridvel.Clear();
    lin_resM_Dphi.Clear();

    // evaluate shape functions and derivatives at integration point
    Base::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    // evaluate shape function derivatives w.r.t. to material coordinates at integration point
    Base::SetupMaterialDerivatives();

    // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ *
    // N_XYZ_^T
    static LINALG::Matrix<Base::nsd_, Base::nsd_> defgrd(false);
    Base::ComputeDefGradient(defgrd, Base::N_XYZ_, Base::xyze_);

    // inverse deformation gradient F^-1
    static LINALG::Matrix<Base::nsd_, Base::nsd_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J and the volume change
    Base::ComputeJacobianDeterminantVolumeChange(
        Base::J_, volchange, defgrd, Base::N_XYZ_, edispnp);

    Base::EvaluateVariablesAtGaussPointOD(params, ebofoaf, evelaf, evelnp, eveln, epreaf, eprenp,
        epren, epressnp_timederiv, epressam_timederiv, epressn_timederiv, eaccam, edispnp, edispn,
        egridv, egridvn, escaaf, emhist, echist, eporositynp);

    //************************************************auxilary variables for computing the porosity

    double dphi_dp = 0.0;
    double dphi_dJ = 0.0;
    double dphi_dJdp = 0.0;
    double dphi_dJJ = 0.0;
    Base::porosity_ = 0.0;

    // compute scalar at n+alpha_F or n+1
    const double scalaraf = Base::funct_.Dot(escaaf);
    params.set<double>("scalar", scalaraf);
    ComputePorosity(params, Base::press_, volchange, *(iquad), Base::funct_, eporositynp,
        Base::porosity_, &dphi_dp, &dphi_dJ, &dphi_dJdp, &dphi_dJJ,
        nullptr,  // dphi_dpp not needed
        false);

    double refporositydot = Base::struct_mat_->RefPorosityTimeDeriv();

    //---------------------------  dJ/dx = dJ/dF : dF/dx = JF^-T : dF/dx at gausspoint
    static LINALG::Matrix<Base::nsd_, 1> gradJ(false);
    // spatial porosity gradient
    static LINALG::Matrix<Base::nsd_, 1> grad_porosity(false);
    //--------------------- linearization of porosity w.r.t. structure displacements
    static LINALG::Matrix<1, Base::nsd_ * Base::nen_> dphi_dus(false);

    //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J *
    // N_x
    static LINALG::Matrix<1, Base::nsd_ * Base::nen_> dJ_dus(false);
    //------------------ d( grad(\phi) ) / du_s = d\phi/(dJ du_s) * dJ/dx+ d\phi/dJ * dJ/(dx*du_s) +
    // d\phi/(dp*du_s) * dp/dx
    static LINALG::Matrix<Base::nsd_, Base::nen_ * Base::nsd_> dgradphi_dus(false);

    //------------------------------------ build F^-T as vector 9x1
    static LINALG::Matrix<Base::nsd_ * Base::nsd_, 1> defgrd_IT_vec(false);
    for (int i = 0; i < Base::nsd_; i++)
      for (int j = 0; j < Base::nsd_; j++) defgrd_IT_vec(i * Base::nsd_ + j) = defgrd_inv(j, i);

    // dF/dx
    static LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nsd_> F_x(false);

    // dF/dX
    static LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nsd_> F_X(false);

    Base::ComputeFDerivative(edispnp, defgrd_inv, F_x, F_X);

    // compute gradients if needed
    Base::ComputeGradients(Base::J_, dphi_dp, dphi_dJ, defgrd_IT_vec, F_x, Base::gradp_,
        eporositynp, gradJ, Base::grad_porosity_, Base::refgrad_porosity_);

    ComputeLinearizationOD(dphi_dJ, dphi_dJJ, dphi_dJdp, defgrd_inv, defgrd_IT_vec, F_x, F_X, gradJ,
        dJ_dus, dphi_dus, dgradphi_dus);

    //----------------------------------------------------------------------
    // potential evaluation of material parameters and/or stabilization
    // parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    Base::GetMaterialParamters(material);

    // set viscous term from previous iteration to zero (required for
    // using routine for evaluation of momentum rhs/residual as given)
    Base::visc_old_.Clear();
    Base::viscs2_.Clear();
    // compute viscous term from previous iteration and viscous operator
    if (Base::is_higher_order_ele_) Base::CalcDivEps(evelaf);

    Base::ComputeSpatialReactionTerms(material, defgrd_inv);

    // compute linearization of spatial reaction tensor w.r.t. structural displacements
    Base::ComputeLinSpatialReactionTerms(material, defgrd_inv, &dJ_dus, nullptr);

    // get stabilization parameters at integration point
    Base::ComputeStabilizationParameters(vol);

    // compute old RHS of momentum equation and subgrid scale velocity
    Base::ComputeOldRHSAndSubgridScaleVelocity();

    // compute old RHS of continuity equation
    Base::ComputeOldRHSConti(dphi_dp);

    // compute strong residual of mixture (structural) equation
    if (Base::porofldpara_->StabBiot() and (not Base::porofldpara_->IsStationaryConti()) and
        Base::struct_mat_->PoroLawType() != INPAR::MAT::m_poro_law_constant)
      Base::ComputeMixtureStrongResidual(params, defgrd, edispnp, edispn, F_X, *iquad, true);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = Base::fldparatimint_->TimeFac() * Base::fac_;
    const double timefacfacpre = Base::fldparatimint_->TimeFacPre() * Base::fac_;

    //***********************************************************************************************
    // 1) coupling terms in momentum balance

    Base::FillMatrixMomentumOD(timefacfac, evelaf, egridv, epreaf, dgradphi_dus, dphi_dp, dphi_dJ,
        dphi_dus, refporositydot, lin_resM_Dus, lin_resM_Dus_gridvel, ecoupl_u);

    //*************************************************************************************************************
    // 2) coupling terms in continuity equation

    Base::FillMatrixContiOD(timefacfacpre, dphi_dp, dphi_dJ, dphi_dJJ, dphi_dJdp, refporositydot,
        dgradphi_dus, dphi_dus, dJ_dus, egridv, lin_resM_Dus, lin_resM_Dus_gridvel, ecoupl_p);

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
      const double porosity_inv = 1.0 / Base::porosity_;
      for (int ui = 0; ui < Base::nen_; ++ui)
      {
        for (int idim = 0; idim < Base::nsd_; ++idim)
        {
          lin_resM_Dphi(idim, ui) +=
              timefacfac * porosity_inv * Base::reac_tensor_vel_(idim) * Base::funct_(ui);
        }
      }
    }


    if (not Base::porofldpara_->IsStationaryMomentum())
    // transient terms
    /*  reaction  */
    /*
      /                           \
     |                             |
  -  |    sigma * v_s D(phi), v    |
     |                             |
      \                           /
     */
    {
      const double porosity_inv = 1.0 / Base::porosity_;
      for (int ui = 0; ui < Base::nen_; ++ui)
      {
        for (int idim = 0; idim < Base::nsd_; ++idim)
        {
          lin_resM_Dphi(idim, ui) +=
              timefacfac * porosity_inv * (-Base::reac_tensor_gridvel_(idim)) * Base::funct_(ui);
        }
      }
    }

    // viscous terms (brinkman terms)
    if (Base::visceff_)
    {
      static LINALG::Matrix<Base::nsd_, Base::nsd_> viscstress(false);
      for (int jdim = 0; jdim < Base::nsd_; ++jdim)
      {
        for (int idim = 0; idim < Base::nsd_; ++idim)
        {
          viscstress(idim, jdim) =
              Base::visceff_ * (Base::vderxy_(jdim, idim) + Base::vderxy_(idim, jdim));
        }
      }

      static LINALG::Matrix<Base::nsd_, 1> viscstress_gradphi(false);
      viscstress_gradphi.Multiply(viscstress, Base::grad_porosity_);

      static LINALG::Matrix<Base::nsd_, Base::nen_> viscstress_derxy(false);
      viscstress_derxy.Multiply(viscstress, Base::derxy_);

      const double porosity_inv = 1.0 / Base::porosity_;

      for (int ui = 0; ui < Base::nen_; ++ui)
      {
        for (int idim = 0; idim < Base::nsd_; ++idim)
        {
          lin_resM_Dphi(idim, ui) += timefacfac * porosity_inv *
                                     (porosity_inv * viscstress_gradphi(idim) * Base::funct_(ui) -
                                         viscstress_derxy(idim, ui));
        }
      }
    }

    for (int ui = 0; ui < Base::nen_; ++ui)
    {
      for (int vi = 0; vi < Base::nen_; ++vi)
      {
        const int fvi = Base::nsd_ * vi;
        for (int idim = 0; idim < Base::nsd_; ++idim)
        {
          ecouplp1_u(fvi + idim, ui) += Base::funct_(vi) * lin_resM_Dphi(idim, ui);
        }
      }
    }

    //*************************************************************************************************************
    // 3.2) Continuity equation

    // transient terms
    /* time derivative*/
    /*
      /                           \
     |                             |
  -  |    D(phi), v                |
     |                             |
      \                           /
     */
    {
      for (int ui = 0; ui < Base::nen_; ++ui)
        for (int vi = 0; vi < Base::nen_; ++vi)
          ecouplp1_p(vi, ui) += Base::fac_ * Base::funct_(vi) * Base::funct_(ui);
    }

    static LINALG::Matrix<Base::nen_, 1> derxy_convel(false);
    derxy_convel.Clear();

    for (int i = 0; i < Base::nen_; i++)
      for (int j = 0; j < Base::nsd_; j++) derxy_convel(i) += Base::derxy_(j, i) * Base::velint_(j);

    if (not Base::porofldpara_->IsStationaryConti())
    {
      for (int i = 0; i < Base::nen_; i++)
        for (int j = 0; j < Base::nsd_; j++)
          derxy_convel(i) += Base::derxy_(j, i) * (-Base::gridvel_int_(j));
    }

    if (!static_cast<DRT::ELEMENTS::FluidEleParameterPoro*>(Base::fldpara_)->PoroContiPartInt())
    {
      /*
        /                           \     /                             \
       |                             |    |                              |
       |    \nabla v_f D(phi), v     | +  |  (v_f-v_s) \nabla  D(phi), v |
       |                             |    |                              |
        \                           /     \                             /
       */
      for (int ui = 0; ui < Base::nen_; ++ui)
      {
        for (int vi = 0; vi < Base::nen_; ++vi)
        {
          ecouplp1_p(vi, ui) += +timefacfacpre * Base::vdiv_ * Base::funct_(vi) * Base::funct_(ui) +
                                timefacfacpre * Base::funct_(vi) * derxy_convel(ui);
        }
      }
    }
    else
    {
      /*
          /                             \
          |                              |
       -  |  (v_f-v_s) \nabla  D(phi), v |
          |                              |
          \                             /
       */
      for (int ui = 0; ui < Base::nen_; ++ui)
      {
        for (int vi = 0; vi < Base::nen_; ++vi)
        {
          ecouplp1_p(vi, ui) += -1.0 * timefacfacpre * derxy_convel(vi) * Base::funct_(ui);
        }
      }
      /*
          /                             \
          |                              |
          |  \nabla v_s D(phi), v        |
          |                              |
          \                             /
       */
      if (not Base::porofldpara_->IsStationaryConti())
      {
        for (int ui = 0; ui < Base::nen_; ++ui)
        {
          for (int vi = 0; vi < Base::nen_; ++vi)
          {
            ecouplp1_p(vi, ui) +=
                timefacfacpre * Base::funct_(vi) * Base::gridvel_div_ * Base::funct_(ui);
          }
        }
      }
    }

    //*************************************************************************************************************
    // PSPG
    if (Base::fldpara_->PSPG())
    {
      double scal_grad_q = 0.0;

      if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
      {
        scal_grad_q = Base::tau_(1);
      }
      else
      {
        scal_grad_q = 0.0;  // Base::fldpara_->AlphaF()*fac3;
      }

      {
        const double v1 = -timefacfacpre * Base::dtau_dphi_(1) / scal_grad_q;
        for (int ui = 0; ui < Base::nen_; ++ui)
        {
          for (int idim = 0; idim < Base::nsd_; ++idim)
          {
            const double v = v1 * Base::sgvelint_(idim) * Base::funct_(ui);

            for (int vi = 0; vi < Base::nen_; ++vi)
            {
              ecouplp1_p(vi, ui) += v * Base::derxy_(idim, vi);
            }
          }
        }
      }

      // linearization of residual in stabilization term w.r.t. porosity
      if (Base::is_higher_order_ele_ || Base::fldpara_->IsNewton())
      {
        static LINALG::Matrix<Base::nen_, Base::nen_> temp(false);
        temp.Clear();

        for (int vi = 0; vi < Base::nen_; ++vi)
        {
          for (int ui = 0; ui < Base::nen_; ++ui)
            for (int idim = 0; idim < Base::nsd_; ++idim)
              temp(vi, ui) += Base::derxy_(idim, vi) * lin_resM_Dphi(idim, ui);
        }

        for (int ui = 0; ui < Base::nen_; ++ui)
        {
          for (int vi = 0; vi < Base::nen_; ++vi)
          {
            ecouplp1_p(vi, ui) += scal_grad_q * temp(vi, ui);
          }
        }
      }
    }

    //*************************************************************************************************************
    if (Base::fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
    {
      double reac_tau;
      if (Base::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
        reac_tau = Base::fldpara_->ViscReaStabFac() * Base::reacoeff_ * Base::tau_(1);
      else
      {
        dserror("Is this factor correct? Check for bugs!");
        reac_tau = 0.0;
      }

      if (Base::is_higher_order_ele_ or Base::fldpara_->IsNewton())
      {
        for (int vi = 0; vi < Base::nen_; ++vi)
        {
          const double v = reac_tau * Base::funct_(vi);

          for (int idim = 0; idim < Base::nsd_; ++idim)
          {
            const int fvi_p_idim = Base::nsd_ * vi + idim;

            for (int ui = 0; ui < Base::nen_; ++ui)
            {
              ecouplp1_u(fvi_p_idim, ui) += v * lin_resM_Dphi(idim, ui);
            }
          }
        }
      }

      {  // linearization of stabilization parameter w.r.t. porosity
        const double v = timefacfac * Base::fldpara_->ViscReaStabFac() *
                         (Base::reacoeff_ * Base::dtau_dphi_(1) / Base::tau_(1) +
                             Base::reacoeff_ / Base::porosity_);
        for (int vi = 0; vi < Base::nen_; ++vi)
        {
          const double w = -1.0 * v * Base::funct_(vi);

          for (int idim = 0; idim < Base::nsd_; ++idim)
          {
            const double w_sgvelint = w * Base::sgvelint_(idim);
            const int fvi = Base::nsd_ * vi + idim;

            for (int ui = 0; ui < Base::nen_; ++ui)
            {
              ecouplp1_u(fvi, ui) += w_sgvelint * Base::funct_(ui);
            }
          }
        }
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputeLinearization(const double& dphi_dp,
    const double& dphi_dpp, const double& dphi_dJdp, const LINALG::Matrix<Base::nsd_, 1>& gradJ,
    LINALG::Matrix<Base::nsd_, Base::nen_>& dgradphi_dp)
{
  // porosity is a primary variable -> d(grad(phi)/d(pressure)) is zero!
  dgradphi_dp.Clear();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputeLinearizationOD(const double& dphi_dJ,
    const double& dphi_dJJ, const double& dphi_dJp,
    const LINALG::Matrix<Base::nsd_, Base::nsd_>& defgrd_inv,
    const LINALG::Matrix<Base::nsd_ * Base::nsd_, 1>& defgrd_IT_vec,
    const LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nsd_>& F_x,
    const LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nsd_>& F_X,
    const LINALG::Matrix<Base::nsd_, 1>& gradJ, LINALG::Matrix<1, Base::nsd_ * Base::nen_>& dJ_dus,
    LINALG::Matrix<1, Base::nsd_ * Base::nen_>& dphi_dus,
    LINALG::Matrix<Base::nsd_, Base::nen_ * Base::nsd_>& dgradphi_dus)
{
  //------------------------------------------------dJ/dus = dJ/dF : dF/dus = J * F^-T . N_X = J *
  // N_x
  for (int i = 0; i < Base::nen_; i++)
    for (int j = 0; j < Base::nsd_; j++) dJ_dus(j + i * Base::nsd_) = Base::J_ * Base::derxy_(j, i);

  //--------------------- linearization of porosity w.r.t. structure displacements
  // porosity is a primary variable -> d(grad(phi)/d(displacement)) is zero!
  dphi_dus.Clear();

  //------------------ d( grad(\phi) ) / du_s is also zero
  dgradphi_dus.Clear();

  if (!static_cast<DRT::ELEMENTS::FluidEleParameterPoro*>(Base::fldpara_)->PoroContiPartInt() or
      Base::visceff_)
  {
    //---------------------d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J *
    // F^-T : N_X_x

    // dF^-T/dus : dF/dx = - (F^-1. dN/dx . u_s)^T  : dF/dx
    static LINALG::Matrix<Base::nsd_, Base::nsd_ * Base::nen_> dFinvdus_dFdx(false);
    dFinvdus_dFdx.Clear();
    for (int i = 0; i < Base::nsd_; i++)
    {
      for (int n = 0; n < Base::nen_; n++)
      {
        for (int j = 0; j < Base::nsd_; j++)
        {
          const int gid = Base::nsd_ * n + j;
          const double defgrd_inv_ij = defgrd_inv(i, j);
          for (int k = 0; k < Base::nsd_; k++)
          {
            const double derxy_kn = Base::derxy_(k, n);
            for (int p = 0; p < Base::nsd_; p++)
              dFinvdus_dFdx(p, gid) += -defgrd_inv_ij * derxy_kn * F_x(k * Base::nsd_ + i, p);
          }
        }
      }
    }

    // F^-T : d(dF/dx)/dus =  F^-T : (N,XX * F^ -1 + dF/dX * F^-1 * N,x)
    static LINALG::Matrix<Base::nsd_, Base::nsd_ * Base::nen_> FinvT_dFx_dus(false);
    FinvT_dFx_dus.Clear();

    for (int n = 0; n < Base::nen_; n++)
    {
      for (int j = 0; j < Base::nsd_; j++)
      {
        const int gid = Base::nsd_ * n + j;
        for (int p = 0; p < Base::nsd_; p++)
        {
          double val = 0.0;
          const double derxy_p_n = Base::derxy_(p, n);
          for (int k = 0; k < Base::nsd_; k++)
          {
            const double defgrd_inv_kj = defgrd_inv(k, j);
            const double defgrd_inv_kp = defgrd_inv(k, p);
            for (int i = 0; i < Base::nsd_; i++)
            {
              val += defgrd_inv(i, j) * Base::N_XYZ2full_(i * Base::nsd_ + k, n) * defgrd_inv_kp;
              for (int l = 0; l < Base::nsd_; l++)
                val += -defgrd_inv(i, l) * F_X(i * Base::nsd_ + l, k) * defgrd_inv_kj * derxy_p_n;
            }
          }
          FinvT_dFx_dus(p, gid) += val;
        }
      }
    }

    //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
    static LINALG::Matrix<1, Base::nsd_> temp;
    temp.MultiplyTN(defgrd_IT_vec, F_x);

    //----d(gradJ)/dus =  dJ/dus * F^-T . : dF/dx + J * dF^-T/dus : dF/dx + J * F^-T : N_X_x
    static LINALG::Matrix<Base::nsd_, Base::nen_ * Base::nsd_> dgradJ_dus;

    dgradJ_dus.MultiplyTN(temp, dJ_dus);

    dgradJ_dus.Update(Base::J_, dFinvdus_dFdx, 1.0);

    dgradJ_dus.Update(Base::J_, FinvT_dFx_dus, 1.0);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::PSPG(
    LINALG::Matrix<Base::nen_, Base::nen_ * Base::nsd_>& estif_q_u,
    LINALG::Matrix<Base::nen_, Base::nen_>& ppmat, LINALG::Matrix<Base::nen_, 1>& preforce,
    const LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nen_>& lin_resM_Du,
    const LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nen_>& lin_resMRea_Du,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& lin_resM_Dp, const double& dphi_dp,
    const double& fac3, const double& timefacfac, const double& timefacfacpre, const double& rhsfac)
{
  Base::PSPG(estif_q_u, ppmat, preforce, lin_resM_Du, lin_resMRea_Du, lin_resM_Dp, dphi_dp, fac3,
      timefacfac, timefacfacpre, rhsfac);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ReacStab(
    LINALG::Matrix<Base::nen_ * Base::nsd_, Base::nen_ * Base::nsd_>& estif_u,
    LINALG::Matrix<Base::nen_ * Base::nsd_, Base::nen_>& estif_p_v,
    LINALG::Matrix<Base::nsd_, Base::nen_>& velforce,
    LINALG::Matrix<Base::nsd_ * Base::nsd_, Base::nen_>& lin_resM_Du,
    const LINALG::Matrix<Base::nsd_, Base::nen_>& lin_resM_Dp, const double& dphi_dp,
    const double& timefacfac, const double& timefacfacpre, const double& rhsfac, const double& fac3)
{
  Base::ReacStab(estif_u, estif_p_v, velforce, lin_resM_Du, lin_resM_Dp, dphi_dp, timefacfac,
      timefacfacpre, rhsfac, fac3);
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::ComputeVolume(Teuchos::ParameterList& params,
    DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1)
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype, Base::nsd_, LINALG::Matrix<Base::nsd_, Base::nen_>>(
      ele, Base::xyze_);
  // set element id
  Base::eid_ = ele->Id();

  LINALG::Matrix<Base::nsd_, Base::nen_> edispnp(true);
  LINALG::Matrix<Base::nen_, 1> eporositynp(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &edispnp, &eporositynp, "dispnp");

  LINALG::Matrix<Base::nsd_, Base::nen_> egridvnp(true);
  Base::ExtractValuesFromGlobalVector(
      discretization, lm, *Base::rotsymmpbc_, &egridvnp, nullptr, "gridv");

  // get new node positions of ALE mesh
  Base::xyze_ += edispnp;

  // integration loop
  for (DRT::UTILS::GaussIntegration::iterator iquad = Base::intpoints_.begin();
       iquad != Base::intpoints_.end(); ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    Base::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    //-----------------------------------computing the porosity
    Base::porosity_ = Base::funct_.Dot(eporositynp);

    // get structure velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    LINALG::Matrix<Base::nsd_, Base::nsd_> gridvelderxy;
    gridvelderxy.MultiplyNT(egridvnp, Base::derxy_);

    Base::gridvel_div_ = 0.0;
    for (int idim = 0; idim < Base::nsd_; ++idim) Base::gridvel_div_ += gridvelderxy(idim, idim);

    elevec1(0) += Base::porosity_ * Base::fac_;
  }

  return 0;
}

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
