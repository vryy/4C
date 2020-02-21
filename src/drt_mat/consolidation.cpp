/*----------------------------------------------------------------------*/
/*! \file
\brief managing common consolidation methods

\level 3
\maintainer Sebastian Proell
 */

/*----------------------------------------------------------------------*
 |  headers                                                                |
 *----------------------------------------------------------------------*/
#include "consolidation.H"
#include "matpar_bundle.H"
#include "../drt_tsi/tsi_defines.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function_library.H"

#include <algorithm>
#include <iterator>

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::Consolidation::Consolidation(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      solidustemp_(matdata->GetDouble("SOLIDUS")),
      liquidustemp_(matdata->GetDouble("LIQUIDUS")),
      latentheat_(matdata->GetDouble("LATENTHEAT"))
{
  const Teuchos::ParameterList& tdynparams = DRT::Problem::Instance()->ThermalDynamicParams();
  heatint_ = DRT::INPUT::IntegralValue<int>(tdynparams, "HEATINTEGRATION") == 1;
  // TODO create different evaluator strategies here
  // only use enthalpy dependent parameters if both isothermal AND heat integration

  if (abs(solidustemp_ - liquidustemp_) < 1e-12)
  {
    enthalpydep_ = true;
    if (!heatint_)
      dserror("Isothermal phase change requires HEATINTEGRATION=yes in thermal dynamic.");
    if (latentheat_ <= 0) dserror("Isothermal phase change requires LATENTHEAT>0 in material.");
  }
  else
    enthalpydep_ = false;
}

Teuchos::RCP<MAT::Material> MAT::PAR::Consolidation::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Consolidation(this));
}

MAT::ConsolidationType MAT::ConsolidationType::instance_;

DRT::ParObject* MAT::ConsolidationType::Create(const std::vector<char>& data)
{
  MAT::Consolidation* mat = new MAT::Consolidation();
  mat->Unpack(data);
  return mat;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public) proell 05/18 |
 *----------------------------------------------------------------------*/
MAT::Consolidation::Consolidation()
    : params_(NULL), cfracn_(Teuchos::null), cfracnp_(Teuchos::null), cfracnode_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)  proell 05/18 |
 *----------------------------------------------------------------------*/
MAT::Consolidation::Consolidation(MAT::PAR::Consolidation* params)
    : params_(params), cfracn_(Teuchos::null), cfracnp_(Teuchos::null), cfracnode_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public) proell 05/18 |
 *----------------------------------------------------------------------*/
void MAT::Consolidation::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack history data
  int histsize;
  // if material is not initialised, i.e. start simulation, nothing to pack
  if (cfracn_ == Teuchos::null)
    histsize = 0;
  else
    histsize = cfracn_->size();

  // TODO this only works for same number of GPs and nodes

  AddtoPack(data, histsize);  // length of history vector(s)
  // add each component of history vectors
  for (int i = 0; i < histsize; i++)
  {
    AddtoPack(data, cfracn_->at(i));
    AddtoPack(data, cfracnp_->at(i));
    AddtoPack(data, cfracnode_->at(i));
    AddtoPack(data, nodallump_->at(i));
    AddtoPack(data, lhrem_->at(i));
  }
  // NOTE: we do not pack internal helpers, just reinitialize them in Unpack()
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public) proell 05/18 |
 *----------------------------------------------------------------------*/
void MAT::Consolidation::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Consolidation*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  // history variables
  int histsize = 0;
  ExtractfromPack(position, data, histsize);

  cfracn_ = Teuchos::rcp(new std::vector<double>);
  cfracnp_ = Teuchos::rcp(new std::vector<double>);
  cfracnode_ = Teuchos::rcp(new std::vector<double>);
  nodallump_ = new std::vector<double>;
  lhrem_ = new std::vector<double>;
  lhremnp_ = new std::vector<double>;

  double temp_d;
  for (int i = 0; i < histsize; i++)
  {
    ExtractfromPack(position, data, temp_d);
    cfracn_->push_back(temp_d);

    ExtractfromPack(position, data, temp_d);
    cfracnp_->push_back(temp_d);

    ExtractfromPack(position, data, temp_d);
    cfracnode_->push_back(temp_d);

    ExtractfromPack(position, data, temp_d);
    nodallump_->push_back(temp_d);

    ExtractfromPack(position, data, temp_d);
    lhrem_->push_back(temp_d);
    // the step updated vector can be intialized to last step values
    lhremnp_->push_back(temp_d);
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  // create the internal helper variables, not(!) history so no communication necessary
  SetupFunctionEvaluateVariables();
}

/*----------------------------------------------------------------------*
 |  setup history-variables                       (public) proell 05/18 |
 *----------------------------------------------------------------------*/
void MAT::Consolidation::Setup(const int numgp)
{
  // create and resize GP-sized history vectors
  cfracn_ = Teuchos::rcp(new std::vector<double>);
  cfracn_->resize(numgp, false);

  cfracnp_ = Teuchos::rcp(new std::vector<double>);
  cfracnp_->resize(numgp, -1.0);

  // TODO implement for general number of nodes and GPs
  // assume there are as man nodes as GPs
  cfracnode_ = Teuchos::rcp(new std::vector<double>);
  cfracnode_->resize(numgp, -1.0);

  nodallump_ = new std::vector<double>(numgp);
  lhrem_ = new std::vector<double>(numgp);
  lhremnp_ = new std::vector<double>(numgp);

  SetupFunctionEvaluateVariables();
}


void MAT::Consolidation::SetupFunctionEvaluateVariables()
{
  frac_ = new std::vector<double>(NUMFUNCT);
  frac_T_ = new std::vector<double>(NUMFUNCT);
  enth_ = 0;

  if (params_->enthalpydep_)
  {
    Ts_mod_ = 0;
    Tl_mod_ = params_->latentheat_;
    fractionInterpolant = Teuchos::rcp(new LinearInterpolant(Ts_mod_, Tl_mod_, 0, 1));
  }
  else
  {
    Ts_mod_ = params_->solidustemp_;
    Tl_mod_ = params_->liquidustemp_;
    fractionInterpolant = Teuchos::rcp(new LinearInterpolant(Ts_mod_, Tl_mod_, 0, 1));
  }

  // create apparent capacity interpolants
  // TODO unnecessary to have these in each material instance, they are always the same
  const double Tmid = (Ts_mod_ + Tl_mod_) / 2;
  const double Cpeak = 2 * params_->latentheat_ / (Tl_mod_ - Ts_mod_);
  const double delta = (Tl_mod_ - Ts_mod_) / 6;
  leftInterpAppCapa_ =
      Teuchos::rcp(new QuadLinQuadInterpolant(Ts_mod_, Tmid, 0, Cpeak, 0, 0, delta));
  rightInterpAppCapa_ =
      Teuchos::rcp(new QuadLinQuadInterpolant(Tmid, Tl_mod_, Cpeak, 0, 0, 0, delta));
}


/*----------------------------------------------------------------------*
 |  update history variables                     (public) proell 05/18  |
 *----------------------------------------------------------------------*/
void MAT::Consolidation::Update(double temperature, int gp)
{
  // TODO remove when heat integration is moved out of consolidation code
  if (params_->enthalpydep_)
    cfracn_->at(gp) = CFracNext(enth_, cfracn_->at(gp));
  else
    cfracn_->at(gp) = CFracNext(temperature, cfracn_->at(gp));
  // careful: copy the actual content and not just the pointer!
  *lhrem_ = *lhremnp_;
}

/*----------------------------------------------------------------------*
 |  function evaluate routine                    (public) proell 05/18  |
 *----------------------------------------------------------------------*/
double MAT::Consolidation::EvaluateTempDependentFunction(
    const double temperature,  ///< temperature at GP
    const unsigned int gp,     ///< current GP to select internal phase information
    const std::vector<int>
        functions  ///< vector of global functions IDs to interpolate, length NUMFUNCT
)
{
  EvaluateCFracnpAtGp(temperature, gp);

  if (temperature > Tl_mod_)
  {
    return FunctionValue(temperature, functions[1]);
  }
  if (temperature < Ts_mod_)
  {
    const double powderVal = FunctionValue(temperature, functions[0]);
    const double solidVal = FunctionValue(temperature, functions[2]);
    return cfracnp_->at(gp) * solidVal + (1 - cfracnp_->at(gp)) * powderVal;
  }
  else
  {
    // quasi-linear interpolation
    EvaluateCFracnpAtGp(temperature, gp);  // also evaluates Tfrac
    EvaluateFractions(temperature, cfracnp_->at(gp));
    const double powderVal = FunctionValue(temperature, functions[0]);
    const double meltVal = FunctionValue(temperature, functions[1]);
    const double solidVal = FunctionValue(temperature, functions[2]);
    return powderVal * frac_->at(0) + meltVal * frac_->at(1) + solidVal * frac_->at(2);
  }
}


/*----------------------------------------------------------------------*
 |  derivative evaluate routine                  (public) proell 05/18  |
 *----------------------------------------------------------------------*/
double MAT::Consolidation::EvaluateTempDependentDerivative(
    const double temperature, const unsigned int gp, const std::vector<int> functions)
{
  EvaluateCFracnpAtGp(temperature, gp);

  if (temperature > Tl_mod_)
  {
    return FunctionDerivative(temperature, functions[1]);
  }
  if (temperature < Ts_mod_)
  {
    const double powderDeriv = FunctionDerivative(temperature, functions[0]);
    const double solidDeriv = FunctionDerivative(temperature, functions[2]);
    return cfracnp_->at(gp) * solidDeriv + (1 - cfracnp_->at(gp)) * powderDeriv;
  }
  else
  {
    // evaluate according to the product rule
    EvaluateCFracnpAtGp(temperature, gp);
    EvaluateFractions(temperature, cfracnp_->at(gp));
    EvaluateFractionDerivatives(temperature, cfracnp_->at(gp), cfracn_->at(gp));
    const double powderDeriv = FunctionDerivative(temperature, functions[0]);
    const double meltDeriv = FunctionDerivative(temperature, functions[1]);
    const double solidDeriv = FunctionDerivative(temperature, functions[2]);
    const double powderVal = FunctionValue(temperature, functions[0]);
    const double meltVal = FunctionValue(temperature, functions[1]);
    const double solidVal = FunctionValue(temperature, functions[2]);
    const double deriv = powderDeriv * frac_->at(0) + meltDeriv * frac_->at(1) +
                         solidDeriv * frac_->at(2) + powderVal * frac_T_->at(0) +
                         meltVal * frac_T_->at(1) + solidVal * frac_T_->at(2);
    return deriv;
  }
}


/*----------------------------------------------------------------------*
 |  helper to determine fracs                   (private) proell 05/18  |
 *----------------------------------------------------------------------*/
void MAT::Consolidation::EvaluateFractions(const double temperature, const double cfracnp)
{
  // all fraction are functions of consolidated fraction and (fractional) temperature
  EvaluateTFrac(temperature);
  // powder fraction = 1-r_c (everything that is not consolidated)
  frac_->at(0) = 1 - cfracnp;
  // melt fraction is independent of history
  frac_->at(1) = Tfrac_;
  // solid fraction = r_c-r_m
  frac_->at(2) = cfracnp - Tfrac_;
}

/*----------------------------------------------------------------------*
 |  helper to determine frac derivatives        (private) proell 05/18  |
 *----------------------------------------------------------------------*/
void MAT::Consolidation::EvaluateFractionDerivatives(
    const double temperature, const double cfracnp, const double cfracn)
{
  const double slope = 1 / (Tl_mod_ - Ts_mod_);
  if (temperature > Ts_mod_ && temperature < Tl_mod_)
  {
    // melt fraction always changes with calculated slope
    (*frac_T_)[1] = slope;
    // depending on an increase/decrease in temperature eitehr solid or powder changes
    // here we always split the slope onto both fractions as a robust approximation
    (*frac_T_)[0] = -0.5 * slope;
    (*frac_T_)[2] = -0.5 * slope;
  }
  else
    std::fill(frac_T_->begin(), frac_T_->end(), 0);
}

/*----------------------------------------------------------------------*
 |  definition of the cfrac update formula      (private) proell 05/18  |
 *----------------------------------------------------------------------*/
double MAT::Consolidation::CFracNext(const double temperature, const double oldcfrac)
{
  EvaluateTFrac(temperature);
  // cfrac_{n+1} = max(cfrac_{n}; T*)
  return std::max(oldcfrac, Tfrac_);
}

/*-----------------------------------------------------------------------*
 |  smoothed fractional temperature              (private) proell 07/18  |
 *-----------------------------------------------------------------------*/
void MAT::Consolidation::EvaluateTFrac(const double temperature)
{
  Tfrac_ = fractionInterpolant->valueAt(temperature);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MAT::Consolidation::EvaluateCFracnpAtGp(const double temperature, const int gp)
{
  cfracnp_->at(gp) = CFracNext(temperature, cfracn_->at(gp));
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
double MAT::Consolidation::FunctionValue(const double temperature, const int function)
{
  return (dynamic_cast<DRT::UTILS::FastPolynomialFunction&>(
              DRT::Problem::Instance()->Funct(function - 1)))
      .Evaluate(temperature);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
double MAT::Consolidation::FunctionDerivative(const double temperature, const int function)
{
  return (dynamic_cast<DRT::UTILS::FastPolynomialFunction&>(
              DRT::Problem::Instance()->Funct(function - 1)))
      .EvaluateDerivative(temperature);
}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
double MAT::Consolidation::ScalarProduct(std::vector<double>* vec1, std::vector<double>* vec2)
{
  return std::inner_product(vec1->begin(), vec1->end(), vec2->begin(), 0.0);
}
