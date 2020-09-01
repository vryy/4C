/*----------------------------------------------------------------------*/
/*! \file
\brief abstract interface for electrode and electrolyte materials carrying concentration and
electric potential as degrees of freedom


\level 2
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include "elchsinglemat.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::ElchSingleMat::ElchSingleMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      diffusion_coefficient_concentration_dependence_funct_(
          matdata->GetInt("DIFF_COEF_CONC_DEP_FUNCT")),
      number_diffusion_coefficent_params_(matdata->GetInt("DIFF_PARA_NUM")),
      diffusion_coefficent_params_(*matdata->Get<std::vector<double>>("DIFF_PARA")),
      conductivity_concentration_dependence_funct_(matdata->GetInt("COND_CONC_DEP_FUNCT")),
      number_conductivity_params_(matdata->GetInt("COND_PARA_NUM")),
      conductivity_params_(*matdata->Get<std::vector<double>>("COND_PARA"))
{
  // safety checks
  if (number_diffusion_coefficent_params_ != static_cast<int>(diffusion_coefficent_params_.size()))
    dserror("Mismatch in number of parameters for diffusion coefficient!");
  if (number_conductivity_params_ != static_cast<int>(conductivity_params_.size()))
    dserror("Mismatch in number of parameters for conductivity!");
  CheckProvidedParams(
      diffusion_coefficient_concentration_dependence_funct_, diffusion_coefficent_params_);
  CheckProvidedParams(conductivity_concentration_dependence_funct_, conductivity_params_);
}


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void MAT::PAR::ElchSingleMat::CheckProvidedParams(
    const int functnr, const std::vector<double>& functparams)
{
  // name of specified curve
  std::string functionname;

  // expected number of parameters for specified curve
  unsigned int nfunctparams = 0;

  // check set of implemented functions with negative curve number
  if (functnr < 0)
  {
    switch (functnr)
    {
      case -1:
      {
        // constant value: functval=functparams[0];
        functionname = "'constant value'";
        nfunctparams = 1;
        break;
      }
      case -2:
      {
        // linear function: functval=functparams[0]+functparams[1]*concentration;
        functionname = "'linear function'";
        nfunctparams = 2;
        break;
      }
      case -3:
      {
        // quadratic function:
        // functval=functparams[0]+functparams[1]*concentration+functparams[2]*concentration*concentration;
        functionname = "'quadratic function'";
        nfunctparams = 3;
        break;
      }
      case -4:
      {
        // power function: functval=functparams[0]*pow(concentration,functparams[1]);
        functionname = "'power function'";
        nfunctparams = 2;
        break;
      }
      case -5:
      {
        // function 1 for conductivity;
        functionname = "'function 1 for conductivity'";
        nfunctparams = 4;
        break;
      }
      case -6:
      {
        // a0*c + a1*c^1.5 + a2*c^3
        functionname = "'a0*c + a1*c^1.5 + a2*c^3'";
        nfunctparams = 3;
        break;
      }
      case -7:
      {
        // a0 + a1*c + a2*c^2 + a3*c^3
        functionname = "'a0 + a1*c + a2*c^2 + a3*c^3'";
        nfunctparams = 4;
        break;
      }
      case -8:
      {
        // thermodynamic factor Nyman 2008
        functionname = "'function thermodynamic factor (Nyman 2008)'";
        nfunctparams = 7;
        break;
      }
      case -9:
      {
        // linear thermodynamic factor including Debye-Hückel theory
        functionname = "'function linear thermodynamic factor (including Debye Hueckel theory)'";
        nfunctparams = 2;
        break;
      }
      case -10:
      {
        // function 1 for conductivity
        functionname = "'function 1 for conductivity: own definition'";
        nfunctparams = 6;
        break;
      }
      case -11:
      {
        // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann,
        // Kee (2012) kappa = a0*c*exp(a1*c^a2)
        functionname =
            "'conductivity as a function of concentration according to Goldin, Colclasure, "
            "Wiedemann, Kee (2012)'";
        nfunctparams = 3;
        break;
      }
      case -12:
      {
        // diffusion coefficient based on a function defined in
        // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion
        // Coefficients in LiPF6 Electrolytic Solutions Journal of The Electrochemical Society,
        // 2008, 155, F13-F16 diff = a0*exp(-a1*c^a2)
        functionname = "'diffusion coefficient as an exponential function: a1*exp(a2*c)'";
        nfunctparams = 2;
        break;
      }
      case -13:
      {
        // TDF based on a function defined in
        // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
        // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
        // sqrt(c)/(1+a2*sqrt(c)) + a2*c
        functionname =
            "'TDF as as a function of concentration according to Landesfeind, Ehrl, Graf, Wall, "
            "Gasteiger (2015)'";
        nfunctparams = 3;
        break;
      }
      default:
      {
        dserror("Curve number %i is not implemented", functnr);
        break;
      }
    }

    // safety check
    if (functparams.size() != nfunctparams)
      dserror(
          "Number of provided parameters does not match number of expected parameters for function "
          "with curve number %i (%s)!",
          functnr, functionname.c_str());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeDiffusionCoefficient(const double concentration) const
{
  double diffusionCoefficient(0.0);

  // evaluate pre implemented concentration dependent diffusion coefficient
  if (DiffusionCoefficientConcentrationDependenceFunct() < 0)
    diffusionCoefficient = EvalFunctValue(DiffusionCoefficientConcentrationDependenceFunct(),
        concentration, DiffusionCoefficientParams());

  // diffusion coefficient is a constant prescribed in input file
  else if (DiffusionCoefficientConcentrationDependenceFunct() == 0)
    diffusionCoefficient = EvalFunctValue(-1, concentration, DiffusionCoefficientParams());

  // diffusion coefficient is a function of the concentration as defined in the input file
  else
    diffusionCoefficient = DRT::Problem::Instance()
                               ->Funct(DiffusionCoefficientConcentrationDependenceFunct() - 1)
                               .EvaluateTime(concentration);

  return diffusionCoefficient;
}


/*------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeFirstDerivDiffCoeff(const double concentration) const
{
  double firstderiv(0.0);

  if (DiffusionCoefficientConcentrationDependenceFunct() < 0)
    firstderiv = EvalFirstDerivFunctValue(DiffusionCoefficientConcentrationDependenceFunct(),
        concentration, DiffusionCoefficientParams());
  else if (DiffusionCoefficientConcentrationDependenceFunct() == 0)
    firstderiv = EvalFirstDerivFunctValue(-1, concentration, DiffusionCoefficientParams());
  else
    firstderiv = (DRT::Problem::Instance()
                      ->Funct(DiffusionCoefficientConcentrationDependenceFunct() - 1)
                      .EvaluateTimeDerivative(concentration, 1))[1];

  return firstderiv;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeConductivity(const double concentration) const
{
  double conductivity(0.0);

  // evaluate pre implemented concentration dependent conductivity
  if (ConductivityConcentrationDependenceFunct() < 0)
    conductivity = EvalFunctValue(
        ConductivityConcentrationDependenceFunct(), concentration, ConductivityParams());

  // conductivitiy is a constant prescribed in input file
  else if (ConductivityConcentrationDependenceFunct() == 0)
    conductivity = EvalFunctValue(-1, concentration, ConductivityParams());

  // conductivitiy is a function of the concentration as defined in the input file
  else
    conductivity = DRT::Problem::Instance()
                       ->Funct(ConductivityConcentrationDependenceFunct() - 1)
                       .EvaluateTime(concentration);

  return conductivity;
}


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeFirstDerivCond(const double concentration) const
{
  double firstderiv(0.0);

  if (ConductivityConcentrationDependenceFunct() < 0)
    firstderiv = EvalFirstDerivFunctValue(
        ConductivityConcentrationDependenceFunct(), concentration, ConductivityParams());
  else if (ConductivityConcentrationDependenceFunct() == 0)
    firstderiv = EvalFirstDerivFunctValue(-1, concentration, ConductivityParams());
  else
    firstderiv = (DRT::Problem::Instance()
                      ->Funct(ConductivityConcentrationDependenceFunct() - 1)
                      .EvaluateTimeDerivative(concentration, 1))[1];

  return firstderiv;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::EvalFunctValue(
    const int functnr, const double concentration, const std::vector<double>& functparams) const
{
  double functval(0.0);

  switch (functnr)
  {
    // a0
    case -1:
      functval = functparams[0];
      break;

    // a0 + a1*c
    case -2:
      functval = functparams[0] + functparams[1] * concentration;
      break;

    // a0 + a1*c + a2*c^2
    case -3:
      functval = functparams[0] + functparams[1] * concentration +
                 functparams[2] * concentration * concentration;
      break;

    // a0*c^a1
    case -4:
      functval = functparams[0] * std::pow(concentration, functparams[1]);
      break;

    // conductivity
    case -5:
    {
      const double nenner =
          (1.0 + functparams[2] * concentration * concentration -
              functparams[3] * concentration * concentration * concentration * concentration);
      // functparams[0]*(functparams[1]*concentration/nenner) + 0.01 -> constant level 0.01 deleted
      // since it does not have a physical meaning (28.04.2014)
      functval = functparams[0] * (functparams[1] * concentration / nenner);
      break;
    }

    // a0*c + a1*c^1.5 + a2*c^3
    case -6:
      functval = functparams[0] * concentration + functparams[1] * std::pow(concentration, 1.5) +
                 functparams[2] * concentration * concentration * concentration;
      break;

    // a0 + a1*c + a2*c^2 + a3*c^3
    case -7:
      functval = functparams[0] + functparams[1] * concentration +
                 functparams[2] * concentration * concentration +
                 functparams[3] * concentration * concentration * concentration;
      break;

    // thermodynamic factor Nyman 2008
    case -8:
    {
      const double num = functparams[0] + functparams[1] * concentration +
                         functparams[2] * concentration * concentration;
      const double denom = functparams[3] + functparams[4] * concentration +
                           functparams[5] * concentration * concentration +
                           functparams[6] * concentration * concentration * concentration;
      functval = num / denom;
      break;
    }

    // linear thermodynamic factor including Debye-Hückel theory
    // 1 + a1*0.5*c^0.5 + a2*c
    case -9:
      functval = 1.0 + functparams[0] * 0.5 * std::pow(concentration, 0.5) +
                 functparams[1] * concentration;
      break;

    // conductivity: own definition which also fulfills the Kohlrausches Square root law
    case -10:
    {
      const double num = functparams[0] * concentration +
                         functparams[1] * std::pow(concentration, 1.5) +
                         functparams[2] * concentration * concentration +
                         functparams[3] * concentration * concentration * concentration;
      const double denom =
          (1.0 + functparams[4] * concentration * concentration +
              functparams[5] * concentration * concentration * concentration * concentration);
      // functparams[0]*(functparams[1]*concentration/nenner) + 0.01 -> constant level 0.01 deleted
      // since it does not have a physical meaning (28.04.2014)
      functval = num / denom;
      break;
    }

    // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann, Kee
    // (2012) kappa = a0*c*exp(a1*c^a2)
    case -11:
    {
      // safety check
      if (concentration < 1.e-12)
        dserror("Concentration value %lf is zero or negative!", concentration);

      const double exponent = functparams[1] * std::pow(concentration, functparams[2]);

      // safety check
      if (exponent > 20.)
        dserror("Overflow detected during conductivity evaluation! Exponent is too large: %lf",
            exponent);

      functval = functparams[0] * concentration * std::exp(exponent);

      break;
    }
    // diffusion coefficient based on a function defined in
    // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion Coefficients in
    // LiPF6 Electrolytic Solutions Journal of The Electrochemical Society, 2008, 155, F13-F16 diff
    // = a0*exp(-a1*c^a2)
    case -12:
    {
      functval = functparams[0] * std::exp(functparams[1] * concentration);
      break;
    }
    case -13:
    {
      // TDF based on a function defined in
      // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
      // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
      // sqrt(c)/(1+a2*sqrt(c)) + a2*c
      functval = 1.0 -
                 (0.5 * functparams[0] * std::pow(concentration, 0.5)) /
                     (std::pow((1 + functparams[1] * std::pow(concentration, 0.5)), 2)) +
                 functparams[2] * concentration;
      break;
    }
    default:
    {
      dserror("Curve number %i is not implemented!", functnr);
      break;
    }
  }

  return functval;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::EvalFirstDerivFunctValue(
    const int functnr, const double concentration, const std::vector<double>& functparams) const
{
  double firstderivfunctval(0.0);

  switch (functnr)
  {
    // d/dc: a0
    case -1:
      firstderivfunctval = 0.0;
      break;

    // d/dc: a0 + a1*c
    case -2:
      firstderivfunctval = functparams[1];
      break;

    // d/dc: a0 + a1*c + a2*c^2
    case -3:
      firstderivfunctval = functparams[1] + 2 * functparams[2] * concentration;
      break;

    // d/dc: a0 + c^a1
    case -4:
      firstderivfunctval =
          functparams[0] * functparams[1] * std::pow(concentration, functparams[1] - 1.0);
      break;

    // d/dc: conductivity
    case -5:
    {
      const double nenner =
          (1.0 + functparams[2] * concentration * concentration -
              functparams[3] * concentration * concentration * concentration * concentration);
      const double nennernenner = nenner * nenner;
      firstderivfunctval =
          functparams[0] * ((functparams[1] * nenner - functparams[1] * concentration *
                                                           (2 * functparams[2] * concentration -
                                                               4 * functparams[3] * concentration *
                                                                   concentration * concentration)) /
                               nennernenner);
      break;
    }

    // d/dc: a0*c + a1*c^1.5 + a2*c^3
    case -6:
      firstderivfunctval = functparams[0] + 1.5 * functparams[1] * std::pow(concentration, 0.5) +
                           3 * functparams[2] * concentration * concentration;
      break;

    // d/dc: a0 + a1*c + a2*c^2 + a3*c^3
    case -7:
      firstderivfunctval = functparams[1] + 2 * functparams[2] * concentration +
                           3 * functparams[3] * concentration * concentration;
      break;

    // d/dc: thermodynamic factor Nyman 2008
    case -8:
    {
      const double num = functparams[0] + functparams[1] * concentration +
                         functparams[2] * concentration * concentration;
      const double denom = functparams[3] + functparams[4] * concentration +
                           functparams[5] * concentration * concentration +
                           functparams[6] * concentration * concentration * concentration;
      const double denomdenom = denom * denom;
      const double derivnum = functparams[1] + 2 * functparams[2] * concentration;
      const double derivdenom = functparams[4] + 2 * functparams[5] * concentration +
                                3 * functparams[6] * concentration * concentration;
      firstderivfunctval = (derivnum * denom - num * derivdenom) / denomdenom;
      break;
    }

    // linear thermodynamic factor including Debye-Hückel theory
    // d/dc: 1 + a1*0.5*c^0.5 + a2*c
    case -9:
      firstderivfunctval =
          functparams[0] * 0.5 * 0.5 * std::pow(concentration, -0.5) + functparams[1];
      break;

    // d/dc: conductivity: own definition which also fulfills the Kohlrausches Square root law
    case -10:
    {
      const double num = functparams[0] * concentration +
                         functparams[1] * std::pow(concentration, 1.5) +
                         functparams[2] * concentration * concentration +
                         functparams[3] * concentration * concentration * concentration;
      const double denom =
          (1.0 + functparams[4] * concentration * concentration +
              functparams[5] * concentration * concentration * concentration * concentration);
      const double denomdenom = denom * denom;
      const double derivnum = functparams[0] + 1.5 * functparams[1] * std::pow(concentration, 0.5) +
                              2.0 * functparams[2] * concentration +
                              3.0 * functparams[3] * concentration * concentration;
      const double derivdenom =
          2.0 * functparams[4] * concentration +
          4.0 * functparams[5] * concentration * concentration * concentration;
      firstderivfunctval = ((derivnum * denom - num * derivdenom) / denomdenom);
      break;
    }

    // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann, Kee
    // (2012) d/dc: kappa = a0*c*exp(a1*c^a2)
    case -11:
    {
      // safety check
      if (concentration < 1.0e-12)
        dserror("Concentration value %lf is zero or negative!", concentration);

      const double exponent = functparams[1] * std::pow(concentration, functparams[2]);

      // safety check
      if (std::abs(exponent) > 20.0)
        dserror(
            "Overflow detected during conductivity evaluation! Absolute value of exponent is too "
            "large: %lf",
            exponent);

      firstderivfunctval =
          functparams[0] * std::exp(exponent) *
          (1 + functparams[1] * functparams[2] * std::pow(concentration, functparams[2]));

      break;
    }
    // diffusion coefficient based on a function defined in
    // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion Coefficients in
    // LiPF6 Electrolytic Solutions Journal of The Electrochemical Society, 2008, 155, F13-F16 diff
    // = a0*exp(-a1*c^a2) deriv (diff) = a0*a1*exp(a1*c^a2)
    case -12:
    {
      firstderivfunctval =
          functparams[0] * functparams[1] * std::exp(functparams[1] * concentration);
      break;
    }
    case -13:
    {
      // TDF based on a function defined in
      // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
      // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
      // sqrt(c)/(1+a2*sqrt(c)) + a2*c
      firstderivfunctval = -(0.25 * functparams[0] * std::pow(concentration, -0.5)) /
                               (std::pow((1 + functparams[1] * std::pow(concentration, 0.5)), 2)) +
                           (0.5 * functparams[0] * functparams[1]) /
                               (std::pow((1 + functparams[1] * std::pow(concentration, 0.5)), 3)) +
                           functparams[2];
      break;
    }
    default:
    {
      dserror("Curve number %i is not implemented!", functnr);
      break;
    }
  }

  return firstderivfunctval;
}
