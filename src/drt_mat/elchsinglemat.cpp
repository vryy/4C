/*!----------------------------------------------------------------------*/
/*!
\file elchsinglemat.cpp

\brief abstract interface for electrode and electrolyte materials carrying concentration and
electric potential as degrees of freedom

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

\level 2
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include "elchsinglemat.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 02/15 |
 *----------------------------------------------------------------------*/
MAT::PAR::ElchSingleMat::ElchSingleMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      diffcurve_(matdata->GetInt("DIFFCOEF")),
      diffparanum_(matdata->GetInt("DIFF_PARA_NUM")),
      diffpara_(*matdata->Get<std::vector<double>>("DIFF_PARA")),
      condcurve_(matdata->GetInt("COND")),
      condparanum_(matdata->GetInt("COND_PARA_NUM")),
      condpara_(*matdata->Get<std::vector<double>>("COND_PARA"))
{
  // safety checks
  if (diffparanum_ != (int)diffpara_.size())
    dserror("Mismatch in number of parameters for diffusion coefficient!");
  if (condparanum_ != (int)condpara_.size())
    dserror("Mismatch in number of parameters for conductivity!");
  CheckProvidedParams(diffcurve_, diffpara_);
  CheckProvidedParams(condcurve_, condpara_);

  return;
}


/*---------------------------------------------------------------------------------*
 | check whether number of parameters is consistent with curve number   fang 02/15 |
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
        // linear function: functval=functparams[0]+functparams[1]*cint;
        functionname = "'linear function'";
        nfunctparams = 2;
        break;
      }
      case -3:
      {
        // quadratic function: functval=functparams[0]+functparams[1]*cint+functparams[2]*cint*cint;
        functionname = "'quadratic function'";
        nfunctparams = 3;
        break;
      }
      case -4:
      {
        // power function: functval=functparams[0]*pow(cint,functparams[1]);
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

  return;
}


/*----------------------------------------------------------------------*
 | compute diffusion coefficient according to curve number   fang 02/15 |
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeDiffusionCoefficient(const double cint) const
{
  double diff(0.);

  if (DiffCurve() < 0)
    diff = EvalFunctValue(DiffCurve(), cint, DiffParams());
  else if (DiffCurve() == 0)
    diff = EvalFunctValue(-1, cint, DiffParams());
  else
    diff = DRT::Problem::Instance()->Funct(DiffCurve() - 1).EvaluateTime(cint);

  return diff;
}


/*------------------------------------------------------------------------------------------*
 | compute first derivative of diffusion coefficient according to curve number   fang 02/15 |
 *------------------------------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeFirstDerivDiffCoeff(const double cint) const
{
  double firstderiv(0.);

  if (DiffCurve() < 0)
    firstderiv = EvalFirstDerivFunctValue(DiffCurve(), cint, DiffParams());
  else if (DiffCurve() == 0)
    firstderiv = EvalFirstDerivFunctValue(-1, cint, DiffParams());
  else
    firstderiv =
        (DRT::Problem::Instance()->Funct(DiffCurve() - 1).EvaluateTimeDerivative(cint, 1))[1];

  return firstderiv;
}


/*----------------------------------------------------------------------*
 | compute conductivity according to curve number            fang 02/15 |
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeConductivity(const double cint) const
{
  double cond(0.);

  if (CondCurve() < 0)
    cond = EvalFunctValue(CondCurve(), cint, CondParams());
  else if (CondCurve() == 0)
    cond = EvalFunctValue(-1, cint, CondParams());
  else
    cond = DRT::Problem::Instance()->Funct(CondCurve() - 1).EvaluateTime(cint);

  return cond;
}


/*---------------------------------------------------------------------------------*
 | compute first derivative of conductivity according to curve number   fang 02/15 |
 *---------------------------------------------------------------------------------*/
double MAT::ElchSingleMat::ComputeFirstDerivCond(const double cint) const
{
  double firstderiv(0.);

  if (CondCurve() < 0)
    firstderiv = EvalFirstDerivFunctValue(CondCurve(), cint, CondParams());
  else if (CondCurve() == 0)
    firstderiv = EvalFirstDerivFunctValue(-1, cint, CondParams());
  else
    firstderiv =
        (DRT::Problem::Instance()->Funct(CondCurve() - 1).EvaluateTimeDerivative(cint, 1))[1];

  return firstderiv;
}


/*----------------------------------------------------------------------*
 | evaluate value of predefined function                     fang 02/15 |
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::EvalFunctValue(
    const int functnr, const double cint, const std::vector<double>& functparams) const
{
  double functval(0.);

  switch (functnr)
  {
    // a0
    case -1:
      functval = functparams[0];
      break;

    // a0 + a1*c
    case -2:
      functval = functparams[0] + functparams[1] * cint;
      break;

    // a0 + a1*c + a2*c^2
    case -3:
      functval = functparams[0] + functparams[1] * cint + functparams[2] * cint * cint;
      break;

    // a0*c^a1
    case -4:
      functval = functparams[0] * pow(cint, functparams[1]);
      break;

    // conductivity
    case -5:
    {
      const double nenner =
          (1.0 + functparams[2] * cint * cint - functparams[3] * cint * cint * cint * cint);
      // functparams[0]*(functparams[1]*cint/nenner) + 0.01 -> constant level 0.01 deleted since it
      // does not have a physical meaning (28.04.2014)
      functval = functparams[0] * (functparams[1] * cint / nenner);
      break;
    }

    // a0*c + a1*c^1.5 + a2*c^3
    case -6:
      functval = functparams[0] * cint + functparams[1] * pow(cint, 1.5) +
                 functparams[2] * cint * cint * cint;
      break;

    // a0 + a1*c + a2*c^2 + a3*c^3
    case -7:
      functval = functparams[0] + functparams[1] * cint + functparams[2] * cint * cint +
                 functparams[3] * cint * cint * cint;
      break;

    // thermodynamic factor Nyman 2008
    case -8:
    {
      const double num = functparams[0] + functparams[1] * cint + functparams[2] * cint * cint;
      const double denom = functparams[3] + functparams[4] * cint + functparams[5] * cint * cint +
                           functparams[6] * cint * cint * cint;
      functval = num / denom;
      break;
    }

    // linear thermodynamic factor including Debye-Hückel theory
    // 1 + a1*0.5*c^0.5 + a2*c
    case -9:
      functval = 1.0 + functparams[0] * 0.5 * pow(cint, 0.5) + functparams[1] * cint;
      break;

    // conductivity: own definition which also fulfills the Kohlrausches Square root law
    case -10:
    {
      const double num = functparams[0] * cint + functparams[1] * pow(cint, 1.5) +
                         functparams[2] * cint * cint + functparams[3] * cint * cint * cint;
      const double denom =
          (1.0 + functparams[4] * cint * cint + functparams[5] * cint * cint * cint * cint);
      // functparams[0]*(functparams[1]*cint/nenner) + 0.01 -> constant level 0.01 deleted since it
      // does not have a physical meaning (28.04.2014)
      functval = num / denom;
      break;
    }

    // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann, Kee
    // (2012) kappa = a0*c*exp(a1*c^a2)
    case -11:
    {
      // safety check
      if (cint < 1.e-12) dserror("Concentration value %lf is zero or negative!", cint);

      const double exponent = functparams[1] * pow(cint, functparams[2]);

      // safety check
      if (exponent > 20.)
        dserror("Overflow detected during conductivity evaluation! Exponent is too large: %lf",
            exponent);

      functval = functparams[0] * cint * exp(exponent);

      break;
    }
    // diffusion coefficient based on a function defined in
    // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion Coefficients in
    // LiPF6 Electrolytic Solutions Journal of The Electrochemical Society, 2008, 155, F13-F16 diff
    // = a0*exp(-a1*c^a2)
    case -12:
    {
      functval = functparams[0] * exp(functparams[1] * cint);
      break;
    }
    case -13:
    {
      // TDF based on a function defined in
      // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
      // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
      // sqrt(c)/(1+a2*sqrt(c)) + a2*c
      functval = 1.0 -
                 (0.5 * functparams[0] * pow(cint, 0.5)) /
                     (pow((1 + functparams[1] * pow(cint, 0.5)), 2)) +
                 functparams[2] * cint;
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
 | evaluate first derivative of predefined function          fang 02/15 |
 *----------------------------------------------------------------------*/
double MAT::ElchSingleMat::EvalFirstDerivFunctValue(
    const int functnr, const double cint, const std::vector<double>& functparams) const
{
  double firstderivfunctval(0.);

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
      firstderivfunctval = functparams[1] + 2 * functparams[2] * cint;
      break;

    // d/dc: a0 + c^a1
    case -4:
      firstderivfunctval = functparams[0] * functparams[1] * pow(cint, functparams[1] - 1.0);
      break;

    // d/dc: conductivity
    case -5:
    {
      const double nenner =
          (1.0 + functparams[2] * cint * cint - functparams[3] * cint * cint * cint * cint);
      const double nennernenner = nenner * nenner;
      firstderivfunctval =
          functparams[0] *
          ((functparams[1] * nenner -
               functparams[1] * cint *
                   (2 * functparams[2] * cint - 4 * functparams[3] * cint * cint * cint)) /
              nennernenner);
      break;
    }

    // d/dc: a0*c + a1*c^1.5 + a2*c^3
    case -6:
      firstderivfunctval =
          functparams[0] + 1.5 * functparams[1] * pow(cint, 0.5) + 3 * functparams[2] * cint * cint;
      break;

    // d/dc: a0 + a1*c + a2*c^2 + a3*c^3
    case -7:
      firstderivfunctval =
          functparams[1] + 2 * functparams[2] * cint + 3 * functparams[3] * cint * cint;
      break;

    // d/dc: thermodynamic factor Nyman 2008
    case -8:
    {
      const double num = functparams[0] + functparams[1] * cint + functparams[2] * cint * cint;
      const double denom = functparams[3] + functparams[4] * cint + functparams[5] * cint * cint +
                           functparams[6] * cint * cint * cint;
      const double denomdenom = denom * denom;
      const double derivnum = functparams[1] + 2 * functparams[2] * cint;
      const double derivdenom =
          functparams[4] + 2 * functparams[5] * cint + 3 * functparams[6] * cint * cint;
      firstderivfunctval = (derivnum * denom - num * derivdenom) / denomdenom;
      break;
    }

    // linear thermodynamic factor including Debye-Hückel theory
    // d/dc: 1 + a1*0.5*c^0.5 + a2*c
    case -9:
      firstderivfunctval = functparams[0] * 0.5 * 0.5 * pow(cint, -0.5) + functparams[1];
      break;

    // d/dc: conductivity: own definition which also fulfills the Kohlrausches Square root law
    case -10:
    {
      const double num = functparams[0] * cint + functparams[1] * pow(cint, 1.5) +
                         functparams[2] * cint * cint + functparams[3] * cint * cint * cint;
      const double denom =
          (1.0 + functparams[4] * cint * cint + functparams[5] * cint * cint * cint * cint);
      const double denomdenom = denom * denom;
      const double derivnum = functparams[0] + 1.5 * functparams[1] * pow(cint, 0.5) +
                              2.0 * functparams[2] * cint + 3.0 * functparams[3] * cint * cint;
      const double derivdenom =
          2.0 * functparams[4] * cint + 4.0 * functparams[5] * cint * cint * cint;
      firstderivfunctval = ((derivnum * denom - num * derivdenom) / denomdenom);
      break;
    }

    // conductivity as a function of concentration according to Goldin, Colclasure, Wiedemann, Kee
    // (2012) d/dc: kappa = a0*c*exp(a1*c^a2)
    case -11:
    {
      // safety check
      if (cint < 1.e-12) dserror("Concentration value %lf is zero or negative!", cint);

      const double exponent = functparams[1] * pow(cint, functparams[2]);

      // safety check
      if (abs(exponent) > 20.)
        dserror(
            "Overflow detected during conductivity evaluation! Absolute value of exponent is too "
            "large: %lf",
            exponent);

      firstderivfunctval = functparams[0] * exp(exponent) *
                           (1 + functparams[1] * functparams[2] * pow(cint, functparams[2]));

      break;
    }
    // diffusion coefficient based on a function defined in
    // Stewart, S. G. & Newman, J. The Use of UV/vis Absorption to Measure Diffusion Coefficients in
    // LiPF6 Electrolytic Solutions Journal of The Electrochemical Society, 2008, 155, F13-F16 diff
    // = a0*exp(-a1*c^a2) deriv (diff) = a0*a1*exp(a1*c^a2)
    case -12:
    {
      firstderivfunctval = functparams[0] * functparams[1] * exp(functparams[1] * cint);
      break;
    }
    case -13:
    {
      // TDF based on a function defined in
      // J. Landesfeind, A. Ehrl, M. Graf, W.A. Wall, H.A. Gasteiger: Direct electrochemical
      // determination of activity coefficients in aprotic binary electrolytes TDF = 1.0 - 0.5 a1
      // sqrt(c)/(1+a2*sqrt(c)) + a2*c
      firstderivfunctval = -(0.25 * functparams[0] * pow(cint, -0.5)) /
                               (pow((1 + functparams[1] * pow(cint, 0.5)), 2)) +
                           (0.5 * functparams[0] * functparams[1]) /
                               (pow((1 + functparams[1] * pow(cint, 0.5)), 3)) +
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
