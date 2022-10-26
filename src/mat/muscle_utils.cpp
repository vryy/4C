/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of utility methods for skeletal muscle


\level 3

*/
/*-----------------------------------------------------------*/

#include "muscle_utils.H"
#include "dserror.H"
#include "standardtypes_cpp.H"

void MAT::UTILS::MUSCLE::EvaluateLambert(
    const double xi, double &W0, const double tol, const int maxiter)
{
  double W0_old =
      std::numeric_limits<double>::infinity();  // s.t. error is infinite in the beginning
  int numiter = 0;                              // number of iterations

  // Halley's method
  while ((std::abs(W0 - W0_old) / std::abs(W0) > tol) && (numiter <= maxiter))
  {
    numiter++;
    W0_old = W0;
    W0 = W0 - (W0 * std::exp(W0) - xi) /
                  ((std::exp(W0) * (W0 + 1) - (W0 + 2) * (W0 * std::exp(W0) - xi) / (2 * W0 + 2)));
  }

  // error handling
  if (numiter >= maxiter)
  {
    dserror(
        "Maximal number of iterations for evaluation of Lambert W function with Halley's method "
        "exceeded for tolerance %E.",
        tol);
  }
}

double MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(
    const double lambdaM, const double lambdaMin, const double lambdaOpt)
{
  double fxi = 0.0;
  double explambda = std::exp(((2 * lambdaMin - lambdaM - lambdaOpt) * (lambdaM - lambdaOpt)) /
                              (2 * std::pow(lambdaMin - lambdaOpt, 2)));  // prefactor
  if (lambdaM > lambdaMin)
  {
    fxi = ((lambdaM - lambdaMin) / (lambdaOpt - lambdaMin)) * explambda;
  }

  return fxi;
}

double MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(
    const double lambdaM, const double lambdaMin, const double lambdaOpt)
{
  double dFxidLamdaM = 0.0;
  double explambda = std::exp(((2 * lambdaMin - lambdaM - lambdaOpt) * (lambdaM - lambdaOpt)) /
                              (2 * std::pow(lambdaMin - lambdaOpt, 2)));  // prefactor

  if (lambdaM > lambdaMin)
  {
    dFxidLamdaM = ((std::pow(lambdaMin - lambdaM, 2) - std::pow(lambdaMin - lambdaOpt, 2)) /
                      std::pow(lambdaMin - lambdaOpt, 3)) *
                  explambda;
  }

  return dFxidLamdaM;
}

double MAT::UTILS::MUSCLE::EvaluateForceVelocityDependencyBoel(const double lambdaM,
    const double lambdaMOld, const double timestepsize, const double dotLambdaMMin, const double de,
    const double dc, const double ke, const double kc)
{
  // dotLambdaM = (lambdaM_n - lambdaM_{n-1})/dt approximated through BW Euler
  double dotLambdaM = (lambdaM - lambdaMOld) / timestepsize;

  // helper variable
  double ratioDotLambdaM = dotLambdaM / dotLambdaMMin;

  // velocity dependency fv
  double fv = 1.0;
  if (dotLambdaM > 0)
  {
    fv = (1 + de) - de * (1 + ratioDotLambdaM) / (1 - ke * kc * ratioDotLambdaM);
  }
  else
  {
    fv = (1 - dc) + dc * (1 - ratioDotLambdaM) / (1 + kc * ratioDotLambdaM);
  }

  return fv;
}

double MAT::UTILS::MUSCLE::EvaluateDerivativeForceVelocityDependencyBoel(const double lambdaM,
    const double lambdaMOld, const double timestepsize, const double dotLambdaMMin, const double de,
    const double dc, const double ke, const double kc)
{
  // dotLambdaM = (lambdaM_n - lambdaM_{n-1})/dt approximated through BW Euler
  double dotLambdaM = (lambdaM - lambdaMOld) / timestepsize;

  // dDotLambdaMdLambdaM = 1/dt approximated through BW Euler
  double dDotLambdaMdLambdaM = 1 / timestepsize;

  // helper variable
  double ratioDotLambdaM = dotLambdaM / dotLambdaMMin;

  // derivative of fv w.r.t. dotLambdaM
  double dFvdDotLambdaM = 1.0;
  if (dotLambdaM > 0)
  {
    dFvdDotLambdaM =
        -de * ((1 + ke * kc) / (dotLambdaMMin * std::pow((1 - ke * kc * ratioDotLambdaM), 2)));
  }
  else
  {
    dFvdDotLambdaM = -dc * ((1 + kc) / (dotLambdaMMin * std::pow((1 + kc * ratioDotLambdaM), 2)));
  }

  // derivative of fv w.r.t. lambdaM
  double dFvdLambdaM = dFvdDotLambdaM * dDotLambdaMdLambdaM;

  return dFvdLambdaM;
}

double MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressEhret(const double Na,
    const int muTypesNum, const std::vector<double> &rho, const std::vector<double> &I,
    const std::vector<double> &F, const std::vector<double> &T, const int actIntervalsNum,
    const std::vector<double> &actTimes, const std::vector<double> &actValues,
    const double currentTime)
{
  // compute twitch force of motor unit (MU) type iMU
  double t_iMU_jstim = 0;                // time of arriving stimulus jstim for MU type iMU
  double t_end = 0;                      // helper variable
  double ratiotime = 1.0;                // helper variable
  std::vector<double> sumg(muTypesNum);  // superposition of single twitches until current time
  std::vector<double> G(muTypesNum);     // gain function for MU type iMU
  std::vector<double> Ft(muTypesNum);    // twitch force for MU type iMU

  // for all motor unit types i
  for (int iMU = 0; iMU < muTypesNum; ++iMU)
  {
    // for all activation intervals
    for (int actinterval = 0; actinterval < actIntervalsNum; ++actinterval)
    {
      // set time of first stimulus jstim=1 to start time of the current activation interval
      t_iMU_jstim = actTimes[actinterval];

      // determine end time of activation interval
      // if inside of actinterval
      if (currentTime < actTimes[actinterval + 1])
      {
        t_end = currentTime;  // set end time to current simulation time
      }
      // if outside of actinterval
      else
      {
        t_end = actTimes[actinterval + 1];  // set end time to end time of actinterval
      }

      // superposition of single twitches
      // for all impulses from start of actinterval until determined end time of actinterval
      while (t_iMU_jstim < t_end)
      {
        ratiotime = (currentTime - t_iMU_jstim) / T[iMU];

        // add single twitch force response for stimulus jstim and motor unit iMU scaled by
        // percentage activation prescribed in actinterval
        sumg[iMU] += actValues[actinterval] * ratiotime * F[iMU] * std::exp(1 - ratiotime);

        // next impulse jstim at time t_iMU_jstim after stimulation interval I
        t_iMU_jstim += I[iMU];
      }
    }  // end actinterval

    // gain function for MU type iMU
    G[iMU] = (1 - std::exp(-2 * std::pow(T[iMU] / I[iMU], 3))) / (T[iMU] / I[iMU]);

    // twitch force for MU type iMU
    Ft[iMU] = G[iMU] * sumg[iMU];
  }  // end iMU

  // compute force-time/stimulation frequency dependency Poptft
  double Poptft = 0.0;
  for (int iMU = 0; iMU < muTypesNum; ++iMU)
  {
    // sum up twitch forces for all MU types weighted by the respective MU type fraction
    Poptft += Na * Ft[iMU] * rho[iMU];
  }

  return Poptft;
}