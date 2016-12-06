/*!----------------------------------------------------------------------
\file combust_time_integration_element.cpp
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and ale displacement.

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "combust_time_integration_element.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> COMBUST::TIMEINT::GetOldPartOfRighthandside(
    const LINALG::Matrix<3,1>&                veln,
    const LINALG::Matrix<3,1>&                velnm,
    const LINALG::Matrix<3,1>&                accn,
    const INPAR::FLUID::TimeIntegrationScheme timealgo,
    const double                              dta,
    const double                              theta
)
{
  LINALG::Matrix<3,1> hist;
  switch (timealgo)
  {
  case INPAR::FLUID::timeint_stationary:
  case INPAR::FLUID::timeint_afgenalpha:
  {
    hist = 0.0;
    break;
  }
  case INPAR::FLUID::timeint_one_step_theta:
  {
    const int nsd = 3;
    for(int isd = 0; isd < nsd; ++isd)
    {
      hist(isd) = veln(isd) + dta*(1.0-theta)*accn(isd);
    }
    break;
  }
  case INPAR::FLUID::timeint_bdf2:
  {
    const int nsd = 3;
    for(int isd = 0; isd < nsd; ++isd)
    {
      hist(isd) = (4.0/3.0)*veln(isd) - (1.0/3.0)*velnm(isd);
    }
    break;
  }
  default:
    dserror("Time integration scheme unknown!");
    exit(1);
  }
  return hist;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> COMBUST::UTILS::GetOldPartOfRighthandside(
    const LINALG::Matrix<3,1>&                 veln,
    const LINALG::Matrix<3,1>&                 velnm,
    const LINALG::Matrix<3,1>&                 accn,
    const INPAR::FLUID::TimeIntegrationScheme  timealgo,
    const double                               dta,
    const double                               theta
)
{
  LINALG::Matrix<3,1> hist;
  switch (timealgo)
  {
    case INPAR::FLUID::timeint_stationary:
    case INPAR::FLUID::timeint_afgenalpha:
    {
      hist = 0.0;
      break;
    }
    case INPAR::FLUID::timeint_one_step_theta:
    {
      const int nsd = 3;
      for(int isd = 0; isd < nsd; ++isd)
      {
        hist(isd) = veln(isd) + dta*(1.0-theta)*accn(isd);
      }
      break;
    }
    case INPAR::FLUID::timeint_bdf2:
    {
      const int nsd = 3;
      for(int isd = 0; isd < nsd; ++isd)
      {
        hist(isd) = (4.0/3.0)*veln(isd) - (1.0/3.0)*velnm(isd);
      }
      break;
    }
    default:
      dserror("Time integration scheme unknown!");
      exit(1);
  }
  return hist;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double COMBUST::TIMEINT::ComputeTimeFac(
    const INPAR::FLUID::TimeIntegrationScheme timealgo,
    const double                              dt,
    const double                              theta,
    const double                              ga_alphaF,
    const double                              ga_alphaM,
    const double                              ga_gamma
  )
{
  double timefac;
  switch (timealgo)
  {
  case INPAR::FLUID::timeint_stationary:
    timefac = 1.0;
    break;

  case INPAR::FLUID::timeint_one_step_theta:
    timefac = theta*dt;
    break;

  case INPAR::FLUID::timeint_bdf2:
    timefac = 2.0/3.0 * dt;
    break;

  case INPAR::FLUID::timeint_afgenalpha:
    timefac = ga_alphaF/ga_alphaM * ga_gamma * dt;
    break;

  default:
    dserror("Time integration scheme unknown!");
    exit(1);
  }
  return timefac;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double COMBUST::UTILS::ComputeTimeFac(
    const INPAR::FLUID::TimeIntegrationScheme timealgo,
    const double                              dt,
    const double                              theta
  )
{
  double timefac;
  switch (timealgo)
  {
    case INPAR::FLUID::timeint_stationary:
      timefac = 1.0;
      break;

    case INPAR::FLUID::timeint_one_step_theta:
    case INPAR::FLUID::timeint_afgenalpha:
      timefac = theta*dt;
      break;

    case INPAR::FLUID::timeint_bdf2:
      timefac = 2.0/3.0 * dt;
      break;

    default:
      dserror("Time integration scheme unknown!");
      exit(1);
  }
  return timefac;
}

