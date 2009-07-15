/*!----------------------------------------------------------------------
\file time_integration_element.cpp
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and ale displacement.

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "time_integration_element.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
    const LINALG::Matrix<3,1>&                 veln,
    const LINALG::Matrix<3,1>&                 velnm,
    const LINALG::Matrix<3,1>&                 accn,
    const FLUID_TIMEINTTYPE             timealgo,
    const double                        dta,
    const double                        theta
)
{
  LINALG::Matrix<3,1> hist;
  switch (timealgo)
  {
  case timeint_stationary:
    hist = 0.0;
    break;

  case timeint_one_step_theta:
  {
    const int nsd = 3;
    for(int isd = 0; isd < nsd; ++isd)
    {
      hist(isd) = veln(isd) + dta*(1.0-theta)*accn(isd);
    }
    break;
  }
  case timeint_bdf2:
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
double FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(
    const FLUID_TIMEINTTYPE    timealgo,
    const double               dt,
    const double               theta
    )
{
  double timefac;
  switch (timealgo)
  {
  case timeint_stationary:
    timefac = 1.0;
    break;

  case timeint_one_step_theta:
    timefac = theta*dt;
    break;

  case timeint_bdf2:
    timefac = 2.0/3.0 * dt;
    break;

  default:
    dserror("Time integration scheme unknown!");
    exit(1);
  }
  return timefac;
}

#endif /* CCADISCRET       */
