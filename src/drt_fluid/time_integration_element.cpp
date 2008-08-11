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
BlitzVec3 FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
    const BlitzVec3&                    veln,
    const BlitzVec3&                    velnm,
    const BlitzVec3&                    accn,
    const FLUID_TIMEINTTYPE             timealgo,
    const double                        dta,
    const double                        theta
)
{
  const int nsd = 3;
  BlitzVec3 hist = 0.0;
  switch (timealgo)
  {
  case timeint_stationary:
    hist = 0.0;
    break;

  case timeint_one_step_theta:
    for(int isd = 0; isd < nsd; ++isd)
    {
      hist(isd) = veln(isd) + dta*(1.0-theta)*accn(isd);
    }
    break;

  case timeint_bdf2:
    for(int isd = 0; isd < nsd; ++isd)
    {
      hist(isd) = (4.0/3.0)*veln(isd) - (1.0/3.0)*velnm(isd);
    }
    break;
  default:
    dserror("Time integration scheme unknown!");
  }
  return hist;
}

#endif /* CCADISCRET       */
