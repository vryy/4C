/*!----------------------------------------------------------------------
\file combust_utils_time_integration.cpp
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

#include "combust_utils_time_integration.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::UTILS::SetOldPartOfRighthandside(
    const Teuchos::RCP<Epetra_Vector>&   veln,
    const Teuchos::RCP<Epetra_Vector>&   velnm,
    const Teuchos::RCP<Epetra_Vector>&   accn,
    const INPAR::FLUID::TimeIntegrationScheme timealgo,
    const double                         dta,
    const double                         theta,
    Teuchos::RCP<Epetra_Vector>&         hist
)
{
  /*
     for low-Mach-number flow: distinguish momentum and continuity part
     (continuity part only meaningful for low-Mach-number flow)

     Stationary/af-generalized-alpha:

                   mom: hist_ = 0.0
                  (con: hist_ = 0.0)

     One-step-Theta:

                   mom: hist_ = veln_  + dt*(1-Theta)*accn_
                  (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

     BDF2: for constant time step:

                   mom: hist_ = 4/3 veln_  - 1/3 velnm_
                  (con: hist_ = 4/3 densn_ - 1/3 densnm_)

  */
  switch (timealgo)
  {
    case INPAR::FLUID::timeint_stationary: /* Stationary algorithm */
    case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
    case INPAR::FLUID::timeint_npgenalpha:
      hist->PutScalar(0.0);
      break;

    case INPAR::FLUID::timeint_one_step_theta: /* One step Theta time integration */
      hist->Update(1.0, *veln, dta*(1.0-theta), *accn, 0.0);
      break;

    case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      hist->Update(4./3., *veln, -1./3., *velnm, 0.0);
      break;

    default:
    {
      dserror("Time integration scheme unknown!");
      break;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::UTILS::CalculateAcceleration(
    const Teuchos::RCP<const Epetra_Vector>    velnp,
    const Teuchos::RCP<const Epetra_Vector>    veln,
    const Teuchos::RCP<const Epetra_Vector>    velnm,
    const Teuchos::RCP<const Epetra_Vector>    accn,
    const INPAR::FLUID::TimeIntegrationScheme  timealgo,
    const int                                  step,
    const double                               theta,
    const double                               dta,
    const double                               dtp,
    const Teuchos::RCP<Epetra_Vector>          accnp
)
{
    /*

    Following formulations are for n+1; acceleration values, however, are
    directly stored in vectors at time n (velocity has not yet been updated).

    One-step-Theta:

     acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


    BDF2:

                   2*dt(n)+dt(n-1)                  dt(n)+dt(n-1)
     acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
                 dt(n)*[dt(n)+dt(n-1)]              dt(n)*dt(n-1)

                         dt(n)
               + ----------------------- vel(n-1)
                 dt(n-1)*[dt(n)+dt(n-1)]

    */

    switch (timealgo)
    {
      case INPAR::FLUID::timeint_stationary: /* no accelerations for stationary problems*/
      {
        accnp->PutScalar(0.0);
        break;
      }
      case INPAR::FLUID::timeint_one_step_theta: /* One-step-theta time integration */
      {
        const double fact1 = 1.0/(theta*dta);
        const double fact2 =-1.0/theta +1.0;   /* = -1/Theta + 1 */

        accnp->Update( fact1,*velnp,0.0);
        accnp->Update(-fact1,*veln ,1.0);
        accnp->Update( fact2,*accn,1.0);
        break;
      }
      case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        if (dta*dtp < EPS15) dserror("Zero time step size!!!!!");
        const double sum = dta + dtp;

        accnp->Update((2.0*dta+dtp)/(dta*sum),*velnp, -sum/(dta*dtp),*veln ,0.0);
        accnp->Update(dta/(dtp*sum),*velnm,1.0);
        break;
      }
      case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
      case INPAR::FLUID::timeint_npgenalpha:
      {
        // do nothing: new acceleration is calculated at beginning of next time step
        break;
      }
      default:
      {
        dserror("Time integration scheme unknown!");
        break;
      }
    }

  return;
}



