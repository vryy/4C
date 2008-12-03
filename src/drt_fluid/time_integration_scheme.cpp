/*!----------------------------------------------------------------------
\file time_integration_scheme.cpp
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

#include "time_integration_scheme.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::TIMEINT_THETA_BDF2::SetOldPartOfRighthandside(
    const Teuchos::RCP<Epetra_Vector>&   veln,
    const Teuchos::RCP<Epetra_Vector>&   velnm,
    const Teuchos::RCP<Epetra_Vector>&   accn,
    const FLUID_TIMEINTTYPE              timealgo,
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
    case timeint_stationary: /* Stationary algorithm */
    case timeint_afgenalpha: /* Af-generalized-alpha time integration */
      hist->PutScalar(0.0);
      break;

    case timeint_one_step_theta: /* One step Theta time integration */
      hist->Update(1.0, *veln, dta*(1.0-theta), *accn, 0.0);
      break;

    case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      hist->Update(4./3., *veln, -1./3., *velnm, 0.0);
      break;

    default:
      dserror("Time integration scheme unknown!");
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::TIMEINT_THETA_BDF2::ExplicitPredictor(
    const Teuchos::RCP<Epetra_Vector>&   veln,
    const Teuchos::RCP<Epetra_Vector>&   velnm,
    const Teuchos::RCP<Epetra_Vector>&   accn,
    const FLUID_TIMEINTTYPE              timealgo,
    const double                         dta,
    const double                         dtp,
    Teuchos::RCP<Epetra_Vector>&         velnp
)
{
  switch (timealgo)
  {
    case timeint_stationary: /* Stationary algorithm */
      // do nothing
      break;
    case timeint_afgenalpha: /* Generalized-alpha time integration */
    {
      // do nothing for the time being, that is, steady-state predictor
      break;
    }
    case timeint_one_step_theta: /* One step Theta time integration */
    case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    {
      const double fact1 = dta*(1.0+dta/dtp);
      const double fact2 = DSQR(dta/dtp);

      velnp->Update( fact1,*accn ,1.0);
      velnp->Update(-fact2,*veln ,1.0);
      velnp->Update( fact2,*velnm,1.0);
      break;
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::TIMEINT_THETA_BDF2::CalculateAcceleration(
    const Teuchos::RCP<Epetra_Vector>&   velnp,
    const Teuchos::RCP<Epetra_Vector>&   veln,
    const Teuchos::RCP<Epetra_Vector>&   velnm,
    const Teuchos::RCP<Epetra_Vector>&   accnp,
    const FLUID_TIMEINTTYPE              timealgo,
    const int                            step,
    const double                         theta,
    const double                         dta,
    const double                         dtp,
    Teuchos::RCP<Epetra_Vector>&         accn
)
{
  if (step == 1)
  {
    switch (timealgo)
    {
      case timeint_stationary: /* no accelerations for stationary problems*/
      {
        accn->PutScalar(0.0);
        break;
      }
      case timeint_one_step_theta: /* One step Theta time integration */
      case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        // do just a linear interpolation within the first timestep
        accn->Update( 1.0/dta,*velnp,-1.0/dta,*veln, 0.0);
        break;
      }
      case timeint_afgenalpha: /* Af-generalized-alpha time integration */
      {
        accn->Update(1.0,*accnp,0.0); //just update acceleration
        break;
      }
      default:
        dserror("Time integration scheme unknown!");
    }
  }
  else
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

    For low-Mach-number flow, the same is done for density values,
    which are located at the "pressure dofs" of "vede"-vectors.

    */

    switch (timealgo)
    {
      case timeint_stationary: /* no accelerations for stationary problems*/
      {
        accn->PutScalar(0.0);
        break;
      }
      case timeint_one_step_theta: /* One-step-theta time integration */
      {
        const double fact1 = 1.0/(theta*dta);
        const double fact2 =-1.0/theta +1.0;   /* = -1/Theta + 1 */

        accn->Update(fact1,*velnp,-fact1,*veln ,fact2);
        break;
      }
      case timeint_afgenalpha: /* Af-generalized-alpha time integration */
      {
        accn->Update(1.0,*accnp,0.0); //just update acceleration
        break;
      }
      case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        if (dta*dtp < EPS15) dserror("Zero time step size!!!!!");
        const double sum = dta + dtp;

        accn->Update((2.0*dta+dtp)/(dta*sum),*velnp, -sum/(dta*dtp),*veln ,0.0);
        accn->Update(dta/(dtp*sum),*velnm,1.0);
        break;
      }
      default:
        dserror("Time integration scheme unknown!");
    }
  }

  return;
}


#endif /* CCADISCRET       */
