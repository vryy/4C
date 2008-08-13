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

  Stationary:

                 hist_ = 0.0

  One-step-Theta:

                 hist_ = veln_ + dt*(1-Theta)*accn_


  BDF2: for constant time step:

                 hist_ = 4/3 veln_ - 1/3 velnm_

  */
  switch (timealgo)
  {
  case timeint_stationary: /* One step Theta time integration */
    hist->Scale(0.0);
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
    const double                         dta,
    const double                         dtp,
    Teuchos::RCP<Epetra_Vector>&         velnp
)
{
  const double fact1 = dta*(1.0+dta/dtp);
  const double fact2 = DSQR(dta/dtp);

  velnp->Update( fact1,*accn ,1.0);
  velnp->Update(-fact2,*veln ,1.0);
  velnp->Update( fact2,*velnm,1.0);

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::TIMEINT_THETA_BDF2::CalculateAcceleration(
    const Teuchos::RCP<Epetra_Vector>&   velnp,
    const Teuchos::RCP<Epetra_Vector>&   veln,
    const Teuchos::RCP<Epetra_Vector>&   velnm,
    const Teuchos::RCP<Epetra_Vector>&   accn,
    const FLUID_TIMEINTTYPE              timealgo,
    const int                            step,
    const double                         theta,
    const double                         dta,
    const double                         dtp,
    Teuchos::RCP<Epetra_Vector>&         accnp
)
{

  // in the first step, we have no old acceleration values
  if (step == 1)
  {
    // do just a linear interpolation within the first timestep
    accnp->Update( 1.0/dta,*velnp,-1.0/dta,*veln, 0.0);
  }
  else
  {
    /*

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
    case timeint_one_step_theta: /* One step Theta time integration */
    {
      const double fact1 = 1.0/(theta*dta);
      const double fact2 =-1.0/theta +1.0;   /* = -1/Theta + 1 */
      
      accnp->Update( fact1,*velnp,0.0);
      accnp->Update(-fact1,*veln ,1.0);
      accnp->Update( fact2,*accn ,1.0);
      
      break;
    }
    case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    {
      if (dta*dtp < EPS15)
        dserror("Zero time step size!!!!!");
      const double sum = dta + dtp;
      
      accnp->Update((2.0*dta+dtp)/(dta*sum),*velnp,
          - sum /(dta*dtp),*veln ,0.0);
      accnp->Update(dta/(dtp*sum),*velnm,1.0);
      break;
    }
    default:
      dserror("Time integration scheme unknown for mass rhs!");
    }
  }
  
  // no accelerations for stationary problems
  if (timealgo == timeint_stationary)
  {
    accnp->PutScalar(0.0);
  }

  return;
}


#endif /* CCADISCRET       */
