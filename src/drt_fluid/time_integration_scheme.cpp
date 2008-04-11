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

#include <stdio.h>

#include "time_integration_scheme.H"

#include "../drt_lib/linalg_ana.H"

#include "../drt_xfem/dof_management.H"
#include "fluid_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
void FLUIDTIMEINTEGRATION::SetOldPartOfRighthandside(
        RCP<Epetra_Vector>&         veln,
        RCP<Epetra_Vector>&         velnm,
        RCP<Epetra_Vector>&         accn,
        const FLUID_TIMEINTTYPE     timealgo,
        const double                dta,
        const double                theta,
        RCP<Epetra_Vector>&         hist
        )
{
  /*

  One-step-Theta:

                 hist_ = veln_ + dt*(1-Theta)*accn_


  BDF2: for constant time step:

                 hist_ = 4/3 veln_ - 1/3 velnm_

  */
  switch (timealgo)
  {
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
 |  do explicit predictor step to start nonlinear iteration from a      |
 |  better value                                             gammi 04/07|
 *----------------------------------------------------------------------*/
void FLUIDTIMEINTEGRATION::ExplicitPredictor(
        RCP<Epetra_Vector>&         veln,
        RCP<Epetra_Vector>&         velnm,
        RCP<Epetra_Vector>&         accn,
        const double                dta,
        const double                dtp,
        RCP<Epetra_Vector>&         velnp
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
 | current solution becomes most recent solution of next timestep       |
 |                                                           gammi 04/07|
 *----------------------------------------------------------------------*/
void FLUIDTIMEINTEGRATION::CalculateAcceleration(
        RCP<Epetra_Vector>&         velnp,
        RCP<Epetra_Vector>&         veln,
        RCP<Epetra_Vector>&         velnm,
        const FLUID_TIMEINTTYPE     timealgo,
        const int                   step,
        const double                theta,
        const double                dta,
        const double                dtp,
        RCP<Epetra_Vector>&         accn,
        RCP<Epetra_Vector>&         accnm
        )
{

  // update acceleration
  if (step == 1)
  {
    accnm->PutScalar(0.0);

    // do just a linear interpolation within the first timestep
    accn->Update( 1.0/dta,*velnp,1.0);

    accn->Update(-1.0/dta,*veln ,1.0);

    // ???
    accnm->Update(1.0,*accn,0.0);

  }
  else
  {
    // prev. acceleration becomes (n-1)-accel. of next time step
    accnm->Update(1.0,*accn,0.0);

    /*

    One-step-Theta:

    acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


    BDF2:

                   2*dt(n)+dt(n-1)          dt(n)+dt(n-1)
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

            accn->Update( fact1,*velnp,0.0);
            accn->Update(-fact1,*veln ,1.0);
            accn->Update( fact2,*accnm ,1.0);

            break;
          }
          case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
          {
            if (dta*dtp < EPS15)
              dserror("Zero time step size!!!!!");
            const double sum = dta + dtp;

            accn->Update((2.0*dta+dtp)/(dta*sum),*velnp,
                          - sum /(dta*dtp),*veln ,0.0);
            accn->Update(dta/(dtp*sum),*velnm,1.0);
          }
          break;
          default:
            dserror("Time integration scheme unknown for mass rhs!");
      }
    }

  return;
}


/*----------------------------------------------------------------------*
 |                                                           chfoe 01/08|
 -----------------------------------------------------------------------*/
void FLUIDTIMEINTEGRATION::UpdateGridv(
        RCP<Epetra_Vector>&         dispnp,
        RCP<Epetra_Vector>&         dispn,
        RCP<Epetra_Vector>&         dispnm,
        const int                   order,
        const int                   step,
        const double                theta,
        const double                dta,
        const double                dtp,
        RCP<Epetra_Vector>&         gridv
        )
{
  // get order of accuracy of grid velocity determination
  // from input file data

  switch (order)
  {
    case 1:
      /* get gridvelocity from BE time discretisation of mesh motion:
           -> cheap
           -> easy
           -> limits FSI algorithm to first order accuracy in time

                  x^n+1 - x^n
             uG = -----------
                    Delta t                        */
      gridv->Update(1/dta, *dispnp, -1/dta, *dispn, 0.0);
    break;
    case 2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacemnt
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridv->Update(1.5/dta, *dispnp, -2.0, *dispn, 0.0);
      gridv->Update(0.5, *dispnm, 1.0);
    break;
  }
}



#endif /* CCADISCRET       */
