/*!----------------------------------------------------------------------
\file scatra_timint_ost.cpp
\brief 

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "scatra_timint_ost.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::TimIntOneStepTheta(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,output)
{
    // -------------------------------------------------------------------
    // get a vector layout from the discretization to construct matching
    // vectors and matrices
    //                 local <-> global dof numbering
    // -------------------------------------------------------------------
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // Vectors passed to the element
    // -----------------------------

    // temporal solution derivative at time n
    phidtn_       = LINALG::CreateVector(*dofrowmap,true);

    return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntOneStepTheta::~TimIntOneStepTheta()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::SetOldPartOfRighthandside()
  {
    // One-step-Theta:   hist_ = phin_ + dt*(1-Theta)*phidtn_

    hist_->Update(1.0, *phin_, dta_*(1.0-theta_), *phidtn_, 0.0);
    return;
  }


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::Update()
{
  // update time derivative of phi
  if (step_ == 1)
  {
    // do just a linear interpolation within the first timestep
    phidtn_->Update( 1.0/dta_,*phinp_,1.0);
    phidtn_->Update(-1.0/dta_,*phin_ ,1.0);
  }
  else
  {
    /*
    One-step-Theta:

      phidt(n+1) = (phi(n+1)-phi(n)) / (Theta * dt(n)) - (1/Theta -1) * phidt(n)

    */

    double fact1 = 1.0/(theta_*dta_);
    double fact2 = (-1.0/theta_) +1.0;

    phidtn_->Update( fact1,*phinp_,-fact1,*phin_ ,fact2);
  }

  // solution of this step becomes most recent solution of the last step
  phin_ ->Update(1.0,*phinp_,0.0);

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::OutputRestart()
{
  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phidtn", phidtn_);
  output_->WriteVector("phin", phin_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepTheta::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed for One-Step-Theta restart
  reader.ReadVector(phinp_,"phinp");
  reader.ReadVector(phin_, "phin");
  reader.ReadVector(phidtn_, "phidtn");

  return;
}

#endif /* CCADISCRET */
