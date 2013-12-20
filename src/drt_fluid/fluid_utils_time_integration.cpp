/*!----------------------------------------------------------------------
\file time_integration_scheme.cpp
\brief routines for fluid (in)stationary time-integration,

     including instationary formulations

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and ale displacement.

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#include "fluid_utils_time_integration.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::ExplicitPredictor(
    const std::string                               predictor,
    const Teuchos::RCP<const Epetra_Vector>    veln,
    const Teuchos::RCP<const Epetra_Vector>    velnm,
    const Teuchos::RCP<const Epetra_Vector>    accn,
    const LINALG::MapExtractor&                velpressplitter,
    const double                               timealgo_constant,
    const double                               dta,
    const double                               dtp,
    const Teuchos::RCP<Epetra_Vector>          velnp,
    const Epetra_Comm&                         comm
)
{

  if(comm.MyPID()==0)
  {
    printf("fluid: using explicit predictor %s",predictor.c_str());
  }

  if (predictor=="steady_state")
  {
    // steady state predictor
    //
    //       n+1    n
    //      u    = u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)

    // this has already been done in TimeUpdate()
  }
  else if(predictor=="zero_acceleration")
  {
    // zero acceleration predictor
    //
    //       n+1    n                   n
    //      u    = u  + (1-gamma)*dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    velnp->Update(1.0,*veln,0.0);

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> inc = velpressplitter.ExtractOtherVector(accn);
    inc->Scale((1.0-timealgo_constant)*dta);

    velpressplitter.AddOtherVector(inc,velnp);
  }
  else if(predictor=="constant_acceleration")
  {
    // constant acceleration predictor
    //
    //       n+1    n         n
    //      u    = u  + dt*acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    velnp->Update(1.0,*veln,0.0);

    Teuchos::RCP<Epetra_Vector> inc = velpressplitter.ExtractOtherVector(accn);
    inc->Scale(dta);

    velpressplitter.AddOtherVector(inc,velnp);
  }
  else if(predictor=="constant_increment")
  {
    // constant increment predictor
    //
    //       n+1      n    n-1
    //      u    = 2*u  - u
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    //
    velnp->Update(1.0,*veln,0.0);

    Teuchos::RCP<Epetra_Vector> un  = velpressplitter.ExtractOtherVector(veln );
    Teuchos::RCP<Epetra_Vector> unm = velpressplitter.ExtractOtherVector(velnm);
    unm->Scale(-1.0);

    velpressplitter.AddOtherVector(un ,velnp);
    velpressplitter.AddOtherVector(unm,velnp);
  }
  else if(predictor=="explicit_second_order_midpoint")
  {
    // the conventional explicit second order predictor (assuming constant dt)
    // also known as leapfrog integration
    /*
    //                        /          n    n-1 \
    //       n+1    n        |      n   u  - u     |
    //      u    = u  + dt * | 2*acc  - ---------  |
    //       (0)             |             dt      |
    //                        \                   /
    // respectively
    //
    //       n+1    n-1               n
    //      u    = u    + 2 * dt * acc
    //       (0)
    //
    //  and
    //
    //       n+1    n
    //      p    = p
    //       (0)
    */
    velnp->Update(1.0,*veln,0.0);

    // split between acceleration and pressure
    Teuchos::RCP<Epetra_Vector> unm = velpressplitter.ExtractOtherVector(velnm);
    Teuchos::RCP<Epetra_Vector> an  = velpressplitter.ExtractOtherVector(accn );

    unm->Update(2.0*dta,*an,1.0);

    velpressplitter.InsertOtherVector(unm,velnp);
  }
  else
    dserror("Unknown fluid predictor %s", predictor.c_str());

  if(comm.MyPID()==0)
  {
    printf("\n");
  }

  return;
}



