/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_stat.cpp
\brief solution algorithm for stationary problems

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "scatra_timint_stat.H"
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntStationary::TimIntStationary(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<ParameterList>            extraparams,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,extraparams,output)
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // fine-scale vector
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    fsphinp_ = LINALG::CreateVector(*dofrowmap,true);

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                   gjb 08/08 |
*----------------------------------------------------------------------*/
SCATRA::TimIntStationary::~TimIntStationary()
{
  return;
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetOldPartOfRighthandside()
{
  hist_->PutScalar(0.0);

  return;
}


/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::SetTimeForNeumannEvaluation(
  ParameterList& params)
{
  params.set("total time",time_);
  return;
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                    vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AddNeumannToResidual()
{
  residual_->Update(1.0,*neumann_loads_,1.0);
  return;
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // AVM3 separation
  Sep_->Multiply(false,*phinp_,*fsphinp_);

  // set fine-scale vector
  discret_->SetState("fsphinp",fsphinp_);

  return;
}


/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AddSpecificTimeIntegrationParameters(
  ParameterList& params)
{
  params.set("using stationary formulation",true);
  params.set("using generalized-alpha time integration",false);
  params.set("total time",time_);

  discret_->SetState("hist",hist_);
  discret_->SetState("phinp",phinp_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 09/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // read state vectors that are needed for restart
  reader.ReadVector(phinp_, "phinp");

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart                 gjb 10/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::OutputRestart()
{
  // This feature enables starting a time-dependent simulation from
  // a non-trivial steady-state solution that was calculated before.
  output_->WriteVector("phin", phinp_);
  output_->WriteVector("phinm", phinp_); // for BDF2
  output_->WriteVector("phidtn", zeros_); // for OST

  return;
}


#endif /* CCADISCRET */
