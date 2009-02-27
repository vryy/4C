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


/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntStationary::TimIntStationary(
  RCP<DRT::Discretization>      actdis,
  RCP<LINALG::Solver>           solver,
  RCP<ParameterList>            params,
  RCP<IO::DiscretizationWriter> output)
: ScaTraTimIntImpl(actdis,solver,params,output)
{
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
 | add parameters specific for time-integration scheme         vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AddSpecificTimeIntegrationParameters(
  ParameterList& params)
{
  params.set("using stationary formulation",true);
  params.set("using generalized-alpha time integration",false);
  params.set("total time",time_);

  discret_->SetState("phinp", phinp_);
  discret_->SetState("densnp",densnp_);
  discret_->SetState("densam",densnp_);

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
  reader.ReadVector(phinp_,"phinp");
  reader.ReadVector(phin_, "phin");

  return;
}


#endif /* CCADISCRET */
