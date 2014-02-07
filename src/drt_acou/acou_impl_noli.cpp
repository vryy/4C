/*!----------------------------------------------------------------------
\file acou_impl_noli.cpp
\brief Main control routine for nonlinear acoustic simulations

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "acou_impl_noli.H"
#include "acou_ele_action.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::TimIntImplNoli::TimIntImplNoli(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      )
:AcouImplicitTimeInt(actdis,solver,params,output),
 ittol_(params->get<double>("CONVTOL")),
 itmax_(params->get<int>("ITEMAX")),
 resnorm_(0.0),
 incnorm_(0.0)
{
  veli_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  intveli_ = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
}

/*----------------------------------------------------------------------*
 |  Nonlinear solver loop                                schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplNoli::Solve()
{
  stopnonliniter_ = false;
  int          itnum = 0;
  const double tcpusolve=Teuchos::Time::wallTime();

  // -------------------------------------------------------------------
  // prepare print out
  // -------------------------------------------------------------------
  if (!myrank_)
  {
    printf("+------------+-------------+-------------+-------------+\n");
    printf("|- step/max -|-    tol    -|-    res    -|-    inc    -|\n");
  }

  if(step_>1)
    AssembleMatAndRHS();

  while(stopnonliniter_ == false)
  {
    // update iter count
    itnum++;

    // solve the system
    solver_->Solve(sysmat_->EpetraOperator(),velnp_,residual_,true,itnum==1,Teuchos::null);

    // check for convergence
    stopnonliniter_ = ConvergenceCheck(itnum);

    // update interior variables, calculate new system matrix and assemble new right hand side vector
    UpdateInteriorVariablesAndAssemebleRHS();


    // update solution increment
    IterTraceUpdate();
  }

  if (!myrank_)
    printf("+------------+-------------+-------------+-------------+\n");

  dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;

  return;
} // Solve

/*----------------------------------------------------------------------*
 |  Update the trace increment                           schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplNoli::IterTraceUpdate()
{
  veli_->Update(1.0,*velnp_,0.0);
  intveli_->Update(1.0,*intvelnp_,0.0);
  return;
} // IterTraceUpdate

/*----------------------------------------------------------------------*
 | update everything                                     schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplNoli::UpdateInteriorVariablesAndAssemebleRHS()
{
  dtele_ = 0.0;

  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::UpdateInteriorVariablesAndAssemebleRHS");

  // get cpu time
  const double tcpu = Teuchos::Time::wallTime();

  // create parameterlist
  Teuchos::ParameterList eleparams;

  // fill in parameters and set states needed by elements
  discret_->SetState(1,"intvel",intvelnp_);
  discret_->SetState(1,"intvelm",intveln_);
  eleparams.set<double>("dt",dtp_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);
  eleparams.set<bool>("converged",stopnonliniter_);

  Teuchos::RCP<std::vector<double> > elevals = Teuchos::rcp(new std::vector<double>(discret_->NumGlobalElements(),0.0));
  eleparams.set<Teuchos::RCP<std::vector<double> > >("elevals",elevals);

  eleparams.set<int>("action",ACOU::update_secondary_solution_and_calc_mat_residual_noli);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);

  discret_->SetState("trace",velnp_);
  discret_->SetState("trace_m",veln_);

  residual_->Scale(0.0);
  sysmat_->Zero();

  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);

  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  // update the error vector
  if(errormaps_)
  {
    std::vector<double> localvals = *(elevals.get());
    for(int el=0; el<discret_->NumMyRowElements(); ++el)
      error_->ReplaceGlobalValue(el,0,localvals[error_->Map().GID(el)]);
  }

  // update internal field for parallel usage
  const Epetra_Vector& intvelnpGhosted = *discret_->GetState(1,"intvel");
  for (int i=0; i<intvelnp_->MyLength(); ++i)
    (*intvelnp_)[i] = intvelnpGhosted[intvelnpGhosted.Map().LID(intvelnp_->Map().GID(i))];

  discret_->ClearState();

  // calculate source term for adjoint simulation
  bool resonly = false;
  eleparams.set<bool>("resonly",resonly);
  if(!resonly || adjoint_)
  {
    // absorbing boundary conditions
    std::string condname = "Absorbing";
    std::vector<DRT::Condition*> absorbingBC;
    discret_->GetCondition(condname,absorbingBC);
    if(absorbingBC.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_abc);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }

  sysmat_->Complete();

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpu;

  return;
} // UpdateInteriorVariablesAndAssemebleRHS

/*----------------------------------------------------------------------*
 |  Convergence Check                                    schoeder 01/14 |
 *----------------------------------------------------------------------*/
bool ACOU::TimIntImplNoli::ConvergenceCheck(int itnum)
{
  // calculate norm of residual and of increment
  Teuchos::RCP<Epetra_Vector> res;
  res = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  sysmat_->Multiply(false,*velnp_,*res);
  res->Update(-1.0,*residual_,1.0);
  res->Norm2(&resnorm_);
  Teuchos::RCP<Epetra_Vector> incvel;
  incvel = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  incvel->Update(1.0,*velnp_,0.0);
  incvel->Update(-1.0,*veli_,1.0);
  incvel->Norm2(&incnorm_);

  // output of iteration information
  if(!myrank_)
  {
    printf("|  %3d/%3d   | %10.3E  | %10.3E  | %10.3E  |\n",itnum,itmax_,ittol_,resnorm_,incnorm_);
  }

  // we fulfill given tolerance
  if(resnorm_ <= ittol_)
    return true;

  // maximal number of iterations reached
  if(itnum == itmax_)
  {
    dserror("        >>>>>> not converged in %d steps!",itmax_);
    return true;
  }

  // neither given tolerance is fulfilled nor the maximal number of iterations is reached
  return false;
} // ConvergenceCheck
