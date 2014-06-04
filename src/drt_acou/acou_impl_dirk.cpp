/*!----------------------------------------------------------------------
\file acou_impl_dirk.cpp
\brief

<pre>
Maintainers: Svenja Schoeder
             schoeder@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "acou_impl_dirk.H"
#include "acou_utils.H"
#include "acou_ele_action.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::TimIntImplDIRK::TimIntImplDIRK(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      )
:AcouImplicitTimeInt(actdis,solver,params,output)
{
  // fill the scheme specific coefficients
  FillDIRKValues(dyna_,dirk_a_,dirk_b_,dirk_q_);

  // locate vectors for stage vectors

  s_.resize(dirk_q_);
  y_.resize(dirk_q_);
  f_.resize(dirk_q_);
  t_.resize(dirk_q_);
  ft_.resize(dirk_q_);

  for(unsigned int i=0; i<dirk_q_; ++i)
  {
    s_[i] = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    y_[i] = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    f_[i] = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    t_[i] = LINALG::CreateVector(*(discret_->DofRowMap(0)),true);
    ft_[i] = LINALG::CreateVector(*(discret_->DofRowMap(0)),true);
  }

  // that's it. the standard constructor did everything else
} // TimIntImplDIRK

/*----------------------------------------------------------------------*
 |  Time loop (public)                                   schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::Integrate(Teuchos::RCP<Epetra_MultiVector> history, Teuchos::RCP<LINALG::MapExtractor> splitter)
{
  // write some information for the curious user
  PrintInformationToScreen();

  // output of initial field (given by function for purely acoustic simulation or given by optics for PAT simulation)
  Output(history,splitter);

  // time loop
  while (step_<stepmax_ and time_<maxtime_)
  {
    // increment time and step
    IncrementTimeAndStep();

    // output to screen
    OutputToScreen();

    // assemble, update and solve all stages of DIRK scheme
    Solve();

    // update solution, current solution becomes old solution of next timestep
    TimeUpdate();

    // output of solution
    Output(history,splitter);

    // evaluate error
    EvaluateErrorComparedToAnalyticalSol();
  }

  if (!myrank_) printf("\n");

  return;
} // Integrate

/*----------------------------------------------------------------------*
 |  Calculate system matrix (public)                     schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::AssembleMatAndRHS()
{
  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual
  residual_->Scale(0.0);

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------

  bool resonly = !(!bool(step_-1) || !bool(step_-restart_-1));
  if(!resonly)
    sysmat_->Zero();

  eleparams.set<bool>("resonly",resonly);
  eleparams.set<double>("dt",dtp_*dirk_a_[0][0]);
  eleparams.set<int>("action",ACOU::calc_systemmat_and_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<int>("step",step_);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);

  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  discret_->ClearState();

  if(!resonly)
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
  if(adjoint_)
  {
    std::string condname = "PressureMonitor";
    std::vector<DRT::Condition*> pressuremon;
    discret_->GetCondition(condname,pressuremon);
    if(pressuremon.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_pressuremon);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }

  // finalize the complete matrix
  if(!resonly)
    sysmat_->Complete();

  return;
} // AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | Update interior field (public)                        schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::UpdateInteriorVariables(int stage)
{
  Teuchos::ParameterList eleparams;

  discret_->SetState(1,"intvel",y_[stage]);

  eleparams.set<double>("dt",dtp_*dirk_a_[0][0]);

  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);
  Teuchos::RCP<std::vector<double> > elevals;
  if(errormaps_)
    elevals = Teuchos::rcp(new std::vector<double>(discret_->NumGlobalElements(),0.0));
  eleparams.set<Teuchos::RCP<std::vector<double> > >("elevals",elevals);

  eleparams.set<int>("action",ACOU::update_secondary_solution);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

  discret_->SetState(0,"trace",t_[stage]);

  residual_->Scale(0.0);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);
  bool resonly = true;
  eleparams.set<bool>("resonly",resonly);

  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // update the error vector
  if(errormaps_)
  {
    std::vector<double> localvals = *(elevals.get());
    for(int el=0; el<discret_->NumMyRowElements(); ++el)
      error_->ReplaceMyValue(el,0,localvals[error_->Map().GID(el)]);
  }

  const Epetra_Vector& intvelnpGhosted = *discret_->GetState(1,"intvel");
  for (int i=0; i<intvelnp_->MyLength(); ++i)
    (*(y_[stage]))[i] = intvelnpGhosted[intvelnpGhosted.Map().LID(intvelnp_->Map().GID(i))];

  discret_->ClearState();

  return;
} // UpdateInteriorVariables

/*----------------------------------------------------------------------*
 |  Loop all DIRK stages (public)                        schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::Solve()
{
  dtsolve_ = 0.0;
  dtele_   = 0.0;

  // initialize some vectors
  intvelnp_->Update(1.0,*intveln_,0.0);
  velnp_->Update(1.0,*veln_,0.0);
  s_[0]->Update(1.0/dirk_a_[0][0],*intvelnp_,0.0);
  y_[0]->Update(1.0,*intvelnp_,0.0);

  // loop over all stages of the DIRK scheme
  for(unsigned int i=0; i<dirk_q_; ++i)
  {
    discret_->SetState("trace",veln_);
    s_[i]->Scale(dirk_a_[i][i]);
    discret_->SetState(1,"intvel",s_[i]);

    // call elements to calculate system matrix/rhs and assemble
    double tcpuele = Teuchos::Time::wallTime();
    AssembleMatAndRHS();
    dtele_ += Teuchos::Time::wallTime()-tcpuele;
    s_[i]->Scale(1.0/dirk_a_[i][i]);

    // apply Dirichlet boundary conditions to system of equations
    ApplyDirichletToSystem();

    // solve the linear equation
    const double tcpusolve = Teuchos::Time::wallTime();
    solver_->Solve(sysmat_->EpetraOperator(),t_[i],residual_,true,false,Teuchos::null);
    dtsolve_ += Teuchos::Time::wallTime()-tcpusolve;

    // update interior variables
    y_[i]->Update(dirk_a_[i][i],*s_[i],0.0);
    tcpuele = Teuchos::Time::wallTime();
    UpdateInteriorVariables(i);
    dtele_ += Teuchos::Time::wallTime()-tcpuele;

    // calculate f[i]
    f_[i]->Update(1.0/dirk_a_[i][i]/dtp_,*y_[i],0.0);
    ft_[i]->Update(1.0/dirk_a_[i][i]/dtp_,*t_[i],0.0);
    f_[i]->Update(-1.0/dirk_a_[i][i]/dtp_,*intveln_,1.0);
    ft_[i]->Update(-1.0/dirk_a_[i][i]/dtp_,*veln_,1.0);
    for(unsigned int j=0; j<i; ++j)
    {
      f_[i]->Update(-dirk_a_[dirk_q_-1][j]/dirk_a_[dirk_q_-1][dirk_q_-1],*f_[j],1.0);
      ft_[i]->Update(-dirk_a_[dirk_q_-1][j]/dirk_a_[dirk_q_-1][dirk_q_-1],*ft_[j],1.0);
    }

    //  now y holds the internal variable values for the i-th stage
    intvelnp_->Update(dtp_*dirk_b_[i],*f_[i],1.0);
    velnp_->Update(dtp_*dirk_b_[i],*ft_[i],1.0);

    // calculate s[i]
    if ((i+1) < dirk_q_)
    {
      s_[i+1]->Update(1.0/dirk_a_[i+1][i+1],*intveln_,0.0);
      for(unsigned int j=0; j<=i; ++j)
      {
        s_[i+1]->Update(dirk_a_[i+1][j]/dirk_a_[i+1][i+1]/dirk_a_[j][j],*y_[j],1.0);
        s_[i+1]->Update(-dirk_a_[i+1][j]/dirk_a_[i+1][i+1],*s_[j],1.0);
      }
    }

  } // for(unsigned int i=0; i<dirk_q_; ++i)

  return;
} // DIRKSolve

/*----------------------------------------------------------------------*
 | Return the name of the time integrator       (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
std::string ACOU::TimIntImplDIRK::Name()
{
  std::string s = DIRKTypeToString(dyna_);
  return s;
} // Name
