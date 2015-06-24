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
#include "acou_ele.H"
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
  FillDIRKValues(dyna_,dirk_a_,dirk_b_,dirk_c_,dirk_q_);
  t_ = LINALG::CreateVector(*(discret_->DofRowMap(0)),true);
  velnp_ = LINALG::CreateVector(*(discret_->DofRowMap(0)),true);
  resonly_ = false;
  // that's it. the standard constructor did everything else
} // TimIntImplDIRK

/*----------------------------------------------------------------------*
 |  Time loop (public)                                   schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::Integrate(Teuchos::RCP<Epetra_MultiVector> history, Teuchos::RCP<LINALG::MapExtractor> splitter)
{
  // output of initial field (given by function for purely acoustic simulation or given by optics for PAT simulation)
  Output(history,splitter);

  // evaluate error
  EvaluateErrorComparedToAnalyticalSol();

  // call elements to calculate system matrix/rhs and assemble
  AssembleMatAndRHS(0);

  // apply Dirichlet boundary conditions to system of equations
  ApplyDirichletToSystem(0);

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

    // p-adaptivity
    UpdatePolyAndState();

    // output of solution
    Output(history,splitter);

    // evaluate error
    EvaluateErrorComparedToAnalyticalSol();
  }

  if (!myrank_) printf("\n");

  return;
} // Integrate


/*----------------------------------------------------------------------*
 |  Update Vectors (public)                              schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::TimeUpdate()
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::TimeUpdate");

  velnp_->Update(1.0,*t_,0.0);

  return;
} // TimeUpdate

/*----------------------------------------------------------------------*
 |  Calculate system matrix (public)                     schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::AssembleMatAndRHS(int stage)
{
  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual
  residual_->Scale(0.0);

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------


  if(padaptivity_) resonly_ = false; // TODO: nur fuer erste stage, bzw. update step
  if(!resonly_)
    sysmat_->Zero();

  eleparams.set<int>("sourcefuncno",sourcefuncno_);
  eleparams.set<bool>("resonly",resonly_);
  eleparams.set<int>("stage",stage);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<double>("dt",dtp_*dirk_a_[0][0]);
  eleparams.set<int>("action",ACOU::calc_systemmat_and_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<int>("step",step_);
  eleparams.set<double>("time",time_);
  eleparams.set<double>("timep",time_+dirk_c_[0]*dtp_);

  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);

  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);
  discret_->ClearState();

  if(!resonly_)
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
  if(!resonly_)
  {
    sysmat_->Complete();
    resonly_ = true;
  }

  std::ofstream file("matrix.dat");
  file << *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
  const double cond = LINALG::Condest(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_),
                                      Ifpack_GMRES, 1000, 1e-10);
  std::cout << "Condition number estimate matrix: " << cond << std::endl;

  return;
} // AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | Update interior field (public)                        schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::UpdateInteriorVariablesAndAssemebleRHS(int stage)
{
  Teuchos::ParameterList eleparams;
  eleparams.set<double>("dt",dtp_*dirk_a_[0][0]);

  eleparams.set<int>("sourcefuncno",sourcefuncno_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);
  eleparams.set<int>("stage",stage);
  eleparams.set<double>("time",time_-dtp_+dirk_c_[stage]*dtp_);

  if(stage==dirk_q_-1)
    eleparams.set<double>("timep",time_+dirk_c_[0]*dtp_);
  else
    eleparams.set<double>("timep",time_-dtp_+dirk_c_[stage+1]*dtp_);

  eleparams.set<bool>("padaptivity",padaptivity_);
  if(stage==dirk_q_-1)eleparams.set<double>("padaptivitytol",padapttol_);
  Teuchos::RCP<std::vector<double> > elevals;
  if(errormaps_)
    elevals = Teuchos::rcp(new std::vector<double>(discret_->NumGlobalElements(),0.0));
  eleparams.set<Teuchos::RCP<std::vector<double> > >("elevals",elevals);

  eleparams.set<int>("action",ACOU::update_secondary_solution_and_calc_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

  discret_->SetState(0,"trace",t_);

  residual_->Scale(0.0);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);
  bool resonly_ = true;
  eleparams.set<bool>("resonly",resonly_);
  eleparams.set<bool>("allelesequal",allelesequal_);

  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  // update the error vector
  if(errormaps_)
  {
    std::vector<double> localvals = *(elevals.get());
    for(int el=0; el<discret_->NumMyRowElements(); ++el)
      error_->ReplaceMyValue(el,0,localvals[error_->Map().GID(el)]);
  }

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
  velnp_->Update(1.0,*veln_,0.0);

  // loop over all stages of the DIRK scheme
  for(int i=0; i<dirk_q_; ++i)
  {
    if(!padaptivity_) discret_->SetState("trace",veln_);

    // solve the linear equation
    const double tcpusolve = Teuchos::Time::wallTime();

    //if(step_<=1)
    //{
       //std::cout<<"LINALG::Condest  "<<LINALG::Condest(dynamic_cast<LINALG::SparseMatrix&>(*sysmat_),Ifpack_GMRES,1000,1e-6)<<std::endl;
      //Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_)->EpetraMatrix()->Print(std::cout);
      //LINALG::PrintMatrixInMatlabFormat("xxx_mat.mat",*(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_)->EpetraMatrix()));
    //}
    solver_->Solve(sysmat_->EpetraOperator(),t_,residual_,true,false,Teuchos::null);
    dtsolve_ += Teuchos::Time::wallTime()-tcpusolve;
    // update interior variables
    double tcpuele = Teuchos::Time::wallTime();
    UpdateInteriorVariablesAndAssemebleRHS(i);
    ApplyDirichletToSystem(i);
    dtele_ += Teuchos::Time::wallTime()-tcpuele;

  } // for(unsigned int i=0; i<dirk_q_; ++i)

  return;
} // DIRKSolve

/*----------------------------------------------------------------------*
 |  Dirichlet function (public)                    schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::ApplyDirichletToSystem(int stage)
{
  TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
  Teuchos::ParameterList params;
  params.set<double>("total time",time_);
  discret_->EvaluateDirichlet(params,zeros_,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  LINALG::ApplyDirichlettoSystem(sysmat_,t_,residual_,Teuchos::null,zeros_,*(dbcmaps_->CondMap()));
  return;
} // ApplyDirichletToSystem


/*----------------------------------------------------------------------*
 | P-Adaptivity                                          schoeder 07/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplDIRK::UpdatePolyAndState()
{
  /* This function serves to supply all required steps for p-adaptivity. What we need
   * is the following:
   * 1.) Do the local postprocessing, calculate delta_k (this is the amount the polynomial shape function needs to change)
   * 2.) Several things
   *     - Update the degree
   *     - Map / project the values
   *     - Rebuild the distributed vectors
   *     - Fill the distributed vectors
   * 3.) Do the next time step
   * UpdateInteriorVariables and ComputeResidual have to be separated (or not, if the element are samrt)
   */

  // 1.) can be fully done by the elements: if p-adaptivity is desired, then the elements
  //     do not store the error in the error vector, but the delta_k!
  if(!padaptivity_) return;

  for(int i=0; i<discret_->NumMyRowElements(); ++i)
    dynamic_cast<DRT::ELEMENTS::Acou*>(discret_->lRowElement(i))->SetDegree(int(error_->operator [](i)));

  // actually we don't want an entire FillComplete call, since nodes and elements remain the same
  // we only want the face and internal dofs and the faces are rebuild
  discret_->BuildFaces();
  discret_->BuildFaceRowMap();
  discret_->BuildFaceColMap();
  discret_->AssignDegreesOfFreedom(0);

  // update maps for global vectors
  velnp_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  residual_.reset(new Epetra_Vector(*(discret_->DofRowMap())));
  t_.reset(new Epetra_Vector(*(discret_->DofRowMap())));

  sysmat_ = Teuchos::null;
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()),108,false,true));

  // now we have to call the calculation of the residual, because we skipped it in
  // UpdateInteriorVariablesAndComputeResidual
  AssembleMatAndRHS(0);

  return;
}

/*----------------------------------------------------------------------*
 | Return the name of the time integrator       (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
std::string ACOU::TimIntImplDIRK::Name()
{
  std::string s = DIRKTypeToString(dyna_);
  return s;
} // Name
