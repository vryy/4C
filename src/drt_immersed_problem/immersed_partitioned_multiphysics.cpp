/*!----------------------------------------------------------------------
\file immersed_partitioned_multiphysics.cpp

\brief base class for all multifield partitioned immersed algorithms

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_multiphysics.H"


IMMERSED::ImmersedPartitionedMultiphysics::ImmersedPartitionedMultiphysics(Teuchos::ParameterList& params, const Epetra_Comm& comm)
{

}// ImmersedPartitionedMultiphysics constructor

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
//  const Teuchos::ParameterList& immerseddyn = DRT::Problem::Instance()->ImmersedMethodParams();
//
//  // Get the top level parameter list
//  Teuchos::ParameterList& nlParams = noxparameterlist_;
//
//  // sublists
//
//  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
//  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
//  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
//
//  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
//
//  // Create printing utilities
//  utils_ = Teuchos::rcp(new NOX::Utils(printParams));
//
//  // ==================================================================
//
//  // log solver iterations
//
//  Teuchos::RCP<std::ofstream> log;
//  if (Comm().MyPID()==0)
//  {
//    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
//    s.append(".iteration");
//    log = Teuchos::rcp(new std::ofstream(s.c_str()));
//    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
//           << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
//           << "# Jacobian       = " << nlParams.get("Jacobian", "None") << "\n"
//           << "# Preconditioner = " << nlParams.get("Preconditioner","None") << "\n"
//           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method","Aitken") << "\n"
//           << "# Predictor      = '" << immerseddyn.sublist("PARTITIONED SOLVER").get<std::string>("PREDICTOR") << "'\n"
//           << "#\n"
//           << "# step | time | time/step | #nliter  |R|  #liter  Residual  Jac  Prec  FD_Res  MF_Res  MF_Jac  User\n"
//      ;
//  }
//
//  Teuchos::Time timer("time step timer");
//
//  // ==================================================================
//
//  while (NotFinished())
//  {
//    // set time step size for this step
//    SetFieldDt();
//
//    // Increment all field counters and predict field values whenever
//    // appropriate.
//    PrepareTimeStep();
//
//    // reset all counters
//    std::fill(counter_.begin(),counter_.end(),0);
//    lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
//    linsolvcount_.resize(0);
//
//    // start time measurement
//    Teuchos::RCP<Teuchos::TimeMonitor> timemonitor = Teuchos::rcp(new Teuchos::TimeMonitor(timer,true));
//
//    /*----------------- CSD - predictor for itnum==0 --------------------*/
//
//    // Begin Nonlinear Solver ************************************
//
//    // Get initial guess.
//    Teuchos::RCP<Epetra_Vector> soln = InitialGuess();
//
//    NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);
//
//    // Create the linear system
//    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
//      CreateLinearSystem(nlParams, interface, noxSoln, utils_);
//
//    // Create the Group
//    Teuchos::RCP<NOX::Epetra::Group> grp =
//      Teuchos::rcp(new NOX::Epetra::Group(printParams, interface, noxSoln, linSys));
//
//    // Convergence Tests
//    Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);
//
//    // Create the solver
//    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,Teuchos::RCP<Teuchos::ParameterList>(&nlParams,false));
//
//    // solve the whole thing
//    NOX::StatusTest::StatusType status = solver->solve();
//
//    if (status != NOX::StatusTest::Converged)
//      dserror("Nonlinear solver failed to converge!");
//
//    // End Nonlinear Solver **************************************
//
//    // Output the parameter list
//    if (utils_->isPrintType(NOX::Utils::Parameters))
//      if (Step()==1 and Comm().MyPID()==0)
//      {
//        utils_->out() << std::endl
//                      << "Final Parameters" << std::endl
//                      << "****************" << std::endl;
//        solver->getList().print(utils_->out());
//        utils_->out() << std::endl;
//      }
//
//    // ==================================================================
//
//    // stop time measurement
//    timemonitor = Teuchos::null;
//
//    if (Comm().MyPID()==0)
//    {
//      (*log) << Step()
//             << "\t" << Time()
//             << "\t" << timer.totalElapsedTime()
//             << "\t" << nlParams.sublist("Output").get("Nonlinear Iterations",0)
//             << "\t" << nlParams.sublist("Output").get("2-Norm of Residual", 0.)
//             << "\t" << lsParams.sublist("Output").get("Total Number of Linear Iterations",0)
//        ;
//      for (std::vector<int>::size_type i=0; i<counter_.size(); ++i)
//      {
//        (*log) << " " << counter_[i];
//      }
//      (*log) << std::endl;
//      log->flush();
//    }
//
//    // ==================================================================
//
//
//    // calculate stresses, strains, energies
//    PrepareOutput();
//
//    // prepare field variables for new time step
//    Update();
//
//    // write current solution
//    Output();
//  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool IMMERSED::ImmersedPartitionedMultiphysics::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
//  const char* flags[] = { "Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL };
//
//  Epetra_Time timer(x.Comm());
//  const double startTime = timer.WallTime();
//
//  if (Comm().MyPID()==0)
//  {
//    utils_->out() << "\n "
//                  << "Global residual calculation"
//                  << ".\n";
//    if (fillFlag!=Residual)
//      utils_->out() << " fillFlag = " << flags[fillFlag] << "\n";
//  }
//
//  // we count the number of times the residuum is build
//  counter_[fillFlag] += 1;
//
//  if (!x.Map().UniqueGIDs())
//    dserror("source map not unique");
//
//
//  // Do the coupling step. The real work is in here.
//  CouplingOp(x,F,fillFlag);
//
//
//  const double endTime = timer.WallTime();
//  if (Comm().MyPID()==0)
//    utils_->out() << "\nTime for residual calculation: " << endTime-startTime << " secs\n\n";
  return true;
}
