/*----------------------------------------------------------------------*/
/*!
\file sti_algorithm.cpp

\brief monolithic algorithm for scatra-thermo interaction

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "sti_algorithm.H"

/*--------------------------------------------------------------------------------*
 | constructor                                                         fang 04/15 |
 *--------------------------------------------------------------------------------*/
STI::Algorithm::Algorithm(
    const Epetra_Comm&              comm,          //! communicator
    const Teuchos::ParameterList&   stidyn,        //! parameter list for scatra-thermo interaction
    const Teuchos::ParameterList&   scatradyn,     //! scalar transport parameter list for scatra and thermo fields
    const Teuchos::ParameterList&   solverparams   //! solver parameter list
    ) :
    // instantiate base class
    AlgorithmBase(comm,scatradyn),

    scatra_(Teuchos::null),
    thermo_(Teuchos::null),
    stiparameters_(Teuchos::rcp(new Teuchos::ParameterList(stidyn))),
    fieldparameters_(Teuchos::rcp(new Teuchos::ParameterList(scatradyn))),
    iter_(0),
    itermax_(fieldparameters_->sublist("NONLINEAR").get<int>("ITEMAX")),
    itertol_(fieldparameters_->sublist("NONLINEAR").get<double>("CONVTOL")),
    restol_(fieldparameters_->sublist("NONLINEAR").get<double>("ABSTOLRES")),

    // initialize global map extractor
    maps_(Teuchos::rcp(new LINALG::MapExtractor(
        *LINALG::MergeMap(
            *scatra_->Discretization()->DofRowMap(),
            *thermo_->Discretization()->DofRowMap(),
            false
            ),
        thermo_->DofRowMap(),
        scatra_->DofRowMap()
        ))),

    // initialize global system matrix
    systemmatrix_(Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
        *maps_,
        *maps_,
        81,
        false,
        true
        ))),

    // initialize scatra-thermo block of global system matrix
    scatrathermoblock_(Teuchos::rcp(new LINALG::SparseMatrix(
        *scatra_->Discretization()->DofRowMap(),
        81,
        true,
        true
        ))),

    // initialize thermo-scatra block of global system matrix
    thermoscatrablock_(Teuchos::rcp(new LINALG::SparseMatrix(
        *thermo_->Discretization()->DofRowMap(),
        81,
        true,
        true
        ))),

    // initialize global increment vector for Newton-Raphson iteration
    increment_(LINALG::CreateVector(
        *DofRowMap(),
        true
        )),

    // initialize timer for Newton-Raphson iteration
    timer_(Teuchos::rcp(new Epetra_Time(comm))),

    // initialize global residual vector
    residual_(LINALG::CreateVector(
        *DofRowMap(),
        true
        )),

    dtsolve_(0.),

    // initialize algebraic solver for global system of equations
    solver_(Teuchos::rcp(new LINALG::Solver(
        solverparams,
        comm,
        DRT::Problem::Instance()->ErrorFile()->Handle()
        ))),

    // initialize L2 norms for Newton-Raphson convergence check
    scatradofnorm_(0.),
    scatraresnorm_(0.),
    scatraincnorm_(0.),
    thermodofnorm_(0.),
    thermoresnorm_(0.),
    thermoincnorm_(0.)
{
  // check maps from scatra and thermo discretizations
  if(scatra_->Discretization()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scatra discretization does not have any degrees of freedom!");
  if(thermo_->Discretization()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Thermo discretization does not have any degrees of freedom!");

  // check global map extractor
  maps_->CheckForValidMapExtractor();

  // additional safety check
  if(!scatra_->IsIncremental())
    dserror("Must have incremental solution approach for scatra-thermo interaction!");

  // initialize scatra time integrator
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(*fieldparameters_,*fieldparameters_,solverparams))->ScaTraField();

  // modify field parameters for thermo field
  ModifyFieldParametersForThermoField();

  // initialize thermo time integrator
  thermo_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(*fieldparameters_,*fieldparameters_,solverparams,"thermo"))->ScaTraField();

  return;
} // STI::Algorithm::Algorithm


/*----------------------------------------------------------------------*
 | assemble global system of equations                       fang 07/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::AssembleMatAndRHS()
{
  // build system matrix and residual for scatra field
  scatra_->PrepareLinearSolve();

  // build system matrix and residual for thermo field
  thermo_->PrepareLinearSolve();

  // build global system matrix
  systemmatrix_->Assign(0,0,LINALG::View,*scatra_->SystemMatrix());
  systemmatrix_->Assign(0,1,LINALG::View,*scatrathermoblock_);
  systemmatrix_->Assign(1,0,LINALG::View,*thermoscatrablock_);
  systemmatrix_->Assign(1,1,LINALG::View,*thermo_->SystemMatrix());
  systemmatrix_->Complete();

  // create full monolithic rhs vector
  maps_->InsertVector(scatra_->Residual(),0,residual_);
  maps_->InsertVector(thermo_->Residual(),1,residual_);

  return;
} // STI::Algorithm::AssembleMatAndRHS


/*----------------------------------------------------------------------*
 | global map of degrees of freedom                          fang 04/15 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& STI::Algorithm::DofRowMap() const
{
  return maps_->FullMap();
}  // STI::Algorithm::DofRowMap()


/*-------------------------------------------------------------------------------------*
 | pass scatra degrees of freedom to thermo discretization and vice versa   fang 04/15 |
 *-------------------------------------------------------------------------------------*/
void STI::Algorithm::ExchangeStateVectors()
{
  // pass scatra degrees of freedom to thermo discretization and vice versa
  thermo_->Discretization()->SetState(1,"scatra",scatra_->Phiafnp());
  scatra_->Discretization()->SetState(1,"thermo",thermo_->Phiafnp());

  return;
} // STI::Algorithm::ExchangeStateVectors()


/*-----------------------------------------------------------------------*
 | check termination criterion for Newton-Raphson iteration   fang 04/15 |
 *-----------------------------------------------------------------------*/
const bool STI::Algorithm::ExitNewtonRaphson()
{
  // initialize exit flag
  bool exit(false);

  // compute vector norms for convergence check
  scatra_->Phinp()->Norm2(&scatradofnorm_);
  maps_->ExtractVector(residual_,0)->Norm2(&scatraresnorm_);
  maps_->ExtractVector(increment_,0)->Norm2(&scatraincnorm_);
  thermo_->Phinp()->Norm2(&thermodofnorm_);
  maps_->ExtractVector(residual_,1)->Norm2(&thermoresnorm_);
  maps_->ExtractVector(increment_,1)->Norm2(&thermoincnorm_);

  // safety checks
  if(std::isnan(scatradofnorm_) or
     std::isnan(scatraresnorm_) or
     std::isnan(scatraincnorm_) or
     std::isnan(thermodofnorm_) or
     std::isnan(thermoresnorm_) or
     std::isnan(thermoincnorm_))
    dserror("Vector norm is not a number!");
  if(std::isinf(scatradofnorm_) or
     std::isinf(scatraresnorm_) or
     std::isinf(scatraincnorm_) or
     std::isinf(thermodofnorm_) or
     std::isinf(thermoresnorm_) or
     std::isinf(thermoincnorm_))
    dserror("Vector norm is infinity!");

  // prevent division by zero
  if(scatradofnorm_ < 1.e-5)
    scatradofnorm_ = 1.e-5;
  if(thermodofnorm_ < 1.e-5)
    scatradofnorm_ = 1.e-5;

  // first Newton-Raphson iteration
  if(iter_ == 1)
  {
    // print first line of convergence table to screen
    // solution increment not yet available during first Newton-Raphson iteration
    if(Comm().MyPID() == 0)
      std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraresnorm_
                << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << thermoresnorm_
                << "   |      --      | "
                << "(       --      , te = "
                << std::setw(10) << std::setprecision(3) << scatra_->DtEle()+thermo_->DtEle() << ")" << std::endl;
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if(Comm().MyPID() == 0)
      std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraresnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraincnorm_/scatradofnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << thermoresnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << thermoincnorm_/thermodofnorm_ << "   | (ts = "
                << std::setw(10) << std::setprecision(3) << dtsolve_ << ", te = "
                << std::setw(10) << std::setprecision(3) << scatra_->DtEle()+thermo_->DtEle() << ")" << std::endl;

    // convergence check
    if(scatraresnorm_ <= itertol_ and
       thermoresnorm_ <= itertol_ and
       scatraincnorm_/scatradofnorm_ <= itertol_ and
       thermoincnorm_/thermodofnorm_ <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional solver calls
  if(scatraresnorm_ < restol_ and thermoresnorm_ < restol_)
    exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without convergence
  if(iter_ == itermax_)
  {
    if(Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;
      std::cout << "|      Newton-Raphson method has not converged after a maximum number of " << std::setw(2) << itermax_ << " iterations!      |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if(exit and Comm().MyPID() == 0)
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl << std::endl;

  return exit;
} // STI::Algorithm::ExitNewtonRaphson()


/*----------------------------------------------------------------------*
 | modify field parameters for thermo field                  fang 06/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::ModifyFieldParametersForThermoField()
{
  // extract parameters for initial temperature field from parameter list for scatra-thermo interaction
  // and overwrite corresponding parameters in parameter list for thermo field
  if(!fieldparameters_->isParameter("INITIALFIELD") or !fieldparameters_->isParameter("INITFUNCNO"))
    dserror("Initial field parameters not properly set in input file section SCALAR TRANSPORT DYNAMIC!");
  if(!stiparameters_->isParameter("THERMO_INITIALFIELD") or !stiparameters_->isParameter("THERMO_INITFUNCNO"))
    dserror("Initial field parameters not properly set in input file section SCALAR TRANSPORT DYNAMIC!");
  fieldparameters_->set<int>("INITIALFIELD",stiparameters_->get<int>("THERMO_INITIALFIELD"));
  fieldparameters_->set<int>("INITFUNCNO",stiparameters_->get<int>("THERMO_INITFUNCNO"));

  return;
} // STI::Algorithm::ModifyFieldParametersForThermoField()


/*----------------------------------------------------------------------*
 | output solution to screen and files                       fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::Output()
{
  // output scatra field
  scatra_->Output();

  // output thermo field
  thermo_->Output();

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // provide scatra and thermo fields with velocities
  scatra_->SetVelocityField(2);
  thermo_->SetVelocityField(2);

  // pass scatra degrees of freedom to thermo discretization and vice versa
  ExchangeStateVectors();

  // prepare time step for scatra and thermo fields
  scatra_->PrepareTimeStep();
  thermo_->PrepareTimeStep();

  // print time step information to screen
  scatra_->PrintTimeStepInfo();

  return;
} // STI::Algorithm::PrepareTimeStep()


/*----------------------------------------------------------------------*
 | read restart data                                         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::ReadRestart(
    int step   //! time step for restart
    )
{
  // read scatra and thermo restart variables
  scatra_->ReadRestart(step);
  thermo_->ReadRestart(step);

  // pass scatra degrees of freedom to thermo discretization and vice versa
  ExchangeStateVectors();

  // set time and time step
  SetTimeStep(scatra_->Time(),step);

  // ToDo: check and remove
  dserror("Restart functionality for scatra-thermo interaction has not been tested yet. Feel free to do it and remove this error.");

  return;
} // STI::Algorithm::ReadRestart


/*----------------------------------------------------------------------*
 | evaluate time step using Newton-Raphson iteration         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::Solve()
{
  // initialize counter for Newton-Raphson iterations
  iter_ = 0;

  // print header of convergence table to screen
  if(Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;
    std::cout << "|- step/max -|- tolerance[norm] -|- scatra-res -|- scatra-inc -|- thermo-res -|- thermo-inc -|" << std::endl;
  }

  // start Newton-Raphson iteration
  while(true)
  {
    // update iteration counter
    iter_ += 1;

    // reset timer
    timer_->ResetStartTime();

    // assemble global system of equations
    AssembleMatAndRHS();

    // safety check
    if(!systemmatrix_->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if(scatra_->FDCheckType() == INPAR::SCATRA::fdcheck_global)
      FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if(ExitNewtonRaphson())
      break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    const double time = timer_->WallTime();

    // solve global system of equations
    // Dirichlet boundary conditions have already been applied to global system of equations
    solver_->Solve(
        systemmatrix_->EpetraOperator(),
        increment_,
        residual_,
        true,
        iter_==1
        );

    // determine time needed for solving global system of equations
    dtsolve_ = timer_->WallTime()-time;

    // update scatra field
    scatra_->UpdateIter(maps_->ExtractVector(increment_,0));
    scatra_->ComputeIntermediateValues();

    // update thermo field
    thermo_->UpdateIter(maps_->ExtractVector(increment_,1));
    thermo_->ComputeIntermediateValues();
  } // Newton-Raphson iteration

  return;
} // STI::Algorithm::Solve


/*----------------------------------------------------------------------*
 | time loop                                                 fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::TimeLoop()
{
  // output initial solution to screen and files
  if(Step() == 0)
    Output();

  // time loop
  while(NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // evaluate time step
    Solve();

    // update scatra and thermo fields
    Update();

    // output solution to screen and files
    Output();
  } // while(NotFinished())

  return;
} // STI::Algorithm::TimeLoop()


/*-------------------------------------------------------------------------*
 | update scatra and thermo fields after time step evaluation   fang 04/15 |
 *-------------------------------------------------------------------------*/
void STI::Algorithm::Update()
{
  // update scatra field
  scatra_->Update();

  // compare scatra field to analytical solution if applicable
  scatra_->EvaluateErrorComparedToAnalyticalSol();

  // update thermo field
  thermo_->Update();

  // compare thermo field to analytical solution if applicable
  thermo_->EvaluateErrorComparedToAnalyticalSol();

  return;
} // STI::Algorithm::Update()


/*---------------------------------------------------------------------------------------------*
 | finite difference check for global system matrix (for debugging only)            fang 07/15 |
 *---------------------------------------------------------------------------------------------*/
void STI::Algorithm::FDCheck()
{
  // initial screen output
  if(Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR STI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(),true));
  maps_->InsertVector(scatra_->Phinp(),0,statenp);
  maps_->InsertVector(thermo_->Phinp(),1,statenp);

  // make a copy of global state vector to undo perturbations later
  Teuchos::RCP<Epetra_Vector> statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if(Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_) != Teuchos::null)
    sysmat_original = (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_)->Merge())))->EpetraMatrix();
  else
    dserror("Global system matrix must be a block sparse matrix!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int col=0; col<sysmat_original->NumGlobalCols(); ++col)
  {
    // fill global state vector with original state variables
    statenp->Update(1.,*statenp_original,0.);

    // impose perturbation
    statenp->SumIntoGlobalValue(col,0,scatra_->FDCheckEps());
    scatra_->Phinp()->Update(1.,*maps_->ExtractVector(statenp,0),0.);
    thermo_->Phinp()->Update(1.,*maps_->ExtractVector(statenp,1),0.);

    // carry perturbation over to state vectors at intermediate time stages if necessary
    scatra_->ComputeIntermediateValues();
    thermo_->ComputeIntermediateValues();

    // calculate element right-hand side vector for perturbed state
    AssembleMatAndRHS();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the system
    // matrix, the second comparison might yield good agreement in spite of the entries being wrong!
    for(int row=0; row<DofRowMap()->NumMyElements(); ++row)
    {
      // get current entry in original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(row);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(row,length,numentries,&values[0],&indices[0]);
      for(int ientry=0; ientry<length; ++ientry)
      {
        if(sysmat_original->ColMap().GID(indices[ientry]) == col)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval = -(*residual_)[row] / scatra_->FDCheckEps() + (*rhs_original)[row] / scatra_->FDCheckEps();

      // confirm accuracy of first comparison
      if(abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        dserror("Finite difference check involves values too close to numerical zero!");

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if(abs(abserr1) > maxabserr)
        maxabserr = abs(abserr1);
      double relerr1(0.);
      if(abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if(abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if(abs(relerr1) > maxrelerr)
        maxrelerr = abs(relerr1);

      // evaluate first comparison
      if(abs(relerr1) > scatra_->FDCheckTol())
      {
        std::cout << "sysmat[" << row << "," << col << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left  = entry - (*rhs_original)[row] / scatra_->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[row] / scatra_->FDCheckEps();

        // confirm accuracy of second comparison
        if(abs(right) > 1.e-17 and abs(right) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if(abs(abserr2) > maxabserr)
          maxabserr = abs(abserr2);
        double relerr2(0.);
        if(abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if(abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if(abs(relerr2) > maxrelerr)
          maxrelerr = abs(relerr2);

        // evaluate second comparison
        if(abs(relerr2) > scatra_->FDCheckTol())
        {
          std::cout << "sysmat[" << row << "," << col << "]-rhs[" << row << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << row << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  Comm().SumAll(&counter,&counterglobal,1);
  double maxabserrglobal(0.);
  Comm().MaxAll(&maxabserr,&maxabserrglobal,1);
  double maxrelerrglobal(0.);
  Comm().MaxAll(&maxrelerr,&maxrelerrglobal,1);

  // final screen output
  if(Comm().MyPID() == 0)
  {
    if(counterglobal)
    {
      printf("--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n",counterglobal);
      dserror("Finite difference check failed for scalar transport system matrix!");
    }
    else
      printf("--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",maxabserrglobal,maxrelerrglobal);
  }

  // undo perturbations of state variables
  scatra_->Phinp()->Update(1.,*maps_->ExtractVector(statenp_original,0),0.);
  scatra_->ComputeIntermediateValues();
  thermo_->Phinp()->Update(1.,*maps_->ExtractVector(statenp_original,1),0.);
  thermo_->ComputeIntermediateValues();

  // recompute system matrix and right-hand side vector based on original state variables
  AssembleMatAndRHS();

  return;
}
