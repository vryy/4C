/*--------------------------------------------------------------------------*/
/*! \file

\brief class for partitioned elastohydrodynamic lubrication (lubrication structure interaction)

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/


#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_lubrication.H"
#include "../drt_adapter/adapter_coupling_ehl_mortar.H"

#include "../drt_lubrication/lubrication_timint_implicit.H"

#include "ehl_partitioned.H"

/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned::Partitioned(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string lubrication_disname)
    : Base(comm, globaltimeparams, lubricationparams, structparams, struct_disname,
          lubrication_disname),
      preincnp_(LINALG::CreateVector(
          *lubrication_->LubricationField()->Discretization()->DofRowMap(0), true)),
      dispincnp_(LINALG::CreateVector(*structure_->DofRowMap(0), true))
{
  // call the EHL parameter lists
  const Teuchos::ParameterList& ehlparams = DRT::Problem::Instance()->ElastoHydroDynamicParams();
  const Teuchos::ParameterList& ehlparamspart =
      DRT::Problem::Instance()->ElastoHydroDynamicParams().sublist("PARTITIONED");

  if (DRT::INPUT::IntegralValue<int>(ehlparams, "DIFFTIMESTEPSIZE"))
  {
    dserror("Different time stepping for two way coupling not implemented yet.");
  }

  // Get the parameters for the ConvergenceCheck
  itmax_ = ehlparams.get<int>("ITEMAX");          // default: =10
  ittol_ = ehlparamspart.get<double>("CONVTOL");  // default: =1e-6

  // no dry contact in partitioned ehl
  if (dry_contact_) dserror("no dry contact model in partitioned ehl");
}

/*----------------------------------------------------------------------*
 | Timeloop for EHL problems                                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::Timeloop()
{
  while (NotFinished())
  {
    PrepareTimeStep();

    OuterLoop();

    UpdateAndOutput();
  }
}


/*----------------------------------------------------------------------*
 | prepare time step                                        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  SetStructSolution(structure_->Dispn());
  structure_->PrepareTimeStep();
  //  SetLubricationSolution(lubrication_->LubricationField()->Quantity()); // todo: what quantity
  lubrication_->LubricationField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*
 | outer Timeloop for EHL without relaxation                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::OuterLoop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  while (stopnonliniter == false)
  {
    itnum++;

    // store pressure from first solution for convergence check (like in
    // elch_algorithm: use current values)
    preincnp_->Update(1.0, *lubrication_->LubricationField()->Prenp(), 0.0);
    dispincnp_->Update(1.0, *structure_->Dispnp(), 0.0);

    // set the external fluid force on the structure, which result from the fluid pressure
    SetLubricationSolution(lubrication_->LubricationField()->Prenp());
    if (itnum != 1) structure_->PreparePartitionStep();
    // solve structural system
    DoStructStep();

    // set mesh displacement, velocity fields and film thickness
    SetStructSolution(structure_->Dispnp());

    // solve lubrication equation and calculate the resulting traction, which will be applied on the
    // solids
    DoLubricationStep();
    // LubricationEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::UpdateAndOutput()
{
  structure_->PrepareOutput();

  structure_->Update();
  lubrication_->LubricationField()->Update();

  lubrication_->LubricationField()->EvaluateErrorComparedToAnalyticalSol();

  structure_->Output();
  lubrication_->LubricationField()->Output();
}


/*----------------------------------------------------------------------*
 | solve structure filed                                    wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_->Solve();
}


/*----------------------------------------------------------------------*
 | solve Lubrication field                                  wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::DoLubricationStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n  LUBRICATION SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                           solve nonlinear
  // -------------------------------------------------------------------
  lubrication_->LubricationField()->Solve();
}


/*----------------------------------------------------------------------*
 | convergence check for both fields (lubrication & structure)          |
 |                                                          wirtz 12/15 |
 *----------------------------------------------------------------------*/
bool EHL::Partitioned::ConvergenceCheck(int itnum)
{
  // convergence check based on the pressure increment
  bool stopnonliniter = false;

  //    | pressure increment |_2
  //  -------------------------------- < Tolerance
  //     | pressure+1 |_2
  //
  // AND
  //
  //    | pressure increment |_2
  //  -------------------------------- < Tolerance
  //             dt * n

  // variables to save different L2 - Norms
  // define L2-norm of incremental pressure and pressure
  double preincnorm_L2(0.0);
  double prenorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current pressure increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  preincnp_->Update(1.0, *(lubrication_->LubricationField()->Prenp()), -1.0);
  dispincnp_->Update(1.0, *(structure_->Dispnp()), -1.0);

  // build the L2-norm of the pressure increment and the pressure
  preincnp_->Norm2(&preincnorm_L2);
  lubrication_->LubricationField()->Prenp()->Norm2(&prenorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero pressure
  if (prenorm_L2 < 1e-6) prenorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout
        << "***********************************************************************************\n";
    std::cout << "    OUTER ITERATION STEP    \n";
    std::cout
        << "***********************************************************************************\n";
    printf(
        "+--------------+---------------------+------------------+-----------------+---------------"
        "-------+------------------+\n");
    printf(
        "|-  step/max  -|-  tol      [norm]  -|-  pressure-inc  -|  disp-inc      -|-  "
        "pressure-rel-inc  -|-  disp-rel-inc  -|\n");
    printf(
        "|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E      |  %10.3E     |  %10.3E          |  "
        "%10.3E      |",
        itnum, itmax_, ittol_, preincnorm_L2 / Dt() / sqrt(preincnp_->GlobalLength()),
        dispincnorm_L2 / Dt() / sqrt(dispincnp_->GlobalLength()), preincnorm_L2 / prenorm_L2,
        dispincnorm_L2 / dispnorm_L2);
    printf("\n");
    printf(
        "+--------------+---------------------+------------------+-----------------+---------------"
        "-------+------------------+\n");
  }

  // converged
  if (((preincnorm_L2 / prenorm_L2) <= ittol_) and ((dispincnorm_L2 / dispnorm_L2) <= ittol_) and
      ((dispincnorm_L2 / Dt() / sqrt(dispincnp_->GlobalLength())) <= ittol_) and
      ((preincnorm_L2 / Dt() / sqrt(preincnp_->GlobalLength())) <= ittol_))
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf(
          "|  Outer Iteration loop converged after iteration %3d/%3d !                             "
          "                            |\n",
          itnum, itmax_);
      printf(
          "+--------------+---------------------+------------------+-----------------+-------------"
          "---------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ((itnum == itmax_) and
      (((preincnorm_L2 / prenorm_L2) > ittol_) or ((dispincnorm_L2 / dispnorm_L2) > ittol_) or
          ((dispincnorm_L2 / Dt() / sqrt(dispincnp_->GlobalLength())) > ittol_) or
          (preincnorm_L2 / Dt() / sqrt(preincnp_->GlobalLength())) > ittol_))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf(
          "|     >>>>>> not converged in itemax steps!                                             "
          "                            |\n");
      printf(
          "+--------------+---------------------+------------------+-----------------+-------------"
          "---------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned EHL solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}
