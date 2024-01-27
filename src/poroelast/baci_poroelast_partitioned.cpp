/*----------------------------------------------------------------------*/
/*! \file

 \brief  Partitioned poroelasticity algorithm

\level 2

 *-----------------------------------------------------------------------*/

#include "baci_poroelast_partitioned.H"

#include "baci_adapter_fld_poro.H"
#include "baci_adapter_str_fpsiwrapper.H"
#include "baci_global_data.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_structure_aux.H"

BACI_NAMESPACE_OPEN

POROELAST::Partitioned::Partitioned(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter)
    : PoroBase(comm, timeparams, porosity_splitter),
      fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(FluidField()->Velnp())))),
      structincnp_(Teuchos::rcp(new Epetra_Vector(*(StructureField()->Dispnp()))))
{
  const Teuchos::ParameterList& porodyn = GLOBAL::Problem::Instance()->PoroelastDynamicParams();
  // Get the parameters for the ConvergenceCheck
  itmax_ = porodyn.get<int>("ITEMAX");     // default: =10
  ittol_ = porodyn.get<double>("INCTOL");  // default: =1e-6

  fluidveln_ = CORE::LINALG::CreateVector(*(FluidField()->DofRowMap()), true);
  fluidveln_->PutScalar(0.0);
}

void POROELAST::Partitioned::DoTimeStep()
{
  PrepareTimeStep();

  Solve();

  UpdateAndOutput();
}

void POROELAST::Partitioned::SetupSystem() {}  // SetupSystem()

void POROELAST::Partitioned::UpdateAndOutput()
{
  constexpr bool force_prepare = false;
  PrepareOutput(force_prepare);

  Update();

  Output();
}

void POROELAST::Partitioned::Solve()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  if (Step() == 1)
  {
    fluidveln_->Update(1.0, *(FluidField()->Veln()), 0.0);
  }

  while (!stopnonliniter)
  {
    itnum++;

    // store increment from first solution for convergence check (like in
    // elch_algorithm: use current values)
    fluidincnp_->Update(1.0, *FluidField()->Velnp(), 0.0);
    structincnp_->Update(1.0, *StructureField()->Dispnp(), 0.0);

    // get current fluid velocities due to solve fluid step, like predictor in FSI
    // 1. iteration: get velocities of old time step (T_n)
    if (itnum == 1)
    {
      fluidveln_->Update(1.0, *(FluidField()->Veln()), 0.0);
    }
    else  // itnum > 1
    {
      // save velocity solution of old iteration step T_{n+1}^i
      fluidveln_->Update(1.0, *(FluidField()->Velnp()), 0.0);
    }

    // set fluid- and structure-based scalar transport values required in FSI
    SetFluidSolution();

    if (itnum != 1) StructureField()->PreparePartitionStep();
    // solve structural system
    DoStructStep();

    // set mesh displacement and velocity fields
    SetStructSolution();

    // solve scalar transport equation
    DoFluidStep();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }
}

void POROELAST::Partitioned::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  StructureField()->Solve();
}

void POROELAST::Partitioned::DoFluidStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n FLUID SOLVER \n***********************\n";
  }

  // FluidField()->PrepareSolve();
  // Newton-Raphson iteration
  FluidField()->Solve();
}

void POROELAST::Partitioned::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  StructureField()->PrepareTimeStep();
  SetStructSolution();
  FluidField()->PrepareTimeStep();
  SetFluidSolution();
}

bool POROELAST::Partitioned::ConvergenceCheck(int itnum)
{
  // convergence check based on the increment
  bool stopnonliniter = false;

  // variables to save different L2 - Norms
  // define L2-norm of increments and solution
  double fluidincnorm_L2(0.0);
  double fluidnorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double structnorm_L2(0.0);

  // build the current increment
  fluidincnp_->Update(1.0, *(FluidField()->Velnp()), -1.0);
  structincnp_->Update(1.0, *(StructureField()->Dispnp()), -1.0);

  // build the L2-norm of the increment and the solution
  fluidincnp_->Norm2(&fluidincnorm_L2);
  FluidField()->Velnp()->Norm2(&fluidnorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  StructureField()->Dispnp()->Norm2(&structnorm_L2);

  // care for the case that there is (almost) zero solution
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;
  if (structnorm_L2 < 1e-6) structnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout
        << "***********************************************************************************\n";
    std::cout << "    OUTER ITERATION STEP    \n";
    std::cout
        << "***********************************************************************************\n";
    printf("+--------------+------------------------+--------------------+--------------------+\n");
    printf(
        "|-  step/max  -|-  tol      [norm]     -|--  fluid-inc      --|--  disp-inc      --|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         |", itnum, itmax_,
        ittol_, fluidincnorm_L2 / fluidnorm_L2, dispincnorm_L2 / structnorm_L2);
    printf("\n");
    printf("+--------------+------------------------+--------------------+--------------------+\n");
  }

  // converged
  if ((fluidincnorm_L2 / fluidnorm_L2 <= ittol_) && (dispincnorm_L2 / structnorm_L2 <= ittol_))
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf(
          "|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n",
          itnum, itmax_);
      printf(
          "+--------------+------------------------+--------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax_) and
      ((fluidincnorm_L2 / fluidnorm_L2 > ittol_) || (dispincnorm_L2 / structnorm_L2 > ittol_)))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))  // and PrintScreenEvry() and (Step()%PrintScreenEvry()==0))
    {
      printf(
          "|     >>>>>> not converged in itemax steps!                                       |\n");
      printf(
          "+--------------+------------------------+--------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}

Teuchos::RCP<const Epetra_Map> POROELAST::Partitioned::DofRowMapStructure()
{
  return StructureField()->DofRowMap();
}

Teuchos::RCP<const Epetra_Map> POROELAST::Partitioned::DofRowMapFluid()
{
  return FluidField()->DofRowMap();
}

BACI_NAMESPACE_CLOSE
