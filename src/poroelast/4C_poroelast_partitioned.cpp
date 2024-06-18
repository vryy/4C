/*----------------------------------------------------------------------*/
/*! \file

 \brief  Partitioned poroelasticity algorithm

\level 2

 *-----------------------------------------------------------------------*/

#include "4C_poroelast_partitioned.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

PoroElast::Partitioned::Partitioned(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter)
    : PoroBase(comm, timeparams, porosity_splitter),
      fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(fluid_field()->Velnp())))),
      structincnp_(Teuchos::rcp(new Epetra_Vector(*(structure_field()->Dispnp()))))
{
  const Teuchos::ParameterList& porodyn = Global::Problem::Instance()->poroelast_dynamic_params();
  // Get the parameters for the convergence_check
  itmax_ = porodyn.get<int>("ITEMAX");     // default: =10
  ittol_ = porodyn.get<double>("INCTOL");  // default: =1e-6

  fluidveln_ = Core::LinAlg::CreateVector(*(fluid_field()->dof_row_map()), true);
  fluidveln_->PutScalar(0.0);
}

void PoroElast::Partitioned::do_time_step()
{
  prepare_time_step();

  Solve();

  update_and_output();
}

void PoroElast::Partitioned::SetupSystem() {}  // SetupSystem()

void PoroElast::Partitioned::update_and_output()
{
  constexpr bool force_prepare = false;
  prepare_output(force_prepare);

  update();

  Output();
}

void PoroElast::Partitioned::Solve()
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
    fluidveln_->Update(1.0, *(fluid_field()->Veln()), 0.0);
  }

  while (!stopnonliniter)
  {
    itnum++;

    // store increment from first solution for convergence check (like in
    // elch_algorithm: use current values)
    fluidincnp_->Update(1.0, *fluid_field()->Velnp(), 0.0);
    structincnp_->Update(1.0, *structure_field()->Dispnp(), 0.0);

    // get current fluid velocities due to solve fluid step, like predictor in FSI
    // 1. iteration: get velocities of old time step (T_n)
    if (itnum == 1)
    {
      fluidveln_->Update(1.0, *(fluid_field()->Veln()), 0.0);
    }
    else  // itnum > 1
    {
      // save velocity solution of old iteration step T_{n+1}^i
      fluidveln_->Update(1.0, *(fluid_field()->Velnp()), 0.0);
    }

    // set fluid- and structure-based scalar transport values required in FSI
    set_fluid_solution();

    if (itnum != 1) structure_field()->prepare_partition_step();
    // solve structural system
    do_struct_step();

    // set mesh displacement and velocity fields
    set_struct_solution();

    // solve scalar transport equation
    do_fluid_step();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = convergence_check(itnum);
  }
}

void PoroElast::Partitioned::do_struct_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_field()->Solve();
}

void PoroElast::Partitioned::do_fluid_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n FLUID SOLVER \n***********************\n";
  }

  // fluid_field()->PrepareSolve();
  // Newton-Raphson iteration
  fluid_field()->Solve();
}

void PoroElast::Partitioned::prepare_time_step()
{
  increment_time_and_step();
  print_header();

  structure_field()->prepare_time_step();
  set_struct_solution();
  fluid_field()->prepare_time_step();
  set_fluid_solution();
}

bool PoroElast::Partitioned::convergence_check(int itnum)
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
  fluidincnp_->Update(1.0, *(fluid_field()->Velnp()), -1.0);
  structincnp_->Update(1.0, *(structure_field()->Dispnp()), -1.0);

  // build the L2-norm of the increment and the solution
  fluidincnp_->Norm2(&fluidincnorm_L2);
  fluid_field()->Velnp()->Norm2(&fluidnorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  structure_field()->Dispnp()->Norm2(&structnorm_L2);

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
    if ((Comm().MyPID() == 0))  // and print_screen_evry() and (Step()%print_screen_evry()==0))
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

Teuchos::RCP<const Epetra_Map> PoroElast::Partitioned::DofRowMapStructure()
{
  return structure_field()->dof_row_map();
}

Teuchos::RCP<const Epetra_Map> PoroElast::Partitioned::DofRowMapFluid()
{
  return fluid_field()->dof_row_map();
}

FOUR_C_NAMESPACE_CLOSE
