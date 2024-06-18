/*----------------------------------------------------------------------*/
/*! \file

 \brief partitioned two way coupled poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/

#include "4C_poroelast_scatra_part_2wc.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
PoroElastScaTra::PoroScatraPart2WC::PoroScatraPart2WC(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraPart(comm, timeparams),
      scaincnp_(Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())))),
      structincnp_(Teuchos::rcp(new Epetra_Vector(*(poro_field()->structure_field()()->Dispnp())))),
      fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(poro_field()->fluid_field()()->Velnp()))))
{
  if (comm.MyPID() == 0) std::cout << "\n Create PoroScatraPart2WC algorithm ... \n" << std::endl;

  const Teuchos::ParameterList& params = Global::Problem::Instance()->poro_scatra_control_params();
  // Get the parameters for the convergence_check
  itmax_ = params.get<int>("ITEMAX");            // default: =10
  ittol_ = params.get<double>("TOLINC_GLOBAL");  // default: =1e-6
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::Timeloop()
{
  // initial_calculations();

  while (NotFinished())
  {
    prepare_time_step();

    Solve();

    prepare_output();

    update();

    output();
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::read_restart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    SetScatraSolution();
    SetPoroSolution();

    poro_field()->read_restart(restart);
    ScaTraField()->read_restart(restart);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (poro_field()->HasSubmeshes())
      replace_dof_sets(structure_field()->discretization(), fluid_field()->discretization(),
          ScaTraField()->discretization());

    // the variables need to be set on other field
    SetScatraSolution();
    SetPoroSolution();

    // second restart needed due to two way coupling.
    ScaTraField()->read_restart(restart);
    poro_field()->read_restart(restart);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (poro_field()->HasSubmeshes())
      replace_dof_sets(structure_field()->discretization(), fluid_field()->discretization(),
          ScaTraField()->discretization());

    SetTimeStep(poro_field()->Time(), restart);

    if (matchinggrid_)
    {
      // Material pointers to other field were deleted during read_restart().
      // They need to be reset.
      PoroElast::UTILS::SetMaterialPointersMatchingGrid(
          poro_field()->structure_field()->discretization(), ScaTraField()->discretization());
      PoroElast::UTILS::SetMaterialPointersMatchingGrid(
          poro_field()->fluid_field()->discretization(), ScaTraField()->discretization());
    }
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::DoPoroStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  poro_field()->Solve();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::do_scatra_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  ScaTraField()->Solve();
}

/*----------------------------------------------------------------------*/
// prepare time step                                        vuong 08/13  |
/*----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::prepare_time_step(bool printheader)
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sync!
  increment_time_and_step();

  if (printheader) print_header();

  SetPoroSolution();
  ScaTraField()->prepare_time_step();
  // set structure-based scalar transport values
  SetScatraSolution();

  poro_field()->prepare_time_step();
  SetPoroSolution();
  // SetScatraSolution();
}


/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::prepare_output()
{
  constexpr bool force_prepare = false;
  poro_field()->prepare_output(force_prepare);
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::update()
{
  poro_field()->update();
  ScaTraField()->update();

  ScaTraField()->evaluate_error_compared_to_analytical_sol();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::output()
{
  poro_field()->Output();
  ScaTraField()->check_and_write_output_and_restart();
}


/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::Solve()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n*******************************************\n Poro-Scatra 2WC OUTER ITERATION "
                 "LOOP \n*******************************************\n";
  }

  while (stopnonliniter == false)
  {
    itnum++;

    // store scalar from first solution for convergence check (like in
    // elch_algorithm: use current values)
    scaincnp_->Update(1.0, *ScaTraField()->Phinp(), 0.0);
    structincnp_->Update(1.0, *poro_field()->structure_field()->Dispnp(), 0.0);
    fluidincnp_->Update(1.0, *poro_field()->fluid_field()->Velnp(), 0.0);

    // set structure-based scalar transport values
    SetScatraSolution();

    // solve structural system
    DoPoroStep();

    // set mesh displacement and velocity fields
    SetPoroSolution();

    // solve scalar transport equation
    do_scatra_step();
    // scatra_evaluate_solve_iter_update();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = convergence_check(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*
 | convergence check for both fields (scatra & poro) (copied from tsi)
 *----------------------------------------------------------------------*/
bool PoroElastScaTra::PoroScatraPart2WC::convergence_check(int itnum)
{
  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //     | scalar+1 |_2

  // variables to save different L2 - Norms
  // define L2-norm of incremental scalar and scalar
  double scaincnorm_L2(0.0);
  double scanorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);
  double fluidincnorm_L2(0.0);
  double fluidnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  scaincnp_->Update(1.0, *(ScaTraField()->Phinp()), -1.0);
  structincnp_->Update(1.0, *(poro_field()->structure_field()->Dispnp()), -1.0);
  fluidincnp_->Update(1.0, *(poro_field()->fluid_field()->Velnp()), -1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  poro_field()->structure_field()->Dispnp()->Norm2(&dispnorm_L2);
  fluidincnp_->Norm2(&fluidincnorm_L2);
  poro_field()->fluid_field()->Velnp()->Norm2(&fluidnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;

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
        "+--------------+------------------------+--------------------+--------------------+-------"
        "-------------+\n");
    printf(
        "|-  step/max  -|-  tol      [norm]     -|-  scalar-inc      -|-  disp-inc        -|-  "
        "fluid-inc       -|\n");
    printf(
        "|   %3d/%3d    |  %10.3E[L_2 ]      |  %10.3E        |  %10.3E        |  %10.3E        |",
        itnum, itmax_, ittol_, scaincnorm_L2 / scanorm_L2, dispincnorm_L2 / dispnorm_L2,
        fluidincnorm_L2 / fluidnorm_L2);
    printf("\n");
    printf(
        "+--------------+------------------------+--------------------+--------------------+-------"
        "-------------+\n");
  }

  // converged
  if ((scaincnorm_L2 / scanorm_L2 <= ittol_) and (dispincnorm_L2 / dispnorm_L2 <= ittol_) and
      (fluidincnorm_L2 / fluidnorm_L2 <= ittol_))
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
      ((scaincnorm_L2 / scanorm_L2 > ittol_) or (dispincnorm_L2 / dispnorm_L2 > ittol_) or
          (fluidincnorm_L2 / fluidnorm_L2 > ittol_)))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
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

FOUR_C_NAMESPACE_CLOSE
