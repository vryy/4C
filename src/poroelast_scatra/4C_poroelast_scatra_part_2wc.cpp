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
      scaincnp_(Teuchos::rcp(new Epetra_Vector(*(sca_tra_field()->phinp())))),
      structincnp_(Teuchos::rcp(new Epetra_Vector(*(poro_field()->structure_field()()->dispnp())))),
      fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(poro_field()->fluid_field()()->velnp()))))
{
  if (comm.MyPID() == 0) std::cout << "\n Create PoroScatraPart2WC algorithm ... \n" << std::endl;

  const Teuchos::ParameterList& params = Global::Problem::instance()->poro_scatra_control_params();
  // Get the parameters for the convergence_check
  itmax_ = params.get<int>("ITEMAX");            // default: =10
  ittol_ = params.get<double>("TOLINC_GLOBAL");  // default: =1e-6
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::timeloop()
{
  // initial_calculations();

  while (not_finished())
  {
    prepare_time_step();

    solve();

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
    set_scatra_solution();
    set_poro_solution();

    poro_field()->read_restart(restart);
    sca_tra_field()->read_restart(restart);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (poro_field()->has_submeshes())
      replace_dof_sets(structure_field()->discretization(), fluid_field()->discretization(),
          sca_tra_field()->discretization());

    // the variables need to be set on other field
    set_scatra_solution();
    set_poro_solution();

    // second restart needed due to two way coupling.
    sca_tra_field()->read_restart(restart);
    poro_field()->read_restart(restart);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (poro_field()->has_submeshes())
      replace_dof_sets(structure_field()->discretization(), fluid_field()->discretization(),
          sca_tra_field()->discretization());

    set_time_step(poro_field()->time(), restart);

    if (matchinggrid_)
    {
      // Material pointers to other field were deleted during read_restart().
      // They need to be reset.
      PoroElast::UTILS::SetMaterialPointersMatchingGrid(
          poro_field()->structure_field()->discretization(), sca_tra_field()->discretization());
      PoroElast::UTILS::SetMaterialPointersMatchingGrid(
          poro_field()->fluid_field()->discretization(), sca_tra_field()->discretization());
    }
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::do_poro_step()
{
  if (get_comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  poro_field()->solve();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::do_scatra_step()
{
  if (get_comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  sca_tra_field()->solve();
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

  set_poro_solution();
  sca_tra_field()->prepare_time_step();
  // set structure-based scalar transport values
  set_scatra_solution();

  poro_field()->prepare_time_step();
  set_poro_solution();
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
  sca_tra_field()->update();

  sca_tra_field()->evaluate_error_compared_to_analytical_sol();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::output()
{
  poro_field()->output();
  sca_tra_field()->check_and_write_output_and_restart();
}


/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart2WC::solve()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (get_comm().MyPID() == 0)
  {
    std::cout << "\n*******************************************\n Poro-Scatra 2WC OUTER ITERATION "
                 "LOOP \n*******************************************\n";
  }

  while (stopnonliniter == false)
  {
    itnum++;

    // store scalar from first solution for convergence check (like in
    // elch_algorithm: use current values)
    scaincnp_->Update(1.0, *sca_tra_field()->phinp(), 0.0);
    structincnp_->Update(1.0, *poro_field()->structure_field()->dispnp(), 0.0);
    fluidincnp_->Update(1.0, *poro_field()->fluid_field()->velnp(), 0.0);

    // set structure-based scalar transport values
    set_scatra_solution();

    // solve structural system
    do_poro_step();

    // set mesh displacement and velocity fields
    set_poro_solution();

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
  scaincnp_->Update(1.0, *(sca_tra_field()->phinp()), -1.0);
  structincnp_->Update(1.0, *(poro_field()->structure_field()->dispnp()), -1.0);
  fluidincnp_->Update(1.0, *(poro_field()->fluid_field()->velnp()), -1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  sca_tra_field()->phinp()->Norm2(&scanorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  poro_field()->structure_field()->dispnp()->Norm2(&dispnorm_L2);
  fluidincnp_->Norm2(&fluidincnorm_L2);
  poro_field()->fluid_field()->velnp()->Norm2(&fluidnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (get_comm().MyPID() == 0)
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
    if (get_comm().MyPID() == 0)
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
    if ((get_comm().MyPID() == 0))
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
