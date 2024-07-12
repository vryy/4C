/*----------------------------------------------------------------------*/
/*! \file

 \brief partitioned one way coupled poroelasticity scalar transport interaction algorithms

\level 2

*/
/*---------------------------------------------------------------------*/


#include "4C_poroelast_scatra_part_1wc.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WC::do_poro_step()
{
  if (get_comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }
  // 1)  solve the step problem. Methods obtained from poroelast->TimeLoop(sdynparams); -->
  // sdynparams
  //      CUIDADO, aqui vuelve a avanzar el paso de tiempo. Hay que corregir eso.
  // 2)  Newton-Raphson iteration
  poro_field()->solve();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WC::do_scatra_step()
{
  if (get_comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }
  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_field()->solve();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WC::prepare_output()
{
  constexpr bool force_prepare = false;
  poro_field()->prepare_output(force_prepare);
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WC::update()
{
  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  poro_field()->update();
  scatra_field()->update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  scatra_field()->evaluate_error_compared_to_analytical_sol();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WC::output()
{
  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  poro_field()->output();
  scatra_field()->check_and_write_output_and_restart();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
PoroElastScaTra::PoroScatraPart1WCPoroToScatra::PoroScatraPart1WCPoroToScatra(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraPart1WC(comm, timeparams)
{
  if (comm.MyPID() == 0)
    std::cout << "\n Create PoroScatraPart1WCPoroToScatra algorithm ... \n" << std::endl;
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCPoroToScatra::timeloop()
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

/*----------------------------------------------------------------------*/
// prepare time step                                  rauch/vuong 04/15  |
/*----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCPoroToScatra::prepare_time_step(bool printheader)
{
  increment_time_and_step();
  if (printheader) print_header();

  poro_field()->prepare_time_step();
  set_poro_solution();
  scatra_field()->prepare_time_step();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCPoroToScatra::solve()
{
  do_poro_step();  // It has its own time and timestep variables, and it increments them by itself.
  set_poro_solution();
  do_scatra_step();  // It has its own time and timestep variables, and it increments them by
                     // itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCPoroToScatra::read_restart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    poro_field()->read_restart(restart);
    set_poro_solution();
    scatra_field()->read_restart(restart);

    set_time_step(poro_field()->time(), restart);

    // Material pointers to other field were deleted during read_restart().
    // They need to be reset.
    PoroElast::UTILS::SetMaterialPointersMatchingGrid(
        poro_field()->structure_field()->discretization(), scatra_field()->discretization());
    PoroElast::UTILS::SetMaterialPointersMatchingGrid(
        poro_field()->fluid_field()->discretization(), scatra_field()->discretization());
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
PoroElastScaTra::PoroScatraPart1WCScatraToPoro::PoroScatraPart1WCScatraToPoro(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraPart1WC(comm, timeparams)
{
  if (comm.MyPID() == 0)
    std::cout << "\n Create PoroScatraPart1WCScatraToPoro algorithm ... \n" << std::endl;

  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset =
      scatra_field()->discretization()->get_dof_set_proxy();

  // check if structure field has 2 discretizations, so that coupling is possible
  if (poro_field()->structure_field()->discretization()->add_dof_set(scatradofset) != 1)
    FOUR_C_THROW("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCScatraToPoro::timeloop()
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

/*----------------------------------------------------------------------*/
// prepare time step                                  rauch/vuong 04/15  |
/*----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCScatraToPoro::prepare_time_step(bool printheader)
{
  increment_time_and_step();
  if (printheader) print_header();

  scatra_field()->prepare_time_step();
  set_scatra_solution();
  poro_field()->prepare_time_step();
}


/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCScatraToPoro::solve()
{
  do_scatra_step();  // It has its own time and timestep variables, and it increments them by
                     // itself.
  set_scatra_solution();
  do_poro_step();  // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraPart1WCScatraToPoro::read_restart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    scatra_field()->read_restart(restart);
    set_scatra_solution();
    poro_field()->read_restart(restart);

    set_time_step(poro_field()->time(), restart);

    // Material pointers to other field were deleted during read_restart().
    // They need to be reset.
    PoroElast::UTILS::SetMaterialPointersMatchingGrid(
        poro_field()->structure_field()->discretization(), scatra_field()->discretization());
    PoroElast::UTILS::SetMaterialPointersMatchingGrid(
        poro_field()->fluid_field()->discretization(), scatra_field()->discretization());
  }
}

FOUR_C_NAMESPACE_CLOSE
