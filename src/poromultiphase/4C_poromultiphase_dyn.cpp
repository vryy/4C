/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for poroelasticity multiphase flow

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_poromultiphase_dyn.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_poromultiphase.hpp"
#include "4C_poromultiphase_base.hpp"
#include "4C_poromultiphase_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | Main control routine                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
void poromultiphase_dyn(int restart)
{
  // define the discretization names
  const std::string struct_disname = "structure";
  const std::string fluid_disname = "porofluid";

  // access the problem
  Global::Problem* problem = Global::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis(struct_disname)->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    POROMULTIPHASE::PrintLogo();
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << problem->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // initialize variables for dof set numbers
  int nds_disp(-1);
  int nds_vel(-1);
  int nds_solidpressure(-1);

  // Setup discretizations and coupling. Assign the dof sets and return the numbers
  std::map<int, std::set<int>> nearbyelepairs =
      POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
          comm, struct_disname, fluid_disname, nds_disp, nds_vel, nds_solidpressure);

  // Parameter reading
  const Teuchos::ParameterList& poroparams = problem->poro_multi_phase_dynamic_params();
  // access scatra params list
  const Teuchos::ParameterList& structdyn = problem->structural_dynamic_params();
  // access poro fluid dynamic params list
  const Teuchos::ParameterList& fluiddyn = problem->poro_fluid_multi_phase_dynamic_params();

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  Inpar::POROMULTIPHASE::SolutionSchemeOverFields solscheme =
      Core::UTILS::IntegralValue<Inpar::POROMULTIPHASE::SolutionSchemeOverFields>(
          poroparams, "COUPALGO");

  Teuchos::RCP<Adapter::PoroMultiPhase> algo =
      POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(solscheme, poroparams, comm);

  // initialize
  algo->Init(poroparams, poroparams, structdyn, fluiddyn, struct_disname, fluid_disname, true,
      nds_disp, nds_vel, nds_solidpressure,
      -1,  // no scalar field
      &nearbyelepairs);

  // read the restart information, set vectors and variables
  if (restart) algo->read_restart(restart);

  // assign poro material for evaluation of porosity
  // note: to be done after potential restart, as in read_restart()
  //       the secondary material is destroyed
  POROMULTIPHASE::UTILS::assign_material_pointers(struct_disname, fluid_disname);

  // Setup the solver (only for the monolithic problem)
  algo->SetupSolver();

  // Run of the actual problem.

  // Some setup needed for the subproblems.
  algo->SetupSystem();

  // Solve the whole problem
  algo->Timeloop();

  // Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test if required
  algo->CreateFieldTest();
  problem->TestAll(comm);

  return;

}  // poromultiphase_dyn

FOUR_C_NAMESPACE_CLOSE
