/*----------------------------------------------------------------------*/
/*! \file
 \brief global access method for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_scatra_dyn.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_poromultiphase_scatra_base.hpp"
#include "4C_poromultiphase_scatra_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void poromultiphasescatra_dyn(int restart)
{
  // define the discretization names
  const std::string struct_disname = "structure";
  const std::string fluid_disname = "porofluid";
  const std::string scatra_disname = "scatra";

  // access the problem
  Global::Problem* problem = Global::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis(struct_disname)->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    PoroMultiPhaseScaTra::PrintLogo();
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << problem->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // Parameter reading
  // access poro multiphase scatra params list
  const Teuchos::ParameterList& poroscatraparams =
      problem->poro_multi_phase_scatra_dynamic_params();
  // access poro multiphase params list
  const Teuchos::ParameterList& poroparams = problem->poro_multi_phase_dynamic_params();
  // access scatra params list
  const Teuchos::ParameterList& structparams = problem->structural_dynamic_params();
  // access poro fluid dynamic params list
  const Teuchos::ParameterList& fluidparams = problem->poro_fluid_multi_phase_dynamic_params();
  // access scatra dynamic params list
  const Teuchos::ParameterList& scatraparams = problem->scalar_transport_dynamic_params();

  // do we perform coupling with 1D artery
  const bool artery_coupl = Core::UTILS::IntegralValue<int>(poroscatraparams, "ARTERY_COUPLING");

  // initialize variables for dof set numbers
  int ndsporo_disp(-1);
  int ndsporo_vel(-1);
  int ndsporo_solidpressure(-1);
  int ndsporofluid_scatra(-1);

  // Setup discretizations and coupling. Assign the dof sets and return the numbers
  std::map<int, std::set<int>> nearbyelepairs =
      PoroMultiPhaseScaTra::UTILS::SetupDiscretizationsAndFieldCoupling(comm, struct_disname,
          fluid_disname, scatra_disname, ndsporo_disp, ndsporo_vel, ndsporo_solidpressure,
          ndsporofluid_scatra, artery_coupl);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields solscheme =
      Core::UTILS::IntegralValue<Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields>(
          poroscatraparams, "COUPALGO");

  Teuchos::RCP<PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase> algo =
      PoroMultiPhaseScaTra::UTILS::CreatePoroMultiPhaseScatraAlgorithm(
          solscheme, poroscatraparams, comm);

  algo->Init(poroscatraparams, poroscatraparams, poroparams, structparams, fluidparams,
      scatraparams, struct_disname, fluid_disname, scatra_disname, true, ndsporo_disp, ndsporo_vel,
      ndsporo_solidpressure, ndsporofluid_scatra, &nearbyelepairs);

  // read the restart information, set vectors and variables
  if (restart) algo->read_restart(restart);

  // assign materials
  // note: to be done after potential restart, as in read_restart()
  //       the secondary material is destroyed
  PoroMultiPhaseScaTra::UTILS::assign_material_pointers(
      struct_disname, fluid_disname, scatra_disname, artery_coupl);

  // Setup Solver (necessary if poro-structure coupling solved monolithically)
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


}  // poromultiphase_dyn

FOUR_C_NAMESPACE_CLOSE
