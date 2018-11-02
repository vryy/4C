/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_dyn.cpp

 \brief entry point (global control routine) for poroelasticity multiphase flow

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/


#include "poromultiphase_dyn.H"
#include "poromultiphase_base.H"

#include "poromultiphase_utils.H"

#include "../drt_inpar/inpar_poromultiphase.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | Main control routine                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
void poromultiphase_dyn(int restart)
{
  // define the discretization names
  const std::string struct_disname = "structure";
  const std::string fluid_disname = "porofluid";

  // access the problem
  DRT::Problem* problem = DRT::Problem::Instance();

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
  POROMULTIPHASE::UTILS::SetupDiscretizationsAndFieldCoupling(
      comm, struct_disname, fluid_disname, nds_disp, nds_vel, nds_solidpressure);

  // Parameter reading
  const Teuchos::ParameterList& poroparams = problem->PoroMultiPhaseDynamicParams();
  // access scatra params list
  const Teuchos::ParameterList& structdyn = problem->StructuralDynamicParams();
  // access poro fluid dynamic params list
  const Teuchos::ParameterList& fluiddyn = problem->PoroFluidMultiPhaseDynamicParams();

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  INPAR::POROMULTIPHASE::SolutionSchemeOverFields solscheme =
      DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::SolutionSchemeOverFields>(
          poroparams, "COUPALGO");

  Teuchos::RCP<ADAPTER::PoroMultiPhase> algo =
      POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(solscheme, poroparams, comm);

  // initialize
  algo->Init(poroparams, poroparams, structdyn, fluiddyn, struct_disname, fluid_disname, true,
      nds_disp, nds_vel, nds_solidpressure,
      -1  // no scalar field
  );

  // read the restart information, set vectors and variables
  if (restart) algo->ReadRestart(restart);

  // assign poro material for evaluation of porosity
  // note: to be done after potential restart, as in ReadRestart()
  //       the secondary material is destroyed
  POROMULTIPHASE::UTILS::AssignMaterialPointers(struct_disname, fluid_disname);

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
