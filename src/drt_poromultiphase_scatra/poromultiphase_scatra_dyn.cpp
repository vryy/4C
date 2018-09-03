/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_dyn.cpp

 \brief global access method for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_dyn.H"
#include "poromultiphase_scatra_base.H"
#include "poromultiphase_scatra_utils.H"

#include "../drt_inpar/inpar_poromultiphase_scatra.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | Main control routine                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
void poromultiphasescatra_dyn(int restart)
{
  // define the discretization names
  const std::string struct_disname = "structure";
  const std::string fluid_disname = "porofluid";
  const std::string scatra_disname = "scatra";

  // access the problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis(struct_disname)->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    POROMULTIPHASESCATRA::PrintLogo();
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << problem->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // Parameter reading
  // access poro multiphase scatra params list
  const Teuchos::ParameterList& poroscatraparams = problem->PoroMultiPhaseScatraDynamicParams();
  // access poro multiphase params list
  const Teuchos::ParameterList& poroparams = problem->PoroMultiPhaseDynamicParams();
  // access scatra params list
  const Teuchos::ParameterList& structparams = problem->StructuralDynamicParams();
  // access poro fluid dynamic params list
  const Teuchos::ParameterList& fluidparams = problem->PoroFluidMultiPhaseDynamicParams();
  // access scatra dynamic params list
  const Teuchos::ParameterList& scatraparams = problem->ScalarTransportDynamicParams();

  // do we perform coupling with 1D artery
  const bool artery_coupl = DRT::INPUT::IntegralValue<int>(poroscatraparams, "ARTERY_COUPLING");

  // initialize variables for dof set numbers
  int ndsporo_disp(-1);
  int ndsporo_vel(-1);
  int ndsporo_solidpressure(-1);
  int ndsporofluid_scatra(-1);

  // Setup discretizations and coupling. Assign the dof sets and return the numbers
  POROMULTIPHASESCATRA::UTILS::SetupDiscretizationsAndFieldCoupling(comm, struct_disname,
      fluid_disname, scatra_disname, ndsporo_disp, ndsporo_vel, ndsporo_solidpressure,
      ndsporofluid_scatra, artery_coupl);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields solscheme =
      DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields>(
          poroscatraparams, "COUPALGO");

  Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase> algo =
      POROMULTIPHASESCATRA::UTILS::CreatePoroMultiPhaseScatraAlgorithm(
          solscheme, poroscatraparams, comm);

  algo->Init(poroscatraparams, poroscatraparams, poroparams, structparams, fluidparams,
      scatraparams, struct_disname, fluid_disname, scatra_disname, true, ndsporo_disp, ndsporo_vel,
      ndsporo_solidpressure, ndsporofluid_scatra);

  // read the restart information, set vectors and variables
  if (restart) algo->ReadRestart(restart);

  // assign materials
  // note: to be done after potential restart, as in ReadRestart()
  //       the secondary material is destroyed
  POROMULTIPHASESCATRA::UTILS::AssignMaterialPointers(
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

  return;

}  // poromultiphase_dyn
