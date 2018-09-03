/*!----------------------------------------------------------------------
\file var_chemdiff_dyn.cpp
\brief Control routine for Chemical diffusion (variational form) module.

\level 2

<pre>
\maintainer Jorge De Anda Salazar
            DeAnda@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#include "../drt_ale/ale_utils_clonestrategy.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include <Epetra_Time.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "var_chemdiff_dyn.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"



/*----------------------------------------------------------------------*/
// entry point for Variational_Forms in DRT
/*----------------------------------------------------------------------*/
void var_chemdiff_dyn(int restart)
{
  // pointer to problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis("fluid")->Comm();

  // print (Var)iational form-Logo to screen
  if (comm.MyPID() == 0) print_var_logo();

  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order; this creates dof numbers with
  //       fluid dof < scatra/elch dof
  fluiddis->FillComplete();
  scatradis->FillComplete();

  // access the scalar transport parameter list
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();
  const INPAR::SCATRA::VelocityField veltype =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");

  // choose algorithm depending on velocity field type
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:      // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // spatial function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes() == 0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      // add proxy of velocity related degrees of freedom to scatra discretization
      Teuchos::RCP<DRT::DofSetInterface> dofsetaux = Teuchos::rcp(
          new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim() + 1, 0, 0, true));
      if (scatradis->AddDofSet(dofsetaux) != 1)
        dserror("Scatra discretization has illegal number of dofsets!");


      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror(
            "no linear solver defined for Variational chemical diffusion problem. Please set "
            "LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly =
          Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

      // now we can call Init() on the base algorithm
      // scatra time integrator is constructed and initialized inside
      scatraonly->Init(
          scatradyn, scatradyn, DRT::Problem::Instance()->SolverParams(linsolvernumber));

      // now me may redistribute or ghost the scatra discretization
      // finalize discretization
      scatradis->FillComplete(true, true, true);

      // all objects relying on the parallel distribution are
      // created and pointers are set.
      scatraonly->Setup();

      // read the restart information, set vectors and variables
      if (restart) (scatraonly->ScaTraField())->ReadRestart(restart);

      // set velocity field
      // note: The order ReadRestart() before SetVelocityField() is important here!!
      // for time-dependent velocity fields, SetVelocityField() is additionally called in each
      // PrepareTimeStep()-call
      (scatraonly->ScaTraField())->SetVelocityField(1);

      // enter time loop to solve problem with given convective velocity
      (scatraonly->ScaTraField())->TimeLoop();

      // perform the result test if required
      DRT::Problem::Instance()->AddFieldTest(scatraonly->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);

      break;
    }
    case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
    {
      dserror("Not implemented yet");
      break;
    }  // case 2
    default:
      dserror("Unknown velocity field type for transport of passive scalar: %d", veltype);
      break;
  }

  return;

}  // var_chemdiff_dyn()


/*----------------------------------------------------------------------*/
// print Variational Forms-Module logo
/*----------------------------------------------------------------------*/

void print_var_logo()
{
  // more at http://www.ascii-art.de under entry "wodka" (a minor but vital modification has been
  // implemented)
  std::cout << "            |HH|  " << std::endl;
  std::cout << "            )  (  " << std::endl;
  std::cout << "            (  )  " << std::endl;
  std::cout << "           _)  (_ "
            << "  __        __     ___         _______          " << std::endl;
  std::cout << "          |      |  "
            << "\\X\\ \\      / /  /X/ _ \\     |X|  __   \\    " << std::endl;
  std::cout << "          ||||||||  "
            << " \\X\\ \\    / /  /X/ /_\\ \\    |X| |__/  /    " << std::endl;
  std::cout << "   ____   | CHEM |  "
            << "  \\X\\ \\  / /  /X/  ___  \\   |X|  _   /      " << std::endl;
  std::cout << "  |~~~~~| |Diffus|  "
            << "   \\X\\ \\/ /  /X/  / \\X\\  \\  |X| |\\X\\ \\ " << std::endl;
  std::cout << "  `--,--' ||||||||  "
            << "    \\X\\__/  /X/__/   \\X\\__\\ |X|_| \\X\\_\\ " << std::endl;
  std::cout << "     |    |      |  " << std::endl;
  std::cout << "   __|__  |______|  " << std::endl;
}
