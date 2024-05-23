/*--------------------------------------------------------------------------*/
/*! \file

\brief entry point for the solution of Lubrication problems

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "4C_lubrication_dyn.hpp"

#include "4C_adapter_lubrication.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_lib_utils_createdis.hpp"
#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | Main control routine for Lubrication problems            wirtz 11/15 |
 *----------------------------------------------------------------------*/
void lubrication_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("lubrication")->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << GLOBAL::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& lubricationdyn =
      GLOBAL::Problem::Instance()->lubrication_dynamic_params();

  // access the lubrication discretization
  Teuchos::RCP<DRT::Discretization> lubricationdis =
      GLOBAL::Problem::Instance()->GetDis("lubrication");

  lubricationdis->FillComplete();

  // we directly use the elements from the Lubrication elements section
  if (lubricationdis->NumGlobalNodes() == 0)
    FOUR_C_THROW("No elements in the ---LUBRICATION ELEMENTS section");

  // add proxy of velocity related degrees of freedom to lubrication discretization
  Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(
          GLOBAL::Problem::Instance()->NDim(), 0, 0, true));
  if (lubricationdis->AddDofSet(dofsetaux) != 1)
    FOUR_C_THROW("lub discretization has illegal number of dofsets!");

  // finalize discretization
  lubricationdis->FillComplete(true, false, false);

  // get linear solver id from LUBRICATION DYNAMIC
  const int linsolvernumber = lubricationdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for LUBRICATION problem. Please set LINEAR_SOLVER in LUBRICATION "
        "DYNAMIC to a valid number!");

  // create instance of Lubrication basis algorithm
  Teuchos::RCP<ADAPTER::LubricationBaseAlgorithm> lubricationonly =
      Teuchos::rcp(new ADAPTER::LubricationBaseAlgorithm());

  // setup Lubrication basis algorithm
  lubricationonly->Setup(
      lubricationdyn, lubricationdyn, GLOBAL::Problem::Instance()->SolverParams(linsolvernumber));

  // read the restart information, set vectors and variables
  if (restart) lubricationonly->LubricationField()->ReadRestart(restart);

  // enter time loop to solve problem
  (lubricationonly->LubricationField())->TimeLoop();

  // perform the result test if required
  GLOBAL::Problem::Instance()->AddFieldTest(lubricationonly->create_lubrication_field_test());
  GLOBAL::Problem::Instance()->TestAll(comm);

  return;

}  // end of lubrication_dyn()

FOUR_C_NAMESPACE_CLOSE
