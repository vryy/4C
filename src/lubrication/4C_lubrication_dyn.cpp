/*--------------------------------------------------------------------------*/
/*! \file

\brief entry point for the solution of Lubrication problems

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "4C_lubrication_dyn.hpp"

#include "4C_adapter_lubrication.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | Main control routine for Lubrication problems            wirtz 11/15 |
 *----------------------------------------------------------------------*/
void lubrication_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("lubrication")->get_comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << Global::Problem::instance()->problem_name()
              << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& lubricationdyn =
      Global::Problem::instance()->lubrication_dynamic_params();

  // access the lubrication discretization
  Teuchos::RCP<Core::FE::Discretization> lubricationdis =
      Global::Problem::instance()->get_dis("lubrication");

  lubricationdis->fill_complete();

  // we directly use the elements from the Lubrication elements section
  if (lubricationdis->num_global_nodes() == 0)
    FOUR_C_THROW("No elements in the ---LUBRICATION ELEMENTS section");

  // add proxy of velocity related degrees of freedom to lubrication discretization
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
          Global::Problem::instance()->n_dim(), 0, 0, true));
  if (lubricationdis->add_dof_set(dofsetaux) != 1)
    FOUR_C_THROW("lub discretization has illegal number of dofsets!");

  // finalize discretization
  lubricationdis->fill_complete(true, false, false);

  // get linear solver id from LUBRICATION DYNAMIC
  const int linsolvernumber = lubricationdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for LUBRICATION problem. Please set LINEAR_SOLVER in LUBRICATION "
        "DYNAMIC to a valid number!");

  // create instance of Lubrication basis algorithm
  Teuchos::RCP<Adapter::LubricationBaseAlgorithm> lubricationonly =
      Teuchos::rcp(new Adapter::LubricationBaseAlgorithm());

  // setup Lubrication basis algorithm
  lubricationonly->setup(
      lubricationdyn, lubricationdyn, Global::Problem::instance()->solver_params(linsolvernumber));

  // read the restart information, set vectors and variables
  if (restart) lubricationonly->lubrication_field()->read_restart(restart);

  // enter time loop to solve problem
  (lubricationonly->lubrication_field())->time_loop();

  // perform the result test if required
  Global::Problem::instance()->add_field_test(lubricationonly->create_lubrication_field_test());
  Global::Problem::instance()->test_all(comm);

  return;

}  // end of lubrication_dyn()

FOUR_C_NAMESPACE_CLOSE
