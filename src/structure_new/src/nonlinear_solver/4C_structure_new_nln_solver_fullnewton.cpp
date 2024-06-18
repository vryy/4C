/*-----------------------------------------------------------*/
/*! \file
\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_nln_solver_fullnewton.hpp"

#include "4C_structure_new_nln_solver_utils.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Nln::SOLVER::FullNewton::FullNewton() : Nox()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Nln::SOLVER::FullNewton::setup()
{
  check_init();

  // setup the nox parameter list for a full Newton solution method
  set_full_newton_params();

  // Call the setup() function of the base class
  // Note, that the issetup_ flag is also updated during this call.
  Nox::setup();

  FOUR_C_ASSERT(is_setup(), "issetup_ should be \"true\" at this point!");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Nln::SOLVER::FullNewton::set_full_newton_params()
{
  check_init();

  // get the nox parameter list and set the necessary parameters for a
  // full Newton solution procedure
  Teuchos::ParameterList& p = data_sdyn().get_nox_params();

  // ---------------------------------------------------------------------------
  // Set-up the full Newton method
  // ---------------------------------------------------------------------------
  // Non-linear solver
  p.set("Nonlinear Solver", "Line Search Based");

  // Direction
  Teuchos::ParameterList& pdir = p.sublist("Direction", true);
  pdir.set("Method", "Newton");

  // Line Search
  Teuchos::ParameterList& plinesearch = p.sublist("Line Search", true);
  plinesearch.set("Method", "Full Step");

  // Line Search/Full Step
  Teuchos::ParameterList& pfullstep = plinesearch.sublist("Full Step", true);
  // check if the default value is set
  double fullstep = pfullstep.get<double>("Full Step");
  if (fullstep != 1.0)
  {
    std::string markerline;
    markerline.assign(40, '!');
    std::cout << markerline << std::endl
              << "WARNING: The Full Step length is " << fullstep << " (default=1.0)" << std::endl
              << markerline << std::endl;
  }

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /* This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the dat-file instead.
   */
  if (not is_xml_status_test_file(data_sdyn().get_nox_params().sublist("Status Test")))
  {
    std::set<enum NOX::Nln::StatusTest::QuantityType> qtypes;
    STR::Nln::SOLVER::create_quantity_types(qtypes, data_sdyn());
    {  // remove the unsupported quantity of status test
      qtypes.erase(NOX::Nln::StatusTest::quantity_eas);  // EAS is removed since it is an element
                                                         // quantity and not nodal dof
    }
    STR::Nln::SOLVER::set_status_test_params(
        data_sdyn().get_nox_params().sublist("Status Test"), data_sdyn(), qtypes);
  }
}

FOUR_C_NAMESPACE_CLOSE
