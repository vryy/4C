/*-----------------------------------------------------------*/
/*! \file
\brief NOX's Newton with full step


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_nln_solver_fullnewton.H"

#include "baci_structure_new_nln_solver_utils.H"
#include "baci_structure_new_timint_base.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::FullNewton::FullNewton() : Nox()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::Setup()
{
  CheckInit();

  // setup the nox parameter list for a full Newton solution method
  SetFullNewtonParams();

  // Call the Setup() function of the base class
  // Note, that the issetup_ flag is also updated during this call.
  Nox::Setup();

  dsassert(IsSetup(), "issetup_ should be \"true\" at this point!");

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetFullNewtonParams()
{
  CheckInit();

  // get the nox parameter list and set the necessary parameters for a
  // full Newton solution procedure
  Teuchos::ParameterList& p = DataSDyn().GetNoxParams();

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
  if (not IsXMLStatusTestFile(DataSDyn().GetNoxParams().sublist("Status Test")))
  {
    std::set<enum NOX::NLN::StatusTest::QuantityType> qtypes;
    STR::NLN::SOLVER::CreateQuantityTypes(qtypes, DataSDyn());
    {  // remove the unsupported quantity of status test
      qtypes.erase(NOX::NLN::StatusTest::quantity_eas);  // EAS is removed since it is an element
                                                         // quantity and not nodal dof
    }
    STR::NLN::SOLVER::SetStatusTestParams(
        DataSDyn().GetNoxParams().sublist("Status Test"), DataSDyn(), qtypes);
  }

  return;
}
