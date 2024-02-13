/*-----------------------------------------------------------*/
/*! \file

\brief Structural single step solver for explicit dynamics

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_nln_solver_singlestep.hpp"

#include "baci_solver_nonlin_nox_globaldata.hpp"
#include "baci_solver_nonlin_nox_group.hpp"
#include "baci_structure_new_integrator.hpp"
#include "baci_structure_new_nln_solver_utils.hpp"
#include "baci_structure_new_timint_base.hpp"

#include <NOX_Solver_Generic.H>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::SingleStep::Setup()
{
  CheckInit();

  // setup the nox parameter list for a full Newton solution method
  SetSingleStepParams();

  // Call the Setup() function of the base class
  // Note, that the issetup_ flag is also updated during this call.
  Nox::Setup();

  dsassert(IsSetup(), "issetup_ should be \"true\" at this point!");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::SingleStep::SetSingleStepParams()
{
  CheckInit();

  // get the nox parameter list and set the necessary parameters for a
  // full Newton solution procedure
  Teuchos::ParameterList& p = DataSDyn().GetNoxParams();
  SetSingleStepParams(p);

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /* This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the dat-file instead.
   */
  if (not IsXMLStatusTestFile(DataSDyn().GetNoxParams().sublist("Status Test")))
  {
    std::set<enum NOX::NLN::StatusTest::QuantityType> qtypes;
    CreateQuantityTypes(qtypes, DataSDyn());
    SetStatusTestParams(DataSDyn().GetNoxParams().sublist("Status Test"), DataSDyn(), qtypes);
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::SingleStep::SetSingleStepParams(Teuchos::ParameterList& p)
{
  // ---------------------------------------------------------------------------
  // Set-up the single step method
  // ---------------------------------------------------------------------------
  // Non-linear solver
  p.set("Nonlinear Solver", "Single Step");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = p.sublist("Printing");
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information",
      ::NOX::Utils::OuterIteration + ::NOX::Utils::OuterIterationStatusTest +
          ::NOX::Utils::InnerIteration + ::NOX::Utils::Details + ::NOX::Utils::Warning);

  // Set NOX solving parameters: "Update Jacobian" is essential to compute the Newton direction,
  // i.e. x=A^-1 b. For details behind how it's handled, see NOX_Solver_SingleStep.C
  Teuchos::ParameterList& sssParams = p.sublist("Single Step Solver");
  sssParams.set("Update Jacobian", true);
  sssParams.set("Print Norms", true);
  sssParams.set("Compute Relative Norm", true);

  // Sublist for linear solver for the single step
  Teuchos::ParameterList& lsParams = sssParams.sublist("Linear Solver");
  lsParams.set("Adaptive Control", false);
  lsParams.set("Number of Nonlinear Iterations", 1);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::SingleStep::ResetParams()
{
  SetSingleStepParams(nlnglobaldata_->GetNlnParameterList());
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ConvergenceStatus STR::NLN::SOLVER::SingleStep::Solve()
{
  CheckInitSetup();

  auto& nln_group = dynamic_cast<NOX::NLN::Group&>(Group());

  const auto& x_epetra = dynamic_cast<const ::NOX::Epetra::Vector&>(Group().getX());

  nln_group.setIsValidNewton(true);  // to circumvent the check in ::NOX::Solver::SingleStep
  nln_group.setIsValidRHS(false);    // force to compute the RHS
  nln_group.resetX();                // to initialize the solution vector to zero

  //// do one non-linear step using solve
  ::NOX::StatusTest::StatusType stepstatus = nlnsolver_->solve();

  // copy the solution group into the class variable
  Group() = nlnsolver_->getSolutionGroup();

  Integrator().SetState(x_epetra.getEpetraVector());

  return ConvertFinalStatus(stepstatus);
}

BACI_NAMESPACE_CLOSE
