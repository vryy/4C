#include "4C_solver_nonlin_nox_meritfunction_factory.hpp"

#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::MeritFunction::Factory::Factory()
{
  // empty;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::MeritFunction::Generic> NOX::Nln::MeritFunction::Factory::build_merit_function(
    const NOX::Nln::GlobalData& noxNlnGlobalData) const
{
  Teuchos::RCP<::NOX::MeritFunction::Generic> mrtFctPtr;

  const Teuchos::ParameterList& solverOptionsList =
      noxNlnGlobalData.get_nln_parameter_list().sublist("Solver Options");

  const auto& mftype =
      solverOptionsList.get<NOX::Nln::MeritFunction::MeritFctName>("Merit Function");

  if (noxNlnGlobalData.is_constrained())
    mrtFctPtr = build_constrained_merit_function(mftype, noxNlnGlobalData);
  else
    mrtFctPtr = build_unconstrained_merit_function(mftype, noxNlnGlobalData);

  return mrtFctPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::MeritFunction::Generic>
NOX::Nln::MeritFunction::Factory::build_unconstrained_merit_function(
    const NOX::Nln::MeritFunction::MeritFctName& mftype,
    const NOX::Nln::GlobalData& noxNlnGlobalData) const
{
  Teuchos::RCP<::NOX::MeritFunction::Generic> mrtFctPtr = Teuchos::null;

  if (mftype == NOX::Nln::MeritFunction::MeritFctName::mrtfct_sum_of_squares)
  {
    // default NOX case, no pointer necessary
    mrtFctPtr = Teuchos::null;
  }
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::Nln::MeritFunction::buildUnconstraintedMeritFunction() - \n"
           "The \"Solver Options > Merit Function\" parameter \""
        << mftype
        << "\" is not a valid UNCONSTRAINTED merit function name.  Please "
           "fix your parameter list!";
    FOUR_C_THROW(msg.str().c_str());
  }

  return mrtFctPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::MeritFunction::Generic>
NOX::Nln::MeritFunction::Factory::build_constrained_merit_function(
    const NOX::Nln::MeritFunction::MeritFctName& mftype,
    const NOX::Nln::GlobalData& noxNlnGlobalData) const
{
  Teuchos::RCP<::NOX::MeritFunction::Generic> mrtFctPtr = Teuchos::null;

  if (mftype == NOX::Nln::MeritFunction::MeritFctName::mrtfct_sum_of_squares)
  {
    // default NOX case, no pointer necessary
    mrtFctPtr = Teuchos::null;
  }
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::Nln::MeritFunction::buildConstrainedMeritFunction() - \n"
           "The \"Solver Options > Merit Function\" parameter \""
        << mftype
        << "\" is not a valid CONSTRAINT merit function name.  Please "
           "fix your parameter list!";
    FOUR_C_THROW(msg.str().c_str());
  }

  return mrtFctPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::MeritFunction::Generic> NOX::Nln::MeritFunction::build_merit_function(
    const NOX::Nln::GlobalData& noxNlnGlobalData)
{
  const Factory factory;
  return factory.build_merit_function(noxNlnGlobalData);
}

FOUR_C_NAMESPACE_CLOSE
