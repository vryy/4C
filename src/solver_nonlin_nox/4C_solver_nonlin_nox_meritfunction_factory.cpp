/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create a merit function evaluation object.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_meritfunction_factory.hpp"

#include "4C_solver_nonlin_nox_meritfunction_lagrangian.hpp"

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

  const std::string& mftype = solverOptionsList.get<std::string>("Merit Function");

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
    const std::string& mftype, const NOX::Nln::GlobalData& noxNlnGlobalData) const
{
  Teuchos::RCP<::NOX::MeritFunction::Generic> mrtFctPtr = Teuchos::null;

  if (mftype == "Sum of Squares")
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
    const std::string& mftype, const NOX::Nln::GlobalData& noxNlnGlobalData) const
{
  Teuchos::RCP<::NOX::MeritFunction::Generic> mrtFctPtr = Teuchos::null;

  if (mftype == "Sum of Squares")
  {
    // default NOX case, no pointer necessary
    mrtFctPtr = Teuchos::null;
  }
  else if (mftype.substr(0, 10) == "Lagrangian")
  {
    mrtFctPtr = Teuchos::rcp(new Lagrangian(mftype, noxNlnGlobalData.get_NOX_utils_ptr()));
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
Teuchos::RCP<::NOX::MeritFunction::Generic> NOX::Nln::MeritFunction::BuildMeritFunction(
    const NOX::Nln::GlobalData& noxNlnGlobalData)
{
  const Factory factory;
  return factory.build_merit_function(noxNlnGlobalData);
}

FOUR_C_NAMESPACE_CLOSE
