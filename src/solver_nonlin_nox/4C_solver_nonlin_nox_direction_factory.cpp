/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_direction_factory.hpp"

#include "4C_solver_nonlin_nox_direction_newton.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Direction::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Direction::Generic> NOX::Nln::Direction::Factory::buildDirection(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
{
  Teuchos::RCP<::NOX::Direction::Generic> direction;
  const std::string method = params.get<std::string>("User Defined Method");

  if (method == "Newton")
    direction = Teuchos::make_rcp<Newton>(gd, params);
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::Nln::Direction::Facotry::buildDirection() - Invalid "
           "choice for \"Method\" or \"User Defined Method\" in \"Direction\" sublist!\n";
    msg << "The \"Direction\"-\"(User Defined) Method\" = " << method
        << " is "
           "not (yet) supported!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  return direction;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Direction::Generic> NOX::Nln::Direction::build_direction(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  Factory factory;
  return factory.buildDirection(gd, params);
}

FOUR_C_NAMESPACE_CLOSE
