/*-----------------------------------------------------------*/
/*!
\file nox_nln_direction_factory.cpp

\maintainer Michael Hiermeier

\date Aug 6, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_direction_factory.H"

// supported directions
#include "nox_nln_direction_newton.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Direction::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic> NOX::NLN::Direction::Factory::BuildDirection(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  Teuchos::RCP<NOX::Direction::Generic> direction;

  std::string method = params.get("Method", "Newton");

  if (method == "Newton")
    direction = Teuchos::rcp(new Newton(gd, params));
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::NLN::Direction::Facotry::buildDirection() - Invalid choice for \"Method\" "
           "in \"Direction\" sublist!\n";
    msg << "The \"Direction\"-\"Method\" = " << method << " is not (yet) supported!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  return direction;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic> NOX::NLN::Direction::BuildDirection(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  Factory factory;
  return factory.BuildDirection(gd, params);
}
