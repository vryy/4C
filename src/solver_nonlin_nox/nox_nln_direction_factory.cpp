/*-----------------------------------------------------------*/
/*!
\file nox_nln_direction_factory.cpp

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_direction_factory.H"

// supported directions
#include "nox_nln_direction_modified_newton.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Direction::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic> NOX::NLN::Direction::Factory::buildDirection(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
{
  Teuchos::RCP<NOX::Direction::Generic> direction;
  const std::string method = params.get<std::string>("User Defined Method");

  if (method == "Newton")
    direction = Teuchos::rcp(new Newton(gd, params));
  else if (method == "Modified Newton")
    direction = Teuchos::rcp(new ModifiedNewton(gd, params));
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::NLN::Direction::Facotry::buildDirection() - Invalid "
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
Teuchos::RCP<NOX::Direction::Generic> NOX::NLN::Direction::BuildDirection(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  Factory factory;
  return factory.buildDirection(gd, params);
}
