/*-----------------------------------------------------------*/
/*!
\file nox_nln_meritfunction_factory.cpp

\brief Factory to create a merit function evaluation object.

\maintainer Michael Hiermeier

\date Jun 12, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_meritfunction_factory.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <NOX_Utils.H>

// All the different merit functions
#include "nox_nln_meritfunction_lagrangian.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::MeritFunction::Factory::Factory()
{
  // empty;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::MeritFunction::Generic> NOX::NLN::MeritFunction::Factory::
  BuildMeritFunction(Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData)
{
  Teuchos::RCP<NOX::MeritFunction::Generic> mrtFctPtr;

  if (noxNlnGlobalData->GetIsConstrained())
    mrtFctPtr = BuildConstrainedMeritFunction(noxNlnGlobalData);
  else
    mrtFctPtr = BuildUnconstrainedMeritFunction(noxNlnGlobalData);

  return mrtFctPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::MeritFunction::Generic> NOX::NLN::MeritFunction::Factory::
  BuildUnconstrainedMeritFunction(Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData)
{
  Teuchos::RCP<NOX::MeritFunction::Generic> mrtFctPtr = Teuchos::null;

  const Teuchos::ParameterList& solverOptionsList =
      noxNlnGlobalData->GetNlnParameterList().sublist("Solver Options");

  std::string mftype = solverOptionsList.get<std::string>("Merit Function");

  if (mftype == "Sum of Squares")
  {
    // default NOX case, no pointer necessary
    mrtFctPtr = Teuchos::null;
  }
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::NLN::MeritFunction::buildUnconstraintedMeritFunction() - \n"
        "The \"Solver Options > Merit Function\" parameter \""
        << mftype << "\" is not a valid UNCONSTRAINTED merit function name.  Please "
            "fix your parameter list!";
    dserror(msg.str().c_str());
  }

  return mrtFctPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::MeritFunction::Generic> NOX::NLN::MeritFunction::Factory::
  BuildConstrainedMeritFunction(Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData)
{
  Teuchos::RCP<NOX::MeritFunction::Generic> mrtFctPtr = Teuchos::null;

  const Teuchos::ParameterList& solverOptionsList =
      noxNlnGlobalData->GetNlnParameterList().sublist("Solver Options");

  std::string mftype = solverOptionsList.get<std::string>("Merit Function");
  if (mftype=="Lagrangian")
  {
    mrtFctPtr = Teuchos::rcp(new Lagrangian(noxNlnGlobalData->GetNoxUtilsPtr()));
  }
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::NLN::MeritFunction::buildConstrainedMeritFunction() - \n"
        "The \"Solver Options > Merit Function\" parameter \""
        << mftype << "\" is not a valid CONSTRAINT merit function name.  Please "
            "fix your parameter list!";
    dserror(msg.str().c_str());
  }

  return mrtFctPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::MeritFunction::Generic> NOX::NLN::MeritFunction::
BuildMeritFunction(Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData)
{
  Factory factory;
  return factory.BuildMeritFunction(noxNlnGlobalData);
}
