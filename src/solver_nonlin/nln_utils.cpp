/*----------------------------------------------------------------------*/
/*!
\file nln_utils.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <string>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

// baci
#include "nln_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Create parameter list from xml-file */
Teuchos::RCP<Teuchos::ParameterList>
NLNSOL::UTILS::CreateParamListFromXML()
{
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());

  std::string xmlfilename = DRT::Problem::Instance()->NonlinearSolverParams().get<std::string>("XML_FILE");

  if(xmlfilename.length())
  {
    Teuchos::updateParametersFromXmlFile(xmlfilename, inoutArg(*params));
  }
  else
  {
    dserror("The file name '%s' is not a valid XML file name.", xmlfilename.c_str());
  }

  return params;
}
