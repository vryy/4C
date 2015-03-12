/*----------------------------------------------------------------------------*/
/*!
\file nln_utils.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <string>

// Teuchos
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

// baci
#include "nln_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::UTILS::StagnationDetection::StagnationDetection()
: isinit_(false),
  stagiter_(1),
  stagitermax_(0),
  stagthreshold_(1.0),
  normprev_(1.0e+12)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::StagnationDetection::Init(const double norminitial)
{
  // Init() may be called only once
  if (IsInit())
    dserror("Init() has already been called. Don't call it again!");

  // initialize
  stagiter_ = 1;
  stagitermax_ = 3;
  stagthreshold_ = 0.95;
  normprev_ = norminitial;

  // Init() has been called
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const bool NLNSOL::UTILS::StagnationDetection::Check(const double norm)
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // ratio of residual norms of two subsequent iterations
  double ratio = norm / normprev_;

  // update
  normprev_ = norm;

  // ---------------------------------------------------------------------------
  // decide whether this is considered as stagnation
  // ---------------------------------------------------------------------------
  if (ratio > stagthreshold_)
    ++stagiter_;
  else
    stagiter_ = 0;

  std::cout << "StagnationCheck: stagiter = " << stagiter_ << " / " << stagitermax_
      << ", ratio = " << ratio << std::endl;

  return Status();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const bool NLNSOL::UTILS::StagnationDetection::Status() const
{
  bool stagnation = false;
  if (stagiter_ >= stagitermax_)
    stagnation = true;

  return stagnation;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> NLNSOL::UTILS::CreateParamListFromXML()
{
  std::string filename =
      DRT::Problem::Instance()->NonlinearSolverParams().get<std::string>(
          "XML_FILE");

  // check for reasonable xml file
  if (filename == "none")
    dserror("Seems like you forgot to set the XML file for configuration of "
        "the nonlinear solver.");

  return CreateParamListFromXML(filename);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> NLNSOL::UTILS::CreateParamListFromXML(
    const std::string filename)
{
  // create a new parameter list to be filled
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList());

  // fill the parameter list with the content of the given XML-file
  CreateParamListFromXML(filename, *params);

  return params;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::CreateParamListFromXML(const std::string filename,
    Teuchos::ParameterList& params)
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::UTILS::CreateParamListFromXML");
  Teuchos::TimeMonitor monitor(*time);

  // check for reasonable filename of xml-file to be read
  if (filename.length() && filename.rfind(".xml"))
  {
    // fill parameter list with content of xml-file
    params = *(Teuchos::getParametersFromXmlFile(filename));
  }
  else
    dserror("The file name '%s' is not a valid XML file name.",
        filename.c_str());

  // loop over all entries (parameters and sublists) of parameter list 'params'
  for (Teuchos::ParameterList::ConstIterator it = params.begin(); it != params.end(); ++it)
  {
    // Is this entry of type "std::string"?
    if (it->second.isType<std::string>())
    {
      // get the string
      std::string str;
      str = it->second.getValue(&str);

      /* Is this string a filename to an xml-file for a sublist?
       * Note: We introduced a naming convention to identify those cases (see
       * header file).
       */
      if (not str.find("sublistfile:") && str.rfind(".xml") != std::string::npos)
      {
        // extract filename from string by removing the flag 'sublistfile: '
        std::string sublistfilename;
        sublistfilename.assign(str.begin()+13, str.end());

        /* remove the place holder Teuchos::ParameterEntry and replace it by a
         * sublist with the same name.
         */
        const std::string sublistname = it->first;
        params.remove(it->first);
        Teuchos::ParameterList& subparams = params.sublist(sublistname);

        // check for reasonable filename of xml-file to be read
        if (filename.length() && filename.rfind(".xml"))
        {
          // fill sublist with content of xml-file
          CreateParamListFromXML(sublistfilename, subparams);
          subparams.print(std::cout);
        }
        else
          dserror("The file name '%s' is not a valid XML file name.",
              filename.c_str());
      }
    }
  }

  return;
}
