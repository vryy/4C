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
#include <Teuchos_FancyOStream.hpp>
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
  params_(Teuchos::null),
  stagiter_(1),
  stagitermax_(0),
  stagratio_(0.0),
  stagthreshold_(1.0),
  normprev_(1.0e+12),
  active_(false)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::StagnationDetection::Init(
    const Teuchos::ParameterList& params,
    const double norminitial)
{
  // Init() may be called only once
  if (IsInit())
    dserror("Init() has already been called. Don't call it again!");

  // initialize
  params_ = Teuchos::rcp(&params, false);

  if (Params()->get<bool>("Stagnation Detection: on off"))
  {
    active_ = true;
    stagitermax_ = Params()->get<int>("Stagnation Detection: max iterations");
    stagthreshold_ = Params()->get<double>("Stagnation Detection: reduction threshold");
    normprev_ = norminitial;
  }
  else
    active_ = false;

  // configure output to screen
  setVerbLevel(
      NLNSOL::UTILS::TranslateVerbosityLevel(
          Params()->get<std::string>("Stagnation Detection: Verbosity")));

  if (getVerbLevel() > Teuchos::VERB_HIGH)
  {
    *getOStream() << "Parmeter List passed to Stagnation Detection: "
        << std::endl;
    Params()->print(*getOStream());
  }

  if (getVerbLevel() > Teuchos::VERB_MEDIUM)
  {
    *getOStream() << "Initialized Stagnation Detection with " << std::endl
        << "  Stagnation Detection: on off = " << IsActive() << std::endl
        << "  Stagnation Detection: max iterations = " << stagiter_ << std::endl
        << "  Stagnation Detection: reduction threshold = " << stagthreshold_
        << std::endl
        << "  Stagnation Detection: Verbosity = "
        << Params()->get<std::string>("Stagnation Detection: Verbosity")
        << std::endl;
  }

  // Init() has been called
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const bool NLNSOL::UTILS::StagnationDetection::Check(const double norm)
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  if (IsActive())
  {
    stagratio_ = norm / normprev_;
    normprev_ = norm;

    // ---------------------------------------------------------------------------
    // decide whether this is considered as stagnation
    // ---------------------------------------------------------------------------
    if (stagratio_ > stagthreshold_)
      ++stagiter_;
    else
      stagiter_ = 0;

    if (getVerbLevel() > Teuchos::VERB_LOW)
    {
      *getOStream() << "StagnationCheck: stagiter = " << stagiter_ << "/"
          << stagitermax_ << ", ratio = " << stagratio_ << std::endl;
    }
    else if (getVerbLevel() == Teuchos::VERB_LOW and stagiter_ > 0)
    {
      *getOStream() << "StagnationCheck: stagiter = " << stagiter_ << "/"
          << stagitermax_ << ", ratio = " << stagratio_ << std::endl;
    }
  }

  return Status();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const bool NLNSOL::UTILS::StagnationDetection::Status(
    Teuchos::RCP<Teuchos::ParameterList> oparams) const
{
  bool stagnation = false;

  if (IsActive() and (stagiter_ >= stagitermax_ or stagratio_ > 1.0))
    stagnation = true;

  // fill output parameter list
  if (not oparams.is_null())
  {
    oparams->set<bool>("Stagnation Detection: status", stagnation);
    oparams->set<double>("Stagnation Detection: ratio", stagratio_);
    oparams->set<int>("Stagnation Detection: iterations", stagiter_);
  }

  if (stagnation)
    *getOStream() << "*** Detected Stagnation ***" << std::endl;

  return stagnation;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList>
NLNSOL::UTILS::StagnationDetection::StatusParams() const
{
  Teuchos::RCP<Teuchos::ParameterList> oparams =
      Teuchos::rcp(new Teuchos::ParameterList());

  Status(oparams);

  return oparams;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList>
NLNSOL::UTILS::StagnationDetection::Params() const
{
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been set, yet.");

  return params_;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// END OF CLASS //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Teuchos::EVerbosityLevel NLNSOL::UTILS::TranslateVerbosityLevel(
    const std::string verblevelstring)
{
  Teuchos::EVerbosityLevel verblevel = Teuchos::VERB_DEFAULT;

  if (verblevelstring == "default")
    verblevel = Teuchos::VERB_DEFAULT;
  else if (verblevelstring == "none")
    verblevel = Teuchos::VERB_NONE;
  else if (verblevelstring == "low")
    verblevel = Teuchos::VERB_LOW;
  else if (verblevelstring == "medium")
    verblevel = Teuchos::VERB_MEDIUM;
  else if (verblevelstring == "high")
    verblevel = Teuchos::VERB_HIGH;
  else if (verblevelstring == "extreme")
    verblevel = Teuchos::VERB_EXTREME;
  else
    dserror("Unknown verbosity level '%s'.", verblevelstring.c_str());

  return verblevel;
}
