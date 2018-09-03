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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// BEGIN OF CLASS /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

NLNSOL::UTILS::NlnConfig::NlnConfig() : params_(Teuchos::null) { return; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::UTILS::NlnConfig::Setup(const std::string filename)
{
  params_ = Teuchos::rcp(new Teuchos::ParameterList());

  NLNSOL::UTILS::CreateParamListFromXML(filename, *params_);

  if (params_.is_null())
  {
    dserror("Filling the parameter list from the XML-file failed.");
    return false;
  }

  return true;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NLNSOL::UTILS::NlnConfig::GetSubList(const std::string listname) const
{
  if (params_.is_null()) dserror("Parameter list 'params_' not properly initialized.");

  if (not params_->isSublist(listname)) dserror("There is no sublist '%s'.", listname.c_str());

  return params_->sublist(listname);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::ParameterList& NLNSOL::UTILS::NlnConfig::GetSubListNonConst(
    const std::string listname) const
{
  if (params_.is_null()) dserror("Parameter list 'params_' not properly initialized.");

  if (not params_->isSublist(listname)) dserror("There is no sublist '%s'.", listname.c_str());

  return params_->sublist(listname);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> NLNSOL::UTILS::NlnConfig::GetAllNonConstRcp() const
{
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList(*params_));
  return params;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// END OF CLASS //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// BEGIN OF CLASS /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::UTILS::StagnationDetection::StagnationDetection()
    : isinit_(false),
      config_(Teuchos::null),
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
void NLNSOL::UTILS::StagnationDetection::Init(Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config,
    const std::string listname, const double norminitial)
{
  // Init() may be called only once
  if (IsInit()) dserror("Init() has already been called. Don't call it again!");

  // initialize
  config_ = config;

  if (config_->GetParameter<bool>(listname, "stagnation detection: on off"))
  {
    active_ = true;
    stagitermax_ = config_->GetParameter<int>(listname, "stagnation detection: max iterations");
    stagthreshold_ =
        config_->GetParameter<double>(listname, "stagnation detection: reduction threshold");
    normprev_ = norminitial;
  }
  else
    active_ = false;

  // configure output to screen
  setVerbLevel(NLNSOL::UTILS::TranslateVerbosityLevelToTeuchos(
      config_->GetParameter<std::string>(listname, "stagnation detection: verbosity")));

  if (getVerbLevel() > Teuchos::VERB_HIGH)
  {
    *getOStream() << "Parmeter List passed to Stagnation Detection: " << std::endl;
  }

  if (getVerbLevel() > Teuchos::VERB_MEDIUM)
  {
    *getOStream() << "Initialized Stagnation Detection with " << std::endl
                  << "  stagnation detection: on off = " << IsActive() << std::endl
                  << "  stagnation detection: max iterations = " << stagiter_ << std::endl
                  << "  stagnation detection: reduction threshold = " << stagthreshold_ << std::endl
                  << "  stagnation detection: verbosity = "
                  << config_->GetParameter<std::string>(listname, "stagnation detection: verbosity")
                  << std::endl;
  }

  // Init() has been called
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::UTILS::StagnationDetection::Check(const double norm)
{
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

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
      *getOStream() << "StagnationCheck: stagiter = " << stagiter_ << "/" << stagitermax_
                    << ", ratio = " << stagratio_ << std::endl;
    }
    else if (getVerbLevel() == Teuchos::VERB_LOW and stagiter_ > 0)
    {
      *getOStream() << "StagnationCheck: stagiter = " << stagiter_ << "/" << stagitermax_
                    << ", ratio = " << stagratio_ << std::endl;
    }
  }

  return Status();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::UTILS::StagnationDetection::Status(Teuchos::RCP<Teuchos::ParameterList> oparams) const
{
  bool stagnation = false;

  if (IsActive() and (stagiter_ >= stagitermax_ or stagratio_ > 1.0)) stagnation = true;

  // fill output parameter list
  if (not oparams.is_null())
  {
    oparams->set<bool>("Stagnation Detection: status", stagnation);
    oparams->set<double>("Stagnation Detection: ratio", stagratio_);
    oparams->set<int>("Stagnation Detection: iterations", stagiter_);
  }

  if (stagnation) *getOStream() << "*** Detected Stagnation ***" << std::endl;

  return stagnation;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> NLNSOL::UTILS::StagnationDetection::StatusParams() const
{
  Teuchos::RCP<Teuchos::ParameterList> oparams = Teuchos::rcp(new Teuchos::ParameterList());

  Status(oparams);

  return oparams;
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
      DRT::Problem::Instance()->NonlinearSolverParams().get<std::string>("XML_FILE");

  // check for reasonable xml file
  if (filename == "none")
    dserror(
        "Seems like you forgot to set the XML file for configuration of "
        "the nonlinear solver.");

  return CreateParamListFromXML(filename);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> NLNSOL::UTILS::CreateParamListFromXML(
    const std::string filename)
{
  // create a new parameter list to be filled
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());

  // fill the parameter list with the content of the given XML-file
  CreateParamListFromXML(filename, *params);

  return params;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::CreateParamListFromXML(
    const std::string filename, Teuchos::ParameterList& params)
{
  // check for reasonable filename of xml-file to be read
  if (filename.length() && filename.rfind(".xml"))
  {
    // fill parameter list with content of xml-file
    params = *(Teuchos::getParametersFromXmlFile(filename));
  }
  else
    dserror("The file name '%s' is not a valid XML file name.", filename.c_str());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::CreateParamListFromXMLOld(
    const std::string filename, Teuchos::ParameterList& params)
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::UTILS::CreateParamListFromXML");
  Teuchos::TimeMonitor monitor(*time);

  // check for reasonable filename of xml-file to be read
  if (filename.length() && filename.rfind(".xml"))
  {
    // fill parameter list with content of xml-file
    params = *(Teuchos::getParametersFromXmlFile(filename));
  }
  else
    dserror("The file name '%s' is not a valid XML file name.", filename.c_str());

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
        sublistfilename.assign(str.begin() + 13, str.end());

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
          dserror("The file name '%s' is not a valid XML file name.", filename.c_str());
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::EVerbosityLevel NLNSOL::UTILS::TranslateVerbosityLevelToTeuchos(
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

#ifdef HAVE_MueLu
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
MueLu::MsgType NLNSOL::UTILS::TranslateVerbosityLevelToMueLu(const std::string verblevelstring)
{
  MueLu::MsgType verblevel = MueLu::Default;

  if (verblevelstring == "default")
    verblevel = MueLu::Default;
  else if (verblevelstring == "none")
    verblevel = MueLu::None;
  else if (verblevelstring == "low")
    verblevel = MueLu::Low;
  else if (verblevelstring == "medium")
    verblevel = MueLu::Medium;
  else if (verblevelstring == "high")
    verblevel = MueLu::High;
  else if (verblevelstring == "extreme")
    verblevel = MueLu::Extreme;
  else
    dserror("Unknown verbosity level '%s'.", verblevelstring.c_str());

  return verblevel;
}
#endif
