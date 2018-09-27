/*----------------------------------------------------------------------------*/
/*!
\file nln_utils_debugwriter.cpp

\brief Write output for debugging purposes

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <string>

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

// Teuchos
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

// baci
#include "nln_utils_debugwriter.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::UTILS::DebugWriterBase::DebugWriterBase()
    : isinit_(false),
      issetup_(false),
      name_("unknown"),
      counter_(0),
      log_(Teuchos::null),
      comm_(Teuchos::null),
      step_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterBase::SetComm(const Epetra_Comm& comm)
{
  comm_ = Teuchos::rcp(&comm, false);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterBase::SetOutputFileName(const std::string& name)
{
  name_ = name;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterBase::SetStep(const unsigned int step)
{
  step_ = step;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterBase::CreateLogFile()
{
  // Create log file
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename.append("-");
  filename.append(GetName());
  filename.append("-step-");
  {
    std::stringstream stepstring;
    stepstring << GetStep();
    filename.append(stepstring.str());
  }

  std::string filelog = filename;
  filelog.append(".dbgout");
  log_ = Teuchos::rcp(new std::ofstream(filelog.c_str()));

  // Write header of log-file (only on proc 0)
  if (Comm()->MyPID() == 0)
  {
    (*log_) << "# " << Label() << std::endl
            << "# Encoding of auto-generated name and user-given description for" << std::endl
            << "# debug output file '" << filename << "'" << std::endl
            << "#" << std::endl
            << "# auto\t\t" << std::setw(16) << std::left << "user" << std::endl;

    (*log_) << "#" << std::endl << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterBase::WriteLogFileEntry(
    const std::string& name, const std::string& description) const
{
  if (Comm()->MyPID() == 0)
  {
    (*log_) << name << "\t\t" << std::setw(16) << std::left << description << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const std::string NLNSOL::UTILS::DebugWriterBase::GenerateName() const
{
  // ToDo (mayr) enable user-provided prefix
  std::stringstream uniquename;
  uniquename << "nlndbg_" << counter_;

  // increase counter
  ++counter_;

  return uniquename.str();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Comm> NLNSOL::UTILS::DebugWriterBase::Comm() const
{
  // check if communicator has already been set
  if (comm_.is_null()) dserror("Communicator 'comm_' has not been set, yet.");

  return comm_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::UTILS::DebugWriterSingleField::DebugWriterSingleField()
    : NLNSOL::UTILS::DebugWriterBase(),
      dis_(Teuchos::null),
      control_(Teuchos::null),
      writer_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterSingleField::Init(const Epetra_Comm& comm,
    Teuchos::RCP<const DRT::Discretization> dis, const std::string& name, const unsigned int step)
{
  if (IsInit()) dserror("Init() has already been called. Yoy may call it only once!");

  // Initialize base class members
  SetComm(comm);
  SetOutputFileName(name);
  SetStep(step);

  // Initialize class members
  dis_ = dis;

  // ToDo (mayr) make this a user input
  setDefaultVerbLevel(Teuchos::VERB_MEDIUM);

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterSingleField::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  CreateLogFile();

  // ---------------------------------------------------------------------------
  // Create separate debug control file
  // ---------------------------------------------------------------------------
  // Create unique filename
  std::string filecontrol = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filecontrol.append("-");
  filecontrol.append(GetName());
  filecontrol.append("-step-");
  {
    std::stringstream stepstring;
    stepstring << GetStep();
    filecontrol.append(stepstring.str());
  }

  // Create the control file
  control_ = Teuchos::rcp(new IO::OutputControl(dis_->Comm(),
      //      "none",                   // we do not have a problem type
      "Structure",                       // ToDo (mayr) provide problem type for post filter
      "Polynomial",                      // this is a FE code ... no nurbs
      "debug-output",                    // no input file either
      filecontrol.c_str(),               // an output file name is needed
      DRT::Problem::Instance()->NDim(),  // spatial dimension of problem
      0,                                 // restart is meaningless here
      1000,                              // we never expect to get 1000 iterations
      DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(), "OUTPUT_BIN")));

  // ---------------------------------------------------------------------------

  // Create and initialize debug discretization writer
  writer_ = Teuchos::rcp(
      new IO::DiscretizationWriter(Teuchos::rcp_const_cast<DRT::Discretization>(dis_)));
  writer_->SetOutput(control_);
  writer_->WriteMesh(0, 0.0);
  writer_->NewStep(0, 0.0);
  writer_->WriteElementData(true);  // write ID of Owner proc

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::UTILS::DebugWriterSingleField::WriteVector(Teuchos::RCP<const Epetra_MultiVector> vec,
    const std::string& description, const IO::VectorType vt) const
{
  // Ensure that naming is unique
  const std::string uniquename = GenerateName();

  if (getVerbLevel() > Teuchos::VERB_NONE)
  {
    *getOStream() << "+++ " << LabelShort() << ": writing vector '" << uniquename
                  << "' to output ... ";
  }

  // write vector to debug output
  writer_->WriteVector(uniquename, Teuchos::rcp_const_cast<Epetra_MultiVector>(vec), vt);

  // create entry in log file
  WriteLogFileEntry(uniquename, description);

  if (getVerbLevel() > Teuchos::VERB_NONE)
  {
    *getOStream() << "done" << std::endl;
  }

  return;
}
