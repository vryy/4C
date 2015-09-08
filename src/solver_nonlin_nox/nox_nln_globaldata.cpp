/*-----------------------------------------------------------*/
/*!
\file nox_nln_globaldata.cpp

\maintainer Michael Hiermeier

\date Jul 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_globaldata.H"   // class definition
#include "nox_nln_constraint_interface_required.H"
#include "nox_nln_meritfunction_factory.H"
#include "nox_nln_prepostoperator.H"

#include <NOX_Utils.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>

#include "../linalg/linalg_solver.H"

#include "../drt_inpar/inpar_structure.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/drt_boolifyparameters.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm,
    const Teuchos::ParameterList& noxparams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::GlobalData::OptimizationProblemType& type,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>& iConstr) :
    comm_(Teuchos::rcp(& comm,false)),
    nlnparams_(Teuchos::rcp(new Teuchos::ParameterList(noxparams))),
    optType_(type),
    linSolvers_(linSolvers),
    iReqPtr_(iReq),
    iJacPtr_(iJac),
    iPrecPtr_(iPrec),
    iConstrPtr_(iConstr),
    mrtFctPtr_(Teuchos::null),
    prePostOpPtr_(Teuchos::null),
    isConstrained_(type!=opt_unconstrained)
{
  // short input check
  if (linSolvers_.size()==0)
    dserror("The linear solver map has the size 0! Required size > 0.");

  typedef std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >::const_iterator CI;
  for (CI iter=linSolvers_.begin(); iter!=linSolvers_.end(); ++iter)
    if (iter->second.is_null())
    {
      std::ostringstream msg;
      msg << "The entry \"" << SolutionType2String(iter->first) << "\" of the linear solver vector is not initialized!";
      dserror(msg.str());
    }

  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::Setup()
{
  // set the nonlinear optimzation problem type
  nlnparams_->set("Optimization Problem Type",optType_);

  // set printing parameters
  SetPrintingParameters();

  // construct the nox utils class
  noxUtils_ = Teuchos::rcp(new NOX::Utils(nlnparams_->sublist("Printing")));

  // set (non-linear) solver option parameters
  SetSolverOptionParameters();

  // set status test parameters
  SetStatusTestParameters();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::SetPrintingParameters()
{
// make all Yes/No integral values to Boolean
  DRT::INPUT::BoolifyValidInputParameters(*nlnparams_);

  // adjust printing parameter list
  Teuchos::ParameterList& printParams = nlnparams_->sublist("Printing");
  printParams.set("MyPID", comm_->MyPID());
  printParams.set("Output Precision", 5);
  printParams.set("Output Processor", 0);
  int outputinformationlevel = NOX::Utils::Error;  // NOX::Utils::Error==0
  if (printParams.get<bool>("Error"))
    outputinformationlevel += NOX::Utils::Error;
  if (printParams.get<bool>("Warning"))
    outputinformationlevel += NOX::Utils::Warning;
  if (printParams.get<bool>("Outer Iteration"))
    outputinformationlevel += NOX::Utils::OuterIteration;
  if (printParams.get<bool>("Inner Iteration"))
    outputinformationlevel += NOX::Utils::InnerIteration;
  if (printParams.get<bool>("Parameters"))
    outputinformationlevel += NOX::Utils::Parameters;
  if (printParams.get<bool>("Details"))
    outputinformationlevel += NOX::Utils::Details;
  if (printParams.get<bool>("Outer Iteration StatusTest"))
    outputinformationlevel += NOX::Utils::OuterIterationStatusTest;
  if (printParams.get<bool>("Linear Solver Details"))
    outputinformationlevel += NOX::Utils::LinearSolverDetails;
  if (printParams.get<bool>("Test Details"))
    outputinformationlevel += NOX::Utils::TestDetails;
  /*  // for LOCA
  if (printParams.get<bool>("Stepper Iteration"))
    outputinformationlevel += NOX::Utils::StepperIteration;
  if (printParams.get<bool>("Stepper Details"))
    outputinformationlevel += NOX::Utils::StepperDetails;
  if (printParams.get<bool>("Stepper Parameters"))
    outputinformationlevel += NOX::Utils::StepperParameters;
  */
  if (printParams.get<bool>("Debug"))
    outputinformationlevel += NOX::Utils::Debug;
  printParams.set("Output Information", outputinformationlevel);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::SetSolverOptionParameters()
{
  // NOX gives you the option to use a user defined merit function by inserting
  // a Teuchos::RCP<NOX::MeritFunction::Generic> pointer into the parameter
  // sublist "Solver Options". Since we are forced to use the
  // NOX_GlobalData class, we have to use this to define our own merit functions
  // by creating the particular class in the NOX::NLN::MeritFunction::Factory
  // and insert it into the parameter list.
  Teuchos::ParameterList& solverOptionsList = nlnparams_->sublist("Solver Options");

  // Pure reading access to the unfinished nox_nln_utils class
  Teuchos::RCP<const NOX::NLN::GlobalData> nlnGlobalDataPtr = Teuchos::rcp(this,false);
  mrtFctPtr_ =
      NOX::NLN::MeritFunction::BuildMeritFunction(nlnGlobalDataPtr);

  // If the mrtFctPtr is Teuchos::null the default "Sum of Squares" NOX internal
  // merit function is used.
  if (not mrtFctPtr_.is_null())
    solverOptionsList.set<Teuchos::RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function",mrtFctPtr_);

  // We use the parameter list to define a PrePostOperator class for the
  // non-linear iteration process.
  prePostOpPtr_ = Teuchos::rcp(new NOX::NLN::PrePostOperator());
  solverOptionsList.set<Teuchos::RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator",prePostOpPtr_);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::SetStatusTestParameters()
{
  Teuchos::ParameterList& statusTestParams = nlnparams_->sublist("Status Test",true);
  // initialize the sublist "Outer Status Test". This one has to be specified.
  Teuchos::ParameterList& outerStatusTestParams = statusTestParams.sublist("Outer Status Test",false);

  Teuchos::ParameterList xmlParams;
  std::string xmlfilename = statusTestParams.get<std::string>("XML File");

  // check the input: path to the "Status Test" xml-file
  if (xmlfilename == "none")
    dserror("Please specify a \"Status Test\"->\"XML File\" for the nox nonlinear solver!");

  if (xmlfilename.length() && xmlfilename.rfind(".xml"))
  {
    xmlParams = *(Teuchos::getParametersFromXmlFile(xmlfilename));
  }
  else
    dserror("The file name '%s' is not a valid XML file name.",
        xmlfilename.c_str());

  // copy the "Outer Status Test" into the nox parameter list
  if (not xmlParams.isSublist("Outer Status Test"))
    dserror("You have to specify a sublist \"Outer Status Test\" in your xml-file.");
  else
    outerStatusTestParams = xmlParams.sublist("Outer Status Test");

  // If  we use a line search based backtracking non-linear solver a
  // a inner-StatusTest list is necessary.
  if ((nlnparams_->get<std::string>("Nonlinear Solver")=="Line Search Based") and
      (nlnparams_->sublist("Line Search").get<std::string>("Method")=="Backtrack"))
  {
    Teuchos::ParameterList& innerStatusTestParams = statusTestParams.sublist("Inner Status Test",false);

    // copy the "Inner Status Test" into the nox parameter list
    if (not xmlParams.isSublist("Inner Status Test") or
        xmlParams.sublist("Inner Status Test").numParams()==0)
      dserror("You have to specify a filled sublist \"Inner Status Test\" in your xml-file to use a "
          "\"Line Search Based\"-\"Backtrack\" method.");
    else
      innerStatusTestParams = xmlParams.sublist("Inner Status Test");
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Utils& NOX::NLN::GlobalData::GetNoxUtils() const
{
  if (noxUtils_.is_null())
    dserror("noxUtils_ was not initialized!");

  return *noxUtils_;
}

const Teuchos::RCP<NOX::Utils>& NOX::NLN::GlobalData::GetNoxUtilsPtr() const
{
  if (noxUtils_.is_null())
    dserror("noxUtils_ was not initialized!");

  return noxUtils_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NOX::NLN::GlobalData::GetNlnParameterList() const
{
  if (nlnparams_.is_null())
    dserror("nlnparams_ was not initialized!");

  return *nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& NOX::NLN::GlobalData::GetNlnParameterList()
{
  if (nlnparams_.is_null())
    dserror("nlnparams_ was not initialized!");

  return *nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<Teuchos::ParameterList>& NOX::NLN::GlobalData::GetNlnParameterListPtr()
{
  if (nlnparams_.is_null())
    dserror("nlnparams_ was not initialized!");

  return nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& NOX::NLN::GlobalData::GetComm() const
{
  if (comm_.is_null())
    dserror("comm_ was not initialized!");

  return *comm_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& NOX::NLN::GlobalData::GetLinSolvers()
{
  return linSolvers_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Interface::Required> NOX::NLN::GlobalData::GetRequiredInterface()
{
  if (iReqPtr_.is_null())
    dserror("Required interface pointer iReqPtr_ was not initialized!");

  return iReqPtr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> NOX::NLN::GlobalData::GetJacobianInterface()
{
  if (iJacPtr_.is_null())
    dserror("Jacobian interface pointer iJacPtr_ was not initialized!");

  return iJacPtr_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> NOX::NLN::GlobalData::GetPreconditionerInterface()
{
  // we allow a return value of Teuchos::NULL, because the preconditioner
  // interface is in many case optional

  return iPrecPtr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> NOX::NLN::GlobalData::GetConstraintInterface()
{
  if (iConstrPtr_.is_null())
    dserror("Constraint interface pointer iConstrPtr_ was not initialized!");

  return iConstrPtr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const bool& NOX::NLN::GlobalData::GetIsConstrained() const
{
  return isConstrained_;
}
