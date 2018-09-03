/*-----------------------------------------------------------*/
/*!
\file nox_nln_globaldata.cpp

\brief Kind of a data container class, which holds many variables
       and objects, which are necessary to setup a NOX::NLN
       solution strategy.

\maintainer Michael Hiermeier

\date Jul 17, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_globaldata.H"  // class definition
#include "nox_nln_meritfunction_factory.H"
#include "nox_nln_solver_prepostop_generic.H"
#include "nox_nln_linearsystem.H"
#include "nox_nln_aux.H"

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
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::OptimizationProblemType& type,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
    const Teuchos::RCP<NOX::Epetra::Scaling>& iScale)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      optType_(type),
      linSolvers_(linSolvers),
      iReqPtr_(iReq),
      iJacPtr_(iJac),
      iPrecPtr_(iPrec),
      iConstr_(iConstr),
      iConstrPrec_(iConstrPrec),
      iScale_(iScale),
      mrtFctPtr_(Teuchos::null),
      prePostOpPtr_(Teuchos::null),
      isConstrained_(type != opt_unconstrained)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, const OptimizationProblemType& type,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      optType_(type),
      linSolvers_(linSolvers),
      iReqPtr_(iReq),
      iJacPtr_(iJac),
      iPrecPtr_(Teuchos::null),  // no pre-conditioner
      iConstr_(iConstr),
      mrtFctPtr_(Teuchos::null),
      prePostOpPtr_(Teuchos::null),
      isConstrained_(type != opt_unconstrained)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      optType_(opt_unconstrained),
      linSolvers_(linSolvers),
      iReqPtr_(iReq),
      iJacPtr_(iJac),
      iPrecPtr_(iPrec),
      mrtFctPtr_(Teuchos::null),
      prePostOpPtr_(Teuchos::null),
      isConstrained_(false)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      optType_(opt_unconstrained),
      linSolvers_(linSolvers),
      iReqPtr_(iReq),
      iJacPtr_(iJac),
      iPrecPtr_(Teuchos::null),  // no pre-conditioner
      mrtFctPtr_(Teuchos::null),
      prePostOpPtr_(Teuchos::null),
      isConstrained_(false)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::CheckInput() const
{
  // short input check
  if (linSolvers_.size() == 0) dserror("The linear solver map has the size 0! Required size > 0.");

  typedef std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>::const_iterator CI;
  for (CI iter = linSolvers_.begin(); iter != linSolvers_.end(); ++iter)
    if (iter->second.is_null())
    {
      std::ostringstream msg;
      msg << "The entry \"" << SolutionType2String(iter->first)
          << "\" of the linear solver vector is not initialized!";
      dserror(msg.str());
    }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::Setup()
{
  // set the nonlinear optimzation problem type
  nlnparams_->set("Optimization Problem Type", optType_);

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
  NOX::NLN::AUX::SetPrintingParameters(*nlnparams_, *comm_);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::SetSolverOptionParameters()
{
  /* NOX gives you the option to use a user defined merit function by inserting
   * a Teuchos::RCP<NOX::MeritFunction::Generic> pointer into the parameter
   * sublist "Solver Options". Since we are forced to use the
   * NOX_GlobalData class, we have to use this to define our own merit functions
   * by creating the particular class in the NOX::NLN::MeritFunction::Factory
   * and insert it into the parameter list. */
  Teuchos::ParameterList& solverOptionsList = nlnparams_->sublist("Solver Options");

  /* If we do a full newton method, we don't need a merit function and can
   * skip the following */
  const std::string& nlnsolver_name = nlnparams_->get<std::string>("Nonlinear Solver");
  const std::string& ls_method_name = nlnparams_->sublist("Line Search").get<std::string>("Method");
  if (((nlnsolver_name == "Line Search Based") or (nlnsolver_name == "Pseudo Transient")) and
      (ls_method_name != "Full Step"))
  {
    // Pure reading access to the unfinished nox_nln_globaldata class
    mrtFctPtr_ = NOX::NLN::MeritFunction::BuildMeritFunction(*this);
  }

  // If the mrtFctPtr is Teuchos::null the default "Sum of Squares" NOX internal
  // merit function is used.
  if (not mrtFctPtr_.is_null())
    solverOptionsList.set<Teuchos::RCP<NOX::MeritFunction::Generic>>(
        "User Defined Merit Function", mrtFctPtr_);

  /* We use the parameter list to define a PrePostOperator class for the
   * non-linear iteration process. */
  prePostOpPtr_ = Teuchos::rcp(new NOX::NLN::Solver::PrePostOp::Generic());
  NOX::NLN::AUX::AddToPrePostOpVector(solverOptionsList, prePostOpPtr_);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::SetStatusTestParameters()
{
  Teuchos::ParameterList& statusTestParams = nlnparams_->sublist("Status Test", true);

  // check if the status test was already set via the dat-file
  if (statusTestParams.isSublist("Outer Status Test") and
      statusTestParams.sublist("Outer Status Test").numParams() != 0)
    return;

  // initialize the sublist "Outer Status Test". This one has to be specified.
  Teuchos::ParameterList& outerStatusTestParams =
      statusTestParams.sublist("Outer Status Test", false);

  Teuchos::ParameterList xmlParams;
  std::string xmlfilename = statusTestParams.get<std::string>("XML File");

  // check the input: path to the "Status Test" xml-file
  if (xmlfilename == "none")
    dserror("Please specify a \"Status Test\"->\"XML File\" for the nox nonlinear solver!");

  if (xmlfilename.length() && xmlfilename.rfind(".xml"))
  {
    try
    {
      xmlParams = *(Teuchos::getParametersFromXmlFile(xmlfilename));
    }
    catch (std::runtime_error& e)
    {
      dserror(
          "The \"Status Test\"->\"XML File\" was not found! "
          "Please check the path in your input file! \n"
          "CURRENT PATH = %s",
          xmlfilename.c_str());
    }
  }
  else
    dserror("The file name '%s' is not a valid XML file name.", xmlfilename.c_str());

  // copy the "Outer Status Test" into the nox parameter list
  if (not xmlParams.isSublist("Outer Status Test"))
    dserror("You have to specify a sublist \"Outer Status Test\" in your xml-file.");
  else
    outerStatusTestParams = xmlParams.sublist("Outer Status Test");

  /* If we use a line search based backtracking non-linear solver and an inner-status test
   * list is necessary. We check it during the nox_nln_linesearch_factory call. */
  if (xmlParams.isSublist("Inner Status Test"))
  {
    Teuchos::ParameterList& innerStatusTestParams =
        statusTestParams.sublist("Inner Status Test", false);
    // copy the "Inner Status Test" into the nox parameter list
    if (xmlParams.sublist("Inner Status Test").numParams())
      innerStatusTestParams = xmlParams.sublist("Inner Status Test");
  }

  // make all Yes/No integral values to Boolean
  DRT::INPUT::BoolifyValidInputParameters(statusTestParams);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Utils& NOX::NLN::GlobalData::GetNoxUtils() const
{
  if (noxUtils_.is_null()) dserror("noxUtils_ was not initialized!");

  return *noxUtils_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Utils>& NOX::NLN::GlobalData::GetNoxUtilsPtr() const
{
  if (noxUtils_.is_null()) dserror("noxUtils_ was not initialized!");

  return noxUtils_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NOX::NLN::GlobalData::GetNlnParameterList() const
{
  if (nlnparams_.is_null()) dserror("nlnparams_ was not initialized!");

  return *nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& NOX::NLN::GlobalData::GetNlnParameterList()
{
  if (nlnparams_.is_null()) dserror("nlnparams_ was not initialized!");

  return *nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<Teuchos::ParameterList>& NOX::NLN::GlobalData::GetNlnParameterListPtr()
{
  if (nlnparams_.is_null()) dserror("nlnparams_ was not initialized!");

  return nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& NOX::NLN::GlobalData::GetComm() const
{
  if (comm_.is_null()) dserror("comm_ was not initialized!");

  return *comm_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::LinearSystem::SolverMap& NOX::NLN::GlobalData::GetLinSolvers()
{
  return linSolvers_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Interface::Required> NOX::NLN::GlobalData::GetRequiredInterface()
{
  if (iReqPtr_.is_null()) dserror("Required interface pointer iReqPtr_ was not initialized!");

  return iReqPtr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> NOX::NLN::GlobalData::GetJacobianInterface()
{
  if (iJacPtr_.is_null()) dserror("Jacobian interface pointer iJacPtr_ was not initialized!");

  return iJacPtr_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>
NOX::NLN::GlobalData::GetPreconditionerInterface()
{
  /* We explicitly allow a return value of Teuchos::NULL, because the
   * preconditioner interface is in many cases optional */
  return iPrecPtr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::CONSTRAINT::ReqInterfaceMap& NOX::NLN::GlobalData::GetConstraintInterfaces()
{
  if (iConstr_.size() == 0) dserror("The constraint interface map is empty!");

  return iConstr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::CONSTRAINT::PrecInterfaceMap& NOX::NLN::GlobalData::GetConstraintPrecInterfaces()
{
  if (iConstrPrec_.size() == 0) dserror("The constraint preconditioner interface map is empty!");

  return iConstrPrec_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const bool& NOX::NLN::GlobalData::GetIsConstrained() const { return isConstrained_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::Epetra::Scaling>& NOX::NLN::GlobalData::GetScalingObject()
{
  return iScale_;
}
