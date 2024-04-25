/*-----------------------------------------------------------*/
/*! \file

\brief Kind of a data container class, which holds many variables
       and objects, which are necessary to setup a NOX::NLN
       solution strategy.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_globaldata.hpp"  // class definition

#include "4C_inpar_boolifyparameters.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_direction_factory.hpp"
#include "4C_solver_nonlin_nox_linearsystem.hpp"
#include "4C_solver_nonlin_nox_meritfunction_factory.hpp"
#include "4C_solver_nonlin_nox_solver_prepostop_generic.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <stdexcept>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::OptimizationProblemType& type,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
    const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
    const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
    const Teuchos::RCP<::NOX::Epetra::Scaling>& iScale)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      opt_type_(type),
      lin_solvers_(linSolvers),
      i_req_ptr_(iReq),
      i_jac_ptr_(iJac),
      i_prec_ptr_(iPrec),
      i_constr_(iConstr),
      i_constr_prec_(iConstrPrec),
      i_scale_(iScale),
      mrt_fct_ptr_(Teuchos::null),
      pre_post_op_ptr_(Teuchos::null),
      is_constrained_(type != opt_unconstrained)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const OptimizationProblemType& type, const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      opt_type_(type),
      lin_solvers_(linSolvers),
      i_req_ptr_(iReq),
      i_jac_ptr_(iJac),
      i_prec_ptr_(Teuchos::null),  // no pre-conditioner
      i_constr_(iConstr),
      mrt_fct_ptr_(Teuchos::null),
      pre_post_op_ptr_(Teuchos::null),
      is_constrained_(type != opt_unconstrained)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      opt_type_(opt_unconstrained),
      lin_solvers_(linSolvers),
      i_req_ptr_(iReq),
      i_jac_ptr_(iJac),
      i_prec_ptr_(iPrec),
      mrt_fct_ptr_(Teuchos::null),
      pre_post_op_ptr_(Teuchos::null),
      is_constrained_(false)
{
  CheckInput();
  // do some setup things
  Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GlobalData::GlobalData(const Epetra_Comm& comm, Teuchos::ParameterList& noxParams,
    const NOX::NLN::LinearSystem::SolverMap& linSolvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac)
    : comm_(Teuchos::rcp(&comm, false)),
      nlnparams_(Teuchos::rcp(&noxParams, false)),
      opt_type_(opt_unconstrained),
      lin_solvers_(linSolvers),
      i_req_ptr_(iReq),
      i_jac_ptr_(iJac),
      i_prec_ptr_(Teuchos::null),  // no pre-conditioner
      mrt_fct_ptr_(Teuchos::null),
      pre_post_op_ptr_(Teuchos::null),
      is_constrained_(false)
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
  if (lin_solvers_.size() == 0)
    FOUR_C_THROW("The linear solver map has the size 0! Required size > 0.");

  typedef std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>::const_iterator CI;
  for (CI iter = lin_solvers_.begin(); iter != lin_solvers_.end(); ++iter)
    if (iter->second.is_null())
    {
      std::ostringstream msg;
      const std::string act_msg("        WARNING: The entry \"" + SolutionType2String(iter->first) +
                                "\" of the linear solver vector is not (yet) initialized!");
      msg << std::setw(109) << std::setfill('!') << "\n" << std::setfill(' ');
      msg << "!!! " << std::left << std::setw(100) << act_msg << " !!!";
      msg << std::setw(109) << std::setfill('!') << "\n" << std::setfill(' ') << std::right;
      if (comm_->MyPID() == 0) std::cout << msg.str() << std::endl;
    }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GlobalData::Setup()
{
  // set the nonlinear optimzation problem type
  nlnparams_->set("Optimization Problem Type", opt_type_);

  // set printing parameters
  SetPrintingParameters();

  // construct the nox utils class
  nox_utils_ = Teuchos::rcp(new ::NOX::Utils(nlnparams_->sublist("Printing")));

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
   * a Teuchos::RCP<::NOX::MeritFunction::Generic> pointer into the parameter
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
    mrt_fct_ptr_ = NOX::NLN::MeritFunction::BuildMeritFunction(*this);
  }

  // If the mrtFctPtr is Teuchos::null the default "Sum of Squares" NOX internal
  // merit function is used.
  if (not mrt_fct_ptr_.is_null())
    solverOptionsList.set<Teuchos::RCP<::NOX::MeritFunction::Generic>>(
        "User Defined Merit Function", mrt_fct_ptr_);

  /* NOX has become quite strict about parameter lists. They may only contain parameters,
   * that are then actually used. Hence, we have to remove some.
   *
   * ToDo: Clean this up by not setting unused parameters in the first place.
   */
  {
    // remove "Merit Function" parameter entry after building the user-defined function
    // (NOX does not accept this parameter)
    solverOptionsList.remove("Merit Function", false);

    // also remove *PrintEqualSign*
    solverOptionsList.remove("*PrintEqualSign*", false);
  }

  // Similar to the merit function we provide our own direction factory, if
  // the direction method is set to "User Defined"
  {
    Teuchos::ParameterList& pdir = nlnparams_->sublist("Direction", true);

    // Set the direction method to user-defined if the method is set to Newton,
    // since we want to use our own NOX::NLN::Newton direction method and not the
    // one provided by NOX.
    if (pdir.get<std::string>("Method") == "Newton")
    {
      pdir.set("Method", "User Defined");
      pdir.set("User Defined Method", "Newton");
    }

    // Set the associated factory
    if (pdir.get<std::string>("Method") == "User Defined")
    {
      direction_factory_ = Teuchos::rcp(new NOX::NLN::Direction::Factory);
      pdir.set<Teuchos::RCP<::NOX::Direction::UserDefinedFactory>>(
          "User Defined Direction Factory", direction_factory_);
    }
  }

  /* We use the parameter list to define a PrePostOperator class for the
   * non-linear iteration process. */
  pre_post_op_ptr_ = Teuchos::rcp(new NOX::NLN::Solver::PrePostOp::Generic());
  NOX::NLN::AUX::AddToPrePostOpVector(solverOptionsList, pre_post_op_ptr_);

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
    FOUR_C_THROW("Please specify a \"Status Test\"->\"XML File\" for the nox nonlinear solver!");

  if (xmlfilename.length() && xmlfilename.rfind(".xml"))
  {
    try
    {
      xmlParams = *(Teuchos::getParametersFromXmlFile(xmlfilename));
    }
    catch (const std::runtime_error&)
    {
      FOUR_C_THROW(
          "The \"Status Test\"->\"XML File\" was not found! "
          "Please check the path in your input file! \n"
          "CURRENT PATH = %s",
          xmlfilename.c_str());
    }
  }
  else
    FOUR_C_THROW("The file name '%s' is not a valid XML file name.", xmlfilename.c_str());

  // copy the "Outer Status Test" into the nox parameter list
  if (not xmlParams.isSublist("Outer Status Test"))
    FOUR_C_THROW("You have to specify a sublist \"Outer Status Test\" in your xml-file.");
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
  INPUT::BoolifyValidInputParameters(statusTestParams);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Utils& NOX::NLN::GlobalData::GetNoxUtils() const
{
  if (nox_utils_.is_null()) FOUR_C_THROW("noxUtils_ was not initialized!");

  return *nox_utils_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<::NOX::Utils>& NOX::NLN::GlobalData::GetNoxUtilsPtr() const
{
  if (nox_utils_.is_null()) FOUR_C_THROW("noxUtils_ was not initialized!");

  return nox_utils_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NOX::NLN::GlobalData::GetNlnParameterList() const
{
  if (nlnparams_.is_null()) FOUR_C_THROW("nlnparams_ was not initialized!");

  return *nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& NOX::NLN::GlobalData::GetNlnParameterList()
{
  if (nlnparams_.is_null()) FOUR_C_THROW("nlnparams_ was not initialized!");

  return *nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<Teuchos::ParameterList>& NOX::NLN::GlobalData::GetNlnParameterListPtr()
{
  if (nlnparams_.is_null()) FOUR_C_THROW("nlnparams_ was not initialized!");

  return nlnparams_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& NOX::NLN::GlobalData::GetComm() const
{
  if (comm_.is_null()) FOUR_C_THROW("comm_ was not initialized!");

  return *comm_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::LinearSystem::SolverMap& NOX::NLN::GlobalData::GetLinSolvers()
{
  return lin_solvers_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::Interface::Required> NOX::NLN::GlobalData::GetRequiredInterface()
{
  if (i_req_ptr_.is_null())
    FOUR_C_THROW("Required interface pointer iReqPtr_ was not initialized!");

  return i_req_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::Interface::Jacobian> NOX::NLN::GlobalData::GetJacobianInterface()
{
  if (i_jac_ptr_.is_null())
    FOUR_C_THROW("Jacobian interface pointer iJacPtr_ was not initialized!");

  return i_jac_ptr_;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>
NOX::NLN::GlobalData::GetPreconditionerInterface()
{
  /* We explicitly allow a return value of Teuchos::nullptr, because the
   * preconditioner interface is in many cases optional */
  return i_prec_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::CONSTRAINT::ReqInterfaceMap& NOX::NLN::GlobalData::GetConstraintInterfaces()
{
  if (i_constr_.size() == 0) FOUR_C_THROW("The constraint interface map is empty!");

  return i_constr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::CONSTRAINT::PrecInterfaceMap& NOX::NLN::GlobalData::GetConstraintPrecInterfaces()
{
  if (i_constr_prec_.size() == 0)
    FOUR_C_THROW("The constraint preconditioner interface map is empty!");

  return i_constr_prec_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const bool& NOX::NLN::GlobalData::GetIsConstrained() const { return is_constrained_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<::NOX::Epetra::Scaling>& NOX::NLN::GlobalData::GetScalingObject()
{
  return i_scale_;
}

FOUR_C_NAMESPACE_CLOSE
