/*-----------------------------------------------------------*/
/*! \file

\brief Library of Continuation Algorithms (LOCA) implementation.

\maintainer Anh-Tu Vuong

\date Nov 18, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_loca_continuation.H"
#include "str_nln_solver_utils.H"
#include "str_timint_locainterface.H"
#include "str_utils.H"
#include "str_impl_generic.H"

#include "../solver_nonlin_nox/nox_nln_globaldata.H"
#include "../loca_continuation/loca_nln_problem.H"

#include <LOCA_MultiContinuation_AbstractGroup.H>
#include <LOCA_Epetra_Factory.H>
#include <LOCA_Parameter_Vector.H>
#include <LOCA_MultiContinuation_ConstraintInterface.H>

#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_LinearSystem.H>
#include "str_nln_linearsystem_scaling.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::LOCAContinuation::LOCAContinuation()
    : p_loca_nox_ptr_(Teuchos::null),
      loca_constraints_ptr_(Teuchos::null),
      loca_stepper_ptr_(Teuchos::null),
      loca_stepper_status_(LOCA::Abstract::Iterator::NotFinished)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::Setup()
{
  CheckInit();
  // call base class setup routine
  Base::Setup();

  CombineParameterLists();
  Teuchos::RCP<LOCA::ParameterVector> loca_param_vec_ptr =
      CreateParameterVector(p_loca_nox_ptr_->sublist("LOCA", false));
  CheckForNecessaryParameters();
  Teuchos::ParameterList& p_loca = p_loca_nox_ptr_->sublist("LOCA");
  Teuchos::ParameterList& p_nox = p_loca_nox_ptr_->sublist("NOX");

  // ---------------------------------------------------------------------------
  // create loca constraint interface
  // ---------------------------------------------------------------------------
  /* Seems to be unnecessary: See LOCA::MultiContinuation::ArcLengthConstraint
   * for an example.
   * loca_constraints_ptr_ =
   *     Teuchos::rcp(new LOCA::STR::MultiContinuation::ArcLengthConstraint()); */
  // ToDo extend this functionality (multi constraint case)
  loca_constraint_param_names_ptr_ = Teuchos::rcp(new std::vector<std::string>(1));
  (*loca_constraint_param_names_ptr_)[0] = "displ_fac";

  // set the constraints in the LOCA parameter list
  Teuchos::ParameterList& p_constraints = p_loca.sublist("Constraints");
  p_constraints.set("Constraint Object", loca_constraints_ptr_);
  p_constraints.set("Constraint Parameter Names", loca_constraint_param_names_ptr_);

  // ---------------------------------------------------------------------------
  // create the user defined factory
  // ToDo extend to a LOCA NLN factory if necessary
  // ---------------------------------------------------------------------------
  Teuchos::RCP<LOCA::Epetra::Factory> loca_epetra_factory_ptr =
      Teuchos::rcp(new LOCA::Epetra::Factory());

  // ---------------------------------------------------------------------------
  // create loca global date object
  // ---------------------------------------------------------------------------
  Teuchos::RCP<LOCA::GlobalData> loca_global_data_ptr =
      LOCA::createGlobalData(p_loca_nox_ptr_, loca_epetra_factory_ptr);

  // ---------------------------------------------------------------------------
  // create the interface for the internal nox and loca compute routines
  // ---------------------------------------------------------------------------
  // cast the base class integrator
  Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr =
      Teuchos::rcp_dynamic_cast<STR::IMPLICIT::Generic>(IntegratorPtr());
  Teuchos::RCP<STR::TIMINT::LocaInterface> loca_interface_ptr =
      Teuchos::rcp(new STR::TIMINT::LocaInterface());
  loca_interface_ptr->Init(DataGlobalStatePtr(), implint_ptr, DBCPtr(), Teuchos::rcp(this, false));
  loca_interface_ptr->Setup();

  // ---------------------------------------------------------------------------
  // create nox nln global data
  // ---------------------------------------------------------------------------
  std::vector<enum NOX::NLN::SolutionType> soltypes(0);
  std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>> linsolvers;
  // convert the INPAR::STR::ModelType to a NOX::NLN::SolType
  STR::NLN::ConvertModelType2SolType(
      soltypes, linsolvers, DataSDyn().GetModelTypes(), DataSDyn().GetLinSolvers());
  // define and initialize the optimization type
  const NOX::NLN::OptimizationProblemType opttype = STR::NLN::OptimizationType(soltypes);

  // set constraint interfaces
  NOX::NLN::CONSTRAINT::ReqInterfaceMap iconstr;
  STR::NLN::CreateConstraintInterfaces(iconstr, Integrator(), soltypes);

  // preconditioner map for constraint problems
  NOX::NLN::CONSTRAINT::PrecInterfaceMap iconstr_prec;
  STR::NLN::CreateConstraintPreconditioner(iconstr_prec, Integrator(), soltypes);

  // create object to scale linear system
  Teuchos::RCP<NOX::Epetra::Scaling> iscale = Teuchos::null;
  STR::NLN::CreateScaling(iscale, DataSDyn(), DataGlobalState());

  Teuchos::RCP<NOX::NLN::GlobalData> nox_nln_global_data_ptr = Teuchos::rcp(
      new NOX::NLN::GlobalData(DataGlobalState().GetComm(), p_nox, linsolvers, loca_interface_ptr,
          loca_interface_ptr, opttype, iconstr, loca_interface_ptr, iconstr_prec, iscale));

  // ---------------------------------------------------------------------------
  // get initial solution vector and jacobian
  // ---------------------------------------------------------------------------
  Teuchos::RCP<NOX::Epetra::Vector> soln_ptr = Teuchos::rcp(new NOX::Epetra::Vector(
      DataGlobalState().GetMutableDisNp(), NOX::Epetra::Vector::CreateView));
  Teuchos::RCP<LINALG::SparseOperator> jac_ptr = DataGlobalState().GetMutableJacobian();

  // ---------------------------------------------------------------------------
  // Build linear system, loca group, outer/inner/loca status tests
  // ---------------------------------------------------------------------------
  LOCA::NLN::Problem loca_problem = LOCA::NLN::Problem(
      nox_nln_global_data_ptr, loca_global_data_ptr, soln_ptr, jac_ptr, loca_param_vec_ptr);

  Teuchos::RCP<NOX::Epetra::LinearSystem> linsys = loca_problem.CreateLinearSystem();
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> loca_grp_ptr =
      Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::AbstractGroup>(
          loca_problem.CreateGroup(linsys));
  if (loca_grp_ptr.is_null()) dserror("Dynamic cast failed!");

  Teuchos::RCP<NOX::StatusTest::Generic> nox_ostatus = Teuchos::null;
  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> nox_istatus = Teuchos::null;
  Teuchos::RCP<LOCA::StatusTest::Abstract> loca_status = Teuchos::null;
  loca_problem.CreateStatusTests(nox_ostatus, nox_istatus, loca_status);

  // ---------------------------------------------------------------------------
  // Build the stepper object
  // ---------------------------------------------------------------------------
  loca_stepper_ptr_ = Teuchos::rcp(new LOCA::NLN::Stepper(loca_global_data_ptr, loca_grp_ptr,
      loca_status, nox_ostatus, nox_istatus, nox_nln_global_data_ptr, p_loca_nox_ptr_));

  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::CheckForNecessaryParameters() const
{
  CheckInit();

  if (p_loca_nox_ptr_->get<double>("Initial Value") == -1.0e+12)
    dserror(
        "The initial value for the continuation parameter has to be "
        "supplied! Adapt your input file accordingly!");
  if (p_loca_nox_ptr_->get<double>("Max Value") == -1.0e+12)
    dserror(
        "The maximum value for the continuation parameter has to be "
        "supplied. Adapt your input file accordingly!");
  if (p_loca_nox_ptr_->get<double>("Min Value") == -1.0e+12)
    dserror(
        "The minimum value for the continuation parameter has to be "
        "supplied. Adapt your input file accordingly!");
  if (p_loca_nox_ptr_->get<double>("Min Value") >= p_loca_nox_ptr_->get<double>("Max Value"))
    dserror(
        "The maximal value for the continuation parameter has to be larger "
        "than the minimal value. Check your input file!");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::CombineParameterLists()
{
  CheckInit();

  Teuchos::ParameterList& p_nox = DataSDyn().GetMutableNoxParams();
  Teuchos::ParameterList& p_loca = DataSDyn().GetMutableLocaParams();

  // The sublists are references to the two separated lists!
  p_loca_nox_ptr_ = Teuchos::rcp(new Teuchos::ParameterList("LOCA and NOX"));
  p_loca_nox_ptr_->sublist("NOX") = p_nox;
  p_loca_nox_ptr_->sublist("LOCA") = p_loca;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LOCA::ParameterVector> STR::TIMINT::LOCAContinuation::CreateParameterVector(
    const Teuchos::ParameterList& p_loca)
{
  CheckInit();

  Teuchos::RCP<LOCA::ParameterVector> loca_param_vec_ptr =
      Teuchos::rcp(new LOCA::ParameterVector());

  const double& ini_fac = p_loca.get<double>("Initial Value");
  loca_param_vec_ptr->addParameter("force_fac", ini_fac);
  loca_param_vec_ptr->addParameter("displ_fac", ini_fac);

  return loca_param_vec_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::LOCAContinuation::Integrate()
{
  CheckInitSetup();
  loca_stepper_status_ = loca_stepper_ptr_->run();
  return loca_stepper_status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::LOCAContinuation::IntegrateStep()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::SetState(const Teuchos::RCP<Epetra_Vector>& x)
{
  dserror("Not supported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::LOCAContinuation::Solve()
{
  CheckInitSetup();
  dserror("Currently unsupported!");
  return INPAR::STR::conv_success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Abstract::Group& STR::TIMINT::LOCAContinuation::GetSolutionGroup() const
{
  CheckInitSetup();
  Teuchos::RCP<const NOX::Abstract::Group> solgrp_ptr = loca_stepper_ptr_->getSolutionGroup();

  if (solgrp_ptr.is_null()) dserror("The solution group is not initialized!");

  return *solgrp_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> STR::TIMINT::LOCAContinuation::SolutionGroupPtr()
{
  CheckInitSetup();
  Teuchos::RCP<const NOX::Abstract::Group> solgrp_ptr = loca_stepper_ptr_->getSolutionGroup();

  Teuchos::RCP<NOX::Abstract::Group> mutable_solgrp_ptr =
      Teuchos::rcp_const_cast<NOX::Abstract::Group>(solgrp_ptr);

  return mutable_solgrp_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::PrepareTimeStep()
{
  dserror(
      "The LOCA integration is currently not supposed to be used without the"
      "LOCA::Stepper object. For coupled problems use the implicit integrator.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::PreparePartitionStep()
{
  dserror(
      "The LOCA integration is currently not supposed to be used without the"
      "LOCA::Stepper object. For coupled problems use the implicit integrator.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::UpdateStateIncrementally(
    Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  dserror(
      "The LOCA integration is currently not supporting coupled monolithic"
      "problems, thus UpdateStateIncrementally() can not be called. Use the implicit integrator "
      "instead.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  dserror(
      "The LOCA integration is currently not supporting coupled monolithic"
      "problems, thus Evaluate() can not be called. Use the implicit integrator instead.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LOCAContinuation::Evaluate()
{
  dserror(
      "The LOCA integration is currently not supporting coupled monolithic"
      "problems, thus Evaluate() can not be called. Use the implicit integrator instead.");
}
