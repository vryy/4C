/*-----------------------------------------------------------*/
/*!
\file str_integrator.cpp

\maintainer Michael Hiermeier

\date Dec 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_integrator.H"
#include "../drt_lib/drt_dserror.H"

#include "str_dbc.H"
#include "str_timint_base.H"
#include "str_timint_noxinterface.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_data.H"
#include "str_model_evaluator_generic.H"
#include "nox_nln_str_linearsystem.H"

#include "../solver_nonlin_nox/nox_nln_aux.H"

#include "../linalg/linalg_sparsematrix.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator::Integrator()
    : isinit_(false),
      issetup_(false),
      modelevaluator_ptr_(Teuchos::null),
      eval_data_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null),
      isequalibriate_initial_state_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& io_ptr,
    const Teuchos::RCP<STR::Dbc>& dbc_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr
    )
{
  issetup_ = false;

  sdyn_ptr_ = sdyn_ptr;
  gstate_ptr_ = gstate_ptr;
  dbc_ptr_ = dbc_ptr;
  timint_ptr_ = timint_ptr;

  // ---------------------------------------------------------------------------
  // build model evaluator
  // ---------------------------------------------------------------------------
  eval_data_ptr_ =
      Teuchos::rcp(new STR::MODELEVALUATOR::Data());
  eval_data_ptr_->Init(sdyn_ptr,io_ptr,timint_ptr,gstate_ptr->GetCommPtr());
  eval_data_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // build model evaluator
  // ---------------------------------------------------------------------------
  modelevaluator_ptr_ =
      Teuchos::rcp(new STR::ModelEvaluator());
  modelevaluator_ptr_->Init(sdyn_ptr->GetModelTypes(),eval_data_ptr_,
      gstate_ptr,io_ptr,Teuchos::rcp(this,false),timint_ptr);
  modelevaluator_ptr_->Setup();

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::EquilibriateInitialState()
{
  CheckInit();

  // temporary right-hand-side
  Teuchos::RCP<Epetra_Vector> rhs_ptr =
      Teuchos::rcp(new Epetra_Vector(*GlobalState().DofRowMapView(),true));
  // wrap the rhs_ptr in a nox_epetra_Vector
  Teuchos::RCP<NOX::Epetra::Vector> nox_rhs_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(rhs_ptr,NOX::Epetra::Vector::CreateView));

  // initialize the structural stiffness matrix
  Teuchos::RCP<LINALG::SparseOperator> jac_ptr =
      GlobalState().GetMutableJacobian();
  Teuchos::RCP<LINALG::SparseOperator> stiff_ptr =
      GlobalState().ExtractDisplBlock(*jac_ptr);

  /* auxiliary vector in order to store accelerations of inhomogeneous
   * Dirichilet-DoFs
   * Meier 2015: This contribution is necessary in order to determine correct
   * initial accelerations in case of inhomogeneous Dirichlet conditions */
  Teuchos::RCP<Epetra_Vector> acc_aux_ptr =
      Teuchos::rcp(new Epetra_Vector(*GlobalState().DofRowMapView(), true));

  // overwrite initial state vectors with Dirichlet BCs
  const double& timen = (*GlobalState().GetMultiTime())[0];
  Teuchos::RCP<Epetra_Vector> disn_ptr = GlobalState().GetMutableDisN();
  Teuchos::RCP<Epetra_Vector> veln_ptr = GlobalState().GetMutableVelN();
  Teuchos::RCP<Epetra_Vector> accn_ptr = GlobalState().GetMutableAccN();
  Dbc().ApplyDirichletBC(timen,disn_ptr,veln_ptr,acc_aux_ptr,false);


  // ---------------------------------------------------------------------------
  // compute the mass matrix
  // ---------------------------------------------------------------------------
  // set the evaluate parameters of the current base class
  ResetEvalParams();
  // !!! evaluate the initial state !!!
  EvalData().SetTotalTime(gstate_ptr_->GetTimeN());

  isequalibriate_initial_state_ = true;
  // evaluate the stiffness only on the structural model evaluator
  ModelEval().Evaluator(INPAR::STR::model_structure).
      ApplyStiff(*disn_ptr,*stiff_ptr);
  // build the entire right-hand-side
  ModelEval().ApplyForce(*disn_ptr,*rhs_ptr);
  isequalibriate_initial_state_ = false;

  /* Meier 2015: Here, we copy the mass matrix in the stiffness block in order to
   * not perform the Dirichlet conditions on the constant mass matrix later on.
   * This is necessary since we need the original mass matrix (without blanked
   * rows) on the Dirichlet DoFs in order to calculate correct reaction
   * forces.*/
  stiff_ptr->Add(*GlobalState().GetMassMatrix(),false,1.0,0.0);
  stiff_ptr->Complete();

  // ---------------------------------------------------------------------------
  // build a NOX::NLN::STR::LinearSystem
  // ---------------------------------------------------------------------------
  // get the structural linear solver
  std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> > str_linsolver;
  str_linsolver[NOX::NLN::sol_structure] =
      TimInt().GetDataSDyn().GetLinSolvers().at(INPAR::STR::model_structure);

  // copy the nox parameter-list
  Teuchos::ParameterList p_nox =
      TimInt().GetDataSDyn().GetNoxParams();
  NOX::NLN::AUX::SetPrintingParameters(p_nox,GlobalState().GetComm());

  // create a copy of the initial displacement vector
  Teuchos::RCP<Epetra_Vector> soln_ptr =
      Teuchos::rcp(new Epetra_Vector(*GlobalState().DofRowMapView(),true));
  // wrap the soln_ptr in a nox_epetra_Vector
  Teuchos::RCP<NOX::Epetra::Vector> nox_soln_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(soln_ptr,NOX::Epetra::Vector::CreateView));

  // Check if we are using a Newton direction
  const std::string& dir_str =
      p_nox.sublist("Direction").get<std::string>("Method");
  if (dir_str != "Newton")
    dserror("The EquilibriateState predictor is currently only working for the "
        "direction-method \"Newton\".");

  // create the linear system
  // printing parameters
  Teuchos::ParameterList& p_print = p_nox.sublist("Printing",true);
  // linear solver parameters
  Teuchos::ParameterList& p_ls    = p_nox.sublist("Direction",true).
        sublist("Newton",true).sublist("Linear Solver",true);

  Teuchos::RCP<NOX::NLN::LinearSystem> linsys_ptr =
      Teuchos::rcp(new NOX::NLN::STR::LinearSystem(p_print,p_ls,
          str_linsolver,Teuchos::null,Teuchos::null,stiff_ptr,*nox_soln_ptr));

  // (re)set the linear solver parameters
  p_ls.set("Number of Nonlinear Iterations",0);
  p_ls.set("Current Time Step",GlobalState().GetTimeNp());
  // ToDo Get the actual tolerance value
  p_ls.set("Wanted Tolerance",1.0e-6);
  // ---------------------------------------------------------------------------
  /* Meier 2015: Due to the Dirichlet conditions applied to the mass matrix, we
   * solely solve for the accelerations at non-Dirichlet DoFs while the
   * resulting accelerations at the Dirichlet-DoFs will be zero. Therefore, the
   * accelerations at DoFs with inhomogeneous Dirichlet conditions will be added
   * below at *). */
  // ---------------------------------------------------------------------------
  // solve the linear system
  linsys_ptr->applyJacobianInverse(p_ls,*nox_rhs_ptr,*nox_soln_ptr);

  // get the solution vector and copy it into the acceleration vector
  accn_ptr->PutScalar(0.0);
  accn_ptr->Update(1.0,nox_soln_ptr->getEpetraVector(),0.0);
  //*) Add contributions of inhomogeneous DBCs
  accn_ptr->Update(1.0,*acc_aux_ptr,1.0);

  /* We need to reset the stiffness matrix because its graph (topology)
   * is not finished yet in case of constraints and possibly other side
   * effects (basically remaining model evaluators). */
  stiff_ptr->Reset();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::DetermineStressStrain()
{
  CheckInitSetup();
  ModelEval().DetermineStressStrain();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::DetermineEnergy()
{
  CheckInitSetup();
  ModelEval().DetermineEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::OutputStepState()
{
  CheckInitSetup();
  ModelEval().OutputStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::GetCondensedUpdateNorm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  CheckInitSetup();

  double myupdatenorm = eval_data_ptr_->GetMyUpdateNorm(qtype);
  const enum NOX::Abstract::Vector::NormType normtype =
      eval_data_ptr_->GetUpdateNormType(qtype);

  return GetCondensedGlobalNorm(qtype,normtype,myupdatenorm);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::GetCondensedPreviousSolNorm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  CheckInitSetup();

  double myprevsolnorm = eval_data_ptr_->GetMyPreviousSolNorm(qtype);
  const enum NOX::Abstract::Vector::NormType normtype =
      eval_data_ptr_->GetUpdateNormType(qtype);

  return GetCondensedGlobalNorm(qtype,normtype,myprevsolnorm);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::GetCondensedSolutionUpdateRMS(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  CheckInitSetup();
  // global relative mean square norm
  double grmsnorm = 0.0;
  // get proc data
  double myrmsnorm = eval_data_ptr_->GetMyRMSNorm(qtype);
  // get total dof number
  int gdofnumber = GetCondensedDofNumber(qtype);
  // sum over all procs
  gstate_ptr_->GetComm().SumAll(&myrmsnorm,&grmsnorm,1);

  // calculate the root mean square and return the value
  return (std::sqrt(grmsnorm)/static_cast<double>(gdofnumber));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::Integrator::GetCondensedDofNumber(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  CheckInitSetup();
  // global dof number of the given quantity
  int gdofnumber = 0.0;
  int mydofnumber = eval_data_ptr_->GetMyDofNumber(qtype);
  gstate_ptr_->GetComm().SumAll(&mydofnumber,&gdofnumber,1);
  return gdofnumber;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::GetCondensedGlobalNorm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const enum NOX::Abstract::Vector::NormType& normtype,
    double& mynorm) const
{
  double gnorm = 0;

  switch (normtype)
  {
    case NOX::Abstract::Vector::OneNorm:
    {
      gstate_ptr_->GetComm().SumAll(&mynorm,&gnorm,1);
      break;
    }
    case NOX::Abstract::Vector::TwoNorm:
    {
      gstate_ptr_->GetComm().SumAll(&mynorm,&gnorm,1);
      gnorm = std::sqrt(gnorm);
      break;
    }
    case NOX::Abstract::Vector::MaxNorm:
    {
      gstate_ptr_->GetComm().MaxAll(&mynorm,&gnorm,1);
      break;
    }
  }
  return gnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ModelEvaluator& STR::Integrator::ModelEval()
{
  CheckInit();
  return *modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::ModelEvaluator& STR::Integrator::ModelEval() const
{
  CheckInit();
  return *modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Data& STR::Integrator::EvalData() const
{
  CheckInit();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Data& STR::Integrator::EvalData()
{
  CheckInit();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataSDyn& STR::Integrator::SDyn()
{
  CheckInit();
  return *sdyn_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataSDyn& STR::Integrator::SDyn() const
{
  CheckInit();
  return *sdyn_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::Integrator::GlobalState() const
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& STR::Integrator::GlobalState()
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Dbc& STR::Integrator::Dbc()
{
  CheckInit();
  return *dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::Dbc& STR::Integrator::Dbc() const
{
  CheckInit();
  return *dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::Base& STR::Integrator::TimInt() const
{
  CheckInit();
  return *timint_ptr_;
}
