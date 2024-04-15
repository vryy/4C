/*-----------------------------------------------------------*/
/*! \file

\brief explicit structural time integration


\level 3

*/
/*-----------------------------------------------------------*/


#include "baci_structure_new_timint_explicit.hpp"

#include "baci_solver_nonlin_nox_group.hpp"
#include "baci_solver_nonlin_nox_linearsystem.hpp"
#include "baci_structure_new_nln_solver_factory.hpp"
#include "baci_structure_new_timint_noxinterface.hpp"

#include <NOX_Abstract_Group.H>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Explicit::Explicit() : STR::TIMINT::Base()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Setup()
{
  // safety check
  CheckInit();
  STR::TIMINT::Base::Setup();
  // ---------------------------------------------------------------------------
  // cast the base class integrator
  // ---------------------------------------------------------------------------
  explint_ptr_ = Teuchos::rcp_dynamic_cast<STR::EXPLICIT::Generic>(IntegratorPtr(), true);
  // ---------------------------------------------------------------------------
  // build NOX interface
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::NoxInterface> noxinterface_ptr =
      Teuchos::rcp(new STR::TIMINT::NoxInterface());
  noxinterface_ptr->Init(DataGlobalStatePtr(), explint_ptr_, DBCPtr(), Teuchos::rcp(this, false));
  noxinterface_ptr->Setup();
  // ---------------------------------------------------------------------------
  // build non-linear solver
  // ---------------------------------------------------------------------------
  enum INPAR::STR::NonlinSolTech nlnSolverType = DataSDyn().GetNlnSolverType();
  if (nlnSolverType != INPAR::STR::soltech_singlestep)
  {
    std::cout << "WARNING!!!Nonlinear solver for explicit dynamics is given (in the dat file) as "
              << INPAR::STR::NonlinSolTechString(nlnSolverType)
              << ". This is not compatible. singlestep solver will be selected." << std::endl;
    nlnSolverType = INPAR::STR::soltech_singlestep;
  }
  nlnsolver_ptr_ = STR::NLN::SOLVER::BuildNlnSolver(nlnSolverType);
  nlnsolver_ptr_->Init(DataGlobalStatePtr(), DataSDynPtr(), noxinterface_ptr, explint_ptr_,
      Teuchos::rcp(this, false));
  nlnsolver_ptr_->Setup();
  // set setup flag
  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::PrepareTimeStep()
{
  CheckInitSetup();
  // things that need to be done before Predict
  PrePredict();

  // ToDo prepare contact for new time step
  // PrepareStepContact();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  CheckInitSetup();
  dserror(
      "All monolithically coupled problems work with implicit time "
      "integration schemes. Thus, calling Evaluate() in an explicit scheme "
      "is not possible.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::DetermineStressStrain() { ExplInt().DetermineStressStrain(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  CheckInitSetup();
  dserror(
      "All monolithically coupled problems work with implicit time "
      "integration schemes. Thus, calling Evaluate() in an explicit scheme "
      "is not possible.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Evaluate()
{
  CheckInitSetup();
  ThrowIfStateNotInSyncWithNOXGroup();
  ::NOX::Abstract::Group& grp = NlnSolver().SolutionGroup();

  auto* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  if (grp_ptr == nullptr) dserror("Dynamic cast failed!");

  // you definitely have to evaluate here. You might be called from a coupled
  // problem and the group might not be aware, that a different state than
  // the internally stored displacements may have changed.
  // This is a hack to get NOX to set IsValid to false.
  grp_ptr->setX(grp_ptr->getX());

  // compute the rhs vector and the stiffness matrix
  grp_ptr->computeFandJacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::SetState(const Teuchos::RCP<Epetra_Vector>& x)
{
  dserror(
      "All coupled problems work with implicit time "
      "integration schemes. Thus, calling SetState() in an explicit scheme "
      "is not considered, yet.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::ResetStep()
{
  // calling the base reset
  STR::TIMINT::Base::ResetStep();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Explicit::Solve()
{
  CheckInitSetup();
  IntegrateStep();
  return INPAR::STR::conv_success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::PreparePartitionStep()
{
  // do nothing for explicit time integrators
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::Update(double endtime)
{
  CheckInitSetup();
  dserror("Not implemented. No time adaptivity available for explicit time integration.");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::PrintStep()
{
  CheckInitSetup();

  if (DataGlobalState().GetMyRank() != 0 or GroupId() != 0) return;

  const int stepmax = DataSDyn().GetStepMax();
  const int stepn = DataGlobalState().GetStepN();
  const double timen = DataGlobalState().GetTimeN();
  const double dt = (*DataGlobalState().GetDeltaTime())[0];
  const double wct = DataGlobalState().GetTimer()->totalElapsedTime(true);

  // open outstd::stringstream
  std::ostringstream oss;

  /* Output of the following quantities
   * time   : total simulated time
   * dt     : used time step
   * wct    : wall clock time */
  oss << "Finalised step " << std::setw(1) << stepn;
  oss << " / " << std::setw(1) << stepmax;
  oss << " | time " << std::setw(9) << std::setprecision(3) << std::scientific << timen;
  oss << " | dt " << std::setw(9) << std::setprecision(3) << std::scientific << dt;
  oss << " | wct " << std::setw(8) << std::setprecision(2) << std::scientific << wct;
  oss << "\n--------------------------------------------------------------------------------\n";

  // print to ofile (could be done differently...)
  fprintf(stdout, "%s\n", oss.str().c_str());

  // print it, now
  fflush(stdout);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::StcScale STR::TIMINT::Explicit::GetSTCAlgo()
{
  CheckInitSetup();
  dserror("GetSTCAlgo() has not been tested for explicit time integration.");
  return INPAR::STR::stc_none;
};


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> STR::TIMINT::Explicit::GetSTCMat()
{
  CheckInitSetup();
  dserror("GetSTCMat() has not been tested for explicit time integration.");
  return Teuchos::null;
};


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::Integrate()
{
  dserror(
      "The function is unused since the ADAPTER::StructureTimeLoop "
      "wrapper gives you all the flexibility you need.");
  return 0;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::IntegrateStep()
{
  CheckInitSetup();
  ThrowIfStateNotInSyncWithNOXGroup();
  // reset the non-linear solver
  NlnSolver().Reset();
  // solve the non-linear problem
  NlnSolver().Solve();
  return 0;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::Explicit::InitialGuess()
{
  dserror("InitialGuess() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::TIMINT::Explicit::GetF() const
{
  dserror("RHS() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::TIMINT::Explicit::Freact()
{
  CheckInitSetup();
  dserror("Not implemented!");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> STR::TIMINT::Explicit::SystemMatrix()
{
  dserror("SystemMatrix() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> STR::TIMINT::Explicit::BlockSystemMatrix()
{
  dserror("BlockSystemMatrix() is not available for explicit time integration");
  return Teuchos::null;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Explicit::UseBlockMatrix(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps)
{
  dserror("UseBlockMatrix() is not available for explicit time integration");
}
///@}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::DynamicType STR::TIMINT::Explicit::MethodName() const
{
  return explint_ptr_->MethodName();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::MethodSteps() const { return explint_ptr_->MethodSteps(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::MethodOrderOfAccuracyDis() const
{
  return explint_ptr_->MethodOrderOfAccuracyDis();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Explicit::MethodOrderOfAccuracyVel() const
{
  return explint_ptr_->MethodOrderOfAccuracyVel();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::Explicit::MethodLinErrCoeffDis() const
{
  return explint_ptr_->MethodLinErrCoeffDis();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::Explicit::MethodLinErrCoeffVel() const
{
  return explint_ptr_->MethodLinErrCoeffVel();
}

FOUR_C_NAMESPACE_CLOSE
