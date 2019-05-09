/*----------------------------------------------------------------------------*/
/**
\brief xcontact level-set one step theta algorithm

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------------*/

#include "xcontact_levelset_timint_ost.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::LEVELSET::TIMINT::OneStepTheta::OneStepTheta(
    const Teuchos::RCP<DRT::Discretization>& actdis, const Teuchos::RCP<LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Teuchos::ParameterList>& sctratimintparams,
    const Teuchos::RCP<Teuchos::ParameterList>& extraparams,
    const Teuchos::RCP<IO::DiscretizationWriter>& output)
    : SCATRA::ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      SCATRA::LevelSetAlgorithm(actdis, solver, params, sctratimintparams, extraparams, output),
      XCONTACT::LEVELSET::Algorithm(actdis, solver, params, sctratimintparams, extraparams, output),
      SCATRA::TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::Init()
{
  // keep the calling order
  SCATRA::TimIntOneStepTheta::Init();
  SCATRA::LevelSetAlgorithm::Init();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::Setup()
{
  // keep the calling order
  SCATRA::TimIntOneStepTheta::Setup();
  XCONTACT::LEVELSET::Algorithm::Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  // do basic restart
  TimIntOneStepTheta::ReadRestart(step, input);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::PrintTimeStepInfo()
{
  if (myrank_) return;

  if (not switchreinit_)
  {
    TimIntOneStepTheta::PrintTimeStepInfo();
    return;
  }

  switch (reinitaction_)
  {
    case INPAR::SCATRA::reinitaction_ellipticeq:
    {
      std::cout << "ELLIPTIC REINITIALIZATION:\n" << std::flush;
      break;
    }
    default:
      /* do nothing */
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::CalcInitialTimeDerivative()
{
  if (not switchreinit_)
  {
    TimIntOneStepTheta::CalcInitialTimeDerivative();
  }
  else
  {
    /* set element parameters with stabilization and artificial diffusivity
     * deactivated */
    SetReinitializationElementParameters(true);

    /* note: time-integration parameter list has not to be overwritten here,
     *       since we rely on incremental solve as already set in
     *       PrepareTimeLoopReinit() */

    // compute time derivative of phi at pseudo-time tau=0
    ScaTraTimIntImpl::CalcInitialTimeDerivative();

    // eventually, undo changes in general parameter list
    SetReinitializationElementParameters();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::SetOldPartOfRighthandside()
{
  if (not switchreinit_)
    SCATRA::TimIntOneStepTheta::SetOldPartOfRighthandside();
  else
    hist_->Update(1.0, *phin_, dtau_ * (1.0 - thetareinit_), *phidtn_, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::UpdateState()
{
  if ((not switchreinit_) and particle_ == Teuchos::null)
  {
    // compute time derivative at time n+1
    ComputeTimeDerivative();

    // after the next command (time shift of solutions) do NOT call
    // ComputeTimeDerivative() anymore within the current time step!!!

    // solution of this step becomes most recent solution of the last step
    phin_->Update(1.0, *phinp_, 0.0);

    // time deriv. of this step becomes most recent time derivative of
    // last step
    phidtn_->Update(1.0, *phidtnp_, 0.0);
  }
  else
  {
    // solution of this step becomes most recent solution of the last step
    phin_->Update(1.0, *phinp_, 0.0);

    // reinitialization is done, reset flag
    if (switchreinit_ == true) switchreinit_ = false;

    // we also have reset the time-integration parameter list, since incremental
    // solver has to be overwritten if used
    SetElementTimeParameter(true);

    // compute time derivative at time n (and n+1)
    ScaTraTimIntImpl::CalcInitialTimeDerivative();

    // reset element time-integration parameters
    SetElementTimeParameter();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::Update(const int num)
{
  // reinitialize
  Reinitialization();

  // update state
  UpdateState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::UpdateReinit()
{
  // TODO: Check if this is really necessary. The whole reinitialization seems to
  //      be one big hack ...
  // compute time derivative at time n+1
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact1 = 1.0 / (thetareinit_ * dtau_);
  const double fact2 = 1.0 - (1.0 / thetareinit_);
  phidtnp_->Update(fact2, *phidtn_, 0.0);
  phidtnp_->Update(fact1, *phinp_, -fact1, *phin_, 1.0);

  // solution of this step becomes most recent solution of the last step
  phin_->Update(1.0, *phinp_, 0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0, *phidtnp_, 0.0);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::TIMINT::OneStepTheta::SetState(const Epetra_Vector& sol_state)
{
  LINALG::AssembleMyVector(0.0, *phinp_, 1.0, sol_state);
}
