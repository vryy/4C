/*-----------------------------------------------------------*/
/*! \file
\brief Parallel redistribution for CONTACT::AUG


\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_contact_aug_parallel_distribution_controller.hpp"

#include "4C_contact_aug_strategy.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ParallelDistributionController::setup(CONTACT::ParamsInterface& cparams)
{
  // time measurement (on each processor)
  global_timer_.setComm(&strat_.Comm());
  global_timer_.reset();
  cparams.SetTimer(&global_timer_, 0);

  acttype_ = cparams.get_action_type();

  // rebuild the slave element evaluation time vector only if the redistribution
  // has been changed
  if (sele_eval_times_.is_null() or not sele_eval_times_->Map().SameAs(data_.GSeleColMap()))
  {
    sele_eval_times_ = Teuchos::null;
    sele_eval_times_ = Teuchos::rcp(new Epetra_Vector(data_.GSeleColMap()));
  }
  // otherwise set its content to zero
  else
    sele_eval_times_->PutScalar(0.0);

  cparams.set(sele_eval_times_.get(), 0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::ParallelDistributionController::check(CONTACT::ParamsInterface& cparams)
{
  cparams.ClearEntry(Core::Gen::AnyDataContainer::DataType::time_monitor, 0);
  cparams.ClearEntry(Core::Gen::AnyDataContainer::DataType::any, 0);

  global_timer_.write(Core::IO::cout.os(Core::IO::verbose));

  if (data_.GSeleEvalTimesPtr().is_null() or
      not sele_eval_times_->Map().SameAs(data_.GSeleEvalTimesPtr()->Map()))
  {
    data_.GSeleEvalTimesPtr() = Teuchos::rcp(new Epetra_Vector(*sele_eval_times_));
  }

  // average over the last force/stiff and pure force evaluation
  switch (acttype_)
  {
    case Mortar::eval_force_stiff:
    {
      // delete too old information
      data_.unbalance_element_factors().clear();
      data_.unbalance_time_factors().clear();

      // (re)set the slave element evaluation times
      data_.GSeleEvalTimesPtr()->Scale(1.0, *sele_eval_times_);
      break;
    }
    case Mortar::eval_force:
      data_.GSeleEvalTimesPtr()->Update(1.0, *sele_eval_times_, 1.0);
      strat_.spread_global_sele_eval_times_to_interfaces();
      break;
    // in all other cases: return
    default:
      return;
  }

  strat_.check_parallel_distribution(global_timer_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::Aug::ParallelDistributionController::redistribute(
    const Teuchos::RCP<const Epetra_Vector>& dis, Teuchos::RCP<const Epetra_Vector> vel,
    const int nlniter)
{
  if (nlniter % interval_ == 0)
  {
    const bool is_redis = strat_.RedistributeContact(dis, vel);
    if (is_redis) sele_eval_times_ = Teuchos::null;

    return is_redis;
  }

  return false;
}

FOUR_C_NAMESPACE_CLOSE
