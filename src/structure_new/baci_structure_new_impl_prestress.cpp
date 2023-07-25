/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a time integrator for prestressing

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_structure_new_impl_prestress.H"
#include "baci_io.H"
#include "baci_structure_new_model_evaluator.H"
#include "baci_structure_new_timint_basedataglobalstate.H"
#include "baci_lib_prestress_service.H"
#include "baci_io_pstream.H"
#include "baci_io.H"
#include "baci_structure_new_timint_basedatasdyn.H"
#include "baci_linalg_utils_sparse_algebra_create.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::PreStress::PreStress() : absoluteDisplacementNorm_(1e9)
{
  // empty constructor
}


void STR::IMPLICIT::PreStress::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  CheckInitSetup();

  const auto zeros = Teuchos::rcp(new Epetra_Vector(*GlobalState().DofRowMapView(), true));

  // write zero dynamic forces (for dynamic restart after  static prestressing)
  iowriter.WriteVector("finert", zeros);
  iowriter.WriteVector("fvisco", zeros);

  ModelEval().WriteRestart(iowriter, forced_writerestart);
}

void STR::IMPLICIT::PreStress::UpdateStepState()
{
  CheckInitSetup();

  // Compute norm of the displacements
  GlobalState().GetDisNp()->NormInf(&absoluteDisplacementNorm_);

  if (!IsMaterialIterativePrestressConverged())
  {
    // Only update prestress if the material iterative prestress is not converged
    // update model specific variables
    ModelEval().UpdateStepState(0.0);
  }
}

void STR::IMPLICIT::PreStress::UpdateStepElement()
{
  CheckInitSetup();

  if (!IsMaterialIterativePrestressConverged())
  {
    // Only update prestress if the material iterative prestress is not converged
    ModelEval().UpdateStepElement();
  }
}

void STR::IMPLICIT::PreStress::PostUpdate()
{
  // Check for prestressing
  if (::UTILS::PRESTRESS::IsMulfActive(GlobalState().GetTimeN()))

  {
    if (GlobalState().GetMyRank() == 0) IO::cout << "====== Resetting Displacements" << IO::endl;
    // This is a MULF step, hence we do not update the displacements at the end of the
    // timestep. This is achieved by resetting the displacements, velocities and
    // accelerations.
    GlobalState().GetMutableDisN()->PutScalar(0.0);
    GlobalState().GetMutableVelN()->PutScalar(0.0);
    GlobalState().GetMutableAccN()->PutScalar(0.0);
  }
  else if (::UTILS::PRESTRESS::IsMaterialIterativeActive(GlobalState().GetTimeN()))
  {
    // Print prestress status update
    if (GlobalState().GetMyRank() == 0)
    {
      IO::cout << "====== Iterative Prestress Status" << IO::endl;
      IO::cout << "abs-dis-inf-norm:                    " << absoluteDisplacementNorm_ << IO::endl;
    }
  }
}

bool STR::IMPLICIT::PreStress::IsMaterialIterativePrestressConverged() const
{
  return ::UTILS::PRESTRESS::IsMaterialIterative() &&
         GlobalState().GetStepN() >= SDyn().GetPreStressMinimumNumberOfLoadSteps() &&
         absoluteDisplacementNorm_ < SDyn().GetPreStressDisplacementTolerance();
}

bool STR::IMPLICIT::PreStress::EarlyStopping() const
{
  CheckInitSetup();

  if (IsMaterialIterativePrestressConverged())
  {
    if (GlobalState().GetMyRank() == 0)
    {
      IO::cout << "Prestress is converged. Stopping simulation." << IO::endl;
      IO::cout << "abs-dis-inf-norm:                    " << absoluteDisplacementNorm_ << IO::endl;
    }
    return true;
  }

  return false;
}

void STR::IMPLICIT::PreStress::PostTimeLoop()
{
  if (::UTILS::PRESTRESS::IsMaterialIterative())
  {
    if (absoluteDisplacementNorm_ > SDyn().GetPreStressDisplacementTolerance())
    {
      dserror(
          "Prestress algorithm did not converged within the given timesteps. "
          "abs-dis-inf-norm is "
          "%f",
          absoluteDisplacementNorm_);
    }
  }
}
