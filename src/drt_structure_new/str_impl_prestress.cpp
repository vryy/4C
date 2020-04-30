/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a time integrator for prestressing

\level 3

*/
/*----------------------------------------------------------------------*/

#include "str_impl_prestress.H"
#include "../drt_io/io.H"
#include "str_model_evaluator.H"
#include "str_timint_basedataglobalstate.H"
#include "../drt_lib/prestress_service.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"
#include "str_timint_basedatasdyn.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

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
    // Compute norm of the displacements before resetting them
    GlobalState().GetDisN()->Norm2(&absoluteDisplacementNorm_);

    if (GlobalState().GetMyRank() == 0)
    {
      IO::cout << "====== Iterative Prestress Status" << IO::endl;
      IO::cout << "abs-dis-norm:                    " << absoluteDisplacementNorm_ << IO::endl;
    }
  }
}

bool STR::IMPLICIT::PreStress::EarlyStopping() const
{
  CheckInitSetup();

  if (::UTILS::PRESTRESS::IsMaterialIterative() && GlobalState().GetStepN() > 0)
  {
    if (absoluteDisplacementNorm_ < SDyn().GetPreStressDisplacementTolerance())
    {
      if (::UTILS::PRESTRESS::IsMaterialIterative() && GlobalState().GetMyRank() == 0)
      {
        IO::cout << "Prestress is converged. Stopping simulation." << IO::endl;
        IO::cout << "abs-dis-norm:                    " << absoluteDisplacementNorm_ << IO::endl;
      }
      return true;
    }
  }

  return false;
}