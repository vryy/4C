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
