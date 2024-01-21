/*-----------------------------------------------------------*/
/*! \file

\brief Fluid time integrator for FS3I-AC problems


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_fluid_timint_ac.H"

#include "baci_global_data.H"
#include "baci_io.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     Thon 12/14 |
 *----------------------------------------------------------------------*/
FLD::TimIntAC::TimIntAC(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid)
{
  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::ReadRestart(int step)
{
  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  const bool restartfrompartfsi = INPUT::IntegralValue<int>(fs3idyn, "RESTART_FROM_PART_FSI");

  if (not restartfrompartfsi)  // standard restart
  {
    IO::DiscretizationReader reader(discret_, step);

    reader.ReadVector(trueresidual_, "trueresidual");
  }

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntAC::Output()
{
  FluidImplicitTimeInt::Output();

  // output of solution
  if (uprestart_ > 0 and step_ % uprestart_ == 0)
  {
    output_->WriteVector("trueresidual", trueresidual_);
  }
  return;
}

BACI_NAMESPACE_CLOSE
