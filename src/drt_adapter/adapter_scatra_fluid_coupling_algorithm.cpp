/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.cpp

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "adapter_scatra_fluid_coupling_algorithm.H"
#include "../drt_fluid/turbulence_statistic_manager.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    bool isale,
    const int disnum
    )
:  AlgorithmBase(comm,prbdyn),
   FluidBaseAlgorithm(prbdyn,isale), // false -> no ALE in fluid algorithm
   ScaTraBaseAlgorithm(prbdyn,isale,disnum), // false -> no ALE in scatra algorithm
   params_(prbdyn)
{
  // transfer the initial convective velocity from initial fluid field to scalar transport field
  // subgrid scales not transferred since they are zero at time t=0.0
  ScaTraField().SetVelocityField(
      FluidField().Velnp(),
      Teuchos::null,
      Teuchos::null,
      FluidField().Discretization()
  );

  // ensure that both single field solvers use the same
  // time integration scheme
  switch (ScaTraField().MethodName())
  {
  case INPAR::SCATRA::timeint_stationary:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  case INPAR::SCATRA::timeint_one_step_theta:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_one_step_theta)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  case INPAR::SCATRA::timeint_bdf2:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_bdf2)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  case INPAR::SCATRA::timeint_gen_alpha:
  {
    if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_gen_alpha and
        FluidField().TimIntScheme() != INPAR::FLUID::timeint_afgenalpha)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  default:
    dserror("Fluid and Scatra time integration schemes do not match");
  }

  // if applicable, provide scatra data to the turbulence statistics
  if (FluidField().TurbulenceStatisticManager() != Teuchos::null)
  {
    // Now, the statistics manager has pointers
    // to ScaTra discretization and result vectors and can access relevant data
    FluidField().TurbulenceStatisticManager()
          ->AddScaTraResults(ScaTraField().Discretization(),ScaTraField().Phinp());
  }

  return;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::~ScaTraFluidCouplingAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidCouplingAlgorithm::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  ScaTraField().ReadRestart(step);
  SetTimeStep(FluidField().Time(),step);

  // read scatra-specific restart data for turbulence statistics
  if (FluidField().TurbulenceStatisticManager() != Teuchos::null)
  {
    IO::DiscretizationReader reader(ScaTraField().Discretization(),step);
    FluidField().TurbulenceStatisticManager()->RestartScaTra(reader,step);
  }

  return;
}


#endif // CCADISCRET
