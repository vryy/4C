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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    bool isale
    )
:  FSI::AlgorithmBase(comm,prbdyn),
   FluidBaseAlgorithm(prbdyn,isale), // false -> no ALE in fluid algorithm
   ScaTraBaseAlgorithm(prbdyn,isale), // false -> no ALE in scatra algorithm
   params_(prbdyn)
{
  // transfer the initial convective velocity from initial fluid field to scalar transport field
  // subgrid scales not transferred since they are zero at time t=0.0
  ScaTraField().SetVelocityField(
      FluidField().Velnp(),
      Teuchos::null,
      FluidField().Discretization()
  );

  // ensure that both single field solvers use the same
  // time integration scheme
  switch (ScaTraField().MethodName())
  {
  case INPAR::SCATRA::timeint_stationary:
  {
    if (FluidField().TimIntScheme() != timeint_stationary)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  case INPAR::SCATRA::timeint_one_step_theta:
  {
    if (FluidField().TimIntScheme() != timeint_one_step_theta)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  case INPAR::SCATRA::timeint_bdf2:
  {
    if (FluidField().TimIntScheme() != timeint_bdf2)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  case INPAR::SCATRA::timeint_gen_alpha:
  {
    if (FluidField().TimIntScheme() != timeint_gen_alpha and
        FluidField().TimIntScheme() != timeint_afgenalpha)
      dserror("Fluid and Scatra time integration schemes do not match");
    break;
  }
  default:
    dserror("Fluid and Scatra time integration schemes do not match");
  }
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
  ScaTraField().ReadRestart(step);
  FluidField().ReadRestart(step);
  SetTimeStep(FluidField().Time(),step);
  return;
}


#endif // CCADISCRET
