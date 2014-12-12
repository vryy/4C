/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.cpp

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/


#include "adapter_scatra_fluid_coupling_algorithm.H"
#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_fluid_xfluid/xfluid.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    bool isale,
    const std::string disname,
    const Teuchos::ParameterList& solverparams
    )
:  AlgorithmBase(comm,prbdyn),
   FluidBaseAlgorithm(prbdyn,DRT::Problem::Instance()->FluidDynamicParams(),"fluid",isale,false), // false -> no immediate initialization of fluid time integration
   ScaTraBaseAlgorithm(prbdyn,isale,disname,solverparams), // false -> no ALE in scatra algorithm
   params_(prbdyn)
{
  // check whether fluid and scatra discret still have the same maps
  // they may change due a modified ghosting required, i.e., for particle level-set methods
  if (ScaTraField()->ScaTraType() == INPAR::SCATRA::scatratype_levelset)
  {
    const Epetra_Map* scatraelecolmap = ScaTraField()->Discretization()->ElementColMap();
    if (not scatraelecolmap->PointSameAs(*FluidField()->Discretization()->ElementColMap()))
    {
      if (comm.MyPID()==0)
        std::cout << "----- adaption of fluid ghosting to scatra ghosting ------" << std::endl;

      // adapt fluid ghosting to scatra ghosting
      if (DRT::Problem::Instance()->ProblemType() != prb_combust)
        FluidField()->Discretization()->ExtendedGhosting(*scatraelecolmap,true,true,true,false);
      else
        FluidField()->Discretization()->ExtendedGhosting(*scatraelecolmap,false,false,false,false);
    }
  }

  // initialize fluid time integration scheme
  FluidField()->Init();

  // create initial state for xfluid
  if(DRT::Problem::Instance()->ProblemType() == prb_fluid_xfem_ls)
  {

    Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true);

    xfluid->SetLevelSetField(ScaTraField()->Phinp(), ScaTraField()->Discretization());
    xfluid->CreateInitialState();
  }

  // set also initial field
  if (DRT::Problem::Instance()->ProblemType() != prb_combust)
  {
    SetInitialFlowField(DRT::Problem::Instance()->FluidDynamicParams());
  }

  if (DRT::Problem::Instance()->ProblemType() != prb_combust and DRT::Problem::Instance()->ProblemType() != prb_fluid_xfem_ls)
  {
    // transfer the initial convective velocity from initial fluid field to scalar transport field
    // subgrid scales not transferred since they are zero at time t=0.0
    ScaTraField()->SetVelocityField(
      FluidField()->ConvectiveVel(),
      Teuchos::null,
      Teuchos::null,
      Teuchos::null,
      Teuchos::null,
      FluidField()->Discretization()
    );
  }

  // ensure that both single field solvers use the same
  // time integration scheme
  switch (ScaTraField()->MethodName())
  {
  case INPAR::SCATRA::timeint_stationary:
  {
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary)
      if (comm.MyPID()==0)
        dserror("Fluid and scatra time integration schemes do not match!");
    break;
  }
  case INPAR::SCATRA::timeint_one_step_theta:
  {
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_one_step_theta)
      if (comm.MyPID()==0)
        std::cout << "WARNING: Fluid and scatra time integration schemes do not match!" << std::endl;
    break;
  }
  case INPAR::SCATRA::timeint_bdf2:
  {
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_bdf2)
      if (comm.MyPID()==0)
        std::cout << "WARNING: Fluid and scatra time integration schemes do not match!" << std::endl;
    break;
  }
  case INPAR::SCATRA::timeint_gen_alpha:
  {
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_npgenalpha and
        FluidField()->TimIntScheme() != INPAR::FLUID::timeint_afgenalpha)
      if (comm.MyPID()==0)
        std::cout << "WARNING: Fluid and scatra time integration schemes do not match!" << std::endl;
    break;
  }
  default:
  {
    dserror("Time integration scheme for scalar transport not recognized!");
    break;
  }
  }

  // if applicable, provide scatra data to the turbulence statistics
  if (FluidField()->TurbulenceStatisticManager() != Teuchos::null and ScaTraField()->MethodName()!= INPAR::SCATRA::timeint_stationary)
  {
    // Now, the statistics manager has access to the scatra time integration
    FluidField()->TurbulenceStatisticManager()->AddScaTraField(ScaTraField());
  }

  // if available, allow scatra field to access dynamic Smagorinsky filter
  if (FluidField()->DynSmagFilter() != Teuchos::null)
    ScaTraField()->AccessDynSmagFilter(FluidField()->DynSmagFilter());

  // if available, allow scatra field to access dynamic Vreman
  if (FluidField()->Vreman() != Teuchos::null)
    ScaTraField()->AccessVreman(FluidField()->Vreman());

  return;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidCouplingAlgorithm::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  ScaTraField()->ReadRestart(step);
  SetTimeStep(FluidField()->Time(),step);

  // read scatra-specific restart data for turbulence statistics
  if (FluidField()->TurbulenceStatisticManager() != Teuchos::null)
  {
    IO::DiscretizationReader reader(ScaTraField()->Discretization(),step);
    FluidField()->TurbulenceStatisticManager()->RestartScaTra(reader,step);
  }

  return;
}


