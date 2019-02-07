/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.cpp

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

\level 1

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/


#include "adapter_scatra_fluid_coupling_algorithm.H"
#include "adapter_coupling_volmortar.H"

#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

// XFEM-specific coupling.
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_levelset/levelset_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn, bool isale, const std::string scatra_disname,
    const Teuchos::ParameterList& solverparams)
    : AlgorithmBase(comm, prbdyn),
      FluidBaseAlgorithm(prbdyn, DRT::Problem::Instance()->FluidDynamicParams(), "fluid", isale,
          false),  // false -> no immediate initialization of fluid time integration
      ScaTraBaseAlgorithm(),
      fieldcoupling_(DRT::INPUT::IntegralValue<INPAR::SCATRA::FieldCoupling>(
          DRT::Problem::Instance()->ScalarTransportDynamicParams(), "FIELDCOUPLING")),
      volcoupl_fluidscatra_(Teuchos::null),
      params_(prbdyn),
      scatra_disname_(scatra_disname),
      issetup_(false),
      isinit_(false)
{
  // keep constructor empty
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidCouplingAlgorithm::Init(
    const Teuchos::ParameterList& prbdyn,        ///< parameter list for global problem
    const Teuchos::ParameterList& scatradyn,     ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
    const std::string& disname,                  ///< name of scalar transport discretization
    const bool isale                             ///< ALE flag
)
{
  SetIsSetup(false);

  ADAPTER::ScaTraBaseAlgorithm::Init(prbdyn, scatradyn, solverparams, disname, isale);

  // perform algorithm specific initialization stuff
  DoAlgorithmSpecificInit();

  // do potential volmortar business
  SetupFieldCoupling("fluid", scatra_disname_);

  SetIsInit(true);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidCouplingAlgorithm::Setup()
{
  CheckIsInit();

  // initialize scatra time integration scheme
  ADAPTER::ScaTraBaseAlgorithm::Setup();

  // initialize fluid time integration scheme
  FluidField()->Init();

  // setup coupling adapter
  if (not volcoupl_fluidscatra_.is_null()) volcoupl_fluidscatra_->Setup();

  // set also initial field
  SetInitialFlowField(DRT::Problem::Instance()->FluidDynamicParams());

  if (DRT::Problem::Instance()->ProblemType() != prb_fluid_xfem_ls)
  {
    // transfer the initial convective velocity from initial fluid field to scalar transport field
    // subgrid scales not transferred since they are zero at time t=0.0
    if (volcoupl_fluidscatra_.is_null())
      ScaTraField()->SetVelocityField(
          FluidField()->ConvectiveVel(), Teuchos::null, Teuchos::null, Teuchos::null, 1);
    else
      ScaTraField()->SetVelocityField(
          volcoupl_fluidscatra_->ApplyVectorMapping21(FluidField()->ConvectiveVel()), Teuchos::null,
          Teuchos::null, Teuchos::null, 1);
  }

  // ensure that both single field solvers use the same
  // time integration scheme
  switch (ScaTraField()->MethodName())
  {
    case INPAR::SCATRA::timeint_stationary:
    {
      if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary)
        if (Comm().MyPID() == 0) dserror("Fluid and scatra time integration schemes do not match!");
      break;
    }
    case INPAR::SCATRA::timeint_one_step_theta:
    {
      if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_one_step_theta)
        if (Comm().MyPID() == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    case INPAR::SCATRA::timeint_bdf2:
    {
      if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_bdf2)
        if (Comm().MyPID() == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    case INPAR::SCATRA::timeint_gen_alpha:
    {
      if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_npgenalpha and
          FluidField()->TimIntScheme() != INPAR::FLUID::timeint_afgenalpha)
        if (Comm().MyPID() == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    default:
    {
      dserror("Time integration scheme for scalar transport not recognized!");
      break;
    }
  }

  // if applicable, provide scatra data to the turbulence statistics
  if (FluidField()->TurbulenceStatisticManager() != Teuchos::null and
      ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
  {
    // Now, the statistics manager has access to the scatra time integration
    FluidField()->TurbulenceStatisticManager()->AddScaTraField(ScaTraField());
  }

  // if available, allow scatra field to access dynamic Smagorinsky filter
  if (FluidField()->DynSmagFilter() != Teuchos::null)
    ScaTraField()->AccessDynSmagFilter(FluidField()->DynSmagFilter());

  // if available, allow scatra field to access dynamic Vreman
  if (FluidField()->Vreman() != Teuchos::null) ScaTraField()->AccessVreman(FluidField()->Vreman());

  // safety check:
  if (volcoupl_fluidscatra_ == Teuchos::null and
      fieldcoupling_ == INPAR::SCATRA::coupling_volmortar)
    dserror("Something went terrible wrong. Sorry about this!");

  SetIsSetup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidCouplingAlgorithm::SetupFieldCoupling(
    const std::string fluid_disname, const std::string scatra_disname)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis(fluid_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

  if (fieldcoupling_ == INPAR::SCATRA::coupling_volmortar)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_fluidscatra_ = Teuchos::rcp(new ADAPTER::MortarVolCoupl());

    // setup projection matrices (use default material strategy)
    volcoupl_fluidscatra_->Init(fluiddis, scatradis, NULL, NULL, NULL, NULL, Teuchos::null, true);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> ADAPTER::ScaTraFluidCouplingAlgorithm::FluidToScatra(
    const Teuchos::RCP<const Epetra_Vector> fluidvector) const
{
  switch (fieldcoupling_)
  {
    case INPAR::SCATRA::coupling_match:
      return fluidvector;
      break;
    case INPAR::SCATRA::coupling_volmortar:
      return volcoupl_fluidscatra_->ApplyVectorMapping21(fluidvector);
      break;
    default:
      dserror("unknown field coupling type");
      return Teuchos::null;
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> ADAPTER::ScaTraFluidCouplingAlgorithm::ScatraToFluid(
    const Teuchos::RCP<const Epetra_Vector> scatravector) const
{
  switch (fieldcoupling_)
  {
    case INPAR::SCATRA::coupling_match:
      return scatravector;
      break;
    case INPAR::SCATRA::coupling_volmortar:
      return volcoupl_fluidscatra_->ApplyVectorMapping12(scatravector);
      break;
    default:
      dserror("unknown field coupling type");
      return Teuchos::null;
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidCouplingAlgorithm::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  ScaTraField()->ReadRestart(step);
  SetTimeStep(FluidField()->Time(), step);

  // read scatra-specific restart data for turbulence statistics
  if (FluidField()->TurbulenceStatisticManager() != Teuchos::null)
  {
    IO::DiscretizationReader reader(ScaTraField()->Discretization(), step);
    FluidField()->TurbulenceStatisticManager()->ReadRestartScaTra(reader, step);
  }

  return;
}
