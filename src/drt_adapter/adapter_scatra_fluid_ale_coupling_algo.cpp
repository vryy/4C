/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and scalar transport equations including deforming meshes
\level 2

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------*/

#include "adapter_scatra_fluid_ale_coupling_algo.H"
#include "adapter_coupling.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidAleCouplingAlgorithm::ScaTraFluidAleCouplingAlgorithm(const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn, const std::string condname,
    const Teuchos::ParameterList& solverparams)
    : ScaTraFluidCouplingAlgorithm(
          comm, prbdyn, true, "scatra", solverparams),  // yes, we need the ALE formulation
      AleBaseAlgorithm(
          prbdyn, DRT::Problem::Instance()->GetDis("ale")),  // construct ale base algorithm as well
      condname_(condname)
{
  // keep constructor empty
  return;
}


/*----------------------------------------------------------------------*
| Setup                                                     rauch 08/16 |
*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidAleCouplingAlgorithm::Init(
    const Teuchos::ParameterList& prbdyn,        ///< parameter list for global problem
    const Teuchos::ParameterList& scatradyn,     ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
    const std::string& disname,                  ///< name of scalar transport discretization
    const bool isale                             ///< ALE flag
)
{
  // call Init() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init(prbdyn, scatradyn, solverparams, disname, isale);

  ale_ = Teuchos::rcp_dynamic_cast<AleFluidWrapper>(AleBaseAlgorithm::AleField(), true);

  return;
}


/*----------------------------------------------------------------------*
| Init                                                      rauch 08/16 |
*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidAleCouplingAlgorithm::Setup()
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  const int ndim = DRT::Problem::Instance()->NDim();

  // set up couplings
  icoupfa_ = Teuchos::rcp(new Coupling());
  icoupfa_->SetupConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), condname_, ndim);

  fscoupfa_ = Teuchos::rcp(new Coupling());
  fscoupfa_->SetupConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->FSCondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSCondMap(), "FREESURFCoupling", ndim);

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

  coupfa_ = Teuchos::rcp(new Coupling());
  coupfa_->SetupCoupling(*FluidField()->Discretization(), *AleField()->Discretization(),
      *fluidnodemap, *alenodemap, ndim);

  FluidField()->SetMeshMap(coupfa_->MasterDofMap());

  // the ale matrix might be build just once!
  AleField()->CreateSystemMatrix(AleField()->Interface());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidAleCouplingAlgorithm::~ScaTraFluidAleCouplingAlgorithm() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidAleCouplingAlgorithm::FluidAleNonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel, bool pseudotransient)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    AleField()->ApplyInterfaceDisplacements(FluidToAle(idisp));
    if (not pseudotransient)
    {
      FluidField()->ApplyInterfaceVelocities(ivel);
    }
  }

  if (FluidField()->Interface()->FSCondRelevant())
  {
    dserror("free surface code in combination with scatra has to be checked");
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField()->Interface()->ExtractFSCondVector(dispnp);
    AleField()->ApplyFreeSurfaceDisplacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  AleField()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField()->WriteAccessDispnp());
  FluidField()->ApplyMeshDisplacement(fluiddisp);

  // no computation of fluid velocities in case only ScaTra and ALE are to compute
  if (not pseudotransient)
  {
    FluidField()->Solve();
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::AleToFluidField(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::AleToFluidField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::FluidToAle(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::FluidToAle(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}
