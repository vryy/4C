/*----------------------------------------------------------------------------*/
/*! \file

\brief Solver for fluid field on a moving ALE mesh


\level 1
*/
/*----------------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_ale.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_dirichletneumann_volcoupl.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::FluidAle::FluidAle(const Teuchos::ParameterList& prbdyn, std::string condname)
    : timeparams_(prbdyn)
{
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(
      prbdyn, GLOBAL::Problem::Instance()->FluidDynamicParams(), "fluid", true, false));
  fluid_ = fluid->FluidField();
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
      new ADAPTER::AleBaseAlgorithm(prbdyn, GLOBAL::Problem::Instance()->GetDis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFluidWrapper>(ale->AleField(), true);

  if (ale_ == Teuchos::null) FOUR_C_THROW("Failed to cast to problem-specific ALE-wrapper");

  const int ndim = GLOBAL::Problem::Instance()->NDim();

  // default parameters for coupling
  double tolerance = 1.e-3;
  int nds_master = 0;
  int nds_slave = 0;

  // set nds_master = 2 in case of HDG discretization
  // (nds = 0 used for trace values, nds = 1 used for interior values)
  if (GLOBAL::Problem::Instance()->spatial_approximation_type() == CORE::FE::ShapeFunctionType::hdg)
  {
    nds_master = 2;
  }

  // check for matching fluid and ale meshes (==true in default case)
  if (CORE::UTILS::IntegralValue<bool>(
          GLOBAL::Problem::Instance()->FSIDynamicParams(), "MATCHGRID_FLUIDALE"))
  {
    // the fluid-ale coupling matches
    const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

    /* Setup coupling adapter
     *
     * Since ALE has been cloned form fluid discretization, nodes reside at the
     * exact same location. Thus, we specify a very tight tolerance for the
     * octree search.
     */
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupfa_matching =
        Teuchos::rcp(new CORE::ADAPTER::Coupling());
    coupfa_matching->SetupCoupling(*FluidField()->Discretization(), *AleField()->Discretization(),
        *fluidnodemap, *alenodemap, ndim,
        CORE::UTILS::IntegralValue<bool>(
            GLOBAL::Problem::Instance()->FSIDynamicParams(), "MATCHALL"),
        tolerance, nds_master, nds_slave);
    coupfa_ = coupfa_matching;
  }
  else
  {
    // non matching volume meshes of fluid and ale
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> coupfa_volmortar =
        Teuchos::rcp(new CORE::ADAPTER::MortarVolCoupl());

    // couple displacement dofs of ale and velocity dofs of fluid

    // projection ale -> fluid : all ndim dofs (displacements)
    std::vector<int> coupleddof12 = std::vector<int>(ndim, 1);

    // projection fluid -> ale : ndim dofs (only velocity, no pressure)
    std::vector<int> coupleddof21 = std::vector<int>(ndim + 1, 1);
    // unmark pressure dof
    coupleddof21[ndim] = 0;

    // define dof sets to be coupled for both projections
    std::pair<int, int> dofsets12(0, 0);
    std::pair<int, int> dofsets21(0, 0);

    // initialize coupling adapter
    coupfa_volmortar->Init(ndim, FluidField()->Discretization(),
        AleField()->write_access_discretization(), &coupleddof12, &coupleddof21, &dofsets12,
        &dofsets21, Teuchos::null, false);

    // setup coupling adapter
    coupfa_volmortar->Setup(GLOBAL::Problem::Instance()->VolmortarParams());

    // set pointer to coupling adapter
    coupfa_ = coupfa_volmortar;
  }

  // Apply initial ALE mesh displacement
  if (CORE::UTILS::IntegralValue<INPAR::ALE::InitialDisp>(
          GLOBAL::Problem::Instance()->AleDynamicParams(), "INITIALDISP") !=
      INPAR::ALE::initdisp_zero_disp)
  {
    FluidField()->SetMeshMap(coupfa_->MasterDofMap(), nds_master);
    Teuchos::RCP<Epetra_Vector> initfluiddisp = AleToFluidField(AleField()->Dispn());
    FluidField()->apply_initial_mesh_displacement(initfluiddisp);
  }

  // initializing the fluid is done later as for xfluids the first cut is done
  // there (coupfa_ cannot be build anymore!!!)
  FluidField()->Init();
  fluid->SetInitialFlowField(
      GLOBAL::Problem::Instance()->FluidDynamicParams());  // call from base algorithm


  if (CORE::UTILS::IntegralValue<bool>(
          GLOBAL::Problem::Instance()->FSIDynamicParams(), "MATCHGRID_STRUCTALE"))
  {
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfa = Teuchos::rcp(new CORE::ADAPTER::Coupling());
    icoupfa->setup_condition_coupling(*FluidField()->Discretization(),
        FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
        AleField()->Interface()->FSICondMap(), condname, ndim, true, nds_master, nds_slave);
    icoupfa_ = icoupfa;
  }
  else
  {
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> icoupfa =
        Teuchos::rcp(new CORE::ADAPTER::MortarVolCoupl());

    // couple displacement dofs of ale and velocity dofs of fluid

    // projection ale -> fluid : all ndim dofs (displacements)
    std::vector<int> coupleddof12 = std::vector<int>(ndim, 1);

    // projection fluid -> ale : ndim dofs (only velocity, no pressure)
    std::vector<int> coupleddof21 = std::vector<int>(ndim + 1, 1);
    // unmark pressure dof
    coupleddof21[ndim] = 0;

    // define dof sets to be coupled for both projections
    std::pair<int, int> dofsets12(0, 0);
    std::pair<int, int> dofsets21(0, 0);

    icoupfa->Init(ndim, GLOBAL::Problem::Instance()->GetDis("fluid"),
        GLOBAL::Problem::Instance()->GetDis("ale"), &coupleddof12, &coupleddof21, &dofsets12,
        &dofsets21, Teuchos::null, false);

    icoupfa->Setup(GLOBAL::Problem::Instance()->VolmortarParams());

    icoupfa_ = icoupfa;
  }

  fscoupfa_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  fscoupfa_->setup_condition_coupling(*FluidField()->Discretization(),
      FluidField()->Interface()->FSCondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSCondMap(), "FREESURFCoupling", ndim, true, nds_master, nds_slave);

  aucoupfa_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  aucoupfa_->setup_condition_coupling(*FluidField()->Discretization(),
      FluidField()->Interface()->AUCondMap(), *AleField()->Discretization(),
      AleField()->Interface()->AUCondMap(), "ALEUPDATECoupling", ndim, true, nds_master, nds_slave);

  FluidField()->SetMeshMap(coupfa_->MasterDofMap(), nds_master);

  // the ale matrix might be build just once
  AleField()->CreateSystemMatrix();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidAle::Discretization()
{
  return FluidField()->Discretization();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::PrepareTimeStep()
{
  FluidField()->PrepareTimeStep();
  AleField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::Update()
{
  FluidField()->Update();
  AleField()->Update();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::Output()
{
  FluidField()->StatisticsAndOutput();

  // Note: We want to write the fsi interface tractions in order to restart
  // monolithically from an partitioned fsi scheme (e.g. fsi prestress simulation).
  // TODO (Thon): this is not the nice way, but fluid-ale and xfem problems may have now FSI
  // interface, so we can not do this in general :(
  if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::fsi)
  {
    // we want to be able to restart monolithically from an partitioned fsi scheme
    const int uprestart = timeparams_.get<int>("RESTARTEVRY");
    const int upres = timeparams_.get<int>("RESULTSEVRY");

    if ((uprestart != 0 && FluidField()->Step() % uprestart == 0) ||
        FluidField()->Step() % upres == 0)
    {
      Teuchos::RCP<Epetra_Vector> lambda = FluidField()->extract_interface_forces();
      Teuchos::RCP<Epetra_Vector> lambdafull =
          FluidField()->Interface()->InsertFSICondVector(lambda);
      FluidField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
    }
  }

  AleField()->Output();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double ADAPTER::FluidAle::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  AleField()->ReadRestart(step);
  return FluidField()->Time();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::NonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    AleField()->apply_interface_displacements(FluidToAle(idisp));
    FluidField()->apply_interface_velocities(ivel);
  }

  // Update the ale update part
  if (FluidField()->Interface()->AUCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField()->Dispnp();
    Teuchos::RCP<Epetra_Vector> audispnp = FluidField()->Interface()->ExtractAUCondVector(dispnp);
    AleField()->apply_ale_update_displacements(aucoupfa_->MasterToSlave(audispnp));
  }

  // Update the free-surface part
  if (FluidField()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField()->Interface()->ExtractFSCondVector(dispnp);
    AleField()->apply_free_surface_displacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  AleField()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField()->Dispnp());
  FluidField()->apply_mesh_displacement(fluiddisp);
  FluidField()->Solve();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::nonlinear_solve_vol_coupl(Teuchos::RCP<Epetra_Vector> idisp,
    Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<FSI::InterfaceCorrector> icorrector)
{
  if (idisp != Teuchos::null)
  {
    AleField()->apply_interface_displacements(AleField()->Interface()->ExtractFSICondVector(idisp));
    FluidField()->apply_interface_velocities(ivel);
  }

  // Update the ale update part
  if (FluidField()->Interface()->AUCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField()->Dispnp();
    Teuchos::RCP<Epetra_Vector> audispnp = FluidField()->Interface()->ExtractAUCondVector(dispnp);
    AleField()->apply_ale_update_displacements(aucoupfa_->MasterToSlave(audispnp));
  }

  // Update the free-surface part
  if (FluidField()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField()->Interface()->ExtractFSCondVector(dispnp);
    AleField()->apply_free_surface_displacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  AleField()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField()->Dispnp());

  icorrector->correct_interface_displacements(fluiddisp, FluidField()->Interface());
  FluidField()->apply_mesh_displacement(fluiddisp);
  FluidField()->Solve();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::apply_interface_values(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    AleField()->apply_interface_displacements(FluidToAle(idisp));
    FluidField()->apply_interface_velocities(ivel);
  }

  if (FluidField()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField()->Interface()->ExtractFSCondVector(dispnp);
    AleField()->apply_free_surface_displacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField()->Dispnp());
  FluidField()->apply_mesh_displacement(fluiddisp);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  // Here we have a mesh position independent of the
  // given trial vector, but still the grid velocity depends on the
  // trial vector only.

  // grid velocity
  AleField()->apply_interface_displacements(FluidToAle(idisp));

  AleField()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField()->Dispnp());
  fluiddisp->Scale(1. / dt);

  FluidField()->ApplyMeshVelocity(fluiddisp);

  // grid position is done inside RelaxationSolve

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return FluidField()->RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::extract_interface_forces()
{
  return FluidField()->extract_interface_forces();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::extract_interface_velnp()
{
  return FluidField()->extract_interface_velnp();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::extract_interface_veln()
{
  return FluidField()->extract_interface_veln();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::integrate_interface_shape()
{
  return FluidField()->integrate_interface_shape();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::UTILS::ResultTest> ADAPTER::FluidAle::CreateFieldTest()
{
  return FluidField()->CreateFieldTest();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::AleToFluidField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::FluidToAle(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

FOUR_C_NAMESPACE_CLOSE
