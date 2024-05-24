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
  fluid_ = fluid->fluid_field();
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
      new ADAPTER::AleBaseAlgorithm(prbdyn, GLOBAL::Problem::Instance()->GetDis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFluidWrapper>(ale->ale_field(), true);

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
    const Epetra_Map* fluidnodemap = fluid_field()->Discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = ale_field()->Discretization()->NodeRowMap();

    /* Setup coupling adapter
     *
     * Since ALE has been cloned form fluid discretization, nodes reside at the
     * exact same location. Thus, we specify a very tight tolerance for the
     * octree search.
     */
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupfa_matching =
        Teuchos::rcp(new CORE::ADAPTER::Coupling());
    coupfa_matching->setup_coupling(*fluid_field()->Discretization(),
        *ale_field()->Discretization(), *fluidnodemap, *alenodemap, ndim,
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
    coupfa_volmortar->Init(ndim, fluid_field()->Discretization(),
        ale_field()->write_access_discretization(), &coupleddof12, &coupleddof21, &dofsets12,
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
    fluid_field()->SetMeshMap(coupfa_->MasterDofMap(), nds_master);
    Teuchos::RCP<Epetra_Vector> initfluiddisp = AleToFluidField(ale_field()->Dispn());
    fluid_field()->apply_initial_mesh_displacement(initfluiddisp);
  }

  // initializing the fluid is done later as for xfluids the first cut is done
  // there (coupfa_ cannot be build anymore!!!)
  fluid_field()->Init();
  fluid->SetInitialFlowField(
      GLOBAL::Problem::Instance()->FluidDynamicParams());  // call from base algorithm


  if (CORE::UTILS::IntegralValue<bool>(
          GLOBAL::Problem::Instance()->FSIDynamicParams(), "MATCHGRID_STRUCTALE"))
  {
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfa = Teuchos::rcp(new CORE::ADAPTER::Coupling());
    icoupfa->setup_condition_coupling(*fluid_field()->Discretization(),
        fluid_field()->Interface()->FSICondMap(), *ale_field()->Discretization(),
        ale_field()->Interface()->FSICondMap(), condname, ndim, true, nds_master, nds_slave);
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
  fscoupfa_->setup_condition_coupling(*fluid_field()->Discretization(),
      fluid_field()->Interface()->FSCondMap(), *ale_field()->Discretization(),
      ale_field()->Interface()->FSCondMap(), "FREESURFCoupling", ndim, true, nds_master, nds_slave);

  aucoupfa_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  aucoupfa_->setup_condition_coupling(*fluid_field()->Discretization(),
      fluid_field()->Interface()->AUCondMap(), *ale_field()->Discretization(),
      ale_field()->Interface()->AUCondMap(), "ALEUPDATECoupling", ndim, true, nds_master,
      nds_slave);

  fluid_field()->SetMeshMap(coupfa_->MasterDofMap(), nds_master);

  // the ale matrix might be build just once
  ale_field()->CreateSystemMatrix();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidAle::Discretization()
{
  return fluid_field()->Discretization();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::prepare_time_step()
{
  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::Update()
{
  fluid_field()->Update();
  ale_field()->Update();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::Output()
{
  fluid_field()->StatisticsAndOutput();

  // Note: We want to write the fsi interface tractions in order to restart
  // monolithically from an partitioned fsi scheme (e.g. fsi prestress simulation).
  // TODO (Thon): this is not the nice way, but fluid-ale and xfem problems may have now FSI
  // interface, so we can not do this in general :(
  if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::fsi)
  {
    // we want to be able to restart monolithically from an partitioned fsi scheme
    const int uprestart = timeparams_.get<int>("RESTARTEVRY");
    const int upres = timeparams_.get<int>("RESULTSEVRY");

    if ((uprestart != 0 && fluid_field()->Step() % uprestart == 0) ||
        fluid_field()->Step() % upres == 0)
    {
      Teuchos::RCP<Epetra_Vector> lambda = fluid_field()->extract_interface_forces();
      Teuchos::RCP<Epetra_Vector> lambdafull =
          fluid_field()->Interface()->InsertFSICondVector(lambda);
      fluid_field()->DiscWriter()->WriteVector("fsilambda", lambdafull);
    }
  }

  ale_field()->Output();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double ADAPTER::FluidAle::read_restart(int step)
{
  fluid_field()->read_restart(step);
  ale_field()->read_restart(step);
  return fluid_field()->Time();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::NonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    ale_field()->apply_interface_displacements(FluidToAle(idisp));
    fluid_field()->apply_interface_velocities(ivel);
  }

  // Update the ale update part
  if (fluid_field()->Interface()->AUCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->Dispnp();
    Teuchos::RCP<Epetra_Vector> audispnp = fluid_field()->Interface()->ExtractAUCondVector(dispnp);
    ale_field()->apply_ale_update_displacements(aucoupfa_->MasterToSlave(audispnp));
  }

  // Update the free-surface part
  if (fluid_field()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = fluid_field()->Interface()->ExtractFSCondVector(dispnp);
    ale_field()->apply_free_surface_displacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  ale_field()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);
  fluid_field()->Solve();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::nonlinear_solve_vol_coupl(Teuchos::RCP<Epetra_Vector> idisp,
    Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<FSI::InterfaceCorrector> icorrector)
{
  if (idisp != Teuchos::null)
  {
    ale_field()->apply_interface_displacements(
        ale_field()->Interface()->ExtractFSICondVector(idisp));
    fluid_field()->apply_interface_velocities(ivel);
  }

  // Update the ale update part
  if (fluid_field()->Interface()->AUCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->Dispnp();
    Teuchos::RCP<Epetra_Vector> audispnp = fluid_field()->Interface()->ExtractAUCondVector(dispnp);
    ale_field()->apply_ale_update_displacements(aucoupfa_->MasterToSlave(audispnp));
  }

  // Update the free-surface part
  if (fluid_field()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = fluid_field()->Interface()->ExtractFSCondVector(dispnp);
    ale_field()->apply_free_surface_displacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  ale_field()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(ale_field()->Dispnp());

  icorrector->correct_interface_displacements(fluiddisp, fluid_field()->Interface());
  fluid_field()->apply_mesh_displacement(fluiddisp);
  fluid_field()->Solve();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::FluidAle::apply_interface_values(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    ale_field()->apply_interface_displacements(FluidToAle(idisp));
    fluid_field()->apply_interface_velocities(ivel);
  }

  if (fluid_field()->Interface()->FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = fluid_field()->Interface()->ExtractFSCondVector(dispnp);
    ale_field()->apply_free_surface_displacements(fscoupfa_->MasterToSlave(fsdispnp));
  }

  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);
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
  ale_field()->apply_interface_displacements(FluidToAle(idisp));

  ale_field()->Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(ale_field()->Dispnp());
  fluiddisp->Scale(1. / dt);

  fluid_field()->ApplyMeshVelocity(fluiddisp);

  // grid position is done inside RelaxationSolve

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::extract_interface_forces()
{
  return fluid_field()->extract_interface_forces();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::extract_interface_velnp()
{
  return fluid_field()->extract_interface_velnp();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::extract_interface_veln()
{
  return fluid_field()->extract_interface_veln();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::integrate_interface_shape()
{
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::UTILS::ResultTest> ADAPTER::FluidAle::CreateFieldTest()
{
  return fluid_field()->CreateFieldTest();
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
