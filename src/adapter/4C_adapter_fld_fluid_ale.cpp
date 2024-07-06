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
Adapter::FluidAle::FluidAle(const Teuchos::ParameterList& prbdyn, std::string condname)
    : timeparams_(prbdyn)
{
  Teuchos::RCP<Adapter::FluidBaseAlgorithm> fluid = Teuchos::rcp(new Adapter::FluidBaseAlgorithm(
      prbdyn, Global::Problem::instance()->fluid_dynamic_params(), "fluid", true, false));
  fluid_ = fluid->fluid_field();
  Teuchos::RCP<Adapter::AleBaseAlgorithm> ale = Teuchos::rcp(
      new Adapter::AleBaseAlgorithm(prbdyn, Global::Problem::instance()->get_dis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<Adapter::AleFluidWrapper>(ale->ale_field(), true);

  if (ale_ == Teuchos::null) FOUR_C_THROW("Failed to cast to problem-specific ALE-wrapper");

  const int ndim = Global::Problem::instance()->n_dim();

  // default parameters for coupling
  double tolerance = 1.e-3;
  int nds_master = 0;
  int nds_slave = 0;

  // set nds_master = 2 in case of HDG discretization
  // (nds = 0 used for trace values, nds = 1 used for interior values)
  if (Global::Problem::instance()->spatial_approximation_type() == Core::FE::ShapeFunctionType::hdg)
  {
    nds_master = 2;
  }

  // check for matching fluid and ale meshes (==true in default case)
  if (Core::UTILS::IntegralValue<bool>(
          Global::Problem::instance()->fsi_dynamic_params(), "MATCHGRID_FLUIDALE"))
  {
    // the fluid-ale coupling matches
    const Epetra_Map* fluidnodemap = fluid_field()->discretization()->node_row_map();
    const Epetra_Map* alenodemap = ale_field()->discretization()->node_row_map();

    /* Setup coupling adapter
     *
     * Since ALE has been cloned form fluid discretization, nodes reside at the
     * exact same location. Thus, we specify a very tight tolerance for the
     * octree search.
     */
    Teuchos::RCP<Core::Adapter::Coupling> coupfa_matching =
        Teuchos::rcp(new Core::Adapter::Coupling());
    coupfa_matching->setup_coupling(*fluid_field()->discretization(),
        *ale_field()->discretization(), *fluidnodemap, *alenodemap, ndim,
        Core::UTILS::IntegralValue<bool>(
            Global::Problem::instance()->fsi_dynamic_params(), "MATCHALL"),
        tolerance, nds_master, nds_slave);
    coupfa_ = coupfa_matching;
  }
  else
  {
    // non matching volume meshes of fluid and ale
    Teuchos::RCP<Core::Adapter::MortarVolCoupl> coupfa_volmortar =
        Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

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
    coupfa_volmortar->init(ndim, fluid_field()->discretization(),
        ale_field()->write_access_discretization(), &coupleddof12, &coupleddof21, &dofsets12,
        &dofsets21, Teuchos::null, false);

    // setup coupling adapter
    coupfa_volmortar->setup(Global::Problem::instance()->volmortar_params(),
        Global::Problem::instance()->cut_general_params());

    // set pointer to coupling adapter
    coupfa_ = coupfa_volmortar;
  }

  // Apply initial ALE mesh displacement
  if (Core::UTILS::IntegralValue<Inpar::ALE::InitialDisp>(
          Global::Problem::instance()->ale_dynamic_params(), "INITIALDISP") !=
      Inpar::ALE::initdisp_zero_disp)
  {
    fluid_field()->set_mesh_map(coupfa_->master_dof_map(), nds_master);
    Teuchos::RCP<Epetra_Vector> initfluiddisp = ale_to_fluid_field(ale_field()->dispn());
    fluid_field()->apply_initial_mesh_displacement(initfluiddisp);
  }

  // initializing the fluid is done later as for xfluids the first cut is done
  // there (coupfa_ cannot be build anymore!!!)
  fluid_field()->init();
  fluid->set_initial_flow_field(
      Global::Problem::instance()->fluid_dynamic_params());  // call from base algorithm


  if (Core::UTILS::IntegralValue<bool>(
          Global::Problem::instance()->fsi_dynamic_params(), "MATCHGRID_STRUCTALE"))
  {
    Teuchos::RCP<Core::Adapter::Coupling> icoupfa = Teuchos::rcp(new Core::Adapter::Coupling());
    icoupfa->setup_condition_coupling(*fluid_field()->discretization(),
        fluid_field()->interface()->fsi_cond_map(), *ale_field()->discretization(),
        ale_field()->interface()->fsi_cond_map(), condname, ndim, true, nds_master, nds_slave);
    icoupfa_ = icoupfa;
  }
  else
  {
    Teuchos::RCP<Core::Adapter::MortarVolCoupl> icoupfa =
        Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

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

    icoupfa->init(ndim, Global::Problem::instance()->get_dis("fluid"),
        Global::Problem::instance()->get_dis("ale"), &coupleddof12, &coupleddof21, &dofsets12,
        &dofsets21, Teuchos::null, false);

    icoupfa->setup(Global::Problem::instance()->volmortar_params(),
        Global::Problem::instance()->cut_general_params());

    icoupfa_ = icoupfa;
  }

  fscoupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  fscoupfa_->setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->interface()->fs_cond_map(), *ale_field()->discretization(),
      ale_field()->interface()->fs_cond_map(), "FREESURFCoupling", ndim, true, nds_master,
      nds_slave);

  aucoupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  aucoupfa_->setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->interface()->au_cond_map(), *ale_field()->discretization(),
      ale_field()->interface()->au_cond_map(), "ALEUPDATECoupling", ndim, true, nds_master,
      nds_slave);

  fluid_field()->set_mesh_map(coupfa_->master_dof_map(), nds_master);

  // the ale matrix might be build just once
  ale_field()->create_system_matrix();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization> Adapter::FluidAle::discretization()
{
  return fluid_field()->discretization();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::FluidAle::prepare_time_step()
{
  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::FluidAle::update()
{
  fluid_field()->update();
  ale_field()->update();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::FluidAle::output()
{
  fluid_field()->statistics_and_output();

  // Note: We want to write the fsi interface tractions in order to restart
  // monolithically from an partitioned fsi scheme (e.g. fsi prestress simulation).
  // TODO (Thon): this is not the nice way, but fluid-ale and xfem problems may have now FSI
  // interface, so we can not do this in general :(
  if (Global::Problem::instance()->get_problem_type() == Core::ProblemType::fsi)
  {
    // we want to be able to restart monolithically from an partitioned fsi scheme
    const int uprestart = timeparams_.get<int>("RESTARTEVRY");
    const int upres = timeparams_.get<int>("RESULTSEVRY");

    if ((uprestart != 0 && fluid_field()->step() % uprestart == 0) ||
        fluid_field()->step() % upres == 0)
    {
      Teuchos::RCP<Epetra_Vector> lambda = fluid_field()->extract_interface_forces();
      Teuchos::RCP<Epetra_Vector> lambdafull =
          fluid_field()->interface()->insert_fsi_cond_vector(lambda);
      fluid_field()->disc_writer()->write_vector("fsilambda", lambdafull);
    }
  }

  ale_field()->output();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double Adapter::FluidAle::read_restart(int step)
{
  fluid_field()->read_restart(step);
  ale_field()->read_restart(step);
  return fluid_field()->time();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::FluidAle::nonlinear_solve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    ale_field()->apply_interface_displacements(fluid_to_ale(idisp));
    fluid_field()->apply_interface_velocities(ivel);
  }

  // Update the ale update part
  if (fluid_field()->interface()->au_cond_relevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->dispnp();
    Teuchos::RCP<Epetra_Vector> audispnp =
        fluid_field()->interface()->extract_au_cond_vector(dispnp);
    ale_field()->apply_ale_update_displacements(aucoupfa_->master_to_slave(audispnp));
  }

  // Update the free-surface part
  if (fluid_field()->interface()->fs_cond_relevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp =
        fluid_field()->interface()->extract_fs_cond_vector(dispnp);
    ale_field()->apply_free_surface_displacements(fscoupfa_->master_to_slave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  ale_field()->solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid_field(ale_field()->dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);
  fluid_field()->solve();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::FluidAle::nonlinear_solve_vol_coupl(Teuchos::RCP<Epetra_Vector> idisp,
    Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<FSI::InterfaceCorrector> icorrector)
{
  if (idisp != Teuchos::null)
  {
    ale_field()->apply_interface_displacements(
        ale_field()->interface()->extract_fsi_cond_vector(idisp));
    fluid_field()->apply_interface_velocities(ivel);
  }

  // Update the ale update part
  if (fluid_field()->interface()->au_cond_relevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->dispnp();
    Teuchos::RCP<Epetra_Vector> audispnp =
        fluid_field()->interface()->extract_au_cond_vector(dispnp);
    ale_field()->apply_ale_update_displacements(aucoupfa_->master_to_slave(audispnp));
  }

  // Update the free-surface part
  if (fluid_field()->interface()->fs_cond_relevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp =
        fluid_field()->interface()->extract_fs_cond_vector(dispnp);
    ale_field()->apply_free_surface_displacements(fscoupfa_->master_to_slave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  ale_field()->solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid_field(ale_field()->dispnp());

  icorrector->correct_interface_displacements(fluiddisp, fluid_field()->interface());
  fluid_field()->apply_mesh_displacement(fluiddisp);
  fluid_field()->solve();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::FluidAle::apply_interface_values(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp != Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    ale_field()->apply_interface_displacements(fluid_to_ale(idisp));
    fluid_field()->apply_interface_velocities(ivel);
  }

  if (fluid_field()->interface()->fs_cond_relevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = fluid_field()->dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp =
        fluid_field()->interface()->extract_fs_cond_vector(dispnp);
    ale_field()->apply_free_surface_displacements(fscoupfa_->master_to_slave(fsdispnp));
  }

  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid_field(ale_field()->dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::relaxation_solve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  // Here we have a mesh position independent of the
  // given trial vector, but still the grid velocity depends on the
  // trial vector only.

  // grid velocity
  ale_field()->apply_interface_displacements(fluid_to_ale(idisp));

  ale_field()->solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid_field(ale_field()->dispnp());
  fluiddisp->Scale(1. / dt);

  fluid_field()->apply_mesh_velocity(fluiddisp);

  // grid position is done inside RelaxationSolve

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->relaxation_solve(idisp);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::extract_interface_forces()
{
  return fluid_field()->extract_interface_forces();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::extract_interface_velnp()
{
  return fluid_field()->extract_interface_velnp();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::extract_interface_veln()
{
  return fluid_field()->extract_interface_veln();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::integrate_interface_shape()
{
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> Adapter::FluidAle::create_field_test()
{
  return fluid_field()->create_field_test();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::ale_to_fluid_field(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->slave_to_master(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::ale_to_fluid_field(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->slave_to_master(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::fluid_to_ale(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->master_to_slave(iv);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidAle::fluid_to_ale(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->master_to_slave(iv);
}

FOUR_C_NAMESPACE_CLOSE
