/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi. Can only be used in conjunction with FLD::FluidImplicitTimeInt

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_fsi.hpp"

#include "4C_adapter_fld_fluid.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::FluidFSI::FluidFSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Core::FE::Discretization> dis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidWrapper(fluid),
      dis_(dis),
      params_(params),
      output_(output),
      dirichletcond_(dirichletcond),
      interface_(Teuchos::rcp(new FLD::UTILS::MapExtractor())),
      meshmap_(Teuchos::rcp(new Core::LinAlg::MapExtractor())),
      locerrvelnp_(Teuchos::null),
      auxintegrator_(Inpar::FSI::timada_fld_none),
      numfsidbcdofs_(0),
      methodadapt_(ada_none)
{
  // make sure
  if (fluid_ == Teuchos::null) FOUR_C_THROW("Failed to create the underlying fluid adapter");
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::init()
{
  // call base class init
  FluidWrapper::init();

  // cast fluid to fluidimplicit
  if (fluidimpl_ == Teuchos::null)
    fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_);

  if (fluidimpl_ == Teuchos::null)
    FOUR_C_THROW("Failed to cast Adapter::Fluid to FLD::FluidImplicitTimeInt.");

  // default dofset for coupling
  int nds_master = 0;

  // set nds_master = 2 in case of HDG discretization
  // (nds = 0 used for trace values, nds = 1 used for interior values)
  if (Global::Problem::instance()->spatial_approximation_type() == Core::FE::ShapeFunctionType::hdg)
  {
    nds_master = 2;
  }

  // create fluid map extractor
  setup_interface(nds_master);

  fluidimpl_->set_surface_splitter(&(*interface_));

  // create map of inner velocity dof (no FSI or Dirichlet conditions)
  build_inner_vel_map();

  if (dirichletcond_)
  {
    // mark all interface velocities as dirichlet values
    fluidimpl_->add_dirich_cond(interface()->fsi_cond_map());
  }

  interfaceforcen_ = Teuchos::rcp(new Epetra_Vector(*(interface()->fsi_cond_map())));

  // time step size adaptivity in monolithic FSI
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const bool timeadapton =
      Core::UTILS::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");
  if (timeadapton)
  {
    // extract the type of auxiliary integrator from the input parameter list
    auxintegrator_ = Core::UTILS::IntegralValue<Inpar::FSI::FluidMethod>(
        fsidyn.sublist("TIMEADAPTIVITY"), "AUXINTEGRATORFLUID");

    if (auxintegrator_ != Inpar::FSI::timada_fld_none)
    {
      // determine type of adaptivity
      if (aux_method_order_of_accuracy() > fluidimpl_->method_order_of_accuracy())
        methodadapt_ = ada_upward;
      else if (aux_method_order_of_accuracy() < fluidimpl_->method_order_of_accuracy())
        methodadapt_ = ada_downward;
      else
        methodadapt_ = ada_orderequal;
    }

    //----------------------------------------------------------------------------
    // Handling of Dirichlet BCs in error estimation
    //----------------------------------------------------------------------------
    // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
    // and are located at the FSI interface.
    std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
    intersectionmaps.push_back(get_dbc_map_extractor()->cond_map());
    intersectionmaps.push_back(interface()->fsi_cond_map());
    Teuchos::RCP<Epetra_Map> intersectionmap =
        Core::LinAlg::MultiMapExtractor::intersect_maps(intersectionmaps);

    // store number of interface DOFs subject to Dirichlet BCs on structure and fluid side of the
    // interface
    numfsidbcdofs_ = intersectionmap->NumGlobalElements();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::FluidFSI::dof_row_map() { return dof_row_map(0); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::FluidFSI::dof_row_map(unsigned nds)
{
  const Epetra_Map* dofrowmap = dis_->dof_row_map(nds);
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::FluidFSI::time_scaling() const
{
  if (params_->get<bool>("interface second order"))
    return 2. / dt();
  else
    return 1. / dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::update()
{
  if (Global::Problem::instance()->spatial_approximation_type() !=
      Core::FE::ShapeFunctionType::hdg)  // TODO als fix this!
  {
    Teuchos::RCP<Epetra_Vector> interfaceforcem =
        interface()->extract_fsi_cond_vector(true_residual());

    interfaceforcen_ = fluidimpl_->extrapolate_end_point(interfaceforcen_, interfaceforcem);
  }

  FluidWrapper::update();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::relaxation_solve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = discretization()->dof_row_map();
  Teuchos::RCP<Epetra_Vector> relax = Core::LinAlg::CreateVector(*dofrowmap, true);
  interface()->insert_fsi_cond_vector(ivel, relax);
  fluidimpl_->linear_relaxation_solve(relax);
  return extract_interface_forces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Adapter::FluidFSI::inner_velocity_row_map() { return innervelmap_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::extract_interface_forces()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem =
      interface()->extract_fsi_cond_vector(true_residual());

  return fluidimpl_->extrapolate_end_point(interfaceforcen_, interfaceforcem);
}

/*----------------------------------------------------------------------*
 | Return interface velocity at new time level n+1                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::extract_interface_velnp()
{
  return interface()->extract_fsi_cond_vector(velnp());
}


/*----------------------------------------------------------------------*
 | Return interface velocity at old time level n                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::extract_interface_veln()
{
  return interface()->extract_fsi_cond_vector(veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::extract_free_surface_veln()
{
  return interface()->extract_fs_cond_vector(veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::apply_interface_velocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  // apply the interface velocities
  interface()->insert_fsi_cond_vector(ivel, fluidimpl_->write_access_velnp());

  const Teuchos::ParameterList& fsipart =
      Global::Problem::instance()->fsi_dynamic_params().sublist("PARTITIONED SOLVER");
  if (Core::UTILS::IntegralValue<int>(fsipart, "DIVPROJECTION"))
  {
    // project the velocity field into a divergence free subspace
    // (might enhance the linear solver, but we are still not sure.)
    proj_vel_to_div_zero();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::apply_initial_mesh_displacement(
    Teuchos::RCP<const Epetra_Vector> initfluiddisp)
{
  // cast fluid to fluidimplicit
  if (fluidimpl_ == Teuchos::null)
    fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_);

  if (fluidimpl_ == Teuchos::null)
    FOUR_C_THROW("Failed to cast Adapter::Fluid to FLD::FluidImplicitTimeInt.");

  meshmap_->insert_cond_vector(initfluiddisp, fluidimpl_->create_dispn());
  meshmap_->insert_cond_vector(initfluiddisp, fluidimpl_->create_dispnp());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::apply_mesh_displacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_->insert_cond_vector(fluiddisp, fluidimpl_->write_access_dispnp());

  // new grid velocity
  fluidimpl_->update_gridv();
}

/*----------------------------------------------------------------------*
 | Update fluid griv velocity via FD approximation           Thon 12/14 |
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::update_gridv()
{
  // new grid velocity via FD approximation
  fluidimpl_->update_gridv();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::apply_mesh_velocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  meshmap_->insert_cond_vector(gridvel, fluidimpl_->write_access_grid_vel());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::set_mesh_map(Teuchos::RCP<const Epetra_Map> mm, const int nds_master)
{
  meshmap_->setup(*dis_->dof_row_map(nds_master), mm,
      Core::LinAlg::SplitMap(*dis_->dof_row_map(nds_master), *mm));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::displacement_to_velocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln_vector =
      interface()->extract_fsi_cond_vector(veln());

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check, whether maps are the same
  if (!fcx->Map().PointSameAs(veln_vector->Map()))
  {
    FOUR_C_THROW("Maps do not match, but they have to.");
  }
#endif

  /*
   * Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double timescale = time_scaling();
  fcx->Update(-timescale * dt(), *veln_vector, timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::velocity_to_displacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln_vector =
      interface()->extract_fsi_cond_vector(veln());

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check, whether maps are the same
  if (!fcx->Map().PointSameAs(veln_vector->Map()))
  {
    FOUR_C_THROW("Maps do not match, but they have to.");
  }
#endif

  /*
   * Delta d(n+1,i+1) = tau * Delta u(n+1,i+1) + dt * u(n)]
   *
   *             / = dt / 2   if interface time integration is second order
   * with tau = |
   *             \ = dt       if interface time integration is first order
   */
  const double tau = 1. / time_scaling();
  fcx->Update(dt(), *veln_vector, tau);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::free_surf_displacement_to_velocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln_vector = interface()->extract_fs_cond_vector(veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = time_scaling();
  fcx->Update(-timescale * dt(), *veln_vector, timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSI::free_surf_velocity_to_displacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln_vector = interface()->extract_fs_cond_vector(veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = 1. / time_scaling();
  fcx->Update(dt(), *veln_vector, timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::integrate_interface_shape()
{
  return interface()->extract_fsi_cond_vector(fluidimpl_->integrate_interface_shape("FSICoupling"));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::use_block_matrix(bool splitmatrix)
{
  Teuchos::RCP<std::set<int>> condelements =
      interface()->conditioned_element_map(*discretization());
  fluidimpl_->use_block_matrix(condelements, *interface(), *interface(), splitmatrix);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Solver> Adapter::FluidFSI::linear_solver()
{
  return FluidWrapper::linear_solver();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::proj_vel_to_div_zero()
{
  // This projection affects also the inner DOFs. Unfortunately, the matrix
  // does not look nice. Hence, the inversion of B^T*B is quite costly and
  // we are not sure yet whether it is worth the effort.

  //   get maps with Dirichlet DOFs and fsi interface DOFs
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcfsimaps;
  dbcfsimaps.push_back(get_dbc_map_extractor()->cond_map());
  dbcfsimaps.push_back(interface()->fsi_cond_map());

  // create a map with all DOFs that have either a Dirichlet boundary condition
  // or are located on the fsi interface
  Teuchos::RCP<Epetra_Map> dbcfsimap = Core::LinAlg::MultiMapExtractor::merge_maps(dbcfsimaps);

  // create an element map with offset
  const int numallele = discretization()->num_global_elements();
  const int mapoffset =
      dbcfsimap->MaxAllGID() + discretization()->element_row_map()->MinAllGID() + 1;
  Teuchos::RCP<Epetra_Map> elemap =
      Teuchos::rcp(new Epetra_Map(numallele, mapoffset, discretization()->get_comm()));

  // create the combination of dbcfsimap and elemap
  std::vector<Teuchos::RCP<const Epetra_Map>> domainmaps;
  domainmaps.push_back(dbcfsimap);
  domainmaps.push_back(elemap);
  Teuchos::RCP<Epetra_Map> domainmap = Core::LinAlg::MultiMapExtractor::merge_maps(domainmaps);

  // build the corresponding map extractor
  Teuchos::RCP<Core::LinAlg::MapExtractor> domainmapex =
      Teuchos::rcp(new Core::LinAlg::MapExtractor(*domainmap, dbcfsimap));

  const int numofrowentries = 82;
  Teuchos::RCP<Core::LinAlg::SparseMatrix> B =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dof_row_map(), numofrowentries, false));

  // define element matrices and vectors
  Core::LinAlg::SerialDenseMatrix elematrix1;
  Core::LinAlg::SerialDenseMatrix elematrix2;
  Core::LinAlg::SerialDenseVector elevector1;
  Core::LinAlg::SerialDenseVector elevector2;
  Core::LinAlg::SerialDenseVector elevector3;

  discretization()->clear_state();
  discretization()->set_state("dispnp", dispnp());

  // loop over all fluid elements
  for (int lid = 0; lid < discretization()->num_my_col_elements(); lid++)
  {
    // get pointer to current element
    Core::Elements::Element* actele = discretization()->l_col_element(lid);

    // get element location vector and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    actele->location_vector(*discretization(), lm, lmowner, lmstride);

    // get dimension of element matrices and vectors
    const int eledim = (int)lm.size();

    // Reshape element matrices and vectors and initialize to zero
    elevector1.size(eledim);

    // set action in order to calculate the integrated divergence operator via an evaluate()-call
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::calc_divop);

    // call the element specific evaluate method
    actele->evaluate(
        params, *discretization(), lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);

    // assembly
    std::vector<int> lmcol(1);
    lmcol[0] = actele->id() + dbcfsimap->MaxAllGID() + 1;
    B->assemble(actele->id(), lmstride, elevector1, lm, lmowner, lmcol);
  }  // end of loop over all fluid elements

  discretization()->clear_state();

  // insert '1's for all DBC and interface DOFs
  for (int i = 0; i < dbcfsimap->NumMyElements(); i++)
  {
    int rowid = dbcfsimap->GID(i);
    int colid = dbcfsimap->GID(i);
    B->assemble(1.0, rowid, colid);
  }

  B->complete(*domainmap, *dof_row_map());

  // Compute the projection operator
  Teuchos::RCP<Core::LinAlg::SparseMatrix> BTB = Core::LinAlg::Multiply(*B, true, *B, false, true);

  Teuchos::RCP<Epetra_Vector> BTvR = Teuchos::rcp(new Epetra_Vector(*domainmap));
  B->multiply(true, *velnp(), *BTvR);
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*dbcfsimap, true));

  domainmapex->insert_cond_vector(zeros, BTvR);

  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*domainmap));

  const Teuchos::ParameterList& fdyn = Global::Problem::instance()->fluid_dynamic_params();
  const int simplersolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  if (simplersolvernumber == (-1))
    FOUR_C_THROW(
        "no simpler solver, that is used to solve this system, defined for fluid pressure problem. "
        "\nPlease set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");

  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::rcp(
      new Core::LinAlg::Solver(Global::Problem::instance()->solver_params(simplersolvernumber),
          discretization()->get_comm(), Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY")));

  if (solver->params().isSublist("ML Parameters"))
  {
    Teuchos::RCP<Epetra_MultiVector> pressure_nullspace =
        Teuchos::rcp(new Epetra_MultiVector(*(dis_->dof_row_map()), 1));
    pressure_nullspace->PutScalar(1.0);

    solver->params().sublist("ML Parameters").set("PDE equations", 1);
    solver->params().sublist("ML Parameters").set("null space: dimension", 1);
    solver->params()
        .sublist("ML Parameters")
        .set("null space: vectors", pressure_nullspace->Values());
    solver->params().sublist("ML Parameters").remove("nullspace", false);  // necessary?
    solver->params()
        .sublist("Michael's secret vault")
        .set<Teuchos::RCP<Epetra_MultiVector>>("pressure nullspace", pressure_nullspace);
  }

  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver->solve(BTB->epetra_operator(), x, BTvR, solver_params);

  Teuchos::RCP<Epetra_Vector> vmod = Teuchos::rcp(new Epetra_Vector(velnp()->Map(), true));
  B->Apply(*x, *vmod);
  write_access_velnp()->Update(-1.0, *vmod, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::reset(bool completeReset, int numsteps, int iter)

{
  FluidWrapper::reset(completeReset, numsteps, iter);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::calculate_error()

{
  fluidimpl_->evaluate_error_compared_to_analytical_sol();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::time_step_auxiliar()
{
  // current state
  Teuchos::RCP<const Epetra_Vector> veln_vector = Teuchos::rcp(new Epetra_Vector(*veln()));
  Teuchos::RCP<const Epetra_Vector> accn_vector = Teuchos::rcp(new Epetra_Vector(*accn()));

  // prepare vector for solution of auxiliary time step
  locerrvelnp_ = Teuchos::rcp(new Epetra_Vector(*fluid_->dof_row_map(), true));

  // ---------------------------------------------------------------------------

  // calculate time step with auxiliary time integrator, i.e. the extrapolated solution
  switch (auxintegrator_)
  {
    case Inpar::FSI::timada_fld_none:
    {
      break;
    }
    case Inpar::FSI::timada_fld_expleuler:
    {
      explicit_euler(*veln_vector, *accn_vector, *locerrvelnp_);

      break;
    }
    case Inpar::FSI::timada_fld_adamsbashforth2:
    {
      if (step() >= 1)  // adams_bashforth2 only if at least second time step
      {
        // Acceleration from previous time step
        Teuchos::RCP<Epetra_Vector> accnm_vector =
            Teuchos::rcp(new Epetra_Vector(*extract_velocity_part(accnm())));

        adams_bashforth2(*veln_vector, *accn_vector, *accnm_vector, *locerrvelnp_);
      }
      else  // explicit_euler as starting algorithm
      {
        explicit_euler(*veln_vector, *accn_vector, *locerrvelnp_);
      }

      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown auxiliary time integration scheme for fluid field.");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::explicit_euler(
    const Epetra_Vector& veln, const Epetra_Vector& accn, Epetra_Vector& velnp) const
{
  // Do a single explicit Euler step
  velnp.Update(1.0, veln, dt(), accn, 0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::adams_bashforth2(const Epetra_Vector& veln, const Epetra_Vector& accn,
    const Epetra_Vector& accnm, Epetra_Vector& velnp) const
{
  // time step sizes of current and previous time step
  const double current_dt = dt();
  const double dto = fluidimpl_->dt_previous();

  // Do a single Adams-Bashforth 2 step
  velnp.Update(1.0, veln, 0.0);
  velnp.Update((2.0 * current_dt * dto + current_dt * current_dt) / (2 * dto), accn,
      -current_dt * current_dt / (2.0 * dto), accnm, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::indicate_error_norms(double& err, double& errcond, double& errother,
    double& errinf, double& errinfcond, double& errinfother)
{
  // compute estimation of local discretization error
  if (methodadapt_ == ada_orderequal)
  {
    const double coeffmarch = fluidimpl_->method_lin_err_coeff_vel();
    const double coeffaux = aux_method_lin_err_coeff_vel();
    locerrvelnp_->Update(-1.0, *velnp(), 1.0);
    locerrvelnp_->Scale(coeffmarch / (coeffaux - coeffmarch));
  }
  else
  {
    // schemes do not have the same order of accuracy
    locerrvelnp_->Update(-1.0, *velnp(), 1.0);
  }

  // set '0' on all pressure DOFs
  auto zeros = Teuchos::rcp(new Epetra_Vector(locerrvelnp_->Map(), true));
  Core::LinAlg::apply_dirichlet_to_system(*locerrvelnp_, *zeros, *pressure_row_map());
  // TODO: Do not misuse apply_dirichlet_to_system()...works for this purpose here: writes zeros
  // into all pressure DoFs

  // set '0' on Dirichlet DOFs
  zeros = Teuchos::rcp(new Epetra_Vector(locerrvelnp_->Map(), true));
  Core::LinAlg::apply_dirichlet_to_system(
      *locerrvelnp_, *zeros, *(get_dbc_map_extractor()->cond_map()));

  // extract the condition part of the full error vector (i.e. only interface velocity DOFs)
  Teuchos::RCP<Epetra_Vector> errorcond =
      Teuchos::rcp(new Epetra_Vector(*interface()->extract_fsi_cond_vector(locerrvelnp_)));

  /* in case of structure split: extract the other part of the full error vector
   * (i.e. interior velocity and all pressure DOFs) */
  Teuchos::RCP<Epetra_Vector> errorother =
      Teuchos::rcp(new Epetra_Vector(*interface()->extract_other_vector(locerrvelnp_)));

  // calculate L2-norms of different subsets of temporal discretization error vector
  // (neglect Dirichlet and pressure DOFs for length scaling)
  err =
      calculate_error_norm(*locerrvelnp_, get_dbc_map_extractor()->cond_map()->NumGlobalElements() +
                                              pressure_row_map()->NumGlobalElements());
  errcond = calculate_error_norm(*errorcond, numfsidbcdofs_);
  errother = calculate_error_norm(
      *errorother, pressure_row_map()->NumGlobalElements() +
                       (get_dbc_map_extractor()->cond_map()->NumGlobalElements() - numfsidbcdofs_));

  // calculate L-inf-norms of temporal discretization errors
  locerrvelnp_->NormInf(&errinf);
  errorcond->NormInf(&errinfcond);
  errorother->NormInf(&errinfother);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidFSI::calculate_wall_shear_stresses()
{
  // get inputs
  Teuchos::RCP<const Epetra_Vector> trueresidual = fluidimpl_->true_residual();
  double dt = fluidimpl_->dt();

  // Get WSSManager
  Teuchos::RCP<FLD::UTILS::StressManager> stressmanager = fluidimpl_->stress_manager();

  // Since the WSS Manager cannot be initialized in the FluidImplicitTimeInt::init()
  // it is not so sure if the WSSManager is jet initialized. So let's be safe here..
  if (stressmanager == Teuchos::null) FOUR_C_THROW("Call of StressManager failed!");
  if (not stressmanager->is_init()) FOUR_C_THROW("StressManager has not been initialized jet!");

  // Call StressManager to calculate WSS from residual
  Teuchos::RCP<Epetra_Vector> wss = stressmanager->get_wall_shear_stresses(trueresidual, dt);

  return wss;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Adapter::FluidFSI::calculate_error_norm(const Epetra_Vector& vec, const int numneglect) const
{
  double norm = 1.0e+12;

  vec.Norm2(&norm);

  if (vec.GlobalLength() - numneglect > 0.0)
    norm /= sqrt((double)(vec.GlobalLength() - numneglect));
  else
    norm = 0.0;

  return norm;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Adapter::FluidFSI::aux_method_order_of_accuracy() const
{
  if (auxintegrator_ == Inpar::FSI::timada_fld_none)
    return 0;
  else if (auxintegrator_ == Inpar::FSI::timada_fld_expleuler)
    return 1;
  else if (auxintegrator_ == Inpar::FSI::timada_fld_adamsbashforth2)
    return 2;
  else
  {
    FOUR_C_THROW("Unknown auxiliary time integration scheme for fluid field.");
    return 0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Adapter::FluidFSI::aux_method_lin_err_coeff_vel() const
{
  if (auxintegrator_ == Inpar::FSI::timada_fld_none)
    return 0.0;
  else if (auxintegrator_ == Inpar::FSI::timada_fld_expleuler)
    return 0.5;
  else if (auxintegrator_ == Inpar::FSI::timada_fld_adamsbashforth2)
  {
    // time step sizes of current and previous time step
    const double dtc = dt();
    const double dto = fluidimpl_->dt_previous();

    // leading error coefficient
    return (2 * dtc + 3 * dto) / (12 * dtc);
  }
  else
  {
    FOUR_C_THROW("Unknown auxiliary time integration scheme for fluid field.");
    return 0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Adapter::FluidFSI::get_tim_ada_err_order() const
{
  if (auxintegrator_ != Inpar::FSI::timada_fld_none)
  {
    if (methodadapt_ == ada_upward)
      return fluidimpl_->method_order_of_accuracy_vel();
    else
      return aux_method_order_of_accuracy();
  }
  else
  {
    FOUR_C_THROW(
        "Cannot return error order for adaptive time integration, since"
        "no auxiliary scheme has been chosen for the fluid field.");
    return 0.0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string Adapter::FluidFSI::get_tim_ada_method_name() const
{
  switch (auxintegrator_)
  {
    case Inpar::FSI::timada_fld_none:
    {
      return "none";
      break;
    }
    case Inpar::FSI::timada_fld_expleuler:
    {
      return "ExplicitEuler";
      break;
    }
    case Inpar::FSI::timada_fld_adamsbashforth2:
    {
      return "AdamsBashforth2";
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown auxiliary time integration scheme for fluid field.");
      return "";
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::setup_interface(const int nds_master)
{
  interface_->setup(*dis_, false, false, nds_master);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::build_inner_vel_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(FluidWrapper::velocity_row_map());
  maps.push_back(interface()->other_map());
  maps.push_back(get_dbc_map_extractor()->other_map());
  innervelmap_ = Core::LinAlg::MultiMapExtractor::intersect_maps(maps);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSI::update_slave_dof(Teuchos::RCP<Epetra_Vector>& f)
{
  fluidimpl_->update_slave_dof(f);
}

FOUR_C_NAMESPACE_CLOSE
