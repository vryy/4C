/*----------------------------------------------------------------------------*/
/*! \file

\brief ALE time integration


\level 1
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale.hpp"

#include "4C_ale_meshsliding.hpp"
#include "4C_ale_resulttest.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_fem_condition_locsys.hpp"
#include "4C_fem_condition_periodic.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ALE::Ale::Ale(Teuchos::RCP<Core::FE::Discretization> actdis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : discret_(actdis),
      solver_(solver),
      params_(params),
      output_(output),
      step_(0),
      numstep_(params_->get<int>("NUMSTEP")),
      time_(0.0),
      maxtime_(params_->get<double>("MAXTIME")),
      dt_(params_->get<double>("TIMESTEP")),
      writerestartevery_(params->get<int>("RESTARTEVRY")),
      writeresultsevery_(params->get<int>("RESULTSEVRY")),
      sysmat_(Teuchos::null),
      residual_(Teuchos::null),
      rhs_(Teuchos::null),
      dispnp_(Teuchos::null),
      dispn_(Teuchos::null),
      disi_(Teuchos::null),
      zeros_(Teuchos::null),
      eledetjac_(Teuchos::null),
      elequality_(Teuchos::null),
      elequalityyesno_(Core::UTILS::IntegralValue<bool>(*params, "ASSESSMESHQUALITY")),
      aletype_(Core::UTILS::IntegralValue<Inpar::ALE::AleDynamic>(*params, "ALE_TYPE")),
      maxiter_(params->get<int>("MAXITER")),
      tolres_(params->get<double>("TOLRES")),
      toldisp_(params->get<double>("TOLDISP")),
      divercont_(Core::UTILS::IntegralValue<Inpar::ALE::DivContAct>(*params, "DIVERCONT")),
      msht_(Core::UTILS::IntegralValue<Inpar::ALE::MeshTying>(*params, "MESHTYING")),
      initialdisp_(Core::UTILS::IntegralValue<Inpar::ALE::InitialDisp>(*params, "INITIALDISP")),
      startfuncno_(params->get<int>("STARTFUNCNO"))
{
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  dispn_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  dispnp_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  disi_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  residual_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  rhs_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  zeros_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  eledetjac_ = Core::LinAlg::CreateVector(*discretization()->element_row_map(), true);
  elequality_ = Core::LinAlg::CreateVector(*discretization()->element_row_map(), true);

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  set_initial_displacement(initialdisp_, startfuncno_);

  setup_dbc_map_ex();

  // ensure that the ALE string was removed from conditions
  {
    Core::Conditions::Condition* cond = discret_->get_condition("ALEDirichlet");
    if (cond) FOUR_C_THROW("Found a ALE Dirichlet condition. Remove ALE string!");
  }

  if (msht_ == Inpar::ALE::meshsliding)
  {
    meshtying_ = Teuchos::rcp(
        new Meshsliding(discret_, *solver_, msht_, Global::Problem::instance()->n_dim(), nullptr));
  }
  else if (msht_ == Inpar::ALE::meshtying)
  {
    meshtying_ = Teuchos::rcp(
        new Meshtying(discret_, *solver_, msht_, Global::Problem::instance()->n_dim(), nullptr));
  }

  // ---------------------------------------------------------------------
  // Create LocSysManager, if needed (used for LocSys-Dirichlet BCs)
  // ---------------------------------------------------------------------
  {
    std::vector<Core::Conditions::Condition*> locsysconditions(0);
    discret_->get_condition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      // Initialize locsys manager
      locsysman_ = Teuchos::rcp(
          new Core::Conditions::LocsysManager(*discret_, Global::Problem::instance()->n_dim()));
    }
  }

  create_system_matrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ALE::Ale::set_initial_displacement(const Inpar::ALE::InitialDisp init, const int startfuncno)
{
  switch (init)
  {
    case Inpar::ALE::initdisp_zero_disp:
    {
      dispn_->PutScalar(0.0);
      dispnp_->PutScalar(0.0);

      break;
    }
    case Inpar::ALE::initdisp_disp_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->num_my_row_nodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->l_row_node(lnodeid);

        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->dof(0, lnode);

        // loop nodal dofs
        for (unsigned int d = 0; d < nodedofset.size(); ++d)
        {
          const int dofgid = nodedofset[d];
          int doflid = dofrowmap->LID(dofgid);

          // evaluate component d of function
          double initialval =
              Global::Problem::instance()
                  ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                  .evaluate(lnode->x().data(), 0, d);

          int err = dispn_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }

      // initialize also the solution vector
      dispnp_->Update(1.0, *dispn_, 0.0);

      break;
    }
    default:
      FOUR_C_THROW("Unknown option for initial displacement: %d", init);
      break;
  }

  return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::create_system_matrix(Teuchos::RCP<const ALE::UTILS::MapExtractor> interface)
{
  if (msht_ != Inpar::ALE::no_meshtying)
  {
    std::vector<int> coupleddof(Global::Problem::instance()->n_dim(), 1);
    sysmat_ = meshtying_->setup(coupleddof, dispnp_);
    meshtying_->dirichlet_on_master(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->cond_map());

    if (interface != Teuchos::null)
    {
      meshtying_->is_multifield(*interface, true);
    }
  }
  else if (interface == Teuchos::null)
  {
    const Epetra_Map* dofrowmap = discret_->dof_row_map();
    sysmat_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap, 81, false, true));
  }
  else
  {
    sysmat_ =
        Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
            *interface, *interface, 81, false, true));
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::evaluate(
    Teuchos::RCP<const Epetra_Vector> stepinc, ALE::UTILS::MapExtractor::AleDBCSetType dbc_type)
{
  // We save the current solution here. This will not change the
  // result of our element call, but the next time somebody asks us we
  // know the displacements.
  //
  // Note: What we get here is the sum of all increments in this time
  // step, not just the latest increment.

  disi_->PutScalar(0.0);

  if (stepinc != Teuchos::null)
  {
    dispnp_->Update(1.0, *stepinc, 1.0, *dispn_, 0.0);
  }

  if (msht_ != Inpar::ALE::no_meshtying)
  {
    meshtying_->msht_split(sysmat_);
  }

  evaluate_elements();
  evaluate_element_quality();

  // prepare meshtying system
  if (msht_ != Inpar::ALE::no_meshtying)
  {
    meshtying_->prepare_meshtying_system(sysmat_, residual_, dispnp_);
    meshtying_->multifield_split(sysmat_);
  }

  // dispnp_ has zeros at the Dirichlet-entries, so we maintain zeros there.
  if (locsys_manager() != Teuchos::null)
  {
    // Transform system matrix and rhs to local coordinate systems
    locsys_manager()->rotate_global_to_local(
        Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_), residual_);

    // When using local systems, a rotated dispnp_ vector needs to be used as dbcval for
    // apply_dirichlet_to_system
    Teuchos::RCP<Epetra_Vector> dispnp_local = Teuchos::rcp(new Epetra_Vector(*(zeros_)));
    locsys_manager()->rotate_global_to_local(dispnp_local);

    if (get_loc_sys_trafo() != Teuchos::null)
    {
      Core::LinAlg::apply_dirichlet_to_system(
          *Core::LinAlg::CastToSparseMatrixAndCheckSuccess(sysmat_), *disi_, *residual_,
          *get_loc_sys_trafo(), *dispnp_local, *(dbcmaps_[dbc_type]->cond_map()));
    }
    else
    {
      Core::LinAlg::apply_dirichlet_to_system(
          *sysmat_, *disi_, *residual_, *dispnp_local, *(dbcmaps_[dbc_type]->cond_map()));
    }
  }
  else
  {
    Core::LinAlg::apply_dirichlet_to_system(
        *sysmat_, *disi_, *residual_, *zeros_, *(dbcmaps_[dbc_type]->cond_map()));
  }

  /* residual_ contains the most recent "mechanical" residual including DBCs.
   * We make this negative and store it in rhs_ for use in Newton-type methods.
   */
  rhs_->Update(-1.0, *residual_, 0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ALE::Ale::solve()
{
  // We need the negative residual here as right hand side of the linear problem
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(*residual_));
  rhs->Scale(-1.0);

  // ToDo (mayr) Why can't we use rhs_ instead of local variable rhs???
  int errorcode = 0;
  if (msht_ == Inpar::ALE::no_meshtying)
  {
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    errorcode = solver_->solve(sysmat_->epetra_operator(), disi_, rhs, solver_params);
  }
  else
    errorcode = meshtying_->solve_meshtying(*solver_, sysmat_, disi_, rhs, dispnp_);
  // calc norm
  disi_->Norm2(&normdisi_);
  normdisi_ /= sqrt(disi_->GlobalLength());

  return errorcode;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::update_iter() { dispnp_->Update(1.0, *disi_, 1.0); }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::update() { dispn_->Update(1.0, *dispnp_, 0.0); }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool ALE::Ale::converged(const int iter)
{
  if (iter == 0) normdisi_ = 0.0;

  bool converged = false;
  // determine norms
  double res_norm;
  residual_->Norm2(&res_norm);
  res_norm /= sqrt(residual_->GlobalLength());
  if (discret_->get_comm().MyPID() == 0)
    std::cout << "ITER: " << iter << "  RES NORM: " << res_norm << " DISP NORM: " << normdisi_
              << std::endl;

  if (res_norm < tolres_ && normdisi_ < toldisp_)
  {
    converged = true;
  }
  return converged;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::evaluate_elements()
{
  sysmat_->zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // set vector values needed by elements
  discret_->clear_state();

  // action for elements
  eleparams.set<std::string>("action", element_action_string(aletype_));
  eleparams.set<bool>("use spatial configuration", update_sys_mat_every_step());

  discret_->set_state("dispnp", dispnp_);

  discret_->evaluate(eleparams, sysmat_, residual_);
  discret_->clear_state();

  sysmat_->complete();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::string ALE::Ale::element_action_string(const enum Inpar::ALE::AleDynamic name)
{
  switch (name)
  {
    case Inpar::ALE::solid:
      return "calc_ale_solid";
      break;
    case Inpar::ALE::solid_linear:
      return "calc_ale_solid_linear";
      break;
    case Inpar::ALE::laplace_material:
      return "calc_ale_laplace_material";
      break;
    case Inpar::ALE::laplace_spatial:
      return "calc_ale_laplace_spatial";
      break;
    case Inpar::ALE::springs_material:
      return "calc_ale_springs_material";
      break;
    case Inpar::ALE::springs_spatial:
      return "calc_ale_springs_spatial";
      break;
    default:
      FOUR_C_THROW("Cannot make std::string for ALE type %d", name);
      return "";
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ALE::Ale::dof_row_map() const
{
  return Teuchos::rcp(discret_->dof_row_map(), false);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> ALE::Ale::system_matrix()
{
  return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> ALE::Ale::block_system_matrix()
{
  return Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ALE::Ale::integrate()
{
  // add eps to prevent stopping one step too early due to memory trash on last digits
  const double eps = 1.0e-12;
  while (step_ < numstep_ and time_ <= maxtime_ + eps)
  {
    prepare_time_step();
    time_step();
    update();
    output();
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::output()
{
  /*  We need ALE output only in case of pure ALE problems. If fluid is present,
   *  the fluid field writes its own displacement field as output.
   *
   *  Though, we might need restart data.
   */

  // Has any output data been written?
  bool datawritten = false;

  // write restart data if necessary
  if (writerestartevery_ != 0 and step_ % writerestartevery_ == 0)
  {
    output_restart(datawritten);
  }

  // write output data if necessary
  if (not datawritten and writeresultsevery_ != 0 and step_ % writeresultsevery_ == 0)
  {
    output_state(datawritten);
  }

  // write domain decomposition for visualization
  if ((step_ == writeresultsevery_ or step_ == 0)) output_->write_element_data(true);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::output_state(bool& datawritten)
{
  // write output data
  output_->new_step(step_, time_);
  output_->write_vector("dispnp", dispnp_);

  if (elequalityyesno_)
  {
    output_->write_vector("det_j", eledetjac_, Core::IO::elementvector);
    output_->write_vector("element_quality", elequality_, Core::IO::elementvector);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::output_restart(bool& datawritten)
{
  // write restart data
  output_->new_step(step_, time_);
  output_->write_vector("dispnp", dispnp_);
  output_->write_vector("dispn", dispn_);

  // restart/output data has been written
  datawritten = true;

  // info dedicated to user's eyes staring at standard out
  // Print restart info only in case of pure ALE problem. Coupled problems
  // print their own restart info.
  if (Global::Problem::instance()->get_problem_type() == Core::ProblemType::ale)
  {
    if (discret_->get_comm().MyPID() == 0)
      Core::IO::cout << "====== Restart written in step " << step_ << Core::IO::endl;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::read_restart(const int step)
{
  Core::IO::DiscretizationReader reader(
      discret_, Global::Problem::instance()->input_control_file(), step);
  time_ = reader.read_double("time");
  step_ = reader.read_int("step");

  reader.read_vector(dispnp_, "dispnp");
  reader.read_vector(dispn_, "dispn");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;

  // Print time step header only in case of pure ALE problem. Coupled problems
  // print their own time step header.
  if (Global::Problem::instance()->get_problem_type() == Core::ProblemType::ale)
    print_time_step_header();

  // Update local coordinate systems (which may be time dependent)
  if (locsysman_ != Teuchos::null)
  {
    discret_->clear_state();
    discret_->set_state("dispnp", dispnp_);
    locsysman_->update(time_, {}, Global::Problem::instance()->function_manager());
    discret_->clear_state();
  }

  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  eleparams.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::instance()->function_manager());

  // Apply Dirichlet boundary conditions on provided state vector
  ALE::Ale::apply_dirichlet_bc(eleparams, dispnp_, Teuchos::null, Teuchos::null, false);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::time_step(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type)
{
  bool is_converged = false;
  int iter = 0;

  // Newton loop to deal with possible non-linearities
  while (!is_converged && iter < maxiter_)
  {
    evaluate(Teuchos::null, dbc_type);

    if (converged(iter))
    {
      is_converged = true;
      continue;
    }
    solve();
    update_iter();
    ++iter;
  }

  if (!is_converged)
  {
    switch (divercont_)
    {
      case Inpar::ALE::divcont_stop:
        FOUR_C_THROW("ALE newton not converged in %i iterations. Abort! ", maxiter_);
        break;
      case Inpar::ALE::divcont_continue:
        if (discret_->get_comm().MyPID() == 0)
        {
          Core::IO::cout << "ALE newton not converged in " << maxiter_ << " iterations. Continue"
                         << Core::IO::endl;
        }
        break;
      default:
        FOUR_C_THROW("Unknown divercont action! ");
        break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::print_time_step_header() const
{
  Core::IO::cout << "TIME: " << time_ << "/" << maxtime_ << "  DT = " << dt_ << "  STEP = " << step_
                 << "/" << numstep_ << Core::IO::endl;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::setup_dbc_map_ex(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type,
    Teuchos::RCP<const ALE::UTILS::MapExtractor> interface,
    Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor> xff_interface)
{
  // set fixed nodes (conditions != 0 are not supported right now). hahn: Why?!
  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  eleparams.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::instance()->function_manager());

  // some consistency checks
  if (interface == Teuchos::null && dbc_type != ALE::UTILS::MapExtractor::dbc_set_std &&
      dbc_type != ALE::UTILS::MapExtractor::dbc_set_x_ff)
    FOUR_C_THROW(
        "For non-standard use of SetupDBCMapEx, please provide a valid ALE::UTILS::MapExtractor.");

  if (xff_interface == Teuchos::null && dbc_type == ALE::UTILS::MapExtractor::dbc_set_x_ff)
    FOUR_C_THROW(
        "For non-standard use of SetupDBCMapEx with fluid-fluid coupling, please provide a "
        "x_fluid_fluid_map_extractor.");
  // REMARK: for all applications, setup of the standard Dirichlet sets is done in the ctor of this
  // class

  switch (dbc_type)
  {
    case ALE::UTILS::MapExtractor::dbc_set_std:
      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std] =
          Teuchos::rcp(new Core::LinAlg::MapExtractor());
      apply_dirichlet_bc(eleparams, dispnp_, Teuchos::null, Teuchos::null, true);
      break;
    case ALE::UTILS::MapExtractor::dbc_set_x_ff:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(xff_interface->xfluid_fluid_cond_map());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->cond_map());
      Teuchos::RCP<Epetra_Map> condmerged = Core::LinAlg::MultiMapExtractor::merge_maps(condmaps);

      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_x_ff] =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), condmerged));
      break;
    }
    case ALE::UTILS::MapExtractor::dbc_set_x_fsi:
    case ALE::UTILS::MapExtractor::dbc_set_biofilm:
    case ALE::UTILS::MapExtractor::dbc_set_part_fsi:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(interface->fsi_cond_map());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->cond_map());
      Teuchos::RCP<Epetra_Map> condmerged = Core::LinAlg::MultiMapExtractor::merge_maps(condmaps);

      dbcmaps_[dbc_type] =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), condmerged));
      break;
    }
    case ALE::UTILS::MapExtractor::dbc_set_wear:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(interface->ale_wear_cond_map());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->cond_map());

      Teuchos::RCP<Epetra_Map> condmerged = Core::LinAlg::MultiMapExtractor::merge_maps(condmaps);
      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_wear] =
          Teuchos::rcp(new Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), condmerged));
      break;
    }
    default:
      FOUR_C_THROW("Undefined type of ALE Dirichlet sets.");
      break;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> ALE::Ale::create_field_test()
{
  return Teuchos::rcp(new ALE::AleResultTest(*this));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::apply_dirichlet_bc(Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> systemvector, Teuchos::RCP<Epetra_Vector> systemvectord,
    Teuchos::RCP<Epetra_Vector> systemvectordd, bool recreatemap)
{
  // In the case of local coordinate systems, we have to rotate forward ...
  // ---------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (systemvector != Teuchos::null) locsysman_->rotate_global_to_local(systemvector);
    if (systemvectord != Teuchos::null) locsysman_->rotate_global_to_local(systemvectord);
    if (systemvectordd != Teuchos::null) locsysman_->rotate_global_to_local(systemvectordd);
  }

  // Apply DBCs
  // ---------------------------------------------------------------------------
  discret_->clear_state();
  if (recreatemap)
  {
    discret_->evaluate_dirichlet(params, systemvector, systemvectord, systemvectordd, Teuchos::null,
        dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]);
  }
  else
  {
    discret_->evaluate_dirichlet(
        params, systemvector, systemvectord, systemvectordd, Teuchos::null, Teuchos::null);
  }
  discret_->clear_state();

  // In the case of local coordinate systems, we have to rotate back into global Cartesian frame
  // ---------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (systemvector != Teuchos::null) locsysman_->rotate_local_to_global(systemvector);
    if (systemvectord != Teuchos::null) locsysman_->rotate_local_to_global(systemvectord);
    if (systemvectordd != Teuchos::null) locsysman_->rotate_local_to_global(systemvectordd);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::reset()
{
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  dispnp_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  dispn_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::reset_step()
{
  dispnp_->Update(1.0, *dispn_, 0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::reset_time(const double dtold)
{
  time_ = time_ - dtold;
  step_ = step_ - 1;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::set_dt(const double dtnew)
{
  dt_ = dtnew;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::SparseMatrix> ALE::Ale::get_loc_sys_trafo() const
{
  if (locsysman_ != Teuchos::null) return locsysman_->trafo();

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::update_slave_dof(Teuchos::RCP<Epetra_Vector>& a)
{
  if (msht_ != Inpar::ALE::no_meshtying)
  {
    meshtying_->update_slave_dof(a, dispnp_);
    meshtying_->recover(a);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool ALE::Ale::evaluate_element_quality()
{
  if (elequalityyesno_)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    discret_->clear_state();
    discret_->set_state("dispnp", dispnp_);

    for (int i = 0; i < discretization()->num_my_row_elements(); ++i)
    {
      Core::Elements::Element* actele;
      actele = discretization()->l_row_element(i);

      // list to define routines at elementlevel
      Teuchos::ParameterList eleparams;
      eleparams.set("action", "calc_jacobian_determinant");

      // initialize element vectors
      Core::Elements::Element::LocationArray la(discretization()->num_dof_sets());
      actele->location_vector(*discretization(), la, false);

      // only two entries per element necessary (detJ and quality measure)
      Core::LinAlg::SerialDenseMatrix elematrix1;
      Core::LinAlg::SerialDenseMatrix elematrix2;
      Core::LinAlg::SerialDenseVector elevector1(2);
      Core::LinAlg::SerialDenseVector elevector2;
      Core::LinAlg::SerialDenseVector elevector3;

      actele->evaluate(
          eleparams, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);

      eledetjac_->ReplaceMyValue(i, 0, elevector1[0]);
      elequality_->ReplaceMyValue(i, 0, elevector1[1]);

    }  // loop elements

    discret_->clear_state();

    // check for non-valid elements
    bool validshapes = true;
    double negdetjac = 0.0;
    eledetjac_->MinValue(&negdetjac);
    if (negdetjac <= 0)
    {
      validshapes = false;
      FOUR_C_THROW("Negative determinant %e in time step %i", negdetjac, step_);
    }

    return validshapes;
  }
  else
  {
    // no assesment of mesh quality. Return true to assume that everything is fine.
    return true;
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// class AleLinear ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*----------------------------------------------------------------------------*/
ALE::AleLinear::AleLinear(Teuchos::RCP<Core::FE::Discretization> actdis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params_in,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : Ale(actdis, solver, params_in, output), validsysmat_(false), updateeverystep_(false)
{
  updateeverystep_ = Core::UTILS::IntegralValue<bool>(params(), "UPDATEMATRIX");
}

/*----------------------------------------------------------------------------*/
void ALE::AleLinear::prepare_time_step()
{
  Ale::prepare_time_step();

  if (updateeverystep_) validsysmat_ = false;

  return;
}

/*----------------------------------------------------------------------------*/
void ALE::AleLinear::time_step(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type)
{
  evaluate(Teuchos::null, dbc_type);
  solve();
  update_iter();

  return;
}

/*----------------------------------------------------------------------------*/
void ALE::AleLinear::evaluate_elements()
{
  if (not validsysmat_)
  {
    Ale::evaluate_elements();

    validsysmat_ = true;
  }
  else if (not system_matrix().is_null())
    system_matrix()->Apply(*dispnp(), *write_access_residual());
  else if (not block_system_matrix().is_null())
    block_system_matrix()->Apply(*dispnp(), *write_access_residual());
  else
    FOUR_C_THROW("Can't compute residual for linear ALE.");

  return;
}

FOUR_C_NAMESPACE_CLOSE
