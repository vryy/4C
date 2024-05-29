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
#include "4C_discretization_condition_locsys.hpp"
#include "4C_discretization_condition_periodic.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
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
ALE::Ale::Ale(Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params, Teuchos::RCP<IO::DiscretizationWriter> output)
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
      elequalityyesno_(CORE::UTILS::IntegralValue<bool>(*params, "ASSESSMESHQUALITY")),
      aletype_(CORE::UTILS::IntegralValue<INPAR::ALE::AleDynamic>(*params, "ALE_TYPE")),
      maxiter_(params->get<int>("MAXITER")),
      tolres_(params->get<double>("TOLRES")),
      toldisp_(params->get<double>("TOLDISP")),
      divercont_(CORE::UTILS::IntegralValue<INPAR::ALE::DivContAct>(*params, "DIVERCONT")),
      msht_(CORE::UTILS::IntegralValue<INPAR::ALE::MeshTying>(*params, "MESHTYING")),
      initialdisp_(CORE::UTILS::IntegralValue<INPAR::ALE::InitialDisp>(*params, "INITIALDISP")),
      startfuncno_(params->get<int>("STARTFUNCNO"))
{
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  dispn_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  dispnp_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  disi_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  residual_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  rhs_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  zeros_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  eledetjac_ = CORE::LINALG::CreateVector(*discretization()->ElementRowMap(), true);
  elequality_ = CORE::LINALG::CreateVector(*discretization()->ElementRowMap(), true);

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  set_initial_displacement(initialdisp_, startfuncno_);

  SetupDBCMapEx();

  // ensure that the ALE string was removed from conditions
  {
    CORE::Conditions::Condition* cond = discret_->GetCondition("ALEDirichlet");
    if (cond) FOUR_C_THROW("Found a ALE Dirichlet condition. Remove ALE string!");
  }

  if (msht_ == INPAR::ALE::meshsliding)
  {
    meshtying_ = Teuchos::rcp(
        new Meshsliding(discret_, *solver_, msht_, GLOBAL::Problem::Instance()->NDim(), nullptr));
  }
  else if (msht_ == INPAR::ALE::meshtying)
  {
    meshtying_ = Teuchos::rcp(
        new Meshtying(discret_, *solver_, msht_, GLOBAL::Problem::Instance()->NDim(), nullptr));
  }

  // ---------------------------------------------------------------------
  // Create LocSysManager, if needed (used for LocSys-Dirichlet BCs)
  // ---------------------------------------------------------------------
  {
    std::vector<CORE::Conditions::Condition*> locsysconditions(0);
    discret_->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      // Initialize locsys manager
      locsysman_ = Teuchos::rcp(
          new CORE::Conditions::LocsysManager(*discret_, GLOBAL::Problem::Instance()->NDim()));
    }
  }

  create_system_matrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ALE::Ale::set_initial_displacement(const INPAR::ALE::InitialDisp init, const int startfuncno)
{
  switch (init)
  {
    case INPAR::ALE::initdisp_zero_disp:
    {
      dispn_->PutScalar(0.0);
      dispnp_->PutScalar(0.0);

      break;
    }
    case INPAR::ALE::initdisp_disp_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        CORE::Nodes::Node* lnode = discret_->lRowNode(lnodeid);

        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // loop nodal dofs
        for (unsigned int d = 0; d < nodedofset.size(); ++d)
        {
          const int dofgid = nodedofset[d];
          int doflid = dofrowmap->LID(dofgid);

          // evaluate component d of function
          double initialval = GLOBAL::Problem::Instance()
                                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                                  .Evaluate(lnode->X().data(), 0, d);

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
  if (msht_ != INPAR::ALE::no_meshtying)
  {
    std::vector<int> coupleddof(GLOBAL::Problem::Instance()->NDim(), 1);
    sysmat_ = meshtying_->Setup(coupleddof, dispnp_);
    meshtying_->DirichletOnMaster(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());

    if (interface != Teuchos::null)
    {
      meshtying_->IsMultifield(*interface, true);
    }
  }
  else if (interface == Teuchos::null)
  {
    const Epetra_Map* dofrowmap = discret_->dof_row_map();
    sysmat_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*dofrowmap, 81, false, true));
  }
  else
  {
    sysmat_ =
        Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
            *interface, *interface, 81, false, true));
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::Evaluate(
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

  if (msht_ != INPAR::ALE::no_meshtying)
  {
    meshtying_->MshtSplit(sysmat_);
  }

  evaluate_elements();
  evaluate_element_quality();

  // prepare meshtying system
  if (msht_ != INPAR::ALE::no_meshtying)
  {
    meshtying_->prepare_meshtying_system(sysmat_, residual_, dispnp_);
    meshtying_->MultifieldSplit(sysmat_);
  }

  // dispnp_ has zeros at the Dirichlet-entries, so we maintain zeros there.
  if (LocsysManager() != Teuchos::null)
  {
    // Transform system matrix and rhs to local coordinate systems
    LocsysManager()->RotateGlobalToLocal(
        Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(sysmat_), residual_);

    // When using local systems, a rotated dispnp_ vector needs to be used as dbcval for
    // apply_dirichlet_to_system
    Teuchos::RCP<Epetra_Vector> dispnp_local = Teuchos::rcp(new Epetra_Vector(*(zeros_)));
    LocsysManager()->RotateGlobalToLocal(dispnp_local);

    if (get_loc_sys_trafo() != Teuchos::null)
    {
      CORE::LINALG::apply_dirichlet_to_system(
          *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(sysmat_), *disi_, *residual_,
          *get_loc_sys_trafo(), *dispnp_local, *(dbcmaps_[dbc_type]->CondMap()));
    }
    else
    {
      CORE::LINALG::apply_dirichlet_to_system(
          *sysmat_, *disi_, *residual_, *dispnp_local, *(dbcmaps_[dbc_type]->CondMap()));
    }
  }
  else
  {
    CORE::LINALG::apply_dirichlet_to_system(
        *sysmat_, *disi_, *residual_, *zeros_, *(dbcmaps_[dbc_type]->CondMap()));
  }

  /* residual_ contains the most recent "mechanical" residual including DBCs.
   * We make this negative and store it in rhs_ for use in Newton-type methods.
   */
  rhs_->Update(-1.0, *residual_, 0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ALE::Ale::Solve()
{
  // We need the negative residual here as right hand side of the linear problem
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(*residual_));
  rhs->Scale(-1.0);

  // ToDo (mayr) Why can't we use rhs_ instead of local variable rhs???
  int errorcode = 0;
  if (msht_ == INPAR::ALE::no_meshtying)
  {
    CORE::LINALG::SolverParams solver_params;
    solver_params.refactor = true;
    errorcode = solver_->Solve(sysmat_->EpetraOperator(), disi_, rhs, solver_params);
  }
  else
    errorcode = meshtying_->SolveMeshtying(*solver_, sysmat_, disi_, rhs, dispnp_);
  // calc norm
  disi_->Norm2(&normdisi_);
  normdisi_ /= sqrt(disi_->GlobalLength());

  return errorcode;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::UpdateIter() { dispnp_->Update(1.0, *disi_, 1.0); }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::Update() { dispn_->Update(1.0, *dispnp_, 0.0); }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool ALE::Ale::Converged(const int iter)
{
  if (iter == 0) normdisi_ = 0.0;

  bool converged = false;
  // determine norms
  double res_norm;
  residual_->Norm2(&res_norm);
  res_norm /= sqrt(residual_->GlobalLength());
  if (discret_->Comm().MyPID() == 0)
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
  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // set vector values needed by elements
  discret_->ClearState();

  // action for elements
  eleparams.set<std::string>("action", ElementActionString(aletype_));
  eleparams.set<bool>("use spatial configuration", update_sys_mat_every_step());

  discret_->set_state("dispnp", dispnp_);

  discret_->Evaluate(eleparams, sysmat_, residual_);
  discret_->ClearState();

  sysmat_->Complete();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::string ALE::Ale::ElementActionString(const enum INPAR::ALE::AleDynamic name)
{
  switch (name)
  {
    case INPAR::ALE::solid:
      return "calc_ale_solid";
      break;
    case INPAR::ALE::solid_linear:
      return "calc_ale_solid_linear";
      break;
    case INPAR::ALE::laplace_material:
      return "calc_ale_laplace_material";
      break;
    case INPAR::ALE::laplace_spatial:
      return "calc_ale_laplace_spatial";
      break;
    case INPAR::ALE::springs_material:
      return "calc_ale_springs_material";
      break;
    case INPAR::ALE::springs_spatial:
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
Teuchos::RCP<CORE::LINALG::SparseMatrix> ALE::Ale::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(sysmat_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ALE::Ale::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(sysmat_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ALE::Ale::Integrate()
{
  // add eps to prevent stopping one step too early due to memory trash on last digits
  const double eps = 1.0e-12;
  while (step_ < numstep_ and time_ <= maxtime_ + eps)
  {
    prepare_time_step();
    TimeStep();
    Update();
    Output();
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::Output()
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
  if ((step_ == writeresultsevery_ or step_ == 0)) output_->WriteElementData(true);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::output_state(bool& datawritten)
{
  // write output data
  output_->NewStep(step_, time_);
  output_->WriteVector("dispnp", dispnp_);

  if (elequalityyesno_)
  {
    output_->WriteVector("det_j", eledetjac_, IO::elementvector);
    output_->WriteVector("element_quality", elequality_, IO::elementvector);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::output_restart(bool& datawritten)
{
  // write restart data
  output_->NewStep(step_, time_);
  output_->WriteVector("dispnp", dispnp_);
  output_->WriteVector("dispn", dispn_);

  // restart/output data has been written
  datawritten = true;

  // info dedicated to user's eyes staring at standard out
  // Print restart info only in case of pure ALE problem. Coupled problems
  // print their own restart info.
  if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::ale)
  {
    if (discret_->Comm().MyPID() == 0)
      IO::cout << "====== Restart written in step " << step_ << IO::endl;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::read_restart(const int step)
{
  IO::DiscretizationReader reader(discret_, GLOBAL::Problem::Instance()->InputControlFile(), step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(dispnp_, "dispnp");
  reader.ReadVector(dispn_, "dispn");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;

  // Print time step header only in case of pure ALE problem. Coupled problems
  // print their own time step header.
  if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::ale)
    print_time_step_header();

  // Update local coordinate systems (which may be time dependent)
  if (locsysman_ != Teuchos::null)
  {
    discret_->ClearState();
    discret_->set_state("dispnp", dispnp_);
    locsysman_->Update(time_, {}, GLOBAL::Problem::Instance()->FunctionManager());
    discret_->ClearState();
  }

  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);

  // Apply Dirichlet boundary conditions on provided state vector
  ALE::Ale::apply_dirichlet_bc(eleparams, dispnp_, Teuchos::null, Teuchos::null, false);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::TimeStep(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type)
{
  bool converged = false;
  int iter = 0;

  // Newton loop to deal with possible non-linearities
  while (!converged && iter < maxiter_)
  {
    Evaluate(Teuchos::null, dbc_type);

    if (Converged(iter))
    {
      converged = true;
      continue;
    }
    Solve();
    UpdateIter();
    ++iter;
  }

  if (!converged)
  {
    switch (divercont_)
    {
      case INPAR::ALE::divcont_stop:
        FOUR_C_THROW("ALE newton not converged in %i iterations. Abort! ", maxiter_);
        break;
      case INPAR::ALE::divcont_continue:
        if (discret_->Comm().MyPID() == 0)
        {
          IO::cout << "ALE newton not converged in " << maxiter_ << " iterations. Continue"
                   << IO::endl;
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
  IO::cout << "TIME: " << time_ << "/" << maxtime_ << "  DT = " << dt_ << "  STEP = " << step_
           << "/" << numstep_ << IO::endl;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::SetupDBCMapEx(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type,
    Teuchos::RCP<const ALE::UTILS::MapExtractor> interface,
    Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor> xff_interface)
{
  // set fixed nodes (conditions != 0 are not supported right now). hahn: Why?!
  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);

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
          Teuchos::rcp(new CORE::LINALG::MapExtractor());
      apply_dirichlet_bc(eleparams, dispnp_, Teuchos::null, Teuchos::null, true);
      break;
    case ALE::UTILS::MapExtractor::dbc_set_x_ff:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(xff_interface->XFluidFluidCondMap());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());
      Teuchos::RCP<Epetra_Map> condmerged = CORE::LINALG::MultiMapExtractor::MergeMaps(condmaps);

      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_x_ff] =
          Teuchos::rcp(new CORE::LINALG::MapExtractor(*(discret_->dof_row_map()), condmerged));
      break;
    }
    case ALE::UTILS::MapExtractor::dbc_set_x_fsi:
    case ALE::UTILS::MapExtractor::dbc_set_biofilm:
    case ALE::UTILS::MapExtractor::dbc_set_part_fsi:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(interface->FSICondMap());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());
      Teuchos::RCP<Epetra_Map> condmerged = CORE::LINALG::MultiMapExtractor::MergeMaps(condmaps);

      dbcmaps_[dbc_type] =
          Teuchos::rcp(new CORE::LINALG::MapExtractor(*(discret_->dof_row_map()), condmerged));
      break;
    }
    case ALE::UTILS::MapExtractor::dbc_set_wear:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(interface->AleWearCondMap());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());

      Teuchos::RCP<Epetra_Map> condmerged = CORE::LINALG::MultiMapExtractor::MergeMaps(condmaps);
      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_wear] =
          Teuchos::rcp(new CORE::LINALG::MapExtractor(*(discret_->dof_row_map()), condmerged));
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
Teuchos::RCP<CORE::UTILS::ResultTest> ALE::Ale::CreateFieldTest()
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
    if (systemvector != Teuchos::null) locsysman_->RotateGlobalToLocal(systemvector);
    if (systemvectord != Teuchos::null) locsysman_->RotateGlobalToLocal(systemvectord);
    if (systemvectordd != Teuchos::null) locsysman_->RotateGlobalToLocal(systemvectordd);
  }

  // Apply DBCs
  // ---------------------------------------------------------------------------
  discret_->ClearState();
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
  discret_->ClearState();

  // In the case of local coordinate systems, we have to rotate back into global Cartesian frame
  // ---------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (systemvector != Teuchos::null) locsysman_->RotateLocalToGlobal(systemvector);
    if (systemvectord != Teuchos::null) locsysman_->RotateLocalToGlobal(systemvectord);
    if (systemvectordd != Teuchos::null) locsysman_->RotateLocalToGlobal(systemvectordd);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::Reset()
{
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  dispnp_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  dispn_ = CORE::LINALG::CreateVector(*dofrowmap, true);

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
Teuchos::RCP<const CORE::LINALG::SparseMatrix> ALE::Ale::get_loc_sys_trafo() const
{
  if (locsysman_ != Teuchos::null) return locsysman_->Trafo();

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& a)
{
  if (msht_ != INPAR::ALE::no_meshtying)
  {
    meshtying_->UpdateSlaveDOF(a, dispnp_);
    meshtying_->Recover(a);
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

    discret_->ClearState();
    discret_->set_state("dispnp", dispnp_);

    for (int i = 0; i < discretization()->NumMyRowElements(); ++i)
    {
      CORE::Elements::Element* actele;
      actele = discretization()->lRowElement(i);

      // list to define routines at elementlevel
      Teuchos::ParameterList eleparams;
      eleparams.set("action", "calc_jacobian_determinant");

      // initialize element vectors
      CORE::Elements::Element::LocationArray la(discretization()->NumDofSets());
      actele->LocationVector(*discretization(), la, false);

      // only two entries per element necessary (detJ and quality measure)
      CORE::LINALG::SerialDenseMatrix elematrix1;
      CORE::LINALG::SerialDenseMatrix elematrix2;
      CORE::LINALG::SerialDenseVector elevector1(2);
      CORE::LINALG::SerialDenseVector elevector2;
      CORE::LINALG::SerialDenseVector elevector3;

      actele->Evaluate(
          eleparams, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);

      eledetjac_->ReplaceMyValue(i, 0, elevector1[0]);
      elequality_->ReplaceMyValue(i, 0, elevector1[1]);

    }  // loop elements

    discret_->ClearState();

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
ALE::AleLinear::AleLinear(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params_in,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : Ale(actdis, solver, params_in, output), validsysmat_(false), updateeverystep_(false)
{
  updateeverystep_ = CORE::UTILS::IntegralValue<bool>(params(), "UPDATEMATRIX");
}

/*----------------------------------------------------------------------------*/
void ALE::AleLinear::prepare_time_step()
{
  Ale::prepare_time_step();

  if (updateeverystep_) validsysmat_ = false;

  return;
}

/*----------------------------------------------------------------------------*/
void ALE::AleLinear::TimeStep(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type)
{
  Evaluate(Teuchos::null, dbc_type);
  Solve();
  UpdateIter();

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
  else if (not SystemMatrix().is_null())
    SystemMatrix()->Apply(*Dispnp(), *write_access_residual());
  else if (not BlockSystemMatrix().is_null())
    BlockSystemMatrix()->Apply(*Dispnp(), *write_access_residual());
  else
    FOUR_C_THROW("Can't compute residual for linear ALE.");

  return;
}

FOUR_C_NAMESPACE_CLOSE
