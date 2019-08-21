/*----------------------------------------------------------------------------*/
/*! \file

\brief ALE time integration

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale.H"
#include "ale_meshsliding.H"
#include "ale_resulttest.H"
#include "ale_utils_mapextractor.H"

// further includes for AleBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_locsys.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_precond.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_inpar/inpar_ale.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_lib/drt_periodicbc.H"

#include "../drt_io/io_pstream.H"

#define SCALING_INFNORM true

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ALE::Ale::Ale(Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
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
      elequalityyesno_(DRT::INPUT::IntegralValue<bool>(*params, "ASSESSMESHQUALITY")),
      precond_(Teuchos::null),
      aletype_(DRT::INPUT::IntegralValue<INPAR::ALE::AleDynamic>(*params, "ALE_TYPE")),
      maxiter_(params->get<int>("MAXITER")),
      tolres_(params->get<double>("TOLRES")),
      toldisp_(params->get<double>("TOLDISP")),
      divercont_(DRT::INPUT::IntegralValue<INPAR::ALE::DivContAct>(*params, "DIVERCONT")),
      msht_(DRT::INPUT::IntegralValue<INPAR::ALE::MeshTying>(*params, "MESHTYING"))
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispn_ = LINALG::CreateVector(*dofrowmap, true);
  dispnp_ = LINALG::CreateVector(*dofrowmap, true);
  disi_ = LINALG::CreateVector(*dofrowmap, true);
  residual_ = LINALG::CreateVector(*dofrowmap, true);
  rhs_ = LINALG::CreateVector(*dofrowmap, true);
  zeros_ = LINALG::CreateVector(*dofrowmap, true);

  eledetjac_ = LINALG::CreateVector(*Discretization()->ElementRowMap(), true);
  elequality_ = LINALG::CreateVector(*Discretization()->ElementRowMap(), true);

  SetupDBCMapEx();

  // ensure that the ALE string was removed from conditions
  {
    DRT::Condition* cond = discret_->GetCondition("ALEDirichlet");
    if (cond) dserror("Found a ALE Dirichlet condition. Remove ALE string!");
  }

  if (msht_ == INPAR::ALE::meshsliding)
  {
    meshtying_ = Teuchos::rcp(
        new Meshsliding(discret_, *solver_, msht_, DRT::Problem::Instance()->NDim(), NULL));
  }
  else if (msht_ == INPAR::ALE::meshtying)
  {
    meshtying_ = Teuchos::rcp(
        new Meshtying(discret_, *solver_, msht_, DRT::Problem::Instance()->NDim(), NULL));
  }

  // ---------------------------------------------------------------------
  // Create LocSysManager, if needed (used for LocSys-Dirichlet BCs)
  // ---------------------------------------------------------------------
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    discret_->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      // Initialize locsys manager
      locsysman_ = Teuchos::rcp(new DRT::UTILS::LocsysManager(*discret_));
    }
  }

  CreateSystemMatrix();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::CreateSystemMatrix(Teuchos::RCP<const ALE::UTILS::MapExtractor> interface)
{
  if (msht_ != INPAR::ALE::no_meshtying)
  {
    std::vector<int> coupleddof(DRT::Problem::Instance()->NDim(), 1);
    sysmat_ = meshtying_->Setup(coupleddof, dispnp_);
    meshtying_->DirichletOnMaster(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());

    if (interface != Teuchos::null)
    {
      meshtying_->IsMultifield(*interface, true);
    }
  }
  else if (interface == Teuchos::null)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 81, false, true));
  }
  else
  {
    sysmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
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

  EvaluateElements();
  EvaluateElementQuality();

  // prepare meshtying system
  if (msht_ != INPAR::ALE::no_meshtying)
  {
    meshtying_->PrepareMeshtyingSystem(sysmat_, residual_, dispnp_);
    meshtying_->MultifieldSplit(sysmat_);
  }

  // dispnp_ has zeros at the Dirichlet-entries, so we maintain zeros there.
  if (LocsysManager() != Teuchos::null)
  {
    // Transform system matrix and rhs to local coordinate systems
    LocsysManager()->RotateGlobalToLocal(
        Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_), residual_);

    // When using local systems, a rotated dispnp_ vector needs to be used as dbcval for
    // ApplyDirichlettoSystem
    Teuchos::RCP<Epetra_Vector> dispnp_local = Teuchos::rcp(new Epetra_Vector(*(zeros_)));
    LocsysManager()->RotateGlobalToLocal(dispnp_local);

    LINALG::ApplyDirichlettoSystem(sysmat_, disi_, residual_, GetLocSysTrafo(), dispnp_local,
        *(dbcmaps_[dbc_type]->CondMap()));
  }
  else
  {
    LINALG::ApplyDirichlettoSystem(
        sysmat_, disi_, residual_, zeros_, *(dbcmaps_[dbc_type]->CondMap()));
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
    errorcode = solver_->Solve(sysmat_->EpetraOperator(), disi_, rhs, true);
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
void ALE::Ale::EvaluateElements()
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
  eleparams.set<bool>("use spatial configuration", UpdateSysMatEveryStep());

  discret_->SetState("dispnp", dispnp_);

  discret_->Evaluate(eleparams, sysmat_, residual_);
  discret_->ClearState();

  sysmat_->Complete();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const std::string ALE::Ale::ElementActionString(const enum INPAR::ALE::AleDynamic name)
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
      dserror("Cannot make std::string for ALE type %d", name);
      return "";
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Preconditioner> ALE::Ale::ConstPreconditioner()
{
  // TODO (mayr) fix const preconditioner stuff

  return precond_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::CreatePreconditioner(bool full)
{
  if (full)
  {
    // partitioned FSI does not use explicit preconditioner objects
  }
  else
  {
    // This is the MFSI case and we need the preconditioner on the inner dofs only
    precond_ = Teuchos::rcp(new LINALG::Preconditioner(solver_));

    Teuchos::RCP<Epetra_CrsMatrix> A = BlockSystemMatrix()->Matrix(0, 0).EpetraMatrix();

    Teuchos::RCP<Epetra_Vector> arowsum;
    Teuchos::RCP<Epetra_Vector> acolsum;

    if (SCALING_INFNORM)
    {
      arowsum = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
      acolsum = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
      A->InvRowSums(*arowsum);
      A->InvColSums(*acolsum);
      if (A->LeftScale(*arowsum) or A->RightScale(*acolsum)) dserror("ale scaling failed");
    }

    precond_->Setup(A);

    if (SCALING_INFNORM)
    {
      arowsum->Reciprocal(*arowsum);
      acolsum->Reciprocal(*acolsum);
      if (A->LeftScale(*arowsum) or A->RightScale(*acolsum)) dserror("ale scaling failed");
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ALE::Ale::DofRowMap() const
{
  return Teuchos::rcp(discret_->DofRowMap(), false);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ALE::Ale::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ALE::Ale::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int ALE::Ale::Integrate()
{
  // add eps to prevent stopping one step too early due to memory trash on last digits
  const double eps = 1.0e-12;
  while (step_ < numstep_ and time_ <= maxtime_ + eps)
  {
    PrepareTimeStep();
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
    OutputRestart(datawritten);
  }

  // write output data if necessary
  if (not datawritten and writeresultsevery_ != 0 and step_ % writeresultsevery_ == 0)
  {
    OutputState(datawritten);
  }

  // write domain decomposition for visualization
  if ((step_ == writeresultsevery_ or step_ == 0)) output_->WriteElementData(true);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::OutputState(bool& datawritten)
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
void ALE::Ale::OutputRestart(bool& datawritten)
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
  if (DRT::Problem::Instance()->ProblemType() == prb_ale)
  {
    if (discret_->Comm().MyPID() == 0)
      IO::cout << "====== Restart written in step " << step_ << IO::endl;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::ReadRestart(const int step)
{
  IO::DiscretizationReader reader(discret_, step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(dispnp_, "dispnp");
  reader.ReadVector(dispn_, "dispn");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  // Print time step header only in case of pure ALE problem. Coupled problems
  // print their own time step header.
  if (DRT::Problem::Instance()->ProblemType() == prb_ale) PrintTimeStepHeader();

  // Update local coordinate systems (which may be time dependent)
  if (locsysman_ != Teuchos::null)
  {
    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
    locsysman_->Setup(time_);
    discret_->ClearState();
  }

  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);

  // Apply Dirichlet boundary conditions on provided state vector
  ALE::Ale::ApplyDirichletBC(eleparams, dispnp_, Teuchos::null, Teuchos::null, false);

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
        dserror("ALE newton not converged in %i iterations. Abort! ", maxiter_);
        break;
      case INPAR::ALE::divcont_continue:
        if (discret_->Comm().MyPID() == 0)
        {
          IO::cout << "ALE newton not converged in " << maxiter_ << " iterations. Continue"
                   << IO::endl;
        }
        break;
      default:
        dserror("Unknown divercont action! ");
        break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::PrintTimeStepHeader() const
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
    dserror(
        "For non-standard use of SetupDBCMapEx, please provide a valid ALE::UTILS::MapExtractor.");

  if (xff_interface == Teuchos::null && dbc_type == ALE::UTILS::MapExtractor::dbc_set_x_ff)
    dserror(
        "For non-standard use of SetupDBCMapEx with fluid-fluid coupling, please provide a "
        "XFluidFluidMapExtractor.");
  // REMARK: for all applications, setup of the standard Dirichlet sets is done in the ctor of this
  // class

  switch (dbc_type)
  {
    case ALE::UTILS::MapExtractor::dbc_set_std:
      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std] = Teuchos::rcp(new LINALG::MapExtractor());
      ApplyDirichletBC(eleparams, dispnp_, Teuchos::null, Teuchos::null, true);
      break;
    case ALE::UTILS::MapExtractor::dbc_set_x_ff:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(xff_interface->XFluidFluidCondMap());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());
      Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);

      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_x_ff] =
          Teuchos::rcp(new LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged));
      break;
    }
    case ALE::UTILS::MapExtractor::dbc_set_x_fsi:
    case ALE::UTILS::MapExtractor::dbc_set_biofilm:
    case ALE::UTILS::MapExtractor::dbc_set_part_fsi:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(interface->FSICondMap());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());
      Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);

      dbcmaps_[dbc_type] =
          Teuchos::rcp(new LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged));
      break;
    }
    case ALE::UTILS::MapExtractor::dbc_set_wear:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
      condmaps.push_back(interface->AleWearCondMap());
      condmaps.push_back(dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]->CondMap());

      Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
      dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_wear] =
          Teuchos::rcp(new LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged));
      break;
    }
    default:
      dserror("Undefined type of ALE Dirichlet sets.");
      break;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ALE::Ale::CreateFieldTest()
{
  return Teuchos::rcp(new ALE::AleResultTest(*this));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::ApplyDirichletBC(Teuchos::ParameterList& params,
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
    discret_->EvaluateDirichlet(params, systemvector, systemvectord, systemvectordd, Teuchos::null,
        dbcmaps_[ALE::UTILS::MapExtractor::dbc_set_std]);
  }
  else
  {
    discret_->EvaluateDirichlet(
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
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispnp_ = LINALG::CreateVector(*dofrowmap, true);
  dispn_ = LINALG::CreateVector(*dofrowmap, true);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::ResetStep()
{
  dispnp_->Update(1.0, *dispn_, 0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::ResetTime(const double dtold)
{
  time_ = time_ - dtold;
  step_ = step_ - 1;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::Ale::SetDt(const double dtnew)
{
  dt_ = dtnew;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix> ALE::Ale::GetLocSysTrafo() const
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
bool ALE::Ale::EvaluateElementQuality()
{
  if (elequalityyesno_)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);

    for (int i = 0; i < Discretization()->NumMyRowElements(); ++i)
    {
      DRT::Element* actele;
      actele = Discretization()->lRowElement(i);

      // list to define routines at elementlevel
      Teuchos::ParameterList eleparams;
      eleparams.set("action", "calc_jacobian_determinant");

      // initialize element vectors
      DRT::Element::LocationArray la(Discretization()->NumDofSets());
      actele->LocationVector(*Discretization(), la, false);

      // only two entries per element necessary (detJ and quality measure)
      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1(2);
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;

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
      dserror("Negative determinant %e in time step %i", negdetjac, step_);
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
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : Ale(actdis, solver, params, output), validsysmat_(false), updateeverystep_(false)
{
  updateeverystep_ = DRT::INPUT::IntegralValue<bool>(Params(), "UPDATEMATRIX");

  return;
}

/*----------------------------------------------------------------------------*/
void ALE::AleLinear::PrepareTimeStep()
{
  Ale::PrepareTimeStep();

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
void ALE::AleLinear::EvaluateElements()
{
  if (not validsysmat_)
  {
    Ale::EvaluateElements();

    validsysmat_ = true;
  }
  else if (not SystemMatrix().is_null())
    SystemMatrix()->Apply(*Dispnp(), *WriteAccessResidual());
  else if (not BlockSystemMatrix().is_null())
    BlockSystemMatrix()->Apply(*Dispnp(), *WriteAccessResidual());
  else
    dserror("Can't compute residual for linear ALE.");

  return;
}
