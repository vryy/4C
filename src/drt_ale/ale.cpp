/*----------------------------------------------------------------------*/
/*!
 * \file ale.cpp
 *
\brief ALE base implementation

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/



#include "ale.H"
#include "ale_lin.H"
#include "ale_laplace.H"
#include "ale_springs.H"
#include "ale_springs_fixed_ref.H"
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
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_inpar/inpar_ale.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_fpsi.H"
#include "../drt_fluid/drt_periodicbc.H"

#include "../drt_io/io_pstream.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ALE::Ale::Ale(Teuchos::RCP<DRT::Discretization> actdis,
              Teuchos::RCP<LINALG::Solver> solver,
              Teuchos::RCP<Teuchos::ParameterList> params,
              Teuchos::RCP<IO::DiscretizationWriter> output,
              bool dirichletcond)
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
    sysmat_(Teuchos::null)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispn_          = LINALG::CreateVector(*dofrowmap,true);
  dispnp_         = LINALG::CreateVector(*dofrowmap,true);
  residual_       = LINALG::CreateVector(*dofrowmap,true);

  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*actdis);

  SetupDBCMapEx(dirichletcond);

  // ensure that the ALE string was removed from conditions
  {
    DRT::Condition* cond = discret_->GetCondition("ALEDirichlet");
    if (cond) dserror("Found a ALE Dirichlet condition. Remove ALE string!");
  }

  // ---------------------------------------------------------------------
  // Create LocSysManager, if needed (used for LocSys-Dirichlet BCs)
  // ---------------------------------------------------------------------
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    discret_->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      // Only 'classic_lin' and 'springs' are supported for use with locsys.
      // If you want to use locsys with another ALE type, just copy from 'classic_lin'
      int aletype = DRT::INPUT::IntegralValue<int>(*params_,"ALE_TYPE");
      if (not ((aletype==INPAR::ALE::classic_lin) || (aletype==INPAR::ALE::springs)))
        dserror("Only ALE types 'classic_lin' and 'springs' are supported for use with locsys conditions.");

      // Initialize locsys manager
      locsysman_ = Teuchos::rcp(new DRT::UTILS::LocsysManager(*discret_));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ALE::Ale::DofRowMap()
{
  return Teuchos::rcp(discret_->DofRowMap(),false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ALE::Ale::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ALE::Ale::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::Integrate()
{
  const double eps = 1.0e-12; // add eps to prevent stopping one step too early due to memory trash on last digits
  while (step_ < numstep_ and time_ <= maxtime_ + eps)
  {
    PrepareTimeStep();
    Solve();
    Update();
    Output();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
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

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::OutputState(bool& datawritten)
{
  // write output data
  output_->NewStep(step_,time_);
  output_->WriteVector("dispnp", dispnp_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::OutputRestart(bool& datawritten)
{
  // write restart data
  output_->NewStep(step_,time_);
  output_->WriteVector("dispnp", dispnp_);
  output_->WriteVector("dispn", dispn_);

  // restart/output data has been written
  datawritten = true;

  // info dedicated to user's eyes staring at standard out
  if (discret_->Comm().MyPID() == 0)
  {
    IO::cout << "====== Restart written in step " << step_ << IO::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::ReadRestart(const int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(dispnp_, "dispnp");
  reader.ReadVector(dispn_,  "dispn");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  // Print time step header only in case of pure ALE problem. Coupled problems
  // print their own time step header.
  if (DRT::Problem::Instance()->ProblemType() == prb_ale)
    PrintTimeStepHeader();

  // Update local coordinate systems (which may be time dependent)
  if (locsysman_ != Teuchos::null)
  {
    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
    locsysman_->Setup(time_);
    discret_->ClearState();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::PrintTimeStepHeader() const
{
  IO::cout << "TIME: " << time_ << "/" << maxtime_
           << "  DT = " << dt_
           << "  STEP = " << step_ << "/" << numstep_
           << IO::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::SetupDBCMapEx(bool dirichletcond)
{
  // set fixed nodes (conditions != 0 are not supported right now). hahn: Why?!
  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  ApplyDirichletBC(eleparams, dispnp_, Teuchos::null, Teuchos::null, true);

  if (dirichletcond)
  {
    // for partitioned FSI the interface becomes a Dirichlet boundary
    // also for structural Lagrangian simulations with contact and wear
    // followed by an Eulerian step to take wear into account, the interface
    // becomes a dirichlet
    std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
    condmaps.push_back(interface_->FSICondMap());
    condmaps.push_back(interface_->AleWearCondMap());
    condmaps.push_back(dbcmaps_->CondMap());
    Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
    *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
  }

  if (dirichletcond and interface_->FSCondRelevant())
  {
    // for partitioned solves the free surface becomes a Dirichlet boundary
    std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
    condmaps.push_back(interface_->FSCondMap());
    condmaps.push_back(dbcmaps_->CondMap());
    Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
    *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ALE::Ale::CreateFieldTest()
{
  return Teuchos::rcp(new ALE::AleResultTest(*this));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp)
{
  interface_->InsertFSCondVector(fsdisp,dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  // applying interface displacements
  if(DRT::Problem::Instance()->ProblemType()!=prb_struct_ale)
    interface_->InsertFSICondVector(idisp,dispnp_);
  else
    interface_->AddAleWearCondVector(idisp,dispnp_);
}

/*---------------------------------------------------------------*/
/* Apply Dirichlet boundary conditions on provided state vectors */
void ALE::Ale::ApplyDirichletBC
(
  Teuchos::ParameterList& params,
  Teuchos::RCP<Epetra_Vector> systemvector,   //!< (may be Teuchos::null)
  Teuchos::RCP<Epetra_Vector> systemvectord,  //!< (may be Teuchos::null)
  Teuchos::RCP<Epetra_Vector> systemvectordd, //!< (may be Teuchos::null)
  bool recreatemap  //!< recreate mapextractor/toggle-vector
)
{
  // In the case of local coordinate systems, we have to rotate forward ...
  // --------------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (systemvector != Teuchos::null)
      locsysman_->RotateGlobalToLocal(systemvector);
    if (systemvectord != Teuchos::null)
      locsysman_->RotateGlobalToLocal(systemvectord);
    if (systemvectordd != Teuchos::null)
      locsysman_->RotateGlobalToLocal(systemvectordd);
  }

  // Apply DBCs
  // --------------------------------------------------------------------------------
  discret_->ClearState();
  if (recreatemap)
  {
    discret_->EvaluateDirichlet(params, systemvector, systemvectord, systemvectordd,
                                Teuchos::null, dbcmaps_);
  }
  else
  {
    discret_->EvaluateDirichlet(params, systemvector, systemvectord, systemvectordd,
                               Teuchos::null, Teuchos::null);
  }
  discret_->ClearState();

  // In the case of local coordinate systems, we have to rotate back into global Cartesian frame
  // --------------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    if (systemvector != Teuchos::null)
      locsysman_->RotateLocalToGlobal(systemvector);
    if (systemvectord != Teuchos::null)
      locsysman_->RotateLocalToGlobal(systemvectord);
    if (systemvectordd != Teuchos::null)
      locsysman_->RotateLocalToGlobal(systemvectordd);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::Reset()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispnp_ = LINALG::CreateVector(*dofrowmap,true);
  dispn_  = LINALG::CreateVector(*dofrowmap,true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::ResetStep()
{
  dispnp_ ->Update(1.0, *dispn_,0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::ResetTime(const double dtold)
{
  time_=time_-dtold;
  step_=step_-1;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::Ale::SetDt(const double dtnew)
{
  dt_ = dtnew;

  return;
}

/*----------------------------------------------------------------------*/
/* Return (rotatory) transformation matrix of local coordinate systems  */
Teuchos::RCP<const LINALG::SparseMatrix> ALE::Ale::GetLocSysTrafo() const
{
  if (locsysman_ != Teuchos::null)
    return locsysman_->Trafo();

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ALE::AleBaseAlgorithm::AleBaseAlgorithm(const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  SetupAle(prbdyn,actdis);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ALE::AleBaseAlgorithm::~AleBaseAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::AleBaseAlgorithm::SetupAle(const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ALE::AleBaseAlgorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // what's the current problem type?
  const PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // connect degrees of freedom for coupled nodes
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions pbc(actdis);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> adyn
    = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->AleDynamicParams()));

  // -------------------------------------------------------------------
  // create a linear solver
  // -------------------------------------------------------------------
  // get the linear solver number
  const int linsolvernumber = adyn->get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for ALE problems. Please set LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // ---------------------------------------------------------------------------
  // overwrite certain parameters when ALE is part of a multi-field problem
  // ---------------------------------------------------------------------------
  adyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  adyn->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  adyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  adyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));

  if (probtype == prb_ale
      or probtype == prb_struct_ale
      or probtype == prb_structure
      or probtype == prb_redairways_tissue
      or probtype == prb_particle)
  {
    adyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));
  }
  else
  {
    adyn->set<int>("RESULTSEVRY", prbdyn.get<int>("UPRES"));
  }


  bool dirichletcond = true;
  if (probtype == prb_fsi or
      probtype == prb_fsi_redmodels or
      probtype == prb_fsi_lung or
      probtype == prb_gas_fsi or
      probtype == prb_thermo_fsi or
      probtype == prb_biofilm_fsi or
      probtype == prb_fluid_fluid_fsi)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicfluidsplit or
        coupling == fsi_iter_constr_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicfluidsplit or
        coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_mortar_monolithicstructuresplit or
        coupling == fsi_iter_mortar_monolithicfluidsplit or
        coupling == fsi_iter_fluidfluid_monolithicstructuresplit or
        coupling == fsi_iter_fluidfluid_monolithicfluidsplit or
        coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nox)
    {
        dirichletcond = false;
    }
  }

  if (probtype == prb_fpsi)
  {
    // FPSI input parameters
    const Teuchos::ParameterList&  fpsidyn = DRT::Problem::Instance()->FPSIDynamicParams();
    int coupling = DRT::INPUT::IntegralValue<int>(fpsidyn,"COUPALGO");
    if (coupling == fpsi_monolithic_plain)
    {
      dirichletcond = false;
    }
    else if (coupling == partitioned)
    {
      dserror("partitioned fpsi solution scheme has not been implemented yet.");
    }
  }

  if (probtype == prb_freesurf)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
         coupling == fsi_iter_monolithicstructuresplit or
         coupling == fsi_iter_constr_monolithicfluidsplit or
         coupling == fsi_iter_constr_monolithicstructuresplit or
         coupling == fsi_iter_lung_monolithicfluidsplit or
         coupling == fsi_iter_lung_monolithicstructuresplit or
         coupling == fsi_iter_mortar_monolithicstructuresplit or
         coupling == fsi_iter_mortar_monolithicfluidsplit)
    {
      dirichletcond = false;
    }
  }

  int aletype = DRT::INPUT::IntegralValue<int>(*adyn,"ALE_TYPE");
  if (aletype == INPAR::ALE::classic_lin)
    ale_ = Teuchos::rcp(new AleLinear(actdis, solver, adyn, output, false, dirichletcond));
  else if (aletype == INPAR::ALE::incr_lin)
    ale_ = Teuchos::rcp(new AleLinear(actdis, solver, adyn, output, true, dirichletcond));
  else if (aletype == INPAR::ALE::laplace)
    ale_ = Teuchos::rcp(new AleLaplace(actdis, solver, adyn, output, true, dirichletcond));
  else if (aletype == INPAR::ALE::springs)
    ale_ = Teuchos::rcp(new AleSprings(actdis, solver, adyn, output, dirichletcond));
  else if (aletype == INPAR::ALE::springs_fixed_ref)
    ale_ = Teuchos::rcp(new AleSpringsFixedRef(actdis, solver, adyn, output, true, dirichletcond));
  else
    dserror("ale type '%s' unsupported", adyn->get<std::string>("ALE_TYPE").c_str());
}
