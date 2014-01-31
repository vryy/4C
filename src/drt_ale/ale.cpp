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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ALE::Ale::Ale(RCP<DRT::Discretization> actdis,
              Teuchos::RCP<LINALG::Solver> solver,
              Teuchos::RCP<Teuchos::ParameterList> params,
              Teuchos::RCP<IO::DiscretizationWriter> output,
              bool dirichletcond)
  : discret_(actdis),
    solver_ (solver),
    params_ (params),
    output_ (output),
    step_(0),
    time_(0.0),
    uprestart_(params->get("write restart every", -1)),
    sysmat_(Teuchos::null)
{
  numstep_ = params_->get<int>("numstep");
  maxtime_ = params_->get<double>("maxtime");
  dt_      = params_->get<double>("dt");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispn_          = LINALG::CreateVector(*dofrowmap,true);
  dispnp_         = LINALG::CreateVector(*dofrowmap,true);
  residual_       = LINALG::CreateVector(*dofrowmap,true);

  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*actdis);

  SetupDBCMapEx(dirichletcond);

  // ensure that the ALE std::string was removed from conditions
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
      const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
      int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");
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
  while (step_ < numstep_-1 and time_ <= maxtime_)
  {
    PrepareTimeStep();
    Solve();
    Update();
    Output();
  }
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

  // Update local coordinate systems (which may be time dependent)
  if (locsysman_ != Teuchos::null) {
    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
    locsysman_->Setup(time_);
    discret_->ClearState();
  }
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

  ApplyDirichletBC(eleparams,dispnp_,Teuchos::null,Teuchos::null,true);

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
  Teuchos::RCP<IO::DiscretizationWriter> output =
    Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the linear solver number
  const int linsolvernumber = adyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for ALE problems. Please set LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList());
  params->set<int>("numstep",    prbdyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", prbdyn.get<double>("MAXTIME"));
  params->set<double>("dt",      prbdyn.get<double>("TIMESTEP"));

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", prbdyn.get<int>("RESTARTEVRY"));

  params->set<int>("ALE_TYPE",DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE"));


  bool dirichletcond = true;
  // what's the current problem type?
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
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

  int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");
  if (aletype==INPAR::ALE::classic_lin)
    ale_ = Teuchos::rcp(new AleLinear(actdis, solver, params, output, false, dirichletcond));
  else if (aletype==INPAR::ALE::incr_lin)
    ale_ = Teuchos::rcp(new AleLinear(actdis, solver, params, output, true , dirichletcond));
  else if (aletype==INPAR::ALE::laplace)
    ale_ = Teuchos::rcp(new AleLaplace(actdis, solver, params, output, true, dirichletcond));
  else if (aletype==INPAR::ALE::springs)
    ale_ = Teuchos::rcp(new AleSprings(actdis, solver, params, output, dirichletcond));
  else if (aletype==INPAR::ALE::springs_fixed_ref)
    ale_ = Teuchos::rcp(new AleSpringsFixedRef(actdis, solver, params, output, true, dirichletcond));
  else
    dserror("ale type '%s' unsupported",adyn.get<std::string>("ALE_TYPE").c_str());
}
