/*!----------------------------------------------------------------------
\file art_net_impl_stationary.cpp
\brief Control routine for arterial network stationary formulation.

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/


#include "art_net_impl_stationary.H"
#include "artery_ele_action.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "artery_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_resulttest.H"
#include <Epetra_Vector.h>



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

ART::ArtNetImplStationary::ArtNetImplStationary(Teuchos::RCP<DRT::Discretization> actdis,
    const int linsolvernumber, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& artparams, FILE* errfile, IO::DiscretizationWriter& output)
    : TimInt(actdis, linsolvernumber, probparams, artparams, errfile, output)
{
  //  exit(1);

}  // ArtNetImplStationary::ArtNetImplStationary



/*----------------------------------------------------------------------*
 | Initialize the time integration.                                     |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::Init(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname)
{
  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  if (coupledTo3D_)
    dserror("this type of coupling is only available for explicit time integration");

  // call base class
  TimInt::Init(globaltimeparams, arteryparams, scatra_disname);

  // ensure that degrees of freedom in the discretization have been set
  if ((not discret_->Filled()) or (not discret_->HaveDofs())) discret_->FillComplete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // -------------------------------------------------------------------
  // create empty system matrix (6 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()), 3, false, true));

  // right hand side vector
  rhs_ = LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap, true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    discret_->EvaluateDirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);  // just in case of change
  }

  // the vector containing body and surface forces
  neumann_loads_ = LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1
  pressurenp_ = LINALG::CreateVector(*dofrowmap, true);
  pressureincnp_ = LINALG::CreateVector(*dofrowmap, true);


  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(
      DRT::INPUT::IntegralValue<INPAR::ARTDYN::InitialField>(arteryparams, "INITIALFIELD"),
      arteryparams.get<int>("INITFUNCNO"));


  if (solvescatra_)
  {
    const Teuchos::ParameterList& myscatraparams =
        DRT::Problem::Instance()->ScalarTransportDynamicParams();
    if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(myscatraparams, "VELOCITYFIELD") !=
        INPAR::SCATRA::velocity_zero)
      dserror("set your velocity field to zero!");
    // scatra problem
    scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

    // initialize the base algo.
    // scatra time integrator is constructed and initialized inside.
    scatra_->Init(globaltimeparams, myscatraparams,
        DRT::Problem::Instance()->SolverParams(linsolvernumber_), scatra_disname, false);

    // only now we must call Setup() on the scatra time integrator.
    // all objects relying on the parallel distribution are
    // created and pointers are set.
    // calls Setup() on the scatra time integrator inside.
    scatra_->ScaTraField()->Setup();
  }
}

/*----------------------------------------------------------------------*
 | Destructor dtor (public)                             kremheller 03/18|
 *----------------------------------------------------------------------*/
ART::ArtNetImplStationary::~ArtNetImplStationary() { return; }

/*----------------------------------------------------------------------*
 | (Linear) Solve.                                                      |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + solve");

  if (coupledTo3D_)
    dserror("this type of coupling is only available for implicit time integration");

  // call elements to calculate system matrix and rhs and assemble
  AssembleMatAndRHS();

  // Prepare Linear Solve (Apply DBC)
  PrepareLinearSolve();

  // solve linear system of equations
  LinearSolve();
}


/*----------------------------------------------------------------------*
 | (Linear) Solve for ScaTra.                                           |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::SolveScatra()
{
  // print user info
  if (Discretization()->Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<      Scalar Transport in 1D Artery Network       >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }

  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + solve scatra");

  if (coupledTo3D_)
    dserror("this type of coupling is only available for explicit time integration");

  // provide scatra discretization with fluid primary variable field
  scatra_->ScaTraField()->Discretization()->SetState(1, "one_d_artery_pressure", pressurenp_);
  scatra_->ScaTraField()->PrepareTimeStep();

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();
}

/*----------------------------------------------------------------------*
 | Prepare Linear Solve (Apply DBC).                                    |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::PrepareLinearSolve()
{
  // apply map: rhs = pressurenp_
  LINALG::ApplyDirichlettoSystem(sysmat_, pressureincnp_, rhs_, zeros_, *(dbcmaps_->CondMap()));

  bool matlab = false;
  if (matlab)
  {
    Teuchos::RCP<const LINALG::SparseMatrix> sparse_matrix =
        Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(sysmat_, true);

    // sparse_matrix
    std::string filename = "../o/mymatrix.dat";
    LINALG::PrintMatrixInMatlabFormat(
        filename, *sparse_matrix->EpetraMatrix());  // *sysmat_->EpetraOperator());
    dserror("exit");
  }
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call elements to calculate system matrix/rhs and assemble            |
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::AssembleMatAndRHS()
{
  dtele_ = 0.0;

  TEUCHOS_FUNC_TIME_MONITOR("      + element calls");


  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // set both system matrix and rhs vector to zero
  sysmat_->Zero();
  rhs_->PutScalar(0.0);

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", ARTERY::calc_sys_matrix_rhs);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0, "pressurenp", pressurenp_);

  // call standard loop over all elements
  discret_->Evaluate(eleparams, sysmat_, rhs_);
  discret_->ClearState();

  // potential addition of Neumann terms
  AddNeumannToResidual();

  // finalize the complete matrix
  sysmat_->Complete();

  // end time measurement for element
  double mydtele = Teuchos::Time::wallTime() - tcpuele;
  discret_->Comm().MaxAll(&mydtele, &dtele_, 1);

}  // ArtNetExplicitTimeInt::AssembleMatAndRHS

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | call linear solver                                                   |
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::LinearSolve()
{
  // time measurement: solver
  TEUCHOS_FUNC_TIME_MONITOR("      + solver");

  // get cpu time
  const double tcpusolve = Teuchos::Time::wallTime();

  // linear solve
  solver_->Solve(sysmat_->EpetraOperator(), pressureincnp_, rhs_, true, 1, Teuchos::null);
  pressurenp_->Update(1.0, *pressureincnp_, 1.0);

  // end time measurement for solver
  double mydtsolve = Teuchos::Time::wallTime() - tcpusolve;
  discret_->Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

}  // ArtNetImplStationary::LinearSolve

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Prepare time step (Apply DBC and Neumann)            kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::PrepareTimeStep()
{
  // call base class
  ART::TimInt::PrepareTimeStep();

  // Apply DBC
  ApplyDirichletBC();

  // Apply Neumann
  ApplyNeumannBC(neumann_loads_);
}

/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}   kremheller 03/18 |
 *----------------------------------------------------------------------*/
void ART::ArtNetImplStationary::ApplyDirichletBC()
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("      + apply dirich cond.");

  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time_);  // actual time t_{n+1}

  // Dirichlet values
  // \c  pressurenp_ then also holds prescribed new Dirichlet values
  discret_->ClearState();
  discret_->EvaluateDirichlet(
      p, pressurenp_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
  discret_->ClearState();
}

/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions                kremheller 03/18 |
 *----------------------------------------------------------------------*/
void ART::ArtNetImplStationary::ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // create parameter list
  Teuchos::ParameterList condparams;
  condparams.set("total time", time_);

  // evaluate Neumann boundary conditions
  discret_->EvaluateNeumann(condparams, *neumann_loads);
  discret_->ClearState();

  return;
}  // ArtNetImplStationary::ApplyNeumannBC

/*----------------------------------------------------------------------*
 | add actual Neumann loads                            kremheller 03/18 |
 *----------------------------------------------------------------------*/
void ART::ArtNetImplStationary::AddNeumannToResidual()
{
  rhs_->Update(1.0, *neumann_loads_, 1.0);
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::TimeUpdate()
{
  // do nothing: stationary

  if (solvescatra_)
  {
    scatra_->ScaTraField()->Update();
    scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
  }

  return;
}  // ArtNetExplicitTimeInt::TimeUpdate

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | prepare the time loop                                kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::PrepareTimeLoop()
{
  // call base class
  ART::TimInt::PrepareTimeLoop();

  // provide information about initial field (do not do for restarts!)
  if (step_ == 0)
  {
    // write out initial state
    Output(false, Teuchos::null);
  }

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of solution vector to binio                   kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::Output(
    bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams)
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("             + output of solution");

  // solution output and potentially restart data
  if (DoOutput())
  {
    // step number and time (only after that data output is possible)
    output_.NewStep(step_, time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    if (step_ == upres_ or step_ == 0) output_.WriteElementData(true);

    // "pressure in the arteries" vector
    output_.WriteVector("one_d_artery_pressure", pressurenp_);

    if (solvescatra_) scatra_->ScaTraField()->Output();
  }

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | test results                                         kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::TestResults()
{
  Teuchos::RCP<DRT::ResultTest> resulttest = CreateFieldTest();
  DRT::Problem::Instance()->AddFieldTest(resulttest);
  if (solvescatra_)
  {
    DRT::Problem::Instance()->AddFieldTest(scatra_->CreateScaTraFieldTest());
  }
  DRT::Problem::Instance()->TestAll(discret_->Comm());
}

/*----------------------------------------------------------------------*
 | create result test for this field                   kremheller 03/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ART::ArtNetImplStationary::CreateFieldTest()
{
  return Teuchos::rcp(new ART::ArteryResultTest(*(this)));
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | ReadRestart (public)                                 kremheller 03/18|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::ReadRestart(int step, bool coupledTo3D)
{
  coupledTo3D_ = coupledTo3D;
  IO::DiscretizationReader reader(discret_, step);

  if (step != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(pressurenp_, "one_d_artery_pressure");
  if (solvescatra_)
    // read restart data for scatra field
    scatra_->ScaTraField()->ReadRestart(step);
}

/*----------------------------------------------------------------------*
 |  set initial field for pressure                     kremheller 04/18 |
 *----------------------------------------------------------------------*/
void ART::ArtNetImplStationary::SetInitialField(
    const INPAR::ARTDYN::InitialField init, const int startfuncno)
{
  switch (init)
  {
    case INPAR::ARTDYN::initfield_zero_field:
    {
      pressurenp_->PutScalar(0.0);
      break;
    }
    case INPAR::ARTDYN::initfield_field_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        int numdofs = nodedofset.size();
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          double initialval =
              DRT::Problem::Instance()->Funct(startfuncno - 1).Evaluate(k, lnode->X(), time_);
          int err = pressurenp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) dserror("dof not on proc");
        }
      }

      break;
    }
    case INPAR::ARTDYN::initfield_field_by_condition:
    {
      // set initial field
      const std::string field = "Artery";
      std::vector<int> localdofs;
      localdofs.push_back(0);
      discret_->EvaluateInitialField(field, pressurenp_, localdofs);

      break;
    }
    default:
      dserror("Unknown option for initial field: %d", init);
      break;
  }  // switch(init)

  return;
}  // ArtNetImplStationary::SetInitialField
