/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for arterial network stationary formulation.


\level 3

*----------------------------------------------------------------------*/


#include "baci_art_net_impl_stationary.H"

#include "baci_adapter_scatra_base_algorithm.H"
#include "baci_art_net_artery_ele_action.H"
#include "baci_art_net_artery_resulttest.H"
#include "baci_global_data.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_linalg_utils_sparse_algebra_print.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_mat_cnst_1d_art.H"
#include "baci_scatra_resulttest.H"
#include "baci_scatra_timint_implicit.H"
#include "baci_utils_function.H"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN



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
    const Teuchos::ParameterList& artparams, IO::DiscretizationWriter& output)
    : TimInt(actdis, linsolvernumber, probparams, artparams, output)
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
  sysmat_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*(discret_->DofRowMap()), 3, false, true));

  // right hand side vector
  rhs_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new CORE::LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    discret_->EvaluateDirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);  // just in case of change
  }

  // the vector containing body and surface forces
  neumann_loads_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1
  pressurenp_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  pressureincnp_ = CORE::LINALG::CreateVector(*dofrowmap, true);

  // for output of volumetric flow
  ele_volflow_ = CORE::LINALG::CreateVector(*discret_->ElementRowMap());

  // for output of element radius
  ele_radius_ = CORE::LINALG::CreateVector(*discret_->ElementRowMap());

  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(INPUT::IntegralValue<INPAR::ARTDYN::InitialField>(arteryparams, "INITIALFIELD"),
      arteryparams.get<int>("INITFUNCNO"));


  if (solvescatra_)
  {
    const Teuchos::ParameterList& myscatraparams =
        GLOBAL::Problem::Instance()->ScalarTransportDynamicParams();
    if (INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(myscatraparams, "VELOCITYFIELD") !=
        INPAR::SCATRA::velocity_zero)
      dserror("set your velocity field to zero!");
    // construct the scatra problem
    scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(globaltimeparams, myscatraparams,
        GLOBAL::Problem::Instance()->SolverParams(linsolvernumber_), scatra_disname, false));

    // initialize the base algo.
    // scatra time integrator is initialized inside.
    scatra_->Init();

    // only now we must call Setup() on the scatra time integrator.
    // all objects relying on the parallel distribution are
    // created and pointers are set.
    // calls Setup() on the scatra time integrator inside.
    scatra_->ScaTraField()->Setup();
  }
}



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
  CORE::LINALG::ApplyDirichletToSystem(
      *sysmat_, *pressureincnp_, *rhs_, *zeros_, *(dbcmaps_->CondMap()));
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
  CORE::LINALG::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver_->Solve(sysmat_->EpetraOperator(), pressureincnp_, rhs_, solver_params);
  // note: incremental form since rhs-coupling with poromultielastscatra-framework might be
  //       nonlinear
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
 | reset artery diameter of previous time step         kremheller 11/20 |
 *----------------------------------------------------------------------*/
void ART::ArtNetImplStationary::ResetArteryDiamPreviousTimeStep()
{
  // set the diameter in material
  for (int i = 0; i < discret_->NumMyColElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = discret_->lColElement(i);

    // get the artery-material
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(actele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");

    const double diam = arterymat->Diam();
    arterymat->SetDiamPreviousTimeStep(diam);
  }
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
  // reset the artery diameter of the previous time step
  ResetArteryDiamPreviousTimeStep();

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
    // set artery diameter of previous time step to intial artery diameter
    ResetArteryDiamPreviousTimeStep();
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
    // and element radius
    if (step_ == upres_ or step_ == 0)
    {
      output_.WriteElementData(true);
    }
    // for variable radius, we need the output of the radius at every time step
    OutputRadius();

    // "pressure in the arteries" vector
    output_.WriteVector("one_d_artery_pressure", pressurenp_);

    // output of flow
    OutputFlow();

    if (solvescatra_) scatra_->ScaTraField()->CheckAndWriteOutputAndRestart();
  }

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of element-based radius                       kremheller 07/19|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::OutputRadius()
{
  // loop over row elements
  const int numrowele = discret_->NumMyRowElements();
  for (int i = 0; i < numrowele; ++i)
  {
    DRT::Element* actele = discret_->lRowElement(i);
    // cast the material to artery material material
    const Teuchos::RCP<const MAT::Cnst_1d_art>& arterymat =
        Teuchos::rcp_dynamic_cast<const MAT::Cnst_1d_art>(actele->Material());
    if (arterymat == Teuchos::null)
      dserror("cast to MAT::Cnst_1d_art failed during output of radius!");
    const double radius = arterymat->Diam() / 2.0;
    ele_radius_->ReplaceGlobalValue(actele->Id(), 0, radius);
  }

  // write the output
  output_.WriteVector("ele_radius", ele_radius_, IO::elementvector);

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | output of element volumetric flow                    kremheller 09/19|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::ArtNetImplStationary::OutputFlow()
{
  CORE::LINALG::SerialDenseMatrix dummyMat;
  CORE::LINALG::SerialDenseVector dummyVec;

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState(0, "pressurenp", pressurenp_);

  // enough to loop over row nodes since element-based quantity
  for (int i = 0; i < discret_->NumMyRowElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = discret_->lRowElement(i);

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action", ARTERY::calc_flow_pressurebased);

    DRT::Element::LocationArray la(discret_->NumDofSets());
    actele->LocationVector(*discret_, la, false);
    CORE::LINALG::SerialDenseVector flowVec(1);

    actele->Evaluate(p, *discret_, la, dummyMat, dummyMat, flowVec, dummyVec, dummyVec);

    int err = ele_volflow_->ReplaceMyValue(i, 0, flowVec(0));
    if (err != 0) dserror("ReplaceMyValue failed with error code %d!", err);
  }

  // write the output
  output_.WriteVector("ele_volflow", ele_volflow_, IO::elementvector);

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
  GLOBAL::Problem::Instance()->AddFieldTest(resulttest);
  if (solvescatra_)
  {
    GLOBAL::Problem::Instance()->AddFieldTest(scatra_->CreateScaTraFieldTest());
  }
  GLOBAL::Problem::Instance()->TestAll(discret_->Comm());
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
  IO::DiscretizationReader reader(discret_, GLOBAL::Problem::Instance()->InputControlFile(), step);

  if (step != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(pressurenp_, "one_d_artery_pressure");

  // read restart for diameter of previous time step
  reader.ReadVector(ele_radius_, "ele_radius");
  Teuchos::RCP<Epetra_Vector> ele_radius_col =
      CORE::LINALG::CreateVector(*discret_->ElementColMap(), true);
  CORE::LINALG::Export(*ele_radius_, *ele_radius_col);

  // set the diameter in material
  for (int i = 0; i < discret_->NumMyColElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = discret_->lColElement(i);

    // get the artery-material
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(actele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");

    const double diam = 2.0 * (*ele_radius_col)[i];

    // reset (if element is collapsed in previous step, set to zero)
    arterymat->SetDiamPreviousTimeStep(diam);
    arterymat->SetDiam(diam);
    if (diam < arterymat->CollapseThreshold()) arterymat->SetDiam(0.0);
  }

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
          double initialval = GLOBAL::Problem::Instance()
                                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                                  .Evaluate(lnode->X().data(), time_, k);
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

}  // ArtNetImplStationary::SetInitialField

BACI_NAMESPACE_CLOSE
