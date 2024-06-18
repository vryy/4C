/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for arterial network stationary formulation.


\level 3

*----------------------------------------------------------------------*/


#include "4C_art_net_impl_stationary.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_art_net_artery_ele_action.hpp"
#include "4C_art_net_artery_resulttest.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_cnst_1d_art.hpp"
#include "4C_scatra_resulttest.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

Arteries::ArtNetImplStationary::ArtNetImplStationary(Teuchos::RCP<Core::FE::Discretization> actdis,
    const int linsolvernumber, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& artparams, Core::IO::DiscretizationWriter& output)
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
void Arteries::ArtNetImplStationary::init(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname)
{
  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + initialization");

  if (coupledTo3D_)
    FOUR_C_THROW("this type of coupling is only available for explicit time integration");

  // call base class
  TimInt::init(globaltimeparams, arteryparams, scatra_disname);

  // ensure that degrees of freedom in the discretization have been set
  if ((not discret_->Filled()) or (not discret_->HaveDofs())) discret_->fill_complete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // -------------------------------------------------------------------
  // create empty system matrix (6 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  sysmat_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(discret_->dof_row_map()), 3, false, true));

  // right hand side vector
  rhs_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    eleparams.set<const Core::UTILS::FunctionManager*>(
        "function_manager", &Global::Problem::Instance()->FunctionManager());
    discret_->evaluate_dirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);  // just in case of change
  }

  // the vector containing body and surface forces
  neumann_loads_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1
  pressurenp_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  pressureincnp_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // for output of volumetric flow
  ele_volflow_ = Core::LinAlg::CreateVector(*discret_->ElementRowMap());

  // for output of element radius
  ele_radius_ = Core::LinAlg::CreateVector(*discret_->ElementRowMap());

  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(
      Core::UTILS::IntegralValue<Inpar::ArtDyn::InitialField>(arteryparams, "INITIALFIELD"),
      arteryparams.get<int>("INITFUNCNO"));


  if (solvescatra_)
  {
    const Teuchos::ParameterList& myscatraparams =
        Global::Problem::Instance()->scalar_transport_dynamic_params();
    if (Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(myscatraparams, "VELOCITYFIELD") !=
        Inpar::ScaTra::velocity_zero)
      FOUR_C_THROW("set your velocity field to zero!");
    // construct the scatra problem
    scatra_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(globaltimeparams, myscatraparams,
        Global::Problem::Instance()->SolverParams(linsolvernumber_), scatra_disname, false));

    // initialize the base algo.
    // scatra time integrator is initialized inside.
    scatra_->init();

    // only now we must call setup() on the scatra time integrator.
    // all objects relying on the parallel distribution are
    // created and pointers are set.
    // calls setup() on the scatra time integrator inside.
    scatra_->ScaTraField()->setup();
  }
}



/*----------------------------------------------------------------------*
 | (Linear) Solve.                                                      |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::ArtNetImplStationary::Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  // time measurement: initialization
  TEUCHOS_FUNC_TIME_MONITOR(" + solve");

  if (coupledTo3D_)
    FOUR_C_THROW("this type of coupling is only available for implicit time integration");

  // call elements to calculate system matrix and rhs and assemble
  assemble_mat_and_rhs();

  // Prepare Linear Solve (Apply DBC)
  PrepareLinearSolve();

  // solve linear system of equations
  linear_solve();
}


/*----------------------------------------------------------------------*
 | (Linear) Solve for ScaTra.                                           |
 |                                                      kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::ArtNetImplStationary::SolveScatra()
{
  // print user info
  if (discretization()->Comm().MyPID() == 0)
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
    FOUR_C_THROW("this type of coupling is only available for explicit time integration");

  // provide scatra discretization with fluid primary variable field
  scatra_->ScaTraField()->discretization()->set_state(1, "one_d_artery_pressure", pressurenp_);
  scatra_->ScaTraField()->prepare_time_step();

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
void Arteries::ArtNetImplStationary::PrepareLinearSolve()
{
  // apply map: rhs = pressurenp_
  Core::LinAlg::apply_dirichlet_to_system(
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
void Arteries::ArtNetImplStationary::assemble_mat_and_rhs()
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
  eleparams.set<int>("action", Arteries::calc_sys_matrix_rhs);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->set_state(0, "pressurenp", pressurenp_);

  // call standard loop over all elements
  discret_->evaluate(eleparams, sysmat_, rhs_);
  discret_->ClearState();

  // potential addition of Neumann terms
  add_neumann_to_residual();

  // finalize the complete matrix
  sysmat_->Complete();

  // end time measurement for element
  double mydtele = Teuchos::Time::wallTime() - tcpuele;
  discret_->Comm().MaxAll(&mydtele, &dtele_, 1);

}  // ArtNetExplicitTimeInt::assemble_mat_and_rhs

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
void Arteries::ArtNetImplStationary::linear_solve()
{
  // time measurement: solver
  TEUCHOS_FUNC_TIME_MONITOR("      + solver");

  // get cpu time
  const double tcpusolve = Teuchos::Time::wallTime();

  // linear solve
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver_->Solve(sysmat_->EpetraOperator(), pressureincnp_, rhs_, solver_params);
  // note: incremental form since rhs-coupling with poromultielastscatra-framework might be
  //       nonlinear
  pressurenp_->Update(1.0, *pressureincnp_, 1.0);

  // end time measurement for solver
  double mydtsolve = Teuchos::Time::wallTime() - tcpusolve;
  discret_->Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

}  // ArtNetImplStationary::linear_solve

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Prepare time step (Apply DBC and Neumann)            kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::ArtNetImplStationary::prepare_time_step()
{
  // call base class
  Arteries::TimInt::prepare_time_step();

  // Apply DBC
  apply_dirichlet_bc();

  // Apply Neumann
  apply_neumann_bc(neumann_loads_);
}

/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}   kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArtNetImplStationary::apply_dirichlet_bc()
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("      + apply dirich cond.");

  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time_);  // actual time t_{n+1}
  p.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::Instance()->FunctionManager());

  // Dirichlet values
  // \c  pressurenp_ then also holds prescribed new Dirichlet values
  discret_->ClearState();
  discret_->evaluate_dirichlet(
      p, pressurenp_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
  discret_->ClearState();
}

/*----------------------------------------------------------------------*
 | reset artery diameter of previous time step         kremheller 11/20 |
 *----------------------------------------------------------------------*/
void Arteries::ArtNetImplStationary::reset_artery_diam_previous_time_step()
{
  // set the diameter in material
  for (int i = 0; i < discret_->NumMyColElements(); ++i)
  {
    // pointer to current element
    Core::Elements::Element* actele = discret_->lColElement(i);

    // get the artery-material
    Teuchos::RCP<Mat::Cnst1dArt> arterymat =
        Teuchos::rcp_dynamic_cast<Mat::Cnst1dArt>(actele->Material());
    if (arterymat == Teuchos::null) FOUR_C_THROW("cast to artery material failed");

    const double diam = arterymat->Diam();
    arterymat->set_diam_previous_time_step(diam);
  }
}

/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions                kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArtNetImplStationary::apply_neumann_bc(
    const Teuchos::RCP<Epetra_Vector>& neumann_loads)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // create parameter list
  Teuchos::ParameterList condparams;
  condparams.set("total time", time_);
  condparams.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::Instance()->FunctionManager());

  // evaluate Neumann boundary conditions
  discret_->evaluate_neumann(condparams, *neumann_loads);
  discret_->ClearState();

  return;
}  // ArtNetImplStationary::apply_neumann_bc

/*----------------------------------------------------------------------*
 | add actual Neumann loads                            kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArtNetImplStationary::add_neumann_to_residual()
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
void Arteries::ArtNetImplStationary::TimeUpdate()
{
  // reset the artery diameter of the previous time step
  reset_artery_diam_previous_time_step();

  if (solvescatra_)
  {
    scatra_->ScaTraField()->Update();
    scatra_->ScaTraField()->evaluate_error_compared_to_analytical_sol();
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
void Arteries::ArtNetImplStationary::prepare_time_loop()
{
  // call base class
  Arteries::TimInt::prepare_time_loop();

  // provide information about initial field (do not do for restarts!)
  if (step_ == 0)
  {
    // set artery diameter of previous time step to intial artery diameter
    reset_artery_diam_previous_time_step();
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
void Arteries::ArtNetImplStationary::Output(
    bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams)
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("             + output of solution");

  // solution output and potentially restart data
  if (DoOutput())
  {
    // step number and time (only after that data output is possible)
    output_.new_step(step_, time_);

    // write domain decomposition for visualization (only once at step "upres"!)
    // and element radius
    if (step_ == upres_ or step_ == 0)
    {
      output_.write_element_data(true);
    }
    // for variable radius, we need the output of the radius at every time step
    OutputRadius();

    // "pressure in the arteries" vector
    output_.write_vector("one_d_artery_pressure", pressurenp_);

    // output of flow
    OutputFlow();

    if (solvescatra_) scatra_->ScaTraField()->check_and_write_output_and_restart();
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
void Arteries::ArtNetImplStationary::OutputRadius()
{
  // loop over row elements
  const int numrowele = discret_->NumMyRowElements();
  for (int i = 0; i < numrowele; ++i)
  {
    Core::Elements::Element* actele = discret_->lRowElement(i);
    // cast the material to artery material material
    const Teuchos::RCP<const Mat::Cnst1dArt>& arterymat =
        Teuchos::rcp_dynamic_cast<const Mat::Cnst1dArt>(actele->Material());
    if (arterymat == Teuchos::null)
      FOUR_C_THROW("cast to Mat::Cnst1dArt failed during output of radius!");
    const double radius = arterymat->Diam() / 2.0;
    ele_radius_->ReplaceGlobalValue(actele->Id(), 0, radius);
  }

  // write the output
  output_.write_vector("ele_radius", ele_radius_, Core::IO::elementvector);

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
void Arteries::ArtNetImplStationary::OutputFlow()
{
  Core::LinAlg::SerialDenseMatrix dummyMat;
  Core::LinAlg::SerialDenseVector dummyVec;

  // set vector values needed by elements
  discret_->ClearState();
  discret_->set_state(0, "pressurenp", pressurenp_);

  // enough to loop over row nodes since element-based quantity
  for (int i = 0; i < discret_->NumMyRowElements(); ++i)
  {
    // pointer to current element
    Core::Elements::Element* actele = discret_->lRowElement(i);

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action", Arteries::calc_flow_pressurebased);

    Core::Elements::Element::LocationArray la(discret_->NumDofSets());
    actele->LocationVector(*discret_, la, false);
    Core::LinAlg::SerialDenseVector flowVec(1);

    actele->evaluate(p, *discret_, la, dummyMat, dummyMat, flowVec, dummyVec, dummyVec);

    int err = ele_volflow_->ReplaceMyValue(i, 0, flowVec(0));
    if (err != 0) FOUR_C_THROW("ReplaceMyValue failed with error code %d!", err);
  }

  // write the output
  output_.write_vector("ele_volflow", ele_volflow_, Core::IO::elementvector);

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
void Arteries::ArtNetImplStationary::TestResults()
{
  Teuchos::RCP<Core::UTILS::ResultTest> resulttest = CreateFieldTest();
  Global::Problem::Instance()->AddFieldTest(resulttest);
  if (solvescatra_)
  {
    Global::Problem::Instance()->AddFieldTest(scatra_->create_sca_tra_field_test());
  }
  Global::Problem::Instance()->TestAll(discret_->Comm());
}

/*----------------------------------------------------------------------*
 | create result test for this field                   kremheller 03/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> Arteries::ArtNetImplStationary::CreateFieldTest()
{
  return Teuchos::rcp(new Arteries::ArteryResultTest(*(this)));
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | read_restart (public)                                 kremheller 03/18|
 -----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::ArtNetImplStationary::read_restart(int step, bool coupledTo3D)
{
  coupledTo3D_ = coupledTo3D;
  Core::IO::DiscretizationReader reader(
      discret_, Global::Problem::Instance()->InputControlFile(), step);

  if (step != reader.read_int("step")) FOUR_C_THROW("Time step on file not equal to given step");

  time_ = reader.read_double("time");
  step_ = reader.read_int("step");

  reader.read_vector(pressurenp_, "one_d_artery_pressure");

  // read restart for diameter of previous time step
  reader.read_vector(ele_radius_, "ele_radius");
  Teuchos::RCP<Epetra_Vector> ele_radius_col =
      Core::LinAlg::CreateVector(*discret_->ElementColMap(), true);
  Core::LinAlg::Export(*ele_radius_, *ele_radius_col);

  // set the diameter in material
  for (int i = 0; i < discret_->NumMyColElements(); ++i)
  {
    // pointer to current element
    Core::Elements::Element* actele = discret_->lColElement(i);

    // get the artery-material
    Teuchos::RCP<Mat::Cnst1dArt> arterymat =
        Teuchos::rcp_dynamic_cast<Mat::Cnst1dArt>(actele->Material());
    if (arterymat == Teuchos::null) FOUR_C_THROW("cast to artery material failed");

    const double diam = 2.0 * (*ele_radius_col)[i];

    // reset (if element is collapsed in previous step, set to zero)
    arterymat->set_diam_previous_time_step(diam);
    arterymat->SetDiam(diam);
    if (diam < arterymat->CollapseThreshold()) arterymat->SetDiam(0.0);
  }

  if (solvescatra_)
    // read restart data for scatra field
    scatra_->ScaTraField()->read_restart(step);
}

/*----------------------------------------------------------------------*
 |  set initial field for pressure                     kremheller 04/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArtNetImplStationary::SetInitialField(
    const Inpar::ArtDyn::InitialField init, const int startfuncno)
{
  switch (init)
  {
    case Inpar::ArtDyn::initfield_zero_field:
    {
      pressurenp_->PutScalar(0.0);
      break;
    }
    case Inpar::ArtDyn::initfield_field_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        int numdofs = nodedofset.size();
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          double initialval = Global::Problem::Instance()
                                  ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                                  .evaluate(lnode->X().data(), time_, k);
          int err = pressurenp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }

      break;
    }
    case Inpar::ArtDyn::initfield_field_by_condition:
    {
      // set initial field
      const std::string field = "Artery";
      std::vector<int> localdofs;
      localdofs.push_back(0);

      discret_->evaluate_initial_field(
          Global::Problem::Instance()->FunctionManager(), field, pressurenp_, localdofs);

      break;
    }
    default:
      FOUR_C_THROW("Unknown option for initial field: %d", init);
      break;
  }  // switch(init)

}  // ArtNetImplStationary::SetInitialField

FOUR_C_NAMESPACE_CLOSE
