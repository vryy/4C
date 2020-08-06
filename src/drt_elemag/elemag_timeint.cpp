/*----------------------------------------------------------------------*/
/*! \file

\brief Base class functions for time integration of electromagnetics

\level 3

*/
/*----------------------------------------------------------------------*/

#include "elemag_timeint.H"
#include "elemag_ele_action.H"
#include "elemag_resulttest.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_solver.H"
#include "../drt_mat/electromagnetic.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
ELEMAG::ElemagTimeInt::ElemagTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG> &actdis,
    const Teuchos::RCP<LINALG::Solver> &solver, const Teuchos::RCP<Teuchos::ParameterList> &params,
    const Teuchos::RCP<IO::DiscretizationWriter> &output)
    : discret_(actdis),
      solver_(solver),
      params_(params),
      output_(output),
      elemagdyna_(DRT::INPUT::IntegralValue<INPAR::ELEMAG::DynamicType>(*params_, "TIMEINT")),
      myrank_(actdis->Comm().MyPID()),
      time_(0.0),
      step_(0),
      restart_(params_->get<int>("restart")),
      maxtime_(params_->get<double>("MAXTIME")),
      stepmax_(params_->get<int>("NUMSTEP")),
      uprestart_(params_->get<int>("RESTARTEVRY", -1)),
      upres_(params_->get<int>("RESULTSEVRY", -1)),
      numdim_(DRT::Problem::Instance()->NDim()),
      dtp_(params_->get<double>("TIMESTEP")),
      tau_(params_->get<double>("TAU")),
      dtele_(0.0),
      dtsolve_(0.0),
      calcerr_(DRT::INPUT::IntegralValue<bool>(*params_, "CALCERR")),
      errfunct_(params_->get<int>("ERRORFUNCNO", -1)),
      sourcefuncno_(params_->get<int>("SOURCEFUNCNO", -1))
{
  // constructor supposed to be empty!

}  // ElemagTimeInt

/*----------------------------------------------------------------------*
 |  Desctructor (public)                               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
ELEMAG::ElemagTimeInt::~ElemagTimeInt() {}


/*----------------------------------------------------------------------*
 |  initialization routine (public)                    berardocco 02/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Init()
{
  // get dof row map
  const Epetra_Map *dofrowmap = discret_->DofRowMap();

  // check time-step length
  if (dtp_ <= 0.0) dserror("Zero or negative time-step length!");

  // Nodevectors for the output
  electric = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  magnetic = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  trace = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  conductivity = Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap()));
  permittivity = Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap()));
  permeability = Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap()));

  // create vector of zeros to be used for enforcing zero Dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap, true);

  trace_ = LINALG::CreateVector(*dofrowmap, true);
  // Map of the dirichlet conditions
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  // Why is this in a new scope?
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);

    // Evaluation of the dirichlet conditions (why is it here and also later?)
    discret_->EvaluateDirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);

    // Initialize elements
    ElementsInit();
  }

  // create system matrix and set to zero
  // the 108 comes from line 282 of /drt_fluid/fluidimplicitintegration.cpp
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 108, false, true));
  // Is it possible to avoid this passage? It is a sparse matrix so it should
  // only contain non-zero entries that have to be initialized
  sysmat_->Zero();

  // create residual vector
  residual_ = LINALG::CreateVector(*dofrowmap, true);

  // write mesh
  output_->WriteMesh(0, 0.0);

  return;
}  // Init


/*----------------------------------------------------------------------*
 |  Print information to screen (public)               berardocco 02/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::PrintInformationToScreen()
{
  if (!myrank_)
  {
    std::cout << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout << "INTEGRATION OF AN ELECTROMAGNETIC PROBLEM USING HDG" << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout << "DISCRETIZATION PARAMETERS:" << std::endl;
    std::cout << "number of DoF sets          " << discret_->NumDofSets() << std::endl;
    std::cout << "number of nodes             " << discret_->NumGlobalNodes() << std::endl;
    std::cout << "number of elements          " << discret_->NumGlobalElements() << std::endl;
    std::cout << "number of faces             " << discret_->NumGlobalFaces() << std::endl;
    std::cout << "number of trace unknowns    " << discret_->DofRowMap(0)->NumGlobalElements()
              << std::endl;
    std::cout << "number of interior unknowns " << discret_->DofRowMap(1)->NumGlobalElements()
              << std::endl;
    std::cout << std::endl;
    std::cout << "SIMULATION PARAMETERS: " << std::endl;
    std::cout << "time step size              " << dtp_ << std::endl;
    std::cout << "time integration scheme     " << Name() << std::endl;
    std::cout << "tau                         " << tau_ << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout << std::endl;
  }
  return;
}  //  PrintInformationToScreen


/*----------------------------------------------------------------------*
 |  Time integration (public)                          berardocco 03/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Integrate()
{
  // Fancy printing
  if (!myrank_)
  {
    std::cout << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout
        << "                              INTEGRATION                                          "
        << std::endl;
  }
  // time measurement: integration
  TEUCHOS_FUNC_TIME_MONITOR("ELEMAG::ElemagTimeInt::Integrate");

  InitializeAlgorithm();

  // call elements to calculate system matrix/rhs and assemble
  AssembleMatAndRHS();

  // Compute matrix for ABC boundary conditions
  ComputeSilverMueller(false);

  // time loop
  while (step_ < stepmax_ and time_ < maxtime_)
  {
    // increment time and step
    IncrementTimeAndStep();

    ApplyDirichletToSystem();
    ComputeSilverMueller(true);
    Solve();

    UpdateInteriorVariablesAndAssembleRHS();

    // The output to file only once in a while
    if (step_ % upres_ == 0)
    {
      Output();
      // Output to screen
      OutputToScreen();
    }
  }  // while (step_<stepmax_ and time_<maxtime_)

  if (!myrank_)
  {
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout << std::endl;
  }

  return;
}  // Integrate

void ELEMAG::ElemagTimeInt::ElementsInit()
{
  // Initializing siome vectors and parameters
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;

  // loop over all elements on the processor
  DRT::Element::LocationArray la(2);
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // Selecting the elements
    DRT::Element *ele = discret_->lColElement(el);

    // This function is a void function and therefore the input goes to the vector "la"
    ele->LocationVector(*discret_, la, true);

    // This is needed to store the dirichlet dofs in the
    if (std::find(la[0].lmdirich_.begin(), la[0].lmdirich_.end(), 1) != la[0].lmdirich_.end())
      initParams.set<std::vector<int> *>("dirichdof", &la[0].lmdirich_);
    else
      initParams.remove("dirichdof", false);

    initParams.set<int>("action", ELEMAG::ele_init);
    initParams.set<INPAR::ELEMAG::DynamicType>("dyna", elemagdyna_);
    Epetra_SerialDenseVector elevec1, elevec2, elevec3;
    Epetra_SerialDenseMatrix elemat1, elemat2;
    ele->Evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Set initial field by given function (public)       berardocco 03/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::SetInitialField(
    const INPAR::ELEMAG::InitialField init, Teuchos::ParameterList &start_params)
{
// time measurement: SetInitialField just in the debug phase
#ifdef DEBUG
  TEUCHOS_FUNC_TIME_MONITOR("ELEMAG::ElemagTimeInt::SetInitialField");
#endif

  // Fancy printing
  if (!myrank_)
  {
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout
        << "                          INITIALIZATION OF THE FIELD                              "
        << std::endl;
  }
  // Core of the routine
  switch (init)
  {
    case INPAR::ELEMAG::initfield_zero_field:
    {
      // Fancy printing to help debugging
      if (!myrank_)
      {
        std::cout << "Initializing a zero field." << std::endl;
      }

      break;
    }
    case INPAR::ELEMAG::initfield_field_by_function:
    {
      if (!myrank_)
      {
        std::cout << "Initializing field as specified by STARTFUNCNO "
                  << start_params.get<int>("startfuncno") << std::endl;
      }

      // Initializing siome vectors and parameters
      Epetra_SerialDenseVector elevec1, elevec2, elevec3;
      Epetra_SerialDenseMatrix elemat1, elemat2;
      start_params.set<int>("action", ELEMAG::project_field);
      start_params.set<double>("time", time_);
      start_params.set<double>("dt", dtp_);
      start_params.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);
      // loop over all elements on the processor
      DRT::Element::LocationArray la(2);
      for (int el = 0; el < discret_->NumMyColElements(); ++el)
      {
        // Selecting the elements
        DRT::Element *ele = discret_->lColElement(el);

        // This function is a void function and therefore the input goes to the vector "la"
        ele->LocationVector(*discret_, la, false);

        // Reshaping the vectors
        if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
          elevec1.Shape(la[0].lm_.size(), 1);
        if (elevec2.M() != discret_->NumDof(1, ele)) elevec2.Shape(discret_->NumDof(1, ele), 1);
        ele->Evaluate(
            start_params, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
      }
      break;
    }  // case INPAR::ELEMAG::initfield_field_by_function
    default:
      dserror("Option for initial field not implemented: %d", init);
      break;
  }  // switch(init)

  // Output of the initial condition
  Output();
  if (!myrank_)
  {
    std::cout << "Initial condition projected." << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
  }
  return;
}  // SetInitialField

/*----------------------------------------------------------------------*
 |  Compute error routine (public)                     berardocco 08/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_SerialDenseVector> ELEMAG::ElemagTimeInt::ComputeError()
{
  // Initializing siome vectors and parameters
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList params;
  params.set<int>("action", ELEMAG::compute_error);
  params.set<int>("funcno", errfunct_);
  params.set<double>("time", time_);
  params.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);

  Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(8));

  // call loop over elements (assemble nothing)
  discret_->EvaluateScalars(params, errors);
  discret_->ClearState(true);

  return errors;
}

void ELEMAG::ElemagTimeInt::PrintErrors(Teuchos::RCP<Epetra_SerialDenseVector> &errors)
{
  if (myrank_ == 0)
  {
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
    std::cout << "---------------------------- Result error wrt FUNCT " << errfunct_
              << " -----------------------------" << std::endl;
    std::cout << "ABSOLUTE ERROR:" << std::endl;
    std::cout << "Electric L2-error: " << std::sqrt((*errors)[0]) << std::endl;
    std::cout << "Magnetic L2-error: " << std::sqrt((*errors)[2]) << std::endl;
    std::cout << "\nRELATIVE ERROR:" << std::endl;
    if ((*errors)[1] == 0 && (*errors)[3] == 0)
      std::cout << "Impossible to compute relative errors. The L2-norm of the analytical "
                   "solution is zero, resulting in a division by zero."
                << std::endl;
    else
    {
      if ((*errors)[1] != 0)
        std::cout << "Electric relative L2-error: " << std::sqrt((*errors)[0] / (*errors)[1])
                  << std::endl;
      else
        std::cout
            << "Impossible to compute the electric relative error. The L2-norm of the analytical "
               "solution is zero, resulting in a division by zero."
            << std::endl;
      if ((*errors)[3] != 0)
        std::cout << "Magnetic relative L2-error: " << std::sqrt((*errors)[2] / (*errors)[3])
                  << std::endl;
      else
        std::cout
            << "Impossible to compute the magnetic relative error. The L2-norm of the analytical "
               "solution is zero, resulting in a division by zero."
            << std::endl;
    }
    std::cout << "\nHDIV ERROR:" << std::endl;
    std::cout << "Electric Hdiv-error: " << std::sqrt((*errors)[4]) << std::endl;
    std::cout << "Magnetic Hdiv-error: " << std::sqrt((*errors)[6]) << std::endl;
    std::cout << "\nHCURL ERROR:" << std::endl;
    std::cout << "Electric Hcurl-error: " << std::sqrt((*errors)[5]) << std::endl;
    std::cout << "Magnetic Hcurl-error: " << std::sqrt((*errors)[7]) << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Project interior variables for testing purposes     berardocco 07/18|
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::ProjectFieldTest(const int startfuncno)
{
  // Initializing siome vectors and parameters
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;
  initParams.set<int>("action", ELEMAG::project_field_test);
  initParams.set("startfuncno", startfuncno);
  initParams.set<double>("time", time_);
  initParams.set<bool>("padaptivity", false);
  initParams.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);

  // loop over all elements on the processor
  DRT::Element::LocationArray la(2);
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // Selecting the elements
    DRT::Element *ele = discret_->lColElement(el);

    // This function is a void function and therefore the input goes to the vector "la"
    ele->LocationVector(*discret_, la, false);

    // Reshaping the vectors
    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1, ele)) elevec2.Shape(discret_->NumDof(1, ele), 1);
    ele->Evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Project trace variable for testing purposes        berardocco 07/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::ProjectFieldTestTrace(const int startfuncno)
{
  // This map contains the trace values
  const Epetra_Map *dofrowmap = discret_->DofRowMap();

  // Initializing siome vectors and parameters
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;
  initParams.set<int>("action", ELEMAG::project_field_test_trace);
  initParams.set("startfuncno", startfuncno);
  initParams.set<double>("time", time_);
  initParams.set<bool>("padaptivity", false);
  initParams.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);
  // loop over all elements on the processor
  DRT::Element::LocationArray la(2);
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // Selecting the elements
    DRT::Element *ele = discret_->lColElement(el);

    // This function is a void function and therefore the input goes to the vector "la"
    ele->LocationVector(*discret_, la, false);

    // Reshaping the vectors
    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1, ele)) elevec2.Shape(discret_->NumDof(1, ele), 1);
    ele->Evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
    // now fill the ele vector into the discretization
    for (unsigned int i = 0; i < la[0].lm_.size(); ++i)
    {
      // In a serial program the global and local ids are the same because there is no need to
      // communicate
      const int lid = dofrowmap->LID(la[0].lm_[i]);
      if (lid >= 0) (*trace_)[lid] = elevec1(i);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Run the first step with BDF1 (public)              berardocco 10/19 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::InitializeAlgorithm()
{
  // In case We don't have BDF1 we initialize the method with BDF1 and then go back
  if (elemagdyna_ == INPAR::ELEMAG::DynamicType::elemag_bdf2 && restart_ == 0)
  {
    INPAR::ELEMAG::DynamicType temp_dyna = elemagdyna_;
    // First step with a BDF1
    elemagdyna_ = INPAR::ELEMAG::DynamicType::elemag_bdf1;

    // call elements to calculate system matrix/rhs and assemble
    AssembleMatAndRHS();

    // Compute matrix for ABC boundary conditions
    ComputeSilverMueller(false);

    // increment time and step
    IncrementTimeAndStep();

    ApplyDirichletToSystem();
    ComputeSilverMueller(true);
    Solve();

    UpdateInteriorVariablesAndAssembleRHS();

    // The output to file only once in a while
    if (step_ % upres_ == 0)
    {
      Output();
      // Output to screen
      OutputToScreen();
    }

    // Set the dynamics back to the original
    elemagdyna_ = temp_dyna;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix and right-hand side (public)       berardocco 04/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::AssembleMatAndRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR("ELEMAG::ElemagTimeInt::AssembleMatAndRHS");

  // Fancy printing to help debugging
  if (!myrank_)
  {
    std::cout << "Creating system matrix elementwise and assembling in the final system matrix."
              << std::endl;
  }

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual and sysmat
  residual_->Scale(0.0);
  sysmat_->Zero();

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------

  // set general vector values needed by elements
  discret_->ClearState(true);

  discret_->SetState("trace", trace_);

  // set time step size
  eleparams.set<double>("dt", dtp_);
  eleparams.set<double>("tau", tau_);

  bool resonly = false;

  // set information needed by the elements
  eleparams.set<int>("sourcefuncno", sourcefuncno_);
  eleparams.set<bool>("resonly", resonly);
  eleparams.set<int>("action", ELEMAG::calc_systemmat_and_residual);
  eleparams.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);
  eleparams.set<double>("time", time_);
  eleparams.set<double>("timep", time_ + dtp_);
  eleparams.set<int>("step", step_);

  // The evaluation of the discretization have to happen before Complete() is
  // called or after an UnComplete() call has been made.
  discret_->Evaluate(eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null);
  discret_->ClearState(true);
  sysmat_->Complete();

  return;
}  // AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | Updates interior variables and compute RHS (public)  berardocco 06/18|
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::UpdateInteriorVariablesAndAssembleRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR("ELEMAG::ElemagTimeInt::UpdateInteriorVariablesAndAssembleRHS");

  // create parameterlist
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("sourcefuncno", sourcefuncno_);
  eleparams.set<double>("dt", dtp_);
  eleparams.set<double>("tau", tau_);
  eleparams.set<double>("time", time_);
  eleparams.set<double>("timep", time_ + dtp_);
  eleparams.set<int>("action", ELEMAG::update_secondary_solution_and_calc_residual);
  eleparams.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);
  eleparams.set<int>("step", step_);
  eleparams.set<bool>("resonly", true);

  residual_->Scale(0.0);
  discret_->SetState("trace", trace_);

  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, residual_, Teuchos::null, Teuchos::null);

  discret_->ClearState(true);


  return;
}  // UpdateInteriorVariablesAndAssembleRHS

/*----------------------------------------------------------------------*
 |  Apply Dirichlet b.c. to system (public)            gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::ApplyDirichletToSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
  Teuchos::ParameterList params;
  params.set<double>("total time", time_);
  discret_->EvaluateDirichlet(
      params, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  LINALG::ApplyDirichlettoSystem(
      sysmat_, trace_, residual_, Teuchos::null, zeros_, *(dbcmaps_->CondMap()));
  return;
}  // ApplyDirichletToSystem

/*----------------------------------------------------------------------*
 |  Compute Silver-Mueller         (public)            berardocco 10/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::ComputeSilverMueller(bool resonly)
{
  TEUCHOS_FUNC_TIME_MONITOR("      + Compute Silver-Mueller BC");

  // absorbing boundary conditions
  std::string condname = "Silver-Mueller";
  std::vector<DRT::Condition *> absorbingBC;
  discret_->GetCondition(condname, absorbingBC);

  // Check if there are Silver-Mueller BC
  if (absorbingBC.size())
  {
    Teuchos::ParameterList eleparams;
    eleparams.set<double>("time", time_);
    eleparams.set<bool>("resonly", resonly);
    eleparams.set<int>("action", ELEMAG::calc_abc);
    // Evaluate the boundary condition
    discret_->EvaluateCondition(
        eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null, condname);
  }

  return;
}  // ApplyDirichletToSystem

/*----------------------------------------------------------------------*
 |  Solve system for trace (public)                    berardocco 06/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Solve()
{
  // This part has only been copied from the fluid part to be able to use algebraic multigrid
  // solvers and has to be checked

  discret_->ComputeNullSpaceIfNecessary(solver_->Params(), true);
  // solve for trace
  solver_->Solve(sysmat_->EpetraOperator(), trace_, residual_, true, false, Teuchos::null);

  return;
}  // Solve


/*----------------------------------------------------------------------*
 |  Output to screen (public)                          berardocco 07/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::OutputToScreen()
{
  if (myrank_ == 0)
  {
    std::cout << "Step: " << step_ << ", time: " << time_ << ", written." << std::endl;
  }
  return;
}  // OutputToScreen

namespace
{
  /*----------------------------------------------------------------------*
  |  Interpolate discontinous values to nodal values     berardocco 03/18 |
  *----------------------------------------------------------------------*/
  // internal helper function for output
  void getNodeVectorsHDG(DRT::Discretization &dis, const Teuchos::RCP<Epetra_Vector> &traceValues,
      const int ndim, Teuchos::RCP<Epetra_MultiVector> &electric,
      Teuchos::RCP<Epetra_MultiVector> &magnetic, Teuchos::RCP<Epetra_MultiVector> &trace,
      Teuchos::RCP<Epetra_Vector> &conductivity, Teuchos::RCP<Epetra_Vector> &permittivity,
      Teuchos::RCP<Epetra_Vector> &permeability)
  {
    // create dofsets for electric and pressure at nodes
    // if there is no pressure vector it means that the vectors have not yet
    // been created and therefor it is needed to create them now.
    if (electric.get() == NULL || electric->GlobalLength() != dis.NumGlobalNodes())
    {
      // The electric is a multivector because it is a vectorial field.
      // The multivector is based on the map of the node
      // owned by the processor. The vectors are zeroed.
      electric.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
      magnetic.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
    }

    // Same for the trace and cell pressure.
    trace.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
    // call element routine for interpolate HDG to elements
    // Here it is used the function that acts in the elements, Evaluate().
    Teuchos::ParameterList params;
    params.set<int>("action", ELEMAG::interpolate_hdg_to_node);
    // Setting a name to the dofs maps
    dis.SetState(0, "trace", traceValues);
    // Declaring all the necessary entry for Evaluate()
    DRT::Element::LocationArray la(2);
    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(dis.NumMyRowNodes());

    // For every element of the processor
    for (int el = 0; el < dis.NumMyColElements(); ++el)
    {
      // Opening the element
      DRT::Element *ele = dis.lColElement(el);

      // Making sure the vector is not a zero dimensional vector and in case
      // resizing it. The interpolVec has to contain all the unknown of the
      // element and therefore has to be carefully sized.
      // The dimensioning has been made as the number of nodes times the number of
      // spatial dimensions (all the unknowns are vectorial fields) times the
      // number of fields that are present:
      // trace
      // Electric field
      // Magnetic field
      ele->LocationVector(dis, la, false);
      if (interpolVec.M() == 0) interpolVec.Resize(ele->NumNode() * 3 * ndim);
      for (int i = 0; i < interpolVec.Length(); i++) interpolVec(i) = 0.0;

      // Interpolating hdg internal values to the node
      ele->Evaluate(params, dis, la[0].lm_, dummyMat, dummyMat, interpolVec, dummyVec, dummyVec);

      // Sum values on nodes into vectors and record the touch count (build average of values)
      // This average is to get a continous inteface out of the discontinous
      // intefaces due to the DG method.
      // Cycling through all the nodes
      for (int i = 0; i < ele->NumNode(); ++i)
      {
        // Get the i-th node of the element
        DRT::Node *node = ele->Nodes()[i];
        // Get the local ID starting from the node's global id
        const int localIndex = dis.NodeRowMap()->LID(node->Id());
        ////If the local index is less than zero skip the for loop
        if (localIndex < 0) continue;
        ////Every entry of the vector gets its own touchCount entry so that we
        ////consider in different ways the position of the various nodes wrt others
        touchCount[localIndex]++;
        for (int d = 0; d < ndim; ++d)
        {
          (*electric)[d][localIndex] += interpolVec(i + d * ele->NumNode());
          (*magnetic)[d][localIndex] +=
              interpolVec(ele->NumNode() * (ndim) + i + d * ele->NumNode());
          (*trace)[d][localIndex] +=
              interpolVec(ele->NumNode() * (2 * ndim) + i + d * ele->NumNode());
        }
      }
    }
    for (int i = 0; i < electric->MyLength(); ++i)
    {
      for (int d = 0; d < ndim; ++d)
      {
        (*electric)[d][i] /= touchCount[i];
        (*magnetic)[d][i] /= touchCount[i];
      }
      for (int d = 0; d < ndim; ++d)
      {
        (*trace)[d][i] /= touchCount[i];
      }
    }
    dis.ClearState(true);
  }

  /*----------------------------------------------------------------------*
  |  Reads material properties from element for output   berardocco 03/18 |
  *----------------------------------------------------------------------*/
  void getElementMaterialProperties(DRT::Discretization &dis,
      Teuchos::RCP<Epetra_Vector> &conductivity, Teuchos::RCP<Epetra_Vector> &permittivity,
      Teuchos::RCP<Epetra_Vector> &permeability)
  {
    // For every element of the processor
    for (int el = 0; el < dis.NumMyRowElements(); ++el)
    {
      // Opening the element
      DRT::Element *ele = dis.lRowElement(el);

      const MAT::ElectromagneticMat *elemagmat =
          static_cast<const MAT::ElectromagneticMat *>(ele->Material().get());
      (*conductivity)[dis.ElementRowMap()->LID(ele->Id())] = elemagmat->sigma(ele->Id());
      (*permittivity)[dis.ElementRowMap()->LID(ele->Id())] = elemagmat->epsilon(ele->Id());
      (*permeability)[dis.ElementRowMap()->LID(ele->Id())] = elemagmat->mu(ele->Id());
    }

    return;
  }
}  // namespace

/*----------------------------------------------------------------------*
 |  Output (public)                                    berardocco 03/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::Output()
{
  TEUCHOS_FUNC_TIME_MONITOR("ELEMAG::ElemagTimeInt::Output");
  // Preparing the vectors that are going to be written in the output file
  electric.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  magnetic.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  trace.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));

  // Get the results from the discretization vectors to the output ones
  getNodeVectorsHDG(*discret_, trace_, numdim_, electric, magnetic, trace, conductivity,
      permittivity, permeability);

  // Create the new step
  output_->NewStep(step_, time_);

  if (step_ == 0)
  {
    getElementMaterialProperties(*discret_, conductivity, permittivity, permeability);
    output_->WriteVector("conductivity", conductivity);
    output_->WriteVector("permittivity", permittivity);
    output_->WriteVector("permeability", permeability);

    output_->WriteElementData(true);

    if (myrank_ == 0) std::cout << "======= Element properties written" << std::endl;
  }

  // Output the reuslts
  output_->WriteVector("magnetic", magnetic, IO::nodevector);
  output_->WriteVector("trace", trace, IO::nodevector);
  output_->WriteVector("electric", electric, IO::nodevector);

  // add restart data

  if (uprestart_ != 0 && step_ % uprestart_ == 0)
  {
    WriteRestart();
  }

  return;
}  // Output


/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                     berardocco 11/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::WriteRestart()
{
  if (myrank_ == 0) std::cout << "======= Restart written in step " << step_ << std::endl;

  output_->WriteVector("traceRestart", trace);

  // write internal field for which we need to create and fill the corresponding vectors
  // since this requires some effort, the WriteRestart method should not be used excessively!
  Teuchos::RCP<Epetra_Vector> intVar = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  Teuchos::RCP<Epetra_Vector> intVarnm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  discret_->SetState(1, "intVar", intVar);
  discret_->SetState(1, "intVarnm", intVarnm);

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", ELEMAG::fill_restart_vecs);
  eleparams.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);

  discret_->Evaluate(eleparams);

  Teuchos::RCP<const Epetra_Vector> matrix_state = discret_->GetState(1, "intVar");
  LINALG::Export(*matrix_state, *intVar);

  matrix_state = discret_->GetState(1, "intVarnm");
  LINALG::Export(*matrix_state, *intVarnm);

  output_->WriteVector("intVar", intVar);
  output_->WriteVector("intVarnm", intVarnm);

  discret_->ClearState(true);

  return;
}  // WriteRestart


/*----------------------------------------------------------------------*
 |  ReadRestart (public)                               berardocco 11/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_, step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");
  Teuchos::RCP<Epetra_Vector> intVar = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  try
  {
    reader.ReadVector(intVar, "intVar");
  }
  catch (...)
  {
    dserror(
        "Impossible to find restart data. Check if the restart step is an existing restart point.");
  }
  discret_->SetState(1, "intVar", intVar);

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", ELEMAG::ele_init_from_restart);
  eleparams.set<INPAR::ELEMAG::DynamicType>("dynamic type", elemagdyna_);

  Teuchos::RCP<Epetra_Vector> intVarnm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  try
  {
    reader.ReadVector(intVarnm, "intVarnm");
  }
  catch (...)
  {
    // if (myrank_ == 0)
    //{
    //  std::cout << "=========== Only one time step was found. Switch to BDF1." << std::endl;
    //}
    // eleparams.set<INPAR::ELEMAG::DynamicType>(
    //    "dynamic type", INPAR::ELEMAG::DynamicType::elemag_bdf1);
    dserror(
        "Impossible to find the additional vector of unknown necessary for the BDF2 integration. "
        "Consider fixing the code or restart a simulation that used BDF2 since the beginning.");
  }
  discret_->SetState(1, "intVarnm", intVarnm);
  reader.ReadMultiVector(trace, "traceRestart");

  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->ClearState(true);

  if (myrank_ == 0)
  {
    std::cout << "======= Restart of a previous simulation" << std::endl;
    std::cout << "Restart time: " << time_ << std::endl;
  }

  return;
}  // ReadRestart

void ELEMAG::ElemagTimeInt::SpySysmat(std::ostream &out)
{
  Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_, true)->EpetraMatrix()->Print(out);
  std::cout << "Routine has to be implemented. In the meanwhile the Print() method from the "
               "Epetra_CsrMatrix is used."
            << std::endl;
  /*
  //Dynamic casting of the sysmat
  Epetra_CrsMatrix *matrix =
  Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_,true)->EpetraMatrix().get(); int r =
  matrix->NumMyRows(); int c = matrix->NumMyCols(); int numentries;
  //double*& values = NULL;
  double* values;
  int* indices;
  //int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values, int*& Indices) const;
  for (unsigned int i = 0; i < r; ++i)
  {
    matrix->ExtractMyRowView(i, numentries, values, indices);
    for (unsigned int j = 0; j < c; ++j)
    {
      for (unsigned int q = 0; q < numentries; ++q)
        if (indices[q] == j && values[q] != 0.0)
          printf("x");
        else
          printf("o");
    }
    std::cout <<std::endl;
  }
  */
}

/*----------------------------------------------------------------------*
 |  Return discretization (public)                     berardocco 08/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ELEMAG::ElemagTimeInt::Discretization()
{
  return Teuchos::rcp_dynamic_cast<DRT::Discretization>(discret_);
}  // Discretization

/*----------------------------------------------------------------------*
 |  Create test field (public)                         berardocco 08/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ELEMAG::ElemagTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new ElemagResultTest(*this));
}  // CreateFieldTest
