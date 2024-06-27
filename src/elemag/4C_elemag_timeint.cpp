/*----------------------------------------------------------------------*/
/*! \file

\brief Base class functions for time integration of electromagnetics

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_elemag_timeint.hpp"

#include "4C_elemag_ele_action.hpp"
#include "4C_elemag_resulttest.hpp"
#include "4C_fem_discretization_hdg.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_electromagnetic.hpp"
#include "4C_utils_result_test.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                               gravemeier 06/17 |
 *----------------------------------------------------------------------*/
EleMag::ElemagTimeInt::ElemagTimeInt(const Teuchos::RCP<Core::FE::DiscretizationHDG> &actdis,
    const Teuchos::RCP<Core::LinAlg::Solver> &solver,
    const Teuchos::RCP<Teuchos::ParameterList> &params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter> &output)
    : discret_(actdis),
      solver_(solver),
      params_(params),
      output_(output),
      elemagdyna_(Core::UTILS::IntegralValue<Inpar::EleMag::DynamicType>(*params_, "TIMEINT")),
      myrank_(actdis->Comm().MyPID()),
      time_(0.0),
      step_(0),
      restart_(params_->get<int>("restart")),
      maxtime_(params_->get<double>("MAXTIME")),
      stepmax_(params_->get<int>("NUMSTEP")),
      uprestart_(params_->get<int>("RESTARTEVRY", -1)),
      upres_(params_->get<int>("RESULTSEVRY", -1)),
      numdim_(Global::Problem::Instance()->NDim()),
      dtp_(params_->get<double>("TIMESTEP")),
      tau_(params_->get<double>("TAU")),
      dtele_(0.0),
      dtsolve_(0.0),
      calcerr_(Core::UTILS::IntegralValue<bool>(*params_, "CALCERR")),
      postprocess_(Core::UTILS::IntegralValue<bool>(*params_, "POSTPROCESS")),
      errfunct_(params_->get<int>("ERRORFUNCNO", -1)),
      sourcefuncno_(params_->get<int>("SOURCEFUNCNO", -1)),
      equilibration_method_(
          Teuchos::getIntegralValue<Core::LinAlg::EquilibrationMethod>(*params, "EQUILIBRATION"))
{
  // constructor supposed to be empty!

}  // ElemagTimeInt



/*----------------------------------------------------------------------*
 |  initialization routine (public)                    berardocco 02/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::init()
{
  // get dof row map
  Teuchos::RCP<const Epetra_Map> dofrowmap = Teuchos::rcp(discret_->dof_row_map(), false);

  // check time-step length
  if (dtp_ <= 0.0) FOUR_C_THROW("Zero or negative time-step length!");

  // Nodevectors for the output
  electric = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  electric_post = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  magnetic = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  trace = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  conductivity = Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap()));
  permittivity = Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap()));
  permeability = Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap()));

  // create vector of zeros to be used for enforcing zero Dirichlet boundary conditions
  zeros_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  trace_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // Map of the dirichlet conditions
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  // Why is this in a new scope?
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    eleparams.set<const Core::UTILS::FunctionManager *>(
        "function_manager", &Global::Problem::Instance()->FunctionManager());
    const Core::ProblemType problem_type = Core::ProblemType::elemag;
    eleparams.set<const Core::ProblemType *>("problem_type", &problem_type);

    // Evaluation of the dirichlet conditions (why is it here and also later?)
    discret_->evaluate_dirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);

    // Initialize elements
    ElementsInit();
  }

  // create system matrix and set to zero
  // the 108 comes from line 282 of /fluid/fluidimplicitintegration.cpp
  sysmat_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap, 108, false, true));
  // Is it possible to avoid this passage? It is a sparse matrix so it should
  // only contain non-zero entries that have to be initialized
  sysmat_->Zero();

  // create residual vector
  residual_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // instantiate equilibration class
  equilibration_ =
      Teuchos::rcp(new Core::LinAlg::EquilibrationSparse(equilibration_method_, dofrowmap));

  // write mesh
  output_->write_mesh(0, 0.0);

  return;
}  // Init


/*----------------------------------------------------------------------*
 |  Print information to screen (public)               berardocco 02/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::print_information_to_screen()
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
    std::cout << "number of trace unknowns    " << discret_->dof_row_map(0)->NumGlobalElements()
              << std::endl;
    std::cout << "number of interior unknowns " << discret_->dof_row_map(1)->NumGlobalElements()
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
}  //  print_information_to_screen


/*----------------------------------------------------------------------*
 |  Time integration (public)                          berardocco 03/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::Integrate()
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
  TEUCHOS_FUNC_TIME_MONITOR("EleMag::ElemagTimeInt::Integrate");

  InitializeAlgorithm();

  // call elements to calculate system matrix/rhs and assemble
  assemble_mat_and_rhs();

  // time loop
  while (step_ < stepmax_ and time_ < maxtime_)
  {
    // increment time and step
    increment_time_and_step();

    // Apply BCS to RHS
    compute_silver_mueller(true);
    apply_dirichlet_to_system(true);

    // Solve for the trace values
    Solve();

    // Update the local solutions
    update_interior_variables_and_assemble_rhs();

    // The output to file only once in a while
    if (step_ % upres_ == 0)
    {
      output();
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

void EleMag::ElemagTimeInt::ElementsInit()
{
  // Initializing siome vectors and parameters
  Core::LinAlg::SerialDenseVector elevec1, elevec2, elevec3;
  Core::LinAlg::SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;

  // loop over all elements on the processor
  Core::Elements::Element::LocationArray la(2);
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // Selecting the elements
    Core::Elements::Element *ele = discret_->lColElement(el);

    // This function is a void function and therefore the input goes to the vector "la"
    ele->LocationVector(*discret_, la, true);

    // This is needed to store the dirichlet dofs in the
    if (std::find(la[0].lmdirich_.begin(), la[0].lmdirich_.end(), 1) != la[0].lmdirich_.end())
      initParams.set<std::vector<int> *>("dirichdof", &la[0].lmdirich_);
    else
      initParams.remove("dirichdof", false);

    initParams.set<int>("action", EleMag::ele_init);
    initParams.set<Inpar::EleMag::DynamicType>("dyna", elemagdyna_);
    Core::LinAlg::SerialDenseVector elevec1, elevec2, elevec3;
    Core::LinAlg::SerialDenseMatrix elemat1, elemat2;
    ele->evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Set initial field by given function (public)       berardocco 03/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::SetInitialField(const Inpar::EleMag::InitialField init, int startfuncno)
{
// time measurement: SetInitialField just in the debug phase
#ifdef FOUR_C_ENABLE_ASSERTIONS
  TEUCHOS_FUNC_TIME_MONITOR("EleMag::ElemagTimeInt::SetInitialField");
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
    case Inpar::EleMag::initfield_zero_field:
    {
      // Fancy printing to help debugging
      if (!myrank_)
      {
        std::cout << "Initializing a zero field." << std::endl;
      }

      break;
    }
    case Inpar::EleMag::initfield_field_by_function:
    {
      if (!myrank_)
      {
        std::cout << "Initializing field as specified by STARTFUNCNO " << startfuncno << std::endl;
      }

      // Initializing siome vectors and parameters
      Core::LinAlg::SerialDenseVector elevec1, elevec2, elevec3;
      Core::LinAlg::SerialDenseMatrix elemat1, elemat2;
      Teuchos::ParameterList initParams;
      initParams.set<int>("action", EleMag::project_field);
      initParams.set("startfuncno", startfuncno);
      initParams.set<double>("time", time_);
      initParams.set<double>("dt", dtp_);
      initParams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);
      // loop over all elements on the processor
      Core::Elements::Element::LocationArray la(2);
      for (int el = 0; el < discret_->NumMyColElements(); ++el)
      {
        // Selecting the elements
        Core::Elements::Element *ele = discret_->lColElement(el);

        // This function is a void function and therefore the input goes to the vector "la"
        ele->LocationVector(*discret_, la, false);

        // Reshaping the vectors
        if (static_cast<std::size_t>(elevec1.numRows()) != la[0].lm_.size())
          elevec1.size(la[0].lm_.size());
        if (elevec2.numRows() != discret_->NumDof(1, ele)) elevec2.size(discret_->NumDof(1, ele));
        ele->evaluate(
            initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
      }
      break;
    }  // case Inpar::EleMag::initfield_field_by_function
    default:
      FOUR_C_THROW("Option for initial field not implemented: %d", init);
      break;
  }  // switch(init)

  // Output of the initial condition
  output();
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
 |  Set initial field by scatra solution (public)      berardocco 05/20 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::set_initial_electric_field(
    Teuchos::RCP<Epetra_Vector> phi, Teuchos::RCP<Core::FE::Discretization> &scatradis)
{
  // we have to call an init for the elements first!
  Teuchos::ParameterList initParams;
  Core::Elements::Element::LocationArray la(2);

  Core::LinAlg::SerialDenseVector elevec1, elevec2;  //, elevec3;
  Core::LinAlg::SerialDenseMatrix elemat;            //, elemat2;

  initParams.set<int>("action", EleMag::project_electric_from_scatra_field);
  initParams.set<Inpar::EleMag::DynamicType>("dyna", elemagdyna_);

  Teuchos::RCP<Epetra_Vector> phicol;
  bool ishdg = false;
  if (Teuchos::rcp_dynamic_cast<Core::FE::DiscretizationHDG>(scatradis) != Teuchos::null)
  {
    phicol = Teuchos::rcp(new Epetra_Vector(*(scatradis->DofColMap(2))));
    ishdg = true;
    initParams.set<bool>("ishdg", ishdg);
  }
  else
    phicol = Teuchos::rcp(new Epetra_Vector(*(scatradis->DofColMap())));

  Core::LinAlg::Export(*phi, *phicol);

  // Loop MyColElements
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // determine owner of the scatra element
    Core::Elements::Element *scatraele = scatradis->lColElement(el);
    Core::Elements::Element *elemagele = discret_->lColElement(el);

    elemagele->LocationVector(*discret_, la, false);
    if (static_cast<std::size_t>(elevec1.numRows()) != la[0].lm_.size())
      elevec1.size(la[0].lm_.size());
    if (elevec2.numRows() != discret_->NumDof(1, elemagele))
      elevec2.size(discret_->NumDof(1, elemagele));

    Teuchos::RCP<Core::LinAlg::SerialDenseVector> nodevals_phi =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector);

    if (ishdg)
    {
      std::vector<int> localDofs = scatradis->Dof(2, scatraele);
      Core::FE::ExtractMyValues(*phicol, (*nodevals_phi), localDofs);
      // Obtain scatra ndofs knowing that the vector contains the values of the transported scalar
      // plus numdim_ components of its gradient
      initParams.set<unsigned int>("ndofs", localDofs.size() / (numdim_ + 1));
    }
    else
    {
      int numscatranode = scatraele->num_node();
      (*nodevals_phi).resize(numscatranode);
      // fill nodevals with node coords and nodebased solution values
      Core::Nodes::Node **scatranodes = scatraele->Nodes();
      for (int i = 0; i < numscatranode; ++i)
      {
        int dof = scatradis->Dof(0, scatranodes[i], 0);
        int lid = phicol->Map().LID(dof);
        if (lid < 0)
          FOUR_C_THROW("given dof is not stored on proc %d although map is colmap", myrank_);
        else
          (*nodevals_phi)[i] = (*(phicol.get()))[lid];
      }
    }

    initParams.set<Teuchos::RCP<Core::LinAlg::SerialDenseVector>>("nodevals_phi", nodevals_phi);

    // evaluate the element
    elemagele->evaluate(
        initParams, *discret_, la[0].lm_, elemat, elemat, elevec1, elevec2, elevec1);
  }
  discret_->Comm().Barrier();  // other procs please wait for the one, who did all the work

  output();

  return;
}  // set_initial_electric_field

/*----------------------------------------------------------------------*
 |  Compute error routine (public)                     berardocco 08/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SerialDenseVector> EleMag::ElemagTimeInt::compute_error()
{
  // Initializing siome vectors and parameters
  Core::LinAlg::SerialDenseVector elevec1, elevec2, elevec3;
  Core::LinAlg::SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList params;
  params.set<int>("action", EleMag::compute_error);
  params.set<int>("funcno", errfunct_);
  params.set<double>("time", time_);
  params.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);
  params.set<bool>("postprocess", postprocess_);

  const int numberOfErrorMeasures = 11;
  Teuchos::RCP<Core::LinAlg::SerialDenseVector> errors =
      Teuchos::rcp(new Core::LinAlg::SerialDenseVector(numberOfErrorMeasures));

  // call loop over elements (assemble nothing)
  discret_->EvaluateScalars(params, errors);
  discret_->ClearState(true);

  return errors;
}

void EleMag::ElemagTimeInt::PrintErrors(Teuchos::RCP<Core::LinAlg::SerialDenseVector> &errors)
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
    if ((*errors)[1] < 1e-14 && (*errors)[3] < 1e-14)
      std::cout << "Impossible to compute relative errors. The L2-norm of the analytical "
                   "solution is zero, resulting in a division by zero."
                << std::endl;
    else
    {
      if ((*errors)[1] > 1e-14)
        std::cout << "Electric relative L2-error: " << std::sqrt((*errors)[0] / (*errors)[1])
                  << std::endl;
      else
        std::cout
            << "Impossible to compute the electric relative error. The L2-norm of the analytical "
               "solution is zero, resulting in a division by zero."
            << std::endl;
      if ((*errors)[3] > 1e-14)
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

    std::cout << "\nPOST-PROCESSED ELECTRIC FIELD:" << std::endl;
    std::cout << "Absolute L2-error: " << std::sqrt((*errors)[8]) << std::endl;
    if ((*errors)[1] > 1e-14)
      std::cout << "Relative L2-error: " << std::sqrt((*errors)[8] / (*errors)[1]) << std::endl;
    else
      std::cout << "Impossible to compute the post-processed electric relative error. The L2-norm "
                   "of the analytical "
                   "solution is zero, resulting in a division by zero."
                << std::endl;
    std::cout << "Hdiv-error: " << std::sqrt((*errors)[9]) << std::endl;
    std::cout << "Hcurl-error: " << std::sqrt((*errors)[10]) << std::endl;
    std::cout
        << "-----------------------------------------------------------------------------------"
        << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Project interior variables for testing purposes     berardocco 07/18|
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::ProjectFieldTest(const int startfuncno)
{
  // Initializing siome vectors and parameters
  Core::LinAlg::SerialDenseVector elevec1, elevec2, elevec3;
  Core::LinAlg::SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;
  initParams.set<int>("action", EleMag::project_field_test);
  initParams.set("startfuncno", startfuncno);
  initParams.set<double>("time", time_);
  initParams.set<bool>("padaptivity", false);
  initParams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);

  // loop over all elements on the processor
  Core::Elements::Element::LocationArray la(2);
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // Selecting the elements
    Core::Elements::Element *ele = discret_->lColElement(el);

    // This function is a void function and therefore the input goes to the vector "la"
    ele->LocationVector(*discret_, la, false);

    // Reshaping the vectors
    if (static_cast<std::size_t>(elevec1.numRows()) != la[0].lm_.size())
      elevec1.size(la[0].lm_.size());
    if (elevec2.numRows() != discret_->NumDof(1, ele)) elevec2.size(discret_->NumDof(1, ele));
    ele->evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Project trace variable for testing purposes        berardocco 07/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::project_field_test_trace(const int startfuncno)
{
  // This map contains the trace values
  const Epetra_Map *dofrowmap = discret_->dof_row_map();

  // Initializing siome vectors and parameters
  Core::LinAlg::SerialDenseVector elevec1, elevec2, elevec3;
  Core::LinAlg::SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;
  initParams.set<int>("action", EleMag::project_field_test_trace);
  initParams.set("startfuncno", startfuncno);
  initParams.set<double>("time", time_);
  initParams.set<bool>("padaptivity", false);
  initParams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);
  // loop over all elements on the processor
  Core::Elements::Element::LocationArray la(2);
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    // Selecting the elements
    Core::Elements::Element *ele = discret_->lColElement(el);

    // This function is a void function and therefore the input goes to the vector "la"
    ele->LocationVector(*discret_, la, false);

    // Reshaping the vectors
    if (static_cast<std::size_t>(elevec1.numRows()) != la[0].lm_.size())
      elevec1.size(la[0].lm_.size());
    if (elevec2.numRows() != discret_->NumDof(1, ele)) elevec2.size(discret_->NumDof(1, ele));
    ele->evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);
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
void EleMag::ElemagTimeInt::InitializeAlgorithm()
{
  // In case We don't have BDF1 we initialize the method with BDF1 and then go back
  if (elemagdyna_ == Inpar::EleMag::DynamicType::elemag_bdf2 && restart_ == 0)
  {
    Inpar::EleMag::DynamicType temp_dyna = elemagdyna_;
    // First step with a BDF1
    elemagdyna_ = Inpar::EleMag::DynamicType::elemag_bdf1;

    // call elements to calculate system matrix/rhs and assemble
    assemble_mat_and_rhs();

    // increment time and step
    increment_time_and_step();

    // Apply BCS to RHS
    compute_silver_mueller(true);
    apply_dirichlet_to_system(true);

    // Solve
    Solve();

    update_interior_variables_and_assemble_rhs();

    // The output to file only once in a while
    if (step_ % upres_ == 0)
    {
      output();
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
void EleMag::ElemagTimeInt::assemble_mat_and_rhs()
{
  TEUCHOS_FUNC_TIME_MONITOR("EleMag::ElemagTimeInt::assemble_mat_and_rhs");

  // Fancy printing to help debugging
  if (!myrank_)
  {
    std::cout << "Creating system matrix elementwise and assembling in the final system matrix."
              << std::endl;
  }

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual and sysmat
  residual_->PutScalar(0.0);
  sysmat_->Zero();

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------

  // set general vector values needed by elements
  discret_->ClearState(true);

  discret_->set_state("trace", trace_);

  // set time step size
  eleparams.set<double>("dt", dtp_);
  eleparams.set<double>("tau", tau_);

  bool resonly = false;

  // set information needed by the elements
  eleparams.set<int>("sourcefuncno", sourcefuncno_);
  eleparams.set<bool>("resonly", resonly);
  eleparams.set<int>("action", EleMag::calc_systemmat_and_residual);
  eleparams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);
  eleparams.set<double>("time", time_);
  eleparams.set<double>("timep", time_ + dtp_);
  eleparams.set<int>("step", step_);

  // The evaluation of the discretization have to happen before Complete() is
  // called or after an UnComplete() call has been made.
  discret_->evaluate(eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null);
  discret_->ClearState(true);
  sysmat_->Complete();

  // Compute matrices and RHS for boundary conditions
  compute_silver_mueller(false);
  apply_dirichlet_to_system(false);

  // Equilibrate matrix
  equilibration_->equilibrate_matrix(sysmat_);

  return;
}  // assemble_mat_and_rhs

/*----------------------------------------------------------------------*
 | Updates interior variables and compute RHS (public)  berardocco 06/18|
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::update_interior_variables_and_assemble_rhs()
{
  TEUCHOS_FUNC_TIME_MONITOR("EleMag::ElemagTimeInt::update_interior_variables_and_assemble_rhs");

  // create parameterlist
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("sourcefuncno", sourcefuncno_);
  eleparams.set<double>("dt", dtp_);
  eleparams.set<double>("tau", tau_);
  eleparams.set<double>("time", time_);
  eleparams.set<double>("timep", time_ + dtp_);
  eleparams.set<int>("action", EleMag::update_secondary_solution_and_calc_residual);
  eleparams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);
  eleparams.set<int>("step", step_);
  eleparams.set<bool>("resonly", true);

  residual_->PutScalar(0.0);
  discret_->set_state("trace", trace_);

  discret_->evaluate(
      eleparams, Teuchos::null, Teuchos::null, residual_, Teuchos::null, Teuchos::null);

  discret_->ClearState(true);


  return;
}  // update_interior_variables_and_assemble_rhs

/*----------------------------------------------------------------------*
 |  Apply Dirichlet b.c. to system (public)            gravemeier 06/17 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::apply_dirichlet_to_system(bool resonly)
{
  TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
  Teuchos::ParameterList params;
  params.set<double>("total time", time_);
  params.set<const Core::UTILS::FunctionManager *>(
      "function_manager", &Global::Problem::Instance()->FunctionManager());
  const Core::ProblemType problem_type = Core::ProblemType::elemag;
  params.set<const Core::ProblemType *>("problem_type", &problem_type);
  discret_->evaluate_dirichlet(
      params, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  if (resonly)
    Core::LinAlg::apply_dirichlet_to_system(*residual_, *zeros_, *(dbcmaps_->cond_map()));
  else
    Core::LinAlg::apply_dirichlet_to_system(
        *sysmat_, *trace_, *residual_, *zeros_, *(dbcmaps_->cond_map()));

  return;
}  // apply_dirichlet_to_system

/*----------------------------------------------------------------------*
 |  Compute Silver-Mueller         (public)            berardocco 10/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::compute_silver_mueller(bool do_rhs)
{
  TEUCHOS_FUNC_TIME_MONITOR("      + Compute Silver-Mueller BC");

  // absorbing boundary conditions
  std::string condname = "Silver-Mueller";
  std::vector<Core::Conditions::Condition *> absorbingBC;
  discret_->GetCondition(condname, absorbingBC);

  // Check if there are Silver-Mueller BC
  if (absorbingBC.size())
  {
    Teuchos::ParameterList eleparams;
    eleparams.set<double>("time", time_);
    eleparams.set<bool>("do_rhs", do_rhs);
    eleparams.set<int>("action", EleMag::calc_abc);
    // Evaluate the boundary condition
    if (do_rhs)
      discret_->evaluate_condition(eleparams, residual_, condname);
    else
      discret_->evaluate_condition(
          eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null, condname);
  }

  return;
}  // apply_dirichlet_to_system

/*----------------------------------------------------------------------*
 |  Solve system for trace (public)                    berardocco 06/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::Solve()
{
  // Equilibrate RHS
  equilibration_->equilibrate_rhs(residual_);

  // This part has only been copied from the fluid part to be able to use algebraic multigrid
  // solvers and has to be checked

  discret_->compute_null_space_if_necessary(solver_->Params(), true);
  // solve for trace
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_->Solve(sysmat_->EpetraOperator(), trace_, residual_, solver_params);


  // Unequilibrate solution vector
  equilibration_->unequilibrate_increment(trace_);

  return;
}  // Solve


/*----------------------------------------------------------------------*
 |  Output to screen (public)                          berardocco 07/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::OutputToScreen()
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
  void getNodeVectorsHDG(Core::FE::Discretization &dis,
      const Teuchos::RCP<Epetra_Vector> &traceValues, const int ndim,
      Teuchos::RCP<Epetra_MultiVector> &electric, Teuchos::RCP<Epetra_MultiVector> &electric_post,
      Teuchos::RCP<Epetra_MultiVector> &magnetic, Teuchos::RCP<Epetra_MultiVector> &trace,
      Teuchos::RCP<Epetra_Vector> &conductivity, Teuchos::RCP<Epetra_Vector> &permittivity,
      Teuchos::RCP<Epetra_Vector> &permeability)
  {
    // create dofsets for electric and pressure at nodes
    // if there is no pressure vector it means that the vectors have not yet
    // been created and therefor it is needed to create them now.
    if (electric.get() == nullptr || electric->GlobalLength() != dis.NumGlobalNodes())
    {
      // The electric is a multivector because it is a vectorial field.
      // The multivector is based on the map of the node
      // owned by the processor. The vectors are zeroed.
      electric.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
      electric_post.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
      magnetic.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
    }

    // Same for the trace and cell pressure.
    trace.reset(new Epetra_MultiVector(*dis.NodeRowMap(), ndim));
    // call element routine for interpolate HDG to elements
    // Here it is used the function that acts in the elements, evaluate().
    Teuchos::ParameterList params;
    params.set<int>("action", EleMag::interpolate_hdg_to_node);
    // Setting a name to the dofs maps
    dis.set_state(0, "trace", traceValues);
    // Declaring all the necessary entry for evaluate()
    Core::Elements::Element::LocationArray la(2);
    Core::LinAlg::SerialDenseMatrix dummyMat;
    Core::LinAlg::SerialDenseVector dummyVec;
    Core::LinAlg::SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(dis.NumMyRowNodes());

    // For every element of the processor
    for (int el = 0; el < dis.NumMyColElements(); ++el)
    {
      // Opening the element
      Core::Elements::Element *ele = dis.lColElement(el);

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
      if (interpolVec.numRows() == 0) interpolVec.resize(ele->num_node() * 4 * ndim);
      for (int i = 0; i < interpolVec.length(); i++) interpolVec(i) = 0.0;

      // Interpolating hdg internal values to the node
      ele->evaluate(params, dis, la[0].lm_, dummyMat, dummyMat, interpolVec, dummyVec, dummyVec);

      // Sum values on nodes into vectors and record the touch count (build average of values)
      // This average is to get a continous inteface out of the discontinous
      // intefaces due to the DG method.
      // Cycling through all the nodes
      for (int i = 0; i < ele->num_node(); ++i)
      {
        // Get the i-th node of the element
        Core::Nodes::Node *node = ele->Nodes()[i];
        // Get the local ID starting from the node's global id
        const int localIndex = dis.NodeRowMap()->LID(node->Id());
        ////If the local index is less than zero skip the for loop
        if (localIndex < 0) continue;
        ////Every entry of the vector gets its own touchCount entry so that we
        ////consider in different ways the position of the various nodes wrt others
        touchCount[localIndex]++;
        for (int d = 0; d < ndim; ++d)
        {
          (*electric)[d][localIndex] += interpolVec(i + d * ele->num_node());
          (*electric_post)[d][localIndex] +=
              interpolVec(ele->num_node() * (2 * ndim) + i + d * ele->num_node());
          (*magnetic)[d][localIndex] +=
              interpolVec(ele->num_node() * (ndim) + i + d * ele->num_node());
          (*trace)[d][localIndex] +=
              interpolVec(ele->num_node() * (3 * ndim) + i + d * ele->num_node());
        }
      }
    }
    for (int i = 0; i < electric->MyLength(); ++i)
    {
      for (int d = 0; d < ndim; ++d)
      {
        (*electric)[d][i] /= touchCount[i];
        (*electric_post)[d][i] /= touchCount[i];
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
  void getElementMaterialProperties(Core::FE::Discretization &dis,
      Teuchos::RCP<Epetra_Vector> &conductivity, Teuchos::RCP<Epetra_Vector> &permittivity,
      Teuchos::RCP<Epetra_Vector> &permeability)
  {
    // For every element of the processor
    for (int el = 0; el < dis.NumMyRowElements(); ++el)
    {
      // Opening the element
      Core::Elements::Element *ele = dis.lRowElement(el);

      const Mat::ElectromagneticMat *elemagmat =
          static_cast<const Mat::ElectromagneticMat *>(ele->Material().get());
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
void EleMag::ElemagTimeInt::output()
{
  TEUCHOS_FUNC_TIME_MONITOR("EleMag::ElemagTimeInt::Output");
  // Preparing the vectors that are going to be written in the output file
  electric.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  electric_post.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  magnetic.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));
  trace.reset(new Epetra_MultiVector(*discret_->NodeRowMap(), numdim_));

  // Get the results from the discretization vectors to the output ones
  getNodeVectorsHDG(*discret_, trace_, numdim_, electric, electric_post, magnetic, trace,
      conductivity, permittivity, permeability);

  // Create the new step
  output_->new_step(step_, time_);

  if (step_ == 0)
  {
    getElementMaterialProperties(*discret_, conductivity, permittivity, permeability);
    output_->write_vector("conductivity", conductivity);
    output_->write_vector("permittivity", permittivity);
    output_->write_vector("permeability", permeability);

    output_->write_element_data(true);

    if (myrank_ == 0) std::cout << "======= Element properties written" << std::endl;
  }

  // Output the reuslts
  output_->write_vector("magnetic", magnetic, Core::IO::nodevector);
  output_->write_vector("trace", trace, Core::IO::nodevector);
  output_->write_vector("electric", electric, Core::IO::nodevector);
  output_->write_vector("electric_post", electric_post, Core::IO::nodevector);

  // add restart data

  if (uprestart_ != 0 && step_ % uprestart_ == 0)
  {
    write_restart();
  }

  return;
}  // Output


/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                     berardocco 11/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::write_restart()
{
  if (myrank_ == 0) std::cout << "======= Restart written in step " << step_ << std::endl;

  output_->write_vector("traceRestart", trace);

  // write internal field for which we need to create and fill the corresponding vectors
  // since this requires some effort, the write_restart method should not be used excessively!
  Teuchos::RCP<Epetra_Vector> intVar = Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map(1))));
  Teuchos::RCP<Epetra_Vector> intVarnm =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map(1))));
  discret_->set_state(1, "intVar", intVar);
  discret_->set_state(1, "intVarnm", intVarnm);

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", EleMag::fill_restart_vecs);
  eleparams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);

  discret_->evaluate(eleparams);

  Teuchos::RCP<const Epetra_Vector> matrix_state = discret_->GetState(1, "intVar");
  Core::LinAlg::Export(*matrix_state, *intVar);

  matrix_state = discret_->GetState(1, "intVarnm");
  Core::LinAlg::Export(*matrix_state, *intVarnm);

  output_->write_vector("intVar", intVar);
  output_->write_vector("intVarnm", intVarnm);

  discret_->ClearState(true);

  return;
}  // write_restart


/*----------------------------------------------------------------------*
 |  read_restart (public)                               berardocco 11/18 |
 *----------------------------------------------------------------------*/
void EleMag::ElemagTimeInt::read_restart(int step)
{
  Core::IO::DiscretizationReader reader(
      discret_, Global::Problem::Instance()->InputControlFile(), step);
  time_ = reader.read_double("time");
  step_ = reader.read_int("step");
  Teuchos::RCP<Epetra_Vector> intVar = Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map(1))));
  try
  {
    reader.read_vector(intVar, "intVar");
  }
  catch (...)
  {
    FOUR_C_THROW(
        "Impossible to find restart data. Check if the restart step is an existing restart point.");
  }
  discret_->set_state(1, "intVar", intVar);

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", EleMag::ele_init_from_restart);
  eleparams.set<Inpar::EleMag::DynamicType>("dynamic type", elemagdyna_);

  Teuchos::RCP<Epetra_Vector> intVarnm =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map(1))));
  try
  {
    reader.read_vector(intVarnm, "intVarnm");
  }
  catch (...)
  {
    // if (myrank_ == 0)
    //{
    //  std::cout << "=========== Only one time step was found. Switch to BDF1." << std::endl;
    //}
    // eleparams.set<Inpar::EleMag::DynamicType>(
    //    "dynamic type", Inpar::EleMag::DynamicType::elemag_bdf1);
    FOUR_C_THROW(
        "Impossible to find the additional vector of unknown necessary for the BDF2 integration. "
        "Consider fixing the code or restart a simulation that used BDF2 since the beginning.");
  }
  discret_->set_state(1, "intVarnm", intVarnm);
  reader.read_multi_vector(trace, "traceRestart");

  discret_->evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->ClearState(true);

  if (myrank_ == 0)
  {
    std::cout << "======= Restart of a previous simulation" << std::endl;
    std::cout << "Restart time: " << time_ << std::endl;
  }

  return;
}  // read_restart

void EleMag::ElemagTimeInt::SpySysmat(std::ostream &out)
{
  Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_, true)->EpetraMatrix()->Print(out);
  std::cout << "Routine has to be implemented. In the meanwhile the print() method from the "
               "Epetra_CsrMatrix is used."
            << std::endl;
  /*
  //Dynamic casting of the sysmat
  Epetra_CrsMatrix *matrix =
  Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_,true)->EpetraMatrix().get(); int r =
  matrix->NumMyRows(); int c = matrix->NumMyCols(); int numentries;
  //double*& values = nullptr;
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
Teuchos::RCP<Core::FE::Discretization> EleMag::ElemagTimeInt::discretization()
{
  return discret_;
}  // discretization

/*----------------------------------------------------------------------*
 |  Create test field (public)                         berardocco 08/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> EleMag::ElemagTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new ElemagResultTest(*this));
}  // CreateFieldTest

FOUR_C_NAMESPACE_CLOSE
