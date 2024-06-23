/*--------------------------------------------------------------------------*/
/*! \file

\brief basis of all monolithic EHL algorithms that perform a coupling between
       the structure field equation and lubrication field equations

\level 3
*/
/*--------------------------------------------------------------------------*/



/*----------------------------------------------------------------------*
 | headers                                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
#include "4C_ehl_monolithic.hpp"

#include "4C_adapter_coupling_ehl_mortar.hpp"
#include "4C_adapter_lubrication.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_node.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_condition_locsys.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_solver.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_lubrication_ele_action.hpp"
#include "4C_lubrication_timint_implicit.hpp"
#include "4C_mat_lubrication_mat.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_node.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.



/*----------------------------------------------------------------------*
 | monolithic                                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
EHL::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string lubrication_disname)
    : Base(comm, globaltimeparams, lubricationparams, structparams, struct_disname,
          lubrication_disname),
      solveradapttol_(
          Core::UTILS::IntegralValue<int>(
              ((Global::Problem::Instance()->elasto_hydro_dynamic_params()).sublist("MONOLITHIC")),
              "ADAPTCONV") == 1),
      solveradaptolbetter_(
          ((Global::Problem::Instance()->elasto_hydro_dynamic_params()).sublist("MONOLITHIC"))
              .get<double>("ADAPTCONV_BETTER")),
      printiter_(true),  // ADD INPUT PARAMETER
      zeros_(Teuchos::null),
      strmethodname_(
          Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(structparams, "DYNAMICTYP")),
      ehldyn_(Global::Problem::Instance()->elasto_hydro_dynamic_params()),
      ehldynmono_(
          (Global::Problem::Instance()->elasto_hydro_dynamic_params()).sublist("MONOLITHIC")),
      blockrowdofmap_(Teuchos::null),
      systemmatrix_(Teuchos::null),
      k_sl_(Teuchos::null),
      k_ls_(Teuchos::null),
      merge_ehl_blockmatrix_(
          Core::UTILS::IntegralValue<bool>(ehldynmono_, "MERGE_EHL_BLOCK_MATRIX")),
      soltech_(Core::UTILS::IntegralValue<Inpar::EHL::NlnSolTech>(ehldynmono_, "NLNSOL")),
      iternorm_(Core::UTILS::IntegralValue<Inpar::EHL::VectorNorm>(ehldynmono_, "ITERNORM")),
      iter_(0),
      sdyn_(structparams),
      timernewton_("EHL_Monolithic_newton", true)
{
  blockrowdofmap_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor);

  // --------------------------------- EHL solver: create a linear solver

  // get iterative solver
  if (merge_ehl_blockmatrix_ == false) create_linear_solver();
  // get direct solver, e.g. UMFPACK
  else  // (merge_ehl_blockmatrix_ == true)
  {
    // get solver parameter list of linear TSI solver
    const int linsolvernumber = ehldynmono_.get<int>("LINEAR_SOLVER");
    const Teuchos::ParameterList& ehlsolverparams =
        Global::Problem::Instance()->SolverParams(linsolvernumber);

    Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);
    *solverparams = ehlsolverparams;

    solver_ = Teuchos::rcp(new Core::LinAlg::Solver(*solverparams, Comm(),
        Global::Problem::Instance()->solver_params_callback(),
        Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
            Global::Problem::Instance()->IOParams(), "VERBOSITY")));
  }  // end BlockMatrixMerge

}  // Monolithic()


/*----------------------------------------------------------------------*
 | prepare time step (public)                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::prepare_time_step()
{
  // counter and print header
  // increment time and step counter
  increment_time_and_step();
  print_header();

  // apply current pressure to structure
  // set the external fluid force on the structure, which result from the fluid pressure
  set_lubrication_solution(lubrication_->LubricationField()->Prenp());

  // call the predictor
  structure_field()->prepare_time_step();
  lubrication_->LubricationField()->prepare_time_step();

}  // prepare_time_step()


/*----------------------------------------------------------------------*
 | create linear solver                                     wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::create_linear_solver()
{
  // get the solver number used for linear EHL solver
  const int linsolvernumber = ehldynmono_.get<int>("LINEAR_SOLVER");
  // check if the EHL solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for monolithic EHL. Please set LINEAR_SOLVER in ELASTO HYDRO "
        "DYNAMIC to a valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  // get parameter list of lubrication dynamics
  const Teuchos::ParameterList& ldyn = Global::Problem::Instance()->lubrication_dynamic_params();
  // use solver blocks for pressure (lubrication field)
  // get the solver number used for lubrication solver
  const int tlinsolvernumber = ldyn.get<int>("LINEAR_SOLVER");
  // check if the EHL solver has a valid solver number
  if (tlinsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for lubrication field. Please set LINEAR_SOLVER in THERMAL "
        "DYNAMIC to a valid number!");

  // get solver parameter list of linear EHL solver
  const Teuchos::ParameterList& ehlsolverparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);

  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(ehlsolverparams, "SOLVER");

  if (solvertype != Core::LinearSolver::SolverType::belos)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and lubrication solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the dat files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    FOUR_C_THROW("Iterative solver expected");
  }
  const auto azprectype =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(ehlsolverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2:
      break;
    case Core::LinearSolver::PreconditionerType::multigrid_muelu:
    case Core::LinearSolver::PreconditionerType::multigrid_nxn:
    case Core::LinearSolver::PreconditionerType::cheap_simple:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      FOUR_C_THROW(
          "Block Gauss-Seidel BGS2x2 preconditioner expected. Alternatively you can define your "
          "own AMG block preconditioner (using an xml file). This is experimental.");
      break;
  }


  // prepare linear solvers and preconditioners
  switch (azprectype)
  {
    case Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2:
    case Core::LinearSolver::PreconditionerType::multigrid_nxn:
    case Core::LinearSolver::PreconditionerType::cheap_simple:
    {
      // This should be the default case (well-tested and used)
      solver_ = Teuchos::rcp(new Core::LinAlg::Solver(ehlsolverparams,
          // ggfs. explizit Comm von STR wie lungscatra
          Comm(), Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY")));

      // use solver blocks for structure and pressure (lubrication field)
      const Teuchos::ParameterList& ssolverparams =
          Global::Problem::Instance()->SolverParams(slinsolvernumber);
      const Teuchos::ParameterList& tsolverparams =
          Global::Problem::Instance()->SolverParams(tlinsolvernumber);

      solver_->put_solver_params_to_sub_params("Inverse1", ssolverparams,
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));
      solver_->put_solver_params_to_sub_params("Inverse2", tsolverparams,
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));

      // prescribe rigid body modes
      structure_field()->discretization()->compute_null_space_if_necessary(
          solver_->Params().sublist("Inverse1"));
      lubrication_->LubricationField()->discretization()->compute_null_space_if_necessary(
          solver_->Params().sublist("Inverse2"));


      if (azprectype == Core::LinearSolver::PreconditionerType::cheap_simple)
      {
        // Tell to the Core::LinAlg::SOLVER::SimplePreconditioner that we use the general
        // implementation
        solver_->Params().set<bool>("GENERAL", true);
      }

      break;
    }
    case Core::LinearSolver::PreconditionerType::multigrid_muelu:
    {
      solver_ = Teuchos::rcp(new Core::LinAlg::Solver(ehlsolverparams,
          // ggfs. explizit Comm von STR wie lungscatra
          Comm(), Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY")));

      // use solver blocks for structure and pressure (lubrication field)
      const Teuchos::ParameterList& ssolverparams =
          Global::Problem::Instance()->SolverParams(slinsolvernumber);
      const Teuchos::ParameterList& tsolverparams =
          Global::Problem::Instance()->SolverParams(tlinsolvernumber);

      // This is not very elegant:
      // first read in solver parameters. These have to contain ML parameters such that...
      solver_->put_solver_params_to_sub_params("Inverse1", ssolverparams,
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));
      solver_->put_solver_params_to_sub_params("Inverse2", tsolverparams,
          Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY"));

      // ... 4C calculates the null space vectors. These are then stored in the sublists
      //     Inverse1 and Inverse2 from where they...
      structure_field()->discretization()->compute_null_space_if_necessary(
          solver_->Params().sublist("Inverse1"));
      lubrication_->LubricationField()->discretization()->compute_null_space_if_necessary(
          solver_->Params().sublist("Inverse2"));

      // ... are copied from here to ...
      const Teuchos::ParameterList& inv1source =
          solver_->Params().sublist("Inverse1").sublist("ML Parameters");
      const Teuchos::ParameterList& inv2source =
          solver_->Params().sublist("Inverse2").sublist("ML Parameters");

      // ... here. The "MueLu Parameters" sublists "Inverse1" and "Inverse2" only contain the basic
      //     information about the corresponding null space vectors, which are actually copied ...
      Teuchos::ParameterList& inv1 =
          solver_->Params().sublist("MueLu Parameters").sublist("Inverse1");
      Teuchos::ParameterList& inv2 =
          solver_->Params().sublist("MueLu Parameters").sublist("Inverse2");

      // ... here.
      inv1.set<int>("PDE equations", inv1source.get<int>("PDE equations"));
      inv2.set<int>("PDE equations", inv2source.get<int>("PDE equations"));
      inv1.set<int>("null space: dimension", inv1source.get<int>("null space: dimension"));
      inv2.set<int>("null space: dimension", inv2source.get<int>("null space: dimension"));
      inv1.set<double*>("null space: vectors", inv1source.get<double*>("null space: vectors"));
      inv2.set<double*>("null space: vectors", inv2source.get<double*>("null space: vectors"));
      inv1.set<Teuchos::RCP<Epetra_MultiVector>>(
          "nullspace", inv1source.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace"));
      inv2.set<Teuchos::RCP<Epetra_MultiVector>>(
          "nullspace", inv2source.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace"));

      solver_->Params().sublist("MueLu Parameters").set("EHL", true);
      break;
    }
    default:
      FOUR_C_THROW("Block Gauss-Seidel BGS2x2 preconditioner expected");
      break;
  }

}  // create_linear_solver()


/*----------------------------------------------------------------------*
 | non-linear solve, i.e. (multiple) corrector (public)     wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::Solve()
{
  // choose solution technique according to input file
  switch (soltech_)
  {
    // Newton-Raphson iteration
    case Inpar::EHL::soltech_newtonfull:
      NewtonFull();
      break;
    // catch problems
    default:
      FOUR_C_THROW("Solution technique \"%s\" is not implemented",
          Inpar::EHL::NlnSolTechString(soltech_).c_str());
      break;
  }  // end switch (soltechnique_)

  return;
}  // Solve()


/*----------------------------------------------------------------------*
 | time loop of the monolithic system                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    // predict solution of both field (call the adapter)
    prepare_time_step();

    // integrate time step
    Solve();

    // calculate stresses, strains, energies
    constexpr bool force_prepare = false;
    prepare_output(force_prepare);

    // update all single field solvers
    update();

    // write output to screen and files
    output();

  }  // not_finished
}  // TimeLoop()


/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration              wirtz 01/16 |
 | in ehl_algorithm: NewtonFull()                                       |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::NewtonFull()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic EHL tangent matrix

  // initialise equilibrium loop
  iter_ = 1;

  // incremental solution vector with length of all EHL dofs
  iterinc_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  zeros_->PutScalar(0.0);

  //------------------------------------------------------ iteration loop

  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    // reset timer
    timernewton_.reset();

    // compute residual forces #rhs_ and tangent #systemmatrix_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) evaluate_force_stiff_residual(),
    evaluate(iterinc_);

    // create the linear system
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    setup_system_matrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled()) FOUR_C_THROW("Effective tangent matrix must be filled here");

    // create full monolithic rhs vector
    // make negative residual not necessary: rhs_ is already negative
    // (STR/LUBRICATION)-RHS is put negative in prepare_system_for_newton_solve()
    setup_rhs();

    if (dry_contact_)
    {
      mortaradapter_->CondenseContact(systemmatrix_, rhs_, structure_field()->Dispnp(), Dt());
      apply_dbc();
    }
    // *********** time measurement ***********
    double dtcpu = timernewton_.wallTime();
    // *********** time measurement ***********
    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in prepare_system_for_newton_solve() within evaluate(iterinc_)
    linear_solve();
    // *********** time measurement ***********
    dtsolve_ = timernewton_.wallTime() - dtcpu;
    // *********** time measurement ***********

    // vector of displacement and pressure increments
    Teuchos::RCP<Epetra_Vector> sx;
    Teuchos::RCP<Epetra_Vector> lx;
    // extract field vectors
    extract_field_vectors(iterinc_, sx, lx);

    if (dry_contact_) mortaradapter_->RecoverCoupled(sx, lx);

    // reset solver tolerance
    solver_->ResetTolerance();

    // --------------------------------------------- build residual norms
    // include all stuff here related with convergence test
    normrhs_ = calculate_vector_norm(iternorm_, rhs_);
    // vector of displacement and pressure residual
    Teuchos::RCP<Epetra_Vector> strrhs;
    Teuchos::RCP<Epetra_Vector> lubricationrhs;
    // extract field vectors
    extract_field_vectors(rhs_, strrhs, lubricationrhs);
    normstrrhs_ = calculate_vector_norm(iternormstr_, strrhs);
    normlubricationrhs_ = calculate_vector_norm(iternormlubrication_, lubricationrhs);

    // --------------------------------- build residual incremental norms
    norminc_ = calculate_vector_norm(iternorm_, iterinc_);
    normdisi_ = calculate_vector_norm(iternormstr_, sx);
    normprei_ = calculate_vector_norm(iternormlubrication_, lx);

    // in case of 'Mix'-convergence criterion: save the norm of the 1st
    // iteration in (norm . iter0_)
    if (iter_ == 1)
    {
      // save initial residual norms
      normrhsiter0_ = normrhs_;
      normstrrhsiter0_ = normstrrhs_;
      normlubricationrhsiter0_ = normlubricationrhs_;
      // save initial incremental norms
      norminciter0_ = norminc_;
      normdisiiter0_ = normdisi_;
      normpreiiter0_ = normprei_;

      // set the minimum of iter0_ and tolrhs_, because we want to prevent the
      // case of a zero characteristic initial norm
      if (normrhsiter0_ == 0.0) normrhsiter0_ = tolrhs_;
      if (normstrrhsiter0_ == 0.0) normstrrhsiter0_ = tolstrrhs_;
      if (normlubricationrhsiter0_ == 0.0) normlubricationrhsiter0_ = tollubricationrhs_;
      if (norminciter0_ == 0.0) norminciter0_ = tolinc_;
      if (normdisiiter0_ == 0.0) normdisiiter0_ = toldisi_;
      if (normpreiiter0_ == 0.0) normpreiiter0_ = tolprei_;
    }

    // print stuff
    print_newton_iter();

    // increment equilibrium loop index
    iter_ += 1;

  }  // end equilibrium loop

  // ----------------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((Converged()) and (Comm().MyPID() == 0))
  {
    print_newton_conv();
  }
  else if (Core::UTILS::IntegralValue<Inpar::STR::DivContAct>(sdyn_, "DIVERCONT") ==
           Inpar::STR::divcont_continue)
    ;
  else if (iter_ >= itermax_)
    FOUR_C_THROW("Newton unconverged in %d iterations", iter_);

}  // NewtonFull()

/*----------------------------------------------------------------------*
 | evaluate the single fields                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::evaluate(Teuchos::RCP<Epetra_Vector> stepinc)
{
  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::Evaluate");

  // displacement and pressure incremental vector
  Teuchos::RCP<Epetra_Vector> sx;
  Teuchos::RCP<Epetra_Vector> lx;

  // if an increment vector exists
  if (stepinc != Teuchos::null)
  {
    // extract displacement sx and pressure lx incremental vector of global
    // unknown incremental vector x
    extract_field_vectors(stepinc, sx, lx);
  }

  // Newton update of the lubrication field
  // update pressure before passed to the structural field
  lubrication_->LubricationField()->UpdateNewton(
      lx);  // The same should be done for the structure field before updating the couplingh...

  // pass the structural values to the lubrication field,
  // i.e. set mesh displacement, velocity fields and film thickness
  // note: the iteration update has not been done yet
  Teuchos::RCP<Epetra_Vector> new_disp = Teuchos::rcp(new Epetra_Vector(*structure_->Dispnp()));
  if (new_disp->Update(1., *sx, 1.)) FOUR_C_THROW("update failed");

  // set interface height, velocity etc to lubrication field
  set_struct_solution(new_disp);

  // apply current pressure to structure
  // set the external fluid force on the structure, which result from the fluid pressure
  set_lubrication_solution(lubrication_->LubricationField()->Prenp());

  // structure Evaluate (builds tangent, residual and applies DBC)
  structure_field()->evaluate(sx);
  structure_field()->discretization()->ClearState(true);

  /// lubrication field

  // lubrication Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  lubrication_->LubricationField()->evaluate();
  lubrication_->LubricationField()->discretization()->ClearState(true);

}  // evaluate()



/*----------------------------------------------------------------------*
 | extract field vectors for calling evaluate() of the      wirtz 01/16 |
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::extract_field_vectors(
    Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector>& sx, Teuchos::RCP<Epetra_Vector>& lx)
{
  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::extract_field_vectors");

  // process structure unknowns of the first field
  sx = extractor()->extract_vector(x, 0);

  // process lubrication unknowns of the second field
  lx = extractor()->extract_vector(x, 1);
}  // extract_field_vectors()


/*----------------------------------------------------------------------*
 | full monolithic dof row map                              wirtz 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> EHL::Monolithic::dof_row_map() const
{
  return extractor()->FullMap();
}  // dof_row_map()


/*----------------------------------------------------------------------*
 | setup system (called in ehl_dyn)                         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetupSystem()
{
  // set parameters that remain the same in the whole calculation
  set_default_parameters();

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

  // use its own dof_row_map, that is the 0th map of the discretization
  vecSpaces.push_back(structure_field()->dof_row_map(0));
  vecSpaces.push_back(lubrication_->LubricationField()->dof_row_map(0));

  if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No pressure equation. Panic.");

  set_dof_row_maps(vecSpaces);

  /*----------------------------------------------------------------------*/
  // initialise EHL-systemmatrix_
  systemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *extractor(), *extractor(), 81, false, true));

  // create empty matrix
  k_sl_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(structure_field()->discretization()->dof_row_map(0)), 81, true, true));

  // create empty matrix
  k_ls_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(lubrication_->LubricationField()->discretization()->dof_row_map(0)), 81, true, true));

}  // SetupSystem()


/*----------------------------------------------------------------------*
 | put the single maps to one full EHL map together         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::merge_maps(maps);

  // full EHL-blockmap
  extractor()->setup(*fullmap, maps);
}  // set_dof_row_maps()


/*----------------------------------------------------------------------*
 | setup system matrix of EHL                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::setup_system_matrix()
{
  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::setup_system_matrix");



  //--------------------------------
  // 0. Prepare the linearizations...
  //--------------------------------

  // Time integration specific parameters..
  double alphaf = -1.;
  switch (Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(sdyn_, "DYNAMICTYP"))
  {
    case Inpar::STR::dyna_genalpha:
    {
      alphaf = sdyn_.sublist("GENALPHA").get<double>("ALPHA_F");
      break;
    }
    case Inpar::STR::dyna_statics:
    {
      alphaf = 0.;
      break;
    }
    default:
      FOUR_C_THROW("unknown time integration strategy for structural problem in coupled EHL");
  }

  //--------------------------------------------------------------------------------------
  // 1. Calculate and assemble k_ss (Linearization of structure residual wrt displacements)
  //--------------------------------------------------------------------------------------

  // k_ss basically consists of the standard effective dynamic stiffness matrix plus the
  // linearization of the lubrication interface force wrt the displacements


  //----------------------------------------------
  // 1.1 Get the effective dynamic stiffness matrix
  //----------------------------------------------

  // Effective dynamic stiffness matrix, Dirichlet already applied.
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ss = structure_field()->system_matrix();
  k_ss->UnComplete();

  //----------------------------------------------------------------------
  // 1.2 Linearization of lubrication interface force wrt the displacements
  //----------------------------------------------------------------------

  // Lubrication interface force: f_lub = B * t = D^T*t - M^T*t
  // Derivative of lubrication interface force: df_lub/dd = dB/dd * t + B * dt/dd

  // 1.2.1 Derivative due to deformation-dependent Mortar matrices
  // d(D^T)/dd * t
  Teuchos::RCP<Epetra_Vector> stritraction_D_col =
      Core::LinAlg::CreateVector(*(mortaradapter_->Interface()->SlaveColDofs()), true);
  Core::LinAlg::Export(*stritraction_D_, *stritraction_D_col);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> slaveiforce_derivd1 =
      mortaradapter_->AssembleEHLLinD(stritraction_D_col);

  // d(-M^T)/dd * t
  Teuchos::RCP<Epetra_Vector> stritraction_M_col =
      Core::LinAlg::CreateVector(*(mortaradapter_->Interface()->SlaveColDofs()), true);
  Core::LinAlg::Export(*stritraction_M_, *stritraction_M_col);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> masteriforce_derivd1 =
      mortaradapter_->AssembleEHLLinM(stritraction_M_col);
  masteriforce_derivd1->Scale(-1.0);

  // Add the negative midpoint values (remember, we are linearizing an extforce-like term)
  k_ss->Add(*slaveiforce_derivd1, false, -(1.0 - alphaf), 1.0);
  k_ss->Add(*masteriforce_derivd1, false, -(1.0 - alphaf), 1.0);

  // linearization of D.t and M.t (the part with D.(dt/dd))
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ds_dd =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->SlaveDofMap(), 81));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dm_dd =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->MasterDofMap(), 81));

  // pressure term
  lin_pressure_force_disp(ds_dd, dm_dd);
  // couette term
  LinCouetteForceDisp(ds_dd, dm_dd);
  // poiseuille term
  lin_poiseuille_force_disp(ds_dd, dm_dd);

  // complete matrices
  ds_dd->Complete(*mortaradapter_->SMdofMap(), *mortaradapter_->SlaveDofMap());
  dm_dd->Complete(*mortaradapter_->SMdofMap(), *mortaradapter_->MasterDofMap());

  // Add the negative midpoint values (remember, we are linearizing an extforce-like term)
  k_ss->Add(*ds_dd, false, -(1.0 - alphaf), 1.0);
  k_ss->Add(*dm_dd, false, -(1.0 - alphaf), 1.0);

  k_ss->Complete();
  k_ss->ApplyDirichlet(*structure_field()->get_dbc_map_extractor()->cond_map(), true);

  // Assign k_ss to system matrix
  systemmatrix_->Assign(0, 0, Core::LinAlg::View, *k_ss);


  //---------------------------------------------------------------------------------
  // 2. Calculate and assemble k_sl (Linearization of structure residual wrt pressure)
  //---------------------------------------------------------------------------------

  // Derivative of the lubrication interface force: df_lub/dp = B * dt/dp

  k_sl_->reset();
  k_sl_->UnComplete();

  // linearization of traction w.r.t. fluid pressure
  Teuchos::RCP<Core::LinAlg::SparseMatrix> slaveiforce_derivp =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->SlaveDofMap(), 81, true, true));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> masteriforce_derivp =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->MasterDofMap(), 81, true, true));

  // pressure force
  lin_pressure_force_pres(slaveiforce_derivp, masteriforce_derivp);

  // poiseuille force
  lin_poiseuille_force_pres(slaveiforce_derivp, masteriforce_derivp);

  // couette force
  LinCouetteForcePres(slaveiforce_derivp, masteriforce_derivp);

  slaveiforce_derivp->Complete(
      *lubrication_->LubricationField()->dof_row_map(), *mortaradapter_->SlaveDofMap());
  masteriforce_derivp->Complete(
      *lubrication_->LubricationField()->dof_row_map(), *mortaradapter_->MasterDofMap());

  // Add the negative midpoint values (remember, we are linearizing an extforce-like term)
  k_sl_->Add(*slaveiforce_derivp, false, -(1.0 - alphaf), 1.0);
  k_sl_->Add(*masteriforce_derivp, false, -(1.0 - alphaf), 1.0);

  k_sl_->Complete(*(extractor()->Map(1)),  // pressue dof map
      *(extractor()->Map(0))               // displacement dof map
  );

  // No DBC need to be applied, since lubrication interface disp-dofs must NOT have dbc conditions!
  k_sl_->ApplyDirichlet(*structure_field()->get_dbc_map_extractor()->cond_map(), false);

  // Assign k_sl to system matrix
  systemmatrix_->Assign(0, 1, Core::LinAlg::View, *(k_sl_));


  //-----------------------------------------------------------------------------------
  // 3. Calculate and assemble k_ll (Linearization of lubrication residual wrt pressure)
  //-----------------------------------------------------------------------------------

  // Pure hydrodynamic stiffness matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ll = lubrication_->LubricationField()->SystemMatrix();

  // Assign k_ll to system matrix
  systemmatrix_->Assign(1, 1, Core::LinAlg::View, *(k_ll));


  //----------------------------------------------------------------------------------------
  // 4. Calculate and assemble k_ls (Linearization of lubrication residual wrt displacements)
  //----------------------------------------------------------------------------------------

  k_ls_->reset();
  k_ls_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *(lubrication_->LubricationField()->discretization()->dof_row_map(0)), 81, false, true));

  /*
   * Three distinct types of derivatives of the lubrication residual wrt to interface displacement
   * exist: -Derivatives due to the deformation dependent film height -Derivatives due to the
   * deformation dependent tangential velocities -Derivatives due to the deformation of the
   * lubrication mesh (not implemented yet)
   */

  /*
   * Instead of assembling the element-wise contributions directly into k_ls, following strategy is
   * followed: The discrete nodal derivatives of the film height and the tangential velocities are
   * calculated and assembled globally. The integrals associated with those derivatives are
   * evaluated element-wise and assembled into global matrices. The product of the discret
   * derivatives with those matrices finally enter the stiffness matrix k_ls
   */

  //-------------------------------------------------------------------------------------------
  // 4.0 Calculate and assemble the matrices, which are associated with the discrete derivatives
  //-------------------------------------------------------------------------------------------

  // Global matrix associated with the discrete film height derivative
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ls_linH =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(extractor()->Map(1)), 81, false, false));

  // Global matrix associated with the discrete tangential velocity derivative
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ls_linV =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(extractor()->Map(1)), 81, false, false));

  // Call elements and assemble
  apply_lubrication_coupl_matrix(k_ls_linH, k_ls_linV);

  //-----------------------------------
  // 4.1 Discrete film height derivative
  //-----------------------------------
  Teuchos::RCP<Core::LinAlg::SparseMatrix> ddgap_dd = mortaradapter_->Nodal_GapDeriv();

  Teuchos::RCP<Core::LinAlg::SparseMatrix> dh_dd = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*ada_strDisp_to_lubDisp_->SlaveDofMap(), 81, true, true));

  Core::LinAlg::MatrixRowTransform()(*ddgap_dd, 1.0,
      Core::Adapter::CouplingMasterConverter(*ada_strDisp_to_lubDisp_), *dh_dd, false);
  dh_dd->Complete(*(extractor()->Map(0)), *ada_strDisp_to_lubDisp_->SlaveDofMap());

  // Multiply with associated matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ls_H =
      Core::LinAlg::MLMultiply(*k_ls_linH, false, *dh_dd, false, false, false, true);

  // Add the contribution of film height derivative to the stiffness matrix
  k_ls_->Add(*k_ls_H, false, 1.0, 1.0);

  //-------------------------------------------
  // 4.2 Discrete tangential velocity derivative
  //-----------------------------------------
  Teuchos::RCP<Core::LinAlg::SparseMatrix> avTangVelDeriv = mortaradapter_->AvTangVelDeriv();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dst = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *lubrication_->LubricationField()->discretization()->dof_row_map(1), 81, true, false));
  Core::LinAlg::MatrixRowTransform().operator()(*avTangVelDeriv, 1.,
      Core::Adapter::CouplingMasterConverter(*ada_strDisp_to_lubDisp_), *dst, false);
  dst->Complete(*structure_field()->dof_row_map(),
      *lubrication_->LubricationField()->discretization()->dof_row_map(1));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp = Core::LinAlg::MLMultiply(*k_ls_linV, *dst, true);
  k_ls_->Add(*tmp, false, -1.0, 1.0);


  k_ls_->Complete(*(extractor()->Map(0)),  // displacement dof map
      *(extractor()->Map(1))               // pressue dof map
  );

  // Apply Dirichet to k_ls
  k_ls_->ApplyDirichlet(
      *lubrication_->LubricationField()->get_dbc_map_extractor()->cond_map(), false);
  if (inf_gap_toggle_lub_ != Teuchos::null) k_ls_->ApplyDirichlet(*inf_gap_toggle_lub_, false);

  // Assign k_ls to system matrix
  systemmatrix_->Assign(1, 0, Core::LinAlg::View, *k_ls_);

  // Finished...
  systemmatrix_->Complete();

}  // setup_system_matrix()


/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::setup_rhs()
{
  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::setup_rhs");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));

  // fill the EHL rhs vector rhs_ with the single field rhs
  setup_vector(*rhs_, structure_field()->RHS(), lubrication_->LubricationField()->RHS());
}  // setup_rhs()


/*----------------------------------------------------------------------*
 | solve linear EHL system                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::linear_solve()
{
  // Solve for inc_ = [disi_,prei_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  Core::LinAlg::SolverParams solver_params;
  if (solveradapttol_ and (iter_ > 1))
  {
    solver_params.nonlin_tolerance = normrhs_;
    solver_params.nonlin_tolerance = tolrhs_;
    solver_params.lin_tol_better = solveradaptolbetter_;
  }

  // Dirichlet boundary conditions are already applied to EHL system, i.e. EHL
  // system is prepared for solve, i.e. EHL systemmatrix, EHL rhs, EHL inc
  // --> in prepare_system_for_newton_solve(): done for rhs and diagonal blocks
  // --> in setup_system_matrix() done for off-diagonal blocks k_sl, k_ls

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // default: use block matrix
  if (merge_ehl_blockmatrix_ == false)
  {
    // Infnormscaling: scale system before solving
    scale_system(*systemmatrix_, *rhs_);

    // solve the problem, work is done here!
    solver_params.refactor = true;
    solver_params.reset = iter_ == 1;
    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, solver_params);

    // Infnormscaling: unscale system after solving
    unscale_solution(*systemmatrix_, *iterinc_, *rhs_);

  }  // use block matrix

  else  // (merge_ehl_blockmatrix_ == true)
  {
    // Infnormscaling: scale system before solving
    scale_system(*systemmatrix_, *rhs_);

    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse = systemmatrix_->Merge();

    // standard solver call
    solver_params.refactor = true;
    solver_params.reset = iter_ == 1;
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, solver_params);

    // Infnormscaling: unscale system after solving
    unscale_solution(*systemmatrix_, *iterinc_, *rhs_);
  }  // MergeBlockMatrix

}  // linear_solve()


/*----------------------------------------------------------------------*
 | setup vector of the structure and lubrication field      wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::setup_vector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> lv)
{
  // extract dofs of the two fields
  // and put the structural/lubrication field vector into the global vector f
  // noticing the block number
  extractor()->insert_vector(*sv, 0, f);
  extractor()->insert_vector(*lv, 1, f);

}  // setup_vector()


/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)           wirtz 01/16 |
 *----------------------------------------------------------------------*/
bool EHL::Monolithic::Converged()
{
  // check for single norms
  bool convrhs = false;
  bool convinc = false;
  bool convstrrhs = false;
  bool convdisp = false;
  bool convlubricationrhs = false;
  bool convpre = false;

  // ----------------------------------------------------------- EHL test
  // residual EHL forces
  switch (normtyperhs_)
  {
    case Inpar::EHL::convnorm_abs:
      convrhs = normrhs_ < tolrhs_;
      break;
    case Inpar::EHL::convnorm_rel:
      convrhs = normrhs_ < std::max(tolrhs_ * normrhsiter0_, 1.0e-15);
      break;
    case Inpar::EHL::convnorm_mix:
      convrhs = ((normrhs_ < tolrhs_) and (normrhs_ < std::max(normrhsiter0_ * tolrhs_, 1.0e-15)));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces!");
      break;
  }

  // residual EHL increments
  switch (normtypeinc_)
  {
    case Inpar::EHL::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    case Inpar::EHL::convnorm_rel:
      convinc = norminc_ < std::max(norminciter0_ * tolinc_, 1e-15);
      break;
    case Inpar::EHL::convnorm_mix:
      convinc = norminc_ < std::max(norminciter0_ * tolinc_, 1e-15);
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of increments!");
      break;
  }  // switch (normtypeinc_)

  // -------------------------------------------------- single field test
  // ---------------------------------------------------------- structure
  // structural residual forces
  switch (normtypestrrhs_)
  {
    case Inpar::STR::convnorm_abs:
      convstrrhs = normstrrhs_ < tolstrrhs_;
      break;
    case Inpar::STR::convnorm_rel:
      convstrrhs = normstrrhs_ < std::max(normstrrhsiter0_ * tolstrrhs_, 1e-15);
      break;
    case Inpar::STR::convnorm_mix:
      convstrrhs = ((normstrrhs_ < tolstrrhs_) or
                    (normstrrhs_ < std::max(normstrrhsiter0_ * tolstrrhs_, 1e-15)));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces!");
      break;
  }  // switch (normtypestrrhs_)

  // residual displacements
  switch (normtypedisi_)
  {
    case Inpar::STR::convnorm_abs:
      convdisp = normdisi_ < toldisi_;
      break;
    case Inpar::STR::convnorm_rel:
      convdisp = normdisi_ < std::max(normdisiiter0_ * toldisi_, 1e-15);
      break;
    case Inpar::STR::convnorm_mix:
      convdisp =
          ((normdisi_ < toldisi_) or (normdisi_ < std::max(normdisiiter0_ * toldisi_, 1e-15)));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of displacements!");
      break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- lubrication
  // lubrication residual forces
  switch (normtypelubricationrhs_)
  {
    case Inpar::LUBRICATION::convnorm_abs:
      convlubricationrhs = normlubricationrhs_ < tollubricationrhs_;
      break;
    case Inpar::LUBRICATION::convnorm_rel:
      convlubricationrhs = normlubricationrhs_ < normlubricationrhsiter0_ * tollubricationrhs_;
      break;
    case Inpar::LUBRICATION::convnorm_mix:
      convlubricationrhs = ((normlubricationrhs_ < tollubricationrhs_) or
                            (normlubricationrhs_ < normlubricationrhsiter0_ * tollubricationrhs_));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces!");
      break;
  }  // switch (normtypelubricationrhs_)

  // residual pressures
  switch (normtypeprei_)
  {
    case Inpar::LUBRICATION::convnorm_abs:
      convpre = normprei_ < tolprei_;
      break;
    case Inpar::LUBRICATION::convnorm_rel:
      convpre = normprei_ < normpreiiter0_ * tolprei_;
      break;
    case Inpar::LUBRICATION::convnorm_mix:
      convpre = ((normprei_ < tolprei_) or (normprei_ < normpreiiter0_ * tolprei_));
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of pressures!");
      break;
  }  // switch (normtypeprei_)

  // -------------------------------------------------------- convergence
  // combine increment-like and force-like residuals, combine EHL and single
  // field values
  bool conv = false;
  if (combincrhs_ == Inpar::EHL::bop_and)
    conv = convinc and convrhs;
  else if (combincrhs_ == Inpar::EHL::bop_or)
    conv = convinc or convrhs;
  else if (combincrhs_ == Inpar::EHL::bop_coupl_and_singl)
    conv = convinc and convrhs and convdisp and convstrrhs and convpre and convlubricationrhs;
  else if (combincrhs_ == Inpar::EHL::bop_coupl_or_singl)
    conv = (convinc and convrhs) or (convdisp and convstrrhs and convpre and convlubricationrhs);
  else if (combincrhs_ == Inpar::EHL::bop_and_singl)
    conv = convdisp and convstrrhs and convpre and convlubricationrhs;
  else if (combincrhs_ == Inpar::EHL::bop_or_singl)
    conv = (convdisp or convstrrhs or convpre or convlubricationrhs);
  else
    FOUR_C_THROW("Something went terribly wrong with binary operator!");

  // return things
  return conv;

}  // Converged()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file  wirtz 01/16 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::print_newton_iter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and print_screen_evry() and (Step() % print_screen_evry() == 0) and
      printiter_)
  {
    if (iter_ == 1) print_newton_iter_header(stdout);
    print_newton_iter_text(stdout);
  }

}  // print_newton_iter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file  wirtz 01/16 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::print_newton_iter_header(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << "numiter";

  // ---------------------------------------------------------------- EHL
  // different style due relative or absolute error checking
  // displacement
  switch (normtyperhs_)
  {
    case Inpar::EHL::convnorm_abs:
      oss << std::setw(15) << "abs-res-norm";
      break;
    case Inpar::EHL::convnorm_rel:
      oss << std::setw(15) << "rel-res-norm";
      break;
    case Inpar::EHL::convnorm_mix:
      oss << std::setw(15) << "mix-res-norm";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::EHL::convnorm_abs:
      oss << std::setw(15) << "abs-inc-norm";
      break;
    case Inpar::EHL::convnorm_rel:
      oss << std::setw(15) << "rel-inc-norm";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  // -------------------------------------------------- single field test
  // ---------------------------------------------------------- structure
  switch (normtypestrrhs_)
  {
    case Inpar::STR::convnorm_rel:
      oss << std::setw(18) << "rel-str-res-norm";
      break;
    case Inpar::STR::convnorm_abs:
      oss << std::setw(18) << "abs-str-res-norm";
      break;
    case Inpar::STR::convnorm_mix:
      oss << std::setw(18) << "mix-str-res-norm";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypestrrhs_)

  switch (normtypedisi_)
  {
    case Inpar::STR::convnorm_rel:
      oss << std::setw(16) << "rel-dis-norm";
      break;
    case Inpar::STR::convnorm_abs:
      oss << std::setw(16) << "abs-dis-norm";
      break;
    case Inpar::STR::convnorm_mix:
      oss << std::setw(16) << "mix-dis-norm";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- lubrication
  switch (normtypelubricationrhs_)
  {
    case Inpar::LUBRICATION::convnorm_rel:
      oss << std::setw(18) << "rel-lub-res-norm";
      break;
    case Inpar::LUBRICATION::convnorm_abs:
      oss << std::setw(18) << "abs-lub-res-norm";
      break;
    case Inpar::LUBRICATION::convnorm_mix:
      oss << std::setw(18) << "mix-lub-res-norm";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypelubricationrhs_)

  switch (normtypeprei_)
  {
    case Inpar::LUBRICATION::convnorm_rel:
      oss << std::setw(16) << "rel-pre-norm";
      break;
    case Inpar::LUBRICATION::convnorm_abs:
      oss << std::setw(16) << "abs-pre-norm";
      break;
    case Inpar::LUBRICATION::convnorm_mix:
      oss << std::setw(16) << "mix-pre-norm";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypeprei_)

  if (mortaradapter_->HasContact())
  {
    oss << std::setw(16) << "L2-Cont-Res";
    oss << std::setw(16) << "L2-Cont-Incr";
    oss << std::setw(11) << "#active";
    oss << std::setw(10) << "#slip";
  }

  // add solution time
  oss << std::setw(12) << "ts";
  oss << std::setw(12) << "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == nullptr) FOUR_C_THROW("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // print_newton_iter_header()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                 wirtz 01/16 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::print_newton_iter_text(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking

  // ----------------------------------------------- test coupled problem
  switch (normtyperhs_)
  {
    case Inpar::EHL::convnorm_abs:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhs_;
      break;
    case Inpar::EHL::convnorm_rel:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhs_ / normrhsiter0_;
      break;
    case Inpar::EHL::convnorm_mix:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << std::min(normrhs_, normrhs_ / normrhsiter0_);
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::EHL::convnorm_abs:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << norminc_;
      break;
    case Inpar::EHL::convnorm_rel:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << norminc_ / norminciter0_;
      break;
    case Inpar::EHL::convnorm_mix:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << std::min(norminc_, norminc_ / norminciter0_);
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypeinc_)

  // ------------------------------------------------- test single fields
  // ---------------------------------------------------------- structure
  // different style due relative or absolute error checking
  // displacement
  switch (normtypestrrhs_)
  {
    case Inpar::STR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normstrrhs_;
      break;
    case Inpar::STR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << normstrrhs_ / normstrrhsiter0_;
      break;
    case Inpar::STR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normstrrhs_, normstrrhs_ / normstrrhsiter0_);
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypestrrhs_)

  switch (normtypedisi_)
  {
    case Inpar::STR::convnorm_abs:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_;
      break;
    case Inpar::STR::convnorm_rel:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_ / normdisiiter0_;
      break;
    case Inpar::STR::convnorm_mix:
      oss << std::setw(16) << std::setprecision(5) << std::scientific
          << std::min(normdisi_, normdisi_ / normdisiiter0_);
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- lubrication
  switch (normtypelubricationrhs_)
  {
    case Inpar::LUBRICATION::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normlubricationrhs_;
      break;
    case Inpar::LUBRICATION::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << normlubricationrhs_ / normlubricationrhsiter0_;
      break;
    case Inpar::LUBRICATION::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normlubricationrhs_, normlubricationrhs_ / normlubricationrhsiter0_);
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypelubricationrhs_)

  switch (normtypeprei_)
  {
    case Inpar::LUBRICATION::convnorm_abs:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normprei_;
      break;
    case Inpar::LUBRICATION::convnorm_rel:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normprei_ / normpreiiter0_;
      break;
    case Inpar::LUBRICATION::convnorm_mix:
      oss << std::setw(16) << std::setprecision(5) << std::scientific
          << std::min(normprei_, normprei_ / normpreiiter0_);
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }  // switch (normtypeprei_)


  if (mortaradapter_->HasContact())
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << mortaradapter_->ContactRes();
    oss << std::setw(16) << std::setprecision(5) << std::scientific
        << mortaradapter_->ContactIncr();
    oss << std::setw(11) << mortaradapter_->ActiveContact();
    oss << std::setw(10) << mortaradapter_->SlipContact();
  }

  // add solution time of to print to screen
  oss << std::setw(12) << std::setprecision(2) << std::scientific << dtsolve_;
  oss << std::setw(12) << std::setprecision(2) << std::scientific
      << timernewton_.totalElapsedTime(true);

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == nullptr) FOUR_C_THROW("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // print_newton_iter_text


/*----------------------------------------------------------------------*
 | print statistics of converged NRI                        wirtz 01/16 |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::print_newton_conv()
{
  // somebody did the door
  return;
}  // print_newton_conv()



/*----------------------------------------------------------------------*
 | evaluate lubrication-mechanical system matrix at state   wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::apply_lubrication_coupl_matrix(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> matheight,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> matvel)
{
  // create the parameters for the discretization
  Teuchos::ParameterList lparams;
  // action for elements
  LUBRICATION::Action action = LUBRICATION::calc_lubrication_coupltang;
  lparams.set<int>("action", action);
  // other parameters that might be needed by the elements
  lparams.set("delta time", Dt());

  // provide bool whether lubrication grid is displaced or not
  lparams.set<bool>("isale", true);

  lubrication_->LubricationField()->discretization()->ClearState(true);
  // set the variables that are needed by the elements
  lubrication_->LubricationField()->discretization()->set_state(
      0, "prenp", lubrication_->LubricationField()->Prenp());

  set_struct_solution(structure_->Dispnp());

  // build specific assemble strategy for the lubrication-mechanical system matrix
  // from the point of view of lubrication_->LubricationField:
  // thermdofset = 0, structdofset = 1
  Core::FE::AssembleStrategy lubricationstrategy(0,  // pressure dofset for row
      1,                                             // displacement dofset for column
      matheight,                                     // lubrication-mechanical matrix
      matvel,                                        // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the lubrication-mechanical system matrix on the lubrication element
  lubrication_->LubricationField()->discretization()->evaluate(lparams, lubricationstrategy);

  matheight->Complete(*(lubrication_->LubricationField()->discretization()->dof_row_map(1)),
      *(lubrication_->LubricationField()->discretization()->dof_row_map(0)));

  matvel->Complete(*(lubrication_->LubricationField()->discretization()->dof_row_map(1)),
      *(lubrication_->LubricationField()->discretization()->dof_row_map(0)));

  lubrication_->LubricationField()->discretization()->ClearState(true);

  return;
}  // apply_lubrication_coupl_matrix()

/*----------------------------------------------------------------------*
 | map containing the dofs with Dirichlet BC                wirtz 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> EHL::Monolithic::combined_dbc_map()
{
  const Teuchos::RCP<const Epetra_Map> scondmap =
      structure_field()->get_dbc_map_extractor()->cond_map();
  const Teuchos::RCP<const Epetra_Map> lcondmap =
      lubrication_->LubricationField()->get_dbc_map_extractor()->cond_map();
  Teuchos::RCP<Epetra_Map> condmap = Core::LinAlg::MergeMap(scondmap, lcondmap, false);
  return condmap;

}  // combined_dbc_map()

/*----------------------------------------------------------------------*
 | scale system, i.e. apply infnorm scaling to linear       wirtz 01/16 |
 | block system before solving system                                   |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::scale_system(Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  // should we scale the system?
  const bool scaling_infnorm = (bool)Core::UTILS::IntegralValue<int>(ehldynmono_, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if ((A->LeftScale(*srowsum_)) or (A->RightScale(*scolsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_)))
      FOUR_C_THROW("structure scaling failed");

    A = mat.Matrix(1, 1).EpetraMatrix();
    lrowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    lcolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*lrowsum_);
    A->InvColSums(*lcolsum_);
    if ((A->LeftScale(*lrowsum_)) or (A->RightScale(*lcolsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->LeftScale(*lrowsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->RightScale(*lcolsum_)))
      FOUR_C_THROW("lubrication scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = extractor()->extract_vector(b, 0);
    Teuchos::RCP<Epetra_Vector> lx = extractor()->extract_vector(b, 1);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (lx->Multiply(1.0, *lrowsum_, *lx, 0.0)) FOUR_C_THROW("lubrication scaling failed");

    extractor()->insert_vector(*sx, 0, b);
    extractor()->insert_vector(*lx, 1, b);
  }
}  // scale_system


/*----------------------------------------------------------------------*
 | unscale solution after solving the linear system         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::unscale_solution(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const bool scaling_infnorm = (bool)Core::UTILS::IntegralValue<int>(ehldynmono_, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = extractor()->extract_vector(x, 0);
    Teuchos::RCP<Epetra_Vector> ly = extractor()->extract_vector(x, 1);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ly->Multiply(1.0, *lcolsum_, *ly, 0.0)) FOUR_C_THROW("lubrication scaling failed");

    extractor()->insert_vector(*sy, 0, x);
    extractor()->insert_vector(*ly, 1, x);

    Teuchos::RCP<Epetra_Vector> sx = extractor()->extract_vector(b, 0);
    Teuchos::RCP<Epetra_Vector> lx = extractor()->extract_vector(b, 1);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (lx->ReciprocalMultiply(1.0, *lrowsum_, *lx, 0.0))
      FOUR_C_THROW("lubrication scaling failed");

    extractor()->insert_vector(*sx, 0, b);
    extractor()->insert_vector(*lx, 1, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if ((A->LeftScale(*srowsum_)) or (A->RightScale(*scolsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_)))
      FOUR_C_THROW("structure scaling failed");

    A = mat.Matrix(1, 1).EpetraMatrix();
    lrowsum_->Reciprocal(*lrowsum_);
    lcolsum_->Reciprocal(*lcolsum_);
    if ((A->LeftScale(*lrowsum_)) or (A->RightScale(*lcolsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->LeftScale(*lrowsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->RightScale(*lcolsum_)))
      FOUR_C_THROW("lubrication scaling failed");

  }  // if (scaling_infnorm)

}  // unscale_solution()


/*----------------------------------------------------------------------*
 | calculate vector norm                                    wirtz 01/16 |
 *----------------------------------------------------------------------*/
double EHL::Monolithic::calculate_vector_norm(
    const enum Inpar::EHL::VectorNorm norm, const Teuchos::RCP<Epetra_Vector> vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == Inpar::EHL::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == Inpar::EHL::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == Inpar::EHL::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == Inpar::EHL::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == Inpar::EHL::norm_l1_scaled)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm / ((double)vect->GlobalLength());
  }
  else
  {
    FOUR_C_THROW("Cannot handle vector norm");
    return 0;
  }
}  // calculate_vector_norm()


/*----------------------------------------------------------------------*
 | set parameters for EHL remaining constant over whole     wirtz 01/16 |
 | simulation                                                           |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::set_default_parameters()
{
  // time parameters
  // call the EHL parameter list
  const Teuchos::ParameterList& ldyn = Global::Problem::Instance()->lubrication_dynamic_params();

  // get the parameters for the Newton iteration
  itermax_ = ehldyn_.get<int>("ITEMAX");
  itermin_ = ehldyn_.get<int>("ITEMIN");

  // what kind of norm do we wanna test for coupled EHL problem
  normtypeinc_ = Core::UTILS::IntegralValue<Inpar::EHL::ConvNorm>(ehldynmono_, "NORM_INC");
  normtyperhs_ = Core::UTILS::IntegralValue<Inpar::EHL::ConvNorm>(ehldynmono_, "NORM_RESF");
  // what kind of norm do we wanna test for the single fields
  normtypedisi_ = Core::UTILS::IntegralValue<Inpar::STR::ConvNorm>(sdyn_, "NORM_DISP");
  normtypestrrhs_ = Core::UTILS::IntegralValue<Inpar::STR::ConvNorm>(sdyn_, "NORM_RESF");
  enum Inpar::STR::VectorNorm striternorm =
      Core::UTILS::IntegralValue<Inpar::STR::VectorNorm>(sdyn_, "ITERNORM");
  normtypeprei_ = Core::UTILS::IntegralValue<Inpar::LUBRICATION::ConvNorm>(ldyn, "NORM_PRE");
  normtypelubricationrhs_ =
      Core::UTILS::IntegralValue<Inpar::LUBRICATION::ConvNorm>(ldyn, "NORM_RESF");
  enum Inpar::LUBRICATION::VectorNorm lubricationiternorm =
      Core::UTILS::IntegralValue<Inpar::LUBRICATION::VectorNorm>(ldyn, "ITERNORM");
  // in total when do we reach a converged state for complete problem
  combincrhs_ = Core::UTILS::IntegralValue<Inpar::EHL::BinaryOp>(ehldynmono_, "NORMCOMBI_RESFINC");

  switch (combincrhs_)
  {
    case Inpar::EHL::bop_and:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of EHL:\n res, inc with 'AND'." << std::endl;
      break;
    }
    case Inpar::EHL::bop_or:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of EHL:\n res, inc with 'OR'." << std::endl;
      break;
    }
    case Inpar::EHL::bop_coupl_and_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of EHL:\n res, inc, str-res, lub-res, dis, pre with 'AND'."
                  << std::endl;
      break;
    }
    case Inpar::EHL::bop_coupl_or_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of EHL:\n (res, inc) or (str-res, lub-res, dis, pre)."
                  << std::endl;
      break;
    }
    case Inpar::EHL::bop_and_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of EHL:\n str-res, lub-res, dis, pre with 'AND'."
                  << std::endl;
      break;
    }
    case Inpar::EHL::bop_or_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of EHL:\n str-res, lub-res, dis, pre with 'OR'."
                  << std::endl;
      break;
    }
    default:
    {
      FOUR_C_THROW("Something went terribly wrong with binary operator!");
      break;
    }
  }  // switch (combincrhs_)

  // convert the single field norms to be used within EHL
  // what norm is used for structure
  switch (striternorm)
  {
    case Inpar::STR::norm_l1:
      iternormstr_ = Inpar::EHL::norm_l1;
      break;
    case Inpar::STR::norm_l2:
      iternormstr_ = Inpar::EHL::norm_l2;
      break;
    case Inpar::STR::norm_rms:
      iternormstr_ = Inpar::EHL::norm_rms;
      break;
    case Inpar::STR::norm_inf:
      iternormstr_ = Inpar::EHL::norm_inf;
      break;
    case Inpar::STR::norm_vague:
    default:
      FOUR_C_THROW("STR norm is not determined");
      break;
  }  // switch (striternorm)

  // what norm is used for lubrication
  switch (lubricationiternorm)
  {
    case Inpar::LUBRICATION::norm_l1:
      iternormlubrication_ = Inpar::EHL::norm_l1;
      break;
    case Inpar::LUBRICATION::norm_l2:
      iternormlubrication_ = Inpar::EHL::norm_l2;
      break;
    case Inpar::LUBRICATION::norm_rms:
      iternormlubrication_ = Inpar::EHL::norm_rms;
      break;
    case Inpar::LUBRICATION::norm_inf:
      iternormlubrication_ = Inpar::EHL::norm_inf;
      break;
    case Inpar::LUBRICATION::norm_vague:
    default:
    {
      FOUR_C_THROW("LUBRICATION norm is not determined.");
      break;
    }
  }  // switch (lubricationiternorm)

  // if scaled L1-norm is wished to be used
  if ((iternorm_ == Inpar::EHL::norm_l1_scaled) and
      ((combincrhs_ == Inpar::EHL::bop_coupl_and_singl) or
          (combincrhs_ == Inpar::EHL::bop_coupl_or_singl)))
  {
    iternormstr_ = Inpar::EHL::norm_l1_scaled;
    iternormlubrication_ = Inpar::EHL::norm_l1_scaled;
  }

  // test the EHL-residual and the EHL-increment
  tolinc_ = ehldynmono_.get<double>("TOLINC");
  tolrhs_ = ehldynmono_.get<double>("CONVTOL");

  // get the single field tolerances from this field itselves
  toldisi_ = sdyn_.get<double>("TOLDISP");
  tolstrrhs_ = sdyn_.get<double>("TOLRES");
  tolprei_ = ldyn.get<double>("TOLPRE");
  tollubricationrhs_ = ldyn.get<double>("TOLRES");

  // initialise norms for coupled EHL problem
  normrhs_ = 0.0;
  normrhsiter0_ = 0.0;
  norminc_ = 0.0;
  norminciter0_ = 0.0;

  // initialise norms for single field tests
  normdisi_ = 0.0;
  normstrrhs_ = 0.0;
  normstrrhsiter0_ = 0.0;
  normprei_ = 0.0;
  normlubricationrhs_ = 0.0;
  normlubricationrhsiter0_ = 0.0;

  return;

}  // SetDefaultParameter()

/*----------------------------------------------------------------------*
 | calculate stresses, strains, energies                    wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::prepare_output(bool force_prepare)
{
  // prepare output (i.e. calculate stresses, strains, energies)
  structure_field()->prepare_output(force_prepare);

  // reset states
  structure_field()->discretization()->ClearState(true);
}
/*----------------------------------------------------------------------*/


void EHL::Monolithic::lin_pressure_force_disp(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dd,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dd)
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> p_deriv_normal =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->NderivMatrix()));
  Teuchos::RCP<Epetra_Vector> p_full =
      Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->dof_row_map(1)));
  if (lubrimaptransform_->Apply(*lubrication_->LubricationField()->Prenp(), *p_full))
    FOUR_C_THROW("apply failed");
  Teuchos::RCP<Epetra_Vector> p_exp =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  p_exp = ada_strDisp_to_lubDisp_->SlaveToMaster(p_full);
  if (p_deriv_normal->LeftScale(*p_exp)) FOUR_C_THROW("leftscale failed");
  if (p_deriv_normal->Scale(-1.)) FOUR_C_THROW("scale failed");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixD(), true, *p_deriv_normal, false, false, false, true);
  if (tmp.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  ds_dd->Add(*tmp, false, +1., 1.);

  tmp = Teuchos::null;
  tmp = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixM(), true, *p_deriv_normal, false, false, false, true);
  if (tmp.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  dm_dd->Add(*tmp, false, -1., 1.);

  return;
}

void EHL::Monolithic::lin_poiseuille_force_disp(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dd,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dd)
{
  Teuchos::RCP<Epetra_Vector> p_int =
      ada_strDisp_to_lubPres_->SlaveToMaster(lubrication_->LubricationField()->Prenp());
  Teuchos::RCP<Epetra_Vector> p_int_full =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  Core::LinAlg::Export(*p_int, *p_int_full);

  Teuchos::RCP<Epetra_Vector> nodal_gap =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *nodal_gap))
    FOUR_C_THROW("multiply failed");

  Teuchos::RCP<Epetra_Vector> grad_p =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (mortaradapter_->SurfGradMatrix()->Apply(*p_int_full, *grad_p)) FOUR_C_THROW("apply failed");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> deriv_Poiseuille = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*mortaradapter_->SlaveDofMap(), 81, false, false));

  Teuchos::RCP<Core::LinAlg::SparseMatrix> derivH_gradP =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->Nodal_GapDeriv()));
  if (derivH_gradP->LeftScale(*grad_p)) FOUR_C_THROW("leftscale failed");
  deriv_Poiseuille->Add(*derivH_gradP, false, -.5, 1.);

  Teuchos::RCP<Epetra_Vector> p_int_full_col =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->Interface()->SlaveColDofs()));
  Core::LinAlg::Export(*p_int_full, *p_int_full_col);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> h_derivGrad_nodalP =
      mortaradapter_->assemble_surf_grad_deriv(p_int_full_col);
  if (h_derivGrad_nodalP->LeftScale(*nodal_gap)) FOUR_C_THROW("leftscale failed");
  deriv_Poiseuille->Add(*h_derivGrad_nodalP, false, -.5, 1.);

  deriv_Poiseuille->Complete(*mortaradapter_->SMdofMap(), *mortaradapter_->SlaveDofMap());

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixD(), true, *deriv_Poiseuille, false, false, false, true);
  if (tmp.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  ds_dd->Add(*tmp, false, +1., 1.);

  tmp = Teuchos::null;
  tmp = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixM(), true, *deriv_Poiseuille, false, false, false, true);
  if (tmp.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  dm_dd->Add(*tmp, false, +1., 1.);
}

void EHL::Monolithic::LinCouetteForceDisp(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dd,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dd)
{
  const int ndim = Global::Problem::Instance()->NDim();
  Core::FE::Discretization& lub_dis = *lubrication_->LubricationField()->discretization();
  Teuchos::RCP<Epetra_Vector> visc_vec =
      Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->dof_row_map(1)));
  for (int i = 0; i < lub_dis.NodeRowMap()->NumMyElements(); ++i)
  {
    Core::Nodes::Node* lnode = lub_dis.lRowNode(i);
    if (!lnode) FOUR_C_THROW("node not found");
    const double p = lubrication_->LubricationField()->Prenp()->operator[](
        lubrication_->LubricationField()->Prenp()->Map().LID(lub_dis.Dof(0, lnode, 0)));

    Teuchos::RCP<Core::Mat::Material> mat = lnode->Elements()[0]->Material(0);
    if (mat.is_null()) FOUR_C_THROW("null pointer");
    Teuchos::RCP<Mat::LubricationMat> lmat =
        Teuchos::rcp_dynamic_cast<Mat::LubricationMat>(mat, true);
    const double visc = lmat->ComputeViscosity(p);

    for (int d = 0; d < ndim; ++d) visc_vec->ReplaceGlobalValue(lub_dis.Dof(1, lnode, d), 0, visc);
  }
  Teuchos::RCP<Epetra_Vector> visc_vec_str = ada_strDisp_to_lubDisp_->SlaveToMaster(visc_vec);

  Teuchos::RCP<Epetra_Vector> height =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *height))
    FOUR_C_THROW("multiply failed");
  Teuchos::RCP<Epetra_Vector> h_inv =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (h_inv->Reciprocal(*height)) FOUR_C_THROW("Reciprocal failed");

  Teuchos::RCP<Epetra_Vector> hinv_visc =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  hinv_visc->Multiply(1., *h_inv, *visc_vec_str, 0.);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> deriv_Couette = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*mortaradapter_->SlaveDofMap(), 81, false, false));

  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->RelTangVelDeriv()));
    tmp->LeftScale(*hinv_visc);
    tmp->Scale(-1.);
    deriv_Couette->Add(*tmp, false, 1., 1.);
  }
  {
    Teuchos::RCP<Epetra_Vector> hinv_hinv_visc =
        Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
    hinv_hinv_visc->Multiply(1., *h_inv, *hinv_visc, 0.);
    Teuchos::RCP<Epetra_Vector> hinv_hinv_visc_vel =
        Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
    hinv_hinv_visc_vel->Multiply(1., *hinv_hinv_visc, *mortaradapter_->RelTangVel(), 0.);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->Nodal_GapDeriv()));
    tmp->LeftScale(*hinv_hinv_visc_vel);
    deriv_Couette->Add(*tmp, false, 1., 1.);
  }
  deriv_Couette->Complete(*mortaradapter_->SMdofMap(), *mortaradapter_->SlaveDofMap());

  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixD(), true, *deriv_Couette, false, false, false, true);
  if (tmp.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  ds_dd->Add(*tmp, false, +1., 1.);

  tmp = Teuchos::null;
  tmp = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixM(), true, *deriv_Couette, false, false, false, true);
  if (tmp.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  dm_dd->Add(*tmp, false, -1., 1.);
}

void EHL::Monolithic::lin_pressure_force_pres(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dp,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dp)
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> tmp = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*mortaradapter_->SlaveDofMap(), 81, false, false));

  Core::LinAlg::MatrixRowTransform().operator()(*lubrimaptransform_, 1.,
      Core::Adapter::CouplingSlaveConverter(*ada_strDisp_to_lubDisp_), *tmp, false);

  tmp->Complete(*lubrication_->LubricationField()->dof_row_map(0), *mortaradapter_->SlaveDofMap());

  if (tmp->LeftScale(*mortaradapter_->Normals())) FOUR_C_THROW("leftscale failed");
  if (tmp->Scale(-1.)) FOUR_C_THROW("scale failed");

  Teuchos::RCP<Core::LinAlg::SparseMatrix> a = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixD(), true, *tmp, false, false, false, true);
  if (a.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  ds_dp->Add(*a, false, +1., 1.);

  a = Teuchos::null;
  a = Core::LinAlg::MLMultiply(
      *mortaradapter_->GetMortarMatrixM(), true, *tmp, false, false, false, true);
  if (a.is_null()) FOUR_C_THROW("MLMULTIPLY failed");
  dm_dp->Add(*a, false, -1., 1.);
}

void EHL::Monolithic::lin_poiseuille_force_pres(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dp,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dp)
{
  Teuchos::RCP<Epetra_Vector> nodal_gap =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *nodal_gap))
    FOUR_C_THROW("multiply failed");

  Core::LinAlg::SparseMatrix m(*mortaradapter_->SurfGradMatrix());
  m.LeftScale(*nodal_gap);
  m.Scale(-.5);

  {
    Teuchos::RCP<const Epetra_Map> r = mortaradapter_->SlaveDofMap();
    Teuchos::RCP<const Epetra_Map> d = lubrication_->LubricationField()->dof_row_map(0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> a = Core::LinAlg::MLMultiply(
        *mortaradapter_->GetMortarMatrixD(), true, m, false, true, false, true);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> b =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(a->RowMap(), 81, false, false));

    Core::LinAlg::MatrixColTransform().operator()(a->RowMap(), a->ColMap(), *a, 1.,
        Core::Adapter::CouplingMasterConverter(*ada_strDisp_to_lubPres_), *b, true, false);
    b->Complete(*d, *r);

    ds_dp->UnComplete();
    ds_dp->Add(*b, false, 1., 1.);
    ds_dp->Complete(*d, *r);
  }

  {
    Teuchos::RCP<const Epetra_Map> r = mortaradapter_->MasterDofMap();
    Teuchos::RCP<const Epetra_Map> d = lubrication_->LubricationField()->dof_row_map(0);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> a = Core::LinAlg::MLMultiply(
        *mortaradapter_->GetMortarMatrixM(), true, m, false, true, false, true);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> b =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(a->RowMap(), 81, false, false));

    Core::LinAlg::MatrixColTransform().operator()(a->RowMap(), a->ColMap(), *a, 1.,
        Core::Adapter::CouplingMasterConverter(*ada_strDisp_to_lubPres_), *b, true, false);
    b->Complete(*d, *r);

    dm_dp->UnComplete();
    dm_dp->Add(*b, false, 1., 1.);
    dm_dp->Complete(*d, *r);
  }
  return;
}

void EHL::Monolithic::LinCouetteForcePres(Teuchos::RCP<Core::LinAlg::SparseMatrix>& ds_dp,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& dm_dp)
{
  const int ndim = Global::Problem::Instance()->NDim();
  const Teuchos::RCP<const Epetra_Vector> relVel = mortaradapter_->RelTangVel();
  Teuchos::RCP<Epetra_Vector> height =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *height))
    FOUR_C_THROW("multiply failed");
  Teuchos::RCP<Epetra_Vector> h_inv =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (h_inv->Reciprocal(*height)) FOUR_C_THROW("Reciprocal failed");
  Teuchos::RCP<Epetra_Vector> hinv_relV =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  hinv_relV->Multiply(1., *h_inv, *relVel, 0.);

  Core::FE::Discretization& lub_dis = *lubrication_->LubricationField()->discretization();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> dVisc_dp =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(lub_dis.dof_row_map(1)), 81));

  for (int i = 0; i < lub_dis.NodeRowMap()->NumMyElements(); ++i)
  {
    Core::Nodes::Node* lnode = lub_dis.lRowNode(i);
    if (!lnode) FOUR_C_THROW("node not found");
    const double p = lubrication_->LubricationField()->Prenp()->operator[](
        lubrication_->LubricationField()->Prenp()->Map().LID(lub_dis.Dof(0, lnode, 0)));

    Teuchos::RCP<Core::Mat::Material> mat = lnode->Elements()[0]->Material(0);
    if (mat.is_null()) FOUR_C_THROW("null pointer");
    Teuchos::RCP<Mat::LubricationMat> lmat =
        Teuchos::rcp_dynamic_cast<Mat::LubricationMat>(mat, true);
    const double visc = lmat->ComputeViscosity(p);
    const double dvisc_dp = lmat->compute_viscosity_deriv(p, visc);

    for (int d = 0; d < ndim; ++d)
      dVisc_dp->Assemble(dvisc_dp, lub_dis.Dof(1, lnode, d), lub_dis.Dof(0, lnode, 0));
  }
  dVisc_dp->Complete(*lub_dis.dof_row_map(0), *lub_dis.dof_row_map(1));

  Teuchos::RCP<Core::LinAlg::SparseMatrix> dVisc_str_dp =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*mortaradapter_->SlaveDofMap(), 81));

  Core::LinAlg::MatrixRowTransform().operator()(*dVisc_dp, 1.,
      Core::Adapter::CouplingSlaveConverter(*ada_strDisp_to_lubDisp_), *dVisc_str_dp, false);

  dVisc_str_dp->Complete(*lub_dis.dof_row_map(0), *mortaradapter_->SlaveDofMap());

  dVisc_str_dp->LeftScale(*hinv_relV);

  {
    Teuchos::RCP<const Epetra_Map> r = mortaradapter_->SlaveDofMap();
    Teuchos::RCP<const Epetra_Map> d = lubrication_->LubricationField()->dof_row_map(0);

    ds_dp->UnComplete();
    ds_dp->Add(*Core::LinAlg::MLMultiply(*mortaradapter_->GetMortarMatrixD(), true, *dVisc_str_dp,
                   false, true, false, true),
        false, 1., 1.);
    ds_dp->Complete(*d, *r);
  }

  {
    Teuchos::RCP<const Epetra_Map> r = mortaradapter_->MasterDofMap();
    Teuchos::RCP<const Epetra_Map> d = lubrication_->LubricationField()->dof_row_map(0);

    dm_dp->UnComplete();
    dm_dp->Add(*Core::LinAlg::MLMultiply(*mortaradapter_->GetMortarMatrixM(), true, *dVisc_str_dp,
                   false, true, false, true),
        false, -1., 1.);
    dm_dp->Complete(*d, *r);
  }
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::apply_dbc()
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ss =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(systemmatrix_->Matrix(0, 0).EpetraMatrix(),
          Core::LinAlg::Copy, true, false, Core::LinAlg::SparseMatrix::CRS_MATRIX));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_sl =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(systemmatrix_->Matrix(0, 1)));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ls =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(systemmatrix_->Matrix(1, 0)));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> k_ll =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(systemmatrix_->Matrix(1, 1)));
  k_ss->ApplyDirichlet(*structure_field()->get_dbc_map_extractor()->cond_map(), true);
  k_sl->ApplyDirichlet(*structure_field()->get_dbc_map_extractor()->cond_map(), false);
  k_ls->ApplyDirichlet(
      *lubrication_->LubricationField()->get_dbc_map_extractor()->cond_map(), false);
  k_ll->ApplyDirichlet(
      *lubrication_->LubricationField()->get_dbc_map_extractor()->cond_map(), true);

  if (inf_gap_toggle_lub_ != Teuchos::null)
  {
    k_ls->ApplyDirichlet(*inf_gap_toggle_lub_, false);
    k_ll->ApplyDirichlet(*inf_gap_toggle_lub_, true);
  }

  systemmatrix_->UnComplete();
  systemmatrix_->Assign(0, 0, Core::LinAlg::View, *k_ss);
  systemmatrix_->Assign(0, 1, Core::LinAlg::View, *k_sl);
  systemmatrix_->Assign(1, 0, Core::LinAlg::View, *k_ls);
  systemmatrix_->Assign(1, 1, Core::LinAlg::View, *k_ll);
  systemmatrix_->Complete();


  Core::LinAlg::apply_dirichlet_to_system(
      *rhs_, *zeros_, *structure_field()->get_dbc_map_extractor()->cond_map());
  Core::LinAlg::apply_dirichlet_to_system(
      *rhs_, *zeros_, *lubrication_->LubricationField()->get_dbc_map_extractor()->cond_map());

  if (inf_gap_toggle_lub_ != Teuchos::null)
    for (int i = 0; i < inf_gap_toggle_lub_->MyLength(); ++i)
      if (abs(inf_gap_toggle_lub_->operator[](i)) > 1.e-12)
        rhs_->ReplaceGlobalValue(inf_gap_toggle_lub_->Map().GID(i), 0, 0.);
}

FOUR_C_NAMESPACE_CLOSE
