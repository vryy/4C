/*----------------------------------------------------------------------*/
/*! \file
\brief Base class for monolithic fluid-fluid-fsi algorithm
 using XFEM (without NOX)

\level 2

*----------------------------------------------------------------------*/


#include "4C_fsi_monolithic_nonox.hpp"

#include "4C_adapter_ale_xffsi.hpp"
#include "4C_adapter_fld_fluid_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ale.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
// constructor (public)
/*----------------------------------------------------------------------*/
FSI::MonolithicNoNOX::MonolithicNoNOX(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicBase(comm, timeparams), zeros_(Teuchos::null)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // use taylored fluid- and ALE-wrappers
  fluid_ = Teuchos::rcp_dynamic_cast<Adapter::FluidFluidFSI>(MonolithicBase::fluid_field());
  ale_ = Teuchos::rcp_dynamic_cast<Adapter::AleXFFsiWrapper>(MonolithicBase::ale_field());

  // enable debugging
  if (Core::UTILS::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") == 1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(structure_field()->discretization()));
    // fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(fluid_field()->discretization()));
  }

  std::string s = Global::Problem::Instance()->OutputControlFile()->FileName();
  s.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(s.c_str()));
  itermax_ = fsimono.get<int>("ITEMAX");
  normtypeinc_ = Core::UTILS::IntegralValue<Inpar::FSI::ConvNorm>(fsimono, "NORM_INC");
  normtypefres_ = Core::UTILS::IntegralValue<Inpar::FSI::ConvNorm>(fsimono, "NORM_RESF");
  combincfres_ = Core::UTILS::IntegralValue<Inpar::FSI::BinaryOp>(fsimono, "NORMCOMBI_RESFINC");
  tolinc_ = fsimono.get<double>("CONVTOL");
  tolfres_ = fsimono.get<double>("CONVTOL");

  tol_dis_res_l2_ = fsimono.get<double>("TOL_DIS_RES_L2");
  tol_dis_res_inf_ = fsimono.get<double>("TOL_DIS_RES_INF");
  tol_dis_inc_l2_ = fsimono.get<double>("TOL_DIS_INC_L2");
  tol_dis_inc_inf_ = fsimono.get<double>("TOL_DIS_INC_INF");
  tol_fsi_res_l2_ = fsimono.get<double>("TOL_FSI_RES_L2");
  tol_fsi_res_inf_ = fsimono.get<double>("TOL_FSI_RES_INF");
  tol_fsi_inc_l2_ = fsimono.get<double>("TOL_FSI_INC_L2");
  tol_fsi_inc_inf_ = fsimono.get<double>("TOL_FSI_INC_INF");
  tol_pre_res_l2_ = fsimono.get<double>("TOL_PRE_RES_L2");
  tol_pre_res_inf_ = fsimono.get<double>("TOL_PRE_RES_INF");
  tol_pre_inc_l2_ = fsimono.get<double>("TOL_PRE_INC_L2");
  tol_pre_inc_inf_ = fsimono.get<double>("TOL_PRE_INC_INF");
  tol_vel_res_l2_ = fsimono.get<double>("TOL_VEL_RES_L2");
  tol_vel_res_inf_ = fsimono.get<double>("TOL_VEL_RES_INF");
  tol_vel_inc_l2_ = fsimono.get<double>("TOL_VEL_INC_L2");
  tol_vel_inc_inf_ = fsimono.get<double>("TOL_VEL_INC_INF");
  // set tolerances for nonlinear solver
}

void FSI::MonolithicNoNOX::SetupSystem()
{
  const int ndim = Global::Problem::Instance()->NDim();

  Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  Core::Adapter::Coupling& coupsa = structure_ale_coupling();
  Core::Adapter::Coupling& coupfa = fluid_ale_coupling();
  Core::Adapter::Coupling& icoupfa = interface_fluid_ale_coupling();

  // structure to fluid

  coupsf.setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->Interface()->FSICondMap(), *fluid_field()->discretization(),
      fluid_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // structure to ale

  coupsa.setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
      ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa.setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
      ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    FOUR_C_THROW("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
    FOUR_C_THROW("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_field()->discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = ale_field()->discretization()->NodeRowMap();

  coupfa.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
      *fluidnodemap, *alenodemap, ndim);

  fluid_field()->SetMeshMap(coupfa.MasterDofMap());
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::Timeloop()
{
  while (NotFinished())
  {
    prepare_time_step();
    newton();
    constexpr bool force_prepare = false;
    prepare_output(force_prepare);
    update();
    output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::newton()
{
  // initialise equilibrium loop
  iter_ = 1;

  x_sum_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  x_sum_->PutScalar(0.0);

  // incremental solution vector with length of all FSI dofs
  iterinc_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  iterinc_->PutScalar(0.0);

  zeros_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  zeros_->PutScalar(0.0);

  // residual vector with length of all FSI dofs
  rhs_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  rhs_->PutScalar(0.0);

  firstcall_ = true;

  // equilibrium iteration loop (loop over k)
  while ((iter_ == 1) or ((not converged()) and (iter_ <= itermax_)))
  {
    // compute residual forces #rhs_ and tangent #tang_
    // build linear system stiffness matrix and rhs/force
    // residual for each field

    evaluate(iterinc_);

    // create the linear system
    // J(x_i) \Delta x_i = - R(x_i)
    // create the systemmatrix
    setup_system_matrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      FOUR_C_THROW("Effective tangent matrix must be filled here");
    }

    setup_rhs(*rhs_, firstcall_);

    linear_solve();

    // reset solver tolerance
    solver_->ResetTolerance();

    // build residual and incremental norms
    // for now use for simplicity only L2/Euclidian norm
    build_convergence_norms();

    // print stuff
    print_newton_iter();

    // increment equilibrium loop index
    iter_ += 1;

    firstcall_ = false;

  }  // end while loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((converged()) and (Comm().MyPID() == 0))
  {
    Core::IO::cout << Core::IO::endl;
    Core::IO::cout << "  Newton Converged! " << Core::IO::endl;
  }
  else if (iter_ >= itermax_)
  {
    Core::IO::cout << Core::IO::endl;
    Core::IO::cout << "  Newton unconverged in " << iter_ << " iterations " << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicNoNOX::converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case Inpar::FSI::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    case Inpar::FSI::convnorm_rel:
      convinc =
          (((normstrincL2_ / ns_) < tol_dis_inc_l2_) and ((normstrincInf_) < tol_dis_inc_inf_) and
              ((norminterfaceincL2_ / ni_) < tol_fsi_inc_l2_) and
              ((norminterfaceincInf_) < tol_fsi_inc_inf_) and
              ((normflvelincL2_ / nfv_) < tol_vel_inc_l2_) and
              ((normflvelincInf_) < tol_vel_inc_inf_) and
              ((normflpresincL2_ / nfp_) < tol_pre_inc_l2_) and
              ((normflpresincInf_) < tol_pre_inc_inf_));
      break;
    case Inpar::FSI::convnorm_mix:
      FOUR_C_THROW("not implemented!");
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual values!");
      break;
  }

  // structural, fluid and ale residual forces
  switch (normtypefres_)
  {
    case Inpar::FSI::convnorm_abs:
      convfres = normrhs_ < tolfres_;
      break;
    case Inpar::FSI::convnorm_rel:
      convfres =
          (((normstrrhsL2_ / ns_) < tol_dis_res_l2_) and ((normstrrhsInf_) < tol_dis_res_inf_) and
              ((norminterfacerhsL2_ / ni_) < tol_fsi_res_l2_) and
              ((norminterfacerhsInf_) < tol_fsi_res_inf_) and
              ((normflvelrhsL2_ / nfv_) < tol_vel_res_l2_) and
              ((normflvelrhsInf_) < tol_vel_res_inf_) and
              ((normflpresrhsL2_ / nfp_) < tol_pre_res_l2_) and
              ((normflpresrhsInf_) < tol_pre_res_inf_));
      break;
    case Inpar::FSI::convnorm_mix:
      FOUR_C_THROW("not implemented!");
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces!");
      break;
  }

  // combined
  bool conv = false;
  if (combincfres_ == Inpar::FSI::bop_and)
    conv = convinc and convfres;
  else
    FOUR_C_THROW("Something went wrong!");

  return conv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::linear_solve()
{
  // merge blockmatrix to SparseMatrix and solve
  Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse = systemmatrix_->Merge();

  // apply Dirichlet BCs to system of equations
  if (firstcall_)
    initial_guess(iterinc_);
  else
    iterinc_->PutScalar(0.0);

  Core::LinAlg::apply_dirichlet_to_system(*sparse, *iterinc_, *rhs_, *zeros_, *combined_dbc_map());

#ifndef moresolvers
  const Teuchos::ParameterList& fdyn = Global::Problem::Instance()->FluidDynamicParams();
  const int fluidsolver = fdyn.get<int>("LINEAR_SOLVER");
  solver_ = Teuchos::rcp(
      new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(fluidsolver), Comm()));
#else
  // get UMFPACK...
  Teuchos::ParameterList solverparams = Global::Problem::Instance()->UMFPACKSolverParams();
  solver_ = Teuchos::rcp(new Core::LinAlg::Solver(solverparams, Comm()));
#endif


  // standard solver call
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = iter_ == 1;
  solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, solver_params);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::evaluate(Teuchos::RCP<const Epetra_Vector> step_increment)
{
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  // Save the inner fluid map that includes the background fluid DOF in order to
  // determine a change.
  const Epetra_BlockMap fluidincrementmap = extractor().ExtractVector(step_increment, 1)->Map();

  if (not firstcall_)
  {
    // structure, ale and fluid fields expects the step increment. So
    // we add all of the increments together to build the step
    // increment.
    //
    // The update of the latest increment with iteration increments:
    // x^n+1_i+1 = x^n+1_i + iterinc
    //
    // The update of the latest increment with step increment:
    // x^n+1_i+1 = x^n     + stepinc

    x_sum_->Update(1.0, *step_increment, 1.0);

    extract_field_vectors(x_sum_, sx, fx, ax);

    if (sdbg_ != Teuchos::null)
    {
      sdbg_->NewIteration();
      sdbg_->WriteVector("x", *structure_field()->Interface()->ExtractFSICondVector(sx));
    }
  }

  // Call all fileds evaluate method and assemble rhs and matrices

  {
    structure_field()->Evaluate(sx);
  }

  {
    // ALE field expects the sum of all increments and not the
    // latest increment. It adds the sum of all increments to the
    // displacement of the last time step. So we need to build the
    // sum of all increments and give it to ALE.

    ale_field()->Evaluate(ax);
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);

  {
    fluid_field()->Evaluate(fx);
  }

  if (has_fluid_dof_map_changed(fluidincrementmap)) handle_fluid_dof_map_change_in_newton();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::set_default_parameters(
    const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list)
{
  // monolithic solver settings
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
  // nlParams.set("Preconditioner", "None");
  // nlParams.set("Norm abs F", fsimono.get<double>("CONVTOL"));

  nlParams.set("Max Iterations", fsimono.get<int>("ITEMAX"));
  // nlParams.set("Max Iterations", 1);

  nlParams.set("Norm abs pres", fsimono.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel", fsimono.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsimono.get<double>("CONVTOL"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  // Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");


  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.set<std::string>("Method", "User Defined");
  //   Teuchos::RCP<::NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this,false);
  //   dirParams.set("User Defined Direction Factory",newtonfactory);


  // status tests are expensive, but instructive
  // solverOptions.set<std::string>("Status Test Check Type","Minimal");
  solverOptions.set<std::string>("Status Test Check Type", "Complete");

  // be explicit about linear solver parameters
  lsParams.set<std::string>("Aztec Solver", "GMRES");
  // lsParams.set<std::string>("BiCGStab","GMRES");
  lsParams.set<std::string>("Orthogonalization", "Modified");

  // "r0", "rhs", "norm", "no scaling", "sol"
  lsParams.set<std::string>("Convergence Test", "r0");

  lsParams.set<int>("Size of Krylov Subspace", fsimono.get<int>("KRYLOV_SIZE"));
  lsParams.set<int>("Max Iterations", fsimono.get<int>("KRYLOV_ITEMAX"));
  lsParams.set<std::string>("Preconditioner", "User Defined");
  lsParams.set<int>("Output Frequency", 10);
  lsParams.set<bool>("Output Solver Details", true);

  // adaptive tolerance settings for linear solver
  lsParams.set<double>("base tolerance", fsimono.get<double>("BASETOL"));  // relative tolerance
  lsParams.set<double>(
      "adaptive distance", fsimono.get<double>("ADAPTIVEDIST"));  // adaptive distance
}

/*----------------------------------------------------------------------*/
/*  print Newton-Raphson iteration to screen and error file             */
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::print_newton_iter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if (Comm().MyPID() == 0)
  {
    if (iter_ == 1) print_newton_iter_header();
    print_newton_iter_text();
  }
}
/*----------------------------------------------------------------------*/
/* print Newton-Raphson iteration to screen and error file              */
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::print_newton_iter_header()
{
  Core::IO::cout << "CONVTOL: " << tolfres_ << Core::IO::endl;

  // open outstringstream
  // std::ostringstream oss;

  Core::IO::cout
      << "===================================================================================="
         "======="
         "========================================================================="
      << Core::IO::endl;

  // enter converged state etc
  Core::IO::cout << "|nit|";

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case Inpar::FSI::convnorm_abs:
      Core::IO::cout << "            "
                     << "abs-res-norm  |";
      break;
    case Inpar::FSI::convnorm_rel:
      Core::IO::cout << "str-rs-l2|"
                     << "fsi-rs-l2|"
                     << "flv-rs-l2|"
                     << "flp-rs-l2|";
      Core::IO::cout << "str-rs-li|"
                     << "fsi-rs-li|"
                     << "flv-rs-li|"
                     << "flp-rs-li|";
      break;
    case Inpar::FSI::convnorm_mix:
      FOUR_C_THROW("not implemented");
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::FSI::convnorm_abs:
      Core::IO::cout << "                  "
                     << "abs-inc-norm";
      break;
    case Inpar::FSI::convnorm_rel:
      Core::IO::cout << "str-in-l2|"
                     << "fsi-in-l2|"
                     << "flv-in-l2|"
                     << "flp-in-l2|";
      Core::IO::cout << "str-in-li|"
                     << "fsi-in-li|"
                     << "flv-in-li|"
                     << "flp-in-li|";
      break;
    case Inpar::FSI::convnorm_mix:
      FOUR_C_THROW("convnorm_mix not implemented");
      break;
    default:
      FOUR_C_THROW("Unknown convergence norm.");
      break;
  }

  // add solution time
  Core::IO::cout << Core::IO::endl;
  Core::IO::cout
      << "===================================================================================="
         "======="
         "========================================================================="
      << Core::IO::endl;
}

/*---------------------------------------------------------------------*/
/*  print Newton-Raphson iteration to screen                           */
/*---------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::print_newton_iter_text()
{
  // enter converged state etc
  Core::IO::cout << " " << iter_ << "/" << itermax_;

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case Inpar::FSI::convnorm_abs:
      Core::IO::cout << "             " << (normrhs_) << Core::IO::endl;
      break;
    case Inpar::FSI::convnorm_rel:
      Core::IO::cout << "|" << (normstrrhsL2_ / ns_) << "|" << (norminterfacerhsL2_ / ni_) << "|"
                     << (normflvelrhsL2_ / nfv_) << "|" << (normflpresrhsL2_ / nfp_) << "|"
                     << (normstrrhsInf_) << "|" << (norminterfacerhsInf_) << "|"
                     << (normflvelrhsInf_) << "|" << (normflpresrhsInf_);
      break;
    case Inpar::FSI::convnorm_mix:
      FOUR_C_THROW("Mixed absolute-relative residual norm not implemented for XFFSI.");
      break;
    default:
      FOUR_C_THROW("Unknown type of residual norm.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::FSI::convnorm_abs:
      Core::IO::cout << "             " << (norminc_) << Core::IO::endl;
      break;
    case Inpar::FSI::convnorm_rel:
      Core::IO::cout << "|" << (normstrincL2_ / ns_) << "|" << (norminterfaceincL2_ / ni_) << "|"
                     << (normflvelincL2_ / nfv_) << "|" << (normflpresincL2_ / nfp_) << "|"
                     << (normstrincInf_) << "|" << (norminterfaceincInf_) << "|"
                     << (normflvelincInf_) << "|" << (normflpresincInf_) << "|" << Core::IO::endl;
      break;
    case Inpar::FSI::convnorm_mix:
      FOUR_C_THROW("Mixed absolute-relative increment norm not implemented for XFFSI.");
      break;
    default:
      FOUR_C_THROW("Unknown type of increment norm.");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::update()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicNoNOX::Update");

  recover_lagrange_multiplier();

  // In case of ALE relaxation
  if (fluid_->monolithic_xffsi_approach() != Inpar::XFEM::XFFSI_Full_Newton and
      fluid_->IsAleRelaxationStep(Step()))
  {
    if (Comm().MyPID() == 0) Core::IO::cout << "Relaxing ALE!" << Core::IO::endl;
    // Set the ALE FSI-DOFs to Dirichlet and solve ALE system again
    // to obtain the true ALE displacement
    ale_field()->Solve();
    // Now apply the ALE-displacement to the (embedded) fluid and update the
    // grid velocity
    fluid_field()->apply_mesh_displacement(ale_to_fluid(ale_field()->Dispnp()));
  }

  // update subsequent fields
  structure_field()->Update();
  fluid_field()->Update();
  ale_field()->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicNoNOX::prepare_time_step()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicNoNOX::prepare_time_step");

  increment_time_and_step();
  print_header();

  structure_field()->prepare_time_step();
  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();

  // no ALE-relaxation or still at the first step? leave!
  if (fluid_->monolithic_xffsi_approach() == Inpar::XFEM::XFFSI_Full_Newton || Step() == 0 ||
      !fluid_->IsAleRelaxationStep(Step() - 1))
    return;

  // recreate the combined dof-map and create a new block system matrix
  // as we have to deal with a new map extrator
  create_combined_dof_row_map();
  systemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          extractor(), extractor(), 81, false, true));
}

FOUR_C_NAMESPACE_CLOSE
