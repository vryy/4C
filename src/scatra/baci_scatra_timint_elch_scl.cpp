/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for elch

\level 2

 *------------------------------------------------------------------------------------------------*/
#include "baci_scatra_timint_elch_scl.H"

#include "baci_adapter_scatra_base_algorithm.H"
#include "baci_coupling_adapter.H"
#include "baci_coupling_adapter_converter.H"
#include "baci_global_data.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_lib_dofset_predefineddofnumber.H"
#include "baci_lib_utils_gid_vector.H"
#include "baci_lib_utils_parallel.H"
#include "baci_lib_utils_parameter_list.H"
#include "baci_linalg_equilibrate.H"
#include "baci_linalg_matrixtransform.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_scatra_ele_action.H"
#include "baci_scatra_resulttest_elch.H"
#include "baci_scatra_timint_elch_service.H"
#include "baci_scatra_timint_meshtying_strategy_s2i_elch.H"
#include "baci_utils_function_of_time.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElchSCL::ScaTraTimIntElchSCL(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntElch(dis, solver, params, sctratimintparams, extraparams, output),
      matrixtype_elch_scl_(
          Teuchos::getIntegralValue<CORE::LINALG::MatrixType>(params->sublist("SCL"), "MATRIXTYPE"))
{
  if (matrixtype_elch_scl_ != CORE::LINALG::MatrixType::sparse and
      matrixtype_elch_scl_ != CORE::LINALG::MatrixType::block_field)
    dserror("Only sparse and block field matrices supported in SCL computations");

  if (INPUT::IntegralValue<int>(*elchparams_, "INITPOTCALC"))
  {
    dserror(
        "Must disable INITPOTCALC for a coupled SCL problem. Use INITPOTCALC in the SCL section "
        "instead.");
  }
  if (!INPUT::IntegralValue<bool>(*params_, "SKIPINITDER"))
  {
    dserror(
        "Must enable SKIPINITDER. Currently, Neumann BCs are not supported in the SCL formulation "
        "and thus, the calculation of the initial time derivative is meaningless.");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("SCL: setup");

  SCATRA::ScaTraTimIntElch::Setup();

  auto* problem = DRT::Problem::Instance();

  auto sdyn_micro =
      Teuchos::rcp(new Teuchos::ParameterList(problem->ScalarTransportDynamicParams()));

  std::string initial_field_type;
  switch (INPUT::IntegralValue<INPAR::SCATRA::InitialField>(
      elchparams_->sublist("SCL"), "INITIALFIELD"))
  {
    case INPAR::SCATRA::initfield_zero_field:
      initial_field_type = "zero_field";
      break;
    case INPAR::SCATRA::initfield_field_by_function:
      initial_field_type = "field_by_function";
      break;
    case INPAR::SCATRA::initfield_field_by_condition:
      initial_field_type = "field_by_condition";
      break;
    default:
      dserror("input type not supported");
      break;
  }
  sdyn_micro->set("INITIALFIELD", initial_field_type);
  sdyn_micro->set("INITFUNCNO", elchparams_->sublist("SCL").get<int>("INITFUNCNO"));

  micro_timint_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(*sdyn_micro, *sdyn_micro,
      problem->SolverParams(sdyn_micro->get<int>("LINEAR_SOLVER")), "scatra_micro", false));

  micro_timint_->Init();

  auto dofset_vel = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(3, 0, 0, true));
  if (micro_timint_->ScaTraField()->Discretization()->AddDofSet(dofset_vel) != 1)
    dserror("unexpected number of dofsets in the scatra micro discretization");
  MicroScaTraField()->SetNumberOfDofSetVelocity(1);

  MicroScaTraField()->Discretization()->FillComplete();

  RedistributeMicroDiscretization();

  MicroScaTraField()->SetVelocityField();

  micro_timint_->Setup();

  // setup coupling between macro and micro field
  SetupCoupling();

  // setup maps for coupled problem
  full_map_elch_scl_ = CORE::LINALG::MergeMap(DofRowMap(), MicroScaTraField()->DofRowMap());
  std::vector<Teuchos::RCP<const Epetra_Map>> block_map_vec_scl;
  switch (matrixtype_elch_scl_)
  {
    case CORE::LINALG::MatrixType::sparse:
      block_map_vec_scl = {full_map_elch_scl_};
      break;
    case CORE::LINALG::MatrixType::block_field:
      block_map_vec_scl = {DofRowMap(), MicroScaTraField()->DofRowMap()};
      break;
    default:
      dserror("Matrix type not supported.");
      break;
  }
  full_block_map_elch_scl_ =
      Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(*full_map_elch_scl_, block_map_vec_scl));

  // setup matrix, rhs, and increment for coupled problem
  increment_elch_scl_ = CORE::LINALG::CreateVector(*full_map_elch_scl_, true);
  residual_elch_scl_ = CORE::LINALG::CreateVector(*full_map_elch_scl_, true);


  switch (matrixtype_elch_scl_)
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      const int expected_entries_per_row = 27;
      const bool explicitdirichlet = false;
      const bool savegraph = true;
      system_matrix_elch_scl_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *full_map_elch_scl_, expected_entries_per_row, explicitdirichlet, savegraph));
      break;
    }
    case CORE::LINALG::MatrixType::block_field:
    {
      const int expected_entries_per_row = 81;
      const bool explicitdirichlet = false;
      const bool savegraph = true;
      system_matrix_elch_scl_ = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *full_block_map_elch_scl_, *full_block_map_elch_scl_, expected_entries_per_row,
              explicitdirichlet, savegraph));
      break;
    }
    default:
      dserror("Matrix type not supported.");
      break;
  }

  // extractor to get micro or macro dofs from global vector
  macro_micro_dofs_ = Teuchos::rcp(
      new CORE::LINALG::MapExtractor(*full_map_elch_scl_, MicroScaTraField()->DofRowMap()));

  dbcmaps_elch_scl_ =
      Teuchos::rcp(new CORE::LINALG::MapExtractor(*full_map_elch_scl_, dbcmaps_->CondMap()));

  // setup solver for coupled problem
  solver_elch_scl_ = Teuchos::rcp(new CORE::LINALG::Solver(
      problem->SolverParams(elchparams_->sublist("SCL").get<int>("SOLVER")), discret_->Comm()));

  switch (matrixtype_elch_scl_)
  {
    case CORE::LINALG::MatrixType::sparse:
      break;
    case CORE::LINALG::MatrixType::block_field:
    {
      std::ostringstream scatrablockstr;
      scatrablockstr << 1;
      Teuchos::ParameterList& blocksmootherparamsscatra =
          solver_elch_scl_->Params().sublist("Inverse" + scatrablockstr.str());
      blocksmootherparamsscatra.sublist("Belos Parameters");
      blocksmootherparamsscatra.sublist("MueLu Parameters");

      Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsscatra);

      std::ostringstream microblockstr;
      microblockstr << 2;
      Teuchos::ParameterList& blocksmootherparamsmicro =
          solver_elch_scl_->Params().sublist("Inverse" + microblockstr.str());
      blocksmootherparamsmicro.sublist("Belos Parameters");
      blocksmootherparamsmicro.sublist("MueLu Parameters");
      MicroScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparamsmicro);

      break;
    }
    default:
      dserror("not supported");
      break;
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::PrepareTimeStep()
{
  if (elchparams_->sublist("SCL").get<int>("ADAPT_TIME_STEP") == Step() + 1)
  {
    const double new_dt = elchparams_->sublist("SCL").get<double>("ADAPTED_TIME_STEP_SIZE");
    if (new_dt <= 0) dserror("new time step size for SCL must be positive.");

    SetDt(new_dt);
    SetTimeStep(Time(), Step());

    MicroScaTraField()->SetDt(new_dt);
    MicroScaTraField()->SetTimeStep(Time(), Step());
    if (discret_->Comm().MyPID() == 0)
      std::cout << "Time step size changed to " << new_dt << std::endl;
  }

  ScaTraTimIntElch::PrepareTimeStep();

  CopySolutionToMicroField();
  MicroScaTraField()->PrepareTimeStep();
  CopySolutionToMicroField();
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::Update()
{
  ScaTraTimIntElch::Update();

  MicroScaTraField()->Update();
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::CheckAndWriteOutputAndRestart()
{
  ScaTraTimIntElch::CheckAndWriteOutputAndRestart();

  MicroScaTraField()->CheckAndWriteOutputAndRestart();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::NonlinearSolve()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:   + nonlin. iteration/lin. solve");

  // out to screen
  PrintTimeStepInfo();

  // prepare Newton-Raphson iteration
  iternum_ = 0;

  CopySolutionToMicroField();

  auto equilibration_method =
      std::vector<CORE::LINALG::EquilibrationMethod>(1, EquilibrationMethod());
  auto equilibration = CORE::LINALG::BuildEquilibration(
      matrixtype_elch_scl_, equilibration_method, full_map_elch_scl_);

  // start Newton-Raphson iteration
  while (true)
  {
    iternum_++;

    // prepare load vector
    neumann_loads_->PutScalar(0.0);

    {
      TEUCHOS_FUNC_TIME_MONITOR("SCL: evaluate");
      // assemble sub problems
      AssembleMatAndRHS();
      MicroScaTraField()->AssembleMatAndRHS();

      // scale micro problem to account for related macro area
      ScaleMicroProblem();

      // couple micro and macro field my nodal mesh tying
      AssembleAndApplyMeshTying();

      system_matrix_elch_scl_->Complete();

      // All DBCs are on the macro scale
      CORE::LINALG::ApplyDirichletToSystem(*system_matrix_elch_scl_, *increment_elch_scl_,
          *residual_elch_scl_, *zeros_, *dbcmaps_elch_scl_->CondMap());

      if (BreakNewtonLoopAndPrintConvergence()) break;
    }

    increment_elch_scl_->PutScalar(0.0);

    {
      TEUCHOS_FUNC_TIME_MONITOR("SCL: solve");

      equilibration->EquilibrateSystem(
          system_matrix_elch_scl_, residual_elch_scl_, full_block_map_elch_scl_);
      CORE::LINALG::SolverParams solver_params;
      solver_params.refactor = true;
      solver_params.reset = iternum_ == 1;
      solver_elch_scl_->Solve(system_matrix_elch_scl_->EpetraOperator(), increment_elch_scl_,
          residual_elch_scl_, solver_params);
      equilibration->UnequilibrateIncrement(increment_elch_scl_);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("SCL: update");

      UpdateIterMicroMacro();

      //-------- update values at intermediate time steps (only for gen.-alpha)
      ComputeIntermediateValues();
      MicroScaTraField()->ComputeIntermediateValues();
      // compute values at the interior of the elements (required for hdg)
      ComputeInteriorValues();
      MicroScaTraField()->ComputeInteriorValues();

      ComputeTimeDerivative();
      MicroScaTraField()->ComputeTimeDerivative();
    }

  }  // nonlinear iteration
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params)
{
  SCATRA::ScaTraTimIntElch::AddProblemSpecificParametersAndVectors(params);

  discret_->SetState("phinp", Phinp());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::CopySolutionToMicroField()
{
  // extract coupled values from macro, copy to micro, and insert into full micro vector
  auto macro_to_micro_coupled_nodes = macro_micro_coupling_adapter_->MasterToSlave(
      macro_coupling_dofs_->ExtractCondVector(Phinp()));
  micro_coupling_dofs_->InsertCondVector(macro_to_micro_coupled_nodes, MicroScaTraField()->Phinp());
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::CreateMeshtyingStrategy()
{
  strategy_ = Teuchos::rcp(new MeshtyingStrategyS2IElchSCL(this, *params_));
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::ReadRestartProblemSpecific(
    int step, IO::DiscretizationReader& reader)
{
  dserror("Restart is not implemented for Elch with SCL.");
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
Teuchos::RCP<SCATRA::ScaTraTimIntImpl> SCATRA::ScaTraTimIntElchSCL::MicroScaTraField()
{
  return micro_timint_->ScaTraField();
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::WriteCouplingToCSV(
    const std::map<int, int>& glob_micro_macro_coupled_node_gids,
    const std::map<int, int>& glob_macro_slave_node_master_node_gids)
{
  std::ofstream file;

  // write GID of coupled nodes to .csv file
  if (myrank_ == 0)
  {
    const std::string file_name_coupling =
        problem_->OutputControlFile()->FileName() + "_micro_macro_coupling.csv";

    file.open(file_name_coupling, std::fstream::trunc);
    file << "macro_slave_node_gid,macro_master_node_gid,micro_slave_node_gid,micro_master_"
            "node_"
            "gid\n";
    file.flush();
    file.close();

    for (const auto& glob_macro_slave_node_master_node_gid : glob_macro_slave_node_master_node_gids)
    {
      const int macro_slave_node_gid = glob_macro_slave_node_master_node_gid.first;
      const int macro_master_node_gid = glob_macro_slave_node_master_node_gid.second;
      const int micro_slave_node_gid =
          glob_micro_macro_coupled_node_gids.find(macro_slave_node_gid)->second;
      const int micro_master_node_gid =
          glob_micro_macro_coupled_node_gids.find(macro_master_node_gid)->second;

      file.open(file_name_coupling, std::fstream::app);
      file << macro_slave_node_gid << "," << macro_master_node_gid << "," << micro_slave_node_gid
           << "," << micro_master_node_gid << "\n";
      file.flush();
      file.close();
    }
  }

  // write node coordinates to .csv file
  const std::string file_name_coords =
      problem_->OutputControlFile()->FileName() + "_micro_macro_coupling_coords.csv";

  if (myrank_ == 0)
  {
    file.open(file_name_coords, std::fstream::trunc);
    file << "node_GID,x,y,z \n";
    file.flush();
    file.close();
  }

  // node coordinates only known by owning proc. Writing of data to file not possible by
  // multiple procs in parallel
  for (int iproc = 0; iproc < discret_->Comm().NumProc(); ++iproc)
  {
    if (iproc == myrank_)
    {
      for (const auto& glob_micro_macro_coupled_node_gid : glob_micro_macro_coupled_node_gids)
      {
        const int macro_node_gid = glob_micro_macro_coupled_node_gid.first;
        const int mirco_node_gid = glob_micro_macro_coupled_node_gid.second;

        if (DRT::UTILS::IsNodeGIDOnThisProc(*discret_, macro_node_gid))
        {
          const auto& macro_coords = discret_->gNode(macro_node_gid)->X();

          file.open(file_name_coords, std::fstream::app);
          file << std::setprecision(16) << std::scientific;
          file << macro_node_gid << "," << macro_coords[0] << "," << macro_coords[1] << ","
               << macro_coords[2] << "\n";
          file.flush();
          file.close();
        }

        if (DRT::UTILS::IsNodeGIDOnThisProc(*MicroScaTraField()->Discretization(), mirco_node_gid))
        {
          const auto& micro_coords =
              MicroScaTraField()->Discretization()->gNode(mirco_node_gid)->X();

          file.open(file_name_coords, std::fstream::app);
          file << std::setprecision(16) << std::scientific;
          file << mirco_node_gid << "," << micro_coords[0] << "," << micro_coords[1] << ","
               << micro_coords[2] << "\n";
          file.flush();
          file.close();
        }
      }
    }
    discret_->Comm().Barrier();
  }
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntElchSCL::BreakNewtonLoopAndPrintConvergence()
{
  // extract processor ID
  const int mypid = discret_->Comm().MyPID();

  const auto& params =
      DRT::Problem::Instance()->ScalarTransportDynamicParams().sublist("NONLINEAR");

  const int itermax = params.get<int>("ITEMAX");
  const double itertol = params.get<double>("CONVTOL");

  auto micro_residual = macro_micro_dofs_->ExtractCondVector(residual_elch_scl_);
  auto macro_residual = macro_micro_dofs_->ExtractOtherVector(residual_elch_scl_);
  auto micro_increment = macro_micro_dofs_->ExtractCondVector(increment_elch_scl_);
  auto macro_increment = macro_micro_dofs_->ExtractOtherVector(increment_elch_scl_);

  double residual_L2, micro_residual_L2, macro_residual_L2, increment_L2, micro_increment_L2,
      macro_increment_L2, micro_state_L2, macro_state_L2;
  residual_elch_scl_->Norm2(&residual_L2);
  micro_residual->Norm2(&micro_residual_L2);
  macro_residual->Norm2(&macro_residual_L2);
  increment_elch_scl_->Norm2(&increment_L2);
  micro_increment->Norm2(&micro_increment_L2);
  macro_increment->Norm2(&macro_increment_L2);
  MicroScaTraField()->Phinp()->Norm2(&micro_state_L2);
  Phinp()->Norm2(&macro_state_L2);

  // safety checks
  if (std::isnan(residual_L2) or std::isnan(micro_residual_L2) or std::isnan(macro_residual_L2) or
      std::isnan(increment_L2) or std::isnan(micro_increment_L2) or
      std::isnan(macro_increment_L2) or std::isnan(micro_state_L2) or std::isnan(macro_state_L2))
    dserror("Calculated vector norm is not a number!");
  if (std::isinf(residual_L2) or std::isinf(micro_residual_L2) or std::isinf(macro_residual_L2) or
      std::isinf(increment_L2) or std::isinf(micro_increment_L2) or
      std::isinf(macro_increment_L2) or std::isinf(micro_state_L2) or std::isinf(macro_state_L2))
    dserror("Calculated vector norm is infinity!");

  micro_state_L2 = micro_state_L2 < 1.0e-10 ? 1.0 : micro_state_L2;
  macro_state_L2 = macro_state_L2 < 1.0e-10 ? 1.0 : micro_state_L2;

  const double state_L2 = std::sqrt(std::pow(micro_state_L2, 2) + std::pow(macro_state_L2, 2));

  micro_increment_L2 = micro_increment_L2 / micro_state_L2;
  macro_increment_L2 = macro_increment_L2 / macro_state_L2;
  increment_L2 = increment_L2 / state_L2;

  const bool finished =
      (residual_L2 < itertol and micro_residual_L2 < itertol and macro_residual_L2 < itertol and
          increment_L2 < itertol and micro_increment_L2 < itertol and
          macro_increment_L2 < itertol and iternum_ > 1) or
      iternum_ == itermax;

  // special case: very first iteration step --> solution increment is not yet available
  if (mypid == 0)
  {
    if (iternum_ == 1)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+-------------+-------------+-------------+---"
                   "----------+-------------+-------------+"
                << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|---  res  ---|---  inc  ---|- micro-res -|- "
                   "micro-inc -|- macro-res -|- macro-inc -| "
                << std::endl;

      // print first line of convergence table to screen
      std::cout << "|  " << std::setw(3) << iternum_ << "/" << std::setw(3) << itermax << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol
                << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << residual_L2 << "  |     --      | " << std::setw(10) << std::setprecision(3)
                << std::scientific << micro_residual_L2 << "  |     --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << macro_residual_L2
                << "  |     --      |  " << std::endl;
    }
    else
    {
      // print current line of convergence table to screen
      std::cout << "|  " << std::setw(3) << iternum_ << "/" << std::setw(3) << itermax << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol
                << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << residual_L2 << "  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << increment_L2 << "  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << micro_residual_L2 << "  | " << std::setw(10)
                << std::setprecision(3) << std::scientific << micro_increment_L2 << "  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << macro_residual_L2
                << "  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << macro_increment_L2 << "  | " << std::endl;


      // convergence check
      if (finished)
      {
        // print finish line of convergence table to screen
        std::cout
            << "+------------+-------------------+-------------+-------------+-------------+---"
               "----------+-------------+-------------+"
            << std::endl;
        if (iternum_ == itermax)
        {
          std::cout << "|      >> Newton-Raphson iteration did not converge! <<                    "
                       "                                          |"
                    << std::endl;
          std::cout
              << "+------------+-------------------+-------------+-------------+-------------+---"
                 "----------+-------------+-------------+"
              << std::endl;
        }
      }
    }
  }
  return finished;
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::SetupCoupling()
{
  TEUCHOS_FUNC_TIME_MONITOR("SCL: setup");

  auto microdis = MicroScaTraField()->Discretization();
  const auto& comm = microdis->Comm();

  // get coupling conditions
  std::vector<DRT::Condition*> macro_coupling_conditions;
  Discretization()->GetCondition("S2ISCLCoupling", macro_coupling_conditions);

  // get all slave and master nodes on this proc from macro coupling condition
  std::vector<int> my_macro_slave_node_gids;
  std::vector<int> my_macro_master_node_gids;
  for (auto* coupling_condition : macro_coupling_conditions)
  {
    for (const int coupling_node_gid : *coupling_condition->Nodes())
    {
      // is this node owned by this proc?
      if (!DRT::UTILS::IsNodeGIDOnThisProc(*discret_, coupling_node_gid)) continue;

      switch (coupling_condition->GetInt("interface side"))
      {
        case INPAR::S2I::side_slave:
          my_macro_slave_node_gids.emplace_back(coupling_node_gid);
          break;
        case INPAR::S2I::side_master:
          my_macro_master_node_gids.emplace_back(coupling_node_gid);
          break;
        default:
          dserror("must be master or slave side");
          break;
      }
    }
  }

  // get master dof(!!) (any proc) to slave nodes(!!) (this proc) from macro coupling adapter
  auto macro_coupling_adapter =
      Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(Strategy())->CouplingAdapter();

  std::map<int, int> my_macro_slave_node_master_dof_gids;
  for (auto my_macro_slave_node_gid : my_macro_slave_node_gids)
  {
    auto* macro_slave_node = discret_->gNode(my_macro_slave_node_gid);
    auto fist_macro_slave_dof_gid = discret_->Dof(0, macro_slave_node)[0];

    for (int slave_dof_lid = 0;
         slave_dof_lid < macro_coupling_adapter->SlaveDofMap()->NumMyElements(); ++slave_dof_lid)
    {
      const int slave_dof_gid = macro_coupling_adapter->SlaveDofMap()->GID(slave_dof_lid);
      if (fist_macro_slave_dof_gid == slave_dof_gid)
      {
        const int first_macro_master_dof_gid =
            macro_coupling_adapter->PermMasterDofMap()->GID(slave_dof_lid);
        my_macro_slave_node_master_dof_gids.insert(
            std::make_pair(my_macro_slave_node_gid, first_macro_master_dof_gid));
        break;
      }
    }
  }
  // distribute all maps to all procs
  const auto glob_macro_slave_node_master_dof_gids =
      DRT::UTILS::BroadcastMap(my_macro_slave_node_master_dof_gids, comm);

  // get master node (this proc) to slave node (any proc)
  std::map<int, int> my_macro_slave_node_master_node_gids;
  for (const auto& glob_macro_slave_node_master_dof_gid : glob_macro_slave_node_master_dof_gids)
  {
    const int master_dof_gid = glob_macro_slave_node_master_dof_gid.second;
    const int slave_node_gid = glob_macro_slave_node_master_dof_gid.first;
    if (DofRowMap()->LID(master_dof_gid) == -1) continue;

    for (const auto my_macro_master_node_gid : my_macro_master_node_gids)
    {
      auto* macro_master_node = discret_->gNode(my_macro_master_node_gid);
      if (discret_->Dof(0, macro_master_node)[0] == master_dof_gid)
      {
        my_macro_slave_node_master_node_gids.insert(
            std::make_pair(slave_node_gid, my_macro_master_node_gid));
        break;
      }
    }
  }
  // distribute all maps to all procs
  const auto glob_macro_slave_node_master_node_gids =
      DRT::UTILS::BroadcastMap(my_macro_slave_node_master_node_gids, comm);

  // we use Dirchlet conditions on micro side to achieve coupling by adapting the DBC value
  std::vector<DRT::Condition*> micro_coupling_conditions;
  microdis->GetCondition("Dirichlet", micro_coupling_conditions);

  if (micro_coupling_conditions.size() != 2) dserror("only 2 DBCs allowed on micro dis");
  if (micro_coupling_conditions[0]->Nodes()->size() !=
      micro_coupling_conditions[1]->Nodes()->size())
    dserror("Number of nodes in micro DBCs are not equal");

  // get all micro coupling nodes
  std::vector<int> my_micro_node_gids;
  for (auto* micro_coupling_condition : micro_coupling_conditions)
  {
    for (const int micro_node_gid : *micro_coupling_condition->Nodes())
    {
      // is this node owned by this proc?
      if (DRT::UTILS::IsNodeGIDOnThisProc(*microdis, micro_node_gid))
        my_micro_node_gids.emplace_back(micro_node_gid);
    }
  }

  // setup coupling between macro and micro problems: find micro problems for this proc (end of last
  // proc)
  int micro_problem_counter = 0;
  int my_micro_problem_counter = 0;
  const unsigned int num_my_macro_slave_node_gids = my_macro_slave_node_gids.size();
  for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
  {
    if (iproc == comm.MyPID())
      micro_problem_counter += static_cast<int>(num_my_macro_slave_node_gids);
    comm.Broadcast(&micro_problem_counter, 1, iproc);

    // start of micro discretization of this proc is end of last proc
    if (iproc == comm.MyPID() - 1) my_micro_problem_counter = micro_problem_counter;
  }

  // global  map between coupled macro nodes and micro nodes
  std::map<int, int> my_macro_micro_coupled_node_gids;
  for (unsigned i = 0; i < num_my_macro_slave_node_gids; ++i)
  {
    const int macro_slave_gid = my_macro_slave_node_gids[i];
    const int macro_master_gid = glob_macro_slave_node_master_node_gids.at(macro_slave_gid);
    const int micro_slave_gid = micro_coupling_conditions[0]->Nodes()->at(my_micro_problem_counter);
    const int micro_master_gid =
        micro_coupling_conditions[1]->Nodes()->at(my_micro_problem_counter);

    my_macro_micro_coupled_node_gids.insert(std::make_pair(macro_slave_gid, micro_slave_gid));
    my_macro_micro_coupled_node_gids.insert(std::make_pair(macro_master_gid, micro_master_gid));
    my_micro_problem_counter++;
  }
  const auto glob_macro_micro_coupled_node_gids =
      DRT::UTILS::BroadcastMap(my_macro_micro_coupled_node_gids, comm);

  // setup macro nodes on this proc and coupled micro nodes (can be other proc)
  std::vector<int> my_micro_permuted_node_gids;
  std::vector<int> my_macro_node_gids;
  for (const auto& glob_macro_micro_coupled_node_gid : glob_macro_micro_coupled_node_gids)
  {
    const int macro_node_gid = glob_macro_micro_coupled_node_gid.first;
    const int mirco_node_gid = glob_macro_micro_coupled_node_gid.second;

    if (!DRT::UTILS::IsNodeGIDOnThisProc(*discret_, macro_node_gid)) continue;

    my_macro_node_gids.emplace_back(macro_node_gid);
    my_micro_permuted_node_gids.emplace_back(mirco_node_gid);
  }

  if (INPUT::IntegralValue<bool>(elchparams_->sublist("SCL"), "COUPLING_OUTPUT"))
    WriteCouplingToCSV(glob_macro_micro_coupled_node_gids, glob_macro_slave_node_master_node_gids);

  // setup Epetra maps for coupled nodes
  auto master_node_map = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, static_cast<int>(my_macro_node_gids.size()), &my_macro_node_gids[0], 0, comm));
  auto slave_node_map = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, static_cast<int>(my_micro_node_gids.size()), &my_micro_node_gids[0], 0, comm));
  auto perm_slave_node_map = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, static_cast<int>(my_micro_permuted_node_gids.size()),
          &my_micro_permuted_node_gids[0], 0, comm));

  // setup coupling adapter between micro (slave) and macro (master) for all dof of the nodes
  auto macro_micro_coupling_adapter_temp = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  macro_micro_coupling_adapter_temp->SetupCoupling(*discret_, *microdis, *master_node_map,
      *slave_node_map, *perm_slave_node_map, NumDofPerNode());

  // setup actual coupling adapter only for dofs for which coupling is activated
  std::vector<int> my_slave_dofs;
  std::vector<int> my_perm_master_dofs;
  for (int slave_lid = 0;
       slave_lid < macro_micro_coupling_adapter_temp->SlaveDofMap()->NumMyElements(); ++slave_lid)
  {
    const int slave_gid = macro_micro_coupling_adapter_temp->SlaveDofMap()->GID(slave_lid);

    for (int dbc_lid = 0; dbc_lid < MicroScaTraField()->DirichMaps()->CondMap()->NumMyElements();
         ++dbc_lid)
    {
      const int dbc_gid = MicroScaTraField()->DirichMaps()->CondMap()->GID(dbc_lid);
      if (slave_gid == dbc_gid)
      {
        my_slave_dofs.emplace_back(slave_gid);
        my_perm_master_dofs.emplace_back(
            macro_micro_coupling_adapter_temp->PermMasterDofMap()->GID(slave_lid));
        break;
      }
    }
  }

  const std::vector<int> glob_slave_dofs = DRT::UTILS::BroadcastVector(my_slave_dofs, comm);

  std::vector<int> my_master_dofs;
  std::vector<int> my_perm_slave_dofs;
  for (int master_lid = 0;
       master_lid < macro_micro_coupling_adapter_temp->MasterDofMap()->NumMyElements();
       ++master_lid)
  {
    const int slave_gid = macro_micro_coupling_adapter_temp->PermSlaveDofMap()->GID(master_lid);
    const int master_gid = macro_micro_coupling_adapter_temp->MasterDofMap()->GID(master_lid);
    if (std::find(glob_slave_dofs.begin(), glob_slave_dofs.end(), slave_gid) !=
        glob_slave_dofs.end())
    {
      my_master_dofs.emplace_back(master_gid);
      my_perm_slave_dofs.emplace_back(slave_gid);
    }
  }

  auto slave_dof_map = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, static_cast<int>(my_slave_dofs.size()), &my_slave_dofs[0], 0, comm));
  auto perm_slave_dof_map = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, static_cast<int>(my_perm_slave_dofs.size()), &my_perm_slave_dofs[0], 0, comm));
  auto master_dof_map = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, static_cast<int>(my_master_dofs.size()), &my_master_dofs[0], 0, comm));
  auto perm_master_dof_map = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, static_cast<int>(my_perm_master_dofs.size()), &my_perm_master_dofs[0], 0, comm));


  macro_micro_coupling_adapter_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  macro_micro_coupling_adapter_->SetupCoupling(
      slave_dof_map, perm_slave_dof_map, master_dof_map, perm_master_dof_map);

  macro_coupling_dofs_ = Teuchos::rcp(
      new CORE::LINALG::MapExtractor(*DofRowMap(), macro_micro_coupling_adapter_->MasterDofMap()));

  micro_coupling_dofs_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(
      *microdis->DofRowMap(), macro_micro_coupling_adapter_->SlaveDofMap()));

  // setup relation between first node of micro sub problem and following nodes. This is required
  // for scaling (see ScaleMicroProblem())
  std::set<int> my_micro_coupling_nodes;
  for (int lid_micro = 0; lid_micro < macro_micro_coupling_adapter_->SlaveDofMap()->NumMyElements();
       ++lid_micro)
  {
    const int gid_micro = macro_micro_coupling_adapter_->SlaveDofMap()->GID(lid_micro);
    my_micro_coupling_nodes.insert(gid_micro);
  }

  const auto glob_micro_coupling_nodes =
      DRT::UTILS::BroadcastSet(my_micro_coupling_nodes, discret_->Comm());

  // by definition, the last node of a micro sub problem is coupled with the macro. Here, all nodes
  // in the sub problem are linked to the coupled node, by looping backwards through all nodes and
  // forwards through the coupled nodes
  for (int lid_micro = MicroScaTraField()->DofRowMap()->NumMyElements() - 1; lid_micro >= 0;
       --lid_micro)
  {
    const int gid_micro = MicroScaTraField()->DofRowMap()->GID(lid_micro);
    for (const auto& coupled_node : glob_micro_coupling_nodes)
    {
      if (coupled_node >= gid_micro)
      {
        coupled_micro_nodes_.insert(std::make_pair(gid_micro, coupled_node));
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::ScaleMicroProblem()
{
  Teuchos::ParameterList condparams;

  // scale micro problem with nodal area of macro discretiaztion
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_nodal_size, condparams);

  auto nodal_size_macro = CORE::LINALG::CreateVector(*DofRowMap(), true);
  discret_->EvaluateCondition(condparams, Teuchos::null, Teuchos::null, nodal_size_macro,
      Teuchos::null, Teuchos::null, "S2ISCLCoupling");

  // extract dof values to node values
  for (int row_lid = 0; row_lid < DofRowMap()->NumMyElements(); row_lid += 2)
  {
    const double row_value = (*nodal_size_macro)[row_lid];
    const double scale_fac = row_value == 0.0 ? 1.0 : row_value;
    for (int dof = 0; dof < NumDofPerNode(); ++dof) (*nodal_size_macro)[row_lid + dof] = scale_fac;
  }

  // transform to micro discretization
  auto nodal_size_micro = macro_micro_coupling_adapter_->MasterToSlave(
      macro_coupling_dofs_->ExtractCondVector(nodal_size_macro));

  // communicate nodal size to all procs to be able to scale all rows in micro discretization
  // attached to a macro node
  std::map<int, double> my_nodal_size_micro;
  for (int lid_micro = 0; lid_micro < macro_micro_coupling_adapter_->SlaveDofMap()->NumMyElements();
       ++lid_micro)
  {
    const int gid_micro = macro_micro_coupling_adapter_->SlaveDofMap()->GID(lid_micro);
    my_nodal_size_micro.insert(std::make_pair(gid_micro, (*nodal_size_micro)[lid_micro]));
  }

  const auto glob_nodal_size_micro =
      DRT::UTILS::BroadcastMap(my_nodal_size_micro, discret_->Comm());

  auto micro_scale = CORE::LINALG::CreateVector(*MicroScaTraField()->DofRowMap(), true);
  for (int lid_micro = MicroScaTraField()->DofRowMap()->NumMyElements() - 1; lid_micro >= 0;
       --lid_micro)
  {
    const int gid_micro = MicroScaTraField()->DofRowMap()->GID(lid_micro);
    const int coupled_node = coupled_micro_nodes_[gid_micro];
    const double scale_val = glob_nodal_size_micro.at(coupled_node);
    (*micro_scale)[lid_micro] = scale_val;
    (*MicroScaTraField()->Residual())[lid_micro] *= scale_val;
  }
  MicroScaTraField()->SystemMatrix()->LeftScale(*micro_scale);
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::AssembleAndApplyMeshTying()
{
  // Meshtying + Assembly RHS
  auto micro_residual = micro_coupling_dofs_->ExtractCondVector(MicroScaTraField()->Residual());
  auto micro_residual_on_macro_side = macro_micro_coupling_adapter_->SlaveToMaster(micro_residual);

  auto full_macro_vector = CORE::LINALG::CreateVector(*DofRowMap(), true);
  macro_coupling_dofs_->InsertCondVector(micro_residual_on_macro_side, full_macro_vector);

  residual_elch_scl_->PutScalar(0.0);
  system_matrix_elch_scl_->Zero();

  macro_micro_dofs_->AddOtherVector(full_macro_vector, residual_elch_scl_);

  macro_micro_dofs_->AddOtherVector(Residual(), residual_elch_scl_);

  // apply pseudo DBC on slave side
  micro_coupling_dofs_->CondPutScalar(*MicroScaTraField()->Residual(), 0.0);

  macro_micro_dofs_->AddCondVector(
      MicroScaTraField()->Residual(), residual_elch_scl_);  // add slave side to total residual

  switch (matrixtype_elch_scl_)
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      auto sparse_systemmatrix =
          CORE::LINALG::CastToSparseMatrixAndCheckSuccess(system_matrix_elch_scl_);

      sparse_systemmatrix->Add(*SystemMatrix(), false, 1.0, 1.0);

      auto micro_side_converter =
          Teuchos::rcp(new CORE::ADAPTER::CouplingSlaveConverter(*macro_micro_coupling_adapter_));

      // micro: interior - interior
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->OtherMap(), *micro_coupling_dofs_->OtherMap(), 1.0, nullptr,
          nullptr, *sparse_systemmatrix, true, true);

      // micro: interior - slave
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->OtherMap(), *micro_coupling_dofs_->CondMap(), 1.0, nullptr,
          &(*micro_side_converter), *sparse_systemmatrix, true, true);

      // micro: slave - interior
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->CondMap(), *micro_coupling_dofs_->OtherMap(), 1.0,
          &(*micro_side_converter), nullptr, *sparse_systemmatrix, true, true);

      // micro: slave - slave
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->CondMap(), *micro_coupling_dofs_->CondMap(), 1.0,
          &(*micro_side_converter), &(*micro_side_converter), *sparse_systemmatrix, true, true);
      break;
    }
    case CORE::LINALG::MatrixType::block_field:
    {
      auto block_systemmatrix =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(system_matrix_elch_scl_);

      block_systemmatrix->Matrix(0, 0).Add(*SystemMatrix(), false, 1.0, 1.0);

      auto micro_side_converter =
          Teuchos::rcp(new CORE::ADAPTER::CouplingSlaveConverter(*macro_micro_coupling_adapter_));

      // micro: interior - interior
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->OtherMap(), *micro_coupling_dofs_->OtherMap(), 1.0, nullptr,
          nullptr, block_systemmatrix->Matrix(1, 1), true, true);

      // micro: interior - slave
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->OtherMap(), *micro_coupling_dofs_->CondMap(), 1.0, nullptr,
          &(*micro_side_converter), block_systemmatrix->Matrix(1, 0), true, true);

      // micro: slave - interior
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->CondMap(), *micro_coupling_dofs_->OtherMap(), 1.0,
          &(*micro_side_converter), nullptr, block_systemmatrix->Matrix(0, 1), true, true);

      // micro: slave - slave
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*MicroScaTraField()->SystemMatrix(),
          *micro_coupling_dofs_->CondMap(), *micro_coupling_dofs_->CondMap(), 1.0,
          &(*micro_side_converter), &(*micro_side_converter), block_systemmatrix->Matrix(0, 0),
          true, true);
      break;
    }
    default:
      dserror("not supported");
      break;
  }

  // pseudo DBCs on slave side
  CORE::LINALG::SparseMatrix& micromatrix =
      matrixtype_elch_scl_ == CORE::LINALG::MatrixType::sparse
          ? *CORE::LINALG::CastToSparseMatrixAndCheckSuccess(system_matrix_elch_scl_)
          : CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(system_matrix_elch_scl_)
                ->Matrix(1, 1);
  auto slavemaps = macro_micro_coupling_adapter_->SlaveDofMap();
  const double one = 1.0;
  for (int doflid_slave = 0; doflid_slave < slavemaps->NumMyElements(); ++doflid_slave)
  {
    // extract global ID of current slave-side row
    const int dofgid_slave = slavemaps->GID(doflid_slave);
    if (dofgid_slave < 0) dserror("Local ID not found!");

    // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column indices
    if (micromatrix.Filled())
    {
      const int rowlid_slave = micromatrix.RowMap().LID(dofgid_slave);
      if (rowlid_slave < 0) dserror("Global ID not found!");
      if (micromatrix.EpetraMatrix()->ReplaceMyValues(rowlid_slave, 1, &one, &rowlid_slave))
        dserror("ReplaceMyValues failed!");
    }

    // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and column
    // indices
    else
    {
      micromatrix.EpetraMatrix()->InsertGlobalValues(dofgid_slave, 1, &one, &dofgid_slave);
    }
  }
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::UpdateIterMicroMacro()
{
  auto increment_macro = macro_micro_dofs_->ExtractOtherVector(increment_elch_scl_);
  auto increment_micro = macro_micro_dofs_->ExtractCondVector(increment_elch_scl_);

  // reconstruct slave result from master side
  auto macro_extract = macro_coupling_dofs_->ExtractCondVector(increment_macro);
  auto macro_extract_to_micro = macro_micro_coupling_adapter_->MasterToSlave(macro_extract);
  micro_coupling_dofs_->InsertCondVector(macro_extract_to_micro, increment_micro);

  UpdateIter(increment_macro);
  MicroScaTraField()->UpdateIter(increment_micro);
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::RedistributeMicroDiscretization()
{
  auto micro_dis = MicroScaTraField()->Discretization();
  const int min_node_gid = micro_dis->NodeRowMap()->MinAllGID();
  const int num_nodes = micro_dis->NodeRowMap()->NumGlobalElements();
  const int num_proc = micro_dis->Comm().NumProc();
  const int myPID = micro_dis->Comm().MyPID();

  const int num_node_per_proc = static_cast<int>(std::floor(num_nodes / num_proc));

  // new node row list: split node list by number of processors
  std::vector<int> my_row_nodes(num_node_per_proc, -1);
  if (myPID == num_proc - 1) my_row_nodes.resize(num_nodes - (num_proc - 1) * num_node_per_proc);
  std::iota(my_row_nodes.begin(), my_row_nodes.end(), min_node_gid + myPID * num_node_per_proc);

  // new node col list: add boundary nodes of other procs (first and last node of list)
  std::vector<int> my_col_nodes = my_row_nodes;
  if (myPID > 0) my_col_nodes.emplace_back(my_row_nodes[0] - 1);
  if (myPID < num_proc - 1) my_col_nodes.emplace_back(my_row_nodes.back() + 1);

  auto new_node_row_map = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      num_nodes, static_cast<int>(my_row_nodes.size()), &my_row_nodes[0], 0, micro_dis->Comm()));

  auto new_node_col_map = Teuchos::rcp<Epetra_Map>(new Epetra_Map(
      -1, static_cast<int>(my_col_nodes.size()), &my_col_nodes[0], 0, micro_dis->Comm()));

  micro_dis->Redistribute(*new_node_row_map, *new_node_col_map, true, true, true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::PrepareTimeLoop()
{
  // call base class routine
  ScaTraTimIntElch::PrepareTimeLoop();

  if (INPUT::IntegralValue<int>(elchparams_->sublist("SCL"), "INITPOTCALC"))
    CalcInitialPotentialField();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::CalcInitialPotentialField()
{
  PreCalcInitialPotentialField();
  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(MicroScaTraField())
      ->PreCalcInitialPotentialField();

  // safety checks
  dsassert(step_ == 0, "Step counter is not zero!");

  if (equpot_ != INPAR::ELCH::equpot_divi)
  {
    dserror(
        "Initial potential field cannot be computed for chosen closing equation for electric "
        "potential!");
  }

  // screen output
  if (myrank_ == 0)
  {
    std::cout << "SCATRA: calculating initial field for electric potential" << std::endl;
    PrintTimeStepInfo();
    std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
    std::cout << "|- step/max -|- tol      [norm] -|--   res   ---|--   inc   ---|" << std::endl;
  }

  // prepare Newton-Raphson iteration
  iternum_ = 0;
  const int itermax = params_->sublist("NONLINEAR").get<int>("ITEMAX");
  const double itertol = params_->sublist("NONLINEAR").get<double>("CONVTOL");

  CopySolutionToMicroField();

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    iternum_++;

    // prepare load vector
    neumann_loads_->PutScalar(0.0);

    // assemble sub problems
    AssembleMatAndRHS();
    MicroScaTraField()->AssembleMatAndRHS();

    // scale micro problem to account for related macro area
    ScaleMicroProblem();

    // couple micro and macro field my nodal mesh tying
    AssembleAndApplyMeshTying();

    system_matrix_elch_scl_->Complete();

    // All DBCs are on the macro scale
    CORE::LINALG::ApplyDirichletToSystem(*system_matrix_elch_scl_, *increment_elch_scl_,
        *residual_elch_scl_, *zeros_, *dbcmaps_elch_scl_->CondMap());

    // apply artificial Dirichlet boundary conditions to system of equations
    // to hold initial concentrations constant when solving for initial potential field
    auto pseudo_dbc_scl =
        CORE::LINALG::MergeMap(splitter_->OtherMap(), MicroScaTraField()->Splitter()->OtherMap());
    auto pseudo_zeros_scl = CORE::LINALG::CreateVector(*pseudo_dbc_scl, true);

    CORE::LINALG::ApplyDirichletToSystem(*system_matrix_elch_scl_, *increment_elch_scl_,
        *residual_elch_scl_, *pseudo_zeros_scl, *pseudo_dbc_scl);

    // compute L2 norm of state vector
    double state_L2_macro, state_L2_micro;
    Phinp()->Norm2(&state_L2_macro);
    MicroScaTraField()->Phinp()->Norm2(&state_L2_micro);
    double state_L2 = std::sqrt(std::pow(state_L2_macro, 2) + std::pow(state_L2_micro, 2));

    // compute L2 residual vector
    double res_L2, inc_L2;
    residual_elch_scl_->Norm2(&res_L2);

    // compute L2 norm of increment vector
    increment_elch_scl_->Norm2(&inc_L2);

    // safety checks
    if (std::isnan(inc_L2) or std::isnan(res_L2)) dserror("calculated vector norm is NaN.");
    if (std::isinf(inc_L2) or std::isinf(res_L2)) dserror("calculated vector norm is INF.");

    // care for the case that nothing really happens
    if (state_L2 < 1.0e-5) state_L2 = 1.0;

    // first iteration step: solution increment is not yet available
    if (iternum_ == 1)
    {
      // print first line of convergence table to screen
      if (myrank_ == 0)
      {
        std::cout << "|  " << std::setw(3) << iternum_ << "/" << std::setw(3) << itermax << "   | "
                  << std::setw(10) << std::setprecision(3) << std::scientific << itertol
                  << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific << 0.0
                  << "   |      --      | " << std::endl;
      }
    }

    // later iteration steps: solution increment can be printed
    else
    {
      // print current line of convergence table to screen
      if (myrank_ == 0)
      {
        std::cout << "|  " << std::setw(3) << iternum_ << "/" << std::setw(3) << itermax << "   | "
                  << std::setw(10) << std::setprecision(3) << std::scientific << itertol
                  << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                  << res_L2 << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                  << inc_L2 / state_L2 << "   | " << std::endl;
      }

      // convergence check
      if (res_L2 <= itertol and inc_L2 / state_L2 <= itertol)
      {
        // print finish line of convergence table to screen
        if (myrank_ == 0)
        {
          std::cout << "+------------+-------------------+--------------+--------------+"
                    << std::endl
                    << std::endl;
        }

        // abort Newton-Raphson iteration
        break;
      }
    }

    // warn if maximum number of iterations is reached without convergence
    if (iternum_ == itermax)
    {
      if (myrank_ == 0)
      {
        std::cout << "+--------------------------------------------------------------+"
                  << std::endl;
        std::cout << "|            >>>>>> not converged!                             |"
                  << std::endl;
        std::cout << "+--------------------------------------------------------------+" << std::endl
                  << std::endl;
      }

      // abort Newton-Raphson iteration
      break;
    }

    // zero out increment vector
    increment_elch_scl_->PutScalar(0.0);

    CORE::LINALG::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = iternum_ == 1;
    solver_elch_scl_->Solve(system_matrix_elch_scl_->EpetraOperator(), increment_elch_scl_,
        residual_elch_scl_, solver_params);

    UpdateIterMicroMacro();

    // copy initial state vector
    phin_->Update(1., *phinp_, 0.);
    MicroScaTraField()->Phin()->Update(1.0, *MicroScaTraField()->Phinp(), 0.0);

    // update state vectors for intermediate time steps (only for generalized alpha)
    ComputeIntermediateValues();
    MicroScaTraField()->ComputeIntermediateValues();
  }  // Newton-Raphson iteration

  // reset global system matrix and its graph, since we solved a very special problem with a
  // special sparsity pattern
  system_matrix_elch_scl_->Reset();

  PostCalcInitialPotentialField();

  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(MicroScaTraField())
      ->PostCalcInitialPotentialField();
}

Teuchos::RCP<DRT::ResultTest> SCATRA::ScaTraTimIntElchSCL::CreateMicroFieldTest()
{
  return Teuchos::rcp(new SCATRA::ElchResultTest(
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(MicroScaTraField())));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElchSCL::TestResults()
{
  DRT::Problem::Instance()->AddFieldTest(CreateScaTraFieldTest());
  DRT::Problem::Instance()->AddFieldTest(CreateMicroFieldTest());
  DRT::Problem::Instance()->TestAll(discret_->Comm());
}

BACI_NAMESPACE_CLOSE
