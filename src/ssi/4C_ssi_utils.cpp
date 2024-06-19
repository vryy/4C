/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 1


 *------------------------------------------------------------------------------------------------*/

#include "4C_ssi_utils.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_comm_broadcast_utils.hpp"
#include "4C_comm_utils_gid_vector.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_coupling_matchingoctree.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_s2i.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_ssi_monolithic.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int SSI::UTILS::CheckTimeStepping(double dt1, double dt2)
{
  const double workdt1 = std::min(dt1, dt2);
  const double workdt2 = std::max(dt1, dt2);
  int i = 0;

  while (true)
  {
    i++;
    const double t1 = i * workdt1;

    if (std::abs(t1 - workdt2) < 10E-10)
      break;
    else if (t1 > workdt2)
      FOUR_C_THROW("Chosen time steps %f and %f are not a multiplicative of each other", dt1, dt2);
  }
  return i;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 10/2014 */
// Modification of time parameter list for problem with different time step size

void SSI::UTILS::ChangeTimeParameter(const Epetra_Comm& comm, Teuchos::ParameterList& ssiparams,
    Teuchos::ParameterList& scatradyn, Teuchos::ParameterList& sdyn)
{
  bool difftimestep = Core::UTILS::IntegralValue<int>(ssiparams, "DIFFTIMESTEPSIZE");

  if (difftimestep)  // Create subproblems with different time steps
  {
    // Check correct choice of time stepping for single fields
    double scatrastep = scatradyn.get<double>("TIMESTEP");
    double solidstep = sdyn.get<double>("TIMESTEP");

    SSI::UTILS::CheckTimeStepping(scatrastep, solidstep);

    // modify global time step size
    ssiparams.set<double>("TIMESTEP", std::min(scatrastep, solidstep));
  }
  else
  {
    // -------------------------------------------------------------------
    // overrule certain parameters for coupled problems
    // -------------------------------------------------------------------
    // the default time step size
    scatradyn.set<double>("TIMESTEP", ssiparams.get<double>("TIMESTEP"));
    sdyn.set<double>("TIMESTEP", ssiparams.get<double>("TIMESTEP"));
    // maximum simulation time
    scatradyn.set<double>("MAXTIME", ssiparams.get<double>("MAXTIME"));
    sdyn.set<double>("MAXTIME", ssiparams.get<double>("MAXTIME"));
    // maximum number of timesteps
    scatradyn.set<int>("NUMSTEP", ssiparams.get<int>("NUMSTEP"));
    sdyn.set<int>("NUMSTEP", ssiparams.get<int>("NUMSTEP"));
  }

  // Check correct input of restart. Code relies that both time value RESTARTEVRYTIME and
  // RESULTSEVRYTIME are given if restart from time is applied
  double restarttime = ssiparams.get<double>("RESTARTEVRYTIME");
  double updatetime = ssiparams.get<double>("RESULTSEVRYTIME");
  if ((updatetime > 0.0) or (restarttime > 0.0))
  {
    if (updatetime <= 0.0 and restarttime <= 0.0)
    {
      FOUR_C_THROW(
          "If time controlled output and restart is desired, both parameters RESTARTEVRYTIME and "
          "RESULTSEVRYTIME has to be set");
    }
  }

  // set restart params
  int scatrarestart;
  int structurerestart;

  if (restarttime > 0.0)
  {
    scatrarestart = SSI::UTILS::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), restarttime);
    structurerestart = SSI::UTILS::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), restarttime);
  }
  else
  {
    int restart = ssiparams.get<int>("RESTARTEVRY");
    scatrarestart = restart;
    structurerestart = restart;
  }

  // set output params
  int scatraupres;
  int structureupres;

  if (updatetime > 0.0)
  {
    scatraupres = SSI::UTILS::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), updatetime);
    structureupres = SSI::UTILS::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), updatetime);
  }
  else
  {
    int update = ssiparams.get<int>("RESULTSEVRY");
    scatraupres = update;
    structureupres = update;
  }

  // restart
  scatradyn.set<int>("RESTARTEVRY", scatrarestart);
  sdyn.set<int>("RESTARTEVRY", structurerestart);
  // solution output
  scatradyn.set<int>("RESULTSEVRY", scatraupres);
  sdyn.set<int>("RESULTSEVRY", structureupres);

  if (comm.MyPID() == 0)
  {
    std::cout << "====================== Overview of chosen time stepping: "
                 "==============================\n"
              << "\t Timestep scatra:           " << scatradyn.get<double>("TIMESTEP") << "\n"
              << "\t Timestep structure:        " << sdyn.get<double>("TIMESTEP") << "\n"
              << "\t Result step scatra:        " << scatradyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Result step structure:     " << sdyn.get<int>("RESULTSEVRY") << "\n"
              << "\t Restart step scatra:       " << scatradyn.get<int>("RESTARTEVRY") << "\n"
              << "\t Restart step structure:    " << sdyn.get<int>("RESTARTEVRY") << "\n"
              << "================================================================================="
                 "=======\n \n";
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSI::UTILS::CloneScaTraManifoldParams(
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& sublist_manifold_params)
{
  Teuchos::ParameterList scatra_manifold_params(scatraparams);

  switch (Core::UTILS::IntegralValue<Inpar::ScaTra::InitialField>(
      sublist_manifold_params, "INITIALFIELD"))
  {
    case Inpar::ScaTra::initfield_zero_field:
    {
      scatra_manifold_params.set<std::string>("INITIALFIELD", "zero_field");
      scatra_manifold_params.set<int>("INITFUNCNO", -1);
      break;
    }
    case Inpar::ScaTra::initfield_field_by_function:
    {
      scatra_manifold_params.set<std::string>("INITIALFIELD", "field_by_function");
      scatra_manifold_params.set<int>("INITFUNCNO", sublist_manifold_params.get<int>("INITFUNCNO"));
      break;
    }
    case Inpar::ScaTra::initfield_field_by_condition:
    {
      scatra_manifold_params.set<std::string>("INITIALFIELD", "field_by_condition");
      scatra_manifold_params.set<int>("INITFUNCNO", -1);
      break;
    }
    default:
      FOUR_C_THROW("Initial field type on manifold not supported.");
      break;
  }

  if (Core::UTILS::IntegralValue<Inpar::ScaTra::OutputScalarType>(scatraparams, "OUTPUTSCALARS") !=
      Inpar::ScaTra::outputscalars_none)
    scatra_manifold_params.set<bool>("output_file_name_discretization", true);

  scatra_manifold_params.set<std::string>("OUTPUTSCALARSMEANGRAD", "No");

  scatra_manifold_params.set<std::string>("ADAPTIVE_TIMESTEPPING", "No");

  return scatra_manifold_params;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList SSI::UTILS::ModifyScaTraParams(const Teuchos::ParameterList& scatraparams)
{
  auto scatraparams_mutable = Teuchos::ParameterList(scatraparams);

  if (Core::UTILS::IntegralValue<Inpar::ScaTra::OutputScalarType>(scatraparams, "OUTPUTSCALARS") !=
      Inpar::ScaTra::outputscalars_none)
    scatraparams_mutable.set<bool>("output_file_name_discretization", true);

  return scatraparams_mutable;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::UTILS::SSIMatrices::SSIMatrices(Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps,
    const Core::LinAlg::MatrixType ssi_matrixtype, const Core::LinAlg::MatrixType scatra_matrixtype,
    const bool is_scatra_manifold)
    : is_scatra_manifold_(is_scatra_manifold),
      scatra_matrixtype_(scatra_matrixtype),
      scatra_dofrowmap_(ssi_maps->ScaTraDofRowMap()),
      structure_dofrowmap_(ssi_maps->StructureDofRowMap())
{
  // fill maps related to scalar transport manifold if relevant
  if (is_scatra_manifold_) scatramanifold_dofrowmap_ = ssi_maps->sca_tra_manifold_dof_row_map();

  initialize_system_matrix(ssi_maps, ssi_matrixtype);

  initialize_main_diag_matrices(ssi_maps);

  initialize_off_diag_matrices(ssi_maps);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::initialize_main_diag_matrices(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
{
  structure_matrix_ = setup_sparse_matrix(structure_dofrowmap_);

  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      scatra_matrix_ = setup_block_matrix(ssi_maps->BlockMapScaTra(), ssi_maps->BlockMapScaTra());
      if (is_scatra_manifold_)
        manifold_matrix_ = setup_block_matrix(
            ssi_maps->block_map_sca_tra_manifold(), ssi_maps->block_map_sca_tra_manifold());

      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      scatra_matrix_ = setup_sparse_matrix(scatra_dofrowmap_);

      if (is_scatra_manifold_) manifold_matrix_ = setup_sparse_matrix(scatramanifold_dofrowmap_);

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::initialize_off_diag_matrices(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps)
{
  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      scatra_structure_matrix_ =
          setup_block_matrix(ssi_maps->BlockMapScaTra(), ssi_maps->BlockMapStructure());

      structure_scatra_matrix_ =
          setup_block_matrix(ssi_maps->BlockMapStructure(), ssi_maps->BlockMapScaTra());

      if (is_scatra_manifold_)
      {
        scatramanifold_structure_matrix_ = setup_block_matrix(
            ssi_maps->block_map_sca_tra_manifold(), ssi_maps->BlockMapStructure());
        scatramanifold_scatra_matrix_ =
            setup_block_matrix(ssi_maps->block_map_sca_tra_manifold(), ssi_maps->BlockMapScaTra());
        scatra_scatramanifold_matrix_ =
            setup_block_matrix(ssi_maps->BlockMapScaTra(), ssi_maps->block_map_sca_tra_manifold());
      }

      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      scatra_structure_matrix_ = setup_sparse_matrix(scatra_dofrowmap_);
      structure_scatra_matrix_ = setup_sparse_matrix(structure_dofrowmap_);

      if (is_scatra_manifold_)
      {
        scatramanifold_structure_matrix_ = setup_sparse_matrix(scatramanifold_dofrowmap_);
        scatramanifold_scatra_matrix_ = setup_sparse_matrix(scatramanifold_dofrowmap_);
        scatra_scatramanifold_matrix_ = setup_sparse_matrix(scatra_dofrowmap_);
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::complete_sca_tra_manifold_sca_tra_matrix()
{
  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
      sca_tra_manifold_sca_tra_matrix()()->Complete(*scatra_dofrowmap_, *scatramanifold_dofrowmap_);
      break;
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
      sca_tra_manifold_sca_tra_matrix()()->Complete();
      break;
    default:
      FOUR_C_THROW("Not supported Core::LinAlg::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::complete_sca_tra_manifold_structure_matrix()
{
  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
      sca_tra_manifold_structure_matrix()->Complete(
          *structure_dofrowmap_, *scatramanifold_dofrowmap_);
      break;
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
      sca_tra_manifold_structure_matrix()->Complete();
      break;
    default:
      FOUR_C_THROW("Not supported Core::LinAlg::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::complete_sca_tra_sca_tra_manifold_matrix()
{
  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
      sca_tra_sca_tra_manifold_matrix()->Complete(*scatramanifold_dofrowmap_, *scatra_dofrowmap_);
      break;
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
      sca_tra_sca_tra_manifold_matrix()->Complete();
      break;
    default:
      FOUR_C_THROW("Not supported Core::LinAlg::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::complete_sca_tra_structure_matrix()
{
  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
      sca_tra_structure_matrix()->Complete(*structure_dofrowmap_, *scatra_dofrowmap_);
      break;
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
      sca_tra_structure_matrix()->Complete();
      break;
    default:
      FOUR_C_THROW("Not supported Core::LinAlg::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::complete_structure_sca_tra_matrix()
{
  switch (scatra_matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
      structure_sca_tra_matrix()->Complete(*scatra_dofrowmap_, *structure_dofrowmap_);
      break;
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
      structure_sca_tra_matrix()->Complete();
      break;
    default:
      FOUR_C_THROW("Not supported Core::LinAlg::MatrixType!");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::ClearMatrices()
{
  system_matrix_->Zero();
  scatra_matrix_->Zero();
  scatra_structure_matrix_->Zero();
  structure_scatra_matrix_->Zero();
  structure_matrix_->Zero();

  if (is_scatra_manifold_)
  {
    scatramanifold_structure_matrix_->Zero();
    manifold_matrix_->Zero();
    scatra_scatramanifold_matrix_->Zero();
    scatramanifold_scatra_matrix_->Zero();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::UTILS::SSIVectors::SSIVectors(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const bool is_scatra_manifold)
    : increment_(Core::LinAlg::CreateVector(*(ssi_maps->MapsSubProblems()->FullMap()), true)),
      is_scatra_manifold_(is_scatra_manifold),
      manifold_residual_(is_scatra_manifold ? Core::LinAlg::CreateVector(
                                                  *(ssi_maps->sca_tra_manifold_dof_row_map()), true)
                                            : Teuchos::null),
      residual_(Core::LinAlg::CreateVector(*(ssi_maps->MapsSubProblems()->FullMap()), true)),
      scatra_residual_(Core::LinAlg::CreateVector(*(ssi_maps->ScaTraDofRowMap()), true)),
      structure_residual_(Core::LinAlg::CreateVector(*(ssi_maps->StructureDofRowMap()), true))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIVectors::ClearIncrement() { increment_->PutScalar(0.0); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::UTILS::SSIVectors::ClearResiduals()
{
  residual_->PutScalar(0.0);
  scatra_residual_->PutScalar(0.0);
  structure_residual_->PutScalar(0.0);
  if (is_scatra_manifold_) manifold_residual_->PutScalar(0.0);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMatrices::initialize_system_matrix(
    Teuchos::RCP<const SSI::UTILS::SSIMaps> ssi_maps, const Core::LinAlg::MatrixType ssi_matrixtype)
{
  switch (ssi_matrixtype)
  {
    case Core::LinAlg::MatrixType::block_field:
    {
      system_matrix_ = setup_block_matrix(
          ssi_maps->block_map_system_matrix(), ssi_maps->block_map_system_matrix());
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      system_matrix_ = setup_sparse_matrix(ssi_maps->MapSystemMatrix());
      break;
    }

    default:
    {
      FOUR_C_THROW("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> SSI::UTILS::SSIMatrices::setup_block_matrix(
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> row_map,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> col_map)
{
  const int expected_entries_per_row = 81;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
      *col_map, *row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> SSI::UTILS::SSIMatrices::setup_sparse_matrix(
    const Teuchos::RCP<const Epetra_Map> row_map)
{
  const int expected_entries_per_row = 27;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      *row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/* ----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::UTILS::SSIMaps::SSIMaps(const SSI::SsiMono& ssi_mono_algorithm)
    : block_maps_sub_problems_(),
      scatra_matrixtype_(ssi_mono_algorithm.ScaTraField()->MatrixType()),
      scatra_manifold_matrixtype_(ssi_mono_algorithm.is_sca_tra_manifold()
                                      ? ssi_mono_algorithm.ScaTraManifold()->MatrixType()
                                      : Core::LinAlg::MatrixType::undefined),
      ssi_matrixtype_(ssi_mono_algorithm.MatrixType())
{
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(
      ssi_mono_algorithm.is_sca_tra_manifold() ? 3 : 2, Teuchos::null);
  Teuchos::RCP<const Epetra_Map> merged_map;

  partial_maps[GetProblemPosition(Subproblem::scalar_transport)] =
      Teuchos::rcp(new Epetra_Map(*ssi_mono_algorithm.ScaTraField()->dof_row_map()));
  partial_maps[GetProblemPosition(Subproblem::structure)] =
      Teuchos::rcp(new Epetra_Map(*ssi_mono_algorithm.structure_field()->dof_row_map()));
  if (ssi_mono_algorithm.is_sca_tra_manifold())
  {
    partial_maps[GetProblemPosition(Subproblem::manifold)] =
        Teuchos::rcp(new Epetra_Map(*ssi_mono_algorithm.ScaTraManifold()->dof_row_map()));
    auto temp_map = Core::LinAlg::MergeMap(partial_maps[0], partial_maps[1], false);
    merged_map = Core::LinAlg::MergeMap(temp_map, partial_maps[2], false);
  }
  else
    merged_map = Core::LinAlg::MergeMap(partial_maps[0], partial_maps[1], false);

  maps_sub_problems_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_sub_problems_->check_for_valid_map_extractor();

  switch (ssi_matrixtype_)
  {
    case Core::LinAlg::MatrixType::block_field:
    {
      auto block_map_structure = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
          *ssi_mono_algorithm.structure_field()->discretization()->dof_row_map(),
          std::vector<Teuchos::RCP<const Epetra_Map>>(
              1, ssi_mono_algorithm.structure_field()->dof_row_map())));

      block_map_structure->check_for_valid_map_extractor();

      block_maps_sub_problems_.insert(std::make_pair(Subproblem::structure, block_map_structure));

      switch (scatra_matrixtype_)
      {
        case Core::LinAlg::MatrixType::sparse:
        {
          auto block_map_scatra = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
              *ssi_mono_algorithm.ScaTraField()->discretization()->dof_row_map(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(
                  1, ssi_mono_algorithm.ScaTraField()->dof_row_map())));

          block_map_scatra->check_for_valid_map_extractor();

          block_maps_sub_problems_.insert(
              std::make_pair(Subproblem::scalar_transport, block_map_scatra));

          if (ssi_mono_algorithm.is_sca_tra_manifold())
          {
            auto block_map_scatra_manifold = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
                *ssi_mono_algorithm.ScaTraManifold()->discretization()->dof_row_map(),
                std::vector<Teuchos::RCP<const Epetra_Map>>(
                    1, ssi_mono_algorithm.ScaTraManifold()->dof_row_map())));

            block_map_scatra_manifold->check_for_valid_map_extractor();

            block_maps_sub_problems_.insert(
                std::make_pair(Subproblem::manifold, block_map_scatra_manifold));
          }
          break;
        }
        case Core::LinAlg::MatrixType::block_condition:
        case Core::LinAlg::MatrixType::block_condition_dof:
        {
          block_maps_sub_problems_.insert(std::make_pair(
              Subproblem::scalar_transport, ssi_mono_algorithm.ScaTraField()->BlockMaps()));

          if (ssi_mono_algorithm.is_sca_tra_manifold())
          {
            block_maps_sub_problems_.insert(std::make_pair(
                Subproblem::manifold, ssi_mono_algorithm.ScaTraManifold()->BlockMaps()));
          }
          else
          {
            block_maps_sub_problems_.insert(std::make_pair(Subproblem::manifold, Teuchos::null));
          }

          break;
        }
        default:
        {
          FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
          break;
        }
      }

      create_and_check_block_maps_sub_problems(ssi_mono_algorithm);

      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      block_maps_sub_problems_.insert(std::make_pair(Subproblem::structure, Teuchos::null));
      block_maps_sub_problems_.insert(std::make_pair(Subproblem::scalar_transport, Teuchos::null));
      block_maps_sub_problems_.insert(std::make_pair(Subproblem::manifold, Teuchos::null));
      break;
    }

    default:
    {
      FOUR_C_THROW("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  map_system_matrix_ = maps_sub_problems_->FullMap();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
std::vector<int> SSI::UTILS::SSIMaps::GetBlockPositions(Subproblem subproblem) const
{
  FOUR_C_ASSERT(
      ssi_matrixtype_ != Core::LinAlg::MatrixType::sparse, "Sparse matrices have just one block");

  std::vector<int> block_position;

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      if (scatra_matrixtype_ == Core::LinAlg::MatrixType::sparse)
        block_position.emplace_back(1);
      else
        block_position.emplace_back(BlockMapScaTra()->NumMaps());
      break;
    }
    case Subproblem::scalar_transport:
    {
      if (scatra_matrixtype_ == Core::LinAlg::MatrixType::sparse)
        block_position.emplace_back(0);
      else
      {
        for (int i = 0; i < BlockMapScaTra()->NumMaps(); ++i) block_position.emplace_back(i);
      }
      break;
    }
    case Subproblem::manifold:
    {
      if (scatra_manifold_matrixtype_ == Core::LinAlg::MatrixType::sparse)
        block_position.emplace_back(2);
      else
      {
        auto scatra_manifold_num_block_maps = block_map_sca_tra_manifold()->NumMaps();

        for (int i = 0; i < scatra_manifold_num_block_maps; ++i)
          block_position.emplace_back(BlockMapScaTra()->NumMaps() + 1 + i);
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of subproblem");
      break;
    }
  }

  return block_position;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
int SSI::UTILS::SSIMaps::GetProblemPosition(Subproblem subproblem)
{
  int position = -1;

  switch (subproblem)
  {
    case Subproblem::structure:
    {
      position = 1;
      break;
    }
    case Subproblem::scalar_transport:
    {
      position = 0;
      break;
    }
    case Subproblem::manifold:
    {
      position = 2;
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of sub problem");
      break;
    }
  }

  return position;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMaps::create_and_check_block_maps_sub_problems(
    const SSI::SsiMono& ssi_mono_algorithm)
{
  const int num_blocks_systemmatrix =
      BlockMapScaTra()->NumMaps() + BlockMapStructure()->NumMaps() +
      (ssi_mono_algorithm.is_sca_tra_manifold() ? block_map_sca_tra_manifold()->NumMaps() : 0);

  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps_system_matrix(
      num_blocks_systemmatrix, Teuchos::null);

  for (int i = 0; i < BlockMapScaTra()->NumMaps(); ++i)

  {
    auto block_positions_scatra = GetBlockPositions(Subproblem::scalar_transport);
    partial_maps_system_matrix[block_positions_scatra.at(i)] = BlockMapScaTra()->Map(i);
  }

  partial_maps_system_matrix.at(GetBlockPositions(Subproblem::structure).at(0)) =
      BlockMapStructure()->FullMap();

  if (ssi_mono_algorithm.is_sca_tra_manifold())
  {
    for (int i = 0; i < block_map_sca_tra_manifold()->NumMaps(); ++i)
    {
      auto block_positions_manifold = GetBlockPositions(Subproblem::manifold);
      partial_maps_system_matrix[block_positions_manifold.at(i)] =
          block_map_sca_tra_manifold()->Map(i);
    }
  }

  block_map_system_matrix_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
      *maps_sub_problems_->FullMap(), partial_maps_system_matrix));

  block_map_system_matrix_->check_for_valid_map_extractor();
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> SSI::UTILS::SSIMaps::BlockMapScaTra() const
{
  return block_maps_sub_problems_.at(Subproblem::scalar_transport);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MultiMapExtractor>
SSI::UTILS::SSIMaps::block_map_sca_tra_manifold() const
{
  return block_maps_sub_problems_.at(Subproblem::manifold);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> SSI::UTILS::SSIMaps::BlockMapStructure() const
{
  return block_maps_sub_problems_.at(Subproblem::structure);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMaps::ScaTraDofRowMap() const
{
  return MapsSubProblems()->Map(SSIMaps::GetProblemPosition(Subproblem::scalar_transport));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMaps::sca_tra_manifold_dof_row_map() const
{
  return MapsSubProblems()->Map(SSIMaps::GetProblemPosition(Subproblem::manifold));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SSI::UTILS::SSIMaps::StructureDofRowMap() const
{
  return MapsSubProblems()->Map(SSIMaps::GetProblemPosition(Subproblem::structure));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::CheckConsistencyOfSSIInterfaceContactCondition(
    const std::vector<Core::Conditions::Condition*>& conditionsToBeTested,
    Teuchos::RCP<Core::FE::Discretization>& structdis)
{
  // get conditions to check against
  std::vector<Core::Conditions::Condition*> s2ikinetics_conditions;
  structdis->GetCondition("S2IKinetics", s2ikinetics_conditions);
  std::vector<Core::Conditions::Condition*> contactconditions;
  structdis->GetCondition("Contact", contactconditions);

  // loop over all ssi conditions and check them
  for (const auto* conditionToBeTested : conditionsToBeTested)
  {
    std::vector<Core::Conditions::Condition*> InterfaceS2IConditions;
    std::vector<Core::Conditions::Condition*> InterfaceContactConditions;

    const int S2IKineticsID = conditionToBeTested->parameters().Get<int>("S2IKineticsID");
    const int contactconditionID = conditionToBeTested->parameters().Get<int>("ContactConditionID");

    if (S2IKineticsID != contactconditionID)
    {
      FOUR_C_THROW(
          "For the 'SSIInterfaceContact' condition we have to demand, that the 'S2ICouplingID' and "
          "the 'ContactConditionID' have the same value as the contact strategy factory relies on "
          "this to set the scatra interface parameters correctly.");
    }

    // loop over all scatra-scatra interface conditions and add them to the vector, if IDs match
    for (auto* s2ikinetics_cond : s2ikinetics_conditions)
    {
      if (s2ikinetics_cond->parameters().Get<int>("ConditionID") != S2IKineticsID) continue;

      InterfaceS2IConditions.push_back(s2ikinetics_cond);
    }

    // loop over all contact conditions and add them to the vector, if IDs match
    for (auto* contactcondition : contactconditions)
    {
      if (contactcondition->parameters().Get<int>("Interface ID") != contactconditionID) continue;

      InterfaceContactConditions.push_back(contactcondition);
    }

    if (InterfaceContactConditions.empty())
      FOUR_C_THROW(
          "Did not find 'Contact' condition as defined in 'SSIInterfaceContact' condition!");

    if (InterfaceS2IConditions.empty())
      FOUR_C_THROW(
          "Did not find 'S2ICoupling' condition as defined in 'SSIInterfaceContact' condition!");

    // now get the nodes
    std::vector<int> InterfaceS2INodes;
    std::vector<int> InterfaceContactNodes;
    Core::Conditions::FindConditionedNodes(*structdis, InterfaceS2IConditions, InterfaceS2INodes);
    Core::Conditions::FindConditionedNodes(
        *structdis, InterfaceContactConditions, InterfaceContactNodes);

    // and compare whether same nodes are defined
    for (const auto InterfaceS2INode : InterfaceS2INodes)
    {
      bool foundit(false);

      for (const auto InterfaceContactNode : InterfaceContactNodes)
      {
        if (InterfaceS2INode == InterfaceContactNode) foundit = true;
      }

      if (!foundit)
      {
        FOUR_C_THROW(
            "The following conditions are defined on a non-consistent set of nodes!\n"
            "The Conditions defined from 'SSIInterfaceContact' condition with the condition-ID: "
            "%i:\n"
            "'S2ICoupling' conditions with ID: %i;\n"
            "'Contact' conditions with ID: %i;\n"
            "The last two conditions are NOT defined on the same node-sets which is not "
            "reasonable. Check your Input-File!",
            conditionToBeTested->parameters().Get<int>("ConditionID"), S2IKineticsID,
            contactconditionID);
      }
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::UTILS::SSIMeshTying::SSIMeshTying(const std::string& conditionname_coupling,
    Teuchos::RCP<Core::FE::Discretization> dis, const bool build_slave_slave_transformation,
    const bool check_over_constrained)
    : comm_(dis->Comm()),
      do_print_(dis->Comm().MyPID() == 0),
      meshtying_handlers_(),
      my_rank_(dis->Comm().MyPID()),
      num_proc_(dis->Comm().NumProc())
{
  setup_mesh_tying_handlers(
      dis, conditionname_coupling, build_slave_slave_transformation, check_over_constrained);

  // construct full slave, master, and interior maps
  std::vector<Teuchos::RCP<const Epetra_Map>> slave_maps;
  std::vector<Teuchos::RCP<const Epetra_Map>> master_maps;
  for (const auto& meshtying : meshtying_handlers_)
  {
    slave_maps.emplace_back(meshtying->SlaveMasterCoupling()->SlaveDofMap());
    master_maps.emplace_back(meshtying->SlaveMasterCoupling()->MasterDofMap());
  }

  full_master_side_map_ = Core::LinAlg::MultiMapExtractor::MergeMaps(master_maps);
  full_slave_side_map_ = Core::LinAlg::MultiMapExtractor::MergeMaps(slave_maps);
  auto interface_map = Core::LinAlg::MergeMap(full_master_side_map_, full_slave_side_map_);

  interior_map_ = Core::LinAlg::SplitMap(*dis->dof_row_map(), *interface_map);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::setup_mesh_tying_handlers(Teuchos::RCP<Core::FE::Discretization> dis,
    const std::string& name_meshtying_condition, const bool build_slave_slave_transformation,
    const bool check_over_constrained)
{
  auto timer = Teuchos::rcp(new Teuchos::Time("MeshtypingHandlers", true));
  const double t0 = timer->wallTime();

  if (do_print_)
  {
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "| Starting to setup mesh tying   |" << std::endl;
  }

  // find pairwise matching nodes
  std::vector<std::pair<int, int>> coupling_pairs;
  find_matching_node_pairs(dis, name_meshtying_condition, coupling_pairs);

  // all nodes on one geometrical point: outer vector defines geometric position, and inner vector
  // all nodes there
  std::vector<std::vector<int>> grouped_matching_nodes;
  group_matching_nodes(coupling_pairs, grouped_matching_nodes);

  if (do_print_) std::cout << "| Finished: group nodes          |" << std::endl;

  // define master gids and pair them with the remaining (slave) nodes
  std::vector<int> master_gids;
  std::map<int, int> slave_master_pair;
  define_master_slave_pairing(
      dis, grouped_matching_nodes, master_gids, slave_master_pair, check_over_constrained);

  if (do_print_) std::cout << "| Finished: master-slave pairing |" << std::endl;

  // get number of slave nodes per master node -> max. number gives number of needed adapters
  std::map<int, int> num_assigned_slave_to_master_nodes;
  int glob_max_adapters = 0;
  get_num_assigned_slave_to_master_nodes(
      slave_master_pair, num_assigned_slave_to_master_nodes, glob_max_adapters);

  // setup mesh tying:
  // - coupling adapter
  // - map extractor
  // - slave_slave_transformation
  std::map<int, int> created_adapters;  // numbers of created coupling adapters per master node
  for (const auto& item : num_assigned_slave_to_master_nodes)
    created_adapters.insert(std::make_pair(item.first, 0));

  for (int iadapter = 0; iadapter < glob_max_adapters; ++iadapter)
  {
    // create vectors of master and slave nodes for this coupling adapter
    std::vector<int> inodegidvec_master;
    std::vector<int> inodegidvec_slave;

    for (auto& pair : slave_master_pair)
    {
      const int master_gid = pair.second;
      if (master_gid != -1)  // this pair is already considered
      {
        const int num_created_adapters = created_adapters[master_gid];
        if (num_assigned_slave_to_master_nodes[master_gid] >= iadapter + 1 and
            num_created_adapters == iadapter)
        {
          const int slave_gid = pair.first;
          Core::Communication::AddOwnedNodeGID(*dis, master_gid, inodegidvec_master);
          Core::Communication::AddOwnedNodeGID(*dis, slave_gid, inodegidvec_slave);

          pair.second = -1;  // do not consider this pair in next iteration
          created_adapters[master_gid] = num_created_adapters + 1;
        }
      }
    }

    // setup coupling adapter
    auto coupling_adapter = Teuchos::rcp(new Core::Adapter::Coupling());
    const int num_dofs =
        static_cast<int>(static_cast<double>(dis->dof_row_map()->NumGlobalElements()) /
                         static_cast<double>(dis->NodeRowMap()->NumGlobalElements()));
    coupling_adapter->setup_coupling(
        *dis, *dis, inodegidvec_master, inodegidvec_slave, num_dofs, true, 1.0e-8);

    // setup multimap extractor for each coupling adapter
    auto slave_map = coupling_adapter->SlaveDofMap();
    auto master_map = coupling_adapter->MasterDofMap();
    auto interior_map =
        Core::LinAlg::SplitMap(*dis->dof_row_map(), *Core::LinAlg::MergeMap(slave_map, master_map));

    std::vector<Teuchos::RCP<const Epetra_Map>> maps(0, Teuchos::null);
    maps.emplace_back(interior_map);
    maps.emplace_back(slave_map);
    maps.emplace_back(master_map);

    auto coupling_map_extractor =
        Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*dis->dof_row_map(), maps));
    coupling_map_extractor->check_for_valid_map_extractor();

    auto slave_slave_transformation = Teuchos::rcp(new Core::Adapter::Coupling());
    if (build_slave_slave_transformation)
    {
      // coupling adapter between new slave nodes (master) and old slave nodes from input file
      // (slave). Needed, to map the lineariaztion dphi/du at the interface to the correct dofs
      std::vector<int> all_coupled_original_slave_gids;
      find_slave_slave_transformation_nodes(
          dis, name_meshtying_condition, inodegidvec_slave, all_coupled_original_slave_gids);

      std::vector<int> my_coupled_original_slave_gids;
      for (const int& slave_gid : all_coupled_original_slave_gids)
        Core::Communication::AddOwnedNodeGID(*dis, slave_gid, my_coupled_original_slave_gids);

      slave_slave_transformation->setup_coupling(*dis, *dis, inodegidvec_slave,
          my_coupled_original_slave_gids, Global::Problem::Instance()->NDim(), true, 1.0e-8);
    }

    // combine coupling adapters and multimap extractor to mesh tying object
    meshtying_handlers_.emplace_back(Teuchos::rcp(new SSIMeshTyingHandler(
        coupling_adapter, coupling_map_extractor, slave_slave_transformation)));
  }

  const double total_time = timer->wallTime() - t0;
  if (do_print_)
  {
    std::cout << "|--------------------------------|" << std::endl;
    std::cout << "|  Mesh tying setup successful   |" << std::endl;
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "|        (" << total_time << " sec.)        |" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl << std::endl;
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
int SSI::UTILS::SSIMeshTying::has_gid(
    const int gid, const std::vector<std::vector<int>>& matching_nodes) const
{
  const int size = static_cast<int>(matching_nodes.size());
  const int load_per_proc = std::floor(static_cast<double>(size) / static_cast<double>(num_proc_));

  // define load of this proc
  const int upper_bound = my_rank_ == num_proc_ - 1 ? size : load_per_proc * (my_rank_ + 1);
  const int lower_bound = load_per_proc * my_rank_;

  int my_return_value = has_gid_partial(gid, lower_bound, upper_bound, matching_nodes);
  int glob_return_value;

  comm_.MaxAll(&my_return_value, &glob_return_value, 1);

  return glob_return_value;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
int SSI::UTILS::SSIMeshTying::has_gid_partial(const int gid, const int start, const int end,
    const std::vector<std::vector<int>>& matching_nodes) const
{
  for (int i = start; i < end; ++i)
  {
    for (const int& j : matching_nodes[i])
      if (gid == j) return i;
  }
  return -1;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::find_matching_node_pairs(Teuchos::RCP<Core::FE::Discretization> dis,
    const std::string& name_meshtying_condition,
    std::vector<std::pair<int, int>>& coupling_pairs) const
{
  // coupled nodes on this proc from input
  std::vector<std::pair<int, int>> my_coupling_pairs;

  // get all mesh tying conditions
  std::vector<Core::Conditions::Condition*> meshtying_conditons(0, nullptr);
  dis->GetCondition(name_meshtying_condition, meshtying_conditons);

  // match nodes between all mesh tying conditons (named with "a" and "b")
  for (int a = 0; a < static_cast<int>(meshtying_conditons.size()); ++a)
  {
    auto* meshtying_condition_a = meshtying_conditons.at(a);

    // nodes of meshtying_condition_a owned by this proc
    std::vector<int> inodegidvec_a;
    Core::Communication::AddOwnedNodeGIDFromList(
        *dis, *meshtying_condition_a->GetNodes(), inodegidvec_a);

    // init node matching octree with nodes from condition a
    auto tree = Core::COUPLING::NodeMatchingOctree();
    tree.Init(*dis, inodegidvec_a, 150, 1.0e-8);
    tree.Setup();

    // find nodes from condition b that match nodes from condition a
    for (int b = a + 1; b < static_cast<int>(meshtying_conditons.size()); ++b)
    {
      auto* meshtying_condition_b = meshtying_conditons.at(b);

      // nodes of meshtying_condition_b owned by this proc
      std::vector<int> inodegidvec_b;
      Core::Communication::AddOwnedNodeGIDFromList(
          *dis, *meshtying_condition_b->GetNodes(), inodegidvec_b);

      // key: master node gid, value: slave node gid and distance
      std::map<int, std::pair<int, double>> coupled_gid_nodes;
      tree.FindMatch(*dis, inodegidvec_b, coupled_gid_nodes);

      // loop over all nodal couplings and find coupled nodes
      for (const auto& pair : coupled_gid_nodes)
      {
        const double distance = pair.second.second;
        if (distance < 1.0e-16)
        {
          const int gid1 = pair.first;
          const int gid2 = pair.second.first;

          my_coupling_pairs.emplace_back(std::make_pair(gid1, gid2));
        }
      }
    }
  }

  // communicate to all other procs
  const auto all_coupling_pairs =
      Core::Communication::BroadcastPairVector(my_coupling_pairs, comm_);

  // remove duplicates (slave node = master node)
  for (const auto& pair : all_coupling_pairs)
    if (pair.first != pair.second) coupling_pairs.emplace_back(pair);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::group_matching_nodes(
    const std::vector<std::pair<int, int>>& coupling_pairs,
    std::vector<std::vector<int>>& grouped_matching_nodes) const
{
  // all nodes on one geometrical point: outer vector defines geometric position, and inner vector
  // all nodes there
  for (const auto& pair : coupling_pairs)
  {
    const int gid1 = pair.first;
    const int gid2 = pair.second;

    std::vector<int> nodes_at_same_point;

    // initial fill of grouped_matching_nodes
    if (grouped_matching_nodes.empty())
    {
      nodes_at_same_point.resize(2);
      nodes_at_same_point[0] = gid1;
      nodes_at_same_point[1] = gid2;
      grouped_matching_nodes.emplace_back(nodes_at_same_point);
    }
    else
    {
      // check if gid1 or gid2 is already part of the matching_nodes vector
      const int index1 = has_gid(gid1, grouped_matching_nodes);
      const int index2 = has_gid(gid2, grouped_matching_nodes);

      // do not add, if both indices are equal, meaning that the different gids are already in the
      // same group
      if (index1 == index2 and index1 != -1)
      {
        continue;
      }
      // if both gids are not part of matching_nodes -> create new entry
      else if (index1 == -1 and index2 == -1)
      {
        nodes_at_same_point.resize(2);
        nodes_at_same_point[0] = gid1;
        nodes_at_same_point[1] = gid2;
        grouped_matching_nodes.emplace_back(nodes_at_same_point);
      }
      // if gid1 is part of matching_nodes and gid2 is not -> add gid2 to entry of gid1
      else if (index1 != -1 and index2 == -1)
      {
        grouped_matching_nodes[index1].emplace_back(gid2);
      }
      // if gid2 is part of matching_nodes and gid1 is not -> add gid1 to entry of gid2
      else if (index1 == -1 and index2 != -1)
      {
        grouped_matching_nodes[index2].emplace_back(gid1);
      }
      // if gid1 and gid2 are part of matching_nodes -> remove nodes from entry of gid1 and move
      // them to entry of gid2
      else if (index1 != -1 and index2 != -1 and index1 != index2)
      {
        for (const int& node : grouped_matching_nodes[index1])
          grouped_matching_nodes[index2].emplace_back(node);

        grouped_matching_nodes.erase(std::next(grouped_matching_nodes.begin(), index1));
      }
      else
        FOUR_C_THROW("This should not be the case");
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::get_num_assigned_slave_to_master_nodes(
    const std::map<int, int>& slave_master_pair,
    std::map<int, int>& num_assigned_slave_to_master_nodes, int& max_assigned_slave_nodes) const
{
  for (const auto& pair : slave_master_pair)
  {
    const int master_node_gid = pair.second;

    // initial fill
    if (num_assigned_slave_to_master_nodes.empty())
    {
      num_assigned_slave_to_master_nodes.insert(std::make_pair(master_node_gid, 1));
      max_assigned_slave_nodes = 1;
    }
    else
    {
      // master node not in map so far -> create entry
      if (num_assigned_slave_to_master_nodes.find(master_node_gid) !=
          num_assigned_slave_to_master_nodes.end())
      {
        num_assigned_slave_to_master_nodes[master_node_gid]++;
        if (max_assigned_slave_nodes < num_assigned_slave_to_master_nodes[master_node_gid])
          max_assigned_slave_nodes = num_assigned_slave_to_master_nodes[master_node_gid];
      }
      // master node already in map-> increase counter by one
      else
        num_assigned_slave_to_master_nodes.insert(std::make_pair(master_node_gid, 1));
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::define_master_slave_pairing(
    Teuchos::RCP<Core::FE::Discretization> dis,
    const std::vector<std::vector<int>>& grouped_matching_nodes, std::vector<int>& master_gids,
    std::map<int, int>& slave_master_pair, const bool check_over_constrained) const
{
  // get Dirichlet nodes -> they define the master side
  std::vector<Core::Conditions::Condition*> dbc_conds;
  dis->GetCondition("Dirichlet", dbc_conds);
  std::set<int> dbc_nodes;
  for (auto* dbc_cond : dbc_conds)
    for (const int& dbc_node : *dbc_cond->GetNodes()) dbc_nodes.insert(dbc_node);

  std::vector<int> my_master_gids;
  std::map<int, int> my_slave_master_pair;

  const int size = static_cast<int>(grouped_matching_nodes.size());
  const int load_per_proc = std::ceil(static_cast<double>(size) / static_cast<double>(num_proc_));

  // define load for this proc in for loop
  const int upper_bound = my_rank_ == num_proc_ - 1 ? size : load_per_proc * (my_rank_ + 1);
  const int lower_bound = load_per_proc * my_rank_;

  for (int i = lower_bound; i < upper_bound; ++i)
  {
    const auto& nodes_at_same_point = grouped_matching_nodes[i];

    // define which of the nodes from nodes_at_same_point is master node
    // in case one node is a Dirichlet node -> use this one
    bool found_dbc_node = false;
    int new_master_gid = -1;
    for (const int& node : nodes_at_same_point)
    {
      for (const int& dbc_node : dbc_nodes)
      {
        if (dbc_node == node)
        {
          if (found_dbc_node and check_over_constrained)
            FOUR_C_THROW("Mesh tying and DBCs are over constrained.");
          else
          {
            found_dbc_node = true;
            new_master_gid = dbc_node;
            break;
          }
#ifndef FOUR_C_ENABLE_ASSERTIONS
          break;  // in release version, we can break here and subsequently skip the safety check
#endif
        }
      }
    }
    // otherwise -> use the first (random) node
    if (!found_dbc_node) new_master_gid = nodes_at_same_point[0];

    // insert new master node to list of master nodes
    my_master_gids.emplace_back(new_master_gid);

    // build slave-master pairs with defined master node
    for (const int& node : nodes_at_same_point)
      if (node != new_master_gid) my_slave_master_pair.insert(std::make_pair(node, new_master_gid));
  }

  master_gids = Core::Communication::BroadcastVector(my_master_gids, comm_);
  slave_master_pair = Core::Communication::BroadcastMap(my_slave_master_pair, comm_);

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check if everything worked fine
  std::set<int> unique_master_gids;
  for (const int& master_gid : master_gids) unique_master_gids.insert(master_gid);
  if (unique_master_gids.size() != master_gids.size()) FOUR_C_THROW("Master gids not unique");
#endif
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::find_slave_slave_transformation_nodes(
    Teuchos::RCP<Core::FE::Discretization> dis, const std::string& name_meshtying_condition,
    const std::vector<int>& inodegidvec_slave,
    std::vector<int>& all_coupled_original_slave_gids) const
{
  // store nodes that are slave nodes from the input
  std::vector<Core::Conditions::Condition*> meshtying_conditons(0, nullptr);
  dis->GetCondition(name_meshtying_condition, meshtying_conditons);

  std::vector<int> original_slave_gids;
  for (auto* meshtying_conditon : meshtying_conditons)
  {
    if (meshtying_conditon->parameters().Get<int>("interface side") == Inpar::S2I::side_slave)
    {
      Core::Communication::AddOwnedNodeGIDFromList(
          *dis, *meshtying_conditon->GetNodes(), original_slave_gids);
    }
  }

  auto tree = Core::COUPLING::NodeMatchingOctree();
  tree.Init(*dis, inodegidvec_slave, 150, 1.0e-8);
  tree.Setup();
  std::map<int, std::pair<int, double>> coupled_gid_nodes;
  tree.FindMatch(*dis, original_slave_gids, coupled_gid_nodes);

  // find matches if distance < 1.0e-16
  std::vector<int> my_coupled_original_slave_gids;
  for (const auto& pair : coupled_gid_nodes)
  {
    const double distance = pair.second.second;
    const int gid = pair.second.first;
    if (distance < 1.0e-16) my_coupled_original_slave_gids.emplace_back(gid);
  }

  // distribute gids from original slave nodes to all procs (matching might be on different proc)
  all_coupled_original_slave_gids =
      Core::Communication::BroadcastVector(my_coupled_original_slave_gids, comm_);
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::UTILS::SSIMeshTying::check_slave_side_has_dirichlet_conditions(
    Teuchos::RCP<const Epetra_Map> struct_dbc_map) const
{
  // check if slave side dofs are part of DBC maps
  std::vector<Teuchos::RCP<const Epetra_Map>> maps(2, Teuchos::null);
  maps[0] = struct_dbc_map;
  for (const auto& meshtying : meshtying_handlers_)
  {
    maps[1] = meshtying->SlaveMasterCoupling()->SlaveDofMap();
    if (Core::LinAlg::MultiMapExtractor::IntersectMaps(maps)->NumGlobalElements() > 0)
      FOUR_C_THROW("Must not apply Dirichlet conditions to slave-side structural displacements!");
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::UTILS::SSIMeshTyingHandler::SSIMeshTyingHandler(
    Teuchos::RCP<Core::Adapter::Coupling> slave_master_coupling,
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> slave_master_extractor,
    Teuchos::RCP<Core::Adapter::Coupling> slave_slave_transformation)
    : slave_master_coupling_(std::move(slave_master_coupling)),
      slave_master_extractor_(std::move(slave_master_extractor)),
      slave_slave_transformation_(std::move(slave_slave_transformation))
{
  slave_side_converter_ =
      Teuchos::rcp(new Core::Adapter::CouplingSlaveConverter(*slave_master_coupling_));
}

FOUR_C_NAMESPACE_CLOSE
