/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSTI

 \level 2


 *------------------------------------------------------------------------------------------------*/


#include "baci_ssti_utils.H"

#include "baci_adapter_str_ssiwrapper.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_par_material.H"
#include "baci_scatra_ele.H"
#include "baci_scatra_timint_implicit.H"
#include "baci_scatra_timint_meshtying_strategy_s2i.H"
#include "baci_so3_nurbs27.H"
#include "baci_ssi_utils.H"
#include "baci_ssti_monolithic.H"

#include <Epetra_Map.h>

BACI_NAMESPACE_OPEN



/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSTI::SSTIMaps::SSTIMaps(const SSTI::SSTIMono& ssti_mono_algorithm)
    : block_map_scatra_(Teuchos::null),
      block_map_structure_(Teuchos::null),
      block_map_thermo_(Teuchos::null),
      maps_subproblems_(Teuchos::null)
{
  // setup maps containing dofs of subproblems
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps(3, Teuchos::null);
  partial_maps[ssti_mono_algorithm.GetProblemPosition(Subproblem::scalar_transport)] =
      Teuchos::rcp(new Epetra_Map(*ssti_mono_algorithm.ScaTraField()->DofRowMap()));
  partial_maps[ssti_mono_algorithm.GetProblemPosition(Subproblem::structure)] =
      Teuchos::rcp(new Epetra_Map(*ssti_mono_algorithm.StructureField()->DofRowMap()));
  partial_maps[ssti_mono_algorithm.GetProblemPosition(Subproblem::thermo)] =
      Teuchos::rcp(new Epetra_Map(*ssti_mono_algorithm.ThermoField()->DofRowMap()));
  Teuchos::RCP<const Epetra_Map> temp_map =
      CORE::LINALG::MergeMap(partial_maps[0], partial_maps[1], false);
  Teuchos::RCP<const Epetra_Map> merged_map =
      CORE::LINALG::MergeMap(temp_map, partial_maps[2], false);
  // initialize global map extractor
  maps_subproblems_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_subproblems_->CheckForValidMapExtractor();

  // initialize map extractors associated with blocks of subproblems
  block_map_structure_ = Teuchos::rcp(
      new CORE::LINALG::MultiMapExtractor(*ssti_mono_algorithm.StructureField()->DofRowMap(),
          std::vector<Teuchos::RCP<const Epetra_Map>>(
              1, ssti_mono_algorithm.StructureField()->DofRowMap())));
  switch (ssti_mono_algorithm.ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      block_map_scatra_ = Teuchos::rcp(
          new CORE::LINALG::MultiMapExtractor(*ssti_mono_algorithm.ScaTraField()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(
                  1, ssti_mono_algorithm.ScaTraField()->DofRowMap())));
      block_map_thermo_ = Teuchos::rcp(
          new CORE::LINALG::MultiMapExtractor(*ssti_mono_algorithm.ThermoField()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(
                  1, ssti_mono_algorithm.ThermoField()->DofRowMap())));
      break;
    }
    case CORE::LINALG::MatrixType::block_condition:
    {
      block_map_scatra_ = ssti_mono_algorithm.ScaTraField()->BlockMaps();
      block_map_thermo_ = ssti_mono_algorithm.ThermoField()->BlockMaps();
      break;
    }
    default:
    {
      dserror("Matrix type not supported");
      break;
    }
  }

  block_map_scatra_->CheckForValidMapExtractor();
  block_map_structure_->CheckForValidMapExtractor();
  block_map_thermo_->CheckForValidMapExtractor();
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> SSTI::SSTIMaps::MapInterface(
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy) const
{
  auto mergedInterfaceMap = CORE::LINALG::MultiMapExtractor::MergeMaps(
      {meshtyingstrategy->CouplingAdapter()->MasterDofMap(),
          meshtyingstrategy->CouplingAdapter()->SlaveDofMap()});
  if (not mergedInterfaceMap->UniqueGIDs()) dserror("Map not unique");
  return mergedInterfaceMap;
}


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::MultiMapExtractor> SSTI::SSTIMaps::MapsInterfaceBlocks(
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
    CORE::LINALG::MatrixType scatramatrixtype, unsigned nummaps) const
{
  Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockmapinterface(Teuchos::null);

  Teuchos::RCP<Epetra_Map> interfacemap = MapInterface(meshtyingstrategy);

  switch (scatramatrixtype)
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      blockmapinterface = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(
          *interfacemap, std::vector<Teuchos::RCP<const Epetra_Map>>(1, interfacemap)));
      break;
    }
    case CORE::LINALG::MatrixType::block_condition:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> partial_blockmapinterface(nummaps, Teuchos::null);
      for (int iblockmap = 0; iblockmap < static_cast<int>(nummaps); ++iblockmap)
      {
        partial_blockmapinterface[iblockmap] = CORE::LINALG::MultiMapExtractor::MergeMaps(
            {meshtyingstrategy->BlockMapsSlave().Map(iblockmap),
                meshtyingstrategy->BlockMapsMaster().Map(iblockmap)});
      }
      blockmapinterface = Teuchos::rcp(
          new CORE::LINALG::MultiMapExtractor(*interfacemap, partial_blockmapinterface));
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  blockmapinterface->CheckForValidMapExtractor();

  return blockmapinterface;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::MultiMapExtractor> SSTI::SSTIMaps::MapsInterfaceBlocksSlave(
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
    CORE::LINALG::MatrixType scatramatrixtype, unsigned nummaps) const
{
  Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockmapinterfaceslave(Teuchos::null);

  switch (scatramatrixtype)
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      const auto slavedofmap = meshtyingstrategy->CouplingAdapter()->SlaveDofMap();
      blockmapinterfaceslave = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(
          *slavedofmap, std::vector<Teuchos::RCP<const Epetra_Map>>(1, slavedofmap)));
      break;
    }
    case CORE::LINALG::MatrixType::block_condition:
    {
      blockmapinterfaceslave =
          Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(meshtyingstrategy->BlockMapsSlave()));
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  blockmapinterfaceslave->CheckForValidMapExtractor();

  return blockmapinterfaceslave;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSTI::SSTIMapsMono::SSTIMapsMono(const SSTI::SSTIMono& ssti_mono_algorithm)
    : SSTIMaps(ssti_mono_algorithm), block_map_system_matrix_(Teuchos::null)
{
  // initialize map extractors associated with blocks of global system matrix
  switch (ssti_mono_algorithm.ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case CORE::LINALG::MatrixType::sparse:
    {
      block_map_system_matrix_ = MapsSubProblems();
      break;
    }
      // many main-diagonal matrix blocks associated with scalar transport field
    case CORE::LINALG::MatrixType::block_condition:
    {
      auto block_positions_scatra =
          ssti_mono_algorithm.GetBlockPositions(Subproblem::scalar_transport);
      auto block_positions_structure = ssti_mono_algorithm.GetBlockPositions(Subproblem::structure);
      auto block_positions_thermo = ssti_mono_algorithm.GetBlockPositions(Subproblem::thermo);

      std::vector<Teuchos::RCP<const Epetra_Map>> maps_systemmatrix(
          block_positions_scatra.size() + block_positions_structure.size() +
          block_positions_thermo.size());
      for (int imap = 0; imap < static_cast<int>(block_positions_scatra.size()); ++imap)
        maps_systemmatrix[block_positions_scatra.at(imap)] = BlockMapScatra()->Map(imap);

      // extract map underlying single main-diagonal matrix block associated with structural
      // field
      maps_systemmatrix[block_positions_structure.at(0)] = BlockMapStructure()->FullMap();

      for (int imap = 0; imap < static_cast<int>(block_positions_thermo.size()); ++imap)
        maps_systemmatrix[block_positions_thermo.at(imap)] = BlockMapThermo()->Map(imap);

      // initialize map extractor associated with blocks of global system matrix
      block_map_system_matrix_ = Teuchos::rcp(
          new CORE::LINALG::MultiMapExtractor(*MapsSubProblems()->FullMap(), maps_systemmatrix));

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSTI::SSTIMatrices::SSTIMatrices(Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono,
    const CORE::LINALG::MatrixType matrixtype_global,
    const CORE::LINALG::MatrixType matrixtype_scatra, bool interfacemeshtying)
    : matrixtype_scatra_(matrixtype_scatra),
      ssti_maps_mono_(ssti_maps_mono),
      systemmatrix_(Teuchos::null),
      scatrastructuredomain_(Teuchos::null),
      scatrastructureinterface_(Teuchos::null),
      scatrathermodomain_(Teuchos::null),
      scatrathermointerface_(Teuchos::null),
      structurescatradomain_(Teuchos::null),
      structurethermodomain_(Teuchos::null),
      thermoscatradomain_(Teuchos::null),
      thermoscatrainterface_(Teuchos::null),
      thermostructuredomain_(Teuchos::null),
      thermostructureinterface_(Teuchos::null),
      interfacemeshtying_(interfacemeshtying)
{
  // perform initializations associated with global system matrix
  switch (matrixtype_global)
  {
    case CORE::LINALG::MatrixType::block_field:
    {
      systemmatrix_ = SetupBlockMatrix(
          ssti_maps_mono->BlockMapSystemMatrix(), ssti_maps_mono->BlockMapSystemMatrix());
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      systemmatrix_ = SetupSparseMatrix(ssti_maps_mono->MapsSubProblems()->FullMap());
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // setup blocks for coupling matrices
  switch (matrixtype_scatra)
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      scatrastructuredomain_ =
          SetupBlockMatrix(ssti_maps_mono->BlockMapScatra(), ssti_maps_mono->BlockMapStructure());
      structurescatradomain_ =
          SetupBlockMatrix(ssti_maps_mono->BlockMapStructure(), ssti_maps_mono->BlockMapScatra());
      structurethermodomain_ =
          SetupBlockMatrix(ssti_maps_mono->BlockMapStructure(), ssti_maps_mono->BlockMapThermo());
      thermostructuredomain_ =
          SetupBlockMatrix(ssti_maps_mono->BlockMapThermo(), ssti_maps_mono->BlockMapStructure());
      scatrathermodomain_ =
          SetupBlockMatrix(ssti_maps_mono->BlockMapScatra(), ssti_maps_mono->BlockMapThermo());
      thermoscatradomain_ =
          SetupBlockMatrix(ssti_maps_mono->BlockMapThermo(), ssti_maps_mono->BlockMapScatra());

      if (interfacemeshtying_)
      {
        scatrastructureinterface_ =
            SetupBlockMatrix(ssti_maps_mono->BlockMapScatra(), ssti_maps_mono->BlockMapStructure());
        thermostructureinterface_ =
            SetupBlockMatrix(ssti_maps_mono->BlockMapThermo(), ssti_maps_mono->BlockMapStructure());
        scatrathermointerface_ =
            SetupBlockMatrix(ssti_maps_mono->BlockMapScatra(), ssti_maps_mono->BlockMapThermo());
        thermoscatrainterface_ =
            SetupBlockMatrix(ssti_maps_mono->BlockMapThermo(), ssti_maps_mono->BlockMapScatra());
      }
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      scatrastructuredomain_ = SetupSparseMatrix(ssti_maps_mono->BlockMapScatra()->FullMap());
      structurescatradomain_ = SetupSparseMatrix(ssti_maps_mono->BlockMapStructure()->FullMap());
      structurethermodomain_ = SetupSparseMatrix(ssti_maps_mono->BlockMapStructure()->FullMap());
      thermostructuredomain_ = SetupSparseMatrix(ssti_maps_mono->BlockMapThermo()->FullMap());
      scatrathermodomain_ = SetupSparseMatrix(ssti_maps_mono->BlockMapScatra()->FullMap());
      thermoscatradomain_ = SetupSparseMatrix(ssti_maps_mono->BlockMapThermo()->FullMap());

      if (interfacemeshtying_)
      {
        scatrastructureinterface_ = SetupSparseMatrix(ssti_maps_mono->BlockMapScatra()->FullMap());
        thermostructureinterface_ = SetupSparseMatrix(ssti_maps_mono->BlockMapThermo()->FullMap());
        scatrathermointerface_ = SetupSparseMatrix(ssti_maps_mono->BlockMapScatra()->FullMap());
        thermoscatrainterface_ = SetupSparseMatrix(ssti_maps_mono->BlockMapThermo()->FullMap());
      }
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSTI::SSTIMatrices::ClearMatrices()
{
  systemmatrix_->Zero();
  scatrastructuredomain_->Zero();
  scatrathermodomain_->Zero();
  structurescatradomain_->Zero();
  structurethermodomain_->Zero();
  thermoscatradomain_->Zero();
  thermostructuredomain_->Zero();

  if (interfacemeshtying_)
  {
    scatrastructureinterface_->Zero();
    scatrathermointerface_->Zero();
    thermoscatrainterface_->Zero();
    thermostructureinterface_->Zero();
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSTI::SSTIMatrices::CompleteCouplingMatrices()
{
  switch (matrixtype_scatra_)
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      scatrastructuredomain_->Complete();
      scatrathermodomain_->Complete();
      structurescatradomain_->Complete();
      structurethermodomain_->Complete();
      thermoscatradomain_->Complete();
      thermostructuredomain_->Complete();

      if (interfacemeshtying_)
      {
        scatrastructureinterface_->Complete();
        scatrathermointerface_->Complete();
        thermoscatrainterface_->Complete();
        thermostructureinterface_->Complete();
      }
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      scatrastructuredomain_->Complete(*ssti_maps_mono_->BlockMapStructure()->FullMap(),
          *ssti_maps_mono_->BlockMapScatra()->FullMap());
      scatrathermodomain_->Complete(*ssti_maps_mono_->BlockMapThermo()->FullMap(),
          *ssti_maps_mono_->BlockMapScatra()->FullMap());
      structurescatradomain_->Complete(*ssti_maps_mono_->BlockMapScatra()->FullMap(),
          *ssti_maps_mono_->BlockMapStructure()->FullMap());
      structurethermodomain_->Complete(*ssti_maps_mono_->BlockMapThermo()->FullMap(),
          *ssti_maps_mono_->BlockMapStructure()->FullMap());
      thermoscatradomain_->Complete(*ssti_maps_mono_->BlockMapScatra()->FullMap(),
          *ssti_maps_mono_->BlockMapThermo()->FullMap());
      thermostructuredomain_->Complete(*ssti_maps_mono_->BlockMapStructure()->FullMap(),
          *ssti_maps_mono_->BlockMapThermo()->FullMap());

      if (interfacemeshtying_)
      {
        scatrastructureinterface_->Complete(*ssti_maps_mono_->BlockMapStructure()->FullMap(),
            *ssti_maps_mono_->BlockMapScatra()->FullMap());
        scatrathermointerface_->Complete(*ssti_maps_mono_->BlockMapThermo()->FullMap(),
            *ssti_maps_mono_->BlockMapScatra()->FullMap());
        thermoscatrainterface_->Complete(*ssti_maps_mono_->BlockMapScatra()->FullMap(),
            *ssti_maps_mono_->BlockMapThermo()->FullMap());
        thermostructureinterface_->Complete(*ssti_maps_mono_->BlockMapStructure()->FullMap(),
            *ssti_maps_mono_->BlockMapThermo()->FullMap());
      }
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSTI::SSTIMatrices::UnCompleteCouplingMatrices()
{
  scatrastructuredomain_->UnComplete();
  scatrathermodomain_->UnComplete();
  structurescatradomain_->UnComplete();
  structurethermodomain_->UnComplete();
  thermoscatradomain_->UnComplete();
  thermostructuredomain_->UnComplete();

  if (interfacemeshtying_)
  {
    scatrastructureinterface_->UnComplete();
    scatrathermointerface_->UnComplete();
    thermoscatrainterface_->UnComplete();
    thermostructureinterface_->UnComplete();
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> SSTI::SSTIMatrices::SetupBlockMatrix(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> row_map,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> col_map)
{
  const int expected_entries_per_row = 81;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
      *col_map, *row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> SSTI::SSTIMatrices::SetupSparseMatrix(
    const Teuchos::RCP<const Epetra_Map> row_map)
{
  const int expected_entries_per_row = 27;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSTI::ConvCheckMono::ConvCheckMono(Teuchos::ParameterList params)
    : itermax_(params.get<int>("ITEMAX")),
      itertol_(params.sublist("MONOLITHIC").get<double>("CONVTOL")),
      restol_(params.sublist("MONOLITHIC").get<double>("ABSTOLRES"))
{
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
bool SSTI::ConvCheckMono::Converged(const SSTI::SSTIMono& ssti_mono)
{
  bool exit(false);

  // compute L2 norm of concentration state vector
  double concdofnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssti_mono.ScaTraField()->Phinp())
      ->Norm2(&concdofnorm);

  // compute L2 norm of concentration increment vector
  double concincnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssti_mono.AllMaps()->MapsSubProblems()->ExtractVector(
          ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&concincnorm);

  // compute L2 norm of concentration residual vector
  double concresnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssti_mono.AllMaps()->MapsSubProblems()->ExtractVector(
          ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&concresnorm);

  // compute L2 norm of potential state vector
  double potdofnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssti_mono.ScaTraField()->Phinp())
      ->Norm2(&potdofnorm);

  // compute L2 norm of potential increment vector
  double potincnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssti_mono.AllMaps()->MapsSubProblems()->ExtractVector(
          ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&potincnorm);

  // compute L2 norm of potential residual vector
  double potresnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssti_mono.AllMaps()->MapsSubProblems()->ExtractVector(
          ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&potresnorm);

  // compute L2 norm of structural state vector
  double structuredofnorm(0.0);
  ssti_mono.StructureField()->Dispnp()->Norm2(&structuredofnorm);

  // compute L2 norm of structural residual vector
  double structureresnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubProblems()
      ->ExtractVector(ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureresnorm);

  // compute L2 norm of structural increment vector
  double structureincnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubProblems()
      ->ExtractVector(ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureincnorm);

  // compute L2 norm of thermo state vector
  double thermodofnorm(0.0);
  ssti_mono.ThermoField()->Phinp()->Norm2(&thermodofnorm);

  // compute L2 norm of thermo residual vector
  double thermoresnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubProblems()
      ->ExtractVector(ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::thermo))
      ->Norm2(&thermoresnorm);

  // compute L2 norm of thermo increment vector
  double thermoincnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubProblems()
      ->ExtractVector(ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::thermo))
      ->Norm2(&thermoincnorm);

  // compute L2 norm of total residual vector
  double totresnorm(0.0);
  ssti_mono.Residual()->Norm2(&totresnorm);

  // safety checks
  if (std::isnan(concdofnorm) or std::isnan(concresnorm) or std::isnan(concincnorm) or
      std::isnan(potdofnorm) or std::isnan(potresnorm) or std::isnan(potincnorm) or
      std::isnan(structuredofnorm) or std::isnan(structureresnorm) or
      std::isnan(structureincnorm) or std::isnan(thermodofnorm) or std::isnan(thermoresnorm) or
      std::isnan(thermoincnorm))
    dserror("Vector norm is not a number!");
  if (std::isinf(concdofnorm) or std::isinf(concresnorm) or std::isinf(concincnorm) or
      std::isinf(potdofnorm) or std::isinf(potresnorm) or std::isinf(potincnorm) or
      std::isinf(structuredofnorm) or std::isinf(structureresnorm) or
      std::isinf(structureincnorm) or std::isnan(thermodofnorm) or std::isnan(thermoresnorm) or
      std::isnan(thermoincnorm))
    dserror("Vector norm is infinity!");

  // prevent division by zero
  if (concdofnorm < 1.e-10) concdofnorm = 1.e-10;
  if (potdofnorm < 1.e-10) potdofnorm = 1.e-10;
  if (structuredofnorm < 1.e-10) structuredofnorm = 1.e-10;
  if (thermodofnorm < 1.e-10) thermodofnorm = 1.e-10;

  // first Newton-Raphson iteration
  if (ssti_mono.NewtonIteration() == 1)
  {
    if (ssti_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+--------------+--------------+----"
                   "----------+"
                << std::endl;
      std::cout << "|- step/max -|- tolerance[norm] -|-- conc-res --|-- conc-inc --|-- pot-res "
                   "---|-- pot-inc ---|- struct-res -|- struct-inc -|- thermo-res -|- thermo-inc "
                   "-|-  tot. res  -|"
                << std::endl;

      // print first line of convergence table to screen
      // solution increment not yet available during first Newton-Raphson iteration
      std::cout << "|  " << std::setw(3) << ssti_mono.NewtonIteration() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << potresnorm << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureresnorm
                << "   |      --      | " << std::setw(10) << std::setprecision(3)
                << std::scientific << thermoresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << totresnorm << "   |    | "
                << "(       --      , te = " << std::setw(10) << std::setprecision(3)
                << ssti_mono.TimeStatistics()[0] << ")" << std::endl;
    }
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if (ssti_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << ssti_mono.NewtonIteration() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << concincnorm / concdofnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << potresnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific
                << potincnorm / potdofnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << structureresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << structureincnorm / structuredofnorm
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << thermoresnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << thermoincnorm / thermodofnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << totresnorm << "   | "
                << "   | (ts = " << std::setw(10) << std::setprecision(3)
                << ssti_mono.TimeStatistics()[1] << ", te = " << std::setw(10)
                << std::setprecision(3) << ssti_mono.TimeStatistics()[0] << ")" << std::endl;
    }

    // convergence check
    if (concresnorm <= itertol_ and potresnorm <= itertol_ and structureresnorm <= itertol_ and
        thermoresnorm <= itertol_ and concincnorm / concdofnorm <= itertol_ and
        potincnorm / potdofnorm <= itertol_ and structureincnorm / structuredofnorm <= itertol_ and
        thermoincnorm / thermodofnorm <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (concresnorm < restol_ and potresnorm < restol_ and structureresnorm < restol_ and
      thermoresnorm < restol_)
    exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence
  if (ssti_mono.NewtonIteration() == itermax_ and !exit)
  {
    if (ssti_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+--------------+--------------+----"
                   "----------+"
                << std::endl;
      std::cout << "|                     Newton-Raphson method has not converged after a maximum "
                   "number of "
                << std::setw(2) << itermax_ << " iterations!                     |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  return exit;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
std::map<std::string, std::string> SSTI::SSTIScatraStructureCloneStrategy::ConditionsToCopy() const
{
  // call base class
  std::map<std::string, std::string> conditions_to_copy =
      SSI::ScatraStructureCloneStrategy::ConditionsToCopy();

  conditions_to_copy.insert({"ThermoDirichlet", "ThermoDirichlet"});
  conditions_to_copy.insert({"ThermoPointNeumann", "ThermoPointNeumann"});
  conditions_to_copy.insert({"ThermoLineNeumann", "ThermoLineNeumann"});
  conditions_to_copy.insert({"ThermoSurfaceNeumann", "ThermoSurfaceNeumann"});
  conditions_to_copy.insert({"ThermoVolumeNeumann", "ThermoVolumeNeumann"});
  conditions_to_copy.insert({"ThermoInitfield", "ThermoInitfield"});
  conditions_to_copy.insert({"SSTIMeshtying3DomainIntersection", "Meshtying3DomainIntersection"});
  conditions_to_copy.insert({"SSTIInterfaceMeshtying", "S2IMeshtying"});

  return conditions_to_copy;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSTI::SSTIScatraStructureCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  auto* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    // set distype as well!
    trans->SetDisType(oldele->Shape());

    // now check whether ImplType is reasonable and if set the ImplType
    INPAR::SCATRA::ImplType impltype = SSI::ScatraStructureCloneStrategy::GetImplType(oldele);

    if (impltype == INPAR::SCATRA::impltype_undefined)
    {
      dserror(
          "ScatraStructureCloneStrategy copies scatra discretization from structure "
          "discretization, but the STRUCTURE elements that are defined in the .dat file are either "
          "not meant to be copied to scatra elements or the ImplType is set 'Undefined' which is "
          "not meaningful for the created scatra discretization! Use SOLIDSCATRA, WALLSCATRA or "
          "SHELLSCATRA elements with meaningful ImplType instead!");
    }
    else
    {
      // find the appropriate thermo type
      if (impltype == INPAR::SCATRA::impltype_elch_electrode)
        trans->SetImplType(INPAR::SCATRA::impltype_elch_electrode_thermo);
      else if (impltype == INPAR::SCATRA::impltype_elch_diffcond)
        trans->SetImplType(INPAR::SCATRA::impltype_elch_diffcond_thermo);
      else
        dserror("Something went wrong");
    }

    // set material
    trans->SetMaterial(matid, oldele);
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
std::map<std::string, std::string> SSTI::SSTIScatraThermoCloneStrategy::ConditionsToCopy() const
{
  // call base class
  std::map<std::string, std::string> conditions_to_copy =
      STI::ScatraThermoCloneStrategy::ConditionsToCopy();

  conditions_to_copy.insert({"Meshtying3DomainIntersection", "Meshtying3DomainIntersection"});
  conditions_to_copy.insert({"SSTIInterfaceMeshtying", "S2IMeshtying"});
  conditions_to_copy.insert({"TotalAndMeanScalar", "TotalAndMeanScalar"});

  return conditions_to_copy;
}
BACI_NAMESPACE_CLOSE
