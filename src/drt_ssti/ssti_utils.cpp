/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSTI

 \level 2


 *------------------------------------------------------------------------------------------------*/


#include <Epetra_Map.h>

#include "ssti_utils.H"

#include "ssti_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSTI::SSTIMaps::SSTIMaps(const SSTI::SSTIMono& ssti_mono_algorithm)
    : map_structure_condensed_(Teuchos::null),
      maps_interface_structure_(Teuchos::null),
      maps_scatra_(Teuchos::null),
      maps_structure_(Teuchos::null),
      maps_subproblems_(Teuchos::null),
      maps_thermo_(Teuchos::null)
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
      LINALG::MergeMap(partial_maps[0], partial_maps[1], false);
  Teuchos::RCP<const Epetra_Map> merged_map = LINALG::MergeMap(temp_map, partial_maps[2], false);
  // initialize global map extractor
  maps_subproblems_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*merged_map, partial_maps));
  // check global map extractor
  maps_subproblems_->CheckForValidMapExtractor();

  // initialize map extractors associated with blocks of subproblems
  maps_structure_ =
      Teuchos::rcp(new LINALG::MultiMapExtractor(*ssti_mono_algorithm.StructureField()->DofRowMap(),
          std::vector<Teuchos::RCP<const Epetra_Map>>(
              1, ssti_mono_algorithm.StructureField()->DofRowMap())));
  switch (ssti_mono_algorithm.ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::sparse:
    {
      maps_scatra_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*ssti_mono_algorithm.ScaTraField()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(
                  1, ssti_mono_algorithm.ScaTraField()->DofRowMap())));
      maps_thermo_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*ssti_mono_algorithm.ThermoField()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(
                  1, ssti_mono_algorithm.ThermoField()->DofRowMap())));
      break;
    }
    case LINALG::MatrixType::block_condition:
    {
      maps_scatra_ = Teuchos::rcpFromRef(ssti_mono_algorithm.ScaTraField()->BlockMaps());
      maps_thermo_ = Teuchos::rcpFromRef(ssti_mono_algorithm.ThermoField()->BlockMaps());
      break;
    }
    default:
    {
      dserror("Matrix type not supported");
      break;
    }
  }

  maps_scatra_->CheckForValidMapExtractor();
  maps_structure_->CheckForValidMapExtractor();
  maps_thermo_->CheckForValidMapExtractor();

  if (ssti_mono_algorithm.InterfaceMeshtying())
  {
    // set up map for interior and master-side structural degrees of freedom
    map_structure_condensed_ = LINALG::SplitMap(*ssti_mono_algorithm.StructureField()->DofRowMap(),
        *ssti_mono_algorithm.CouplingAdapterStructure()->SlaveDofMap());

    // set up structural map extractor holding interface maps of dofs
    std::vector<Teuchos::RCP<const Epetra_Map>> maps_interface(0, Teuchos::null);
    maps_interface.emplace_back(ssti_mono_algorithm.CouplingAdapterStructure()->SlaveDofMap());
    maps_interface.emplace_back(ssti_mono_algorithm.CouplingAdapterStructure()->MasterDofMap());
    maps_interface.emplace_back(LINALG::SplitMap(*map_structure_condensed_,
        *ssti_mono_algorithm.CouplingAdapterStructure()->MasterDofMap()));
    maps_interface_structure_ = Teuchos::rcp(new LINALG::MultiMapExtractor(
        *ssti_mono_algorithm.StructureField()->DofRowMap(), maps_interface));
    maps_interface_structure_->CheckForValidMapExtractor();
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Map> SSTI::SSTIMaps::MapInterface(
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy) const
{
  // only relevant in case of meshtying
  if (meshtyingstrategy == Teuchos::null)
    return Teuchos::null;
  else
  {
    auto mergedInterfaceMap =
        LINALG::MultiMapExtractor::MergeMaps({meshtyingstrategy->CouplingAdapter()->MasterDofMap(),
            meshtyingstrategy->CouplingAdapter()->SlaveDofMap()});
    if (not mergedInterfaceMap->UniqueGIDs()) dserror("Map not unique");
    return mergedInterfaceMap;
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
const Teuchos::RCP<LINALG::MultiMapExtractor> SSTI::SSTIMaps::MapsInterfaceBlocks(
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
    LINALG::MatrixType scatramatrixtype, unsigned nummaps) const
{
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapinterface(Teuchos::null);

  Teuchos::RCP<Epetra_Map> interfacemap = MapInterface(meshtyingstrategy);

  switch (scatramatrixtype)
  {
    case LINALG::MatrixType::sparse:
    {
      blockmapinterface = Teuchos::rcp(new LINALG::MultiMapExtractor(
          *interfacemap, std::vector<Teuchos::RCP<const Epetra_Map>>(1, interfacemap)));
      break;
    }
    case LINALG::MatrixType::block_condition:
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> partial_blockmapinterface(nummaps, Teuchos::null);
      for (int iblockmap = 0; iblockmap < static_cast<int>(nummaps); ++iblockmap)
      {
        partial_blockmapinterface[iblockmap] = LINALG::MultiMapExtractor::MergeMaps(
            {meshtyingstrategy->BlockMapsSlave().Map(iblockmap),
                meshtyingstrategy->BlockMapsMaster().Map(iblockmap)});
      }
      blockmapinterface =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*interfacemap, partial_blockmapinterface));
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
const Teuchos::RCP<LINALG::MultiMapExtractor> SSTI::SSTIMaps::MapsInterfaceBlocksSlave(
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
    LINALG::MatrixType scatramatrixtype, unsigned nummaps) const
{
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapinterfaceslave(Teuchos::null);

  switch (scatramatrixtype)
  {
    case LINALG::MatrixType::sparse:
    {
      const auto slavedofmap = meshtyingstrategy->CouplingAdapter()->SlaveDofMap();
      blockmapinterfaceslave = Teuchos::rcp(new LINALG::MultiMapExtractor(
          *slavedofmap, std::vector<Teuchos::RCP<const Epetra_Map>>(1, slavedofmap)));
      break;
    }
    case LINALG::MatrixType::block_condition:
    {
      blockmapinterfaceslave =
          Teuchos::rcp(new LINALG::MultiMapExtractor(meshtyingstrategy->BlockMapsSlave()));
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
    : SSTIMaps(ssti_mono_algorithm), maps_systemmatrix_subblocks_(Teuchos::null)
{
  // initialize map extractors associated with blocks of global system matrix
  switch (ssti_mono_algorithm.ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case LINALG::MatrixType::sparse:
    {
      maps_systemmatrix_subblocks_ = MapsSubproblems();
      break;
    }
      // many main-diagonal matrix blocks associated with scalar transport field
    case LINALG::MatrixType::block_condition:
    {
      Teuchos::RCP<std::vector<int>> block_positions_scatra =
          ssti_mono_algorithm.GetBlockPositions(Subproblem::scalar_transport);
      Teuchos::RCP<std::vector<int>> block_positions_structure =
          ssti_mono_algorithm.GetBlockPositions(Subproblem::structure);
      Teuchos::RCP<std::vector<int>> block_positions_thermo =
          ssti_mono_algorithm.GetBlockPositions(Subproblem::thermo);

      std::vector<Teuchos::RCP<const Epetra_Map>> maps_systemmatrix(
          block_positions_scatra->size() + block_positions_structure->size() +
          block_positions_thermo->size());
      for (int imap = 0; imap < static_cast<int>(block_positions_scatra->size()); ++imap)
        maps_systemmatrix[block_positions_scatra->at(imap)] = MapsScatra()->Map(imap);

      // extract map underlying single main-diagonal matrix block associated with structural
      // field
      maps_systemmatrix[block_positions_structure->at(0)] = MapsStructure()->FullMap();

      for (int imap = 0; imap < static_cast<int>(block_positions_thermo->size()); ++imap)
        maps_systemmatrix[block_positions_thermo->at(imap)] = MapsThermo()->Map(imap);

      // initialize map extractor associated with blocks of global system matrix
      maps_systemmatrix_subblocks_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*MapsSubproblems()->FullMap(), maps_systemmatrix));

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
    const LINALG::MatrixType matrixtype_global, const LINALG::MatrixType matrixtype_scatra,
    Teuchos::RCP<Epetra_Map> interface_map_scatra, Teuchos::RCP<Epetra_Map> interface_map_thermo,
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface,
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface, bool isinterfacemeshtying)
    : systemmatrix_(Teuchos::null),
      scatrastructuredomain_(Teuchos::null),
      scatrastructureinterface_(Teuchos::null),
      scatrathermodomain_(Teuchos::null),
      scatrathermointerface_(Teuchos::null),
      structurescatradomain_(Teuchos::null),
      structurethermodomain_(Teuchos::null),
      thermoscatradomain_(Teuchos::null),
      thermoscatrainterface_(Teuchos::null),
      thermostructuredomain_(Teuchos::null),
      thermostructureinterface_(Teuchos::null)
{
  // perform initializations associated with global system matrix
  switch (matrixtype_global)
  {
    case LINALG::MatrixType::block_field:
    {
      systemmatrix_ = SetupBlockMatrix(
          ssti_maps_mono->MapsSystemMatrixSubblocks(), ssti_maps_mono->MapsSystemMatrixSubblocks());
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      systemmatrix_ = SetupSparseMatrix(ssti_maps_mono->MapsSubproblems()->FullMap());
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
    case LINALG::MatrixType::block_condition:
    {
      scatrastructuredomain_ =
          SetupBlockMatrix(ssti_maps_mono->MapsScatra(), ssti_maps_mono->MapsStructure());
      structurescatradomain_ =
          SetupBlockMatrix(ssti_maps_mono->MapsStructure(), ssti_maps_mono->MapsScatra());
      structurethermodomain_ =
          SetupBlockMatrix(ssti_maps_mono->MapsStructure(), ssti_maps_mono->MapsThermo());
      thermostructuredomain_ =
          SetupBlockMatrix(ssti_maps_mono->MapsThermo(), ssti_maps_mono->MapsStructure());
      scatrathermodomain_ =
          SetupBlockMatrix(ssti_maps_mono->MapsScatra(), ssti_maps_mono->MapsThermo());
      thermoscatradomain_ =
          SetupBlockMatrix(ssti_maps_mono->MapsThermo(), ssti_maps_mono->MapsScatra());

      if (isinterfacemeshtying)
      {
        scatrastructureinterface_ =
            SetupBlockMatrix(blockmapscatrainterface, ssti_maps_mono->MapsStructure());
        thermostructureinterface_ =
            SetupBlockMatrix(blockmapthermointerface, ssti_maps_mono->MapsStructure());
        scatrathermointerface_ =
            SetupBlockMatrix(blockmapscatrainterface, ssti_maps_mono->MapsThermo());
        thermoscatrainterface_ =
            SetupBlockMatrix(blockmapthermointerface, ssti_maps_mono->MapsScatra());
      }
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      scatrastructuredomain_ = SetupSparseMatrix(ssti_maps_mono->MapsScatra()->FullMap());
      structurescatradomain_ = SetupSparseMatrix(ssti_maps_mono->MapsStructure()->FullMap());
      structurethermodomain_ = SetupSparseMatrix(ssti_maps_mono->MapsStructure()->FullMap());
      thermostructuredomain_ = SetupSparseMatrix(ssti_maps_mono->MapsThermo()->FullMap());
      scatrathermodomain_ = SetupSparseMatrix(ssti_maps_mono->MapsScatra()->FullMap());
      thermoscatradomain_ = SetupSparseMatrix(ssti_maps_mono->MapsThermo()->FullMap());

      if (isinterfacemeshtying)
      {
        scatrastructureinterface_ = SetupSparseMatrix(interface_map_scatra);
        thermostructureinterface_ = SetupSparseMatrix(interface_map_thermo);
        scatrathermointerface_ = SetupSparseMatrix(interface_map_scatra);
        thermoscatrainterface_ = SetupSparseMatrix(interface_map_thermo);
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
Teuchos::RCP<LINALG::BlockSparseMatrixBase> SSTI::SSTIMatrices::SetupBlockMatrix(
    Teuchos::RCP<const LINALG::MultiMapExtractor> row_map,
    Teuchos::RCP<const LINALG::MultiMapExtractor> col_map)
{
  const int expected_entries_per_row = 81;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *col_map, *row_map, expected_entries_per_row, explicitdirichlet, savegraph));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> SSTI::SSTIMatrices::SetupSparseMatrix(
    const Teuchos::RCP<const Epetra_Map> row_map)
{
  const int expected_entries_per_row = 27;
  const bool explicitdirichlet = false;
  const bool savegraph = true;

  return Teuchos::rcp(
      new LINALG::SparseMatrix(*row_map, expected_entries_per_row, explicitdirichlet, savegraph));
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
      ->ExtractOtherVector(ssti_mono.AllMaps()->MapsSubproblems()->ExtractVector(
          ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&concincnorm);

  // compute L2 norm of concentration residual vector
  double concresnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssti_mono.AllMaps()->MapsSubproblems()->ExtractVector(
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
      ->ExtractCondVector(ssti_mono.AllMaps()->MapsSubproblems()->ExtractVector(
          ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&potincnorm);

  // compute L2 norm of potential residual vector
  double potresnorm(0.0);
  ssti_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssti_mono.AllMaps()->MapsSubproblems()->ExtractVector(
          ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&potresnorm);

  // compute L2 norm of structural state vector
  double structuredofnorm(0.0);
  ssti_mono.StructureField()->Dispnp()->Norm2(&structuredofnorm);

  // compute L2 norm of structural residual vector
  double structureresnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubproblems()
      ->ExtractVector(ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureresnorm);

  // compute L2 norm of structural increment vector
  double structureincnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubproblems()
      ->ExtractVector(ssti_mono.Increment(), ssti_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&structureincnorm);

  // compute L2 norm of thermo state vector
  double thermodofnorm(0.0);
  ssti_mono.ThermoField()->Phinp()->Norm2(&thermodofnorm);

  // compute L2 norm of thermo residual vector
  double thermoresnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubproblems()
      ->ExtractVector(ssti_mono.Residual(), ssti_mono.GetProblemPosition(Subproblem::thermo))
      ->Norm2(&thermoresnorm);

  // compute L2 norm of thermo increment vector
  double thermoincnorm(0.0);
  ssti_mono.AllMaps()
      ->MapsSubproblems()
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
  if (ssti_mono.NewtonIteration() == itermax_)
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
std::map<std::string, std::string> SSTI::SSTIScatraStructureCloneStrategy::ConditionsToCopy()
{
  // call base class
  std::map<std::string, std::string> conditions_to_copy =
      SSI::ScatraStructureCloneStrategy::ConditionsToCopy();

  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ThermoDirichlet", "ThermoDirichlet"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ThermoPointNeumann", "ThermoPointNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ThermoLineNeumann", "ThermoLineNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ThermoSurfaceNeumann", "ThermoSurfaceNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ThermoVolumeNeumann", "ThermoVolumeNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ThermoInitfield", "ThermoInitfield"));

  return conditions_to_copy;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSTI::PrintSSTILogo(int pid)
{
  if (pid == 0)
  {
    std::cout << "    ██████   ██████ ▄▄▄█████▓ ██▓" << std::endl;
    std::cout << "  ▒██    ▒ ▒██    ▒ ▓  ██▒ ▓▒▓██▒" << std::endl;
    std::cout << "  ░ ▓██▄   ░ ▓██▄   ▒ ▓██░ ▒░▒██▒" << std::endl;
    std::cout << "    ▒   ██▒  ▒   ██▒░ ▓██▓ ░ ░██░" << std::endl;
    std::cout << "  ▒██████▒▒▒██████▒▒  ▒██▒ ░ ░██░" << std::endl;
    std::cout << "  ▒ ▒▓▒ ▒ ░▒ ▒▓▒ ▒ ░  ▒ ░░   ░▓" << std::endl;
    std::cout << "  ░ ░▒  ░ ░░ ░▒  ░ ░    ░     ▒ ░" << std::endl;
    std::cout << "  ░  ░  ░  ░  ░  ░    ░       ▒ ░" << std::endl;
    std::cout << "        ░        ░            ░" << std::endl;
  }
}
