/*----------------------------------------------------------------------*/
/*! \file

\brief monolithic coupling algorithm for scatra-thermo interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#include "sti_monolithic.H"
#include "sti_monolithic_evaluate_OffDiag.H"
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_coupling.H"
#include "../linalg/linalg_matrixtransform.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_utils_densematrix_communication.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

/*--------------------------------------------------------------------------------*
 | constructor                                                         fang 09/17 |
 *--------------------------------------------------------------------------------*/
STI::Monolithic::Monolithic(const Epetra_Comm& comm,  //! communicator
    const Teuchos::ParameterList& stidyn,  //! parameter list for scatra-thermo interaction
    const Teuchos::ParameterList&
        scatradyn,  //! scalar transport parameter list for scatra and thermo fields
    const Teuchos::ParameterList&
        solverparams,  //! solver parameter list for scatra-thermo interaction
    const Teuchos::ParameterList& solverparams_scatra,  //! solver parameter list for scatra field
    const Teuchos::ParameterList& solverparams_thermo   //! solver parameter list for thermo field
    )
    :  // instantiate base class
      Algorithm(comm, stidyn, scatradyn, solverparams_scatra, solverparams_thermo),

      restol_(fieldparameters_->sublist("NONLINEAR").get<double>("ABSTOLRES")),
      maps_(Teuchos::null),
      condensationthermo_(DRT::INPUT::IntegralValue<bool>(stidyn, "THERMO_CONDENSATION")),
      systemmatrix_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<LINALG::MatrixType>(
          stidyn.sublist("MONOLITHIC"), "MATRIXTYPE")),
      scatrathermoblockdomain_(Teuchos::null),
      scatrathermoblockinterface_(Teuchos::null),
      thermoscatrablockdomain_(Teuchos::null),
      thermoscatrablockinterface_(Teuchos::null),
      blockmaps_(Teuchos::null),
      blockmapthermo_(Teuchos::null),
      increment_(Teuchos::null),
      residual_(Teuchos::null),
      dtele_(0.),
      dtsolve_(0.),

      // initialize algebraic solver for global system of equations
      solver_(Teuchos::rcp(
          new LINALG::Solver(solverparams, comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      invrowsums_(Teuchos::null),

      icoupscatra_(Teuchos::null),
      icoupthermo_(Teuchos::null),
      islavetomasterrowtransformscatraod_(Teuchos::null),
      islavetomastercoltransformthermood_(Teuchos::null),
      islavetomasterrowtransformthermood_(Teuchos::null),
      equilibration_(Teuchos::null)
{
  // safety checks
  if (!ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for scatra-thermo interaction!");
  if (ThermoField()->SystemMatrix() == Teuchos::null)
    dserror("System matrix associated with temperature field must be a sparse matrix!");

  // set control parameters for Newton-Raphson iteration loop
  itermax_ = fieldparameters_->sublist("NONLINEAR").get<int>("ITEMAX");
  itertol_ = fieldparameters_->sublist("NONLINEAR").get<double>("CONVTOL");


  // extract coupling adapters for scatra-scatra interface coupling
  if (ScaTraField()->S2ICoupling() and
      strategyscatra_->CouplingType() == INPAR::S2I::coupling_matching_nodes)
  {
    icoupscatra_ = strategyscatra_->CouplingAdapter();
    icoupthermo_ = strategythermo_->CouplingAdapter();
  }

  // initialize map associated with single thermo block of global system matrix
  Teuchos::RCP<const Epetra_Map> mapthermo(Teuchos::null);
  if (condensationthermo_)
    mapthermo = LINALG::MergeMap(
        *strategythermo_->InterfaceMaps()->Map(0), *strategythermo_->InterfaceMaps()->Map(2));
  else
    mapthermo = ThermoField()->DofRowMap();

  // initialize global map extractor
  maps_ = Teuchos::rcp(new LINALG::MapExtractor(
      *LINALG::MergeMap(*ScaTraField()->Discretization()->DofRowMap(), *mapthermo, false),
      mapthermo, ScaTraField()->DofRowMap()));

  // check global map extractor
  maps_->CheckForValidMapExtractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize global residual vector
  residual_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize transformation operators
  islavetomasterrowtransformscatraod_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
  islavetomastercoltransformthermood_ = Teuchos::rcp(new LINALG::MatrixColTransform);
  islavetomasterrowtransformthermood_ = Teuchos::rcp(new LINALG::MatrixRowTransform);

  // merge slave and master side block maps for interface matrix for thermo and scatra
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);
  Teuchos::RCP<Epetra_Map> interface_map_thermo(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface(Teuchos::null);
  if (ScaTraField()->S2ICoupling())
  {
    // merge slave and master side full maps for interface matrix for thermo and scatra
    interface_map_scatra = LINALG::MultiMapExtractor::MergeMaps(
        {strategyscatra_->InterfaceMaps()->Map(1), strategyscatra_->InterfaceMaps()->Map(2)});
    interface_map_thermo = LINALG::MultiMapExtractor::MergeMaps(
        {strategythermo_->InterfaceMaps()->Map(1), strategythermo_->InterfaceMaps()->Map(2)});
    // build block map for thermo interface by using full thermo map
    blockmapthermointerface = Teuchos::rcp(new LINALG::MultiMapExtractor(*interface_map_thermo,
        std::vector<Teuchos::RCP<const Epetra_Map>>(1, interface_map_thermo)));
    blockmapthermointerface->CheckForValidMapExtractor();
  }

  // initialize map extractors associated with blocks of global system matrix
  switch (ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case LINALG::MatrixType::sparse:
    {
      blockmaps_ = maps_;

      if (ScaTraField()->S2ICoupling())
      {
        blockmapscatrainterface = Teuchos::rcp(new LINALG::MultiMapExtractor(*interface_map_scatra,
            std::vector<Teuchos::RCP<const Epetra_Map>>(1, interface_map_scatra)));
        blockmapscatrainterface->CheckForValidMapExtractor();
      }
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case LINALG::MatrixType::block_condition:
    {
      // extract maps underlying main-diagonal matrix blocks associated with scalar transport field
      const unsigned nblockmapsscatra = ScaTraField()->BlockMaps().NumMaps();
      std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps(nblockmapsscatra + 1);
      for (unsigned iblockmap = 0; iblockmap < nblockmapsscatra; ++iblockmap)
        blockmaps[iblockmap] = ScaTraField()->BlockMaps().Map(iblockmap);

      // extract map underlying single main-diagonal matrix block associated with temperature field
      blockmaps[nblockmapsscatra] = mapthermo;

      // initialize map extractor associated with blocks of global system matrix
      blockmaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*DofRowMap(), blockmaps));

      // initialize map extractor associated with all degrees of freedom inside temperature field
      blockmapthermo_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*ThermoField()->Discretization()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(1, ThermoField()->DofRowMap())));

      // safety check
      blockmapthermo_->CheckForValidMapExtractor();
      if (ScaTraField()->S2ICoupling())
      {
        // build block map for scatra interface by merging slave and master side for each block
        std::vector<Teuchos::RCP<const Epetra_Map>> partial_blockmapscatrainterface(
            nblockmapsscatra, Teuchos::null);
        for (unsigned iblockmap = 0; iblockmap < nblockmapsscatra; ++iblockmap)
          partial_blockmapscatrainterface.at(iblockmap) = LINALG::MultiMapExtractor::MergeMaps(
              {strategyscatra_->BlockMapsSlave().Map(iblockmap),
                  strategyscatra_->BlockMapsMaster().Map(iblockmap)});
        blockmapscatrainterface = Teuchos::rcp(
            new LINALG::MultiMapExtractor(*interface_map_scatra, partial_blockmapscatrainterface));
        blockmapscatrainterface->CheckForValidMapExtractor();
      }

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // safety check
  blockmaps_->CheckForValidMapExtractor();

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case LINALG::MatrixType::block_condition:
    {
      // safety check
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize global system matrix
      systemmatrix_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *blockmaps_, *blockmaps_, 81, false, true));

      // feed AMGnxn block preconditioner with null space information for each block of global block
      // system matrix
      BuildNullSpaces();

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // safety check
      if (ScaTraField()->MatrixType() != LINALG::MatrixType::sparse)
        dserror("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(), 27, false, true));

      // feed AMG preconditioner with null space information associated with global system matrix if
      // applicable
      ComputeNullSpaceIfNecessary(solver_->Params());

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scatra-thermo interaction not recognized!");
      break;
    }
  }

  // initialize scatra-thermo block and thermo-scatra block of global system matrix
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      // initialize scatra-thermo blocks
      scatrathermoblockdomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, ScaTraField()->BlockMaps(), 81, false, true));
      scatrathermoblockinterface_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, *blockmapscatrainterface, 81, false, true));

      // initialize thermo-scatra blocks
      thermoscatrablockdomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ScaTraField()->BlockMaps(), *blockmapthermo_, 81, false, true));
      thermoscatrablockinterface_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ScaTraField()->BlockMaps(), *blockmapthermointerface, 81, false, true));

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // initialize scatra-thermo blocks
      scatrathermoblockdomain_ = Teuchos::rcp(
          new LINALG::SparseMatrix(*ScaTraField()->Discretization()->DofRowMap(), 27, false, true));
      scatrathermoblockinterface_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_scatra, 27, false, true));

      // initialize thermo-scatra blocks
      thermoscatrablockdomain_ = Teuchos::rcp(
          new LINALG::SparseMatrix(*ThermoField()->Discretization()->DofRowMap(), 27, false, true));
      thermoscatrablockinterface_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_thermo, 27, false, true));

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // initialize OD evaluation strategy
  scatrathermooffdiagcoupling_ =
      STI::BuildScatraThermoOffDiagCoupling(strategyscatra_->CouplingType(), blockmapthermo_,
          blockmapthermointerface, maps_->Map(0), maps_->Map(1), interface_map_scatra,
          interface_map_thermo, strategyscatra_, strategythermo_, scatra_, thermo_);

  // instantiate appropriate equilibration class
  equilibration_ = LINALG::BuildEquilibration(
      matrixtype_, ScaTraField()->EquilibrationMethod(), maps_->FullMap());
}

/*---------------------------------------------------------------------------------------------*
 | finite difference check for global system matrix (for debugging only)            fang 07/15 |
 *---------------------------------------------------------------------------------------------*/
void STI::Monolithic::FDCheck()
{
  // initial screen output
  if (Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR STI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(), true));
  maps_->InsertVector(ScaTraField()->Phinp(), 0, statenp);
  maps_->InsertVector(ThermoField()->Phinp(), 1, statenp);

  // make a copy of global state vector to undo perturbations later
  auto statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if (Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_) != Teuchos::null)
    sysmat_original =
        (new LINALG::SparseMatrix(
             *(Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_)->Merge())))
            ->EpetraMatrix();
  else
    dserror("Global system matrix must be a block sparse matrix!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  auto rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid = 0; colgid <= sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // check whether current column index is a valid global column index and continue loop if not
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    Comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // fill global state vector with original state variables
    statenp->Update(1., *statenp_original, 0.);

    // impose perturbation
    if (statenp->Map().MyGID(colgid))
      if (statenp->SumIntoGlobalValue(colgid, 0, ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    ScaTraField()->Phinp()->Update(1., *maps_->ExtractVector(statenp, 0), 0.);
    ThermoField()->Phinp()->Update(1., *maps_->ExtractVector(statenp, 1), 0.);

    // carry perturbation over to state vectors at intermediate time stages if necessary
    ScaTraField()->ComputeIntermediateValues();
    ThermoField()->ComputeIntermediateValues();

    // calculate element right-hand side vector for perturbed state
    AssembleMatAndRHS();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the
    // system matrix, the second comparison might yield good agreement in spite of the entries being
    // wrong!
    for (int rowlid = 0; rowlid < DofRowMap()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if (rowgid < 0) dserror("Invalid global ID of matrix row!");

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid, length, numentries, &values[0], &indices[0]);
      for (int ientry = 0; ientry < length; ++ientry)
      {
        if (sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval = -(*residual_)[rowlid] / ScaTraField()->FDCheckEps() +
                           (*rhs_original)[rowlid] / ScaTraField()->FDCheckEps();

      // confirm accuracy of first comparison
      if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        dserror("Finite difference check involves values too close to numerical zero!");

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if (abs(abserr1) > maxabserr) maxabserr = abs(abserr1);
      double relerr1(0.);
      if (abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if (abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if (abs(relerr1) > maxrelerr) maxrelerr = abs(relerr1);

      // evaluate first comparison
      if (abs(relerr1) > ScaTraField()->FDCheckTol())
      {
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left = entry - (*rhs_original)[rowlid] / ScaTraField()->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / ScaTraField()->FDCheckEps();

        // confirm accuracy of second comparison
        if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if (abs(abserr2) > maxabserr) maxabserr = abs(abserr2);
        double relerr2(0.);
        if (abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if (abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if (abs(relerr2) > maxrelerr) maxrelerr = abs(relerr2);

        // evaluate second comparison
        if (abs(relerr2) > ScaTraField()->FDCheckTol())
        {
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid
                    << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  Comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  Comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  Comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (Comm().MyPID() == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      dserror("Finite difference check failed for STI system matrix!");
    }
    else
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
  }

  // undo perturbations of state variables
  ScaTraField()->Phinp()->Update(1., *maps_->ExtractVector(statenp_original, 0), 0.);
  ScaTraField()->ComputeIntermediateValues();
  ThermoField()->Phinp()->Update(1., *maps_->ExtractVector(statenp_original, 1), 0.);
  ThermoField()->ComputeIntermediateValues();

  // recompute system matrix and right-hand side vector based on original state variables
  AssembleMatAndRHS();

  return;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------------------------*
 | output matrix to *.csv file for debugging purposes, with global row and column IDs of matrix
 components in ascending order across all processors   fang 01/17 |
 *---------------------------------------------------------------------------------------------------------------------------------------------------------------*/
void STI::Monolithic::OutputMatrixToFile(
    const Teuchos::RCP<const LINALG::SparseOperator>
        sparseoperator,     //!< sparse or block sparse matrix to be output
    const int precision,    //!< output precision
    const double tolerance  //!< output omission tolerance
)
{
  // safety check
  if (!sparseoperator->Filled()) dserror("Sparse operator must be filled for output!");

  // extract communicator
  const Epetra_Comm& comm = sparseoperator->Comm();

  // determine whether sparse matrix or block sparse matrix should be output
  const auto sparsematrix = Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(sparseoperator);
  const auto blocksparsematrix =
      Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(sparseoperator);
  if (sparsematrix == Teuchos::null and blocksparsematrix == Teuchos::null)
    dserror("Unknown type of sparse operator!");

  // extract row map
  const Epetra_Map& rowmap =
      sparsematrix != Teuchos::null ? sparsematrix->RowMap() : blocksparsematrix->FullRowMap();

  // safety check
  if (!rowmap.UniqueGIDs()) dserror("Row map of matrix must be non-overlapping!");

  // copy global IDs of matrix rows stored on current processor into vector
  std::vector<int> myrowgids(rowmap.NumMyElements(), 0);
  int* myglobalelements = rowmap.MyGlobalElements();
  std::copy(myglobalelements, myglobalelements + rowmap.NumMyElements(), &myrowgids[0]);

  // communicate global IDs
  std::vector<int> rowgids(0, 0);
  LINALG::AllreduceVector(myrowgids, rowgids, comm);

  // retain communicated global IDs only on processor with ID 0
  if (comm.MyPID()) rowgids.clear();

  // create full row map on processor with ID 0
  const Epetra_Map fullrowmap(-1, rowgids.size(), rowgids.size() ? &rowgids[0] : NULL, 0, comm);

  // import matrix to processor with ID 0
  Epetra_CrsMatrix crsmatrix(Copy, fullrowmap, 0);
  if (sparsematrix != Teuchos::null)
  {
    if (crsmatrix.Import(*sparsematrix->EpetraMatrix(), Epetra_Import(fullrowmap, rowmap), Insert))
      dserror("Matrix import failed!");
  }
  else
  {
    for (int i = 0; i < blocksparsematrix->Rows(); ++i)
      for (int j = 0; j < blocksparsematrix->Cols(); ++j)
        if (crsmatrix.Import(*blocksparsematrix->Matrix(i, j).EpetraMatrix(),
                Epetra_Import(fullrowmap, blocksparsematrix->RangeMap(i)), Insert))
          dserror("Matrix import failed!");
  }

  // let processor with ID 0 output matrix to file
  if (comm.MyPID() == 0)
  {
    // set file name
    std::ostringstream nproc;
    nproc << comm.NumProc();
    const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName() +
                               ".matrix_" + nproc.str() + "proc.csv");

    // open file and write header at beginning
    std::ofstream file;
    file.open(filename.c_str(), std::fstream::trunc);
    file << std::setprecision(precision) << std::scientific << "RowGIDs,ColumnGIDs,Values"
         << std::endl;

    // write matrix to file
    for (int rowlid = 0; rowlid < crsmatrix.NumMyRows(); ++rowlid)
    {
      // extract global ID of current matrix row
      const int rowgid = fullrowmap.GID(rowlid);

      // extract current matrix row
      int numentries;
      double* values;
      int* indices;
      if (crsmatrix.ExtractGlobalRowView(rowgid, numentries, values, indices))
        dserror("Cannot extract matrix row with global ID %d!", rowgid);

      // sort entries in current matrix row in ascending order of column global ID via map
      std::map<int, double> entries;

      // loop over all entries in current matrix row
      for (int j = 0; j < numentries; ++j)
        // add current matrix entry to map
        entries[indices[j]] = values[j];

      // loop over all sorted entries in current matrix row
      for (std::map<int, double>::iterator j = entries.begin(); j != entries.end(); ++j)
        // write current matrix entry to file
        if (std::abs(j->second) > tolerance)
          file << rowgid << "," << j->first << "," << j->second << std::endl;
    }

    // close file
    file.close();
  }

  // wait until output is complete
  comm.Barrier();

  // throw error to abort simulation for debugging
  dserror("Matrix was output to *.csv file!");

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------------------------*
 | output vector to *.csv file for debugging purposes, with global IDs of vector components in
 ascending order across all processors   fang 01/17 |
 *------------------------------------------------------------------------------------------------------------------------------------------------*/
void STI::Monolithic::OutputVectorToFile(const Epetra_MultiVector& vector,  //!< vector to be output
    const int precision,                                                    //!< output precision
    const double tolerance  //!< output omission tolerance
)
{
  // extract communicator
  const Epetra_Comm& comm = vector.Comm();

  // extract vector map
  const Epetra_BlockMap& map = vector.Map();

  // safety check
  if (!map.UniqueGIDs())
    dserror("Vector output to *.csv file currently only works for non-overlapping vector maps!");

  // copy global IDs of vector components stored on current processor into vector
  std::vector<int> mygids(map.NumMyElements(), 0);
  int* myglobalelements = map.MyGlobalElements();
  std::copy(myglobalelements, myglobalelements + map.NumMyElements(), &mygids[0]);

  // communicate global IDs
  std::vector<int> gids(0, 0);
  LINALG::AllreduceVector(mygids, gids, comm);

  // retain communicated global IDs only on processor with ID 0
  if (comm.MyPID()) gids.clear();

  // create full vector map on processor with ID 0
  const Epetra_Map fullmap(-1, gids.size(), gids.size() ? &gids[0] : NULL, 0, comm);

  // export vector to processor with ID 0
  Epetra_MultiVector fullvector(fullmap, vector.NumVectors(), true);
  LINALG::Export(vector, fullvector);

  // let processor with ID 0 output vector to file
  if (comm.MyPID() == 0)
  {
    // set file name
    std::ostringstream nproc;
    nproc << comm.NumProc();
    const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName() +
                               ".vector_" + nproc.str() + "proc.csv");

    // open file and write header at beginning
    std::ofstream file;
    file.open(filename.c_str(), std::fstream::trunc);
    file << std::setprecision(precision) << std::scientific << "GIDs,Values" << std::endl;

    // write vector to file
    for (int lid = 0; lid < fullvector.MyLength(); ++lid)
    {
      // inner loop index
      int j(-1);

      // check output omission tolerance
      for (j = 0; j < fullvector.NumVectors(); ++j)
        if (std::abs(fullvector[j][lid]) > tolerance) break;

      // perform output if applicable
      if (j < fullvector.NumVectors())
      {
        // write global ID of current vector component
        file << fullmap.GID(lid);

        // loop over all subvectors
        for (j = 0; j < fullvector.NumVectors(); ++j)
          // write current vector component to file
          file << "," << fullvector[j][lid];

        // output line break
        file << std::endl;
      }
    }

    // close file
    file.close();
  }

  // wait until output is complete
  comm.Barrier();

  // throw error to abort simulation for debugging
  dserror("Vector was output to *.csv file!");

  return;
}


/*----------------------------------------------------------------------*
 | assemble global system of equations                       fang 07/15 |
 *----------------------------------------------------------------------*/
void STI::Monolithic::AssembleMatAndRHS()
{
  // pass thermo degrees of freedom to scatra discretization
  TransferThermoToScatra(ThermoField()->Phiafnp());

  // build system matrix and residual for scatra field
  ScaTraField()->PrepareLinearSolve();

  // pass scatra degrees of freedom to thermo discretization
  TransferScatraToThermo(ScaTraField()->Phiafnp());

  // build system matrix and residual for thermo field
  ThermoField()->PrepareLinearSolve();

  // evaluate scatra-thermo-OD coupling. Contributions from domain
  scatrathermooffdiagcoupling_->EvaluateOffDiagBlockScatraThermoDomain(scatrathermoblockdomain_);

  // evaluate scatra-thermo-OD coupling. Contributions from interface
  if (ScaTraField()->S2ICoupling())
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockScatraThermoInterface(
        scatrathermoblockinterface_);

  // evaluate thermo-scatra-OD coupling. Contributions from domain
  scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraDomain(thermoscatrablockdomain_);

  // evaluate thermo-scatra-OD coupling. Contributions from interface
  if (ScaTraField()->S2ICoupling())
    scatrathermooffdiagcoupling_->EvaluateOffDiagBlockThermoScatraInterface(
        thermoscatrablockinterface_);

  // OD blocks containing assembled domain and interface contributions
  Teuchos::RCP<LINALG::SparseOperator> scatrathermo_domain_interface(Teuchos::null);
  Teuchos::RCP<LINALG::SparseOperator> thermoscatra_domain_interface(Teuchos::null);

  // assemble interface and domain contributions of OD blocks
  AssembleDomainInterfaceOffDiag(scatrathermo_domain_interface, thermoscatra_domain_interface);

  // apply Dirichlet BCs to OD blocks
  ApplyDirichletOffDiag(scatrathermo_domain_interface, thermoscatra_domain_interface);

  // build global system matrix
  switch (matrixtype_)
  {
    case LINALG::MatrixType::block_condition:
    {
      // check global system matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_);
      if (blocksystemmatrix == Teuchos::null) dserror("System matrix is not a block matrix!");

      switch (ScaTraField()->MatrixType())
      {
        case LINALG::MatrixType::block_condition:
        {
          // extract number of matrix row or column blocks associated with scalar transport field
          const unsigned nblockmapsscatra = ScaTraField()->BlockMaps().NumMaps();

          // construct global system matrix by assigning matrix blocks
          for (unsigned iblock = 0; iblock < nblockmapsscatra; ++iblock)
          {
            for (unsigned jblock = 0; jblock < nblockmapsscatra; ++jblock)
              blocksystemmatrix->Assign(iblock, jblock, LINALG::View,
                  ScaTraField()->BlockSystemMatrix()->Matrix(iblock, jblock));

            // perform second condensation before assigning matrix blocks
            if (condensationthermo_)
            {
              const auto& scatrathermoblock =
                  Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(
                      scatrathermo_domain_interface)
                      ->Matrix(iblock, 0);
              LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                  scatrathermoblock.RangeMap(), *maps_->Map(1), 1.0, NULL, NULL,
                  blocksystemmatrix->Matrix(iblock, nblockmapsscatra));

              switch (strategyscatra_->CouplingType())
              {
                case INPAR::S2I::coupling_matching_nodes:
                {
                  ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
                  LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                      scatrathermoblock.RangeMap(), *strategythermo_->InterfaceMaps()->Map(1), 1.0,
                      NULL, &converter, blocksystemmatrix->Matrix(iblock, nblockmapsscatra), true,
                      true);
                  break;
                }
                case INPAR::S2I::coupling_mortar_standard:
                {
                  // initialize temporary matrix for slave-side columns of current matrix block
                  LINALG::SparseMatrix scatrathermocolsslave(scatrathermoblock.RowMap(), 81);

                  // fill temporary matrix for slave-side columns of current matrix block
                  LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                      scatrathermoblock.RangeMap(), *strategythermo_->InterfaceMaps()->Map(1), 1.0,
                      NULL, NULL, scatrathermocolsslave);

                  // finalize temporary matrix for slave-side columns of current matrix block
                  scatrathermocolsslave.Complete(
                      *strategythermo_->InterfaceMaps()->Map(1), scatrathermoblock.RowMap());

                  // transform and assemble temporary matrix for slave-side columns of current
                  // matrix block
                  blocksystemmatrix->Matrix(iblock, nblockmapsscatra)
                      .Add(*LINALG::MLMultiply(scatrathermocolsslave, *strategythermo_->P(), true),
                          false, 1.0, 1.0);

                  break;
                }
                default:
                {
                  dserror("Invalid type of scatra-scatra interface coupling!");
                  break;
                }
              }

              const auto& thermoscatrablock =
                  Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(
                      thermoscatra_domain_interface)
                      ->Matrix(0, iblock);
              LINALG::MatrixLogicalSplitAndTransform()(thermoscatrablock, *maps_->Map(1),
                  thermoscatrablock.DomainMap(), 1.0, NULL, NULL,
                  blocksystemmatrix->Matrix(nblockmapsscatra, iblock));
            }

            // assign matrix blocks directly
            else
            {
              blocksystemmatrix->Assign(iblock, nblockmapsscatra, LINALG::View,
                  Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(
                      scatrathermo_domain_interface)
                      ->Matrix(iblock, 0));

              blocksystemmatrix->Assign(nblockmapsscatra, iblock, LINALG::View,
                  Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(
                      thermoscatra_domain_interface)
                      ->Matrix(0, iblock));
            }
          }

          // perform second condensation before assigning thermo-thermo matrix block
          if (condensationthermo_)
          {
            LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(), *maps_->Map(1),
                *maps_->Map(1), 1.0, NULL, NULL,
                blocksystemmatrix->Matrix(nblockmapsscatra, nblockmapsscatra));

            switch (strategyscatra_->CouplingType())
            {
              case INPAR::S2I::coupling_matching_nodes:
              {
                ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
                LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(),
                    *maps_->Map(1), *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL,
                    &converter, blocksystemmatrix->Matrix(nblockmapsscatra, nblockmapsscatra), true,
                    true);
                break;
              }
              case INPAR::S2I::coupling_mortar_standard:
              {
                // initialize temporary matrix for slave-side columns of thermo-thermo matrix block
                LINALG::SparseMatrix thermothermocolsslave(*maps_->Map(1), 81);

                // fill temporary matrix for slave-side columns of thermo-thermo matrix block
                LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(),
                    *maps_->Map(1), *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL, NULL,
                    thermothermocolsslave);

                // finalize temporary matrix for slave-side columns of thermo-thermo matrix block
                thermothermocolsslave.Complete(
                    *strategythermo_->InterfaceMaps()->Map(1), *maps_->Map(1));

                // transform and assemble temporary matrix for slave-side columns of thermo-thermo
                // matrix block
                blocksystemmatrix->Matrix(nblockmapsscatra, nblockmapsscatra)
                    .Add(*LINALG::MLMultiply(thermothermocolsslave, *strategythermo_->P(), true),
                        false, 1.0, 1.0);

                break;
              }
              default:
              {
                dserror("Invalid type of scatra-scatra interface coupling!");
                break;
              }
            }
          }

          // assign thermo-thermo matrix block directly
          else
            blocksystemmatrix->Assign(
                nblockmapsscatra, nblockmapsscatra, LINALG::View, *ThermoField()->SystemMatrix());

          break;
        }

        case LINALG::MatrixType::sparse:
        {
          // construct global system matrix by assigning matrix blocks
          blocksystemmatrix->Assign(0, 0, LINALG::View, *ScaTraField()->SystemMatrix());

          // perform second condensation before assigning matrix blocks
          if (condensationthermo_)
          {
            const auto& scatrathermoblock = *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(
                scatrathermo_domain_interface);
            LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                scatrathermoblock.RangeMap(), *maps_->Map(1), 1.0, NULL, NULL,
                blocksystemmatrix->Matrix(0, 1));

            LINALG::MatrixLogicalSplitAndTransform()(
                *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(
                    thermoscatra_domain_interface),
                *maps_->Map(1), thermoscatra_domain_interface->DomainMap(), 1.0, NULL, NULL,
                blocksystemmatrix->Matrix(1, 0));

            LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(), *maps_->Map(1),
                *maps_->Map(1), 1.0, NULL, NULL, blocksystemmatrix->Matrix(1, 1));

            switch (strategyscatra_->CouplingType())
            {
              case INPAR::S2I::coupling_matching_nodes:
              {
                ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
                LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                    scatrathermoblock.RangeMap(), *strategythermo_->InterfaceMaps()->Map(1), 1.0,
                    NULL, &converter, blocksystemmatrix->Matrix(0, 1), true, true);

                LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(),
                    *maps_->Map(1), *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL,
                    &converter, blocksystemmatrix->Matrix(1, 1), true, true);

                break;
              }

              case INPAR::S2I::coupling_mortar_standard:
              {
                // initialize temporary matrix for slave-side columns of scatra-thermo matrix block
                LINALG::SparseMatrix scatrathermocolsslave(
                    *ScaTraField()->Discretization()->DofRowMap(), 81);

                // fill temporary matrix for slave-side columns of scatra-thermo matrix block
                LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                    scatrathermoblock.RangeMap(), *strategythermo_->InterfaceMaps()->Map(1), 1.0,
                    NULL, NULL, scatrathermocolsslave);

                // finalize temporary matrix for slave-side columns of scatra-thermo matrix block
                scatrathermocolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1),
                    *ScaTraField()->Discretization()->DofRowMap());

                // transform and assemble temporary matrix for slave-side columns of scatra-thermo
                // matrix block
                blocksystemmatrix->Matrix(0, 1).Add(
                    *LINALG::MLMultiply(scatrathermocolsslave, *strategythermo_->P(), true), false,
                    1.0, 1.0);

                // initialize temporary matrix for slave-side columns of thermo-thermo matrix block
                LINALG::SparseMatrix thermothermocolsslave(*maps_->Map(1), 81);

                // fill temporary matrix for slave-side columns of thermo-thermo matrix block
                LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(),
                    *maps_->Map(1), *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL, NULL,
                    thermothermocolsslave);

                // finalize temporary matrix for slave-side columns of thermo-thermo matrix block
                thermothermocolsslave.Complete(
                    *strategythermo_->InterfaceMaps()->Map(1), *maps_->Map(1));

                // transform and assemble temporary matrix for slave-side columns of thermo-thermo
                // matrix block
                blocksystemmatrix->Matrix(1, 1).Add(
                    *LINALG::MLMultiply(thermothermocolsslave, *strategythermo_->P(), true), false,
                    1.0, 1.0);

                break;
              }

              default:
              {
                dserror("Invalid type of scatra-scatra interface coupling!");
                break;
              }
            }
          }

          // assign matrix blocks directly
          else
          {
            blocksystemmatrix->Assign(0, 1, LINALG::View,
                *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatrathermo_domain_interface));
            blocksystemmatrix->Assign(1, 0, LINALG::View,
                *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatra_domain_interface));
            blocksystemmatrix->Assign(1, 1, LINALG::View, *ThermoField()->SystemMatrix());
          }

          break;
        }

        default:
        {
          dserror("Invalid matrix type associated with scalar transport field!");
          break;
        }
      }

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // check global system matrix
      auto systemmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
      if (systemmatrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

      // construct global system matrix by adding matrix blocks
      systemmatrix->Add(*ScaTraField()->SystemMatrix(), false, 1.0, 0.0);

      // perform second condensation before adding matrix blocks
      if (condensationthermo_)
      {
        const auto& scatrathermoblock =
            *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrathermo_domain_interface);
        LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock, scatrathermoblock.RangeMap(),
            *maps_->Map(1), 1.0, NULL, NULL, *systemmatrix, true, true);

        LINALG::MatrixLogicalSplitAndTransform()(
            *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatra_domain_interface),
            *maps_->Map(1), thermoscatra_domain_interface->DomainMap(), 1.0, NULL, NULL,
            *systemmatrix, true, true);

        LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(), *maps_->Map(1),
            *maps_->Map(1), 1.0, NULL, NULL, *systemmatrix, true, true);

        switch (strategyscatra_->CouplingType())
        {
          case INPAR::S2I::coupling_matching_nodes:
          {
            ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
            LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                scatrathermoblock.RangeMap(), *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL,
                &converter, *systemmatrix, true, true);

            LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(), *maps_->Map(1),
                *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL, &converter, *systemmatrix,
                true, true);

            break;
          }

          case INPAR::S2I::coupling_mortar_standard:
          {
            // initialize temporary matrix for slave-side columns of global system matrix
            LINALG::SparseMatrix systemmatrixcolsslave(*DofRowMap(), 81);

            // fill temporary matrix for slave-side columns of global system matrix
            LINALG::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                scatrathermoblock.RangeMap(), *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL,
                NULL, systemmatrixcolsslave);

            LINALG::MatrixLogicalSplitAndTransform()(*ThermoField()->SystemMatrix(), *maps_->Map(1),
                *strategythermo_->InterfaceMaps()->Map(1), 1.0, NULL, NULL, systemmatrixcolsslave,
                true, true);

            // finalize temporary matrix for slave-side columns of global system matrix
            systemmatrixcolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1), *DofRowMap());

            // transform and assemble temporary matrix for slave-side columns of global system
            // matrix
            systemmatrix->Add(
                *LINALG::MLMultiply(systemmatrixcolsslave, *strategythermo_->P(), true), false, 1.0,
                1.0);

            break;
          }

          default:
          {
            dserror("Invalid type of scatra-scatra interface coupling!");
            break;
          }
        }
      }

      // add matrix blocks directly
      else
      {
        systemmatrix->Add(*scatrathermo_domain_interface, false, 1.0, 1.0);
        systemmatrix->Add(*thermoscatra_domain_interface, false, 1.0, 1.0);
        systemmatrix->Add(*ThermoField()->SystemMatrix(), false, 1.0, 1.0);
      }

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scatra-thermo interaction not recognized!");
      break;
    }
  }

  // finalize global system matrix
  systemmatrix_->Complete();

  // create full monolithic right-hand side vector
  maps_->InsertVector(ScaTraField()->Residual(), 0, residual_);
  Teuchos::RCP<Epetra_Vector> thermoresidual(Teuchos::null);
  if (condensationthermo_)
  {
    thermoresidual = Teuchos::rcp(new Epetra_Vector(*maps_->Map(1)));
    LINALG::Export(*ThermoField()->Residual(), *thermoresidual);
  }
  else
    thermoresidual = ThermoField()->Residual();
  maps_->InsertVector(thermoresidual, 1, residual_);

  return;
}  // STI::Monolithic::AssembleMatAndRHS()


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void STI::Monolithic::BuildNullSpaces() const
{
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      ScaTraField()->BuildBlockNullSpaces(solver_, 0);
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sublists to trigger null space
      // computation
      Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse1");
      blocksmootherparams.sublist("Aztec Parameters");
      blocksmootherparams.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with temperature field as string
  std::stringstream iblockstr;
  iblockstr << blockmaps_->NumMaps();

  // equip smoother for thermo matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for thermo matrix block with null space associated with all degrees of freedom
  // on thermo discretization
  ThermoField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

  // reduce full null space to match degrees of freedom associated with thermo matrix block if
  // necessary
  if (condensationthermo_)
    LINALG::Solver::FixMLNullspace("Block " + iblockstr.str(),
        *ThermoField()->Discretization()->DofRowMap(), *maps_->Map(1), blocksmootherparams);

  return;
}  // STI::Monolithic::BuildBlockNullSpaces

/*------------------------------------------------------------------------------------------------*
 | compute null space information associated with global system matrix if applicable   fang 06/17 |
 *------------------------------------------------------------------------------------------------*/
void STI::Monolithic::ComputeNullSpaceIfNecessary(
    Teuchos::ParameterList& solverparams  //! solver parameter list for scatra-thermo interaction
    ) const
{
  // compute vector-based null space information for ML preconditioner
  if (solverparams.isSublist("ML Parameters"))
  {
    // extract parameter list for ML preconditioner
    Teuchos::ParameterList& mllist = solverparams.sublist("ML Parameters", true);

    // determine null space dimension
    const int numdofpernode_scatra = ScaTraField()->NumDofPerNode();
    const int numdofpernode_thermo = ThermoField()->NumDofPerNode();
    const int dimns = numdofpernode_scatra + numdofpernode_thermo;

    // allocate vector for null space information
    const auto ns = Teuchos::rcp(new std::vector<double>(dimns * DofRowMap()->NumMyElements(), 0.));

    // compute null space modes associated with scatra field
    const DRT::Discretization& scatradis = *ScaTraField()->Discretization();
    double* modes_scatra[numdofpernode_scatra];
    for (int i = 0; i < numdofpernode_scatra; ++i)
      modes_scatra[i] = &((*ns)[i * DofRowMap()->NumMyElements()]);
    for (int i = 0; i < scatradis.NumMyRowNodes(); ++i)
    {
      const int lid = DofRowMap()->LID(scatradis.Dof(0, scatradis.lRowNode(i), 0));
      if (lid < 0) dserror("Cannot find scatra degree of freedom!");
      for (int j = 0; j < numdofpernode_scatra; ++j) modes_scatra[j][lid + j] = 1.;
    }

    // compute null space modes associated with thermo field
    const DRT::Discretization& thermodis = *ThermoField()->Discretization();
    double* modes_thermo[numdofpernode_thermo];
    for (int i = 0; i < numdofpernode_thermo; ++i)
      modes_thermo[i] = &((*ns)[(numdofpernode_scatra + i) * DofRowMap()->NumMyElements()]);
    for (int i = 0; i < thermodis.NumMyRowNodes(); ++i)
    {
      const int lid = DofRowMap()->LID(thermodis.Dof(0, thermodis.lRowNode(i), 0));
      if (lid < 0) dserror("Cannot find thermo degree of freedom!");
      for (int j = 0; j < numdofpernode_thermo; ++j) modes_thermo[j][lid + j] = 1.;
    }

    // fill parameter list
    mllist.set("PDE equations", dimns);
    mllist.set("null space: dimension", dimns);
    mllist.set("null space: type", "pre-computed");
    mllist.set("null space: add default vectors", false);
    mllist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
    mllist.set("null space: vectors", &((*ns)[0]));
    mllist.set<bool>("ML validate parameter list", false);
  }

  // compute point-based null space information for MueLu preconditioner
  else if (solverparams.isSublist("MueLu Parameters"))
  {
    // extract and fill parameter list for MueLu preconditioner
    Teuchos::ParameterList& mllist = solverparams.sublist("MueLu Parameters", true);
    mllist.set("PDE equations", 1);
    mllist.set("null space: dimension", 1);
    mllist.set("null space: type", "pre-computed");
    mllist.set("null space: add default vectors", false);
    const Teuchos::RCP<std::vector<double>> ns =
        Teuchos::rcp(new std::vector<double>(DofRowMap()->NumMyElements(), 1.));
    mllist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
    mllist.set("null space: vectors", &((*ns)[0]));
    mllist.set<bool>("ML validate parameter list", false);
  }

  return;
}  // STI::Monolithic::ComputeNullSpaceIfNecessary


/*----------------------------------------------------------------------*
 | global map of degrees of freedom                          fang 04/15 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& STI::Monolithic::DofRowMap() const
{
  return maps_->FullMap();
}  // STI::Monolithic::DofRowMap()



/*-----------------------------------------------------------------------*
 | check termination criterion for Newton-Raphson iteration   fang 04/15 |
 *-----------------------------------------------------------------------*/
bool STI::Monolithic::ExitNewtonRaphson()
{
  // initialize exit flag
  bool exit(false);

  // perform Newton-Raphson convergence check depending on type of scalar transport
  switch (
      Teuchos::getIntegralValue<INPAR::STI::ScaTraTimIntType>(*stiparameters_, "SCATRATIMINTTYPE"))
  {
    case INPAR::STI::ScaTraTimIntType::elch:
    {
      // compute L2 norm of concentration state vector
      double concdofnorm(0.);
      ScaTraField()->Splitter()->ExtractOtherVector(*ScaTraField()->Phinp())->Norm2(&concdofnorm);

      // compute L2 norm of concentration residual vector
      double concresnorm(0.);
      ScaTraField()
          ->Splitter()
          ->ExtractOtherVector(*maps_->ExtractVector(residual_, 0))
          ->Norm2(&concresnorm);

      // compute L2 norm of concentration increment vector
      double concincnorm(0.);
      ScaTraField()
          ->Splitter()
          ->ExtractOtherVector(*maps_->ExtractVector(increment_, 0))
          ->Norm2(&concincnorm);

      // compute L2 norm of potential state vector
      double potdofnorm(0.);
      ScaTraField()->Splitter()->ExtractCondVector(*ScaTraField()->Phinp())->Norm2(&potdofnorm);

      // compute L2 norm of potential residual vector
      double potresnorm(0.);
      ScaTraField()
          ->Splitter()
          ->ExtractCondVector(*maps_->ExtractVector(residual_, 0))
          ->Norm2(&potresnorm);

      // compute L2 norm of potential increment vector
      double potincnorm(0.);
      ScaTraField()
          ->Splitter()
          ->ExtractCondVector(*maps_->ExtractVector(increment_, 0))
          ->Norm2(&potincnorm);

      // compute L2 norm of thermo state vector
      double thermodofnorm(0.);
      ThermoField()->Phinp()->Norm2(&thermodofnorm);

      // compute L2 norm of thermo residual vector
      double thermoresnorm(0.);
      maps_->ExtractVector(residual_, 1)->Norm2(&thermoresnorm);

      // compute L2 norm of thermo increment vector
      double thermoincnorm(0.);
      maps_->ExtractVector(increment_, 1)->Norm2(&thermoincnorm);

      // safety checks
      if (std::isnan(concdofnorm) or std::isnan(concresnorm) or std::isnan(concincnorm) or
          std::isnan(potdofnorm) or std::isnan(potresnorm) or std::isnan(potincnorm) or
          std::isnan(thermodofnorm) or std::isnan(thermoresnorm) or std::isnan(thermoincnorm))
        dserror("Vector norm is not a number!");
      if (std::isinf(concdofnorm) or std::isinf(concresnorm) or std::isinf(concincnorm) or
          std::isinf(potdofnorm) or std::isinf(potresnorm) or std::isinf(potincnorm) or
          std::isinf(thermodofnorm) or std::isinf(thermoresnorm) or std::isinf(thermoincnorm))
        dserror("Vector norm is infinity!");

      // prevent division by zero
      if (concdofnorm < 1.e-10) concdofnorm = 1.e-10;
      if (potdofnorm < 1.e-10) potdofnorm = 1.e-10;
      if (thermodofnorm < 1.e-10) thermodofnorm = 1.e-10;

      // first Newton-Raphson iteration
      if (iter_ == 1)
      {
        if (Comm().MyPID() == 0)
        {
          // print header of convergence table to screen
          std::cout << "+------------+-------------------+--------------+--------------+-----------"
                       "---+--------------+--------------+--------------+"
                    << std::endl;
          std::cout << "|- step/max -|- tolerance[norm] -|-- conc-res --|-- conc-inc --|-- pot-res "
                       "---|-- pot-inc ---|- thermo-res -|- thermo-inc -|"
                    << std::endl;

          // print first line of convergence table to screen
          // solution increment not yet available during first Newton-Raphson iteration
          std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                    << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
                    << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                    << concresnorm << "   |      --      | " << std::setw(10)
                    << std::setprecision(3) << std::scientific << potresnorm
                    << "   |      --      | " << std::setw(10) << std::setprecision(3)
                    << std::scientific << thermoresnorm << "   |      --      | "
                    << "(       --      , te = " << std::setw(10) << std::setprecision(3) << dtele_
                    << ")" << std::endl;
        }
      }

      // subsequent Newton-Raphson iterations
      else
      {
        // print current line of convergence table to screen
        if (Comm().MyPID() == 0)
          std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                    << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
                    << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                    << concresnorm << "   | " << std::setw(10) << std::setprecision(3)
                    << std::scientific << concincnorm / concdofnorm << "   | " << std::setw(10)
                    << std::setprecision(3) << std::scientific << potresnorm << "   | "
                    << std::setw(10) << std::setprecision(3) << std::scientific
                    << potincnorm / potdofnorm << "   | " << std::setw(10) << std::setprecision(3)
                    << std::scientific << thermoresnorm << "   | " << std::setw(10)
                    << std::setprecision(3) << std::scientific << thermoincnorm / thermodofnorm
                    << "   | (ts = " << std::setw(10) << std::setprecision(3) << dtsolve_
                    << ", te = " << std::setw(10) << std::setprecision(3) << dtele_ << ")"
                    << std::endl;

        // convergence check
        if (concresnorm <= itertol_ and potresnorm <= itertol_ and thermoresnorm <= itertol_ and
            concincnorm / concdofnorm <= itertol_ and potincnorm / potdofnorm <= itertol_ and
            thermoincnorm / thermodofnorm <= itertol_)
          // exit Newton-Raphson iteration upon convergence
          exit = true;
      }

      // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary
      // additional solver calls
      if (concresnorm < restol_ and potresnorm < restol_ and thermoresnorm < restol_) exit = true;

      // print warning to screen if maximum number of Newton-Raphson iterations is reached without
      // convergence
      if (iter_ == itermax_ and !exit)
      {
        if (Comm().MyPID() == 0)
        {
          std::cout << "+------------+-------------------+--------------+--------------+-----------"
                       "---+--------------+--------------+--------------+"
                    << std::endl;
          std::cout << "|                     Newton-Raphson method has not converged after a "
                       "maximum number of "
                    << std::setw(2) << itermax_ << " iterations!                     |"
                    << std::endl;
        }

        // proceed to next time step
        exit = true;
      }

      // print finish line of convergence table to screen
      if (exit and Comm().MyPID() == 0)
        std::cout << "+------------+-------------------+--------------+--------------+-------------"
                     "-+--------------+--------------+--------------+"
                  << std::endl
                  << std::endl;

      break;
    }

    default:
    {
      dserror("Type of scalar transport not yet available!");
      break;
    }
  }

  return exit;
}  // STI::Monolithic::ExitNewtonRaphson()


/*----------------------------------------------------------------------*
 | prepare time step                                         fang 09/17 |
 *----------------------------------------------------------------------*/
void STI::Monolithic::PrepareTimeStep()
{
  // call base class routine
  Algorithm::PrepareTimeStep();

  // print time step information to screen
  ScaTraField()->PrintTimeStepInfo();

  return;
}  // STI::Monolithic::PrepareTimeStep()


/*----------------------------------------------------------------------*
 | evaluate time step using Newton-Raphson iteration         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Monolithic::Solve()
{
  // initialize counter for Newton-Raphson iterations
  iter_ = 0;

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    ++iter_;

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // assemble global system of equations
    AssembleMatAndRHS();

    // determine time needed for evaluating elements and assembling global system of equations,
    // and take maximum over all processors via communication
    double mydtele = timer_->WallTime() - time;
    Comm().MaxAll(&mydtele, &dtele_, 1);

    // safety check
    if (!systemmatrix_->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if (ScaTraField()->FDCheckType() == INPAR::SCATRA::fdcheck_global) FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if (ExitNewtonRaphson()) break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    time = timer_->WallTime();

    // equilibrate global system of equations if necessary
    equilibration_->EquilibrateSystem(systemmatrix_, residual_, *blockmaps_);

    // solve global system of equations
    // Dirichlet boundary conditions have already been applied to global system of equations
    solver_->Solve(systemmatrix_->EpetraOperator(), increment_, residual_, true, iter_ == 1);

    equilibration_->UnequilibrateIncrement(increment_);

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime() - time;
    Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

    // output performance statistics associated with linear solver into text file if applicable
    if (DRT::INPUT::IntegralValue<int>(*fieldparameters_, "OUTPUTLINSOLVERSTATS"))
      ScaTraField()->OutputLinSolverStats(
          *solver_, dtsolve_, Step(), iter_, residual_->Map().NumGlobalElements());

    // update scatra field
    ScaTraField()->UpdateIter(maps_->ExtractVector(increment_, 0));
    ScaTraField()->ComputeIntermediateValues();

    // update thermo field
    Teuchos::RCP<Epetra_Vector> thermoincrement(Teuchos::null);
    if (condensationthermo_)
    {
      thermoincrement =
          Teuchos::rcp(new Epetra_Vector(*ThermoField()->Discretization()->DofRowMap()));
      LINALG::Export(*maps_->ExtractVector(increment_, 1), *thermoincrement);
      const Teuchos::RCP<const Epetra_Vector> masterincrement =
          strategythermo_->InterfaceMaps()->ExtractVector(*thermoincrement, 2);
      const Teuchos::RCP<Epetra_Vector> slaveincrement =
          LINALG::CreateVector(*strategythermo_->InterfaceMaps()->Map(1));
      switch (strategyscatra_->CouplingType())
      {
        case INPAR::S2I::coupling_matching_nodes:
        {
          icoupthermo_->MasterToSlave(masterincrement, slaveincrement);
          break;
        }
        case INPAR::S2I::coupling_mortar_standard:
        {
          strategythermo_->P()->Multiply(false, *masterincrement, *slaveincrement);
          break;
        }
        default:
        {
          dserror("Invalid type of scatra-scatra interface coupling!");
          break;
        }
      }
      strategythermo_->InterfaceMaps()->InsertVector(slaveincrement, 1, thermoincrement);
    }
    else
      thermoincrement = maps_->ExtractVector(increment_, 1);
    ThermoField()->UpdateIter(thermoincrement);
    ThermoField()->ComputeIntermediateValues();
  }  // Newton-Raphson iteration

  return;
}  // STI::Monolithic::Solve

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::Monolithic::ApplyDirichletOffDiag(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermo_domain_interface,
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatra_domain_interface)
{
  // apply Dirichlet boundary conditions to scatra-thermo matrix block
  scatrathermo_domain_interface->ApplyDirichlet(*ScaTraField()->DirichMaps()->CondMap(), false);

  // apply Dirichlet boundary conditions to scatra-thermo matrix block
  thermoscatra_domain_interface->ApplyDirichlet(*ThermoField()->DirichMaps()->CondMap(), false);

  // zero out slave-side rows of thermo-scatra matrix block after having added them to the
  // corresponding master-side rows to finalize condensation of slave-side thermo dofs
  if (ThermoField()->S2ICoupling())
  {
    switch (strategythermo_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        if (!ThermoField()->Discretization()->GetCondition("PointCoupling"))
          thermoscatra_domain_interface->ApplyDirichlet(*icoupthermo_->SlaveDofMap(), false);
        break;
      }

      case INPAR::S2I::coupling_mortar_condensed_bubnov:
      {
        thermoscatra_domain_interface->ApplyDirichlet(
            *(strategythermo_->InterfaceMaps()->Map(1)), false);
        break;
      }

      default:
      {
        dserror("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }

  return;
}  // STI::Monolithic::ApplyDirichletOD

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::Monolithic::AssembleDomainInterfaceOffDiag(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermo_domain_interface,
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatra_domain_interface)
{
  // initialize scatra-thermo blocks
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      scatrathermo_domain_interface =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, ScaTraField()->BlockMaps(), 81, false, true));
      thermoscatra_domain_interface =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              ScaTraField()->BlockMaps(), *blockmapthermo_, 81, false, true));

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      scatrathermo_domain_interface = Teuchos::rcp(
          new LINALG::SparseMatrix(*ScaTraField()->Discretization()->DofRowMap(), 27, false, true));
      thermoscatra_domain_interface = Teuchos::rcp(
          new LINALG::SparseMatrix(*ThermoField()->Discretization()->DofRowMap(), 27, false, true));
      break;
    }
    default:
    {
      dserror("Unknown matrix type");
      break;
    }
  }

  // add domain and interface contibutions
  scatrathermo_domain_interface->Add(*scatrathermoblockdomain_, false, 1.0, 1.0);
  scatrathermo_domain_interface->Add(*scatrathermoblockinterface_, false, 1.0, 1.0);
  thermoscatra_domain_interface->Add(*thermoscatrablockdomain_, false, 1.0, 1.0);
  thermoscatra_domain_interface->Add(*thermoscatrablockinterface_, false, 1.0, 1.0);

  if (strategythermo_->CouplingType() == INPAR::S2I::coupling_matching_nodes)
  {  // standard meshtying algorithm with Lagrange multipliers condensed out
    if (!ThermoField()->Discretization()->GetCondition("PointCoupling"))
    {
      // during the very first run of the following code, Complete() has not yet been called
      // on the thermo-scatra block of the global system matrix experiments have shown that
      // the ExtractMatrixRows routine called in the following does not properly work in this
      // case if the thermo-scatra block exhibits a block structure in particular, some matrix
      // entries are simply omitted during extraction, although the underlying routine
      // ExtractGlobalRowCopy returns a zero error code as a consequence, the final
      // thermo-scatra block lacks some entries, and so does the final matrix graph after
      // calling Complete() this leads to errors involving unknown global indices of matrix
      // columns during the second run of the following code, with the thermo-scatra block now
      // being Complete() to circumvent this problem, we call Complete() on the thermo-scatra
      // block before the very first run of the following code
      if (ScaTraField()->MatrixType() == LINALG::MatrixType::block_condition and Step() == 1 and
          iter_ == 1)
        thermoscatra_domain_interface->Complete();

      // loop over all thermo-scatra matrix blocks
      for (int iblock = 0; iblock < blockmaps_->NumMaps() - 1; ++iblock)
      {
        // initialize temporary matrix for slave-side rows of current thermo-scatra matrix
        // block
        LINALG::SparseMatrix thermoscatrarowsslave(*icoupthermo_->SlaveDofMap(), 27, false, true);

        // extract current thermo-scatra matrix block
        auto& thermoscatrablock =
            ScaTraField()->MatrixType() == LINALG::MatrixType::block_condition
                ? Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(
                      thermoscatra_domain_interface)
                      ->Matrix(0, iblock)
                : *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatra_domain_interface);

        // extract slave-side rows of thermo-scatra matrix block into temporary matrix
        SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(
            thermoscatrablock, thermoscatrarowsslave, *icoupthermo_->SlaveDofMap());

        // finalize temporary matrix with slave-side rows of thermo-scatra matrix block
        thermoscatrarowsslave.Complete(*maps_->Map(0), *icoupthermo_->SlaveDofMap());

        // undo Complete() from above before performing subsequent matrix row transformation
        if (ScaTraField()->MatrixType() == LINALG::MatrixType::block_condition and Step() == 1 and
            iter_ == 1)
          thermoscatrablock.UnComplete();

        // add slave-side rows of thermo-scatra matrix block to corresponding slave-side rows
        const Teuchos::RCP<LINALG::MatrixRowTransform> islavetomasterrowtransformthermood =
            ScaTraField()->MatrixType() == LINALG::MatrixType::block_condition
                ? Teuchos::rcp(new LINALG::MatrixRowTransform())
                : islavetomasterrowtransformthermood_;
        (*islavetomasterrowtransformthermood)(thermoscatrarowsslave, 1.,
            ADAPTER::CouplingSlaveConverter(*icoupthermo_), thermoscatrablock, true);
      }
    }
  }
  else if (strategythermo_->CouplingType() == INPAR::S2I::coupling_mortar_condensed_bubnov)
  {
    // standard meshtying algorithm with Lagrange multipliers condensed out
    // loop over all thermo-scatra matrix blocks
    for (int iblock = 0; iblock < blockmaps_->NumMaps() - 1; ++iblock)
    {
      // initialize temporary matrix for slave-side rows of current thermo-scatra matrix block
      LINALG::SparseMatrix thermoscatrarowsslave(
          *strategythermo_->InterfaceMaps()->Map(1), 27, false, true);

      // extract current thermo-scatra matrix block
      auto& thermoscatrablock =
          ScaTraField()->MatrixType() == LINALG::MatrixType::block_condition
              ? Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(
                    thermoscatra_domain_interface)
                    ->Matrix(0, iblock)
              : *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatra_domain_interface);

      // extract slave-side rows of thermo-scatra matrix block into temporary matrix
      SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(
          thermoscatrablock, thermoscatrarowsslave, *strategythermo_->InterfaceMaps()->Map(1));

      // finalize temporary matrix with slave-side rows of thermo-scatra matrix block
      thermoscatrarowsslave.Complete(*maps_->Map(0), *strategythermo_->InterfaceMaps()->Map(1));

      // add projected slave-side rows of thermo-scatra matrix block to corresponding
      // master-side rows
      thermoscatrablock.Add(*LINALG::MLMultiply(*strategythermo_->P(), true, thermoscatrarowsslave,
                                false, false, false, true),
          false, 1.0, 1.0);
    }
  }

  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      scatrathermo_domain_interface->Complete();
      thermoscatra_domain_interface->Complete();

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      scatrathermo_domain_interface->Complete(
          *ThermoField()->Discretization()->DofRowMap(), *maps_->Map(0));
      thermoscatra_domain_interface->Complete(*maps_->Map(0), *maps_->Map(1));

      break;
    }
    default:
    {
      dserror("Unknown matrix type");
      break;
    }
  }

  return;
}  // AssembleDomainInterfaceOD
