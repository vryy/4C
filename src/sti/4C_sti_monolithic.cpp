/*----------------------------------------------------------------------*/
/*! \file

\brief monolithic coupling algorithm for scatra-thermo interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_sti_monolithic.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_sti_monolithic_evaluate_OffDiag.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
STI::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& stidyn,
    const Teuchos::ParameterList& scatradyn, const Teuchos::ParameterList& solverparams,
    const Teuchos::ParameterList& solverparams_scatra,
    const Teuchos::ParameterList& solverparams_thermo)
    : Algorithm(comm, stidyn, scatradyn, solverparams_scatra, solverparams_thermo),
      restol_(fieldparameters_->sublist("NONLINEAR").get<double>("ABSTOLRES")),
      maps_(Teuchos::null),
      condensationthermo_(Core::UTILS::IntegralValue<bool>(stidyn, "THERMO_CONDENSATION")),
      systemmatrix_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<Core::LinAlg::MatrixType>(
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
      solver_(Teuchos::rcp(new Core::LinAlg::Solver(solverparams, comm,
          Global::Problem::instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::instance()->io_params(), "VERBOSITY")))),
      invrowsums_(Teuchos::null),
      icoupscatra_(Teuchos::null),
      icoupthermo_(Teuchos::null),
      islavetomasterrowtransformscatraod_(Teuchos::null),
      islavetomastercoltransformthermood_(Teuchos::null),
      islavetomasterrowtransformthermood_(Teuchos::null),
      equilibration_(Teuchos::null)
{
  // safety checks
  if (!scatra_field()->is_incremental())
    FOUR_C_THROW("Must have incremental solution approach for scatra-thermo interaction!");
  if (thermo_field()->system_matrix() == Teuchos::null)
    FOUR_C_THROW("System matrix associated with temperature field must be a sparse matrix!");

  // set control parameters for Newton-Raphson iteration loop
  itermax_ = fieldparameters_->sublist("NONLINEAR").get<int>("ITEMAX");
  itertol_ = fieldparameters_->sublist("NONLINEAR").get<double>("CONVTOL");


  // extract coupling adapters for scatra-scatra interface mesh tying
  if (scatra_field()->s2_i_meshtying() and
      strategyscatra_->coupling_type() == Inpar::S2I::coupling_matching_nodes)
  {
    icoupscatra_ = strategyscatra_->coupling_adapter();
    icoupthermo_ = strategythermo_->coupling_adapter();
  }

  // initialize map associated with single thermo block of global system matrix
  Teuchos::RCP<const Epetra_Map> mapthermo(Teuchos::null);
  if (condensationthermo_)
  {
    mapthermo = Core::LinAlg::MergeMap(
        *strategythermo_->interface_maps()->Map(0), *strategythermo_->interface_maps()->Map(2));
  }
  else
    mapthermo = thermo_field()->dof_row_map();

  // initialize global map extractor
  maps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(
      *Core::LinAlg::MergeMap(*scatra_field()->discretization()->dof_row_map(), *mapthermo, false),
      mapthermo, scatra_field()->dof_row_map()));

  // check global map extractor
  maps_->check_for_valid_map_extractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = Core::LinAlg::CreateVector(*dof_row_map(), true);

  // initialize global residual vector
  residual_ = Core::LinAlg::CreateVector(*dof_row_map(), true);

  // initialize transformation operators
  islavetomasterrowtransformscatraod_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
  islavetomastercoltransformthermood_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  islavetomasterrowtransformthermood_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);

  // merge slave and master side block maps for interface matrix for thermo and scatra
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);
  Teuchos::RCP<Epetra_Map> interface_map_thermo(Teuchos::null);
  Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockmapscatrainterface(Teuchos::null);
  Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockmapthermointerface(Teuchos::null);
  Teuchos::RCP<Core::LinAlg::MultiMapExtractor> blockmapthermointerfaceslave(Teuchos::null);

  if (scatra_field()->s2_i_meshtying())
  {
    // merge slave and master side full maps for interface matrix for thermo and scatra
    interface_map_scatra = Core::LinAlg::MultiMapExtractor::merge_maps(
        {strategyscatra_->interface_maps()->Map(1), strategyscatra_->interface_maps()->Map(2)});
    interface_map_thermo = Core::LinAlg::MultiMapExtractor::merge_maps(
        {strategythermo_->interface_maps()->Map(1), strategythermo_->interface_maps()->Map(2)});
    // build block map for thermo interface by using full thermo map
    blockmapthermointerface =
        Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*interface_map_thermo,
            std::vector<Teuchos::RCP<const Epetra_Map>>(1, interface_map_thermo)));
    blockmapthermointerface->check_for_valid_map_extractor();
    blockmapthermointerfaceslave = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
        *strategythermo_->interface_maps()->Map(1),
        std::vector<Teuchos::RCP<const Epetra_Map>>(1, strategythermo_->interface_maps()->Map(1))));
    blockmapthermointerfaceslave->check_for_valid_map_extractor();
  }

  // initialize map extractors associated with blocks of global system matrix
  switch (scatra_field()->matrix_type())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case Core::LinAlg::MatrixType::sparse:
    {
      blockmaps_ = maps_;

      if (scatra_field()->s2_i_meshtying())
      {
        blockmapscatrainterface =
            Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*interface_map_scatra,
                std::vector<Teuchos::RCP<const Epetra_Map>>(1, interface_map_scatra)));
        blockmapscatrainterface->check_for_valid_map_extractor();
      }
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case Core::LinAlg::MatrixType::block_condition:
    {
      // extract maps underlying main-diagonal matrix blocks associated with scalar transport field
      const int nblockmapsscatra = static_cast<int>(scatra_field()->block_maps()->num_maps());
      std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps(nblockmapsscatra + 1);
      for (int iblockmap = 0; iblockmap < nblockmapsscatra; ++iblockmap)
        blockmaps[iblockmap] = scatra_field()->block_maps()->Map(iblockmap);

      // extract map underlying single main-diagonal matrix block associated with temperature field
      blockmaps[nblockmapsscatra] = mapthermo;

      // initialize map extractor associated with blocks of global system matrix
      blockmaps_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*dof_row_map(), blockmaps));

      // initialize map extractor associated with all degrees of freedom inside temperature field
      blockmapthermo_ = Teuchos::rcp(
          new Core::LinAlg::MultiMapExtractor(*thermo_field()->discretization()->dof_row_map(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(1, thermo_field()->dof_row_map())));

      // safety check
      blockmapthermo_->check_for_valid_map_extractor();
      if (scatra_field()->s2_i_meshtying())
      {
        // build block map for scatra interface by merging slave and master side for each block
        std::vector<Teuchos::RCP<const Epetra_Map>> partial_blockmapscatrainterface(
            nblockmapsscatra, Teuchos::null);
        for (int iblockmap = 0; iblockmap < nblockmapsscatra; ++iblockmap)
        {
          partial_blockmapscatrainterface.at(iblockmap) =
              Core::LinAlg::MultiMapExtractor::merge_maps(
                  {strategyscatra_->block_maps_slave().Map(iblockmap),
                      strategyscatra_->block_maps_master().Map(iblockmap)});
        }
        blockmapscatrainterface = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
            *interface_map_scatra, partial_blockmapscatrainterface));
        blockmapscatrainterface->check_for_valid_map_extractor();
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // safety check
  blockmaps_->check_for_valid_map_extractor();

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      // safety check
      if (!solver_->params().isSublist("AMGnxn Parameters"))
        FOUR_C_THROW(
            "Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *blockmaps_, *blockmaps_, 81, false, true));

      // feed AMGnxn block preconditioner with null space information for each block of global block
      // system matrix
      build_null_spaces();

      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      // safety check
      if (scatra_field()->matrix_type() != Core::LinAlg::MatrixType::sparse)
        FOUR_C_THROW("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dof_row_map(), 27, false, true));

      // feed AMG preconditioner with null space information associated with global system matrix if
      // applicable
      compute_null_space_if_necessary(solver_->params());

      break;
    }

    default:
    {
      FOUR_C_THROW("Type of global system matrix for scatra-thermo interaction not recognized!");
      break;
    }
  }

  // initialize scatra-thermo block and thermo-scatra block of global system matrix
  switch (scatra_field()->matrix_type())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      // initialize scatra-thermo blocks
      scatrathermoblockdomain_ = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, *scatra_field()->block_maps(), 81, false, true));
      scatrathermoblockinterface_ = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, *blockmapscatrainterface, 81, false, true));

      // initialize thermo-scatra blocks
      thermoscatrablockdomain_ = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *scatra_field()->block_maps(), *blockmapthermo_, 81, false, true));
      thermoscatrablockinterface_ = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *scatra_field()->block_maps(), *blockmapthermointerface, 81, false, true));

      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      // initialize scatra-thermo blocks
      scatrathermoblockdomain_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *scatra_field()->discretization()->dof_row_map(), 27, false, true));
      scatrathermoblockinterface_ =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*interface_map_scatra, 27, false, true));

      // initialize thermo-scatra blocks
      thermoscatrablockdomain_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *thermo_field()->discretization()->dof_row_map(), 27, false, true));
      thermoscatrablockinterface_ =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*interface_map_thermo, 27, false, true));

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // STI algorithm has non deforming mesh
  static constexpr bool isAle = false;

  // initialize OD evaluation strategy
  scatrathermooffdiagcoupling_ = STI::BuildScatraThermoOffDiagCoupling(
      strategyscatra_->coupling_type(), blockmapthermo_, blockmapthermointerface,
      blockmapthermointerfaceslave, maps_->Map(0), maps_->Map(1), interface_map_scatra,
      interface_map_thermo, isAle, strategyscatra_, strategythermo_, scatra_, thermo_);

  // instantiate appropriate equilibration class
  auto equilibration_method =
      std::vector<Core::LinAlg::EquilibrationMethod>(1, scatra_field()->equilibration_method());
  equilibration_ =
      Core::LinAlg::BuildEquilibration(matrixtype_, equilibration_method, maps_->full_map());
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::fd_check()
{
  // initial screen output
  if (get_comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR STI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(Core::LinAlg::CreateVector(*dof_row_map(), true));
  maps_->insert_vector(scatra_field()->phinp(), 0, statenp);
  maps_->insert_vector(thermo_field()->phinp(), 1, statenp);

  // make a copy of global state vector to undo perturbations later
  auto statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if (Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(systemmatrix_) !=
      Teuchos::null)
  {
    sysmat_original =
        (new Core::LinAlg::SparseMatrix(
             *(Teuchos::rcp_static_cast<Core::LinAlg::BlockSparseMatrixBase>(systemmatrix_)
                     ->merge())))
            ->epetra_matrix();
  }
  else
    FOUR_C_THROW("Global system matrix must be a block sparse matrix!");
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
    get_comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // fill global state vector with original state variables
    statenp->Update(1., *statenp_original, 0.);

    // impose perturbation
    if (statenp->Map().MyGID(colgid))
      if (statenp->SumIntoGlobalValue(colgid, 0, scatra_field()->fd_check_eps()))
        FOUR_C_THROW(
            "Perturbation could not be imposed on state vector for finite difference check!");
    scatra_field()->phinp()->Update(1., *maps_->extract_vector(statenp, 0), 0.);
    thermo_field()->phinp()->Update(1., *maps_->extract_vector(statenp, 1), 0.);

    // carry perturbation over to state vectors at intermediate time stages if necessary
    scatra_field()->compute_intermediate_values();
    thermo_field()->compute_intermediate_values();

    // calculate element right-hand side vector for perturbed state
    assemble_mat_and_rhs();

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
    for (int rowlid = 0; rowlid < dof_row_map()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if (rowgid < 0) FOUR_C_THROW("Invalid global ID of matrix row!");

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid, length, numentries, values.data(), indices.data());
      for (int ientry = 0; ientry < length; ++ientry)
      {
        if (sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval = -(*residual_)[rowlid] / scatra_field()->fd_check_eps() +
                           (*rhs_original)[rowlid] / scatra_field()->fd_check_eps();

      // confirm accuracy of first comparison
      if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        FOUR_C_THROW("Finite difference check involves values too close to numerical zero!");

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
      if (abs(relerr1) > scatra_field()->fd_check_tol())
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
        const double left = entry - (*rhs_original)[rowlid] / scatra_field()->fd_check_eps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / scatra_field()->fd_check_eps();

        // confirm accuracy of second comparison
        if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
          FOUR_C_THROW("Finite difference check involves values too close to numerical zero!");

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
        if (abs(relerr2) > scatra_field()->fd_check_tol())
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
  get_comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  get_comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  get_comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (get_comm().MyPID() == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      FOUR_C_THROW("Finite difference check failed for STI system matrix!");
    }
    else
    {
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
    }
  }

  // undo perturbations of state variables
  scatra_field()->phinp()->Update(1., *maps_->extract_vector(statenp_original, 0), 0.);
  scatra_field()->compute_intermediate_values();
  thermo_field()->phinp()->Update(1., *maps_->extract_vector(statenp_original, 1), 0.);
  thermo_field()->compute_intermediate_values();

  // recompute system matrix and right-hand side vector based on original state variables
  assemble_mat_and_rhs();
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::output_matrix_to_file(
    const Teuchos::RCP<const Core::LinAlg::SparseOperator> sparseoperator, const int precision,
    const double tolerance)
{
  // safety check
  if (!sparseoperator->filled()) FOUR_C_THROW("Sparse operator must be filled for output!");

  // extract communicator
  const Epetra_Comm& comm = sparseoperator->Comm();

  // determine whether sparse matrix or block sparse matrix should be output
  const auto sparsematrix =
      Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(sparseoperator);
  const auto blocksparsematrix =
      Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(sparseoperator);
  if (sparsematrix == Teuchos::null and blocksparsematrix == Teuchos::null)
    FOUR_C_THROW("Unknown type of sparse operator!");

  // extract row map
  const Epetra_Map& rowmap =
      sparsematrix != Teuchos::null ? sparsematrix->row_map() : blocksparsematrix->full_row_map();

  // safety check
  if (!rowmap.UniqueGIDs()) FOUR_C_THROW("Row map of matrix must be non-overlapping!");

  // copy global IDs of matrix rows stored on current processor into vector
  std::vector<int> myrowgids(rowmap.NumMyElements(), 0);
  int* myglobalelements = rowmap.MyGlobalElements();
  std::copy(myglobalelements, myglobalelements + rowmap.NumMyElements(), myrowgids.data());

  // communicate global IDs
  std::vector<int> rowgids(0, 0);
  Core::LinAlg::AllreduceVector(myrowgids, rowgids, comm);

  // retain communicated global IDs only on processor with ID 0
  if (comm.MyPID()) rowgids.clear();

  // create full row map on processor with ID 0
  const Epetra_Map fullrowmap(
      -1, static_cast<int>(rowgids.size()), rowgids.size() ? rowgids.data() : nullptr, 0, comm);

  // import matrix to processor with ID 0
  Epetra_CrsMatrix crsmatrix(Copy, fullrowmap, 0);
  if (sparsematrix != Teuchos::null)
  {
    if (crsmatrix.Import(*sparsematrix->epetra_matrix(), Epetra_Import(fullrowmap, rowmap), Insert))
      FOUR_C_THROW("Matrix import failed!");
  }
  else
  {
    for (int i = 0; i < blocksparsematrix->rows(); ++i)
    {
      for (int j = 0; j < blocksparsematrix->cols(); ++j)
      {
        if (crsmatrix.Import(*blocksparsematrix->matrix(i, j).epetra_matrix(),
                Epetra_Import(fullrowmap, blocksparsematrix->range_map(i)), Insert))
          FOUR_C_THROW("Matrix import failed!");
      }
    }
  }

  // let processor with ID 0 output matrix to file
  if (comm.MyPID() == 0)
  {
    // set file name
    std::ostringstream nproc;
    nproc << comm.NumProc();
    const std::string filename(Global::Problem::instance()->output_control_file()->file_name() +
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
        FOUR_C_THROW("Cannot extract matrix row with global ID %d!", rowgid);

      // sort entries in current matrix row in ascending order of column global ID via map
      std::map<int, double> entries;

      // loop over all entries in current matrix row
      for (int j = 0; j < numentries; ++j)
        // add current matrix entry to map
        entries[indices[j]] = values[j];

      // loop over all sorted entries in current matrix row
      for (auto& entrie : entries)
      {
        // write current matrix entry to file
        if (std::abs(entrie.second) > tolerance)
          file << rowgid << "," << entrie.first << "," << entrie.second << std::endl;
      }
    }

    // close file
    file.close();
  }

  // wait until output is complete
  comm.Barrier();

  // throw error to abort simulation for debugging
  FOUR_C_THROW("Matrix was output to *.csv file!");
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::output_vector_to_file(
    const Epetra_MultiVector& vector, const int precision, const double tolerance)
{
  // extract communicator
  const Epetra_Comm& comm = vector.Comm();

  // extract vector map
  const Epetra_BlockMap& map = vector.Map();

  // safety check
  if (!map.UniqueGIDs())
    FOUR_C_THROW(
        "Vector output to *.csv file currently only works for non-overlapping vector maps!");

  // copy global IDs of vector components stored on current processor into vector
  std::vector<int> mygids(map.NumMyElements(), 0);
  int* myglobalelements = map.MyGlobalElements();
  std::copy(myglobalelements, myglobalelements + map.NumMyElements(), mygids.data());

  // communicate global IDs
  std::vector<int> gids(0, 0);
  Core::LinAlg::AllreduceVector(mygids, gids, comm);

  // retain communicated global IDs only on processor with ID 0
  if (comm.MyPID()) gids.clear();

  // create full vector map on processor with ID 0
  const Epetra_Map fullmap(
      -1, static_cast<int>(gids.size()), gids.size() ? gids.data() : nullptr, 0, comm);

  // export vector to processor with ID 0
  Epetra_MultiVector fullvector(fullmap, vector.NumVectors(), true);
  Core::LinAlg::export_to(vector, fullvector);

  // let processor with ID 0 output vector to file
  if (comm.MyPID() == 0)
  {
    // set file name
    std::ostringstream nproc;
    nproc << comm.NumProc();
    const std::string filename(Global::Problem::instance()->output_control_file()->file_name() +
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
  FOUR_C_THROW("Vector was output to *.csv file!");
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::assemble_mat_and_rhs()
{
  // pass thermo degrees of freedom to scatra discretization
  transfer_thermo_to_scatra(thermo_field()->phiafnp());

  // build system matrix and residual for scatra field
  scatra_field()->prepare_linear_solve();

  // pass scatra degrees of freedom to thermo discretization
  transfer_scatra_to_thermo(scatra_field()->phiafnp());

  // build system matrix and residual for thermo field
  thermo_field()->prepare_linear_solve();

  // evaluate scatra-thermo-OD coupling. Contributions from domain
  scatrathermooffdiagcoupling_->evaluate_off_diag_block_scatra_thermo_domain(
      scatrathermoblockdomain_);

  // evaluate scatra-thermo-OD coupling. Contributions from interface
  if (scatra_field()->s2_i_meshtying())
    scatrathermooffdiagcoupling_->evaluate_off_diag_block_scatra_thermo_interface(
        scatrathermoblockinterface_);

  // evaluate thermo-scatra-OD coupling. Contributions from domain
  scatrathermooffdiagcoupling_->evaluate_off_diag_block_thermo_scatra_domain(
      thermoscatrablockdomain_);

  // evaluate thermo-scatra-OD coupling. Contributions from interface
  if (scatra_field()->s2_i_meshtying())
    scatrathermooffdiagcoupling_->evaluate_off_diag_block_thermo_scatra_interface(
        thermoscatrablockinterface_);

  // OD blocks containing assembled domain and interface contributions
  Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermo_domain_interface(Teuchos::null);
  Teuchos::RCP<Core::LinAlg::SparseOperator> thermoscatra_domain_interface(Teuchos::null);

  // assemble interface and domain contributions of OD blocks
  assemble_domain_interface_off_diag(scatrathermo_domain_interface, thermoscatra_domain_interface);

  // apply Dirichlet BCs to OD blocks
  apply_dirichlet_off_diag(scatrathermo_domain_interface, thermoscatra_domain_interface);

  // build global system matrix
  switch (matrixtype_)
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      // check global system matrix
      Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blocksystemmatrix =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(systemmatrix_);
      if (blocksystemmatrix == Teuchos::null) FOUR_C_THROW("System matrix is not a block matrix!");

      switch (scatra_field()->matrix_type())
      {
        case Core::LinAlg::MatrixType::block_condition:
        {
          // extract number of matrix row or column blocks associated with scalar transport field
          const int nblockmapsscatra = static_cast<int>(scatra_field()->block_maps()->num_maps());

          // construct global system matrix by assigning matrix blocks
          for (int iblock = 0; iblock < nblockmapsscatra; ++iblock)
          {
            for (int jblock = 0; jblock < nblockmapsscatra; ++jblock)
              blocksystemmatrix->assign(iblock, jblock, Core::LinAlg::View,
                  scatra_field()->block_system_matrix()->matrix(iblock, jblock));

            // perform second condensation before assigning matrix blocks
            if (condensationthermo_)
            {
              const auto& scatrathermoblock =
                  Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(
                      scatrathermo_domain_interface)
                      ->matrix(iblock, 0);
              Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                  scatrathermoblock.range_map(), *maps_->Map(1), 1.0, nullptr, nullptr,
                  blocksystemmatrix->matrix(iblock, nblockmapsscatra));

              switch (strategyscatra_->coupling_type())
              {
                case Inpar::S2I::coupling_matching_nodes:
                {
                  Core::Adapter::CouplingSlaveConverter converter(*icoupthermo_);
                  Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                      scatrathermoblock.range_map(), *strategythermo_->interface_maps()->Map(1),
                      1.0, nullptr, &converter, blocksystemmatrix->matrix(iblock, nblockmapsscatra),
                      true, true);
                  break;
                }
                case Inpar::S2I::coupling_mortar_standard:
                {
                  // initialize temporary matrix for slave-side columns of current matrix block
                  Core::LinAlg::SparseMatrix scatrathermocolsslave(scatrathermoblock.row_map(), 81);

                  // fill temporary matrix for slave-side columns of current matrix block
                  Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                      scatrathermoblock.range_map(), *strategythermo_->interface_maps()->Map(1),
                      1.0, nullptr, nullptr, scatrathermocolsslave);

                  // finalize temporary matrix for slave-side columns of current matrix block
                  scatrathermocolsslave.complete(
                      *strategythermo_->interface_maps()->Map(1), scatrathermoblock.row_map());

                  // transform and assemble temporary matrix for slave-side columns of current
                  // matrix block
                  blocksystemmatrix->matrix(iblock, nblockmapsscatra)
                      .add(*Core::LinAlg::MLMultiply(
                               scatrathermocolsslave, *strategythermo_->p(), true),
                          false, 1.0, 1.0);

                  break;
                }
                default:
                {
                  FOUR_C_THROW("Invalid type of scatra-scatra interface coupling!");
                  break;
                }
              }

              const auto& thermoscatrablock =
                  Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(
                      thermoscatra_domain_interface)
                      ->matrix(0, iblock);
              Core::LinAlg::MatrixLogicalSplitAndTransform()(thermoscatrablock, *maps_->Map(1),
                  thermoscatrablock.domain_map(), 1.0, nullptr, nullptr,
                  blocksystemmatrix->matrix(nblockmapsscatra, iblock));
            }

            // assign matrix blocks directly
            else
            {
              blocksystemmatrix->assign(iblock, nblockmapsscatra, Core::LinAlg::View,
                  Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(
                      scatrathermo_domain_interface)
                      ->matrix(iblock, 0));

              blocksystemmatrix->assign(nblockmapsscatra, iblock, Core::LinAlg::View,
                  Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(
                      thermoscatra_domain_interface)
                      ->matrix(0, iblock));
            }
          }

          // perform second condensation before assigning thermo-thermo matrix block
          if (condensationthermo_)
          {
            Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                *maps_->Map(1), *maps_->Map(1), 1.0, nullptr, nullptr,
                blocksystemmatrix->matrix(nblockmapsscatra, nblockmapsscatra));

            switch (strategyscatra_->coupling_type())
            {
              case Inpar::S2I::coupling_matching_nodes:
              {
                Core::Adapter::CouplingSlaveConverter converter(*icoupthermo_);
                Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                    *maps_->Map(1), *strategythermo_->interface_maps()->Map(1), 1.0, nullptr,
                    &converter, blocksystemmatrix->matrix(nblockmapsscatra, nblockmapsscatra), true,
                    true);
                break;
              }
              case Inpar::S2I::coupling_mortar_standard:
              {
                // initialize temporary matrix for slave-side columns of thermo-thermo matrix block
                Core::LinAlg::SparseMatrix thermothermocolsslave(*maps_->Map(1), 81);

                // fill temporary matrix for slave-side columns of thermo-thermo matrix block
                Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                    *maps_->Map(1), *strategythermo_->interface_maps()->Map(1), 1.0, nullptr,
                    nullptr, thermothermocolsslave);

                // finalize temporary matrix for slave-side columns of thermo-thermo matrix block
                thermothermocolsslave.complete(
                    *strategythermo_->interface_maps()->Map(1), *maps_->Map(1));

                // transform and assemble temporary matrix for slave-side columns of thermo-thermo
                // matrix block
                blocksystemmatrix->matrix(nblockmapsscatra, nblockmapsscatra)
                    .add(*Core::LinAlg::MLMultiply(
                             thermothermocolsslave, *strategythermo_->p(), true),
                        false, 1.0, 1.0);

                break;
              }
              default:
              {
                FOUR_C_THROW("Invalid type of scatra-scatra interface coupling!");
                break;
              }
            }
          }

          // assign thermo-thermo matrix block directly
          else
            blocksystemmatrix->assign(nblockmapsscatra, nblockmapsscatra, Core::LinAlg::View,
                *thermo_field()->system_matrix());

          break;
        }

        case Core::LinAlg::MatrixType::sparse:
        {
          // construct global system matrix by assigning matrix blocks
          blocksystemmatrix->assign(0, 0, Core::LinAlg::View, *scatra_field()->system_matrix());

          // perform second condensation before assigning matrix blocks
          if (condensationthermo_)
          {
            const auto& scatrathermoblock =
                *Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(
                    scatrathermo_domain_interface);
            Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                scatrathermoblock.range_map(), *maps_->Map(1), 1.0, nullptr, nullptr,
                blocksystemmatrix->matrix(0, 1));

            Core::LinAlg::MatrixLogicalSplitAndTransform()(
                *Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(
                    thermoscatra_domain_interface),
                *maps_->Map(1), thermoscatra_domain_interface->domain_map(), 1.0, nullptr, nullptr,
                blocksystemmatrix->matrix(1, 0));

            Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                *maps_->Map(1), *maps_->Map(1), 1.0, nullptr, nullptr,
                blocksystemmatrix->matrix(1, 1));

            switch (strategyscatra_->coupling_type())
            {
              case Inpar::S2I::coupling_matching_nodes:
              {
                Core::Adapter::CouplingSlaveConverter converter(*icoupthermo_);
                Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                    scatrathermoblock.range_map(), *strategythermo_->interface_maps()->Map(1), 1.0,
                    nullptr, &converter, blocksystemmatrix->matrix(0, 1), true, true);

                Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                    *maps_->Map(1), *strategythermo_->interface_maps()->Map(1), 1.0, nullptr,
                    &converter, blocksystemmatrix->matrix(1, 1), true, true);

                break;
              }

              case Inpar::S2I::coupling_mortar_standard:
              {
                // initialize temporary matrix for slave-side columns of scatra-thermo matrix block
                Core::LinAlg::SparseMatrix scatrathermocolsslave(
                    *scatra_field()->discretization()->dof_row_map(), 81);

                // fill temporary matrix for slave-side columns of scatra-thermo matrix block
                Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                    scatrathermoblock.range_map(), *strategythermo_->interface_maps()->Map(1), 1.0,
                    nullptr, nullptr, scatrathermocolsslave);

                // finalize temporary matrix for slave-side columns of scatra-thermo matrix block
                scatrathermocolsslave.complete(*strategythermo_->interface_maps()->Map(1),
                    *scatra_field()->discretization()->dof_row_map());

                // transform and assemble temporary matrix for slave-side columns of scatra-thermo
                // matrix block
                blocksystemmatrix->matrix(0, 1).add(
                    *Core::LinAlg::MLMultiply(scatrathermocolsslave, *strategythermo_->p(), true),
                    false, 1.0, 1.0);

                // initialize temporary matrix for slave-side columns of thermo-thermo matrix block
                Core::LinAlg::SparseMatrix thermothermocolsslave(*maps_->Map(1), 81);

                // fill temporary matrix for slave-side columns of thermo-thermo matrix block
                Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                    *maps_->Map(1), *strategythermo_->interface_maps()->Map(1), 1.0, nullptr,
                    nullptr, thermothermocolsslave);

                // finalize temporary matrix for slave-side columns of thermo-thermo matrix block
                thermothermocolsslave.complete(
                    *strategythermo_->interface_maps()->Map(1), *maps_->Map(1));

                // transform and assemble temporary matrix for slave-side columns of thermo-thermo
                // matrix block
                blocksystemmatrix->matrix(1, 1).add(
                    *Core::LinAlg::MLMultiply(thermothermocolsslave, *strategythermo_->p(), true),
                    false, 1.0, 1.0);

                break;
              }

              default:
              {
                FOUR_C_THROW("Invalid type of scatra-scatra interface coupling!");
                break;
              }
            }
          }

          // assign matrix blocks directly
          else
          {
            blocksystemmatrix->assign(0, 1, Core::LinAlg::View,
                *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(
                    scatrathermo_domain_interface));
            blocksystemmatrix->assign(1, 0, Core::LinAlg::View,
                *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(
                    thermoscatra_domain_interface));
            blocksystemmatrix->assign(1, 1, Core::LinAlg::View, *thermo_field()->system_matrix());
          }

          break;
        }

        default:
        {
          FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
          break;
        }
        case Core::LinAlg::MatrixType::undefined:
        case Core::LinAlg::MatrixType::block_field:
        case Core::LinAlg::MatrixType::block_condition_dof:
          break;
      }

      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      // check global system matrix
      auto systemmatrix = Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(systemmatrix_);
      if (systemmatrix == Teuchos::null) FOUR_C_THROW("System matrix is not a sparse matrix!");

      // construct global system matrix by adding matrix blocks
      systemmatrix->add(*scatra_field()->system_matrix(), false, 1.0, 0.0);

      // perform second condensation before adding matrix blocks
      if (condensationthermo_)
      {
        const auto& scatrathermoblock =
            *Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(
                scatrathermo_domain_interface);
        Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
            scatrathermoblock.range_map(), *maps_->Map(1), 1.0, nullptr, nullptr, *systemmatrix,
            true, true);

        Core::LinAlg::MatrixLogicalSplitAndTransform()(
            *Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(
                thermoscatra_domain_interface),
            *maps_->Map(1), thermoscatra_domain_interface->domain_map(), 1.0, nullptr, nullptr,
            *systemmatrix, true, true);

        Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
            *maps_->Map(1), *maps_->Map(1), 1.0, nullptr, nullptr, *systemmatrix, true, true);

        switch (strategyscatra_->coupling_type())
        {
          case Inpar::S2I::coupling_matching_nodes:
          {
            Core::Adapter::CouplingSlaveConverter converter(*icoupthermo_);
            Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                scatrathermoblock.range_map(), *strategythermo_->interface_maps()->Map(1), 1.0,
                nullptr, &converter, *systemmatrix, true, true);

            Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                *maps_->Map(1), *strategythermo_->interface_maps()->Map(1), 1.0, nullptr,
                &converter, *systemmatrix, true, true);

            break;
          }

          case Inpar::S2I::coupling_mortar_standard:
          {
            // initialize temporary matrix for slave-side columns of global system matrix
            Core::LinAlg::SparseMatrix systemmatrixcolsslave(*dof_row_map(), 81);

            // fill temporary matrix for slave-side columns of global system matrix
            Core::LinAlg::MatrixLogicalSplitAndTransform()(scatrathermoblock,
                scatrathermoblock.range_map(), *strategythermo_->interface_maps()->Map(1), 1.0,
                nullptr, nullptr, systemmatrixcolsslave);

            Core::LinAlg::MatrixLogicalSplitAndTransform()(*thermo_field()->system_matrix(),
                *maps_->Map(1), *strategythermo_->interface_maps()->Map(1), 1.0, nullptr, nullptr,
                systemmatrixcolsslave, true, true);

            // finalize temporary matrix for slave-side columns of global system matrix
            systemmatrixcolsslave.complete(
                *strategythermo_->interface_maps()->Map(1), *dof_row_map());

            // transform and assemble temporary matrix for slave-side columns of global system
            // matrix
            systemmatrix->add(
                *Core::LinAlg::MLMultiply(systemmatrixcolsslave, *strategythermo_->p(), true),
                false, 1.0, 1.0);

            break;
          }

          default:
          {
            FOUR_C_THROW("Invalid type of scatra-scatra interface coupling!");
            break;
          }
        }
      }

      // add matrix blocks directly
      else
      {
        systemmatrix->add(*scatrathermo_domain_interface, false, 1.0, 1.0);
        systemmatrix->add(*thermoscatra_domain_interface, false, 1.0, 1.0);
        systemmatrix->add(*thermo_field()->system_matrix(), false, 1.0, 1.0);
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Type of global system matrix for scatra-thermo interaction not recognized!");
      break;
    }
  }

  // finalize global system matrix
  systemmatrix_->complete();

  // create full monolithic right-hand side vector
  maps_->insert_vector(scatra_field()->residual(), 0, residual_);
  Teuchos::RCP<Epetra_Vector> thermoresidual(Teuchos::null);
  if (condensationthermo_)
  {
    thermoresidual = Teuchos::rcp(new Epetra_Vector(*maps_->Map(1)));
    Core::LinAlg::export_to(*thermo_field()->residual(), *thermoresidual);
  }
  else
    thermoresidual = thermo_field()->residual();
  maps_->insert_vector(thermoresidual, 1, residual_);
}  // STI::Monolithic::assemble_mat_and_rhs()


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void STI::Monolithic::build_null_spaces() const
{
  switch (scatra_field()->matrix_type())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      scatra_field()->build_block_null_spaces(solver_, 0);
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sublists to trigger null space
      // computation
      Teuchos::ParameterList& blocksmootherparams = solver_->params().sublist("Inverse1");
      blocksmootherparams.sublist("Belos Parameters");
      blocksmootherparams.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      scatra_field()->discretization()->compute_null_space_if_necessary(blocksmootherparams);

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with temperature field as string
  std::stringstream iblockstr;
  iblockstr << blockmaps_->num_maps();

  // equip smoother for thermo matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Belos Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for thermo matrix block with null space associated with all degrees of freedom
  // on thermo discretization
  thermo_field()->discretization()->compute_null_space_if_necessary(blocksmootherparams);

  // reduce full null space to match degrees of freedom associated with thermo matrix block if
  // necessary
  if (condensationthermo_)
    Core::LinearSolver::Parameters::fix_null_space("Block " + iblockstr.str(),
        *thermo_field()->discretization()->dof_row_map(), *maps_->Map(1), blocksmootherparams);
}  // STI::Monolithic::build_block_null_spaces

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::compute_null_space_if_necessary(Teuchos::ParameterList& solverparams) const
{
  // compute vector-based null space information for ML preconditioner
  if (solverparams.isSublist("ML Parameters"))
  {
    // extract parameter list for ML preconditioner
    Teuchos::ParameterList& mllist = solverparams.sublist("ML Parameters", true);

    // determine null space dimension
    const int numdofpernode_scatra = scatra_field()->num_dof_per_node();
    const int numdofpernode_thermo = thermo_field()->num_dof_per_node();
    const int dimns = numdofpernode_scatra + numdofpernode_thermo;

    // allocate vector for null space information
    const auto ns =
        Teuchos::rcp(new std::vector<double>(dimns * dof_row_map()->NumMyElements(), 0.));

    // compute null space modes associated with scatra field
    const Core::FE::Discretization& scatradis = *scatra_field()->discretization();
    std::vector<double*> modes_scatra(numdofpernode_scatra);
    for (int i = 0; i < numdofpernode_scatra; ++i)
      modes_scatra[i] = &((*ns)[i * dof_row_map()->NumMyElements()]);
    for (int i = 0; i < scatradis.num_my_row_nodes(); ++i)
    {
      const int lid = dof_row_map()->LID(scatradis.dof(0, scatradis.l_row_node(i), 0));
      if (lid < 0) FOUR_C_THROW("Cannot find scatra degree of freedom!");
      for (int j = 0; j < numdofpernode_scatra; ++j) modes_scatra[j][lid + j] = 1.;
    }

    // compute null space modes associated with thermo field
    const Core::FE::Discretization& thermodis = *thermo_field()->discretization();
    std::vector<double*> modes_thermo(numdofpernode_thermo);
    for (int i = 0; i < numdofpernode_thermo; ++i)
      modes_thermo[i] = &((*ns)[(numdofpernode_scatra + i) * dof_row_map()->NumMyElements()]);
    for (int i = 0; i < thermodis.num_my_row_nodes(); ++i)
    {
      const int lid = dof_row_map()->LID(thermodis.dof(0, thermodis.l_row_node(i), 0));
      if (lid < 0) FOUR_C_THROW("Cannot find thermo degree of freedom!");
      for (int j = 0; j < numdofpernode_thermo; ++j) modes_thermo[j][lid + j] = 1.;
    }

    // fill parameter list
    mllist.set("PDE equations", dimns);
    mllist.set("null space: dimension", dimns);
    mllist.set("null space: type", "pre-computed");
    mllist.set("null space: add default vectors", false);

    Teuchos::RCP<Epetra_MultiVector> nullspace =
        Teuchos::rcp(new Epetra_MultiVector(dof_row_map().operator*(), dimns, true));
    Core::LinAlg::StdVectorToEpetraMultiVector(*ns, nullspace, dimns);

    mllist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
    mllist.set("null space: vectors", nullspace->Values());
    mllist.set("ML validate parameter list", false);
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

    Teuchos::RCP<Epetra_MultiVector> nullspace =
        Teuchos::rcp(new Epetra_MultiVector(dof_row_map().operator*(), 1, true));
    nullspace->PutScalar(1.0);

    mllist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
    mllist.set("null space: vectors", nullspace->Values());
    mllist.set("ML validate parameter list", false);

    Teuchos::RCP<Epetra_MultiVector> coordinates =
        scatra_field()->discretization()->build_node_coordinates();

    mllist.set<Teuchos::RCP<Epetra_MultiVector>>("Coordinates", coordinates);
  }
}  // STI::Monolithic::compute_null_space_if_necessary

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& STI::Monolithic::dof_row_map() const
{
  return maps_->full_map();
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
bool STI::Monolithic::exit_newton_raphson()
{
  // initialize exit flag
  bool exit(false);

  // perform Newton-Raphson convergence check depending on type of scalar transport
  switch (
      Teuchos::getIntegralValue<Inpar::STI::ScaTraTimIntType>(*stiparameters_, "SCATRATIMINTTYPE"))
  {
    case Inpar::STI::ScaTraTimIntType::elch:
    {
      // compute L2 norm of concentration state vector
      double concdofnorm(0.);
      scatra_field()
          ->splitter()
          ->extract_other_vector(*scatra_field()->phinp())
          ->Norm2(&concdofnorm);

      // compute L2 norm of concentration residual vector
      double concresnorm(0.);
      scatra_field()
          ->splitter()
          ->extract_other_vector(*maps_->extract_vector(residual_, 0))
          ->Norm2(&concresnorm);

      // compute L2 norm of concentration increment vector
      double concincnorm(0.);
      scatra_field()
          ->splitter()
          ->extract_other_vector(*maps_->extract_vector(increment_, 0))
          ->Norm2(&concincnorm);

      // compute L2 norm of potential state vector
      double potdofnorm(0.);
      scatra_field()->splitter()->extract_cond_vector(*scatra_field()->phinp())->Norm2(&potdofnorm);

      // compute L2 norm of potential residual vector
      double potresnorm(0.);
      scatra_field()
          ->splitter()
          ->extract_cond_vector(*maps_->extract_vector(residual_, 0))
          ->Norm2(&potresnorm);

      // compute L2 norm of potential increment vector
      double potincnorm(0.);
      scatra_field()
          ->splitter()
          ->extract_cond_vector(*maps_->extract_vector(increment_, 0))
          ->Norm2(&potincnorm);

      // compute L2 norm of thermo state vector
      double thermodofnorm(0.);
      thermo_field()->phinp()->Norm2(&thermodofnorm);

      // compute L2 norm of thermo residual vector
      double thermoresnorm(0.);
      maps_->extract_vector(residual_, 1)->Norm2(&thermoresnorm);

      // compute L2 norm of thermo increment vector
      double thermoincnorm(0.);
      maps_->extract_vector(increment_, 1)->Norm2(&thermoincnorm);

      // safety checks
      if (std::isnan(concdofnorm) or std::isnan(concresnorm) or std::isnan(concincnorm) or
          std::isnan(potdofnorm) or std::isnan(potresnorm) or std::isnan(potincnorm) or
          std::isnan(thermodofnorm) or std::isnan(thermoresnorm) or std::isnan(thermoincnorm))
        FOUR_C_THROW("Vector norm is not a number!");
      if (std::isinf(concdofnorm) or std::isinf(concresnorm) or std::isinf(concincnorm) or
          std::isinf(potdofnorm) or std::isinf(potresnorm) or std::isinf(potincnorm) or
          std::isinf(thermodofnorm) or std::isinf(thermoresnorm) or std::isinf(thermoincnorm))
        FOUR_C_THROW("Vector norm is infinity!");

      // prevent division by zero
      if (concdofnorm < 1.e-10) concdofnorm = 1.e-10;
      if (potdofnorm < 1.e-10) potdofnorm = 1.e-10;
      if (thermodofnorm < 1.e-10) thermodofnorm = 1.e-10;

      // first Newton-Raphson iteration
      if (iter_ == 1)
      {
        if (get_comm().MyPID() == 0)
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
        if (get_comm().MyPID() == 0)
        {
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
        }

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
        if (get_comm().MyPID() == 0)
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
      if (exit and get_comm().MyPID() == 0)
      {
        std::cout << "+------------+-------------------+--------------+--------------+-------------"
                     "-+--------------+--------------+--------------+"
                  << std::endl
                  << std::endl;
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Type of scalar transport not yet available!");
      break;
    }
  }

  return exit;
}  // STI::Monolithic::exit_newton_raphson()

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::prepare_time_step()
{
  // call base class routine
  Algorithm::prepare_time_step();

  // print time step information to screen
  scatra_field()->print_time_step_info();
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::Monolithic::solve()
{
  // initialize counter for Newton-Raphson iterations
  iter_ = 0;

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    ++iter_;

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->wallTime();

    // assemble global system of equations
    assemble_mat_and_rhs();

    // determine time needed for evaluating elements and assembling global system of equations,
    // and take maximum over all processors via communication
    double mydtele = timer_->wallTime() - time;
    get_comm().MaxAll(&mydtele, &dtele_, 1);

    // safety check
    if (!systemmatrix_->filled())
      FOUR_C_THROW("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if (scatra_field()->fd_check_type() == Inpar::ScaTra::fdcheck_global) fd_check();

    // check termination criterion for Newton-Raphson iteration
    if (exit_newton_raphson()) break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    time = timer_->wallTime();

    // equilibrate global system of equations if necessary
    equilibration_->equilibrate_system(systemmatrix_, residual_, blockmaps_);

    // solve global system of equations
    // Dirichlet boundary conditions have already been applied to global system of equations
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = iter_ == 1;
    solver_->solve(systemmatrix_->epetra_operator(), increment_, residual_, solver_params);

    equilibration_->unequilibrate_increment(increment_);

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->wallTime() - time;
    get_comm().MaxAll(&mydtsolve, &dtsolve_, 1);

    // output performance statistics associated with linear solver into text file if applicable
    if (Core::UTILS::IntegralValue<int>(*fieldparameters_, "OUTPUTLINSOLVERSTATS"))
      scatra_field()->output_lin_solver_stats(*solver_, dtsolve_, step(), static_cast<int>(iter_),
          residual_->Map().NumGlobalElements());

    // update scatra field
    scatra_field()->update_iter(maps_->extract_vector(increment_, 0));
    scatra_field()->compute_intermediate_values();

    // update thermo field
    Teuchos::RCP<Epetra_Vector> thermoincrement(Teuchos::null);
    if (condensationthermo_)
    {
      thermoincrement =
          Teuchos::rcp(new Epetra_Vector(*thermo_field()->discretization()->dof_row_map()));
      Core::LinAlg::export_to(*maps_->extract_vector(increment_, 1), *thermoincrement);
      const Teuchos::RCP<const Epetra_Vector> masterincrement =
          strategythermo_->interface_maps()->extract_vector(*thermoincrement, 2);
      const Teuchos::RCP<Epetra_Vector> slaveincrement =
          Core::LinAlg::CreateVector(*strategythermo_->interface_maps()->Map(1));
      switch (strategyscatra_->coupling_type())
      {
        case Inpar::S2I::coupling_matching_nodes:
        {
          icoupthermo_->master_to_slave(masterincrement, slaveincrement);
          break;
        }
        case Inpar::S2I::coupling_mortar_standard:
        {
          strategythermo_->p()->multiply(false, *masterincrement, *slaveincrement);
          break;
        }
        default:
        {
          FOUR_C_THROW("Invalid type of scatra-scatra interface coupling!");
          break;
        }
      }
      strategythermo_->interface_maps()->insert_vector(slaveincrement, 1, thermoincrement);
    }
    else
      thermoincrement = maps_->extract_vector(increment_, 1);
    thermo_field()->update_iter(thermoincrement);
    thermo_field()->compute_intermediate_values();
  }  // Newton-Raphson iteration
}  // STI::Monolithic::Solve

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::Monolithic::apply_dirichlet_off_diag(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& scatrathermo_domain_interface,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& thermoscatra_domain_interface)
{
  // apply Dirichlet boundary conditions to scatra-thermo matrix block
  scatrathermo_domain_interface->apply_dirichlet(*scatra_field()->dirich_maps()->cond_map(), false);

  // apply Dirichlet boundary conditions to scatra-thermo matrix block
  thermoscatra_domain_interface->apply_dirichlet(*thermo_field()->dirich_maps()->cond_map(), false);

  // zero out slave-side rows of thermo-scatra matrix block after having added them to the
  // corresponding master-side rows to finalize condensation of slave-side thermo dofs
  if (thermo_field()->s2_i_meshtying())
  {
    switch (strategythermo_->coupling_type())
    {
      case Inpar::S2I::coupling_matching_nodes:
      {
        if (!thermo_field()->discretization()->get_condition("PointCoupling"))
          thermoscatra_domain_interface->apply_dirichlet(*icoupthermo_->slave_dof_map(), false);
        break;
      }

      case Inpar::S2I::coupling_mortar_condensed_bubnov:
      {
        thermoscatra_domain_interface->apply_dirichlet(
            *(strategythermo_->interface_maps()->Map(1)), false);
        break;
      }

      default:
      {
        FOUR_C_THROW("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }
}  // STI::Monolithic::ApplyDirichletOD

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::Monolithic::assemble_domain_interface_off_diag(
    Teuchos::RCP<Core::LinAlg::SparseOperator>& scatrathermo_domain_interface,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& thermoscatra_domain_interface)
{
  // initialize scatra-thermo blocks
  switch (scatra_field()->matrix_type())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      scatrathermo_domain_interface = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, *scatra_field()->block_maps(), 81, false, true));
      thermoscatra_domain_interface = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *scatra_field()->block_maps(), *blockmapthermo_, 81, false, true));

      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      scatrathermo_domain_interface = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *scatra_field()->discretization()->dof_row_map(), 27, false, true));
      thermoscatra_domain_interface = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *thermo_field()->discretization()->dof_row_map(), 27, false, true));
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown matrix type");
      break;
    }
  }

  // add domain and interface contibutions
  scatrathermo_domain_interface->add(*scatrathermoblockdomain_, false, 1.0, 1.0);
  scatrathermo_domain_interface->add(*scatrathermoblockinterface_, false, 1.0, 1.0);
  thermoscatra_domain_interface->add(*thermoscatrablockdomain_, false, 1.0, 1.0);
  thermoscatra_domain_interface->add(*thermoscatrablockinterface_, false, 1.0, 1.0);

  if (strategythermo_->coupling_type() == Inpar::S2I::coupling_matching_nodes)
  {  // standard meshtying algorithm with Lagrange multipliers condensed out
    if (!thermo_field()->discretization()->get_condition("PointCoupling"))
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
      if (scatra_field()->matrix_type() == Core::LinAlg::MatrixType::block_condition and
          step() == 1 and iter_ == 1)
        thermoscatra_domain_interface->complete();

      // loop over all thermo-scatra matrix blocks
      for (int iblock = 0; iblock < blockmaps_->num_maps() - 1; ++iblock)
      {
        // initialize temporary matrix for slave-side rows of current thermo-scatra matrix
        // block
        Core::LinAlg::SparseMatrix thermoscatrarowsslave(
            *icoupthermo_->slave_dof_map(), 27, false, true);

        // extract current thermo-scatra matrix block
        auto& thermoscatrablock =
            scatra_field()->matrix_type() == Core::LinAlg::MatrixType::block_condition
                ? Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(
                      thermoscatra_domain_interface)
                      ->matrix(0, iblock)
                : *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(
                      thermoscatra_domain_interface);

        // extract slave-side rows of thermo-scatra matrix block into temporary matrix
        ScaTra::MeshtyingStrategyS2I::extract_matrix_rows(
            thermoscatrablock, thermoscatrarowsslave, *icoupthermo_->slave_dof_map());

        // finalize temporary matrix with slave-side rows of thermo-scatra matrix block
        thermoscatrarowsslave.complete(*maps_->Map(0), *icoupthermo_->slave_dof_map());

        // undo Complete() from above before performing subsequent matrix row transformation
        if (scatra_field()->matrix_type() == Core::LinAlg::MatrixType::block_condition and
            step() == 1 and iter_ == 1)
          thermoscatrablock.un_complete();

        // add slave-side rows of thermo-scatra matrix block to corresponding slave-side rows
        const Teuchos::RCP<Core::LinAlg::MatrixRowTransform> islavetomasterrowtransformthermood =
            scatra_field()->matrix_type() == Core::LinAlg::MatrixType::block_condition
                ? Teuchos::rcp(new Core::LinAlg::MatrixRowTransform())
                : islavetomasterrowtransformthermood_;
        (*islavetomasterrowtransformthermood)(thermoscatrarowsslave, 1.,
            Core::Adapter::CouplingSlaveConverter(*icoupthermo_), thermoscatrablock, true);
      }
    }
  }
  else if (strategythermo_->coupling_type() == Inpar::S2I::coupling_mortar_condensed_bubnov)
  {
    // standard meshtying algorithm with Lagrange multipliers condensed out
    // loop over all thermo-scatra matrix blocks
    for (int iblock = 0; iblock < blockmaps_->num_maps() - 1; ++iblock)
    {
      // initialize temporary matrix for slave-side rows of current thermo-scatra matrix block
      Core::LinAlg::SparseMatrix thermoscatrarowsslave(
          *strategythermo_->interface_maps()->Map(1), 27, false, true);

      // extract current thermo-scatra matrix block
      auto& thermoscatrablock =
          scatra_field()->matrix_type() == Core::LinAlg::MatrixType::block_condition
              ? Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(
                    thermoscatra_domain_interface)
                    ->matrix(0, iblock)
              : *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(
                    thermoscatra_domain_interface);

      // extract slave-side rows of thermo-scatra matrix block into temporary matrix
      ScaTra::MeshtyingStrategyS2I::extract_matrix_rows(
          thermoscatrablock, thermoscatrarowsslave, *strategythermo_->interface_maps()->Map(1));

      // finalize temporary matrix with slave-side rows of thermo-scatra matrix block
      thermoscatrarowsslave.complete(*maps_->Map(0), *strategythermo_->interface_maps()->Map(1));

      // add projected slave-side rows of thermo-scatra matrix block to corresponding
      // master-side rows
      thermoscatrablock.add(*Core::LinAlg::MLMultiply(*strategythermo_->p(), true,
                                thermoscatrarowsslave, false, false, false, true),
          false, 1.0, 1.0);
    }
  }

  switch (scatra_field()->matrix_type())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      scatrathermo_domain_interface->complete();
      thermoscatra_domain_interface->complete();

      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      scatrathermo_domain_interface->complete(
          *thermo_field()->discretization()->dof_row_map(), *maps_->Map(0));
      thermoscatra_domain_interface->complete(*maps_->Map(0), *maps_->Map(1));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown matrix type");
      break;
    }
  }
}  // AssembleDomainInterfaceOD

FOUR_C_NAMESPACE_CLOSE
