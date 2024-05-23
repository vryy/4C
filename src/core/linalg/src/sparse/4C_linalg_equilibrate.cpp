/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for matrix equilibration

\level 1

*/
/*---------------------------------------------------------------------*/
#include "4C_linalg_equilibrate.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
CORE::LINALG::Equilibration::Equilibration(Teuchos::RCP<const Epetra_Map> dofrowmap)
    : invcolsums_(CORE::LINALG::CreateVector(*dofrowmap, false)),
      invrowsums_(CORE::LINALG::CreateVector(*dofrowmap, false))

{
}

CORE::LINALG::EquilibrationUniversal::EquilibrationUniversal(
    EquilibrationMethod method, Teuchos::RCP<const Epetra_Map> dofrowmap)
    : Equilibration(dofrowmap), method_(method)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::EquilibrationSparse::EquilibrationSparse(
    EquilibrationMethod method, Teuchos::RCP<const Epetra_Map> dofrowmap)
    : EquilibrationUniversal(method, dofrowmap)
{
  if (method == EquilibrationMethod::symmetry)
    FOUR_C_THROW("symmetric equilibration not implemented for sparse matrices");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::EquilibrationBlock::EquilibrationBlock(
    EquilibrationMethod method, Teuchos::RCP<const Epetra_Map> dofrowmap)
    : EquilibrationUniversal(method, dofrowmap)
{
  if (method == EquilibrationMethod::symmetry)
    FOUR_C_THROW("symmetric equilibration not implemented for block matrices");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::EquilibrationBlockSpecific::EquilibrationBlockSpecific(
    const std::vector<EquilibrationMethod>& method, Teuchos::RCP<const Epetra_Map> dofrowmap)
    : Equilibration(dofrowmap), method_blocks_(method)
{
  for (const auto& method_block : method_blocks_)
  {
    if (method_block == EquilibrationMethod::columns_full or
        method_block == EquilibrationMethod::rows_full or
        method_block == EquilibrationMethod::rowsandcolumns_full)
      FOUR_C_THROW("full matrix equilibration not reasonable for block based equilibration");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::LINALG::Equilibration::ComputeInvRowSums(const CORE::LINALG::SparseMatrix& matrix,
    Teuchos::RCP<Epetra_Vector> invrowsums, const EquilibrationMethod method) const
{
  // compute inverse row sums of matrix
  if (matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    FOUR_C_THROW("Inverse row sums of matrix could not be successfully computed!");

  // take square root of inverse row sums if matrix is scaled from left and right
  if (method == EquilibrationMethod::rowsandcolumns_full or
      method == EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invrowsums->MyLength(); ++i) (*invrowsums)[i] = std::sqrt((*invrowsums)[i]);
}


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void CORE::LINALG::Equilibration::ComputeInvColSums(const CORE::LINALG::SparseMatrix& matrix,
    Teuchos::RCP<Epetra_Vector> invcolsums, const EquilibrationMethod method) const
{
  // compute inverse column sums of matrix
  if (matrix.EpetraMatrix()->InvColSums(*invcolsums))
    FOUR_C_THROW("Inverse column sums of matrix could not be successfully computed!");

  // take square root of inverse column sums if matrix is scaled from left and right
  if (method == EquilibrationMethod::rowsandcolumns_full or
      method == EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invcolsums->MyLength(); ++i) (*invcolsums)[i] = std::sqrt((*invcolsums)[i]);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void CORE::LINALG::Equilibration::ComputeInvSymmetry(
    const CORE::LINALG::SparseMatrix& matrix, Teuchos::RCP<Epetra_Vector> invsymmetry) const
{
  Teuchos::RCP<Epetra_Vector> diag = CORE::LINALG::CreateVector(matrix.RangeMap(), true);
  matrix.ExtractDiagonalCopy(*diag);

  for (int my_row = 0; my_row < diag->Map().NumMyElements(); ++my_row)
  {
    (*invsymmetry)[my_row] = 1.0 / std::sqrt((*diag)[my_row]);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Equilibration::equilibrate_matrix_rows(
    CORE::LINALG::SparseMatrix& matrix, Teuchos::RCP<const Epetra_Vector> invrowsums) const
{
  if (matrix.LeftScale(*invrowsums)) FOUR_C_THROW("Row equilibration of matrix failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Equilibration::equilibrate_matrix_columns(
    CORE::LINALG::SparseMatrix& matrix, Teuchos::RCP<const Epetra_Vector> invcolsums) const
{
  if (matrix.RightScale(*invcolsums)) FOUR_C_THROW("Column equilibration of matrix failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Equilibration::EquilibrateSystem(
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> residual,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const
{
  // Equilibrate the matrix given the chosen method and matrix type
  EquilibrateMatrix(systemmatrix, blockmaps);

  // Equilibrate the RHS
  EquilibrateRHS(residual);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationUniversal::unequilibrate_increment(
    Teuchos::RCP<Epetra_Vector> increment) const
{
  // unequilibrate global increment vector if necessary
  if (Method() == EquilibrationMethod::columns_full or
      Method() == EquilibrationMethod::columns_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    if (increment->Multiply(1.0, *invcolsums_, *increment, 0.0))
      FOUR_C_THROW("Unequilibration of global increment vector failed!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationUniversal::EquilibrateRHS(
    Teuchos::RCP<Epetra_Vector> residual) const
{
  // perform equilibration of global residual vector
  if (Method() == EquilibrationMethod::rows_full or
      Method() == EquilibrationMethod::rows_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    if (residual->Multiply(1.0, *invrowsums_, *residual, 0.0))
      FOUR_C_THROW("Equilibration of global residual vector failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationBlockSpecific::unequilibrate_increment(
    Teuchos::RCP<Epetra_Vector> increment) const
{
  if (increment->Multiply(1.0, *invcolsums_, *increment, 0.0))
    FOUR_C_THROW("Unequilibration of global increment vector failed!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationBlockSpecific::EquilibrateRHS(
    Teuchos::RCP<Epetra_Vector> residual) const
{
  if (residual->Multiply(1.0, *invrowsums_, *residual, 0.0))
    FOUR_C_THROW("Equilibration of global residual vector failed!");
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationSparse::EquilibrateMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const
{
  EquilibrateMatrix(systemmatrix);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationSparse::EquilibrateMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sparsematrix =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(systemmatrix);

  // perform row equilibration
  if (Method() == EquilibrationMethod::rows_full or
      Method() == EquilibrationMethod::rows_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    // compute inverse row sums of global system matrix
    ComputeInvRowSums(*sparsematrix, invrowsums_, Method());

    // perform row equilibration of global system matrix
    equilibrate_matrix_rows(*sparsematrix, invrowsums_);
  }

  // perform column equilibration
  if (Method() == EquilibrationMethod::columns_full or
      Method() == EquilibrationMethod::columns_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    // compute inverse column sums of global system matrix
    ComputeInvColSums(*sparsematrix, invcolsums_, Method());

    // perform column equilibration of global system matrix
    equilibrate_matrix_columns(*sparsematrix, invcolsums_);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationBlock::EquilibrateMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const
{
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blocksparsematrix =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  // perform row equilibration
  if (Method() == EquilibrationMethod::rows_full or
      Method() == EquilibrationMethod::rows_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    for (int i = 0; i < blocksparsematrix->Rows(); ++i)
    {
      // initialize vector for inverse row sums
      auto invrowsums = Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i, i).RowMap()));

      // compute inverse row sums of current main diagonal matrix block
      if (Method() == EquilibrationMethod::rows_maindiag or
          Method() == EquilibrationMethod::rowsandcolumns_maindiag)
      {
        ComputeInvRowSums(blocksparsematrix->Matrix(i, i), invrowsums, Method());
      }
      // compute inverse row sums of current row block of global system matrix
      else
      {
        // loop over all column blocks of global system matrix
        for (int j = 0; j < blocksparsematrix->Cols(); ++j)
        {
          // extract current block of global system matrix
          const CORE::LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i, j);

          // loop over all rows of current matrix block
          for (int irow = 0; irow < matrix.RowMap().NumMyElements(); ++irow)
          {
            // determine length of current matrix row
            const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

            if (length > 0)
            {
              // extract current matrix row from matrix block
              int numentries(0);
              std::vector<double> values(length, 0.);
              if (matrix.EpetraMatrix()->ExtractMyRowCopy(irow, length, numentries, values.data()))
                FOUR_C_THROW("Cannot extract matrix row with local ID %d from matrix block!", irow);

              // compute and store current row sum
              double rowsum(0.);
              for (int ientry = 0; ientry < numentries; ++ientry)
                rowsum += std::abs(values[ientry]);
              (*invrowsums)[irow] += rowsum;
            }
          }
        }

        // invert row sums
        if (invrowsums->Reciprocal(*invrowsums)) FOUR_C_THROW("Vector could not be inverted!");

        // take square root of inverse row sums if matrix is scaled from left and right
        if (Method() == EquilibrationMethod::rowsandcolumns_full or
            Method() == EquilibrationMethod::rowsandcolumns_maindiag)
          for (int j = 0; j < invrowsums->MyLength(); ++j)
            (*invrowsums)[j] = std::sqrt((*invrowsums)[j]);
      }

      // perform row equilibration of matrix blocks in current row block of global system
      // matrix
      for (int j = 0; j < blocksparsematrix->Cols(); ++j)
        equilibrate_matrix_rows(blocksparsematrix->Matrix(i, j), invrowsums);

      // insert inverse row sums of current main diagonal matrix block into global vector
      blockmaps->InsertVector(invrowsums, i, invrowsums_);
    }
  }

  // perform column equilibration
  if (Method() == EquilibrationMethod::columns_full or
      Method() == EquilibrationMethod::columns_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    for (int j = 0; j < blocksparsematrix->Cols(); ++j)
    {
      // initialize vector for inverse column sums
      Teuchos::RCP<Epetra_Vector> invcolsums(
          Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(j, j).DomainMap())));

      // compute inverse column sums of current main diagonal matrix block
      if (Method() == EquilibrationMethod::columns_maindiag or
          Method() == EquilibrationMethod::rowsandcolumns_maindiag)
      {
        ComputeInvColSums(blocksparsematrix->Matrix(j, j), invcolsums, Method());
      }
      // compute inverse column sums of current column block of global system matrix
      else
      {
        // loop over all row blocks of global system matrix
        for (int i = 0; i < blocksparsematrix->Rows(); ++i)
        {
          // extract current block of global system matrix
          const CORE::LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i, j);

          // loop over all rows of current matrix block
          for (int irow = 0; irow < matrix.RowMap().NumMyElements(); ++irow)
          {
            // determine length of current matrix row
            const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

            if (length > 0)
            {
              // extract current matrix row from matrix block
              int numentries(0);
              std::vector<double> values(length, 0.);
              std::vector<int> indices(length, 0);
              if (matrix.EpetraMatrix()->ExtractMyRowCopy(
                      irow, length, numentries, values.data(), indices.data()))
                FOUR_C_THROW("Cannot extract matrix row with local ID %d from matrix block!", irow);

              // add entries of current matrix row to column sums
              for (int ientry = 0; ientry < numentries; ++ientry)
                invcolsums->SumIntoGlobalValue(
                    matrix.ColMap().GID(indices[ientry]), 0, std::abs(values[ientry]));
            }
          }
        }

        // invert column sums
        if (invcolsums->Reciprocal(*invcolsums)) FOUR_C_THROW("Vector could not be inverted!");

        // take square root of inverse column sums if matrix is scaled from left and right
        if (Method() == EquilibrationMethod::rowsandcolumns_full or
            Method() == EquilibrationMethod::rowsandcolumns_maindiag)
          for (int i = 0; i < invcolsums->MyLength(); ++i)
            (*invcolsums)[i] = std::sqrt((*invcolsums)[i]);
      }

      // perform column equilibration of matrix blocks in current column block of global
      // system matrix
      for (int i = 0; i < blocksparsematrix->Rows(); ++i)
        equilibrate_matrix_columns(blocksparsematrix->Matrix(i, j), invcolsums);

      // insert inverse column sums of current main diagonal matrix block into global vector
      blockmaps->InsertVector(invcolsums, j, invcolsums_);
    }
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void CORE::LINALG::EquilibrationBlockSpecific::EquilibrateMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const
{
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blocksparsematrix =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(systemmatrix);

  if (blocksparsematrix->Rows() != static_cast<int>(method_blocks_.size()))
    FOUR_C_THROW("No match between number of equilibration methods and Matrix blocks");

  // init: no scaling
  invcolsums_->PutScalar(1.0);
  invrowsums_->PutScalar(1.0);

  // loop over all blocks of matrix and apply equilibration for each block
  for (int i = 0; i < blocksparsematrix->Rows(); ++i)
  {
    const EquilibrationMethod method = method_blocks_.at(i);
    if (method == EquilibrationMethod::rows_maindiag or
        method == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      auto invrowsums = Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i, i).RowMap()));
      ComputeInvRowSums(blocksparsematrix->Matrix(i, i), invrowsums, method);

      // perform row equilibration of matrix blocks in current row block of global system
      // matrix
      for (int j = 0; j < blocksparsematrix->Cols(); ++j)
        equilibrate_matrix_rows(blocksparsematrix->Matrix(i, j), invrowsums);

      // insert inverse row sums of current main diagonal matrix block into global vector
      blockmaps->InsertVector(invrowsums, i, invrowsums_);
    }
    if (method == EquilibrationMethod::columns_maindiag or
        method == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      auto invcolsums =
          Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i, i).DomainMap()));
      ComputeInvColSums(blocksparsematrix->Matrix(i, i), invcolsums, method);

      // perform column equilibration of matrix blocks in current column block of global
      // system matrix
      for (int j = 0; j < blocksparsematrix->Cols(); ++j)
        equilibrate_matrix_columns(blocksparsematrix->Matrix(j, i), invcolsums);

      // insert inverse column sums of current main diagonal matrix block into global vector
      blockmaps->InsertVector(invcolsums, i, invcolsums_);
    }
    if (method == EquilibrationMethod::symmetry)
    {
      auto invsymmetry =
          CORE::LINALG::CreateVector(blocksparsematrix->Matrix(i, i).RangeMap(), true);

      ComputeInvSymmetry(blocksparsematrix->Matrix(i, i), invsymmetry);

      for (int j = 0; j < blocksparsematrix->Cols(); ++j)
      {
        equilibrate_matrix_rows(blocksparsematrix->Matrix(i, j), invsymmetry);
        equilibrate_matrix_columns(blocksparsematrix->Matrix(j, i), invsymmetry);
      }

      // insert inverse row sums of current main diagonal matrix block into global vector
      blockmaps->InsertVector(invsymmetry, i, invcolsums_);
      blockmaps->InsertVector(invsymmetry, i, invrowsums_);
    }
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::Equilibration> CORE::LINALG::BuildEquilibration(MatrixType type,
    const std::vector<EquilibrationMethod>& method, Teuchos::RCP<const Epetra_Map> dofrowmap)
{
  Teuchos::RCP<CORE::LINALG::Equilibration> equilibration = Teuchos::null;

  if (method.size() == 1)
  {
    EquilibrationMethod method_global = method.at(0);

    if (method_global == CORE::LINALG::EquilibrationMethod::none)
      equilibration = Teuchos::rcp(new CORE::LINALG::EquilibrationNone(dofrowmap));
    else
    {
      switch (type)
      {
        case CORE::LINALG::MatrixType::sparse:
        {
          equilibration =
              Teuchos::rcp(new CORE::LINALG::EquilibrationSparse(method_global, dofrowmap));
          break;
        }
        case CORE::LINALG::MatrixType::block_field:
        case CORE::LINALG::MatrixType::block_condition:
        case CORE::LINALG::MatrixType::block_condition_dof:
        {
          equilibration =
              Teuchos::rcp(new CORE::LINALG::EquilibrationBlock(method_global, dofrowmap));
          break;
        }
        default:
        {
          FOUR_C_THROW("Unknown type of global system matrix");
          break;
        }
      }
    }
  }
  else
    equilibration = Teuchos::rcp(new CORE::LINALG::EquilibrationBlockSpecific(method, dofrowmap));

  return equilibration;
}

FOUR_C_NAMESPACE_CLOSE
