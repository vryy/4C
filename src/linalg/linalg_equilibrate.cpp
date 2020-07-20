/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for matrix equilibration

\level 1

*/
/*---------------------------------------------------------------------*/
#include "linalg_equilibrate.H"

#include "linalg_utils_sparse_algebra_create.H"

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
LINALG::Equilibration::Equilibration(EquilibrationMethod method, Epetra_Map dofrowmap)
    : invcolsums_(Teuchos::rcp(new Epetra_Vector(dofrowmap, false))),
      invrowsums_(Teuchos::rcp(new Epetra_Vector(dofrowmap, false))),
      method_(method)

{
  return;
};

LINALG::EquilibrationSparse::EquilibrationSparse(EquilibrationMethod method, Epetra_Map dofrowmap)
    : Equilibration(method, dofrowmap)
{
  return;
};

LINALG::EquilibrationBlock::EquilibrationBlock(EquilibrationMethod method, Epetra_Map dofrowmap)
    : Equilibration(method, dofrowmap)
{
  return;
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LINALG::Equilibration::ComputeInvRowSums(const LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invrowsums  //!< inverse sums of absolute values of row entries in matrix
    ) const
{
  // compute inverse row sums of matrix
  if (matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    dserror("Inverse row sums of matrix could not be successfully computed!");

  // take square root of inverse row sums if matrix is scaled from left and right
  if (Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invrowsums->MyLength(); ++i) (*invrowsums)[i] = sqrt((*invrowsums)[i]);
}


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void LINALG::Equilibration::ComputeInvColSums(const LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invcolsums  //!< inverse sums of absolute values of column entries in matrix
    ) const
{
  // compute inverse column sums of matrix
  if (matrix.EpetraMatrix()->InvColSums(*invcolsums))
    dserror("Inverse column sums of matrix could not be successfully computed!");

  // take square root of inverse column sums if matrix is scaled from left and right
  if (Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invcolsums->MyLength(); ++i) (*invcolsums)[i] = sqrt((*invcolsums)[i]);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Equilibration::EquilibrateMatrixRows(LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invrowsums  //!< sums of absolute values of row entries in matrix
    ) const
{
  if (matrix.LeftScale(*invrowsums)) dserror("Row equilibration of matrix failed!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Equilibration::EquilibrateMatrixColumns(LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invcolsums  //!< sums of absolute values of column entries in matrix
    ) const
{
  if (matrix.RightScale(*invcolsums)) dserror("Column equilibration of matrix failed!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Equilibration::UnequilibrateIncrement(
    const Teuchos::RCP<Epetra_Vector>& increment  //!< increment vector
    ) const
{
  // unequilibrate global increment vector if necessary
  if (Method() == EquilibrationMethod::columns_full or
      Method() == EquilibrationMethod::columns_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
  {
    if (increment->Multiply(1., *invcolsums_, *increment, 0.))
      dserror("Unequilibration of global increment vector failed!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Equilibration::EquilibrateRHS(const Teuchos::RCP<Epetra_Vector>& residual) const
{
  // perform equilibration of global residual vector
  if (Method() == EquilibrationMethod::rows_full or
      Method() == EquilibrationMethod::rows_maindiag or
      Method() == EquilibrationMethod::rowsandcolumns_full or
      Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    if (residual->Multiply(1., *invrowsums_, *residual, 0.))
      dserror("Equilibration of global residual vector failed!");
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void LINALG::EquilibrationSparse::EquilibrateSystem(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& residual, const LINALG::MultiMapExtractor& blockmaps) const
{
  if (Method() != EquilibrationMethod::none)
  {
    // check matrix
    Teuchos::RCP<LINALG::SparseMatrix> sparsematrix =
        Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
    if (sparsematrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

    // perform row equilibration
    if (Method() == EquilibrationMethod::rows_full or
        Method() == EquilibrationMethod::rows_maindiag or
        Method() == EquilibrationMethod::rowsandcolumns_full or
        Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      // compute inverse row sums of global system matrix
      ComputeInvRowSums(*sparsematrix, invrowsums_);

      // perform row equilibration of global system matrix
      EquilibrateMatrixRows(*sparsematrix, invrowsums_);
    }

    // perform column equilibration
    if (Method() == EquilibrationMethod::columns_full or
        Method() == EquilibrationMethod::columns_maindiag or
        Method() == EquilibrationMethod::rowsandcolumns_full or
        Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      // compute inverse column sums of global system matrix
      ComputeInvColSums(*sparsematrix, invcolsums_);

      // perform column equilibration of global system matrix
      EquilibrateMatrixColumns(*sparsematrix, invcolsums_);
    }
    // equilibrate residual
    EquilibrateRHS(residual);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void LINALG::EquilibrationBlock::EquilibrateSystem(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& residual, const LINALG::MultiMapExtractor& blockmaps) const
{
  if (Method() != EquilibrationMethod::none)
  {
    // check matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksparsematrix =
        Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
    if (blocksparsematrix == Teuchos::null) dserror("System matrix is not a block sparse matrix!");

    // perform row equilibration
    if (Method() == EquilibrationMethod::rows_full or
        Method() == EquilibrationMethod::rows_maindiag or
        Method() == EquilibrationMethod::rowsandcolumns_full or
        Method() == EquilibrationMethod::rowsandcolumns_maindiag)
    {
      for (int i = 0; i < blocksparsematrix->Rows(); ++i)
      {
        // initialize vector for inverse row sums
        Teuchos::RCP<Epetra_Vector> invrowsums(
            Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i, i).RowMap())));

        // compute inverse row sums of current main diagonal matrix block
        if (Method() == EquilibrationMethod::rows_maindiag or
            Method() == EquilibrationMethod::rowsandcolumns_maindiag)
          ComputeInvRowSums(blocksparsematrix->Matrix(i, i), invrowsums);

        // compute inverse row sums of current row block of global system matrix
        else
        {
          // loop over all column blocks of global system matrix
          for (int j = 0; j < blocksparsematrix->Cols(); ++j)
          {
            // extract current block of global system matrix
            const LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i, j);

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
                if (matrix.EpetraMatrix()->ExtractMyRowCopy(irow, length, numentries, &values[0]))
                  dserror("Cannot extract matrix row with local ID %d from matrix block!", irow);

                // compute and store current row sum
                double rowsum(0.);
                for (int ientry = 0; ientry < numentries; ++ientry)
                  rowsum += std::abs(values[ientry]);
                (*invrowsums)[irow] += rowsum;
              }
            }
          }

          // invert row sums
          if (invrowsums->Reciprocal(*invrowsums)) dserror("Vector could not be inverted!");

          // take square root of inverse row sums if matrix is scaled from left and right
          if (Method() == EquilibrationMethod::rowsandcolumns_full or
              Method() == EquilibrationMethod::rowsandcolumns_maindiag)
            for (int i = 0; i < invrowsums->MyLength(); ++i)
              (*invrowsums)[i] = sqrt((*invrowsums)[i]);
        }

        // perform row equilibration of matrix blocks in current row block of global system
        // matrix
        for (int j = 0; j < blocksparsematrix->Cols(); ++j)
          EquilibrateMatrixRows(blocksparsematrix->Matrix(i, j), invrowsums);

        // insert inverse row sums of current main diagonal matrix block into global vector
        blockmaps.InsertVector(invrowsums, i, invrowsums_);
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
          ComputeInvColSums(blocksparsematrix->Matrix(j, j), invcolsums);

        // compute inverse column sums of current column block of global system matrix
        else
        {
          // loop over all row blocks of global system matrix
          for (int i = 0; i < blocksparsematrix->Rows(); ++i)
          {
            // extract current block of global system matrix
            const LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i, j);

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
                        irow, length, numentries, &values[0], &indices[0]))
                  dserror("Cannot extract matrix row with local ID %d from matrix block!", irow);

                // add entries of current matrix row to column sums
                for (int ientry = 0; ientry < numentries; ++ientry)
                  invcolsums->SumIntoGlobalValue(
                      matrix.ColMap().GID(indices[ientry]), 0, std::abs(values[ientry]));
              }
            }
          }

          // invert column sums
          if (invcolsums->Reciprocal(*invcolsums)) dserror("Vector could not be inverted!");

          // take square root of inverse column sums if matrix is scaled from left and right
          if (Method() == EquilibrationMethod::rowsandcolumns_full or
              Method() == EquilibrationMethod::rowsandcolumns_maindiag)
            for (int i = 0; i < invcolsums->MyLength(); ++i)
              (*invcolsums)[i] = sqrt((*invcolsums)[i]);
        }

        // perform column equilibration of matrix blocks in current column block of global
        // system matrix
        for (int i = 0; i < blocksparsematrix->Rows(); ++i)
          EquilibrateMatrixColumns(blocksparsematrix->Matrix(i, j), invcolsums);

        // insert inverse column sums of current main diagonal matrix block into global vector
        blockmaps.InsertVector(invcolsums, j, invcolsums_);
      }
    }
    // equilibrate residual
    EquilibrateRHS(residual);
  }
}
