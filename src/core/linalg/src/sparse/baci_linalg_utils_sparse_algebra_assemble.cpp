/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic assemble methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_linalg_utils_sparse_algebra_assemble.hpp"

#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  assemble a matrix                                         popp 01/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::Assemble(Epetra_CrsMatrix& A, const CORE::LINALG::SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();
  // allow Aele to provide entries past the end of lmrow and lmcol that are
  // not used here, therefore check only for ">" rather than "!="
  if (lrowdim != (int)lmrowowner.size() || lrowdim > Aele.numRows() || lcoldim > Aele.numCols())
    dserror("Mismatch in dimensions");

  const int myrank = A.Comm().MyPID();
  const Epetra_Map& rowmap = A.RowMap();

  // this 'Assemble' is not implemented for a Filled() matrix A
  if (A.Filled())
    dserror("Sparse matrix A is already Filled()");

  else
  {
    // loop rows of local matrix
    for (int lrow = 0; lrow < lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lmrow[lrow];
      if (!(rowmap.MyGID(rgid))) dserror("Sparse matrix A does not have global row %d", rgid);

      for (int lcol = 0; lcol < lcoldim; ++lcol)
      {
        double val = Aele(lrow, lcol);
        int cgid = lmcol[lcol];

        // Now that we do not rebuild the sparse mask in each step, we
        // are bound to assemble the whole thing. Zeros included.
        int errone = A.SumIntoGlobalValues(rgid, 1, &val, &cgid);
        if (errone > 0)
        {
          int errtwo = A.InsertGlobalValues(rgid, 1, &val, &cgid);
          if (errtwo < 0)
            dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code %d", errtwo);
        }
        else if (errone)
          dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d", errone);
      }  // for (int lcol=0; lcol<lcoldim; ++lcol)
    }    // for (int lrow=0; lrow<lrowdim; ++lrow)
  }
}

/*----------------------------------------------------------------------*
 |  assemble a vector                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void CORE::LINALG::Assemble(Epetra_Vector& V, const CORE::LINALG::SerialDenseVector& Vele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  const int ldim = (int)lm.size();
  // allow Vele to provide entries past the end of lm that are not used here,
  // therefore check only for ">" rather than "!="
  if (ldim != (int)lmowner.size() || ldim > Vele.length()) dserror("Mismatch in dimensions");

  const int myrank = V.Comm().MyPID();

  for (int lrow = 0; lrow < ldim; ++lrow)
  {
    if (lmowner[lrow] != myrank) continue;
    int rgid = lm[lrow];
    if (!V.Map().MyGID(rgid)) dserror("Sparse vector V does not have global row %d", rgid);
    int rlid = V.Map().LID(rgid);
    V[rlid] += Vele[lrow];
  }  // for (int lrow=0; lrow<ldim; ++lrow)
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::LINALG::AssembleMyVector(
    double scalar_target, Epetra_Vector& target, double scalar_source, const Epetra_Vector& source)
{
  for (int slid = 0; slid < source.Map().NumMyElements(); ++slid)
  {
    const int sgid = source.Map().GID(slid);
    const int tlid = target.Map().LID(sgid);
    if (tlid == -1)
      dserror(
          "The target vector has no global row %i"
          " on processor %i!",
          sgid, target.Comm().MyPID());

    // update the vector row
    target[tlid] = scalar_target * target[tlid] + scalar_source * source[slid];
  }
}

/*----------------------------------------------------------------------*
 |  assemble a vector  (wrapper for CORE::LINALG::Matrix<3,1>)     katta 10/16|
 *----------------------------------------------------------------------*/
void CORE::LINALG::Assemble(Epetra_Vector& V, CORE::LINALG::Matrix<3, 1>& Vele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  const CORE::LINALG::SerialDenseVector VeleNew(Teuchos::View, &(Vele(0)), 3);
  CORE::LINALG::Assemble(V, VeleNew, lm, lmowner);
}

/*----------------------------------------------------------------------*
 |  assemble a vector  (wrapper for 1 owner)                 katta 10/16|
 *----------------------------------------------------------------------*/
void CORE::LINALG::Assemble(Epetra_Vector& V, CORE::LINALG::Matrix<3, 1>& Vele,
    const std::vector<int>& lm, const int& lmowner)
{
  const std::vector<int> lmownerNew(3, lmowner);
  CORE::LINALG::Assemble(V, Vele, lm, lmownerNew);
}

/*----------------------------------------------------------------------*
 |  assemble a vector  (wrapper, node-based)                 katta 10/16|
 *----------------------------------------------------------------------*/
void CORE::LINALG::Assemble(Epetra_Vector& V, double& Vele, const int& lm, const int& lmowner)
{
  const CORE::LINALG::SerialDenseVector VeleNew(Teuchos::View, &Vele, 1);
  const std::vector<int> lmNew(1, lm);
  const std::vector<int> lmownerNew(1, lmowner);
  CORE::LINALG::Assemble(V, VeleNew, lmNew, lmownerNew);
}

/*----------------------------------------------------------------------*
 |  assemble a vector into MultiVector (public)              mwgee 01/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::Assemble(Epetra_MultiVector& V, const int n,
    const CORE::LINALG::SerialDenseVector& Vele, const std::vector<int>& lm,
    const std::vector<int>& lmowner)
{
  CORE::LINALG::Assemble(*(V(n)), Vele, lm, lmowner);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::ApplyDirichletToSystem(
    Epetra_Vector& x, Epetra_Vector& b, const Epetra_Vector& dbcval, const Epetra_Vector& dbctoggle)
{
  // set the prescribed value in x and b
  const int mylength = dbcval.MyLength();
  for (int i = 0; i < mylength; ++i)
  {
    if (dbctoggle[i] == 1.0)
    {
      x[i] = dbcval[i];
      b[i] = dbcval[i];
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::ApplyDirichletToSystem(
    Epetra_Vector& x, Epetra_Vector& b, const Epetra_Vector& dbcval, const Epetra_Map& dbcmap)
{
  if (not dbcmap.UniqueGIDs()) dserror("unique map required");

  // We use two maps since we want to allow dbcv and X to be independent of
  // each other. So we are slow and flexible...
  const Epetra_BlockMap& xmap = x.Map();
  const Epetra_BlockMap& dbcvmap = dbcval.Map();

  const int mylength = dbcmap.NumMyElements();
  const int* mygids = dbcmap.MyGlobalElements();
  for (int i = 0; i < mylength; ++i)
  {
    int gid = mygids[i];

    int dbcvlid = dbcvmap.LID(gid);
    if (dbcvlid < 0) dserror("illegal Dirichlet map");

    int xlid = xmap.LID(gid);
    if (xlid < 0) dserror("illegal Dirichlet map");

    x[xlid] = dbcval[dbcvlid];
    b[xlid] = dbcval[dbcvlid];
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::ApplyDirichletToSystem(
    Epetra_Vector& b, const Epetra_Vector& dbcval, const Epetra_Map& dbcmap)
{
  if (not dbcmap.UniqueGIDs()) dserror("unique map required");

  const int mylength = dbcmap.NumMyElements();
  const int* mygids = dbcmap.MyGlobalElements();
  for (int i = 0; i < mylength; ++i)
  {
    const int gid = mygids[i];

    const int dbcvlid = dbcval.Map().LID(gid);

    const int blid = b.Map().LID(gid);
    // Note:
    // if gid is not found in vector b, just continue
    // b might only be a subset of a larger field vector
    if (blid >= 0)
    {
      if (dbcvlid < 0)
        dserror("illegal Dirichlet map");
      else
        b[blid] = dbcval[dbcvlid];
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::ApplyDirichletToSystem(CORE::LINALG::SparseOperator& A, Epetra_Vector& x,
    Epetra_Vector& b, const Epetra_Vector& dbcval, const Epetra_Vector& dbctoggle)
{
  A.ApplyDirichlet(dbctoggle);
  ApplyDirichletToSystem(x, b, dbcval, dbctoggle);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::ApplyDirichletToSystem(CORE::LINALG::SparseOperator& A, Epetra_Vector& x,
    Epetra_Vector& b, const Epetra_Vector& dbcval, const Epetra_Map& dbcmap)
{
  A.ApplyDirichlet(dbcmap);
  ApplyDirichletToSystem(x, b, dbcval, dbcmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::ApplyDirichletToSystem(CORE::LINALG::SparseMatrix& A, Epetra_Vector& x,
    Epetra_Vector& b, const CORE::LINALG::SparseMatrix& trafo, const Epetra_Vector& dbcval,
    const Epetra_Map& dbcmap)
{
  A.ApplyDirichletWithTrafo(trafo, dbcmap);
  ApplyDirichletToSystem(x, b, dbcval, dbcmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::MapExtractor> CORE::LINALG::ConvertDirichletToggleVectorToMaps(
    const Teuchos::RCP<const Epetra_Vector>& dbctoggle)
{
  const Epetra_BlockMap& fullblockmap = dbctoggle->Map();
  // this copy is needed because the constructor of CORE::LINALG::MapExtractor
  // accepts only Epetra_Map and not Epetra_BlockMap
  const Epetra_Map fullmap =
      Epetra_Map(fullblockmap.NumGlobalElements(), fullblockmap.NumMyElements(),
          fullblockmap.MyGlobalElements(), fullblockmap.IndexBase(), fullblockmap.Comm());
  const int mylength = dbctoggle->MyLength();
  const int* fullgids = fullmap.MyGlobalElements();
  // build sets containing the DBC or free global IDs, respectively
  std::vector<int> dbcgids;
  std::vector<int> freegids;
  for (int i = 0; i < mylength; ++i)
  {
    const int gid = fullgids[i];
    const int compo = (int)round((*dbctoggle)[i]);
    if (compo == 0)
      freegids.push_back(gid);
    else if (compo == 1)
      dbcgids.push_back(gid);
    else
      dserror("Unexpected component %f. It is neither 1.0 nor 0.0.", (*dbctoggle)[i]);
  }
  // build map of Dirichlet DOFs
  Teuchos::RCP<Epetra_Map> dbcmap = Teuchos::null;
  {
    int nummyelements = 0;
    int* myglobalelements = nullptr;
    if (dbcgids.size() > 0)
    {
      nummyelements = dbcgids.size();
      myglobalelements = dbcgids.data();
    }
    dbcmap = Teuchos::rcp(
        new Epetra_Map(-1, nummyelements, myglobalelements, fullmap.IndexBase(), fullmap.Comm()));
  }
  // build map of free DOFs
  Teuchos::RCP<Epetra_Map> freemap = Teuchos::null;
  {
    int nummyelements = 0;
    int* myglobalelements = nullptr;
    if (freegids.size() > 0)
    {
      nummyelements = freegids.size();
      myglobalelements = freegids.data();
    }
    freemap = Teuchos::rcp(
        new Epetra_Map(-1, nummyelements, myglobalelements, fullmap.IndexBase(), fullmap.Comm()));
  }

  // build and return the map extractor of Dirichlet-conditioned and free DOFs
  return Teuchos::rcp(new CORE::LINALG::MapExtractor(fullmap, dbcmap, freemap));
}

BACI_NAMESPACE_CLOSE
