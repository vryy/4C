/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic assemble methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_sparse_algebra_assemble.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  assemble a matrix                                         popp 01/08|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();
  // allow Aele to provide entries past the end of lmrow and lmcol that are
  // not used here, therefore check only for ">" rather than "!="
  if (lrowdim != (int)lmrowowner.size() || lrowdim > Aele.M() || lcoldim > Aele.N())
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
  return;
}

/*----------------------------------------------------------------------*
 |  assemble a vector                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  const int ldim = (int)lm.size();
  // allow Vele to provide entries past the end of lm that are not used here,
  // therefore check only for ">" rather than "!="
  if (ldim != (int)lmowner.size() || ldim > Vele.Length()) dserror("Mismatch in dimensions");

  const int myrank = V.Comm().MyPID();

  for (int lrow = 0; lrow < ldim; ++lrow)
  {
    if (lmowner[lrow] != myrank) continue;
    int rgid = lm[lrow];
    if (!V.Map().MyGID(rgid)) dserror("Sparse vector V does not have global row %d", rgid);
    int rlid = V.Map().LID(rgid);
    V[rlid] += Vele[lrow];
  }  // for (int lrow=0; lrow<ldim; ++lrow)

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LINALG::AssembleMyVector(
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
 |  assemble a vector  (wrapper for LINALG::Matrix<3,1>)     katta 10/16|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_Vector& V, LINALG::Matrix<3, 1>& Vele, const std::vector<int>& lm,
    const std::vector<int>& lmowner)
{
  const Epetra_SerialDenseVector VeleNew(::View, &(Vele(0)), 3);
  LINALG::Assemble(V, VeleNew, lm, lmowner);
  return;
}

/*----------------------------------------------------------------------*
 |  assemble a vector  (wrapper for 1 owner)                 katta 10/16|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(
    Epetra_Vector& V, LINALG::Matrix<3, 1>& Vele, const std::vector<int>& lm, const int& lmowner)
{
  const std::vector<int> lmownerNew(3, lmowner);
  LINALG::Assemble(V, Vele, lm, lmownerNew);
  return;
}

/*----------------------------------------------------------------------*
 |  assemble a vector  (wrapper, node-based)                 katta 10/16|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_Vector& V, double& Vele, const int& lm, const int& lmowner)
{
  const Epetra_SerialDenseVector VeleNew(::View, &Vele, 1);
  const std::vector<int> lmNew(1, lm);
  const std::vector<int> lmownerNew(1, lmowner);
  LINALG::Assemble(V, VeleNew, lmNew, lmownerNew);
  return;
}

/*----------------------------------------------------------------------*
 |  assemble a vector into MultiVector (public)              mwgee 01/08|
 *----------------------------------------------------------------------*/
void LINALG::Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  LINALG::Assemble(*(V(n)), Vele, lm, lmowner);
  return;
}

/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
    const Teuchos::RCP<const Epetra_Vector> dbcval,
    const Teuchos::RCP<const Epetra_Vector> dbctoggle)
{
  const Epetra_Vector& dbct = *dbctoggle;
  if (x != Teuchos::null && b != Teuchos::null)
  {
    Epetra_Vector& X = *x;
    Epetra_Vector& B = *b;
    const Epetra_Vector& dbcv = *dbcval;
    // set the prescribed value in x and b
    const int mylength = dbcv.MyLength();
    for (int i = 0; i < mylength; ++i)
      if (dbct[i] == 1.0)
      {
        X[i] = dbcv[i];
        B[i] = dbcv[i];
      }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
    const Teuchos::RCP<const Epetra_Vector> dbcval, const Epetra_Map& dbcmap)
{
  if (not dbcmap.UniqueGIDs()) dserror("unique map required");

  if (x != Teuchos::null and b != Teuchos::null)
  {
    Epetra_Vector& X = *x;
    Epetra_Vector& B = *b;
    const Epetra_Vector& dbcv = *dbcval;

    // We use two maps since we want to allow dbcv and X to be independent of
    // each other. So we are slow and flexible...
    const Epetra_BlockMap& xmap = X.Map();
    const Epetra_BlockMap& dbcvmap = dbcv.Map();

    const int mylength = dbcmap.NumMyElements();
    const int* mygids = dbcmap.MyGlobalElements();
    for (int i = 0; i < mylength; ++i)
    {
      int gid = mygids[i];

      int dbcvlid = dbcvmap.LID(gid);
      if (dbcvlid < 0) dserror("illegal Dirichlet map");

      int xlid = xmap.LID(gid);
      if (xlid < 0) dserror("illegal Dirichlet map");

      X[xlid] = dbcv[dbcvlid];
      B[xlid] = dbcv[dbcvlid];
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& b,
    const Teuchos::RCP<const Epetra_Vector> dbcval, const Epetra_Map& dbcmap)
{
  if (not dbcmap.UniqueGIDs()) dserror("unique map required");

  if (b != Teuchos::null)
  {
    Epetra_Vector& B = *b;
    const Epetra_Vector& dbcv = *dbcval;

    // We use two maps since we want to allow dbcv and X to be independent of
    // each other. So we are slow and flexible...
    const Epetra_BlockMap& bmap = B.Map();
    const Epetra_BlockMap& dbcvmap = dbcv.Map();

    const int mylength = dbcmap.NumMyElements();
    const int* mygids = dbcmap.MyGlobalElements();
    for (int i = 0; i < mylength; ++i)
    {
      int gid = mygids[i];

      int dbcvlid = dbcvmap.LID(gid);

      int blid = bmap.LID(gid);
      // Note:
      // if gid is not found in vector b, just continue
      // b might only be a subset of a larger field vector
      if (blid >= 0)
      {
        if (dbcvlid < 0)
          dserror("illegal Dirichlet map");
        else
          B[blid] = dbcv[dbcvlid];
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
    Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
    const Teuchos::RCP<const Epetra_Vector> dbcval,
    const Teuchos::RCP<const Epetra_Vector> dbctoggle)
{
  A->ApplyDirichlet(dbctoggle);
  ApplyDirichlettoSystem(x, b, dbcval, dbctoggle);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
    Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
    const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap)
{
  A->ApplyDirichlet(dbcmap);
  ApplyDirichlettoSystem(x, b, dbcval, dbcmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
    Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
    Teuchos::RCP<const LINALG::SparseMatrix> trafo, const Teuchos::RCP<const Epetra_Vector>& dbcval,
    const Epetra_Map& dbcmap)
{
  if (trafo != Teuchos::null)
    Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(A, true)->ApplyDirichletWithTrafo(
        trafo, dbcmap);
  else
    // trafo==Teuchos::null
    A->ApplyDirichlet(dbcmap);
  ApplyDirichlettoSystem(x, b, dbcval, dbcmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
    Teuchos::RCP<Epetra_Vector>& b, Teuchos::RCP<const LINALG::SparseMatrix> trafo,
    const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap)
{
  if (trafo != Teuchos::null)
    Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(A, true)->ApplyDirichletWithTrafo(
        trafo, dbcmap);
  else
    // trafo==Teuchos::null
    A->ApplyDirichlet(dbcmap);
  ApplyDirichlettoSystem(b, dbcval, dbcmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::MapExtractor> LINALG::ConvertDirichletToggleVectorToMaps(
    const Teuchos::RCP<const Epetra_Vector>& dbctoggle)
{
  const Epetra_BlockMap& fullblockmap = dbctoggle->Map();
  // this copy is needed because the constructor of LINALG::MapExtractor
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
    int* myglobalelements = NULL;
    if (dbcgids.size() > 0)
    {
      nummyelements = dbcgids.size();
      myglobalelements = &(dbcgids[0]);
    }
    dbcmap = Teuchos::rcp(
        new Epetra_Map(-1, nummyelements, myglobalelements, fullmap.IndexBase(), fullmap.Comm()));
  }
  // build map of free DOFs
  Teuchos::RCP<Epetra_Map> freemap = Teuchos::null;
  {
    int nummyelements = 0;
    int* myglobalelements = NULL;
    if (freegids.size() > 0)
    {
      nummyelements = freegids.size();
      myglobalelements = &(freegids[0]);
    }
    freemap = Teuchos::rcp(
        new Epetra_Map(-1, nummyelements, myglobalelements, fullmap.IndexBase(), fullmap.Comm()));
  }

  // build and return the map extractor of Dirichlet-conditioned and free DOFs
  return Teuchos::rcp(new LINALG::MapExtractor(fullmap, dbcmap, freemap));
}
