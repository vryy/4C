/*----------------------------------------------------------------------*/
/*!
\file drt_assemblestrategy.cpp

\brief Routines for the handing a collection of element matrices and vectors to
       the actual assembly calls into one global sparse matrix and global load vector

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#include "drt_assemblestrategy.H"
#include "drt_discret.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

DRT::AssembleStrategy::AssembleStrategy(int firstdofset, int seconddofset,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
    : firstdofset_(firstdofset),
      seconddofset_(seconddofset),
      systemmatrix1_(systemmatrix1),
      systemmatrix2_(systemmatrix2),
      systemvector1_(systemvector1),
      systemvector2_(systemvector2),
      systemvector3_(systemvector3)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::AssembleStrategy::~AssembleStrategy() {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> DRT::AssembleStrategy::MatrixGraph(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Map> dbcmap)
{
  if (!dis.Filled()) dserror("FillComplete() was not called on this discretization");

  const Epetra_Map& dofrowmap = *dis.DofRowMap();
  const int myrank = dis.Comm().MyPID();

  // build graph. insert directly into CrsGraph. there used to be a temporary
  // field of type std::map<int,std::set<int> > to set up the dofs, but std::set
  // is usually much slower for inserting a few ints than the internal routines
  // in CrsGraph.
  Teuchos::RCP<Epetra_CrsGraph> crsgraph(new Epetra_CrsGraph(Copy, dofrowmap, 108, false));

  Element::LocationArray la(dis.NumDofSets());
  int row = FirstDofSet();
  int col = SecondDofSet();

  const int numcolele = dis.NumMyColElements();
  for (int i = 0; i < numcolele; ++i)
  {
    DRT::Element* actele = dis.lColElement(i);
    actele->LocationVector(dis, la, false);

    const std::vector<int>& lmrow = la[row].lm_;
    const std::vector<int>& lmrowowner = la[row].lmowner_;
    const std::vector<int>& lmcol = la[col].lm_;

    const int lrowdim = (int)lmrow.size();
    const int lcoldim = (int)lmcol.size();

    for (int lrow = 0; lrow < lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lmrow[lrow];
      //#ifdef DEBUG
      if (!dofrowmap.MyGID(rgid)) dserror("Proc %d does not have global row %d", myrank, rgid);
      //#endif

      // if we have a Dirichlet map check if this row is a Dirichlet row
      int err;
      if (dbcmap != Teuchos::null and dbcmap->MyGID(rgid))
        err = crsgraph->InsertGlobalIndices(rgid, 1, &rgid);
      else
        err = crsgraph->InsertGlobalIndices(rgid, lcoldim, const_cast<int*>(&lmcol[0]));
      if (err < 0) dserror("graph->InsertGlobalIndices returned err=%d", err);
    }
  }

  int err = crsgraph->FillComplete();
  if (err) dserror("graph->FillComplete() returned err=%d", err);
  err = crsgraph->OptimizeStorage();
  if (err) dserror("graph->OptimizeStorage() returned err=%d", err);

  return crsgraph;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Zero()
{
  if (Assemblemat1())
  {
    systemmatrix1_->Zero();
  }
  if (Assemblemat2())
  {
    systemmatrix1_->Zero();
  }
  if (Assemblevec1())
  {
    systemvector1_->PutScalar(0.0);
  }
  if (Assemblevec2())
  {
    systemvector2_->PutScalar(0.0);
  }
  if (Assemblevec3())
  {
    systemvector3_->PutScalar(0.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Complete()
{
  if (Assemblemat1())
  {
    systemmatrix1_->Complete();
  }
  if (Assemblemat2())
  {
    systemmatrix2_->Complete();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::ClearElementStorage(int rdim, int cdim)
{
  if (Assemblemat1())
  {
    if (elematrix1_.M() != rdim or elematrix1_.N() != cdim)
      elematrix1_.Shape(rdim, cdim);
    else
      memset(elematrix1_.A(), 0, rdim * cdim * sizeof(double));
  }
  if (Assemblemat2())
  {
    if (elematrix2_.M() != rdim or elematrix2_.N() != cdim)
      elematrix2_.Shape(rdim, cdim);
    else
      memset(elematrix2_.A(), 0, rdim * cdim * sizeof(double));
  }
  if (Assemblevec1())
  {
    if (elevector1_.Length() != rdim)
      elevector1_.Size(rdim);
    else
      memset(elevector1_.Values(), 0, rdim * sizeof(double));
  }
  if (Assemblevec2())
  {
    if (elevector2_.Length() != rdim)
      elevector2_.Size(rdim);
    else
      memset(elevector2_.Values(), 0, rdim * sizeof(double));
  }
  if (Assemblevec3())
  {
    if (elevector3_.Length() != rdim)
      elevector3_.Size(rdim);
    else
      memset(elevector3_.Values(), 0, rdim * sizeof(double));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat, int eid,
    const std::vector<int>& lmstride, const Epetra_SerialDenseMatrix& Aele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  sysmat.Assemble(eid, lmstride, Aele, lm, lmowner);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat, int eid,
    const std::vector<int>& lmstride, const Epetra_SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  sysmat.Assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid)
{
  sysmat.Assemble(val, rgid, cgid);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  LINALG::Assemble(V, Vele, lm, lmowner);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(Epetra_MultiVector& V, const int n,
    const Epetra_SerialDenseVector& Vele, const std::vector<int>& lm,
    const std::vector<int>& lmowner)
{
  LINALG::Assemble(V, n, Vele, lm, lmowner);
}
