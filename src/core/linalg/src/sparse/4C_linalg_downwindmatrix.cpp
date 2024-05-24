/*----------------------------------------------------------------------*/
/*! \file

\brief An approximate block factorization preconditioner based on the
       SIMPLE family of methods

\level 2

*----------------------------------------------------------------------*/
#include "4C_linalg_downwindmatrix.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <ml_utils.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/08|
 *----------------------------------------------------------------------*/
CORE::LINALG::DownwindMatrix::DownwindMatrix(Teuchos::RCP<Epetra_CrsMatrix> A, const int nv,
    const int np, const double tau, const int outlevel)
    : outlevel_(outlevel), nv_(nv), np_(np), bs_(nv + np), tau_(tau)
{
  Setup(*A);
  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 03/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::DownwindMatrix::Setup(const Epetra_CrsMatrix& A)
{
  Teuchos::Time time("", true);
  if (!A.Filled()) FOUR_C_THROW("Input matrix has to be fill_complete");
  const int numdofrows = A.RowMap().NumMyElements();
  const int bsize = bs_;
  const int vsize = nv_;
  if (numdofrows % bsize) FOUR_C_THROW("Local number of matrix rows does not divide by block size");
  const int numnoderows = numdofrows / bsize;

  // compute a nodal row map
  Teuchos::RCP<Epetra_Map> onoderowmap;
  {
    std::vector<int> gnodeids(numnoderows);
    int count = 0;
    int i;
    for (i = 0; i < numdofrows; ++i)
    {
      const int gdofid = A.RowMap().GID(i);
      //      if (gdofid%bsize)
      //        FOUR_C_THROW("dof id is not a multiple of block size");
      //      const int gnodeid = gdofid/bsize;
      const int gnodeid = gdofid;
      gnodeids[count++] = gnodeid;
      i += (bsize - 1);
    }
    if (count != numnoderows) FOUR_C_THROW("# nodes wrong: %d != %d", count, numnoderows);
    onoderowmap = Teuchos::rcp(new Epetra_Map(-1, numnoderows, gnodeids.data(), 0, A.Comm()));
  }


  // compute a nodal block weighted graph
  Teuchos::RCP<Epetra_CrsMatrix> onodegraph;
  {
    const int maxnumentries = A.MaxNumEntries();
    Teuchos::RCP<SparseMatrix> tmp = Teuchos::rcp(new SparseMatrix(*onoderowmap, maxnumentries));
    std::vector<int> indices(maxnumentries);
    std::vector<double> values(maxnumentries);
    for (int i = 0; i < numdofrows; ++i)
    {
      const int gdofrow = A.RowMap().GID(i);
      //      if (gdofrow%bsize) FOUR_C_THROW("Row map is not a multiple of bsize");
      //      const int gnoderow = gdofrow / bsize;
      const int gnoderow = gdofrow;

      for (int ii = 0; ii < vsize; ++ii)
      {
        int iii = i + ii;
        int numentries;
        int err = A.ExtractMyRowCopy(iii, maxnumentries, numentries, values.data(), indices.data());
        if (err) FOUR_C_THROW("Epetra_CrsMatrix::ExtractMyRowCopy returned err=%d", err);
        for (int j = 0; j < numentries; ++j) indices[j] = A.ColMap().GID(indices[j]);
        ML_az_sort(indices.data(), numentries, nullptr, values.data());
        for (int j = 0; j < numentries; ++j)
        {
          const int gdofcol = indices[j];
          //          if (gdofcol%bsize)
          //            FOUR_C_THROW("Col map is not a multiple of bsize indices[j]=%d
          //            bsize=%d",indices[j],bsize);

          //          const int gnodecol = gdofcol / bsize;
          const int gnodecol = gdofcol;
          if (gnodecol == gnoderow)
          {
            j += (bsize - 1);
            continue;
          }

          double sum = 0.0;
          for (int jj = 0; jj < vsize; ++jj)
          {
            int jjj = jj + j;
            sum += values[jjj] * values[jjj];
          }
          if (sum != 0.0) tmp->Assemble(sum, gnoderow, gnodecol);
          j += (bsize - 1);
        }
      }
      i += (bsize - 1);
    }
    tmp->Complete(*onoderowmap, *onoderowmap);
    onodegraph = tmp->EpetraMatrix();
  }

  // scale nodal block weighted graph by sqrt and compute number of inflows for each node
  Epetra_Vector oaverweight(onodegraph->RowMap(), false);
  Epetra_IntVector oninflow(onodegraph->RowMap(), false);
  {
    for (int i = 0; i < onodegraph->NumMyRows(); ++i)
    {
      int numentries;
      double* values;
      onodegraph->ExtractMyRowView(i, numentries, values);
      if (numentries == 0)
      {
        oaverweight[i] = 0.0;
        oninflow[i] = 0;
        continue;
      }
      double sum = 0.0;
      for (int j = 0; j < numentries; ++j)
      {
        values[j] = sqrt(values[j]);
        sum += values[j];
      }
      oaverweight[i] = sum / numentries;
      oninflow[i] = numentries;
    }
  }

  // create directed graph
  Teuchos::RCP<Epetra_CrsMatrix> nnodegraph;
  {
    // create a transposed of the full graph
    Teuchos::RCP<Epetra_CrsMatrix> onodegrapht = CORE::LINALG::Transpose(onodegraph);
    // create a new graph that will store the directed graph
    Teuchos::RCP<SparseMatrix> tmp =
        Teuchos::rcp(new SparseMatrix(*onoderowmap, (int)(onodegraph->MaxNumEntries())));
    const Epetra_Map& rowmap = onodegraph->RowMap();
    const Epetra_Map& colmap = onodegraph->ColMap();
    const Epetra_Map& colmapt = onodegrapht->ColMap();
    for (int i = 0; i < rowmap.NumMyElements(); ++i)
    {
      // Dirichlet BC
      if (oninflow[i] == 0) continue;
      const int grnode = rowmap.GID(i);

      int numentries;
      int* indices;
      double* values;
      onodegraph->ExtractMyRowView(i, numentries, values, indices);

      int numentriest;
      int* indicest;
      double* valuest;
      onodegrapht->ExtractMyRowView(i, numentriest, valuest, indicest);

      for (int j = 0; j < numentries; ++j)
      {
        const int gcnode = colmap.GID(indices[j]);
        if (!rowmap.MyGID(gcnode)) continue;
        bool foundit = false;
        int k;
        for (k = 0; k < numentriest; ++k)
        {
          const int gcnodet = colmapt.GID(indicest[k]);
          if (gcnodet == gcnode)
          {
            foundit = true;
            break;
          }
        }
        if (!foundit || (values[j] >= valuest[k])) tmp->Assemble(values[j], grnode, gcnode);
      }
    }  // for (int i=0; i<rowmap.NumMyElements(); ++i)
    tmp->Complete(rowmap, rowmap);
    nnodegraph = tmp->EpetraMatrix();
  }
  // coarsen directed graph
  {
    // create a new graph that will store the directed graph
    Teuchos::RCP<SparseMatrix> tmp =
        Teuchos::rcp(new SparseMatrix(nnodegraph->RowMap(), (int)(nnodegraph->MaxNumEntries())));
    const Epetra_Map& rowmap = nnodegraph->RowMap();
    const Epetra_Map& colmap = onodegraph->ColMap();
    for (int i = 0; i < rowmap.NumMyElements(); ++i)
    {
      // Dirichlet BC
      if (oninflow[i] == 0) continue;
      const int grnode = rowmap.GID(i);
      const double average = oaverweight[i];
      int numentries;
      int* indices;
      double* values;
      onodegraph->ExtractMyRowView(i, numentries, values, indices);

      for (int j = 0; j < numentries; ++j)
      {
        const int gcnode = colmap.GID(indices[j]);
        if (!rowmap.MyGID(gcnode)) continue;
        if (values[j] >= tau_ * average) tmp->Assemble(values[j], grnode, gcnode);
      }
    }  // for (int i=0; i<rowmap.NumMyElements(); ++i)
    tmp->Complete(rowmap, rowmap);
    nnodegraph = tmp->EpetraMatrix();
  }

  Epetra_IntVector index(nnodegraph->RowMap(), false);

  downwind_bey_wittum(*nnodegraph, index, oninflow);
  // downwind_hackbusch(*nnodegraph,index,oninflow);

  // index now specifies the local order in which the rows should be processed
  Teuchos::RCP<Epetra_Map> nnoderowmap;
  {
    Epetra_IntVector gindices(nnodegraph->RowMap(), false);
    const int length = index.MyLength();
    for (int i = 0; i < length; ++i)
    {
      const int gid = nnodegraph->RowMap().GID(index[i]);
      if (gid < 0) FOUR_C_THROW("nnodegraph->RowMap().GID(index[i]) failed");
      gindices[i] = gid;
    }
    nnoderowmap = Teuchos::rcp(new Epetra_Map(-1, length, gindices.Values(), 0, A.Comm()));
  }

  // create a new dof row map which is the final result
  // Also, create a sorted variant which then will be used to build the matrix
  {
    const int mynodelength = nnoderowmap->NumMyElements();
    std::vector<int> gindices(mynodelength * bs_);
    for (int i = 0; i < mynodelength; ++i)
      for (int j = 0; j < bs_; ++j) gindices[i * bs_ + j] = nnoderowmap->GID(i) * bs_ + j;
    ndofrowmap_ =
        Teuchos::rcp(new Epetra_Map(-1, mynodelength * bs_, gindices.data(), 0, A.Comm()));
    ML_az_sort(gindices.data(), mynodelength * bs_, nullptr, nullptr);
    sndofrowmap_ =
        Teuchos::rcp(new Epetra_Map(-1, mynodelength * bs_, gindices.data(), 0, A.Comm()));
  }


  // Allocate an exporter from ndofrowmap_ to sndofrowmap_ and back
  sexporter_ = Teuchos::rcp(new Epetra_Export(*ndofrowmap_, *sndofrowmap_));
  rexporter_ = Teuchos::rcp(new Epetra_Export(*sndofrowmap_, *ndofrowmap_));


  if (!A.Comm().MyPID() && outlevel_)
    std::cout << "                Downwinding Setup time " << time.totalElapsedTime(true) << " s\n"
              << "                nv " << nv_ << " np " << np_ << " bs " << bs_ << " tau " << tau_
              << std::endl;

  return;
}



/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 03/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::DownwindMatrix::downwind_bey_wittum(
    const Epetra_CrsMatrix& nnodegraph, Epetra_IntVector& index, const Epetra_IntVector& oninflow)
{
  const int myrank = nnodegraph.Comm().MyPID();
  index.PutValue(-1);
  int nf = 0;

  // number Dirichlet BCs first
  for (int i = 0; i < nnodegraph.NumMyRows(); ++i)
    if (oninflow[i] == 0)
    {
      index[i] = nf;
      ++nf;
    }

  // do downwind numbering
  for (int i = 0; i < nnodegraph.NumMyRows(); ++i)
    if (index[i] < 0) set_f(i, nf, index, nnodegraph, 0);
  if (outlevel_)
  {
    nnodegraph.Comm().Barrier();
    if (!myrank) std::cout << "                Downwinding:" << std::endl;
    nnodegraph.Comm().Barrier();
    std::cout << "                Proc " << myrank << " lastindex " << index.MyLength()
              << " lastdownwind " << nf << std::endl;
  }

  // number everything that's left over in old order
  for (int i = 0; i < nnodegraph.NumMyRows(); ++i)
    if (index[i] < 0)
    {
      index[i] = nf;
      ++nf;
    }
  if (nf != index.MyLength()) FOUR_C_THROW("Local number of nodes wrong");
  return;
}

/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 03/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::DownwindMatrix::downwind_hackbusch(
    const Epetra_CrsMatrix& nnodegraph, Epetra_IntVector& index, const Epetra_IntVector& oninflow)
{
  const int myrank = nnodegraph.Comm().MyPID();
  index.PutValue(-1);
  int nf = 0;
  int nl = index.MyLength() - 1;
  // number Dirichlet BCs first
  for (int i = 0; i < nnodegraph.NumMyRows(); ++i)
    if (oninflow[i] == 0)
    {
      index[i] = nf;
      ++nf;
    }

  // do down- and upwind numbering
  for (int i = 0; i < nnodegraph.NumMyRows(); ++i)
  {
    if (index[i] < 0) set_f(i, nf, index, nnodegraph, 0);
    if (index[i] < 0) set_l(i, nl, index, nnodegraph, 0);
  }
  if (outlevel_)
  {
    nnodegraph.Comm().Barrier();
    if (!myrank) std::cout << "                Downwinding:" << std::endl;
    nnodegraph.Comm().Barrier();
    std::cout << "                Proc " << myrank << " lastindex " << index.MyLength()
              << " lastdownwind " << nf - 1 << " lastupwind " << nl + 1 << std::endl;
  }
  // number everything that's left over in old order
  for (int i = 0; i < nnodegraph.NumMyRows(); ++i)
    if (index[i] < 0)
    {
      index[i] = nf;
      ++nf;
    }
  return;
}


/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 03/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::DownwindMatrix::set_f(
    const int i, int& nf, Epetra_IntVector& index, const Epetra_CrsMatrix& graph, int rec)
{
  // std::cout << "set_f::Recursion " << rec << "\n"; fflush(stdout);
  const Epetra_Map& rowmap = graph.RowMap();
  const Epetra_Map& colmap = graph.ColMap();
  const int gi = rowmap.GID(i);
  int numentries;
  int* indices;
  double* values;
  graph.ExtractMyRowView(i, numentries, values, indices);
  bool allhaveindex = true;
  for (int j = 0; j < numentries; ++j)
  {
    const int gpre = colmap.GID(indices[j]);
    if (!rowmap.MyGID(gpre)) continue;
    const int lpre = rowmap.LID(gpre);
    if (index[lpre] < 0)
    {
      allhaveindex = false;
      break;
    }
  }
  if (allhaveindex)
  {
    // std::cout << "F: Setting new index " << nf << " to old index " << i << std::endl;
    // std::cout.flush();
    index[i] = nf;
    ++nf;
    for (int k = 0; k < graph.NumMyRows(); ++k)
    {
      if (index[k] >= 0) continue;
      if (is_successor(k, gi, graph)) set_f(k, nf, index, graph, rec + 1);
    }
  }
  return;
}



/*----------------------------------------------------------------------*
 |  (private)                                                mwgee 03/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::DownwindMatrix::set_l(
    const int i, int& nl, Epetra_IntVector& index, const Epetra_CrsMatrix& graph, int rec)
{
  // std::cout << "set_l::Recursion " << rec << "\n"; fflush(stdout);
  const Epetra_Map& rowmap = graph.RowMap();
  const Epetra_Map& colmap = graph.ColMap();
  int numentries;
  int* indices;
  double* values;
  graph.ExtractMyRowView(i, numentries, values, indices);
  bool allhaveindex = true;
  const int gi = rowmap.GID(i);
  for (int k = 0; k < graph.NumMyRows(); ++k)
  {
    if (!is_successor(k, gi, graph)) continue;
    if (index[k] < 0)
    {
      allhaveindex = false;
      break;
    }
  }
  if (allhaveindex)
  {
    // std::cout << "L: Setting new index " << nl << " to old index " << i << std::endl;
    // std::cout.flush();
    index[i] = nl;
    --nl;
    int numentries;
    int* indices;
    double* values;
    graph.ExtractMyRowView(i, numentries, values, indices);
    for (int j = 0; j < numentries; ++j)
    {
      int gj = colmap.GID(indices[j]);
      if (!rowmap.MyGID(gj)) continue;
      int lj = rowmap.LID(gj);
      if (index[lj] >= 0) continue;
      set_l(lj, nl, index, graph, rec + 1);
    }
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
