/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of utility functions for mortar methods

\level 1

\maintainer Matthias Mayr
*/
/*-----------------------------------------------------------------------*/

#include "mortar_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_multiply.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_nurbs_discret/drt_knotvector.H"

/*!
\brief Sort vector in ascending order

This routine is taken from Trilinos MOERTEL package.

\param dlist (in): vector to be sorted (unsorted on input, sorted on output)
\param N (in):     length of vector to be sorted
\param list2 (in): another vector which is sorted accordingly

*/
void MORTAR::Sort(double* dlist, int N, int* list2)
{
  int l, r, j, i, flag;
  int RR2;
  double dRR, dK;

  if (N <= 1) return;

  l = N / 2 + 1;
  r = N - 1;
  l = l - 1;
  dRR = dlist[l - 1];
  dK = dlist[l - 1];

  if (list2 != NULL)
  {
    RR2 = list2[l - 1];
    while (r != 0)
    {
      j = l;
      flag = 1;

      while (flag == 1)
      {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else
        {
          if (j < r + 1)
            if (dlist[j] > dlist[j - 1]) j = j + 1;

          if (dlist[j - 1] > dK)
          {
            dlist[i - 1] = dlist[j - 1];
            list2[i - 1] = list2[j - 1];
          }
          else
          {
            flag = 0;
          }
        }
      }
      dlist[i - 1] = dRR;
      list2[i - 1] = RR2;

      if (l == 1)
      {
        dRR = dlist[r];
        RR2 = list2[r];
        dK = dlist[r];
        dlist[r] = dlist[0];
        list2[r] = list2[0];
        r = r - 1;
      }
      else
      {
        l = l - 1;
        dRR = dlist[l - 1];
        RR2 = list2[l - 1];
        dK = dlist[l - 1];
      }
    }
    dlist[0] = dRR;
    list2[0] = RR2;
  }
  else
  {
    while (r != 0)
    {
      j = l;
      flag = 1;
      while (flag == 1)
      {
        i = j;
        j = j + j;
        if (j > r + 1)
          flag = 0;
        else
        {
          if (j < r + 1)
            if (dlist[j] > dlist[j - 1]) j = j + 1;
          if (dlist[j - 1] > dK)
          {
            dlist[i - 1] = dlist[j - 1];
          }
          else
          {
            flag = 0;
          }
        }
      }
      dlist[i - 1] = dRR;
      if (l == 1)
      {
        dRR = dlist[r];
        dK = dlist[r];
        dlist[r] = dlist[0];
        r = r - 1;
      }
      else
      {
        l = l - 1;
        dRR = dlist[l - 1];
        dK = dlist[l - 1];
      }
    }
    dlist[0] = dRR;
  }

  return;
}

/*----------------------------------------------------------------------*
 | transform the row map of a matrix (GIDs)                   popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowTransformGIDs(
    Teuchos::RCP<const LINALG::SparseMatrix> inmat, Teuchos::RCP<const Epetra_Map> newrowmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat =
      Teuchos::rcp(new LINALG::SparseMatrix(*newrowmap, 100, false, true));

  // transform input matrix to newrowmap
  for (int i = 0; i < (inmat->EpetraMatrix())->NumMyRows(); ++i)
  {
    int NumEntries = 0;
    double* Values;
    int* Indices;
    int err = (inmat->EpetraMatrix())->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err != 0) dserror("ExtractMyRowView error: %d", err);

    // pull indices back to global
    std::vector<int> idx(NumEntries);
    for (int j = 0; j < NumEntries; ++j)
    {
      idx[j] = (inmat->ColMap()).GID(Indices[j]);
    }

    err = (outmat->EpetraMatrix())
              ->InsertGlobalValues(
                  newrowmap->GID(i), NumEntries, const_cast<double*>(Values), &idx[0]);
    if (err < 0) dserror("InsertGlobalValues error: %d", err);
  }

  // complete output matrix
  outmat->Complete(inmat->DomainMap(), *newrowmap);

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the column map of a matrix (GIDs)                popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixColTransformGIDs(
    Teuchos::RCP<const LINALG::SparseMatrix> inmat, Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat =
      Teuchos::rcp(new LINALG::SparseMatrix(inmat->RowMap(), 100, false, true));

  // mapping of column gids
  std::map<int, int> gidmap;
  DRT::Exporter ex(inmat->DomainMap(), inmat->ColMap(), inmat->Comm());
  for (int i = 0; i < inmat->DomainMap().NumMyElements(); ++i)
    gidmap[inmat->DomainMap().GID(i)] = newdomainmap->GID(i);
  ex.Export(gidmap);

  // transform input matrix to newdomainmap
  for (int i = 0; i < (inmat->EpetraMatrix())->NumMyRows(); ++i)
  {
    int NumEntries = 0;
    double* Values;
    int* Indices;
    int err = (inmat->EpetraMatrix())->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err != 0) dserror("ExtractMyRowView error: %d", err);
    std::vector<int> idx;
    std::vector<double> vals;
    idx.reserve(NumEntries);
    vals.reserve(NumEntries);

    for (int j = 0; j < NumEntries; ++j)
    {
      int gid = (inmat->ColMap()).GID(Indices[j]);
      std::map<int, int>::const_iterator iter = gidmap.find(gid);
      if (iter != gidmap.end())
      {
        idx.push_back(iter->second);
        vals.push_back(Values[j]);
      }
      else
        dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
    }

    Values = &vals[0];
    NumEntries = vals.size();
    err = (outmat->EpetraMatrix())
              ->InsertGlobalValues(
                  inmat->RowMap().GID(i), NumEntries, const_cast<double*>(Values), &idx[0]);
    if (err < 0) dserror("InsertGlobalValues error: %d", err);
  }

  // complete output matrix
  outmat->Complete(*newdomainmap, inmat->RowMap());

  return outmat;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::CreateNewColMap(const LINALG::SparseMatrix& mat, const Epetra_Map& newdomainmap,
    Teuchos::RCP<Epetra_Map>& newcolmap)
{
  if (not mat.Filled()) dserror("Matrix must be filled!");

  if (not newcolmap.is_null() and mat.ColMap().SameAs(*newcolmap)) return;

  // reset old no longer correct column map
  newcolmap = Teuchos::null;

  // mapping of column gids
  std::map<int, int> gidmap;
  DRT::Exporter exDomain2Col(mat.DomainMap(), mat.ColMap(), mat.Comm());

  const int nummyelements = mat.DomainMap().NumMyElements();
  if (nummyelements != newdomainmap.NumMyElements())
    dserror("NumMyElements must be the same on each proc!");

  const int* old_gids = mat.DomainMap().MyGlobalElements();
  const int* new_gids = newdomainmap.MyGlobalElements();

  for (int i = 0; i < nummyelements; ++i) gidmap[old_gids[i]] = new_gids[i];

  exDomain2Col.Export(gidmap);

  std::vector<int> my_col_gids(gidmap.size(), -1);
  for (std::map<int, int>::const_iterator cit = gidmap.begin(); cit != gidmap.end(); ++cit)
  {
    const int lid = mat.ColMap().LID(cit->first);
    if (lid == -1)
      dserror("Couldn't find the GID %d in the old column map on proc %d.", cit->first,
          mat.Comm().MyPID());

    my_col_gids[lid] = cit->second;
  }

  newcolmap = Teuchos::rcp(new Epetra_Map(mat.ColMap().NumGlobalElements(),
      static_cast<int>(my_col_gids.size()), &my_col_gids[0], 0, mat.Comm()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::ReplaceColumnAndDomainMap(LINALG::SparseMatrix& mat, const Epetra_Map& newdomainmap,
    Teuchos::RCP<Epetra_Map>* const newcolmap_ptr)
{
  if (not mat.Filled()) dserror("Matrix must be filled!");

  Teuchos::RCP<Epetra_Map> newcolmap = Teuchos::null;
  if (newcolmap_ptr)
  {
    CreateNewColMap(mat, newdomainmap, *newcolmap_ptr);
    newcolmap = *newcolmap_ptr;
  }
  else
    CreateNewColMap(mat, newdomainmap, newcolmap);

  int err = mat.EpetraMatrix()->ReplaceColMap(*newcolmap);
  if (err) dserror("ReplaceColMap failed! ( err = %d )", err);

  Epetra_Import importer(*newcolmap, newdomainmap);

  err = mat.EpetraMatrix()->ReplaceDomainMapAndImporter(newdomainmap, &importer);
  if (err) dserror("ReplaceDomainMapAndImporter failed! ( err = %d )", err);
}

/*----------------------------------------------------------------------*
 | transform the row and column maps of a matrix (GIDs)       popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowColTransformGIDs(
    Teuchos::RCP<const LINALG::SparseMatrix> inmat, Teuchos::RCP<const Epetra_Map> newrowmap,
    Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat =
      Teuchos::rcp(new LINALG::SparseMatrix(*newrowmap, 100, true, true));

  // mapping of column gids
  std::map<int, int> gidmap;
  DRT::Exporter ex(inmat->DomainMap(), inmat->ColMap(), inmat->Comm());
  for (int i = 0; i < inmat->DomainMap().NumMyElements(); ++i)
    gidmap[inmat->DomainMap().GID(i)] = newdomainmap->GID(i);
  ex.Export(gidmap);

  // transform input matrix to newrowmap and newdomainmap
  for (int i = 0; i < (inmat->EpetraMatrix())->NumMyRows(); ++i)
  {
    int NumEntries = 0;
    double* Values;
    int* Indices;
    int err = (inmat->EpetraMatrix())->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err != 0) dserror("ExtractMyRowView error: %d", err);
    std::vector<int> idx;
    std::vector<double> vals;
    idx.reserve(NumEntries);
    vals.reserve(NumEntries);

    for (int j = 0; j < NumEntries; ++j)
    {
      int gid = (inmat->ColMap()).GID(Indices[j]);
      std::map<int, int>::const_iterator iter = gidmap.find(gid);
      if (iter != gidmap.end())
      {
        idx.push_back(iter->second);
        vals.push_back(Values[j]);
      }
      else
        dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
    }

    Values = &vals[0];
    NumEntries = vals.size();
    err = (outmat->EpetraMatrix())
              ->InsertGlobalValues(
                  newrowmap->GID(i), NumEntries, const_cast<double*>(Values), &idx[0]);
    if (err < 0) dserror("InsertGlobalValues error: %d", err);
  }

  // complete output matrix
  outmat->Complete(*newdomainmap, *newrowmap);

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the row map of a matrix                          popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowTransform(
    Teuchos::RCP<const LINALG::SparseMatrix> inmat, Teuchos::RCP<const Epetra_Map> newrowmap)
{
  // redistribute input matrix
  Teuchos::RCP<Epetra_CrsMatrix> permmat = Redistribute(*inmat, *newrowmap, inmat->DomainMap());

  // output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat =
      Teuchos::rcp(new LINALG::SparseMatrix(permmat, LINALG::Copy, true));

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the column map of a matrix                       popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixColTransform(
    Teuchos::RCP<const LINALG::SparseMatrix> inmat, Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(*inmat));

  // complete output matrix
  outmat->UnComplete();
  outmat->Complete(*newdomainmap, inmat->RowMap());

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the row and column maps of a matrix              popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowColTransform(
    Teuchos::RCP<const LINALG::SparseMatrix> inmat, Teuchos::RCP<const Epetra_Map> newrowmap,
    Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // redistribute input matrix
  Teuchos::RCP<Epetra_CrsMatrix> permmat = Redistribute(*inmat, *newrowmap, *newdomainmap);

  // output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat =
      Teuchos::rcp(new LINALG::SparseMatrix(permmat, LINALG::Copy, false));

  return outmat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> MORTAR::Redistribute(
    const LINALG::SparseMatrix& src, const Epetra_Map& permrowmap, const Epetra_Map& permdomainmap)
{
  Teuchos::RCP<Epetra_Export> exporter = Teuchos::rcp(new Epetra_Export(permrowmap, src.RowMap()));

  Teuchos::RCP<Epetra_CrsMatrix> permsrc =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, permrowmap, src.MaxNumEntries()));
  int err = permsrc->Import(*src.EpetraMatrix(), *exporter, Insert);
  if (err) dserror("Import failed with err=%d", err);

  permsrc->FillComplete(permdomainmap, permrowmap);
  return permsrc;
}

/*----------------------------------------------------------------------*
 |  Sort points to obtain final clip polygon                  popp 11/08|
 *----------------------------------------------------------------------*/
int MORTAR::SortConvexHullPoints(bool out, Epetra_SerialDenseMatrix& transformed,
    std::vector<Vertex>& collconvexhull, std::vector<Vertex>& respoly, double& tol)
{
  //**********************************************************************
  // - this yields the final clip polygon
  // - sanity of the generated output is checked
  //**********************************************************************

  // (1) Find point with smallest x-value
  // (if more than 1 point with identical x-value exists, choose the one with the smallest y-value)

  // initialize starting point
  int startindex = 0;
  double startpoint[2] = {transformed(0, 0), transformed(1, 0)};

  int np = (int)collconvexhull.size();
  for (int i = 1; i < np; ++i)
  {
    if (transformed(0, i) < startpoint[0])
    {
      startpoint[0] = transformed(0, i);
      startpoint[1] = transformed(1, i);
      startindex = i;
    }
    else if (transformed(0, i) == startpoint[0])
    {
      if (transformed(1, i) < startpoint[1])
      {
        startpoint[1] = transformed(1, i);
        startindex = i;
      }
    }
    else
    {
      // do nothing: starting point did not change
    }
  }

  if (out)
    std::cout << "Start of convex hull: Index " << startindex << "\t" << startpoint[0] << "\t"
              << startpoint[1] << std::endl;

  // (2) Sort remaining points ascending w.r.t their angle with the y-axis
  // (if more than 1 point with identical angle exists, sort ascending w.r.t. their y-value)
  std::vector<double> cotangle(0);
  std::vector<double> yvalues(0);
  std::vector<int> sorted(0);
  std::vector<int> onxline(0);

  for (int i = 0; i < np; ++i)
  {
    // do nothing for starting point
    if (i == startindex) continue;

    // compute angle and store
    double xdiff = transformed(0, i) - startpoint[0];
    double ydiff = transformed(1, i) - startpoint[1];

    if (xdiff < 0) dserror("ERROR: Found point with x < x_start for convex hull!");
    if (xdiff >= tol)
    {
      cotangle.push_back(ydiff / xdiff);
      sorted.push_back(i);
    }
    else
    {
      // these points need further investigation
      onxline.push_back(i);
    }
  }

  // check points on x-line with starting point and only add
  // those with min and max value in y-direction
  {
    double y_max = std::numeric_limits<double>::min();
    double y_min = std::numeric_limits<double>::max();
    int i_max = -1;
    int i_min = -1;
    for (size_t i = 0; i < onxline.size(); ++i)
    {
      const double yval = transformed(1, onxline[i]) - startpoint[1];
      if (yval < y_min && yval < 0.0)
      {
        y_min = yval;
        i_min = onxline[i];
      }
      else if (yval > y_max && yval > 0.0)
      {
        y_max = yval;
        i_max = onxline[i];
      }
    }
    if (i_max > -1)
    {
      cotangle.push_back(std::numeric_limits<double>::max());
      sorted.push_back(i_max);
    }
    if (i_min > -1)
    {
      cotangle.push_back(-std::numeric_limits<double>::max());
      sorted.push_back(i_min);
    }
  }

  // start index not yet included
  np = (int)sorted.size() + 1;

  if (out)
  {
    std::cout << "Unsorted convex hull:\n";
    std::cout << "Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1]
              << std::endl;
    for (int i = 0; i < np - 1; ++i)
      std::cout << "Index " << sorted[i] << "\t" << transformed(0, sorted[i]) << "\t"
                << transformed(1, sorted[i]) << "\t" << cotangle[i] << std::endl;
  }

  // check if sizes are correct
  if ((int)cotangle.size() != np - 1) dserror("ERROR: Size went wrong for cot angle!");

  // now sort descending w.r.t cotangle = ascending w.r.t angle
  MORTAR::Sort(&cotangle[0], np - 1, &sorted[0]);
  std::reverse(cotangle.begin(), cotangle.end());
  std::reverse(sorted.begin(), sorted.end());

  // get associated y-values
  for (int i = 0; i < np - 1; ++i) yvalues.push_back(transformed(1, sorted[i]));

  // now sort ascending w.r.t value wherever angles are identical
  // (bubblesort: we might need np-2 rounds if all np-1 angles identical)
  for (int round = 0; round < np - 2; ++round)
    for (int i = 0; i < np - 2; ++i)
      if (cotangle[i] == cotangle[i + 1])
        if (yvalues[i] > yvalues[i + 1])
        {
          MORTAR::Swap(cotangle[i], cotangle[i + 1]);
          MORTAR::Swap(yvalues[i], yvalues[i + 1]);
          MORTAR::Swap(sorted[i], sorted[i + 1]);
        }

  if (out)
  {
    std::cout << "Sorted convex hull:\n";
    std::cout << "Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1]
              << std::endl;
    for (int i = 0; i < np - 1; ++i)
      std::cout << "Index " << sorted[i] << "\t" << transformed(0, sorted[i]) << "\t"
                << transformed(1, sorted[i]) << "\t" << cotangle[i] << std::endl;
  }

  // (3) Go through sorted list of points
  // (keep adding points as long as the last 3 points rotate clockwise)
  // (if 3 points rotate counter-clockwise, do NOT add current point and continue)

  // always push pack starting point
  Vertex* current = &collconvexhull[startindex];
  respoly.push_back(Vertex(current->Coord(), current->VType(), current->Nodeids(), NULL, NULL,
      false, false, NULL, -1.0));

  // number of points removed from convex hull
  int removed = (int)collconvexhull.size() - np;

  // go through sorted list and check for clockwise rotation
  std::vector<bool> haveremovedthis(np - 1);
  for (int i = 0; i < np - 1; ++i) haveremovedthis[i] = false;

  for (int i = 0; i < np - 1; ++i)
  {
    double edge1[2] = {0.0, 0.0};
    double edge2[2] = {0.0, 0.0};

    // first triple
    if (i == 0)
    {
      edge1[0] = transformed(0, sorted[0]) - startpoint[0];
      edge1[1] = transformed(1, sorted[0]) - startpoint[1];
      edge2[0] = transformed(0, sorted[1]) - transformed(0, sorted[0]);
      edge2[1] = transformed(1, sorted[1]) - transformed(1, sorted[0]);
    }

    // standard case
    else if (i < np - 2)
    {
      // go back and find first non-removed partner
      bool foundpartner = false;
      int k = i - 1;

      while (foundpartner == false)
      {
        // found non-removed partner
        if (haveremovedthis[k] == false)
        {
          edge1[0] = transformed(0, sorted[i]) - transformed(0, sorted[k]);
          edge1[1] = transformed(1, sorted[i]) - transformed(1, sorted[k]);
          edge2[0] = transformed(0, sorted[i + 1]) - transformed(0, sorted[i]);
          edge2[1] = transformed(1, sorted[i + 1]) - transformed(1, sorted[i]);
          foundpartner = true;
        }
        else
        {
          // decrease counter
          k -= 1;

          // use starting point if all in between removed
          if (k < 0)
          {
            edge1[0] = transformed(0, sorted[i]) - startpoint[0];
            edge1[1] = transformed(1, sorted[i]) - startpoint[1];
            edge2[0] = transformed(0, sorted[i + 1]) - transformed(0, sorted[i]);
            edge2[1] = transformed(1, sorted[i + 1]) - transformed(1, sorted[i]);
            foundpartner = true;
          }
        }
      }
    }

    // last triple
    else /* if i = np-1 */
    {
      // go back and find first non-removed partner
      bool foundpartner = false;
      int k = i - 1;

      while (foundpartner == false)
      {
        // found non-removed partner
        if (haveremovedthis[k] == false)
        {
          edge1[0] = transformed(0, sorted[i]) - transformed(0, sorted[k]);
          edge1[1] = transformed(1, sorted[i]) - transformed(1, sorted[k]);
          edge2[0] = startpoint[0] - transformed(0, sorted[i]);
          edge2[1] = startpoint[1] - transformed(1, sorted[i]);
          foundpartner = true;
        }
        else
        {
          // decrease counter
          k -= 1;

          // use starting point if all in between removed
          if (k < 0)
          {
            edge1[0] = transformed(0, sorted[i]) - startpoint[0];
            edge1[1] = transformed(1, sorted[i]) - startpoint[1];
            edge2[0] = startpoint[0] - transformed(0, sorted[i]);
            edge2[1] = startpoint[1] - transformed(1, sorted[i]);
            foundpartner = true;
          }
        }
      }
    }

    // check for clockwise rotation
    double cw = edge1[0] * edge2[1] - edge1[1] * edge2[0];

    // add point to convex hull if clockwise triple
    // (use tolerance to remove almost straight lines of 3 points)
    if (cw <= -tol)
    {
      Vertex* current = &collconvexhull[sorted[i]];
      respoly.push_back(Vertex(current->Coord(), current->VType(), current->Nodeids(), NULL, NULL,
          false, false, NULL, -1.0));
    }
    // mark vertex as "removed" if counter-clockwise triple
    else
    {
      removed++;
      haveremovedthis[i] = true;
    }
  }

  return removed;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MORTAR::UTILS::CreateVolumeGhosting(const DRT::Discretization& dis_src,
    const std::vector<std::string> dis_tar, std::vector<std::pair<int, int>> material_links,
    bool check_on_in, bool check_on_exit)
{
  if (dis_tar.size() == 0) return;

  DRT::Problem* problem = DRT::Problem::Instance();
  std::vector<Teuchos::RCP<DRT::Discretization>> voldis;
  for (int name = 0; name < (int)dis_tar.size(); ++name)
    voldis.push_back(problem->GetDis(dis_tar.at(name)));

  if (check_on_in)
    for (int c = 1; c < (int)voldis.size(); ++c)
      if (voldis.at(c)->ElementRowMap()->SameAs(*voldis.at(0)->ElementRowMap()) == false)
        dserror("row maps on input do not coincide");

  const Epetra_Map* ielecolmap = dis_src.ElementColMap();

  // 1 Ghost all Volume Element + Nodes,for all col elements in dis_src
  for (unsigned disidx = 0; disidx < voldis.size(); ++disidx)
  {
    std::vector<int> rdata;

    // Fill rdata with existing colmap

    const Epetra_Map* elecolmap = voldis[disidx]->ElementColMap();
    const Teuchos::RCP<Epetra_Map> allredelecolmap =
        LINALG::AllreduceEMap(*voldis[disidx]->ElementRowMap());

    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      rdata.push_back(gid);
    }

    // Find elements, which are ghosted on the interface but not in the volume discretization
    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      DRT::Element* ele = dis_src.gElement(gid);
      if (!ele) dserror("ERROR: Cannot find element with gid %", gid);
      DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(ele);
      if (!faceele) dserror("source element is not a face element");
      int volgid = faceele->ParentElementId();
      // Ghost the parent element additionally
      if (elecolmap->LID(volgid) == -1 &&
          allredelecolmap->LID(volgid) !=
              -1)  // Volume Discretization has not Element on this proc but on another
        rdata.push_back(volgid);
    }

    // re-build element column map
    Teuchos::RCP<Epetra_Map> newelecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, voldis[disidx]->Comm()));
    rdata.clear();

    // redistribute the volume discretization according to the
    // new (=old) element column layout & and ghost also nodes!
    voldis[disidx]->ExtendedGhosting(*newelecolmap, true, true, true, false);  // no check!!!
  }

  // 2 Reconnect Face Element -- Parent Element Pointers to first dis in dis_tar
  {
    const Epetra_Map* elecolmap = voldis[0]->ElementColMap();

    for (int i = 0; i < ielecolmap->NumMyElements(); ++i)
    {
      int gid = ielecolmap->GID(i);

      DRT::Element* ele = dis_src.gElement(gid);
      if (!ele) dserror("ERROR: Cannot find element with gid %", gid);
      DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(ele);
      if (!faceele) dserror("source element is not a face element");
      int volgid = faceele->ParentElementId();

      if (elecolmap->LID(volgid) == -1)  // Volume Discretization has not Element
        dserror("CreateVolumeGhosting: Element %d does not exist on this Proc!", volgid);

      DRT::Element* vele = voldis[0]->gElement(volgid);
      if (!vele) dserror("ERROR: Cannot find element with gid %", volgid);

      faceele->SetParentMasterElement(vele, faceele->FaceParentNumber());
    }
  }

  if (check_on_exit)
    for (int c = 1; c < (int)voldis.size(); ++c)
    {
      if (voldis.at(c)->ElementRowMap()->SameAs(*voldis.at(0)->ElementRowMap()) == false)
        dserror("row maps on exit do not coincide");
      if (voldis.at(c)->ElementColMap()->SameAs(*voldis.at(0)->ElementColMap()) == false)
        dserror("col maps on exit do not coincide");
    }

  // 3 setup material pointers between newly ghosted elements
  for (std::vector<std::pair<int, int>>::const_iterator m = material_links.begin();
       m != material_links.end(); ++m)
  {
    Teuchos::RCP<DRT::Discretization> dis_src_mat = voldis.at(m->first);
    Teuchos::RCP<DRT::Discretization> dis_tar_mat = voldis.at(m->second);

    for (int i = 0; i < dis_tar_mat->NumMyColElements(); ++i)
    {
      DRT::Element* targetele = dis_tar_mat->lColElement(i);
      const int gid = targetele->Id();

      DRT::Element* sourceele = dis_src_mat->gElement(gid);

      targetele->AddMaterial(sourceele->Material());
    }
  }
}



/*----------------------------------------------------------------------*
 |  Prepare mortar element for nurbs-case                    farah 11/14|
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::PrepareNURBSElement(DRT::Discretization& discret,
    Teuchos::RCP<DRT::Element> ele, Teuchos::RCP<MORTAR::MortarElement> cele, int dim)
{
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discret));

  Teuchos::RCP<DRT::NURBS::Knotvector> knots = (*nurbsdis).GetKnotVector();
  std::vector<Epetra_SerialDenseVector> parentknots(dim);
  std::vector<Epetra_SerialDenseVector> mortarknots(dim - 1);

  double normalfac = 0.0;
  Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele, true);
  bool zero_size = knots->GetBoundaryEleAndParentKnots(parentknots, mortarknots, normalfac,
      faceele->ParentMasterElement()->Id(), faceele->FaceMasterNumber());

  // store nurbs specific data to node
  cele->ZeroSized() = zero_size;
  cele->Knots() = mortarknots;
  cele->NormalFac() = normalfac;

  return;
}


/*----------------------------------------------------------------------*
 |  Prepare mortar node for nurbs-case                       farah 11/14|
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::PrepareNURBSNode(DRT::Node* node, Teuchos::RCP<MORTAR::MortarNode> mnode)
{
  DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(node);

  mnode->NurbsW() = cp->W();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::MortarMatrixCondensation(Teuchos::RCP<LINALG::SparseMatrix>& k,
    const Teuchos::RCP<const LINALG::SparseMatrix>& p_row,
    const Teuchos::RCP<const LINALG::SparseMatrix>& p_col)
{
  // prepare maps
  Teuchos::RCP<Epetra_Map> gsrow =
      Teuchos::rcp_const_cast<Epetra_Map>(Teuchos::rcpFromRef<const Epetra_Map>(p_row->RangeMap()));
  Teuchos::RCP<Epetra_Map> gmrow = Teuchos::rcp_const_cast<Epetra_Map>(
      Teuchos::rcpFromRef<const Epetra_Map>(p_row->DomainMap()));
  Teuchos::RCP<Epetra_Map> gsmrow = LINALG::MergeMap(gsrow, gmrow, false);
  Teuchos::RCP<Epetra_Map> gnrow = LINALG::SplitMap(k->RangeMap(), *gsmrow);

  Teuchos::RCP<Epetra_Map> gscol =
      Teuchos::rcp_const_cast<Epetra_Map>(Teuchos::rcpFromRef<const Epetra_Map>(p_col->RangeMap()));
  Teuchos::RCP<Epetra_Map> gmcol = Teuchos::rcp_const_cast<Epetra_Map>(
      Teuchos::rcpFromRef<const Epetra_Map>(p_col->DomainMap()));
  Teuchos::RCP<Epetra_Map> gsmcol = LINALG::MergeMap(gscol, gmcol, false);
  Teuchos::RCP<Epetra_Map> gncol = LINALG::SplitMap(k->DomainMap(), *gsmcol);

  /*--------------------------------------------------------------------*/
  /* Split kteff into 3x3 block matrix                                  */
  /*--------------------------------------------------------------------*/
  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kss = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> ksm = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> ksn = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> kms = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> kmm = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> kmn = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> kns = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> knm = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> knn = Teuchos::null;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<LINALG::SparseMatrix> ksmsm = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> ksmn = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> knsm = Teuchos::null;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1 = Teuchos::null;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2 = Teuchos::null;

  // split
  LINALG::SplitMatrix2x2(k, gsmrow, gnrow, gsmcol, gncol, ksmsm, ksmn, knsm, knn);
  LINALG::SplitMatrix2x2(ksmsm, gsrow, gmrow, gscol, gmcol, kss, ksm, kms, kmm);
  LINALG::SplitMatrix2x2(ksmn, gsrow, gmrow, gncol, tempmap, ksn, tempmtx1, kmn, tempmtx2);
  LINALG::SplitMatrix2x2(knsm, gnrow, tempmap, gscol, gmcol, kns, knm, tempmtx1, tempmtx2);

  Teuchos::RCP<LINALG::SparseMatrix> kteffnew =
      Teuchos::rcp(new LINALG::SparseMatrix(k->RowMap(), 81, true, false, k->GetMatrixtype()));

  // build new stiffness matrix
  kteffnew->Add(*knn, false, 1.0, 1.0);
  kteffnew->Add(*knm, false, 1.0, 1.0);
  kteffnew->Add(*kmn, false, 1.0, 1.0);
  kteffnew->Add(*kmm, false, 1.0, 1.0);
  kteffnew->Add(*LINALG::Multiply(*kns, false, *p_col, false, true, false, true), false, 1., 1.);
  kteffnew->Add(*LINALG::Multiply(*p_row, true, *ksn, false, true, false, true), false, 1., 1.);
  kteffnew->Add(*LINALG::Multiply(*kms, false, *p_col, false, true, false, true), false, 1., 1.);
  kteffnew->Add(*LINALG::Multiply(*p_row, true, *ksm, false, true, false, true), false, 1., 1.);
  kteffnew->Add(*LINALG::Multiply(*p_row, true,
                    *LINALG::Multiply(*kss, false, *p_col, false, true, false, true), false, true,
                    false, true),
      false, 1., 1.);
  if (p_row == p_col) kteffnew->Add(*LINALG::Eye(*gsrow), false, 1., 1.);

  kteffnew->Complete(k->DomainMap(), k->RangeMap());

  // return new matrix
  k = kteffnew;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::MortarRhsCondensation(
    Teuchos::RCP<Epetra_Vector>& rhs, const Teuchos::RCP<LINALG::SparseMatrix>& p)
{
  // prepare maps
  Teuchos::RCP<Epetra_Map> gsdofrowmap =
      Teuchos::rcp_const_cast<Epetra_Map>(Teuchos::rcpFromRef<const Epetra_Map>(p->RangeMap()));
  Teuchos::RCP<Epetra_Map> gmdofrowmap =
      Teuchos::rcp_const_cast<Epetra_Map>(Teuchos::rcpFromRef<const Epetra_Map>(p->DomainMap()));

  Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap));
  Teuchos::RCP<Epetra_Vector> fm_cond = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap));
  LINALG::Export(*rhs, *fs);
  Teuchos::RCP<Epetra_Vector> fs_full = Teuchos::rcp(new Epetra_Vector(rhs->Map()));
  LINALG::Export(*fs, *fs_full);
  if (rhs->Update(-1., *fs_full, 1.)) dserror("update failed");

  if (p->Multiply(true, *fs, *fm_cond)) dserror("multiply failed");

  Teuchos::RCP<Epetra_Vector> fm_cond_full = Teuchos::rcp(new Epetra_Vector(rhs->Map()));
  LINALG::Export(*fm_cond, *fm_cond_full);
  if (rhs->Update(1., *fm_cond_full, 1.)) dserror("update failed");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::MortarRecover(
    Teuchos::RCP<Epetra_Vector>& inc, const Teuchos::RCP<LINALG::SparseMatrix>& p)
{
  // prepare maps
  Teuchos::RCP<Epetra_Map> gsdofrowmap =
      Teuchos::rcp_const_cast<Epetra_Map>(Teuchos::rcpFromRef<const Epetra_Map>(p->RangeMap()));
  Teuchos::RCP<Epetra_Map> gmdofrowmap =
      Teuchos::rcp_const_cast<Epetra_Map>(Teuchos::rcpFromRef<const Epetra_Map>(p->DomainMap()));

  Teuchos::RCP<Epetra_Vector> m_inc = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap));
  LINALG::Export(*inc, *m_inc);

  Teuchos::RCP<Epetra_Vector> s_inc = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap));
  if (p->Multiply(false, *m_inc, *s_inc)) dserror("multiply failed");
  Teuchos::RCP<Epetra_Vector> s_inc_full = Teuchos::rcp(new Epetra_Vector(inc->Map()));
  LINALG::Export(*s_inc, *s_inc_full);
  if (inc->Update(1., *s_inc_full, 1.)) dserror("update failed");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::MortarMatrixCondensation(Teuchos::RCP<LINALG::BlockSparseMatrixBase>& k,
    const std::vector<Teuchos::RCP<LINALG::SparseMatrix>>& p)
{
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> cond_mat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          k->DomainExtractor(), k->RangeExtractor(), 81, false, true));

  for (int row = 0; row < k->Rows(); ++row)
    for (int col = 0; col < k->Cols(); ++col)
    {
      Teuchos::RCP<LINALG::SparseMatrix> new_matrix =
          Teuchos::rcp(new LINALG::SparseMatrix(k->Matrix(row, col)));
      MortarMatrixCondensation(new_matrix, p.at(row), p.at(col) /*,row!=col*/);
      cond_mat->Assign(row, col, LINALG::Copy, *new_matrix);
    }

  cond_mat->Complete();

  k = cond_mat;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::MortarRhsCondensation(
    Teuchos::RCP<Epetra_Vector>& rhs, const std::vector<Teuchos::RCP<LINALG::SparseMatrix>>& p)
{
  for (unsigned i = 0; i < p.size(); MortarRhsCondensation(rhs, p[i++]))
    ;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::UTILS::MortarRecover(
    Teuchos::RCP<Epetra_Vector>& inc, const std::vector<Teuchos::RCP<LINALG::SparseMatrix>>& p)
{
  for (unsigned i = 0; i < p.size(); MortarRecover(inc, p[i++]))
    ;
}
