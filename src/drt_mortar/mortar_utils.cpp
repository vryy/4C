/*!----------------------------------------------------------------------
\file mortar_utils.cpp

\brief A set of utility functions for mortar methods

\level 1

\maintainer Alexander Popp

*----------------------------------------------------------------------*/

#include "mortar_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"
#include "../linalg/linalg_sparsematrix.H"

/*!
\brief Sort vector in ascending order

This routine is taken from Trilinos MOERTEL package.

\param dlist (in): vector to be sorted (unsorted on input, sorted on output)
\param N (in):     length of vector to be sorted
\param list2 (in): another vector which is sorted accordingly

*/
void MORTAR::Sort(double* dlist, int N, int* list2)
{
  int    l, r, j, i, flag;
  int    RR2;
  double dRR, dK;

  if (N <= 1) return;

  l    = N / 2 + 1;
  r    = N - 1;
  l    = l - 1;
  dRR  = dlist[l - 1];
  dK   = dlist[l - 1];

  if (list2 != NULL) {
     RR2 = list2[l - 1];
     while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
           i = j;
           j = j + j;

           if (j > r + 1)
              flag = 0;
           else {
              if (j < r + 1)
                 if (dlist[j] > dlist[j - 1]) j = j + 1;

              if (dlist[j - 1] > dK) {
                 dlist[ i - 1] = dlist[ j - 1];
                 list2[i - 1] = list2[j - 1];
              }
              else {
                 flag = 0;
              }
           }
        }
        dlist[ i - 1] = dRR;
        list2[i - 1] = RR2;

        if (l == 1) {
           dRR  = dlist [r];
           RR2 = list2[r];
           dK = dlist[r];
           dlist[r ] = dlist[0];
           list2[r] = list2[0];
           r = r - 1;
         }
         else {
            l   = l - 1;
            dRR  = dlist[ l - 1];
            RR2 = list2[l - 1];
            dK   = dlist[l - 1];
         }
      }
      dlist[ 0] = dRR;
      list2[0] = RR2;
   }
   else {
      while (r != 0) {
         j = l;
         flag = 1;
         while (flag == 1) {
            i = j;
            j = j + j;
            if (j > r + 1)
               flag = 0;
            else {
               if (j < r + 1)
                  if (dlist[j] > dlist[j - 1]) j = j + 1;
               if (dlist[j - 1] > dK) {
                  dlist[ i - 1] = dlist[ j - 1];
               }
               else {
                  flag = 0;
               }
            }
         }
         dlist[ i - 1] = dRR;
         if (l == 1) {
            dRR  = dlist [r];
            dK = dlist[r];
            dlist[r ] = dlist[0];
            r = r - 1;
         }
         else {
            l   = l - 1;
            dRR  = dlist[ l - 1];
            dK   = dlist[l - 1];
         }
      }
      dlist[ 0] = dRR;
   }

  return;
}

/*----------------------------------------------------------------------*
 | transform the row map of a matrix (GIDs)                   popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowTransformGIDs(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                                  Teuchos::RCP<const Epetra_Map> newrowmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(*newrowmap,100,false,true));

  // transform input matrix to newrowmap
  for (int i=0; i<(inmat->EpetraMatrix())->NumMyRows(); ++i)
  {
    int NumEntries = 0;
    double *Values;
    int *Indices;
    int err = (inmat->EpetraMatrix())->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0) dserror("ExtractMyRowView error: %d", err);

    // pull indices back to global
    std::vector<int> idx(NumEntries);
    for (int j=0; j<NumEntries; ++j)
    {
      idx[j] = (inmat->ColMap()).GID(Indices[j]);
    }

    err = (outmat->EpetraMatrix())->InsertGlobalValues(newrowmap->GID(i), NumEntries, const_cast<double*>(Values),&idx[0]);
    if (err<0) dserror("InsertGlobalValues error: %d", err);
  }

  // complete output matrix
  outmat->Complete(inmat->DomainMap(),*newrowmap);

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the column map of a matrix (GIDs)                popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixColTransformGIDs(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                                  Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(inmat->RowMap(),100,false,true));

  // mapping of column gids
  std::map<int,int> gidmap;
  DRT::Exporter ex(inmat->DomainMap(),inmat->ColMap(),inmat->Comm());
  for (int i=0; i<inmat->DomainMap().NumMyElements(); ++i)
    gidmap[inmat->DomainMap().GID(i)] = newdomainmap->GID(i);
  ex.Export(gidmap);

  // transform input matrix to newdomainmap
  for (int i=0;i<(inmat->EpetraMatrix())->NumMyRows();++i)
  {
    int NumEntries = 0;
    double *Values;
    int *Indices;
    int err = (inmat->EpetraMatrix())->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0) dserror("ExtractMyRowView error: %d", err);
    std::vector<int> idx;
    std::vector<double> vals;
    idx.reserve(NumEntries);
    vals.reserve(NumEntries);

    for (int j=0;j<NumEntries;++j)
    {
      int gid = (inmat->ColMap()).GID(Indices[j]);
      std::map<int,int>::const_iterator iter = gidmap.find(gid);
      if (iter!=gidmap.end())
      {
        idx.push_back(iter->second);
        vals.push_back(Values[j]);
      }
      else
        dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
    }

    Values = &vals[0];
    NumEntries = vals.size();
    err = (outmat->EpetraMatrix())->InsertGlobalValues(inmat->RowMap().GID(i), NumEntries, const_cast<double*>(Values),&idx[0]);
    if (err<0) dserror("InsertGlobalValues error: %d", err);
  }

  // complete output matrix
  outmat->Complete(*newdomainmap,inmat->RowMap());

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the row and column maps of a matrix (GIDs)       popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowColTransformGIDs(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                                     Teuchos::RCP<const Epetra_Map> newrowmap,
                                                                     Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(*newrowmap,100,true,true));

  // mapping of column gids
  std::map<int,int> gidmap;
  DRT::Exporter ex(inmat->DomainMap(),inmat->ColMap(),inmat->Comm());
  for (int i=0; i<inmat->DomainMap().NumMyElements(); ++i)
    gidmap[inmat->DomainMap().GID(i)] = newdomainmap->GID(i);
  ex.Export(gidmap);

  // transform input matrix to newrowmap and newdomainmap
  for (int i=0;i<(inmat->EpetraMatrix())->NumMyRows();++i)
  {
    int NumEntries = 0;
    double *Values;
    int *Indices;
    int err = (inmat->EpetraMatrix())->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0) dserror("ExtractMyRowView error: %d", err);
    std::vector<int> idx;
    std::vector<double> vals;
    idx.reserve(NumEntries);
    vals.reserve(NumEntries);

    for (int j=0;j<NumEntries;++j)
    {
      int gid = (inmat->ColMap()).GID(Indices[j]);
      std::map<int,int>::const_iterator iter = gidmap.find(gid);
      if (iter!=gidmap.end())
      {
        idx.push_back(iter->second);
        vals.push_back(Values[j]);
      }
      else
        dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
    }

    Values = &vals[0];
    NumEntries = vals.size();
    err = (outmat->EpetraMatrix())->InsertGlobalValues(newrowmap->GID(i), NumEntries, const_cast<double*>(Values),&idx[0]);
    if (err<0) dserror("InsertGlobalValues error: %d", err);
  }

  // complete output matrix
  outmat->Complete(*newdomainmap,*newrowmap);

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the row map of a matrix                          popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowTransform(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                              Teuchos::RCP<const Epetra_Map> newrowmap)
{
  // redistribute input matrix
  Teuchos::RCP<Epetra_CrsMatrix> permmat = Redistribute(*inmat,*newrowmap,inmat->DomainMap());

  // output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(permmat,LINALG::Copy,true));

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the column map of a matrix                       popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixColTransform(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                              Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // initialize output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(*inmat));

  // complete output matrix
  outmat->UnComplete();
  outmat->Complete(*newdomainmap,inmat->RowMap());

  return outmat;
}

/*----------------------------------------------------------------------*
 | transform the row and column maps of a matrix              popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> MORTAR::MatrixRowColTransform(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                                 Teuchos::RCP<const Epetra_Map> newrowmap,
                                                                 Teuchos::RCP<const Epetra_Map> newdomainmap)
{
  // redistribute input matrix
  Teuchos::RCP<Epetra_CrsMatrix> permmat = Redistribute(*inmat,*newrowmap,*newdomainmap);

  // output matrix
  Teuchos::RCP<LINALG::SparseMatrix> outmat = Teuchos::rcp(new LINALG::SparseMatrix(permmat,LINALG::Copy,false));

  return outmat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> MORTAR::Redistribute(const LINALG::SparseMatrix& src,
                                                    const Epetra_Map& permrowmap,
                                                    const Epetra_Map& permdomainmap)
{
  Teuchos::RCP<Epetra_Export> exporter = Teuchos::rcp(new Epetra_Export(permrowmap,src.RowMap()));

  Teuchos::RCP<Epetra_CrsMatrix> permsrc = Teuchos::rcp(new Epetra_CrsMatrix(Copy,permrowmap,src.MaxNumEntries()));
  int err = permsrc->Import(*src.EpetraMatrix(),*exporter,Insert);
  if (err) dserror("Import failed with err=%d",err);

  permsrc->FillComplete(permdomainmap,permrowmap);
  return permsrc;
}

/*----------------------------------------------------------------------*
 |  Sort points to obtain final clip polygon                  popp 11/08|
 *----------------------------------------------------------------------*/
int MORTAR::SortConvexHullPoints(bool out,
                                 Epetra_SerialDenseMatrix& transformed,
                                 std::vector<Vertex>& collconvexhull,
                                 std::vector<Vertex>& respoly,
                                 double& tol)
{
  //**********************************************************************
  // - this yields the final clip polygon
  // - sanity of the generated output is checked
  //**********************************************************************

  // (1) Find point with smallest x-value
  // (if more than 1 point with identical x-value exists, choose the one with the smallest y-value)

  // initialize starting point
  int startindex = 0;
  double startpoint[2] = {transformed(0,0), transformed(1,0)};

  int np = (int)collconvexhull.size();
  for (int i=1;i<np;++i)
  {
    if (transformed(0,i) < startpoint[0])
    {
      startpoint[0] = transformed(0,i);
      startpoint[1] = transformed(1,i);
      startindex = i;
    }
    else if (transformed(0,i) == startpoint[0])
    {
      if (transformed(1,i) < startpoint[1])
      {
        startpoint[1] = transformed(1,i);
        startindex = i;
      }
    }
    else
    {
      // do nothing: starting point did not change
    }
  }

  if (out) std::cout << "Start of convex hull: Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1] << std::endl;

  // (2) Sort remaining points ascending w.r.t their angle with the y-axis
  // (if more than 1 point with identical angle exists, sort ascending w.r.t. their y-value)
  std::vector<double> cotangle(0);
  std::vector<double> yvalues(0);
  std::vector<int> sorted(0);
  std::vector<int> onxline(0);

  for (int i=0;i<np;++i)
  {
    // do nothing for starting point
    if (i==startindex) continue;

    // compute angle and store
    double xdiff = transformed(0,i) - startpoint[0];
    double ydiff = transformed(1,i) - startpoint[1];

    if (xdiff < 0) dserror("ERROR: Found point with x < x_start for convex hull!");
    if (xdiff >= tol)
    {
      cotangle.push_back(ydiff/xdiff);
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
    int i_max=-1;
    int i_min=-1;
    for (size_t i=0;i<onxline.size();++i)
    {
      const double yval = transformed(1,onxline[i]) - startpoint[1];
      if(yval < y_min && yval < 0.0)
      {
        y_min = yval;
        i_min = onxline[i];
      }
      else if(yval > y_max && yval > 0.0)
      {
        y_max = yval;
        i_max = onxline[i];
      }
    }
    if(i_max>-1)
    {
      cotangle.push_back(std::numeric_limits<double>::max());
      sorted.push_back(i_max);
    }
    if(i_min>-1)
    {
      cotangle.push_back(-std::numeric_limits<double>::max());
      sorted.push_back(i_min);
    }
  }

  // start index not yet included
  np = (int)sorted.size()+1;

  if (out)
  {
    std::cout << "Unsorted convex hull:\n";
    std::cout << "Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1] << std::endl;
    for (int i=0;i<np-1;++i)
      std::cout << "Index " << sorted[i] << "\t" << transformed(0,sorted[i]) << "\t" << transformed(1,sorted[i]) << "\t" << cotangle[i] << std::endl;
  }

  // check if sizes are correct
  if ((int)cotangle.size() != np-1) dserror("ERROR: Size went wrong for cot angle!");

  // now sort descending w.r.t cotangle = ascending w.r.t angle
  MORTAR::Sort(&cotangle[0],np-1,&sorted[0]);
  std::reverse(cotangle.begin(), cotangle.end());
  std::reverse(sorted.begin(), sorted.end());

  // get associated y-values
  for (int i=0;i<np-1;++i)
    yvalues.push_back(transformed(1,sorted[i]));

  // now sort ascending w.r.t value wherever angles are identical
  // (bubblesort: we might need np-2 rounds if all np-1 angles identical)
  for (int round=0;round<np-2;++round)
    for (int i=0;i<np-2;++i)
      if (cotangle[i]==cotangle[i+1])
        if (yvalues[i]>yvalues[i+1])
        {
          MORTAR::Swap(cotangle[i],cotangle[i+1]);
          MORTAR::Swap(yvalues[i],yvalues[i+1]);
          MORTAR::Swap(sorted[i],sorted[i+1]);
        }

  if (out)
  {
    std::cout << "Sorted convex hull:\n";
    std::cout << "Index " << startindex << "\t" << startpoint[0] << "\t" << startpoint[1] << std::endl;
    for (int i=0;i<np-1;++i)
      std::cout << "Index " << sorted[i] << "\t" << transformed(0,sorted[i]) << "\t" << transformed(1,sorted[i]) << "\t" << cotangle[i] << std::endl;
  }

  // (3) Go through sorted list of points
  // (keep adding points as long as the last 3 points rotate clockwise)
  // (if 3 points rotate counter-clockwise, do NOT add current point and continue)

  // always push pack starting point
  Vertex* current = &collconvexhull[startindex];
  respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));

  // number of points removed from convex hull
  int removed = (int)collconvexhull.size() - np;

  // go through sorted list and check for clockwise rotation
  std::vector<bool> haveremovedthis(np-1);
  for (int i=0;i<np-1;++i) haveremovedthis[i] = false;

  for (int i=0;i<np-1;++i)
  {
    double edge1[2] = {0.0, 0.0};
    double edge2[2] = {0.0, 0.0};

    // first triple
    if (i==0)
    {
      edge1[0] = transformed(0,sorted[0]) - startpoint[0];
      edge1[1] = transformed(1,sorted[0]) - startpoint[1];
      edge2[0] = transformed(0,sorted[1]) - transformed(0,sorted[0]);
      edge2[1] = transformed(1,sorted[1]) - transformed(1,sorted[0]);
    }

    // standard case
    else if (i<np-2)
    {
      // go back and find first non-removed partner
      bool foundpartner = false;
      int k=i-1;

      while (foundpartner==false)
      {
        // found non-removed partner
        if (haveremovedthis[k]==false)
        {
          edge1[0] = transformed(0,sorted[i]) - transformed(0,sorted[k]);
          edge1[1] = transformed(1,sorted[i]) - transformed(1,sorted[k]);
          edge2[0] = transformed(0,sorted[i+1]) - transformed(0,sorted[i]);
          edge2[1] = transformed(1,sorted[i+1]) - transformed(1,sorted[i]);
          foundpartner=true;
        }
        else
        {
          // decrease counter
          k-=1;

          // use starting point if all in between removed
          if (k<0)
          {
            edge1[0] = transformed(0,sorted[i]) - startpoint[0];
            edge1[1] = transformed(1,sorted[i]) - startpoint[1];
            edge2[0] = transformed(0,sorted[i+1]) - transformed(0,sorted[i]);
            edge2[1] = transformed(1,sorted[i+1]) - transformed(1,sorted[i]);
            foundpartner=true;
          }
        }
      }
    }

    // last triple
    else /* if i = np-1 */
    {
      // go back and find first non-removed partner
      bool foundpartner = false;
      int k=i-1;

      while (foundpartner==false)
      {
        // found non-removed partner
        if (haveremovedthis[k]==false)
        {
          edge1[0] = transformed(0,sorted[i]) - transformed(0,sorted[k]);
          edge1[1] = transformed(1,sorted[i]) - transformed(1,sorted[k]);
          edge2[0] = startpoint[0] - transformed(0,sorted[i]);
          edge2[1] = startpoint[1] - transformed(1,sorted[i]);
          foundpartner=true;
        }
        else
        {
          // decrease counter
          k-=1;

          // use starting point if all in between removed
          if (k<0)
          {
            edge1[0] = transformed(0,sorted[i]) - startpoint[0];
            edge1[1] = transformed(1,sorted[i]) - startpoint[1];
            edge2[0] = startpoint[0] - transformed(0,sorted[i]);
            edge2[1] = startpoint[1] - transformed(1,sorted[i]);
            foundpartner=true;
          }
        }
      }
    }

    // check for clockwise rotation
    double cw = edge1[0]*edge2[1]-edge1[1]*edge2[0];

    // add point to convex hull if clockwise triple
    // (use tolerance to remove almost straight lines of 3 points)
    if (cw <= -tol)
    {
      Vertex* current = &collconvexhull[sorted[i]];
      respoly.push_back(Vertex(current->Coord(),current->VType(),current->Nodeids(),NULL,NULL,false,false,NULL,-1.0));
    }
    // mark vertex as "removed" if counter-clockwise triple
    else
    {
      removed++;
      haveremovedthis[i]=true;
    }
  }

  return removed;
}

