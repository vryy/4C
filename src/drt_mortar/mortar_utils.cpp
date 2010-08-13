/*!----------------------------------------------------------------------
\file mortar_utils.cpp
\brief A set of utility functions for mortar methods

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "mortar_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"

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
 | transform the row map of a matrix                          popp 08/10|
 *----------------------------------------------------------------------*/
RCP<LINALG::SparseMatrix> MORTAR::MatrixRowTransform(RCP<LINALG::SparseMatrix> inmat,
		                                                 RCP<Epetra_Map> newrowmap)
{
  // initialize output matrix
  RCP<LINALG::SparseMatrix> outmat = rcp(new LINALG::SparseMatrix(*newrowmap,100,false,true));

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

  // complete transformed constraint matrix kzd
  outmat->Complete(inmat->DomainMap(),*newrowmap);

	return outmat;
}

/*----------------------------------------------------------------------*
 | transform the column map of a matrix                       popp 08/10|
 *----------------------------------------------------------------------*/
RCP<LINALG::SparseMatrix> MORTAR::MatrixColTransform(RCP<LINALG::SparseMatrix> inmat,
		                                                 RCP<Epetra_Map> newdomainmap)
{
	//kdz->Complete(*slavemap,*dispmap);

	// initialize output matrix
	RCP<LINALG::SparseMatrix> outmat = rcp(new LINALG::SparseMatrix(inmat->RowMap(),100,false,true));

	// mapping of gids
	map<int,int> gidmap;
	DRT::Exporter ex(inmat->RowMap(),inmat->ColMap(),inmat->Comm());
	for (int i=0; i<inmat->DomainMap().NumMyElements(); ++i) gidmap[inmat->DomainMap().GID(i)] = newdomainmap->GID(i);
	ex.Export(gidmap);

	// transform constraint matrix inmat to lmdofmap
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

	// complete transformed constraint matrix inmat
	outmat->Complete(*newdomainmap,inmat->RowMap());

	return outmat;
}

/*----------------------------------------------------------------------*
 | transform the row and column maps of a matrix              popp 08/10|
 *----------------------------------------------------------------------*/
RCP<LINALG::SparseMatrix> MORTAR::MatrixRowColTransform(RCP<LINALG::SparseMatrix> inmat,
		                                                    RCP<Epetra_Map> newrowmap,
		                                                    RCP<Epetra_Map> newdomainmap)
{
	// initialize output matrix
	RCP<LINALG::SparseMatrix> outmat = rcp(new LINALG::SparseMatrix(*newrowmap,100,false,true));

	// mapping of gids
	map<int,int> gidmap;
	DRT::Exporter ex(inmat->RowMap(),inmat->ColMap(),inmat->Comm());
	for (int i=0; i<inmat->DomainMap().NumMyElements(); ++i) gidmap[inmat->DomainMap().GID(i)] = newdomainmap->GID(i);
	ex.Export(gidmap);

	// transform constraint matrix inmat to lmdofmap
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

	// complete transformed constraint matrix inmat
	outmat->Complete(*newdomainmap,*newrowmap);

	return outmat;
}

#endif  // #ifdef CCADISCRET
