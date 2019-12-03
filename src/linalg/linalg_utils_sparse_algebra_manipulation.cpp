/*----------------------------------------------------------------------*/
/*! \file
 * \file linalg_utils_sparse_algebra_manipulation.cpp

\brief A collection of dense matrix manipulation methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*-----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_lib/drt_dserror.H"
#include <Epetra_IntVector.h>

/*----------------------------------------------------------------------*
 |  export a Epetra_Vector  (public)                         mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::Export(const Epetra_MultiVector& source, Epetra_MultiVector& target)
{
  try
  {
    const bool sourceunique = source.Map().UniqueGIDs();
    const bool targetunique = target.Map().UniqueGIDs();

    // both are unique, does not matter whether ex- or import
    if (sourceunique && targetunique && source.Comm().NumProc() == 1 &&
        target.Comm().NumProc() == 1)
    {
      if (source.NumVectors() != target.NumVectors())
        dserror("number of vectors in source and target not the same!");
      for (int k = 0; k < source.NumVectors(); ++k)
        for (int i = 0; i < target.Map().NumMyElements(); ++i)
        {
          const int gid = target.Map().GID(i);
          if (gid < 0) dserror("No gid for i");
          const int lid = source.Map().LID(gid);
          if (lid < 0) continue;
          // dserror("No source for target");
          (*target(k))[i] = (*source(k))[lid];
        }
      return;
    }
    else if (sourceunique && targetunique)
    {
      Epetra_Export exporter(source.Map(), target.Map());
      int err = target.Export(source, exporter, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
      return;
    }
    else if (sourceunique && !targetunique)
    {
      Epetra_Import importer(target.Map(), source.Map());
      int err = target.Import(source, importer, Insert);
      if (err) dserror("Export using importer returned err=%d", err);
      return;
    }
    else if (!sourceunique && targetunique)
    {
      // copy locally data from source to target
      // do not allow for inter-processor communication to obtain source
      // as this may give a non-unique answer depending on the proc which is asked
      const Epetra_BlockMap& sourcemap = source.Map();
      const Epetra_BlockMap& targetmap = target.Map();
      for (int targetlid = 0; targetlid < targetmap.NumMyElements(); ++targetlid)
      {
        const int sourcelid = sourcemap.LID(targetmap.GID(targetlid));
        if (sourcelid < 0)
          dserror("Export of non-unique source failed. Source data not available on target proc");

        for (int k = 0; k < source.NumVectors(); ++k)
          (*target(k))[targetlid] = (*source(k))[sourcelid];
      }
      return;
    }
    else if (!sourceunique && !targetunique)
    {
      // Neither target nor source are unique - this is a problem.
      // We need a unique in between stage which we have to create artificially.
      // That's nasty.
      // As it is unclear whether this will ever be needed - do it later.
      dserror("Neither target nor source maps are unique - cannot export");
    }
    else
      dserror("VERY strange");
  }
  catch (int error)
  {
    dserror("Caught an Epetra exception %d", error);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  export a Epetra_IntVector  (public)                      mwgee 01/13|
 *----------------------------------------------------------------------*/
void LINALG::Export(const Epetra_IntVector& source, Epetra_IntVector& target)
{
  try
  {
    const bool sourceunique = source.Map().UniqueGIDs();
    const bool targetunique = target.Map().UniqueGIDs();

    // both are unique, does not matter whether ex- or import
    if (sourceunique && targetunique && source.Comm().NumProc() == 1 &&
        target.Comm().NumProc() == 1)
    {
      for (int i = 0; i < target.Map().NumMyElements(); ++i)
      {
        const int gid = target.Map().GID(i);
        if (gid < 0) dserror("No gid for i");
        const int lid = source.Map().LID(gid);
        if (lid < 0) continue;
        target[i] = source[lid];
      }
      return;
    }
    else if (sourceunique && targetunique)
    {
      Epetra_Export exporter(source.Map(), target.Map());
      int err = target.Export(source, exporter, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
      return;
    }
    else if (sourceunique && !targetunique)
    {
      Epetra_Import importer(target.Map(), source.Map());
      int err = target.Import(source, importer, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
      return;
    }
    else if (!sourceunique && targetunique)
    {
      Epetra_Export exporter(source.Map(), target.Map());
      int err = target.Export(source, exporter, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
      return;
    }
    else if (!sourceunique && !targetunique)
    {
      // Neither target nor source are unique - this is a problem.
      // We need a unique in between stage which we have to create artificially.
      // That's nasty.
      // As it is unclear whether this will ever be needed - do it later.
      dserror("Neither target nor source maps are unique - cannot export");
    }
    else
      dserror("VERY strange");
  }
  catch (int error)
  {
    dserror("Caught an Epetra exception %d", error);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::ExtractMyVector(
    const Epetra_Vector& source, const Epetra_Map& target_map)
{
  Teuchos::RCP<Epetra_Vector> target = Teuchos::rcp(new Epetra_Vector(target_map));

  ExtractMyVector(source, *target);

  return target;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LINALG::ExtractMyVector(const Epetra_Vector& source, Epetra_Vector& target)
{
  const int my_num_target_gids = target.Map().NumMyElements();
  const int* my_target_gids = target.Map().MyGlobalElements();

  double* target_values = target.Values();

  const double* src_values = source.Values();

  for (int tar_lid = 0; tar_lid < my_num_target_gids; ++tar_lid)
  {
    const int target_gid = my_target_gids[tar_lid];

    const int src_lid = source.Map().LID(target_gid);
    // check if the target_map is a local sub-set of the source map on each proc
    if (src_lid == -1)
      dserror("Couldn't find the target GID %d in the source map on proc %d.", target_gid,
          source.Comm().MyPID());

    target_values[tar_lid] = src_values[src_lid];
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LINALG::ExtractMyVector(
    double scalar_source, const Epetra_Vector& source, double scalar_target, Epetra_Vector& target)
{
  const int my_num_target_gids = target.Map().NumMyElements();
  const int* my_target_gids = target.Map().MyGlobalElements();

  double* target_values = target.Values();

  const double* src_values = source.Values();

  for (int tar_lid = 0; tar_lid < my_num_target_gids; ++tar_lid)
  {
    const int target_gid = my_target_gids[tar_lid];

    const int src_lid = source.Map().LID(target_gid);
    // check if the target_map is a local sub-set of the source map on each proc
    if (src_lid == -1)
      dserror("Couldn't find the target GID %d in the source map on proc %d.", target_gid,
          source.Comm().MyPID());

    target_values[tar_lid] *= scalar_target;
    target_values[tar_lid] += scalar_source * src_values[src_lid];
  }
}

/*----------------------------------------------------------------------*
 | split matrix into 2x2 block system                              06/06|
 *----------------------------------------------------------------------*/
bool LINALG::SplitMatrix2x2(Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::RCP<Epetra_Map>& A11rowmap,
    Teuchos::RCP<Epetra_Map>& A22rowmap, Teuchos::RCP<Epetra_CrsMatrix>& A11,
    Teuchos::RCP<Epetra_CrsMatrix>& A12, Teuchos::RCP<Epetra_CrsMatrix>& A21,
    Teuchos::RCP<Epetra_CrsMatrix>& A22)
{
  if (A == Teuchos::null) dserror("LINALG::SplitMatrix2x2: A==null on entry");

  if (A11rowmap == Teuchos::null && A22rowmap != Teuchos::null)
    A11rowmap = LINALG::SplitMap(A->RowMap(), *A22rowmap);
  else if (A11rowmap != Teuchos::null && A22rowmap == Teuchos::null)
    A22rowmap = LINALG::SplitMap(A->RowMap(), *A11rowmap);
  else if (A11rowmap == Teuchos::null && A22rowmap == Teuchos::null)
    dserror("LINALG::SplitMatrix2x2: Both A11rowmap and A22rowmap == null on entry");

  std::vector<Teuchos::RCP<const Epetra_Map>> maps(2);
  maps[0] = Teuchos::rcp(new Epetra_Map(*A11rowmap));
  maps[1] = Teuchos::rcp(new Epetra_Map(*A22rowmap));
  LINALG::MultiMapExtractor extractor(A->RowMap(), maps);

  // create SparseMatrix view to input matrix A
  SparseMatrix a(A, View);

  // split matrix into pieces, where main diagonal blocks are square
  Teuchos::RCP<BlockSparseMatrix<DefaultBlockMatrixStrategy>> Ablock =
      a.Split<DefaultBlockMatrixStrategy>(extractor, extractor);
  Ablock->Complete();

  // get Epetra objects out of the block matrix (prevents them from dying)
  A11 = (*Ablock)(0, 0).EpetraMatrix();
  A12 = (*Ablock)(0, 1).EpetraMatrix();
  A21 = (*Ablock)(1, 0).EpetraMatrix();
  A22 = (*Ablock)(1, 1).EpetraMatrix();

  return true;
}

/*----------------------------------------------------------------------*
 | split matrix into 2x2 block system                          gee 02/08|
 *----------------------------------------------------------------------*/
bool LINALG::SplitMatrix2x2(Teuchos::RCP<LINALG::SparseMatrix> A,
    Teuchos::RCP<Epetra_Map>& A11rowmap, Teuchos::RCP<Epetra_Map>& A22rowmap,
    Teuchos::RCP<Epetra_Map>& A11domainmap, Teuchos::RCP<Epetra_Map>& A22domainmap,
    Teuchos::RCP<LINALG::SparseMatrix>& A11, Teuchos::RCP<LINALG::SparseMatrix>& A12,
    Teuchos::RCP<LINALG::SparseMatrix>& A21, Teuchos::RCP<LINALG::SparseMatrix>& A22)
{
  if (A == Teuchos::null) dserror("LINALG::SplitMatrix2x2: A==null on entry");

  // check and complete input row maps
  if (A11rowmap == Teuchos::null && A22rowmap != Teuchos::null)
    A11rowmap = LINALG::SplitMap(A->RowMap(), *A22rowmap);
  else if (A11rowmap != Teuchos::null && A22rowmap == Teuchos::null)
    A22rowmap = LINALG::SplitMap(A->RowMap(), *A11rowmap);
  else if (A11rowmap == Teuchos::null && A22rowmap == Teuchos::null)
    dserror("LINALG::SplitMatrix2x2: Both A11rowmap and A22rowmap == null on entry");

  // check and complete input domain maps
  if (A11domainmap == Teuchos::null && A22domainmap != Teuchos::null)
    A11domainmap = LINALG::SplitMap(A->DomainMap(), *A22domainmap);
  else if (A11domainmap != Teuchos::null && A22domainmap == Teuchos::null)
    A22domainmap = LINALG::SplitMap(A->DomainMap(), *A11domainmap);
  else if (A11rowmap == Teuchos::null && A22rowmap == Teuchos::null)
    dserror("LINALG::SplitMatrix2x2: Both A11domainmap and A22domainmap == null on entry");

  // local variables
  std::vector<Teuchos::RCP<const Epetra_Map>> rangemaps(2);
  std::vector<Teuchos::RCP<const Epetra_Map>> domainmaps(2);
  rangemaps[0] = Teuchos::rcp(new Epetra_Map(*A11rowmap));
  rangemaps[1] = Teuchos::rcp(new Epetra_Map(*A22rowmap));
  domainmaps[0] = Teuchos::rcp(new Epetra_Map(*A11domainmap));
  domainmaps[1] = Teuchos::rcp(new Epetra_Map(*A22domainmap));
  LINALG::MultiMapExtractor range(A->RangeMap(), rangemaps);
  LINALG::MultiMapExtractor domain(A->DomainMap(), domainmaps);

  Teuchos::RCP<BlockSparseMatrix<DefaultBlockMatrixStrategy>> Ablock =
      A->Split<DefaultBlockMatrixStrategy>(domain, range);

#if 0  // debugging
  std::cout << "A00\n" << (*Ablock)(0,0);
  std::cout << "A10\n" << (*Ablock)(1,0);
  std::cout << "A01\n" << (*Ablock)(0,1);
  std::cout << "A11\n" << (*Ablock)(1,1);
  std::cout << "A->Range\n" << A->RangeMap();
  std::cout << "A->Domain\n" << A->DomainMap();
  std::cout << "A11domainmap\n" << *A11domainmap;
  std::cout << "A22domainmap\n" << *A22domainmap;
#endif

  Ablock->Complete();
  // extract internal data from Ablock in Teuchos::RCP form and let Ablock die
  // (this way, internal data from Ablock will live)
  A11 = Teuchos::rcp(new SparseMatrix((*Ablock)(0, 0), View));
  A12 = Teuchos::rcp(new SparseMatrix((*Ablock)(0, 1), View));
  A21 = Teuchos::rcp(new SparseMatrix((*Ablock)(1, 0), View));
  A22 = Teuchos::rcp(new SparseMatrix((*Ablock)(1, 1), View));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int LINALG::InsertMyRowDiagonalIntoUnfilledMatrix(
    LINALG::SparseMatrix& mat, const Epetra_Vector& diag)
{
  if (mat.Filled()) return -1;

  Teuchos::RCP<Epetra_CrsMatrix> dst_mat_ptr = mat.EpetraMatrix();
  Epetra_CrsMatrix& dst_mat = *dst_mat_ptr;

  const int my_num_entries = diag.Map().NumMyElements();
  const int* my_gids = diag.Map().MyGlobalElements();

  double* diag_values = diag.Values();

  for (int lid = 0; lid < my_num_entries; ++lid)
  {
    const int rgid = my_gids[lid];

    // skip rows which are not part of the matrix
    if (not dst_mat.RangeMap().MyGID(rgid))
      dserror(
          "Could not find the row GID %d in the destination matrix RowMap"
          " on proc %d.",
          rgid, dst_mat.Comm().MyPID());

    if (dst_mat.NumAllocatedGlobalEntries(rgid))
    {
      // add all values, including zeros, as we need a proper matrix graph
      int err = dst_mat.SumIntoGlobalValues(rgid, 1, (diag_values + lid), &rgid);
      if (err > 0)
      {
        err = dst_mat.InsertGlobalValues(rgid, 1, (diag_values + lid), &rgid);
        if (err < 0) dserror("InsertGlobalValues error: %d", err);
      }
      else if (err < 0)
        dserror("SumIntoGlobalValues error: %d", err);
    }
    else
    {
      const int err = dst_mat.InsertGlobalValues(rgid, 1, (diag_values + lid), &rgid);
      if (err < 0) dserror("InsertGlobalValues error: %d", err);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 | split a map into 2 pieces with given Agiven                     06/06|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::SplitMap(const Epetra_Map& Amap, const Epetra_Map& Agiven)
{
  const Epetra_Comm& Comm = Amap.Comm();
  const Epetra_Map& Ag = Agiven;

  int count = 0;
  std::vector<int> myaugids(Amap.NumMyElements());
  for (int i = 0; i < Amap.NumMyElements(); ++i)
  {
    const int gid = Amap.GID(i);
    if (Ag.MyGID(gid)) continue;
    myaugids[count] = gid;
    ++count;
  }
  myaugids.resize(count);
  int gcount;
  Comm.SumAll(&count, &gcount, 1);
  Teuchos::RCP<Epetra_Map> Aunknown =
      Teuchos::rcp(new Epetra_Map(gcount, count, &myaugids[0], 0, Comm));

  return Aunknown;
}

/*----------------------------------------------------------------------*
 | fill matrix row and check for success                     farah 06/14|
 *----------------------------------------------------------------------*/
void LINALG::InsertGlobalValues(
    Teuchos::RCP<Epetra_CrsMatrix> mat, int GlobalRow, int NumEntries, double* Values, int* Indices)
{
  int err;

  if (NumEntries > 0)
    err = mat->InsertGlobalValues(GlobalRow, NumEntries, Values, Indices);
  else
    err = mat->InsertGlobalValues(GlobalRow, NumEntries, Values, 0);

  if (err) dserror("InsertGlobalValues err=%d", err);

  return;
}

/*----------------------------------------------------------------------*
 | merge two given maps to one map                            popp 01/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::MergeMap(
    const Epetra_Map& map1, const Epetra_Map& map2, bool overlap)
{
  // check for unique GIDs and for identity
  // if ((!map1.UniqueGIDs()) || (!map2.UniqueGIDs()))
  //  dserror("LINALG::MergeMap: One or both input maps are not unique");
  if (map1.SameAs(map2))
  {
    if ((overlap == false) && map1.NumGlobalElements() > 0)
      dserror("LINALG::MergeMap: Result map is overlapping");
    else
      return Teuchos::rcp(new Epetra_Map(map1));
  }

  std::vector<int> mygids(map1.NumMyElements() + map2.NumMyElements());
  int count = map1.NumMyElements();

  // get GIDs of input map1
  for (int i = 0; i < count; ++i) mygids[i] = map1.GID(i);

  // add GIDs of input map2 (only new ones)
  for (int i = 0; i < map2.NumMyElements(); ++i)
  {
    // check for overlap
    if (map1.MyGID(map2.GID(i)))
    {
      if (overlap == false) dserror("LINALG::MergeMap: Result map is overlapping");
    }
    // add new GIDs to mygids
    else
    {
      mygids[count] = map2.GID(i);
      ++count;
    }
  }
  mygids.resize(count);

  // sort merged map
  sort(mygids.begin(), mygids.end());

  return Teuchos::rcp(new Epetra_Map(-1, (int)mygids.size(), &mygids[0], 0, map1.Comm()));
}

/*----------------------------------------------------------------------*
 | merge two given maps to one map                            popp 01/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::MergeMap(const Teuchos::RCP<const Epetra_Map>& map1,
    const Teuchos::RCP<const Epetra_Map>& map2, bool overlap)
{
  // check for cases with null Teuchos::RCPs
  if (map1 == Teuchos::null && map2 == Teuchos::null)
    return Teuchos::null;
  else if (map1 == Teuchos::null)
    return Teuchos::rcp(new Epetra_Map(*map2));
  else if (map2 == Teuchos::null)
    return Teuchos::rcp(new Epetra_Map(*map1));

  // wrapped call to non-Teuchos::RCP version of MergeMap
  return LINALG::MergeMap(*map1, *map2, overlap);
}

/*----------------------------------------------------------------------*
 | Find the intersection of two maps                     hiermeier 10/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::IntersectMap(const Epetra_Map& map1, const Epetra_Map& map2)
{
  // check if the maps are identical
  if (map1.SameAs(map2))
  {
    return Teuchos::rcp(new Epetra_Map(map1));
  }

  std::vector<int> mygids(std::min(map1.NumMyElements(), map2.NumMyElements()), -1);
  int count = 0;

  for (int i = 0; i < map1.NumMyElements(); ++i)
  {
    // check for intersecting gids
    if (map2.MyGID(map1.GID(i)))
    {
      mygids[count] = map1.GID(i);
      ++count;
    }
  }
  mygids.resize(count);

  // sort merged map
  sort(mygids.begin(), mygids.end());

  return Teuchos::rcp(new Epetra_Map(-1, (int)mygids.size(), &mygids[0], 0, map1.Comm()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::ExtractMyOverlappingSubMap(
    const Epetra_BlockMap& src_map, const Epetra_BlockMap& ref_map)
{
  std::vector<int> my_overlapping_gids;
  my_overlapping_gids.reserve(src_map.NumMyElements());

  const int num_my_src_entries = src_map.NumMyElements();
  const int* my_src_gids = src_map.MyGlobalElements();

  for (int i = 0; i < num_my_src_entries; ++i)
  {
    const int my_src_gid = my_src_gids[i];
    if (ref_map.MyGID(my_src_gid)) my_overlapping_gids.push_back(my_src_gid);
  }

  return Teuchos::rcp(new Epetra_Map(
      -1, my_overlapping_gids.size(), my_overlapping_gids.data(), 0, src_map.Comm()));
}



/*----------------------------------------------------------------------*
 | split a vector into 2 pieces with given submaps            popp 02/08|
 *----------------------------------------------------------------------*/
bool LINALG::SplitVector(const Epetra_Map& xmap, const Epetra_Vector& x,
    Teuchos::RCP<Epetra_Map>& x1map, Teuchos::RCP<Epetra_Vector>& x1,
    Teuchos::RCP<Epetra_Map>& x2map, Teuchos::RCP<Epetra_Vector>& x2)
{
  // map extractor with fullmap(xmap) and two other maps (x1map and x2map)
  LINALG::MapExtractor extractor(xmap, x1map, x2map);

  // extract subvectors from fullvector
  x1 = extractor.ExtractVector(x, 1);
  x2 = extractor.ExtractVector(x, 0);

  return true;
}

/*----------------------------------------------------------------------*
 | split a vector into 2 pieces with given submaps           farah 02/16|
 *----------------------------------------------------------------------*/
bool LINALG::SplitVector(const Epetra_Map& xmap, const Epetra_Vector& x,
    Teuchos::RCP<const Epetra_Map>& x1map, Teuchos::RCP<Epetra_Vector>& x1,
    Teuchos::RCP<const Epetra_Map>& x2map, Teuchos::RCP<Epetra_Vector>& x2)
{
  // map extractor with fullmap(xmap) and two other maps (x1map and x2map)
  LINALG::MapExtractor extractor(xmap, x1map, x2map);

  // extract subvectors from fullvector
  x1 = extractor.ExtractVector(x, 1);
  x2 = extractor.ExtractVector(x, 0);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ComputeDofMapsFromNodeMaps(const int dofset_id,
    const std::vector<Teuchos::RCP<Epetra_Map>>& node_maps,
    const DRT::DiscretizationInterface& discret, std::vector<Teuchos::RCP<Epetra_Map>>& dof_maps)
{
  dof_maps.reserve(dof_maps.size() + node_maps.size());
  for (const Teuchos::RCP<Epetra_Map>& node_map : node_maps)
    dof_maps.push_back(ComputeDofMapFromNodeMap(dofset_id, *node_map, discret));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::ComputeDofMapFromNodeMap(
    const int dofset_id, const Epetra_Map& node_map, const DRT::DiscretizationInterface& discret)
{
  std::set<int> dof_set;
  std::vector<int> dof_vec;

  const int my_num_nodes = node_map.NumMyElements();
  const int* my_ngids = node_map.MyGlobalElements();

  for (int nlid = 0; nlid < my_num_nodes; ++nlid)
  {
    const DRT::Node* node = discret.gNode(my_ngids[nlid]);

    const int numdofs = discret.NumDof(dofset_id, node);
    for (int d = 0; d < numdofs; ++d) dof_set.insert(discret.Dof(dofset_id, node, d));
  }

  dof_vec.resize(dof_set.size());
  std::copy(dof_set.begin(), dof_set.end(), dof_vec.begin());

  return Teuchos::rcp(new Epetra_Map(-1, dof_vec.size(), dof_vec.data(), 0, discret.Comm()));
}
