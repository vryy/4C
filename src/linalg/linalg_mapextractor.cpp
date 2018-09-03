/*----------------------------------------------------------------------*/
/*!
\file linalg_mapextractor.cpp

\brief Implementation

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "linalg_mapextractor.H"
#include "linalg_utils.H"

#include <Teuchos_getConst.hpp>

#include <numeric>
#include <cmath>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MultiMapExtractor::MultiMapExtractor() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MultiMapExtractor::MultiMapExtractor(
    const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Setup(fullmap, maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::Setup(
    const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  fullmap_ = Teuchos::rcp(new Epetra_Map(fullmap));
  maps_ = maps;

  importer_.resize(maps_.size());
  for (unsigned i = 0; i < importer_.size(); ++i)
  {
    if (maps_[i] != Teuchos::null)
    {
      importer_[i] = Teuchos::rcp(new Epetra_Import(*maps_[i], *fullmap_));
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::CheckForValidMapExtractor() const
{
  if (maps_.size() == 0)
  {
    dserror("no maps_ available");
  }

  for (unsigned i = 0; i < maps_.size(); ++i)
  {
    if (maps_[i] != Teuchos::null)
    {
      if (maps_[i]->DataPtr() == NULL)
      {
        dserror("Got zero data pointer on setup of block %d of maps_\n", i);
      }
      if (not maps_[i]->UniqueGIDs())
      {
        dserror("map %d not unique", i);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::MultiMapExtractor::MergeMaps(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  if (maps.size() == 0) dserror("no maps to merge");
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == Teuchos::null) dserror("can not merge extractor with null maps");
    if (not maps[i]->UniqueGIDs()) dserror("map %d not unique", i);
  }
  std::set<int> mapentries;
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    std::copy(map.MyGlobalElements(), map.MyGlobalElements() + map.NumMyElements(),
        std::inserter(mapentries, mapentries.begin()));
  }
  return LINALG::CreateMap(mapentries, maps[0]->Comm());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::MultiMapExtractor::MergeMapsKeepOrder(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  if (maps.empty()) dserror("no maps to merge");

  // sanity checks
  for (std::size_t i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == Teuchos::null) dserror("can not merge extractor with null maps");
    if (not maps[i]->UniqueGIDs()) dserror("map %d not unique", i);
  }

  // collect gids
  std::vector<int> gids;
  for (std::size_t i = 0; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    for (int j = 0; j < map.NumMyElements(); ++j) gids.push_back(map.GID(j));
  }

  // build combined map
  Teuchos::RCP<Epetra_Map> fullmap =
      Teuchos::rcp(new Epetra_Map(-1, gids.size(), &gids[0], 0, maps[0]->Comm()));
  return fullmap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> LINALG::MultiMapExtractor::IntersectMaps(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  if (maps.size() == 0) dserror("no maps to intersect");
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == Teuchos::null) dserror("can not intersect extractor with null maps");
    if (not maps[i]->UniqueGIDs()) dserror("map %d not unique", i);
  }
  std::set<int> mapentries(
      maps[0]->MyGlobalElements(), maps[0]->MyGlobalElements() + maps[0]->NumMyElements());
  for (unsigned i = 1; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    std::set<int> newset;
    int numele = map.NumMyElements();
    int* ele = map.MyGlobalElements();
    for (int j = 0; j < numele; ++j)
    {
      if (mapentries.find(ele[j]) != mapentries.end())
      {
        newset.insert(ele[j]);
      }
    }
    std::swap(mapentries, newset);
  }
  return LINALG::CreateMap(mapentries, maps[0]->Comm());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MultiMapExtractor::ExtractVector(
    const Epetra_Vector& full, int block) const
{
  if (maps_[block] == Teuchos::null) dserror("null map at block %d", block);
  Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*maps_[block]));
  ExtractVector(full, block, *vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::MultiMapExtractor::ExtractVector(
    const Epetra_MultiVector& full, int block) const
{
  if (maps_[block] == Teuchos::null) dserror("null map at block %d", block);
  Teuchos::RCP<Epetra_MultiVector> vec =
      Teuchos::rcp(new Epetra_MultiVector(*maps_[block], full.NumVectors()));
  ExtractVector(full, block, *vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::ExtractVector(
    const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const
{
  if (maps_[block] == Teuchos::null) dserror("null map at block %d", block);
  int err = partial.Import(full, *importer_[block], Insert);
  if (err) dserror("Import using importer returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> LINALG::MultiMapExtractor::InsertVector(
    const Epetra_Vector& partial, int block) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*fullmap_));
  InsertVector(partial, block, *full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::MultiMapExtractor::InsertVector(
    const Epetra_MultiVector& partial, int block) const
{
  Teuchos::RCP<Epetra_MultiVector> full =
      Teuchos::rcp(new Epetra_MultiVector(*fullmap_, partial.NumVectors()));
  InsertVector(partial, block, *full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::InsertVector(
    const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const
{
  if (maps_[block] == Teuchos::null) dserror("null map at block %d", block);
  int err = full.Export(partial, *importer_[block], Insert);
  if (err) dserror("Export using importer returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::AddVector(
    const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full, double scale) const
{
  Teuchos::RCP<Epetra_MultiVector> v = ExtractVector(full, block);
  if (not v->Map().SameAs(partial.Map())) dserror("The maps of the vectors must be the same!");
  v->Update(scale, partial, 1.0);
  InsertVector(*v, block, full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::PutScalar(Epetra_Vector& full, int block, double scalar) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *FullMap();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) dserror("maps do not match");
    full[lid] = scalar;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double LINALG::MultiMapExtractor::Norm2(const Epetra_Vector& full, int block) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *FullMap();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  double local_norm = 0;

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) dserror("maps do not match");
    double value = full[lid];
    local_norm += value * value;
  }

  double global_norm = 0;
  fm.Comm().SumAll(&local_norm, &global_norm, 1);
  return std::sqrt(global_norm);
}


/*----------------------------------------------------------------------*
 | Scale one block only                                      fang 08/16 |
 *----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::Scale(Epetra_Vector& full, int block, double scalar) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *FullMap();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) dserror("maps do not match");
    full[lid] *= scalar;
  }
}


/*----------------------------------------------------------------------*
 | Scale one block only                                      fang 08/16 |
 *----------------------------------------------------------------------*/
void LINALG::MultiMapExtractor::Scale(Epetra_MultiVector& full, int block, double scalar) const
{
  for (int i = 0; i < full.NumVectors(); ++i) Scale(*full(i), block, scalar);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MapExtractor::MapExtractor() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MapExtractor::MapExtractor(const Epetra_Map& fullmap,
    Teuchos::RCP<const Epetra_Map> condmap, Teuchos::RCP<const Epetra_Map> othermap)
{
  Setup(fullmap, condmap, othermap);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::MapExtractor::MapExtractor(
    const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> partialmap, bool iscondmap)
{
  // initialise other DOFs by inserting all DOFs of full map
  std::set<int> othergids;
  const int* fullgids = fullmap.MyGlobalElements();
  copy(fullgids, fullgids + fullmap.NumMyElements(), inserter(othergids, othergids.begin()));

  // throw away all DOFs which are in condmap
  if (partialmap->NumMyElements() > 0)
  {
    const int* condgids = partialmap->MyGlobalElements();
    for (int lid = 0; lid < partialmap->NumMyElements(); ++lid) othergids.erase(condgids[lid]);
  }

  // create (non-overlapping) othermap for non-condmap DOFs
  Teuchos::RCP<Epetra_Map> othermap = LINALG::CreateMap(othergids, fullmap.Comm());

  // create the extractor based on choice 'iscondmap'
  if (iscondmap)
    Setup(fullmap, partialmap, othermap);
  else
    Setup(fullmap, othermap, partialmap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::Setup(const Epetra_Map& fullmap,
    const Teuchos::RCP<const Epetra_Map>& condmap, const Teuchos::RCP<const Epetra_Map>& othermap)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(othermap);
  maps.push_back(condmap);
  MultiMapExtractor::Setup(fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MapExtractor::Setup(
    const Epetra_Map& fullmap, const Teuchos::RCP<const Epetra_Map>& partialmap, bool iscondmap)
{
  // initialise other DOFs by inserting all DOFs of full map
  std::set<int> othergids;
  const int* fullgids = fullmap.MyGlobalElements();
  copy(fullgids, fullgids + fullmap.NumMyElements(), inserter(othergids, othergids.begin()));

  // throw away all DOFs which are in condmap
  if (partialmap->NumMyElements() > 0)
  {
    const int* condgids = partialmap->MyGlobalElements();
    for (int lid = 0; lid < partialmap->NumMyElements(); ++lid) othergids.erase(condgids[lid]);
  }

  // create (non-overlapping) othermap for non-condmap DOFs
  Teuchos::RCP<Epetra_Map> othermap = LINALG::CreateMap(othergids, fullmap.Comm());

  // create the extractor based on choice 'iscondmap'
  if (iscondmap)
    Setup(fullmap, partialmap, othermap);
  else
    Setup(fullmap, othermap, partialmap);
}
