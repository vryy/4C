/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 0

*/
/*----------------------------------------------------------------------*/

#include "4C_linalg_mapextractor.hpp"

#include "4C_linalg_utils_sparse_algebra_create.hpp"

#include <Teuchos_getConst.hpp>

#include <cmath>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MultiMapExtractor::MultiMapExtractor() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MultiMapExtractor::MultiMapExtractor(
    const Epetra_Map& fullmap, const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  setup(fullmap, maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::setup(
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
void Core::LinAlg::MultiMapExtractor::check_for_valid_map_extractor() const
{
  if (maps_.size() == 0)
  {
    FOUR_C_THROW("no maps_ available");
  }

  for (unsigned i = 0; i < maps_.size(); ++i)
  {
    if (maps_[i] != Teuchos::null)
    {
      if (maps_[i]->DataPtr() == nullptr)
      {
        FOUR_C_THROW("Got zero data pointer on setup of block %d of maps_\n", i);
      }
      if (not maps_[i]->UniqueGIDs())
      {
        FOUR_C_THROW("map %d not unique", i);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> Core::LinAlg::MultiMapExtractor::MergeMaps(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  if (maps.size() == 0) FOUR_C_THROW("no maps to merge");
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == Teuchos::null) FOUR_C_THROW("can not merge extractor with null maps");
    if (not maps[i]->UniqueGIDs()) FOUR_C_THROW("map %d not unique", i);
  }
  std::set<int> mapentries;
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    const Epetra_Map& map = *maps[i];
    std::copy(map.MyGlobalElements(), map.MyGlobalElements() + map.NumMyElements(),
        std::inserter(mapentries, mapentries.begin()));
  }
  return Core::LinAlg::CreateMap(mapentries, maps[0]->Comm());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> Core::LinAlg::MultiMapExtractor::MergeMapsKeepOrder(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  if (maps.empty()) FOUR_C_THROW("no maps to merge");

  // sanity checks
  for (std::size_t i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == Teuchos::null) FOUR_C_THROW("can not merge extractor with null maps");
    if (not maps[i]->UniqueGIDs()) FOUR_C_THROW("map %d not unique", i);
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
      Teuchos::rcp(new Epetra_Map(-1, gids.size(), gids.data(), 0, maps[0]->Comm()));
  return fullmap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> Core::LinAlg::MultiMapExtractor::IntersectMaps(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  if (maps.size() == 0) FOUR_C_THROW("no maps to intersect");
  for (unsigned i = 0; i < maps.size(); ++i)
  {
    if (maps[i] == Teuchos::null) FOUR_C_THROW("can not intersect extractor with null maps");
    if (not maps[i]->UniqueGIDs()) FOUR_C_THROW("map %d not unique", i);
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
  return Core::LinAlg::CreateMap(mapentries, maps[0]->Comm());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Core::LinAlg::MultiMapExtractor::ExtractVector(
    const Epetra_Vector& full, int block) const
{
  if (maps_[block] == Teuchos::null) FOUR_C_THROW("null map at block %d", block);
  Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp(new Epetra_Vector(*maps_[block]));
  ExtractVector(full, block, *vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> Core::LinAlg::MultiMapExtractor::ExtractVector(
    const Epetra_MultiVector& full, int block) const
{
  if (maps_[block] == Teuchos::null) FOUR_C_THROW("null map at block %d", block);
  Teuchos::RCP<Epetra_MultiVector> vec =
      Teuchos::rcp(new Epetra_MultiVector(*maps_[block], full.NumVectors()));
  ExtractVector(full, block, *vec);
  return vec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::ExtractVector(
    const Epetra_MultiVector& full, int block, Epetra_MultiVector& partial) const
{
  if (maps_[block] == Teuchos::null) FOUR_C_THROW("null map at block %d", block);
  int err = partial.Import(full, *importer_[block], Insert);
  if (err) FOUR_C_THROW("Import using importer returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Core::LinAlg::MultiMapExtractor::InsertVector(
    const Epetra_Vector& partial, int block) const
{
  Teuchos::RCP<Epetra_Vector> full = Teuchos::rcp(new Epetra_Vector(*fullmap_));
  InsertVector(partial, block, *full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> Core::LinAlg::MultiMapExtractor::InsertVector(
    const Epetra_MultiVector& partial, int block) const
{
  Teuchos::RCP<Epetra_MultiVector> full =
      Teuchos::rcp(new Epetra_MultiVector(*fullmap_, partial.NumVectors()));
  InsertVector(partial, block, *full);
  return full;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::InsertVector(
    const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full) const
{
  if (maps_[block] == Teuchos::null) FOUR_C_THROW("null map at block %d", block);
  int err = full.Export(partial, *importer_[block], Insert);
  if (err) FOUR_C_THROW("Export using importer returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::AddVector(
    const Epetra_MultiVector& partial, int block, Epetra_MultiVector& full, double scale) const
{
  Teuchos::RCP<Epetra_MultiVector> v = ExtractVector(full, block);
  if (not v->Map().SameAs(partial.Map())) FOUR_C_THROW("The maps of the vectors must be the same!");
  v->Update(scale, partial, 1.0);
  InsertVector(*v, block, full);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::PutScalar(Epetra_Vector& full, int block, double scalar) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *FullMap();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) FOUR_C_THROW("maps do not match");
    full[lid] = scalar;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Core::LinAlg::MultiMapExtractor::Norm2(const Epetra_Vector& full, int block) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *FullMap();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  double local_norm = 0;

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) FOUR_C_THROW("maps do not match");
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
void Core::LinAlg::MultiMapExtractor::Scale(Epetra_Vector& full, int block, double scalar) const
{
  const Epetra_Map& bm = *Map(block);
  const Epetra_Map& fm = *FullMap();

  int numv = bm.NumMyElements();
  int* v = bm.MyGlobalElements();

  for (int i = 0; i < numv; ++i)
  {
    int lid = fm.LID(v[i]);
    if (lid == -1) FOUR_C_THROW("maps do not match");
    full[lid] *= scalar;
  }
}


/*----------------------------------------------------------------------*
 | Scale one block only                                      fang 08/16 |
 *----------------------------------------------------------------------*/
void Core::LinAlg::MultiMapExtractor::Scale(
    Epetra_MultiVector& full, int block, double scalar) const
{
  for (int i = 0; i < full.NumVectors(); ++i) Scale(*full(i), block, scalar);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MapExtractor::MapExtractor() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MapExtractor::MapExtractor(const Epetra_Map& fullmap,
    Teuchos::RCP<const Epetra_Map> condmap, Teuchos::RCP<const Epetra_Map> othermap)
{
  setup(fullmap, condmap, othermap);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::LinAlg::MapExtractor::MapExtractor(
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
  Teuchos::RCP<Epetra_Map> othermap = Core::LinAlg::CreateMap(othergids, fullmap.Comm());

  // create the extractor based on choice 'iscondmap'
  if (iscondmap)
    setup(fullmap, partialmap, othermap);
  else
    setup(fullmap, othermap, partialmap);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MapExtractor::setup(const Epetra_Map& fullmap,
    const Teuchos::RCP<const Epetra_Map>& condmap, const Teuchos::RCP<const Epetra_Map>& othermap)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(othermap);
  maps.push_back(condmap);
  MultiMapExtractor::setup(fullmap, maps);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::MapExtractor::setup(
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
  Teuchos::RCP<Epetra_Map> othermap = Core::LinAlg::CreateMap(othergids, fullmap.Comm());

  // create the extractor based on choice 'iscondmap'
  if (iscondmap)
    setup(fullmap, partialmap, othermap);
  else
    setup(fullmap, othermap, partialmap);
}

FOUR_C_NAMESPACE_CLOSE
