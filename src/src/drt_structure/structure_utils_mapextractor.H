#ifndef STRUCTURE_UTILS_MAPEXTRACTOR_H
#define STRUCTURE_UTILS_MAPEXTRACTOR_H

#include "../drt_lib/linalg_mapextractor.H"

namespace DRT
{
  class Discretization;
}

namespace STR
{
namespace UTILS
{

  /// specific MultiMapExtractor to handle the fluid field
  class MapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_other     = 0,
      cond_fsi       = 1,
      cond_lung_asi  = 2
    };

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis,  Epetra_Map fullmap);

    /// get all element gids those nodes are touched by any condition
    Teuchos::RCP<std::set<int> > ConditionedElementMap(const DRT::Discretization& dis) const;

    MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(FSICond,cond_fsi)
    MAP_EXTRACTOR_VECTOR_METHODS(LungASICond,cond_lung_asi)
  };
}
}

#endif
