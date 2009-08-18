#ifndef ALE_UTILS_MAPEXTRACTOR_H
#define ALE_UTILS_MAPEXTRACTOR_H

#include "../drt_lib/linalg_mapextractor.H"

namespace DRT
{
  class Discretization;
}

namespace ALE
{
namespace UTILS
{

  /// specific MultiMapExtractor to handle the ale field
  class MapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_other = 0,
      cond_fsi   = 1,
      cond_fs    = 2
    };

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis);

    /// get all element gids those nodes are touched by any condition
    Teuchos::RCP<std::set<int> > ConditionedElementMap(const DRT::Discretization& dis) const;

    MAP_EXTRACTOR_VECTOR_METHODS(Other,  cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(FSICond,cond_fsi)
    MAP_EXTRACTOR_VECTOR_METHODS(FSCond, cond_fs)
  };
}
}

#endif
