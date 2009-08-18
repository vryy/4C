#ifndef FLUID_UTILS_MAPEXTRACTOR_H
#define FLUID_UTILS_MAPEXTRACTOR_H

#include "../drt_lib/linalg_mapextractor.H"

namespace DRT
{
  class Discretization;
}

namespace FLD
{
namespace UTILS
{

  /// specific MultiMapExtractor to handle the fluid field
  class MapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_other = 0,
      cond_fsi   = 1,
      cond_fs    = 2
    };

    void Setup(const DRT::Discretization& dis);

    MAP_EXTRACTOR_VECTOR_METHODS(Other,  cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(FSICond,cond_fsi)
    MAP_EXTRACTOR_VECTOR_METHODS(FSCond, cond_fs)
  };
}
}

#endif
