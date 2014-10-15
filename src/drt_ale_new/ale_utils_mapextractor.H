/*----------------------------------------------------------------------------*/
/*!

\file ale_utils_mapextractor.H

\brief map extractor functionalities for interface problems for ALE field

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALENEW_UTILS_MAPEXTRACTOR_H
#define ALENEW_UTILS_MAPEXTRACTOR_H

#include "../linalg/linalg_mapextractor.H"

namespace DRT
{
  class Discretization;
}

namespace ALENEW
{
namespace UTILS
{

  /// specific MultiMapExtractor to handle the ale field
  class MapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_other     = 0,
      cond_fsi       = 1,
      cond_fs        = 2,
      cond_lung_asi  = 3,
      cond_ale_wear  = 4,
      cond_bio_gr    = 5,
      cond_au        = 6,
      cond_fpsi      = 7
    };

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis, bool overlapping = false);

    /// get all element gids those nodes are touched by any condition
    Teuchos::RCP<std::set<int> > ConditionedElementMap(const DRT::Discretization& dis) const;

    MAP_EXTRACTOR_VECTOR_METHODS(Other,  cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(FSICond,cond_fsi)
    MAP_EXTRACTOR_VECTOR_METHODS(FSCond, cond_fs)
    MAP_EXTRACTOR_VECTOR_METHODS(LungASICond,cond_lung_asi)
    MAP_EXTRACTOR_VECTOR_METHODS(AleWearCond,cond_ale_wear)
    MAP_EXTRACTOR_VECTOR_METHODS(BioGrCond,cond_bio_gr)
    MAP_EXTRACTOR_VECTOR_METHODS(AUCond, cond_au)
    MAP_EXTRACTOR_VECTOR_METHODS(FPSICond,cond_fpsi)
  };

  /// specific MultiMapExtractor to handle the fluid_fluid_Coupling
  class XFluidFluidMapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_other       = 0,
      cond_xfluidfluid = 1
    };

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis);

    MAP_EXTRACTOR_VECTOR_METHODS(Other,cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(XFluidFluidCond,cond_xfluidfluid)
  };
}
}

#endif
