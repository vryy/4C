#ifndef FLUID_UTILS_MAPEXTRACTOR_H
#define FLUID_UTILS_MAPEXTRACTOR_H

#include "../linalg/linalg_mapextractor.H"
#include "../drt_xfem/dof_management.H"

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
      cond_other                    = 0,
      cond_fsi                      = 1,
      cond_fs                       = 2,
      cond_lung_asi                 = 3,
      cond_vol_surf_flow            = 4,
      cond_total_traction_corrector = 5,
      cond_bio_gr                   = 6,
      cond_ll                       = 7,
      cond_fpsi                     = 8
    };

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis, bool withpressure = false, bool overlapping = false);

    /// get all element gids those nodes are touched by any condition
    Teuchos::RCP<std::set<int> > ConditionedElementMap(const DRT::Discretization& dis) const;

    MAP_EXTRACTOR_VECTOR_METHODS(Other,  cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(FSICond,cond_fsi)
    MAP_EXTRACTOR_VECTOR_METHODS(FSCond, cond_fs)
    MAP_EXTRACTOR_VECTOR_METHODS(LungASICond,cond_lung_asi)
    MAP_EXTRACTOR_VECTOR_METHODS(VolumetricSurfaceFlowCond,cond_vol_surf_flow)
    MAP_EXTRACTOR_VECTOR_METHODS(TotalTractionCorrectionCond,cond_total_traction_corrector)
    MAP_EXTRACTOR_VECTOR_METHODS(BioGrCond,cond_bio_gr)
    MAP_EXTRACTOR_VECTOR_METHODS(LLCond, cond_ll)
    MAP_EXTRACTOR_VECTOR_METHODS(FPSICond,cond_fpsi)
  };

  /// specific MultiMapExtractor to handle the part of fluid with Krylov space projection
  class KSPMapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_other     = 0,
      cond_ksp       = 1
    };

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis);

    /// get all element gids those nodes are touched by any condition
    Teuchos::RCP<std::set<int> > ConditionedElementMap(const DRT::Discretization& dis) const;

    MAP_EXTRACTOR_VECTOR_METHODS(Other,  cond_other)
    MAP_EXTRACTOR_VECTOR_METHODS(KSPCond,cond_ksp)
  };

  /// specific MultiMapExtractor to handle the velocity-pressure split
  class VelPressExtractor : public LINALG::MultiMapExtractor
  {
  public:

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis);

    MAP_EXTRACTOR_VECTOR_METHODS(Velocity, 0)
    MAP_EXTRACTOR_VECTOR_METHODS(Pressure, 1)
  };

  /// specific MultiMapExtractor to handle the fluid field
  class FluidXFluidMapExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_fluid     = 0,
      cond_xfluid    = 1,
    };

    /// setup the whole thing
    void Setup(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> fluidmap, Teuchos::RCP<const Epetra_Map> xfluidmap);

    MAP_EXTRACTOR_VECTOR_METHODS(Fluid, cond_fluid)
    MAP_EXTRACTOR_VECTOR_METHODS(XFluid, cond_xfluid)
  };

}
}

#endif
