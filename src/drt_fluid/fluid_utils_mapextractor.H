/*-----------------------------------------------------------*/
/*! \file

\brief extracting maps of fluid discretizations

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_UTILS_MAPEXTRACTOR_H
#define FLUID_UTILS_MAPEXTRACTOR_H


#include "../linalg/linalg_mapextractor.H"


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
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_fs = 2,
        cond_lung_asi = 3,
        cond_mortar = 4,
        cond_au = 5
      };

      /// setup the whole thing
      void Setup(const DRT::Discretization& dis, bool withpressure = false,
          bool overlapping = false, const int nds_master = 0);

      /*!
       * \brief setup from an existing extractor
       * By calling this setup version we create a map extractor from
       * (1) an existing map extractor and
       * (2) a DOF-map from another discretization, which is appended to othermap.
       * We need this in the context of XFFSI.
       * \param (in) additionalothermap : map of additional unconditioned DOF
       * \param (in) extractor : extractor, from which the conditions are cloned
       * \author kruse
       * \date 05/2014
       */
      void Setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
          const FLD::UTILS::MapExtractor& extractor);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> ConditionedElementMap(const DRT::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(FSICond, cond_fsi)
      MAP_EXTRACTOR_VECTOR_METHODS(FSCond, cond_fs)
      MAP_EXTRACTOR_VECTOR_METHODS(LungASICond, cond_lung_asi)
      MAP_EXTRACTOR_VECTOR_METHODS(MortarCond, cond_mortar)
      MAP_EXTRACTOR_VECTOR_METHODS(AUCond, cond_au)
    };

    /// specific MultiMapExtractor to handle the part of fluid with volumetric surface flow
    /// condition
    class VolumetricFlowMapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_vol_surf_flow = 1
      };

      /// setup the whole thing
      void Setup(const DRT::Discretization& dis);

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(VolumetricSurfaceFlowCond, cond_vol_surf_flow)
    };

    /// specific MultiMapExtractor to handle the part of fluid with Krylov space projection
    class KSPMapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_ksp = 1
      };

      /// setup the whole thing
      void Setup(const DRT::Discretization& dis);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> ConditionedElementMap(const DRT::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(KSPCond, cond_ksp)
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

    /// specific MultiMapExtractor to handle the fsi and ale meshtying at the same time
    class FsiMapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1
      };

      /// setup the whole thing
      void Setup(const DRT::Discretization& dis);

      void Setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
          const FLD::UTILS::FsiMapExtractor& extractor);

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(FSI, cond_fsi)
    };

    /// specific MultiMapExtractor to handle the fluid field
    class XFluidFluidMapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_fluid = 0,
        cond_xfluid = 1,
      };

      /// setup the whole thing
      void Setup(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> fluidmap,
          Teuchos::RCP<const Epetra_Map> xfluidmap);

      MAP_EXTRACTOR_VECTOR_METHODS(Fluid, cond_fluid)
      MAP_EXTRACTOR_VECTOR_METHODS(XFluid, cond_xfluid)
    };

  }  // namespace UTILS
}  // namespace FLD

#endif
