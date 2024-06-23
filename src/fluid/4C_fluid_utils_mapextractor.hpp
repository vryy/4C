/*-----------------------------------------------------------*/
/*! \file

\brief extracting maps of fluid discretizations


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_UTILS_MAPEXTRACTOR_HPP
#define FOUR_C_FLUID_UTILS_MAPEXTRACTOR_HPP


#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FLD
{
  namespace UTILS
  {
    /// specific MultiMapExtractor to handle the fluid field
    class MapExtractor : public Core::LinAlg::MultiMapExtractor
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
      void setup(const Core::FE::Discretization& dis, bool withpressure = false,
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
      void setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
          const FLD::UTILS::MapExtractor& extractor);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> conditioned_element_map(
          const Core::FE::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(fsi_cond, cond_fsi)
      MAP_EXTRACTOR_VECTOR_METHODS(fs_cond, cond_fs)
      MAP_EXTRACTOR_VECTOR_METHODS(lung_asi_cond, cond_lung_asi)
      MAP_EXTRACTOR_VECTOR_METHODS(mortar_cond, cond_mortar)
      MAP_EXTRACTOR_VECTOR_METHODS(au_cond, cond_au)
    };

    /// specific MultiMapExtractor to handle the part of fluid with volumetric surface flow
    /// condition
    class VolumetricFlowMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_vol_surf_flow = 1
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      MAP_EXTRACTOR_VECTOR_METHODS(other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(volumetric_surface_flow_cond, cond_vol_surf_flow)
    };

    /// specific MultiMapExtractor to handle the part of fluid with Krylov space projection
    class KSPMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_ksp = 1
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> conditioned_element_map(
          const Core::FE::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(ksp_cond, cond_ksp)
    };

    /// specific MultiMapExtractor to handle the velocity-pressure split
    class VelPressExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      MAP_EXTRACTOR_VECTOR_METHODS(velocity, 0)
      MAP_EXTRACTOR_VECTOR_METHODS(pressure, 1)
    };

    /// specific MultiMapExtractor to handle the fsi and ale meshtying at the same time
    class FsiMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1
      };

      /// setup the whole thing
      void setup(const Core::FE::Discretization& dis);

      void setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
          const FLD::UTILS::FsiMapExtractor& extractor);

      MAP_EXTRACTOR_VECTOR_METHODS(other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(fsi, cond_fsi)
    };

    /// specific MultiMapExtractor to handle the fluid field
    class XFluidFluidMapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_fluid = 0,
        cond_xfluid = 1,
      };

      /// setup the whole thing
      void setup(const Epetra_Map& fullmap, Teuchos::RCP<const Epetra_Map> fluidmap,
          Teuchos::RCP<const Epetra_Map> xfluidmap);

      MAP_EXTRACTOR_VECTOR_METHODS(fluid, cond_fluid)
      MAP_EXTRACTOR_VECTOR_METHODS(x_fluid, cond_xfluid)
    };

  }  // namespace UTILS
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
