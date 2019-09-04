/*----------------------------------------------------------------------------*/
/*! \file


\brief map extractor functionalities for interface problems for ALE field

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE_UTILS_MAPEXTRACTOR_H
#define ALE_UTILS_MAPEXTRACTOR_H

#include "../linalg/linalg_mapextractor.H"

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
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_fs = 2,
        cond_lung_asi = 3,
        cond_ale_wear = 4,
        cond_bio_gr = 5,
        cond_au = 6,
        cond_fpsi = 7,
        cond_mortar = 8
      };

      /// application-specific types of Dirichlet sets
      enum AleDBCSetType
      {
        dbc_set_std = 0,       ///< type of Dirichlet set in standard ALE-FSI
        dbc_set_x_ff = 1,      ///< Dirichlet sets include fluid-fluid-interface DOF (during Newton
                               ///< iteration of XFFSI)
        dbc_set_x_fsi = 2,     ///< Dirichlet sets include ALE-sided FSI interface DOF (in
                               ///< ALE-relaxation step of XFFSI)
        dbc_set_biofilm = 3,   ///< type of Dirichlet set for biofilm applications
        dbc_set_part_fsi = 4,  ///< type of Dirichlet set for partitioned FSI
        dbc_set_wear = 5,      ///< type of Dirichlet set for wear application
        dbc_set_count = 6      ///< total number of types
      };

      /// setup the whole thing
      void Setup(const DRT::Discretization& dis, bool overlapping = false);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> ConditionedElementMap(const DRT::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(FSICond, cond_fsi)
      MAP_EXTRACTOR_VECTOR_METHODS(FSCond, cond_fs)
      MAP_EXTRACTOR_VECTOR_METHODS(LungASICond, cond_lung_asi)
      MAP_EXTRACTOR_VECTOR_METHODS(AleWearCond, cond_ale_wear)
      MAP_EXTRACTOR_VECTOR_METHODS(BioGrCond, cond_bio_gr)
      MAP_EXTRACTOR_VECTOR_METHODS(AUCond, cond_au)
      MAP_EXTRACTOR_VECTOR_METHODS(FPSICond, cond_fpsi)
      MAP_EXTRACTOR_VECTOR_METHODS(Mortar, cond_mortar)
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

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(FSI, cond_fsi)
    };

    /// specific MultiMapExtractor to handle the fluid_fluid_Coupling
    class XFluidFluidMapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_xfluidfluid = 1
      };

      /// setup the whole thing
      void Setup(const DRT::Discretization& dis);

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(XFluidFluidCond, cond_xfluidfluid)
    };
  }  // namespace UTILS
}  // namespace ALE

#endif
