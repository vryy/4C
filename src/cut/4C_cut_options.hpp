/*---------------------------------------------------------------------*/
/*! \file

\brief options to set up the cut

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_OPTIONS_HPP
#define FOUR_C_CUT_OPTIONS_HPP

#include "4C_config.hpp"

#include "4C_cut_enum.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    /*!
     \brief Options defined to configure the cutting behavior
     */
    class Options
    {
     public:
      Options(NodalDofSetStrategy nodal_dofset_strategy = NDS_Strategy_full, bool positions = true,
          bool simpleshapes = true, bool genline2 = true, bool genquad4 = true, bool genhex8 = true,
          bool genwedge6 = false, bool genpyramid5 = false)
          : geomintersect_floattype_(floattype_double),
            geomdistance_floattype_(floattype_double),
            general_position_dist_floattype_(floattype_none),
            general_position_pos_floattype_(floattype_none),
            direct_divergence_refplane_(DirDiv_refplane_none),
            nodal_dofset_strategy_(nodal_dofset_strategy),
            positions_(positions),
            simpleshapes_(simpleshapes),
            genquad4_(genquad4),
            genhex8_(genhex8),
            genwedge6_(genwedge6),
            genpyramid5_(genpyramid5),
            do_selfcut_(true),
            selfcut_do_meshcorrection_(true),
            selfcut_island_geom_multiplicator_(2),
            gen_bcell_position_(bcells_on_cut_side),
            split_cutsides_(true),
            bc_cubaturedegree_(20)
      {
      }

      /// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT
      /// GENERAL)
      void Init_by_Paramlist();

      /// Initializes Cut Parameters by Parameterlist (typically from *.dat-file section CUT
      /// GENERAL)
      void Init_by_Paramlist(const Teuchos::ParameterList& cutparams);

      /// Initializes Cut Parameters for Cuttests (use full cln)
      void Init_for_Cuttests();

      void SetFindPositions(bool positions) { positions_ = positions; }

      void set_nodal_dof_set_strategy(Core::Geo::Cut::NodalDofSetStrategy nodal_dofset_strategy)
      {
        nodal_dofset_strategy_ = nodal_dofset_strategy;
      }

      Core::Geo::Cut::NodalDofSetStrategy get_nodal_dof_set_strategy()
      {
        return nodal_dofset_strategy_;
      }

      bool FindPositions() const { return positions_; }

      bool SimpleShapes() const { return simpleshapes_; }

      void SetSimpleShapes(bool simpleshapes) { simpleshapes_ = simpleshapes; }

      /** \brief Set the position for the boundary cell creation
       *
       *  \author hiermeier \date 01/17 */
      void set_gen_boundary_cell_position(BoundaryCellPosition gen_bcell_position)
      {
        gen_bcell_position_ = gen_bcell_position;
      };

      /** \brief Float_type for geometric intersection computation */
      enum Core::Geo::Cut::CutFloatType geom_intersect_floattype() const
      {
        return geomintersect_floattype_;
      }

      /** \brief Float_type for geometric distance computation */
      enum Core::Geo::Cut::CutFloatType geom_distance_floattype() const
      {
        return geomdistance_floattype_;
      }

      /** \brief Which Referenceplanes are used in DirectDivergence */
      enum CutDirectDivergenceRefplane direct_divergence_refplane() const
      {
        return direct_divergence_refplane_;
      }

      enum BoundaryCellPosition gen_boundary_cell_position() const { return gen_bcell_position_; }

      /** \brief get the quad4 integration cell generation indicator
       *  (TRUE: QUAD4, FALSE: TRI3's are used) */
      bool GenQuad4() const { return genquad4_; }

      /** \brief get the hex8 integration cell generation indicator
       *  (TRUE: HEX8, FALSE: TRI4's are used) */
      bool GenHex8() const { return genhex8_; }

      /** \brief get the wedge6 integration cell generation indicator
       *  (TRUE: WEDGE6, FALSE: TRI4's are used) */
      bool GenWedge6() const { return genwedge6_; }

      /** \brief get the pyramid5 integration cell generation indicator
       *  (TRUE: PYRAMID5, FALSE: TRI4's are used) */
      bool GenPyramid5() const { return genpyramid5_; }

      /** \brief Indicator if quad4 cutsides are split into tri3 cutsides */
      bool SplitCutSides() const { return split_cutsides_; }

      /** \brief Cubaturedegree for creating of integrationpoints on boundarycells */
      int BC_Cubaturedegree() const { return bc_cubaturedegree_; }

      /** \brief Run the SelfCut Algorithm*/
      bool Do_SelfCut() { return do_selfcut_; }

      /** \brief perform the mesh correction in the selfcut algorithm*/
      bool self_cut_do_mesh_correction() { return selfcut_do_meshcorrection_; }

      /** \brief multiplicator of the solid element size to specify the maximal size of an island in
       * the selfcut*/
      int self_cut_island_geom_multiplicator() { return selfcut_island_geom_multiplicator_; }

     private:
      /** \brief Float_type for geometric intersection computation */
      Core::Geo::Cut::CutFloatType geomintersect_floattype_;

      /** \brief Float_type for geometric distance computation */
      Core::Geo::Cut::CutFloatType geomdistance_floattype_;

      /** \brief Float_type used in Core::Geo::Cut::POSITION for compute_distance*/
      Core::Geo::Cut::CutFloatType general_position_dist_floattype_;

      /** \brief Float_type used in Core::Geo::Cut::POSITION for ComputePosition*/
      Core::Geo::Cut::CutFloatType general_position_pos_floattype_;

      /** \brief Specifies which Referenceplanes are used in DirectDivergence*/
      CutDirectDivergenceRefplane direct_divergence_refplane_;

      Core::Geo::Cut::NodalDofSetStrategy nodal_dofset_strategy_;

      bool positions_;

      bool simpleshapes_;

      /** \brief Indicator if QUAD4 integrations cells are going to be used or not
       *
       *  If TRUE a QUAD4 element is used as integration cell or otherwise the
       *  QUAD4 is splitted into 2 TRI3 elements. */
      bool genquad4_;

      /** \brief Indicator if HEX8 integrations cells are going to be used or not
       *
       *  If TRUE a HEX8 element is used as integration cell or otherwise the
       *  HEX8 is splitted into 5 TRI4 elements. */
      bool genhex8_;

      /** \brief Indicator if WEDGE6 integrations cells are going to be used or not
       *
       *  If TRUE a WEDGE6 element is used as integration cell or otherwise the
       *  WEDGE6 is splitted into 3 TRI4 elements. */
      bool genwedge6_;

      /** \brief Indicator if PYRAMID5 integrations cells are going to be used or not
       *
       *  If TRUE a PYRAMID5 element is used as integration cell or otherwise the
       *  PYRAMID5 is splitted into 2 TRI4 elements. */
      bool genpyramid5_;

      /** \brief Run the SelfCut Algorithm*/
      bool do_selfcut_;

      /** \brief perform the mesh correction in the selfcut algorithm*/
      bool selfcut_do_meshcorrection_;

      /** \brief multiplicator of the solid element size to specify the maximal size of an island in
       * the selfcut*/
      int selfcut_island_geom_multiplicator_;


      /** \brief Where to create boundary cells */
      enum BoundaryCellPosition gen_bcell_position_;

      /** \brief Indicator if quad4 cutsides are split into tri3 cutsides */
      bool split_cutsides_;

      /** \brief Cubaturedegree for creating of integrationpoints on boundarycells */
      int bc_cubaturedegree_;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
