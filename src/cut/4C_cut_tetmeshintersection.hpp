/*---------------------------------------------------------------------*/
/*! \file

\brief Cut tets used for the Tessellation routine

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_TETMESHINTERSECTION_HPP
#define FOUR_C_CUT_TETMESHINTERSECTION_HPP

#include "4C_config.hpp"

#include "4C_cut_mesh.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class Side;
    class Mesh;
    class Options;
    class TetMesh;
    class VolumeCell;

    /// interface class for scary recursive cuts.
    class TetMeshIntersection
    {
     public:
      TetMeshIntersection(Options& options, Element* element,
          const std::vector<std::vector<int>>& tets, const std::vector<int>& accept_tets,
          const std::vector<Point*>& points, const plain_side_set& cut_sides);

      void Cut(Mesh& parent_mesh, Element* element, const plain_volumecell_set& parent_cells,
          int count, bool tetcellsonly = false);

      Mesh& NormalMesh() { return mesh_; }

     private:
      struct ChildCell
      {
        ChildCell() : done_(false), parent_(nullptr) {}

        bool ContainsChild(VolumeCell* vc) { return cells_.count(vc) > 0; }

        bool done_;
        VolumeCell* parent_;
        plain_volumecell_set cells_;
        std::map<Side*, std::vector<Facet*>> facetsonsurface_;
      };

      struct FacetMesh
      {
        void Add(Facet* f)
        {
          const std::vector<Point*>& points = f->Points();
          for (unsigned i = 0; i != points.size(); ++i)
          {
            Point* p1 = points[i];
            Point* p2 = points[(i + 1) % points.size()];

            if (p1 > p2) std::swap(p1, p2);

            facet_mesh_[std::make_pair(p1, p2)].push_back(f);
          }
        }

        void Erase(Facet* f)
        {
          const std::vector<Point*>& points = f->Points();
          for (unsigned i = 0; i != points.size(); ++i)
          {
            Point* p1 = points[i];
            Point* p2 = points[(i + 1) % points.size()];

            if (p1 > p2) std::swap(p1, p2);

            std::vector<Facet*>& facets = facet_mesh_[std::make_pair(p1, p2)];
            std::vector<Facet*>::iterator j = std::find(facets.begin(), facets.end(), f);
            if (j != facets.end()) facets.erase(j);
          }
        }

        std::map<std::pair<Point*, Point*>, std::vector<Facet*>> facet_mesh_;
      };

      void find_edge_cuts();

      /// find the mapping between child VolumeCells and parent VolumeCells.
      void map_volume_cells(Mesh& parent_mesh, Element* element,
          const plain_volumecell_set& parent_cells, std::map<VolumeCell*, ChildCell>& cellmap);

      /// Generate IntegrationCells and BoundaryCells within the parent VolumeCell
      void Fill(Mesh& parent_mesh, Element* element, const plain_volumecell_set& parent_cells,
          std::map<VolumeCell*, ChildCell>& cellmap);

      /// Fill a parent cell with its child cells by means of the child cell topology.
      /*!
        - needs some seed child cells to start with
        - fails if there a flat rats that isolate different regions of the
          parent cell
       */
      void Fill(VolumeCell* parent_cell, ChildCell& childcell);

      /// find some (most) of the child cells for each parent cell
      void seed_cells(Mesh& parent_mesh, const plain_volumecell_set& parent_cells,
          std::map<VolumeCell*, ChildCell>& cellmap, plain_volumecell_set& done_child_cells);

      void build_surface_cell_map(VolumeCell* vc, ChildCell& cc);

      void register_new_points(Mesh& parent_mesh, const plain_volumecell_set& childset);

      /// put all volume cells at point to cell set
      void find_volume_cell(Point* p, plain_volumecell_set& childset);

      void to_parent(std::vector<Point*>& points) { swap_points(child_to_parent_, points); }
      void to_child(std::vector<Point*>& points) { swap_points(parent_to_child_, points); }

      void to_parent(Mesh& mesh, std::vector<Point*>& points)
      {
        swap_points(mesh, child_to_parent_, points);
      }
      void to_child(Mesh& mesh, std::vector<Point*>& points)
      {
        swap_points(mesh, parent_to_child_, points);
      }

      void to_parent(PointSet& points) { swap_points(child_to_parent_, points); }
      void to_child(PointSet& points) { swap_points(parent_to_child_, points); }

      Point* to_parent(Point* point) { return swap_point(child_to_parent_, point); }
      Point* to_child(Point* point) { return swap_point(parent_to_child_, point); }

      /// convert points between meshes and create points if not found
      void swap_points(
          Mesh& mesh, const std::map<Point*, Point*>& pointmap, std::vector<Point*>& points);

      /// convert points between meshes and raise error if point not found
      void swap_points(const std::map<Point*, Point*>& pointmap, std::vector<Point*>& points);

      /// convert points between meshes and raise error if point not found
      void swap_points(const std::map<Point*, Point*>& pointmap, PointSet& points);

      Point* swap_point(const std::map<Point*, Point*>& pointmap, Point* point);

      void Register(Point* parent_point, Point* child_point);

      void copy_cut_side(Side* s, Facet* f);

      Teuchos::RCP<PointPool> pp_;

      Mesh mesh_;
      Mesh cut_mesh_;

      std::map<Point*, Point*> parent_to_child_;
      std::map<Point*, Point*> child_to_parent_;

      std::map<Side*, std::vector<Side*>> side_parent_to_child_;
      // std::map<Side*, Side*> side_child_to_parent_;
    };
  }  // namespace CUT
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif
