/*---------------------------------------------------------------------*/
/*!
\file cut_tetmeshintersection.H

\brief Cut tets used for the Tessellation routine

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/

#ifndef CUT_TETMESHINTERSECTION_H
#define CUT_TETMESHINTERSECTION_H

#include "cut_mesh.H"

namespace GEO
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

      void Status();

      Mesh& NormalMesh() { return mesh_; }

     private:
      struct ChildCell
      {
        ChildCell() : done_(false), parent_(NULL) {}

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

      void FindEdgeCuts();

      /// find the mapping between child VolumeCells and parent VolumeCells.
      void MapVolumeCells(Mesh& parent_mesh, Element* element,
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
      void SeedCells(Mesh& parent_mesh, const plain_volumecell_set& parent_cells,
          std::map<VolumeCell*, ChildCell>& cellmap, plain_volumecell_set& done_child_cells);

      void BuildSurfaceCellMap(VolumeCell* vc, ChildCell& cc);

      void RegisterNewPoints(Mesh& parent_mesh, const plain_volumecell_set& childset);

      /// put all volume cells at point to cell set
      void FindVolumeCell(Point* p, plain_volumecell_set& childset);

      void ToParent(std::vector<Point*>& points) { SwapPoints(child_to_parent_, points); }
      void ToChild(std::vector<Point*>& points) { SwapPoints(parent_to_child_, points); }

      void ToParent(Mesh& mesh, std::vector<Point*>& points)
      {
        SwapPoints(mesh, child_to_parent_, points);
      }
      void ToChild(Mesh& mesh, std::vector<Point*>& points)
      {
        SwapPoints(mesh, parent_to_child_, points);
      }

      void ToParent(PointSet& points) { SwapPoints(child_to_parent_, points); }
      void ToChild(PointSet& points) { SwapPoints(parent_to_child_, points); }

      Point* ToParent(Point* point) { return SwapPoint(child_to_parent_, point); }
      Point* ToChild(Point* point) { return SwapPoint(parent_to_child_, point); }

      /// convert points between meshes and create points if not found
      void SwapPoints(
          Mesh& mesh, const std::map<Point*, Point*>& pointmap, std::vector<Point*>& points);

      /// convert points between meshes and raise error if point not found
      void SwapPoints(const std::map<Point*, Point*>& pointmap, std::vector<Point*>& points);

      /// convert points between meshes and raise error if point not found
      void SwapPoints(const std::map<Point*, Point*>& pointmap, PointSet& points);

      Point* SwapPoint(const std::map<Point*, Point*>& pointmap, Point* point);

      void Register(Point* parent_point, Point* child_point);

      void CopyCutSide(Side* s, Facet* f);

      Teuchos::RCP<PointPool> pp_;

      Mesh mesh_;
      Mesh cut_mesh_;

      std::map<Point*, Point*> parent_to_child_;
      std::map<Point*, Point*> child_to_parent_;

      std::map<Side*, std::vector<Side*>> side_parent_to_child_;
      // std::map<Side*, Side*> side_child_to_parent_;
    };
  }  // namespace CUT
}  // namespace GEO

#endif
