#ifndef CUT_TETMESHINTERSECTION_H
#define CUT_TETMESHINTERSECTION_H

#include <vector>
#include <set>

#include "cut_mesh.H"
#include "cut_levelsetside.H"

namespace GEO
{
  namespace CUT
  {
    class Side;
    class Mesh;
    class Options;
    class TetMesh;
    class VolumeCell;

    class TetMeshIntersection
    {
    public:

      TetMeshIntersection( const Options & options,
                           Element * element,
                           const std::vector<std::vector<int> > & tets,
                           const std::vector<int> & accept_tets,
                           const std::vector<Point*> & points,
                           const std::set<Side*> & cut_sides,
                           bool levelset );

      void Cut( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, bool levelset );

      void Status();

      Mesh & NormalMesh() { return mesh_; }

    private:

      struct ChildCell
      {
        ChildCell()
          : done_( false ),
            parent_( NULL )
        {
        }

        bool ContainsChild( VolumeCell * vc )
        {
          return cells_.count( vc ) > 0;
        }

        bool done_;
        VolumeCell * parent_;
        std::set<VolumeCell*> cells_;
        std::map<Side*, std::vector<Facet*> > facetsonsurface_;
      };

      void MapVolumeCells( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, std::map<VolumeCell*, ChildCell> & cellmap );

      void Fill( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, std::map<VolumeCell*, ChildCell> & cellmap, bool levelset );

      void Fill( VolumeCell * parent_cell, ChildCell & childcell );

      void RegisterNewPoints( Mesh & parent_mesh, const std::set<VolumeCell*> & childset );

      void FindVolumeCell( Point * p, std::set<VolumeCell*> & childset );

      void ToParent( std::vector<Point*> & points ) { SwapPoints( child_to_parent_, points ); }
      void ToChild ( std::vector<Point*> & points ) { SwapPoints( parent_to_child_, points ); }

      void ToParent( std::set<Point*> & points ) { SwapPoints( child_to_parent_, points ); }
      void ToChild ( std::set<Point*> & points ) { SwapPoints( parent_to_child_, points ); }

      Point * ToParent( Point * point ) { return SwapPoint( child_to_parent_, point ); }
      Point * ToChild ( Point * point ) { return SwapPoint( parent_to_child_, point ); }

      void SwapPoints( const std::map<Point*, Point*> & pointmap, std::vector<Point*> & points );
      void SwapPoints( const std::map<Point*, Point*> & pointmap, std::set<Point*> & points );

      Point * SwapPoint( const std::map<Point*, Point*> & pointmap, Point * point );

      void Register( Point * parent_point, Point * child_point );

      Teuchos::RCP<PointPool> pp_;

      Mesh mesh_;
      Mesh cut_mesh_;

      std::map<Point*, Point*> parent_to_child_;
      std::map<Point*, Point*> child_to_parent_;

      std::map<Side*, std::vector<Side*> > side_parent_to_child_;
      //std::map<Side*, Side*> side_child_to_parent_;

      LevelSetSide ls_side_;
    };
  }
}

#endif
