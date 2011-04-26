#ifndef CUT_VOLUMECELL_H
#define CUT_VOLUMECELL_H

#include <map>
#include <vector>
#include <set>
#include <iostream>

#include "cut_point.H"
#include "cut_facet.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Facet;
  class Element;
  class BoundaryCell;
  class IntegrationCell;
  class IntegrationCellCreator;
  class Mesh;

/// a partial volume of a cutted element of any shape
/*!
  A volume cell is a partial (or whole) element that is defined by its
  facets. It knows

  - any number of integration cells (that have simple shapes)
  - any number of boundary cells (at the cut facets with outward normals)
  - its position to the cut surface
  - the dof set number of all element nodes

 */
class VolumeCell
{
public:

  VolumeCell( const std::set<Facet*> & facets,
              const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
              Element * element );

  bool Contains( Point * p );

  void Neighbors( Point * p,
                  const std::set<VolumeCell*> & cells,
                  const std::set<VolumeCell*> & done,
                  std::set<VolumeCell*> & connected,
                  std::set<Element*> & elements );

  Element * ParentElement() { return element_; }

  void GetAllPoints( Mesh & mesh, std::set<Point*> & cut_points );

  void CreateIntegrationCells( Mesh & mesh );

  void CreateTet4IntegrationCells( Mesh & mesh, Point::PointPosition position, const std::vector<Point*> & points, const std::set<Facet*> & facets, bool project );

  void CreateTet4IntegrationCells( Mesh & mesh, const std::vector<std::vector<Point*> > & tets, const std::map<Facet*, std::vector<Point*> > & sides_xyz );

  void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  void GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells );

  void ConnectNodalDOFSets( bool include_inner );

  GEO::CUT::Point::PointPosition Position() const { return position_; }

  void Position( GEO::CUT::Point::PointPosition position );

  void Print( std::ostream & stream );

  const std::set<Facet*> & Facets() const { return facets_; }

  const std::set<IntegrationCell*> & IntegrationCells() const { return integrationcells_; }

  const std::set<BoundaryCell*> & BoundaryCells() const { return bcells_; }

  const std::vector<int> & NodalDofSet() const { return nodaldofset_; }

  /// create integration cell of specified shape
  void NewIntegrationCell( Mesh & mesh, DRT::Element::DiscretizationType shape, const std::vector<Point*> & x );

  void NewHex8Cell( Mesh & mesh, const std::vector<Point*> & x );

  IntegrationCell * NewTet4Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewWedge6Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewPyramid5Cell( Mesh & mesh, const std::vector<Point*> & x );

  /// create boundary cell of specified shape
  void NewBoundaryCell( Mesh & mesh, DRT::Element::DiscretizationType shape, Facet * f, const std::vector<Point*> & x );

  void NewTri3Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  void NewQuad4Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  /// calculate cell volume from integration cells
  double Volume();

  int NumGaussPoints( DRT::Element::DiscretizationType shape );

private:

  void SetTetPoints( const int * totet4,
                     const std::vector<Point*> & points,
                     std::vector<Point*> & tet4_points ) const
  {
    for ( int i=0; i<4; ++i )
    {
      tet4_points[i] = points[totet4[i]];
    }
  }

  /// the element this is a part of
  Element * element_;

  /// position relative to cut surface
  GEO::CUT::Point::PointPosition position_;

  /// defining facets
  std::set<Facet*> facets_;

  /// all points of the volume
  std::vector<Point*> points_;

  /// volume preserving simple shaped cells
  std::set<IntegrationCell*> integrationcells_;

  /// boundary cells at cut surface with outward normals
  std::set<BoundaryCell*> bcells_;

  /// dof set number of all element nodes
  std::vector<int> nodaldofset_;

  static int hex8totet4[5][4];
  static int wedge6totet4[3][4];
  static int pyramid5totet4[2][4];
};


}
}

#endif
