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
  class Mesh;


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

  bool CreateSimpleShapeIntegrationCell( Mesh & mesh );

  bool CreateSingleElementIntegrationCell( Mesh & mesh );

  void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  const std::set<IntegrationCell*> & IntegrationCells() const { return integrationcells_; }

  const std::set<BoundaryCell*> & BoundaryCells() const { return bcells_; }

  void GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells );

  void ConnectNodalDOFSets( bool include_inner );

  const std::vector<int> & NodalDofSet() const { return nodaldofset_; }

  GEO::CUT::Point::PointPosition Position() const { return position_; }

  void Position( GEO::CUT::Point::PointPosition position );

  void Print( std::ostream & stream );

  const std::set<Facet*> & Facets() const { return facets_; }

  void NewIntegrationCell( Mesh & mesh, DRT::Element::DiscretizationType shape, const std::vector<Point*> & x );

  void NewHex8Cell( Mesh & mesh, const std::vector<Point*> & x );

  IntegrationCell * NewTet4Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewWedge6Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewPyramid5Cell( Mesh & mesh, const std::vector<Point*> & x );

  void NewBoundaryCell( Mesh & mesh, DRT::Element::DiscretizationType shape, Facet * f, const std::vector<Point*> & x );

  void NewTri3Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  void NewQuad4Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  double Volume();

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

  Element * element_;

  GEO::CUT::Point::PointPosition position_;

  std::set<Facet*> facets_;

  std::vector<Point*> points_;

  std::set<IntegrationCell*> integrationcells_;

  std::set<BoundaryCell*> bcells_;

  std::vector<int> nodaldofset_;

  static int hex8totet4[5][4];
  static int wedge6totet4[3][4];
  static int pyramid5totet4[2][4];
};

}
}

#endif
