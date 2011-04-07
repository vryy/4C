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

  void CreateTet4IntegrationCells( Mesh & mesh, Point::PointPosition position, const std::vector<Point*> & points, const std::set<Facet*> & facets );

  void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  const std::set<IntegrationCell*> & IntegrationCells() const { return integrationcells_; }

  const std::set<BoundaryCell*> & BoundaryCells() const { return bcells_; }

  void GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells );

  void ConnectNodalDOFSets( bool include_inner );

  const std::vector<int> & NodalDofSet() const { return nodaldofset_; }

  GEO::CUT::Point::PointPosition Position();

  void Print( std::ostream & stream );

  const std::set<Facet*> & Facets() const { return facets_; }

  void NewTri3Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

//  void NewTri3Cells( Mesh & mesh, Facet * f, const std::vector<Epetra_SerialDenseMatrix> & xyz );

  void NewQuad4Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x );

  double Volume();

private:

  Element * element_;

  std::set<Facet*> facets_;

  std::vector<Point*> points_;

  std::set<IntegrationCell*> integrationcells_;

  std::set<BoundaryCell*> bcells_;

  std::vector<int> nodaldofset_;
};

}
}

#endif
