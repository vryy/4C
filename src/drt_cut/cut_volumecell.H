#ifndef CUT_VOLUMECELL_H
#define CUT_VOLUMECELL_H

#include <map>
#include <vector>
#include <set>
#include <iostream>

#include "cut_point.H"

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

  void GetAllPoints( std::set<Point*> & cut_points );

  void CreateIntegrationCells( Mesh & mesh );

  void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  void SetIntegrationCells( Mesh & mesh,
                            const std::vector<std::vector<std::vector<Point*> >::iterator> & tets,
                            const std::map<std::vector<std::vector<Point*> >::iterator, std::vector<std::vector<Point*> > > & cut_sides );

  void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells );

  void ConnectNodalDOFSets();

  const std::vector<int> & NodalDofSet() const { return nodaldofset_; }

  GEO::CUT::Point::PointPosition Position();

  void Print( std::ostream & stream );

  bool OwnsFacet( Facet * f );

  bool OwnsTet( const std::vector<Point*> & tet );

  void NewTri3Cell( Mesh & mesh, Facet * f, const Epetra_SerialDenseMatrix & x );

  void NewTri3Cells( Mesh & mesh, Facet * f, const std::vector<Epetra_SerialDenseMatrix> & xyz );

  void NewQuad4Cell( Mesh & mesh, Facet * f, const Epetra_SerialDenseMatrix & x );

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
