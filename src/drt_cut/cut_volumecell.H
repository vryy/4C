#ifndef CUT_VOLUMECELL_H
#define CUT_VOLUMECELL_H

#include <map>
#include <vector>
#include <set>
#include <iostream>

#include "cut_point.H"
#include "cut_facet.H"

#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

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

  VolumeCell( const plain_facet_set & facets,
              const std::map<std::pair<Point*, Point*>, plain_facet_set > & volume_lines,
              Element * element );

  bool Contains( Point * p );

  void Neighbors( Point * p,
                  const plain_volumecell_set & cells,
                  const plain_volumecell_set & done,
                  plain_volumecell_set & connected,
                  plain_element_set & elements );

  Element * ParentElement() { return element_; }

  void GetAllPoints( Mesh & mesh, PointSet & cut_points );

  void CreateTet4IntegrationCells( Mesh & mesh, const std::vector<std::vector<Point*> > & tets, const std::map<Facet*, std::vector<Point*> > & sides_xyz );

  void GetIntegrationCells( plain_integrationcell_set & cells );

  void MomentFitGaussWeights(Element *elem);

  void GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells );

  void ConnectNodalDOFSets( bool include_inner );

  GEO::CUT::Point::PointPosition Position() const { return position_; }

  void Position( GEO::CUT::Point::PointPosition position );

  void Print( std::ostream & stream );

  const plain_facet_set & Facets() const { return facets_; }

  const plain_integrationcell_set & IntegrationCells() const { return integrationcells_; }

  const plain_boundarycell_set & BoundaryCells() const { return bcells_; }

  bool Empty() const { return integrationcells_.size()==0 and bcells_.size()==0; }

  void Disconnect();

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

  void SimplifyIntegrationCells( Mesh & mesh );

  /// calculate cell volume from integration cells
  double Volume();

  void TestSurface();

  int NumGaussPoints( DRT::Element::DiscretizationType shape );

  void DumpGmsh(const std::vector<std::vector<double> >&gauspts);

  //return the Gauss points computed using moment fitting equations
  Teuchos::RCP<DRT::UTILS::GaussPoints> GaussPointsFitting();
 // void GaussPointsFitting();
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
  plain_facet_set facets_;

  /// volume preserving simple shaped cells
  plain_integrationcell_set integrationcells_;

  /// boundary cells at cut surface with outward normals
  plain_boundarycell_set bcells_;

  /// dof set number of all element nodes
  std::vector<int> nodaldofset_;

  static int hex8totet4[5][4];
  static int wedge6totet4[3][4];
  static int pyramid5totet4[2][4];

 //store the Gauss point weights and location when moment fitting equations are
 //used to arrive at them
  Epetra_SerialDenseVector weights_;
  std::vector<std::vector<double> > gausPts_;
};


}
}

#endif
