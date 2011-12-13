#ifndef CUT_FACET_H
#define CUT_FACET_H

#include <fstream>
#include <iostream>

#include <map>
#include <set>
#include <vector>
#include <stdexcept>

#include "../linalg/linalg_fixedsizematrix.H"

#include "../drt_lib/drt_element.H"
#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

// #include "../drt_lib/drt_element.H"

#include "cut_point.H"
#include "cut_line.H"

namespace GEO
{

namespace CUT
{
class BoundaryCell;
class VolumeCell;
class Element;
class Mesh;
class Facet;

//typedef EntityIdLess<Facet> FacetPidLess;

/// Surface of arbitrary shape, defined by its points.
class Facet
{
public:

  Facet( Mesh & mesh, const std::vector<Point*> & points, Side * side, bool cutsurface );

  void Register( VolumeCell * cell );

  void DisconnectVolume( VolumeCell * cell );

  void Print()
  {
    for ( std::vector<Point*>::iterator i=points_.begin();
          i!=points_.end();
          ++i )
    {
      Point & p = **i;
      std::cout << p.Pid() << " ";
    }
    std::cout << "\n";
  }

  void Print( std::ostream & stream );

  bool OnCutSide() { return SideId() > -1; }

  int SideId();

  int PositionSideId();

  Side* ParentSide() const { return parentside_; }

  void Coordinates( double * x );

  void CornerCoordinates( double * x );

  void GetAllPoints( Mesh & mesh, PointSet & cut_points, bool dotriangulate=false );

  void AddHole( Facet * hole );

  void ExchangeSide( Side * side, bool cutsurface )
  {
    parentside_ = side;
    if ( cutsurface )
    {
      Position( Point::oncutsurface );
      for ( std::vector<Point*>::const_iterator i=points_.begin(); i!=points_.end(); ++i )
      {
        Point * p = *i;
        p->Position( Point::oncutsurface );
      }
    }
  }

  bool Equals( const std::vector<Point*> & facet_points ) { return Equals( points_, facet_points ); }

  bool Equals( DRT::Element::DiscretizationType distype );

  bool CornerEquals( const std::vector<Point*> & facet_points ) { return Equals( corner_points_, facet_points ); }

  bool IsCutSide( Side * side );

  Point::PointPosition Position() const { return position_; }

  void Position( Point::PointPosition p );

  void GetLines( std::map<std::pair<Point*, Point*>, plain_facet_set> & lines );

  bool IsLine( Point * p1, Point * p2 );

  bool Contains( Point * p ) const;

  bool Contains( const std::vector<Point*> & side ) const;

  bool ContainsSome( const std::vector<Point*> & side ) const;

  bool Touches( Facet * f );

  VolumeCell * Neighbor( VolumeCell * cell );

  void Neighbors( Point * p,
                  const plain_volumecell_set & cells,
                  const plain_volumecell_set & done,
                  plain_volumecell_set & connected,
                  plain_element_set & elements );

  void Neighbors( Point * p,
                  const plain_volumecell_set & cells,
                  const plain_volumecell_set & done,
                  plain_volumecell_set & connected);

  const std::vector<Point*> & Points() const { return points_; }

  const std::vector<Point*> CornerPoints() const { return corner_points_; }

  const std::vector<std::vector<double> > CornerPointsLocal(Element *elem1,bool print);

  const std::vector<std::vector<Point*> > & Triangulation() const { return triangulation_; }

  void TriangulationPoints( PointSet & points );

  void AllPoints( PointSet & points )
  {
    if ( IsTriangulated() )
    {
      TriangulationPoints( points );
    }
    else
    {
      std::copy( points_.begin(), points_.end(), std::inserter( points, points.begin() ) );
    }
  }

  void NewTri3Cell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, plain_boundarycell_set & bcells );

//  void NewTri3Cells( Mesh & mesh, VolumeCell * volume, const std::vector<Epetra_SerialDenseMatrix> & xyz, plain_boundarycell_set & bcells );

  void NewQuad4Cell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, plain_boundarycell_set & bcells );

  void NewArbitraryCell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, plain_boundarycell_set & bcells,
      const DRT::UTILS::GaussIntegration& gp );

  void GetBoundaryCells( plain_boundarycell_set & bcells );

  void TestFacetArea( double tolerance );

  bool IsTriangle( const std::vector<Point*> & tri ) const;

  bool IsTriangulated() const { return triangulation_.size() > 0; }

  bool IsTriangulatedSide( const std::vector<Point*> & tri ) const;

  bool HasHoles() const { return holes_.size() > 0; }

  const plain_facet_set & Holes() const { return holes_; }

  unsigned NumPoints();

  const plain_volumecell_set & Cells() const { return cells_; }

  Point * OtherPoint( Point * p1, Point * p2 );

  void DoTriangulation(Mesh & mesh, const std::vector<Point*> & points ){CreateTriangulation(mesh,points );}

private:

  Facet( const Facet & );
  Facet & operator=( const Facet & );

  bool IsPlanar( Mesh & mesh, bool dotriangulate );

  bool IsPlanar( Mesh & mesh, const std::vector<Point*> & points );

  void CreateTriangulation( Mesh & mesh, const std::vector<Point*> & points );

  void GetNodalIds( Mesh & mesh, const std::vector<Point*> & points, std::vector<int> & nids );

  unsigned Normal( const std::vector<Point*> & points,
                   LINALG::Matrix<3,1> & x1,
                   LINALG::Matrix<3,1> & x2,
                   LINALG::Matrix<3,1> & x3,
                   LINALG::Matrix<3,1> & b1,
                   LINALG::Matrix<3,1> & b2,
                   LINALG::Matrix<3,1> & b3 );

  void FindCornerPoints();

  bool IsLine( const std::vector<Point*> & points, Point * p1, Point * p2 );

  void GetLines( const std::vector<Point*> & points,
                 std::map<std::pair<Point*, Point*>, plain_facet_set> & lines );

  bool Equals( const std::vector<Point*> & my_points, const std::vector<Point*> & facet_points );

  std::vector<Point*> points_;

  std::vector<Point*> corner_points_;

  plain_facet_set holes_;

  std::vector<std::vector<Point*> > triangulation_;

  Side* parentside_;

  bool planar_;

  bool planar_known_;

  Point::PointPosition position_;

  plain_volumecell_set cells_;
};

  inline std::ostream & operator<<( std::ostream & stream, Facet & f )
  {
    stream << "facet: {";
    if ( f.IsTriangulated() )
    {
      const std::vector<std::vector<Point*> > & triangulation = f.Triangulation();
      for ( std::vector<std::vector<Point*> >::const_iterator i=triangulation.begin();
            i!=triangulation.end();
            ++i )
      {
        const std::vector<Point*> & tri = *i;
        stream << "{";
        for ( std::vector<Point*>::const_iterator i=tri.begin(); i!=tri.end(); ++i )
        {
          Point * p = *i;
          p->Print( stream );
          stream << ",";
        }
        stream << "},";
      }
    }
    else
    {
      const std::vector<Point*> & points = f.Points();
      for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
      {
        Point * p = *i;
        p->Print( stream );
        stream << ",";
      }
      if ( f.HasHoles() )
      {
        const plain_facet_set & holes = f.Holes();
        for ( plain_facet_set::const_iterator i=holes.begin(); i!=holes.end(); ++i )
        {
          Facet & h = **i;
          stream << h;
        }
      }
    }
    stream << "}";
    return stream;
  }

  inline Facet * FindFacet( const plain_facet_set & facets, const std::vector<Point*> & side )
  {
    Facet * found = NULL;
    for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->CornerEquals( side ) )
      {
        if ( found==NULL )
        {
          found = f;
        }
        else
        {
          throw std::runtime_error( "not unique" );
        }
      }
    }
    return found;
  }

//inline int EntityId( const Facet & f ) { return f.InternalId(); }
//inline int EntityId( const Facet & f ) { return reinterpret_cast<int>( &f ); }

inline void RemoveNonmatchingTriangulatedFacets( const std::vector<Point*> & side, plain_facet_set & facets )
{
  if ( side.size()==3 )
  {
    for ( plain_facet_set::iterator i=facets.begin(); i!=facets.end(); )
    {
      Facet * f = *i;
      if ( f->IsTriangulated() )
      {
        if ( not f->IsTriangulatedSide( side ) )
        {
          set_erase( facets, i );
        }
        else
        {
          ++i;
        }
      }
      else
      {
        ++i;
      }
    }
  }
}

inline void FindCommonFacets( const std::vector<Point*> & side, plain_facet_set & facets )
{
  std::vector<Point*>::const_iterator is = side.begin();
  facets = ( *is )->Facets();
  for ( ++is; is!=side.end(); ++is )
  {
    Point * p = *is;
    p->Intersection( facets );
    if ( facets.size()==0 )
    {
      break;
    }
  }
  RemoveNonmatchingTriangulatedFacets( side, facets );
}

inline void FindCommonFacets( Point * p1, Point * p2, Point * p3, plain_facet_set & facets )
{
  facets = p1->Facets();
  p2->Intersection( facets );
  p3->Intersection( facets );

  std::vector<Point*> side( 3 );
  side[0] = p1;
  side[1] = p2;
  side[2] = p3;
  RemoveNonmatchingTriangulatedFacets( side, facets );
}

inline void FindCommonFacets( Point * p1, Point * p2, Point * p3, Point * p4, plain_facet_set & facets )
{
  facets = p1->Facets();
  p2->Intersection( facets );
  p3->Intersection( facets );
  p4->Intersection( facets );
}

}
}

#endif
