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

class Facet
{
public:

  Facet( Mesh & mesh, const std::vector<Point*> & points, Side * side, bool cutsurface );

  void Register( VolumeCell * cell ) { cells_.insert( cell ); }

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

  void GetAllPoints( Mesh & mesh, std::set<Point*> & cut_points );

#if 0
  void CreateLinearElements( Mesh & mesh );
#endif

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

  bool Equals( const std::vector<Point*> & facet_points );

  bool Equals( DRT::Element::DiscretizationType distype );

  bool IsCutSide( Side * side );

  Point::PointPosition Position() const { return position_; }

  void Position( Point::PointPosition p );

  void GetLines( std::map<std::pair<Point*, Point*>, std::set<Facet*> > & lines );

  bool IsLine( Point * p1, Point * p2 );

  bool Contains( Point * p ) const;

  bool Contains( const std::vector<Point*> & side ) const;

  bool Touches( Facet * f );

  VolumeCell * Neighbor( VolumeCell * cell );

  void Neighbors( Point * p,
                  const std::set<VolumeCell*> & cells,
                  const std::set<VolumeCell*> & done,
                  std::set<VolumeCell*> & connected,
                  std::set<Element*> & elements );

  const std::vector<Point*> & Points() const { return points_; }

  const std::vector<Point*> & CornerPoints() const { return corner_points_; }

  const std::vector<std::vector<Point*> > & Triangulation() const { return triangulation_; }

  void TriangulationPoints( std::set<Point*> & points );

  void AllPoints( std::set<Point*> & points )
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

  void NewTri3Cell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, std::set<BoundaryCell*> & bcells );

//  void NewTri3Cells( Mesh & mesh, VolumeCell * volume, const std::vector<Epetra_SerialDenseMatrix> & xyz, std::set<BoundaryCell*> & bcells );

  void NewQuad4Cell( Mesh & mesh, VolumeCell * volume, const std::vector<Point*> & points, std::set<BoundaryCell*> & bcells );

  void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells );

  bool IsTriangle( const std::vector<Point*> & tri ) const;

  bool IsTriangulated() const { return triangulation_.size() > 0; }

  bool IsTriangulatedSide( const std::vector<Point*> & tri ) const;

  bool HasHoles() const { return holes_.size() > 0; }

  const std::set<Facet*> & Holes() const { return holes_; }

  unsigned NumPoints();

  const std::set<VolumeCell*> & Cells() const { return cells_; }

private:

  Facet( const Facet & );
  Facet & operator=( const Facet & );

  bool IsPlanar( Mesh & mesh );

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

  unsigned FindNextCornerPoint( const std::vector<Point*> & points,
                                LINALG::Matrix<3,1> & x1,
                                LINALG::Matrix<3,1> & x2,
                                LINALG::Matrix<3,1> & x3,
                                LINALG::Matrix<3,1> & b1,
                                LINALG::Matrix<3,1> & b2,
                                LINALG::Matrix<3,1> & b3,
                                unsigned i );

  bool IsLine( const std::vector<Point*> & points, Point * p1, Point * p2 );

  void GetLines( const std::vector<Point*> & points,
                 std::map<std::pair<Point*, Point*>, std::set<Facet*> > & lines );


  std::vector<Point*> points_;

  std::vector<Point*> corner_points_;

  std::set<Facet*> holes_;

  std::vector<std::vector<Point*> > triangulation_;

  Side* parentside_;

  bool planar_;

  bool planar_known_;

  Point::PointPosition position_;

  std::set<VolumeCell*> cells_;
};


  inline Facet * FindFacet( const std::set<Facet*> & facets, const std::vector<Point*> & side )
  {
    Facet * found = NULL;
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Equals( side ) )
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

inline void RemoveNonmatchingTriangulatedFacets( const std::vector<Point*> & side, std::set<Facet*> & facets )
{
  if ( side.size()==3 )
  {
    for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); )
    {
      Facet * f = *i;
      if ( f->IsTriangulated() )
      {
        if ( not f->IsTriangulatedSide( side ) )
        {
          facets.erase( i++ );
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

inline void FindCommonFacets( const std::vector<Point*> & side, std::set<Facet*> & facets )
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

inline void FindCommonFacets( Point * p1, Point * p2, Point * p3, std::set<Facet*> & facets )
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

inline void FindCommonFacets( Point * p1, Point * p2, Point * p3, Point * p4, std::set<Facet*> & facets )
{
  facets = p1->Facets();
  p2->Intersection( facets );
  p3->Intersection( facets );
  p4->Intersection( facets );
}

}
}

#endif
