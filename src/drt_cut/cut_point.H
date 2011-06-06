#ifndef CUT_POINT_H
#define CUT_POINT_H

#include <fstream>
#include <iostream>

#include <map>
#include <set>
#include <vector>
#include <stdexcept>

#include <Teuchos_RCP.hpp>

#include "../linalg/linalg_fixedsizematrix.H"

namespace GEO
{
namespace CUT
{

class Node;
class Edge;
class Side;
class Element;

class Line;
class Facet;

class Mesh;

class PointLineFilter
{
public:

  virtual ~PointLineFilter() {}

  virtual bool operator()( Line * line ) const = 0;
};

class Point
{
public:

  enum PointPosition
  {
    undecided,
    oncutsurface = -1,
    inside = -2,
    outside = -3
  };

  static Point * NewPoint( Mesh & mesh, const double * x, double t, Edge * cut_edge, Side * cut_side );

  static Point * InsertCut( Edge * cut_edge, Side * cut_side, Node * n );

  Point( unsigned pid, const double * x, Edge * cut_edge, Side * cut_side );

  int Id() const { return pid_; }

  const double * X() const { return x_; }

  PointPosition Position() const { return position_; }

  void Position( PointPosition p );

  void AddEdge( Edge* e );

  void AddSide( Side* s );

  void AddElement( Element* e ) { if ( e!=NULL ) { cut_elements_.insert( e ); } }

  bool IsCut( Edge* e )
  {
    return cut_edges_.count( e )>0;
  }

  bool IsCut( Facet* f )
  {
    return facets_.count( f )>0;
  }

  bool IsCut( Side* s )
  {
    return cut_sides_.count( s )>0;
  }

  bool IsCut( Element* s )
  {
    return cut_elements_.count( s )>0;
  }

  void Coordinates( double * x )
  {
    std::copy( x_, x_+3, x );
  }

  void Register( Line * line )
  {
    lines_.insert( line );
  }

  void Register( Facet * facet )
  {
    facets_.insert( facet );
  }

  void CommonEdge( Point * other, std::set<Edge *> & edges );

  void CommonSide( Point * other, std::set<Side *> & sides );

  Line * CommonLine( Point * other );

  Line * CutLine( Side * side, bool unique=true );

  Line * CutLine( Line * line, Side * side, Element * element );

  Line * CutLine( const PointLineFilter & filter, bool unique=true );

  Line * CutLine( Line * line, const PointLineFilter & filter, bool unique=true );

  void CutLines( const PointLineFilter & filter, std::set<Line*> & cut_lines );

  void CutLines( Side * side, std::set<Line*> & cut_lines );

  //std::vector<Edge*> CutEdges( Point * other );

  const std::set<Edge*> & CutEdges() { return cut_edges_; }

  const std::set<Side*> & CutSides() { return cut_sides_; }

  Side * CutSide( Side * side, Point * other );

  void CutEdge( Side * side, Line * other_line, std::vector<Edge*> & matches );

  void Print( std::ostream & stream=std::cout )
  {
    stream << "("
           << pid_ << "; "
           << x_[0] << ","
           << x_[1] << ","
           << x_[2] << ")";
  }

  void Plot( std::ostream & f )
  {
    f << std::setprecision( 10 ) << x_[0] << " "
      << std::setprecision( 10 ) << x_[1] << " "
      << std::setprecision( 10 ) << x_[2] << " "
      << pid_
      << "\n";
  }

  double t( Edge* edge );

  //void t( Edge* edge, double pos ) { t_[edge] = pos; }

  unsigned Pid() const { return pid_; }

  bool NodalPoint( const std::vector<Node*> & nodes ) const;

  Node * CutNode();

  void Intersection( std::set<Side*> & sides );

  void Intersection( std::set<Facet*> & facets );

  const std::set<Element*> & Elements() const { return cut_elements_; }

  const std::set<Facet*> & Facets() const { return facets_; }

  Edge * CommonCutEdge( Side * side );

private:

  // find the unique line that matches the filter criteria
  template <class Filter>
  Line * Find( Filter & filter, bool unique=true );

  template <class Filter>
  void Find( Filter & filter, std::set<Line*> & cut_lines );

  unsigned pid_;
  PointPosition position_;

  double x_[3];
  std::map<Edge*, double> t_;

  std::set<Edge*> cut_edges_;
  std::set<Side*> cut_sides_;
  std::set<Element*> cut_elements_;

  std::set<Line*> lines_;
  std::set<Facet*> facets_;
};

inline std::ostream & operator<<( std::ostream & stream, Point & point )
{
  point.Print( stream );
  return stream;
}

inline std::ostream & operator<<( std::ostream & stream, Point * point )
{
  point->Print( stream );
  return stream;
}

inline int EntityId( const Point & p ) { return p.Pid(); }

template <class T>
class EntityIdLess
{
public:

  bool operator()( const T& p1, const T& p2 ) const
  {
    return EntityId( p1 ) < EntityId( p2 );
  }

  bool operator()( const T* p1, const T* p2 ) const
  {
    return EntityId( *p1 ) < EntityId( *p2 );
  }

  bool operator()( const Teuchos::RCP<T> p1, const Teuchos::RCP<T> p2 ) const
  {
    return EntityId( *p1 ) < EntityId( *p2 );
  }

};

typedef EntityIdLess<Point> PointPidLess;
typedef std::set<Point*, PointPidLess> PointSet;

class PointPositionLess
{
public:

  explicit PointPositionLess( Edge* edge ) : edge_( edge ) {}

  bool operator()( Point& p1, Point& p2 ) const
  {
    if ( not p1.IsCut( edge_ ) or not p2.IsCut( edge_ ) )
      throw std::runtime_error( "point position compare only on cut edges" );
    return p1.t( edge_ ) < p2.t( edge_ );
  }

  bool operator()( Point* p1, Point* p2 ) const
  {
    if ( not p1->IsCut( edge_ ) or not p2->IsCut( edge_ ) )
      throw std::runtime_error( "point position compare only on cut edges" );
    return p1->t( edge_ ) < p2->t( edge_ );
  }

private:
  Edge* edge_;
};

typedef std::set<Point*, PointPositionLess> PointPositionSet;

class ExcludeLineFilter : public PointLineFilter
{
public:

  ExcludeLineFilter( Line * line, const PointLineFilter & filter )
    : line_( line ),
      filter_( filter )
  {
  }

  virtual bool operator()( Line * line ) const
  {
    return line!=line_ and filter_( line );
  }

private:
  Line * line_;
  const PointLineFilter & filter_;
};

class SideCutFilter : public PointLineFilter
{
public:

  SideCutFilter( Side * side )
    : side_( side )
  {
  }

  virtual bool operator()( Line * line ) const;

private:
  Side * side_;
};

class SideElementCutFilter : public PointLineFilter
{
public:

  SideElementCutFilter( Side * side, Element * element )
    : side_( side ),
      element_( element )
  {
  }

  virtual bool operator()( Line * line ) const;

private:
  Side * side_;
  Element * element_;
};

class SideSideCutFilter : public PointLineFilter
{
public:

  SideSideCutFilter( Side * side1, Side * side2 )
    : side1_( side1 ),
      side2_( side2 )
  {
  }

  virtual bool operator()( Line * line ) const;

private:
  Side * side1_;
  Side * side2_;
};

inline void FindCommonSides( const std::vector<Point*> & side, std::set<Side*> & sides )
{
  std::vector<Point*>::const_iterator is = side.begin();
  sides = ( *is )->CutSides();
  for ( ++is; is!=side.end(); ++is )
  {
    Point * p = *is;
    p->Intersection( sides );
    if ( sides.size()==0 )
    {
      break;
    }
  }
}

inline void FindCommonSides( Point * p1, Point * p2, Point * p3, std::set<Side*> & sides )
{
  sides = p1->CutSides();
  p2->Intersection( sides );
  p3->Intersection( sides );
}

inline void FindCommonSides( Point * p1, Point * p2, Point * p3, Point * p4, std::set<Side*> & sides )
{
  sides = p1->CutSides();
  p2->Intersection( sides );
  p3->Intersection( sides );
  p4->Intersection( sides );
}

}
}

#endif
