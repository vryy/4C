#ifndef CUT_SIDE_H
#define CUT_SIDE_H

#include "cut_node.H"
#include "cut_edge.H"
#include "cut_tolerance.H"

namespace GEO
{
namespace CUT
{
  class Facet;
  class LineSegment;

class Side
{
public:

  Side( int sid,
        const std::vector<Node*> & nodes,
        const std::vector<Edge*> & edges )
    : sid_( sid ),
      nodes_( nodes ),
      edges_( edges )
  {
  }

  virtual ~Side() {}

  int Id() const { return sid_; }

  // debug only
  void SetId( int sid ) { sid_ = sid; }

  void Register( Element * element ) { elements_.insert( element ); }

  void CollectElements( std::set<Element*> & elements )
  {
    std::copy( elements_.begin(), elements_.end(), std::inserter( elements, elements.begin() ) );
  }

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points ) = 0;

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  virtual void FillComplete( Mesh & mesh );

  virtual bool Cut( Mesh & mesh, Side & side, Element * element ) = 0;

  void EdgeCuts( Mesh & mesh, Side & side, std::set<Point*, PointPidLess> & cuts );

  void Coordinates( double * xyze )
  {
    double * x = xyze;
    for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
    {
      Node & n = **i;
      n.Coordinates( x );
      x += 3;
    }
  }

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets ) = 0;

  virtual void MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets ) = 0;

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets ) = 0;

  virtual bool IsCut() = 0;

  //bool FullSideCut() { return cut_lines_.size()==edges_.size() and facets_.size()==1; }

  bool OnSide( const std::set<Point*, PointPidLess> & points );

  bool OnEdge( Point * point );

  bool OnEdge( Line * line );

  virtual void ExchangeFacetSide( Side * side ) = 0;

  virtual void AddLine( Line* cut_line ) = 0;

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points ) = 0;

  void Print();

  Node * OnNode( const LINALG::Matrix<3,1> & x );

  const std::vector<Edge*> & Edges() const { return edges_; }

  Edge * FindEdge( Point * begin, Point * end );

  const std::vector<Node*> & Nodes() const { return nodes_; }

protected:

private:

  int sid_;

  std::vector<Node*> nodes_;

  std::vector<Edge*> edges_;

  std::set<Element*> elements_;
};


class LinearSide : public Side
{
public:

  LinearSide( int sid,
              const std::vector<Node*> & nodes,
              const std::vector<Edge*> & edges )
    : Side( sid, nodes, edges )
  {
    for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      e->Register( this );
    }
  }

  virtual bool Cut( Mesh & mesh, Side & side, Element * element );

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual bool IsCut();

  virtual void ExchangeFacetSide( Side * side );

  // add cut_line to cut_lines_ if not yet included
  virtual void AddLine( Line* cut_line );

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points );

  Teuchos::RCP<LineSegment> CreateLineSegment( Mesh & mesh, Element * element );

protected:

  void CreateLineSegmentList( Mesh & mesh,
                              Element * element,
                              std::vector<Teuchos::RCP<LineSegment> > & segments,
                              bool inner=false );

  std::vector<Line*> cut_lines_;

  std::vector<Facet*> facets_;
};


class QuadraticSide : public Side
{
public:

  QuadraticSide( int sid,
                 const std::vector<Node*> & nodes,
                 const std::vector<Edge*> & edges )
    : Side( sid, nodes, edges )
  {}

  //virtual void FillComplete( Mesh & mesh );

  virtual bool Cut( Mesh & mesh, Side & side, Element * element );

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual bool IsCut();

  virtual void ExchangeFacetSide( Side * side )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void AddLine( Line* cut_line )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points )
  { throw std::runtime_error( "not supposed to end up here" ); }

protected:

  std::vector<Side*> subsides_;
};


template <>
class ConcreteSide<DRT::Element::tri3> : public LinearSide
{
public:

  //const int numNodesSurface = 3;

  ConcreteSide<DRT::Element::tri3>( int sid,
                                    const std::vector<Node*> & nodes,
                                    const std::vector<Edge*> & edges )
    : LinearSide( sid, nodes, edges )
  {}

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  {
    return edge.Cut( mesh, *this, cut_points );
  }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  {
    if ( fabs( s     ) < TOLERANCE ) edges.push_back( Edges()[0] );
    if ( fabs( r+s-1 ) < TOLERANCE ) edges.push_back( Edges()[1] );
    if ( fabs( r     ) < TOLERANCE ) edges.push_back( Edges()[2] );
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 3> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


template <>
class ConcreteSide<DRT::Element::quad4> : public LinearSide
{
public:

  //const int numNodesSurface = 4;

  ConcreteSide<DRT::Element::quad4>( int sid,
                                     const std::vector<Node*> & nodes,
                                     const std::vector<Edge*> & edges )
    : LinearSide( sid, nodes, edges )
  {}

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  {
    return edge.Cut( mesh, *this, cut_points );
  }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  {
    if ( fabs( s+1 ) < TOLERANCE ) edges.push_back( Edges()[0] );
    if ( fabs( r-1 ) < TOLERANCE ) edges.push_back( Edges()[1] );
    if ( fabs( s-1 ) < TOLERANCE ) edges.push_back( Edges()[2] );
    if ( fabs( r+1 ) < TOLERANCE ) edges.push_back( Edges()[3] );
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 4> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


template <>
class ConcreteSide<DRT::Element::tri6> : public QuadraticSide
{
public:

  //const int numNodesSurface = 6;

  ConcreteSide<DRT::Element::tri6>( int sid,
                                    const std::vector<Node*> & nodes,
                                    const std::vector<Edge*> & edges )
    : QuadraticSide( sid, nodes, edges )
  {}

  virtual void FillComplete( Mesh & mesh );

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 6> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


template <>
class ConcreteSide<DRT::Element::quad8> : public QuadraticSide
{
public:

  //const int numNodesSurface = 8;

  ConcreteSide<DRT::Element::quad8>( int sid,
                                     const std::vector<Node*> & nodes,
                                     const std::vector<Edge*> & edges )
    : QuadraticSide( sid, nodes, edges )
  {}

  virtual void FillComplete( Mesh & mesh );

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 8> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


template <>
class ConcreteSide<DRT::Element::quad9> : public QuadraticSide
{
public:

  //const int numNodesSurface = 9;

  ConcreteSide<DRT::Element::quad9>( int sid,
                                     const std::vector<Node*> & nodes,
                                     const std::vector<Edge*> & edges )
    : QuadraticSide( sid, nodes, edges )
  {}

  virtual void FillComplete( Mesh & mesh );

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  { throw std::runtime_error( "not supposed to end up here" ); }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 9> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


}
}

#endif
