#ifndef CUT_SIDE_H
#define CUT_SIDE_H

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "cut_node.H"
#include "cut_edge.H"
#include "cut_tolerance.H"

namespace GEO
{
namespace CUT
{
  class Facet;
  class LineSegment;
  class BoundaryCell;
  class LineSegmentList;
  class Creator;

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
    for ( std::vector<Edge*>::const_iterator i=edges.begin(); i!=edges.end(); ++i )
    {
      Edge * e = *i;
      e->Register( this );
    }
  }

  virtual ~Side() {}

  int Id() const { return sid_; }

  // debug only
  void SetId( int sid ) { sid_ = sid; }

  bool IsCutSide() { return sid_ > -1; }

  void Register( Element * element ) { elements_.insert( element ); }

  void CollectElements( std::set<Element*> & elements )
  {
    std::copy( elements_.begin(), elements_.end(), std::inserter( elements, elements.begin() ) );
  }

  virtual DRT::Element::DiscretizationType Shape() = 0;

  virtual const CellTopologyData * Topology() = 0;

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points ) = 0;

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges ) = 0;

  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

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

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

//   virtual void MakeSideCutFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, std::set<Facet*> & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, const std::vector<Point*> & points, std::set<Facet*> & facets );

  virtual bool IsCut();

//   virtual bool DoTriangulation() { return true; }

  //bool FullSideCut() { return cut_lines_.size()==edges_.size() and facets_.size()==1; }

  bool OnSide( const std::set<Point*, PointPidLess> & points );

  bool OnEdge( Point * point );

  bool OnEdge( Line * line );

  bool HaveCommonEdge( Side & side );

//   virtual void ExchangeFacetSide( Side * side, bool cutsurface ) = 0;

  virtual void AddLine( Line* cut_line );

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points );

  virtual bool FindAmbiguousCutLines( Mesh & mesh, Element * element, Side & side, const std::set<Point*> & cut );

  virtual void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells );

  void CreateMissingLines( Creator & creator, Element * element );

  void Print();

  Node * OnNode( const LINALG::Matrix<3,1> & x );

  const std::vector<Edge*> & Edges() const { return edges_; }

  Edge * FindEdge( Point * begin, Point * end );

  const std::vector<Node*> & Nodes() const { return nodes_; }

  bool FindCutPoints( Mesh & mesh, Element * element, Side & other );

  bool FindCutLines( Mesh & mesh, Element * element, Side & other );

  void GetCutPoints( Element * element, Side & other, std::set<Point*> & cuts );

  const std::vector<Line*> & CutLines() const { return cut_lines_; }

  const std::set<Element*> & Elements() const { return elements_; }

  const std::vector<Facet*> & Facets() const { return facets_; }

protected:

  bool AllOnNodes( const std::set<Point*> & points );

private:

  int sid_;

  std::vector<Node*> nodes_;

  std::vector<Edge*> edges_;

  std::set<Element*> elements_;

  std::vector<Line*> cut_lines_;

  std::vector<Facet*> facets_;
};

  inline std::ostream & operator<<( std::ostream & stream, Side & s )
  {
    stream << "side: {";
    const std::vector<Node*> & nodes = s.Nodes();
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      n->point()->Print( stream );
      stream << ",";
    }
    stream << "}";
    return stream;
  }

template <>
class ConcreteSide<DRT::Element::tri3> : public Side
{
public:

  //const int numNodesSurface = 3;

  ConcreteSide<DRT::Element::tri3>( int sid,
                                    const std::vector<Node*> & nodes,
                                    const std::vector<Edge*> & edges )
    : Side( sid, nodes, edges )
  {}

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::tri3;
  }

  virtual const CellTopologyData * Topology()
  {
    return shards::getCellTopologyData< shards::Triangle<3> >();
  }

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  {
    return edge.Cut( mesh, *this, cut_points );
  }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  {
    if ( fabs( s     ) < MINIMALTOL ) edges.push_back( Edges()[0] );
    if ( fabs( r+s-1 ) < MINIMALTOL ) edges.push_back( Edges()[1] );
    if ( fabs( r     ) < MINIMALTOL ) edges.push_back( Edges()[2] );
  }

  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 3> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


template <>
class ConcreteSide<DRT::Element::quad4> : public Side
{
public:

  //const int numNodesSurface = 4;

  ConcreteSide<DRT::Element::quad4>( int sid,
                                     const std::vector<Node*> & nodes,
                                     const std::vector<Edge*> & edges )
    : Side( sid, nodes, edges )
  {}

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::quad4;
  }

  virtual const CellTopologyData * Topology()
  {
    return shards::getCellTopologyData< shards::Quadrilateral<4> >();
  }

  virtual void Cut( Mesh & mesh, Edge & edge, std::set<Point*, PointPidLess> & cut_points )
  {
    return edge.Cut( mesh, *this, cut_points );
  }

  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  {
    if ( fabs( s+1 ) < MINIMALTOL ) edges.push_back( Edges()[0] );
    if ( fabs( r-1 ) < MINIMALTOL ) edges.push_back( Edges()[1] );
    if ( fabs( s-1 ) < MINIMALTOL ) edges.push_back( Edges()[2] );
    if ( fabs( r+1 ) < MINIMALTOL ) edges.push_back( Edges()[3] );
  }

  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 4> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }
};


}
}

#endif
