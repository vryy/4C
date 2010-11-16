#ifndef CUT_EDGE_H
#define CUT_EDGE_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "cut_node.H"

namespace GEO
{
namespace CUT
{
  class Mesh;
  class LinearElement;
  class LinearSide;
  class LevelSetSide;

template <DRT::Element::DiscretizationType edgetype>
class ConcreteEdge;

template <>
class ConcreteEdge<DRT::Element::line2>;

template <>
class ConcreteEdge<DRT::Element::line3>;


template <DRT::Element::DiscretizationType sidetype>
class ConcreteSide;

template <>
class ConcreteSide<DRT::Element::tri3>;

template <>
class ConcreteSide<DRT::Element::tri6>;

template <>
class ConcreteSide<DRT::Element::quad4>;

template <>
class ConcreteSide<DRT::Element::quad8>;

template <>
class ConcreteSide<DRT::Element::quad9>;


template <DRT::Element::DiscretizationType elementtype>
class ConcreteElement;

template <>
class ConcreteElement<DRT::Element::hex8>;

template <>
class ConcreteElement<DRT::Element::hex20>;

template <>
class ConcreteElement<DRT::Element::hex27>;

template <>
class ConcreteElement<DRT::Element::tet4>;

template <>
class ConcreteElement<DRT::Element::tet10>;

template <>
class ConcreteElement<DRT::Element::wedge6>;

template <>
class ConcreteElement<DRT::Element::wedge15>;

template <>
class ConcreteElement<DRT::Element::pyramid5>;


class Edge
{
public:

  explicit Edge( const std::vector<Node*> & nodes )
    : nodes_( nodes )
  {
  }

  virtual ~Edge() {}

  void Register( Side* side )
  {
    sides_.insert( side );
    for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
    {
      Node * n = *i;
      n->Register( side );
    }
  }

  bool AtSide( Side* side ) { return sides_.count( side )>0; }

  bool Matches( Point * begin, Point * end )
  { return ( ( BeginNode()->point()==begin and EndNode()->point()==end ) or
             ( BeginNode()->point()==end   and EndNode()->point()==begin ) ); }

  //virtual bool HasPoint( Point * p ) = 0;

  const std::set<Side*> & Sides() { return sides_; }

  virtual void FillComplete( Mesh & mesh ) = 0;

  virtual void Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, std::set<Point*, PointPidLess> & cuts ) = 0;

  virtual void Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, std::set<Point*, PointPidLess> & cuts ) = 0;

  virtual void LevelSetCut( Mesh & mesh, LevelSetSide & side, std::set<Point*, PointPidLess> & cuts ) = 0;

  virtual void AddPoint( Point* cut_point ) = 0;

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

  void Print()
  {
    nodes_[0]->Print();
    for ( unsigned i=1; i<nodes_.size(); ++i )
    {
      std::cout << "--";
      nodes_[i]->Print();
    }
  }

  void Plot( std::ofstream & f )
  {
    f << "# edge\n";
    BeginNode()->Plot( f );
    if ( nodes_.size()==3 )
      MiddleNode()->Plot( f );
    EndNode()->Plot( f );
    f << "\n\n";
  }

  virtual void CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points ) = 0;

  virtual void CutPoints( Side * side, std::set<Point*, PointPidLess> & cut_points ) = 0;

  virtual void CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line ) = 0;

  virtual bool IsCut( Side * side ) = 0;

  Node* BeginNode() { return nodes_.front(); }

  Node* MiddleNode()
  {
    if ( nodes_.size()!=3 )
      throw std::runtime_error( "middle node in line3 only" );
    return nodes_[2];
  }

  Node* EndNode() { return nodes_[1]; }

  Point* NodeInElement( Element * element, Point * other );

  const std::vector<Node*> & Nodes() const { return nodes_; }

private:

  std::vector<Node*> nodes_;

  std::set<Side*> sides_;
};


template <>
class ConcreteEdge<DRT::Element::line2> : public Edge
{
public:

  //const int numNodesLine = 2;

  explicit ConcreteEdge<DRT::Element::line2>( const std::vector<Node*> & nodes )
    : Edge( nodes ),
      cut_points_( PointPositionLess( this ) )
  {
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      n->Register( this );
      n->point()->AddEdge( this );
    }
  }

  virtual void FillComplete( Mesh & mesh );

  //virtual bool HasPoint( Point * p ) { return p->IsCut( this ) and cut_points_.count( p )>0; }

  virtual void Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, std::set<Point*, PointPidLess> & cuts );
  virtual void Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, std::set<Point*, PointPidLess> & cuts );
  virtual void LevelSetCut( Mesh & mesh, LevelSetSide & side, std::set<Point*, PointPidLess> & cuts );

  virtual void AddPoint( Point* cut_point );

  virtual void CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points );

  virtual void CutPoints( Side * side, std::set<Point*, PointPidLess> & cut_points );

  virtual void CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line );

  virtual bool IsCut( Side * side );

  void Coordinates( LINALG::Matrix<3, 2> & xyze_lineElement )
  { Edge::Coordinates( xyze_lineElement.A() ); }

  bool FindCutPoints( Mesh & mesh, LinearElement * element, LinearSide & side, LinearSide & other );

  void GetCutPoints( Element * element, Side & side, Side & other, std::set<Point*> & cuts );

private:

  std::set<Point*, PointPositionLess> cut_points_;

};


template <>
class ConcreteEdge<DRT::Element::line3> : public Edge
{
public:

  //const int numNodesLine = 3;

  explicit ConcreteEdge<DRT::Element::line3>( const std::vector<Node*> & nodes )
    : Edge( nodes ),
      subedge1_( NULL ),
      subedge2_( NULL )
  {
  }

  virtual void FillComplete( Mesh & mesh );

  //virtual bool HasPoint( Point * p );

  virtual void Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, std::set<Point*, PointPidLess> & cuts );
  virtual void Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, std::set<Point*, PointPidLess> & cuts );
  virtual void LevelSetCut( Mesh & mesh, LevelSetSide & side, std::set<Point*, PointPidLess> & cuts );

  virtual void AddPoint( Point* cut_point );

  virtual void CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points );

  virtual void CutPoints( Side * side, std::set<Point*, PointPidLess> & cut_points );

  virtual void CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line );

  virtual bool IsCut( Side * side );

  void Coordinates( LINALG::Matrix<3, 3> & xyze_lineElement )
  { Edge::Coordinates( xyze_lineElement.A() ); }

private:

  Edge* subedge1_;
  Edge* subedge2_;
};


}
}

#endif
