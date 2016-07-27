/*!-----------------------------------------------------------------------------------------------*
 \file cut_edge.H

 \brief class representing a geometrical edge

 <pre>
\level 3
\maintainer Andy Wirtz
 wirtz@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15270
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_EDGE_H
#define CUT_EDGE_H

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "cut_node.H"

namespace GEO
{
namespace CUT
{
  class Mesh;
  class Element;
  class LevelSetSide;

template <DRT::Element::DiscretizationType sidetype>
class ConcreteSide;

template <>
class ConcreteSide<DRT::Element::tri3>;

template <>
class ConcreteSide<DRT::Element::quad4>;

template <DRT::Element::DiscretizationType elementtype>
class ConcreteElement;

template <>
class ConcreteElement<DRT::Element::hex8>;

template <>
class ConcreteElement<DRT::Element::tet4>;

template <>
class ConcreteElement<DRT::Element::wedge6>;

template <>
class ConcreteElement<DRT::Element::pyramid5>;

/*!
\brief Linear edge between two nodes. The edge nodes are always cut points. There can be further cut points on the edge.
 */
class Edge
{
public:

  explicit Edge( const std::vector<Node*> & nodes )
    : nodes_( nodes ),
      cut_points_( PointPositionLess( this ) )
  {
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      n->Register( this );
      n->point()->AddEdge( this );
    }
    selfcutposition_ = Point::undecided;
  }

  virtual ~Edge() {}

  /*!
  \brief Add the side to the list of sides cut by this edge
   */
  void Register( Side* side )
  {
    sides_.insert( side );
    for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
    {
      Node * n = *i;
      n->Register( side );
    }
  }

  /*!
  \brief Check whether the edge is part of this side
   */
  bool AtSide( Side* side ) { return sides_.count( side )>0; }

  /*!
  \brief Check whether the edges defined by these two points coincide with this edge
   */
  bool Matches( Point * begin, Point * end )
  { return ( ( BeginNode()->point()==begin and EndNode()->point()==end ) or
             ( BeginNode()->point()==end   and EndNode()->point()==begin ) ); }

  /*!
  \brief In some cases, we need to check whether a cut line coincide with any edge.
  Even if they coincide the point ids are different, only coordinates match
   */
  bool MatchWithCoordinates( Point *begin, Point *end );

  //virtual bool HasPoint( Point * p ) = 0;

  /*!
  \brief Get all the sides on which this edge is a part of
   */
  const plain_side_set & Sides() { return sides_; }

  /*!
  \brief Get the intersection points of this edge with the tri3 side and store the cut points in cuts
   */
  void Cut( Mesh & mesh, ConcreteSide<DRT::Element::tri3> & side, PointSet & cuts );

  /*!
  \brief Get the intersection points of this edge with the quad4 side and store the cut points in cuts
   */
  void Cut( Mesh & mesh, ConcreteSide<DRT::Element::quad4> & side, PointSet & cuts );

  /*!
  \brief Get the intersection points of this edge with the levelset side and store the cut points in cuts
   */
  void LevelSetCut( Mesh & mesh, LevelSetSide & side, PointSet & cuts );

  /*!
  \brief Add cut_point to the list of points at which this edge cuts some other elements
   */
  void AddPoint( Point* cut_point );

  /*!
  \brief Get the coordinates of the nodes of the side
   */
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

  /*!
  \brief Print the coordinates of the nodes in screen
   */
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

  /*!
  \brief Get the cut points on the edge defined by the edge_start and edge_end nodes
   */
  void CutPoint( Node* edge_start, Node* edge_end, std::vector<Point*> & edge_points );

  /*!
  \brief Unused
   */
  void CutPoints( Side * side, PointSet & cut_points );

  /*!
  \brief Unused
   */
  void CutPointsBetween( Point* begin, Point* end, std::vector<Point*> & line );

  void CutPointsIncluding( Point* begin, Point* end, std::vector<Point*> & line );

  void CutPointsInside( Element * element, std::vector<Point*> & line );

  bool IsCut( Side * side );

  bool ComputeCut( Edge * other, double & pos, LINALG::Matrix<3,1> & x, double& tolerance );

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

  void Coordinates( LINALG::Matrix<3, 2> & xyze_lineElement )
  { Edge::Coordinates( xyze_lineElement.A() ); }

  /*!
  \brief Computes the points at which both the sides intersect
   */
  bool FindCutPoints( Mesh & mesh, Element * element, Side & side, Side & other, int recursion );

  /*!
  \brief Cut points falling on this edge that are common to the two given sides are extracted
  */
  void GetCutPoints( Element * element, Side & side, Side & other, PointSet & cuts );

  /*!
  \brief Cut points falling on this edge that are common to the given edge are extracted
  */
  void GetCutPoints( Edge * other, PointSet & cuts );

  const PointPositionSet & CutPoints() const { return cut_points_; }
  //const std::vector<Point*> & CutPoints() const { return cut_points_; }

  void RectifyCutNumerics();

  /// Returns the selfcutposition of this edge
  Point::PointPosition SelfCutPosition(){ return selfcutposition_; }

  /// Gives this edge a selfcutposition and spreads the positional information
  void SelfCutPosition( Point::PointPosition p );

  /// Changes the selfcutposition of this edge and spreads the positional information
  void ChangeSelfCutPosition( Point::PointPosition p );

  /// Erase the cutside from this edge because it is deleted in the selfcut
  void EraseCutSide( Side * cutside ){ sides_.erase( cutside ); }

  /*!
  \brief Replace the node "nod" of this edge with given node "replwith"
   */
  void replaceNode( Node* nod, Node* replwith );


private:

  std::vector<Node*> nodes_;

  plain_side_set sides_;

  //! sorted vector contains all points (end points and cut points) on this edge
  PointPositionSet cut_points_;

  //! every cutsideedge knows its selfcutposition
  Point::PointPosition selfcutposition_;

};


}
}

#endif
