#ifndef CUT_SIDE_H
#define CUT_SIDE_H

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

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
  class Cycle;

/*!
\brief Base class for dealing sides in cut algorithm
 */
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

  /*!
  \brief Returns the ID of the side
   */
  int Id() const { return sid_; }

  /*!
  \brief Set the side ID to the input value
   */
  void SetId( int sid ) { sid_ = sid; }

  /*!
  \brief Returns true if this is a cut side
   */
  bool IsCutSide() { return sid_ > -1; }

  /*!
  \brief Register this element in which the side is a part of
   */
  void Register( Element * element ) { elements_.insert( element ); }

  /*void CollectElements( plain_element_set & elements )
  {
    std::copy( elements_.begin(), elements_.end(), std::inserter( elements, elements.begin() ) );
  }*/

  /*!
  \brief is this side closer to the startpoint than the other side? set is_closer, return if check was successful
   */
  bool IsCloserSide( LINALG::Matrix<3,1>& startpoint_xyz, GEO::CUT::Side* other,  bool& is_closer);

  /*!
  \brief Returns the geometrical shape of this side
   */
  virtual DRT::Element::DiscretizationType Shape() = 0;

  /*!
  \brief Returns the topology data for the side from Shards library
   */
  virtual const CellTopologyData * Topology() = 0;

  /*!
  \brief Calculates the points at which the side is cut by this edge
   */
  virtual void Cut( Mesh & mesh, Edge & edge, PointSet & cut_points ) = 0;

  /*!
   \brief get all edges adjacent to given local coordinates
   */
  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges ) = 0;

  /*!
   \brief get the global coordinates on side at given local coordinates
   */
  virtual void PointAt( double r, double s, LINALG::Matrix<3,1> & xyz) = 0;

  /*!
   \brief get global coordinates of the center of the side
   */
  virtual void SideCenter( LINALG::Matrix<3,1> & midpoint ) = 0;

  /*!
  \brief Calculates the local coordinates (rst) with respect to the element shape from its global coordinates (xyz), return if successful
   */
  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist = false ) = 0;

  /*!
   \brief lies point with given coordinates within this side?
   */
  virtual bool WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist, const double & Tol) = 0;

  /*
   \brief compute the cut of a ray through two points with the 2D space defined by the side
   */
  virtual bool RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi) = 0;

  /*!
  \brief Calculates the normal vector with respect to the element shape at local coordinates xsi
   */
  virtual void Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) = 0;

  /*!
  \brief Returns the global coordinates of the nodes of this side
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

  virtual void MakeOwnedSideFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

//   virtual void MakeSideCutFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, plain_facet_set & facets );

  virtual void MakeInternalFacets( Mesh & mesh, Element * element, const Cycle & points, plain_facet_set & facets );

  virtual bool IsCut();

//   virtual bool DoTriangulation() { return true; }

  //bool FullSideCut() { return cut_lines_.size()==edges_.size() and facets_.size()==1; }

  bool OnSide( const PointSet & points );

  bool OnEdge( Point * point );

  bool OnEdge( Line * line );

  bool HaveCommonNode( Side & side );

  bool HaveCommonEdge( Side & side );

  Element * CommonElement( Side * other );

//   virtual void ExchangeFacetSide( Side * side, bool cutsurface ) = 0;

  virtual void AddPoint( Point * cut_point );

  virtual void AddLine( Line* cut_line );

  virtual Facet * FindFacet( const std::vector<Point*> & facet_points );

  virtual bool FindAmbiguousCutLines( Mesh & mesh, Element * element, Side & side, const PointSet & cut );

  virtual void GetBoundaryCells( plain_boundarycell_set & bcells );

  virtual void CreateMissingLines( Creator & creator, Element * element );

  void Print();

  Node * OnNode( const LINALG::Matrix<3,1> & x );

  const std::vector<Edge*> & Edges() const { return edges_; }

  /*!
  \
  brief Returns the edge of this side with given begin and end points
   */
  Edge * FindEdge( Point * begin, Point * end );

  const std::vector<Node*> & Nodes() const { return nodes_; }

  /*!
  \brief Calculate the points at which the other side intersects with this considered side
   */
  bool FindCutPoints( Mesh & mesh, Element * element, Side & other, int recursion );

  /*!
  \brief Draw cut lines between the cut points of this edge and "other"
   */
  bool FindCutLines( Mesh & mesh, Element * element, Side & other );

  /*!
  \brief Get (not calculate) all the cut points between this edge and "other"
   */
  void GetCutPoints( Element * element, Side & other, PointSet & cuts );

  /*!
  \brief Get all the cut points that are produced by this edge
   */
  const PointSet & CutPoints() const { return cut_points_; }

  const std::vector<Line*> & CutLines() const { return cut_lines_; }

  const plain_element_set & Elements() const { return elements_; }

  const std::vector<Facet*> & Facets() const { return facets_; }

protected:

  bool AllOnNodes( const PointSet & points );

private:

  int sid_;

  std::vector<Node*> nodes_;

  std::vector<Edge*> edges_;

  plain_element_set elements_;

  std::vector<Line*> cut_lines_;

  PointSet cut_points_;

  std::vector<Facet*> facets_;
};

/// tri3 side
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

  /// Returns the geometrical shape of this side
  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::tri3;
  }

  /// Returns the topology data for the side from Shards library
  virtual const CellTopologyData * Topology()
  {
    return shards::getCellTopologyData< shards::Triangle<3> >();
  }

  /// Calculates the points at which the side is cut by this edge
  virtual void Cut( Mesh & mesh, Edge & edge, PointSet & cut_points )
  {
    return edge.Cut( mesh, *this, cut_points );
  }

  /// get all edges adjacent to given local coordinates
  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  {
    if ( fabs( s     ) < MINIMALTOL ) edges.push_back( Edges()[0] );
    if ( fabs( r+s-1 ) < MINIMALTOL ) edges.push_back( Edges()[1] );
    if ( fabs( r     ) < MINIMALTOL ) edges.push_back( Edges()[2] );
  }

  /// get the global coordinates on side at given local coordinates
  virtual void PointAt( double r, double s, LINALG::Matrix<3,1> & xyz);

  /// get global coordinates of the center of the side
  virtual void SideCenter( LINALG::Matrix<3,1> & midpoint );

  ///  lies point with given coordinates within this side?
  virtual bool WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist, const double & Tol );

  /// compute the cut of a ray through two points with the 2D space defined by the side
  virtual bool RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi);

  /// Calculates the local coordinates (rst) with respect to the element shape from its global coordinates (xyz), return if successful
  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist=false );

  /// Calculates the normal vector with respect to the element shape at local coordinates xsi
  virtual void Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal );

  /// get coordinates of side
  void Coordinates( LINALG::Matrix<3, 3> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }

};


/// quad4 side
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

  /// Returns the geometrical shape of this side
  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::quad4;
  }

  /// Returns the topology data for the side from Shards library
  virtual const CellTopologyData * Topology()
  {
    return shards::getCellTopologyData< shards::Quadrilateral<4> >();
  }

  /// Calculates the points at which the side is cut by this edge
  virtual void Cut( Mesh & mesh, Edge & edge, PointSet & cut_points )
  {
    return edge.Cut( mesh, *this, cut_points );
  }

  /// get all edges adjacent to given local coordinates
  virtual void EdgeAt( double r, double s, std::vector<Edge*> & edges )
  {
    if ( fabs( s+1 ) < MINIMALTOL ) edges.push_back( Edges()[0] );
    if ( fabs( r-1 ) < MINIMALTOL ) edges.push_back( Edges()[1] );
    if ( fabs( s-1 ) < MINIMALTOL ) edges.push_back( Edges()[2] );
    if ( fabs( r+1 ) < MINIMALTOL ) edges.push_back( Edges()[3] );
  }

  /// get the global coordinates on side at given local coordinates
  virtual void PointAt( double r, double s, LINALG::Matrix<3,1> & xyz);

  /// get global coordinates of the center of the side
  virtual void SideCenter( LINALG::Matrix<3,1> & midpoint );

  ///  lies point with given coordinates within this side?
  virtual bool WithinSide( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs, double & dist, const double & Tol );

  /// compute the cut of a ray through two points with the 2D space defined by the side
  virtual bool RayCut( const LINALG::Matrix<3,1> & p1_xyz, const LINALG::Matrix<3,1> & p2_xyz, LINALG::Matrix<2,1> & rs, double & line_xi);

  /// Calculates the local coordinates (rst) with respect to the element shape from its global coordinates (xyz), return if successful
  virtual bool LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst, bool allow_dist = false );

  /// Calculates the normal vector with respect to the element shape at local coordinates xsi
  virtual void Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal );

  /// get coordinates of side
  void Coordinates( LINALG::Matrix<3, 4> & xyze_surfaceElement )
  { Side::Coordinates( xyze_surfaceElement.A() ); }

};


}
}

std::ostream & operator<<( std::ostream & stream, GEO::CUT::Side & s );

#endif
