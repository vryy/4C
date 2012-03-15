#ifndef CUT_POINTPOOL_H
#define CUT_POINTPOOL_H

#include <set>

#include <Teuchos_RCP.hpp>

#include "cut_boundingbox.H"
#include "cut_point.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Edge;
  class Side;

/// internal point pool class
/*!
 * Might need some speedup.
 */
class OctTreeNode
{
public:

  OctTreeNode( double norm ) : norm_( norm ) {}

  /*!
  \brief If a point with the coordinates "x" does not exists, it creates a new point correspondingly
   */
  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance );

  /*!
  \brief Get the point with the specified coordinates "x" from the pointpool
   */
  Point* GetPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance );

  unsigned size() const { return points_.size(); }

  bool IsLeaf() const { return nodes_[0]==Teuchos::null; }

  /*!
  \brief Simply insert p into the pointpool and correspondingly modify the boundingbox size
   */
  void AddPoint( const double * x, Teuchos::RCP<Point> p );

  void CollectEdges( const BoundingBox & edgebox, plain_edge_set & edges );

  void CollectSides( const BoundingBox & sidebox, plain_side_set & sides );

  void CollectElements( const BoundingBox & sidebox, plain_element_set & elements );

  void ResetOutsidePoints();

  void Print( int level, std::ostream & stream );

private:

  /*!
  \brief Create a point with the specified ID
   */
  Teuchos::RCP<Point> CreatePoint( unsigned newid, const double * x, Edge * cut_edge, Side * cut_side );

  OctTreeNode* Leaf( const double * x );

  void Split( int level );

  Teuchos::RCP<OctTreeNode> nodes_[8];

  RCPPointSet points_;

  BoundingBox bb_;

  LINALG::Matrix<3,1> splitpoint_;

  double norm_;
};

/*!
\brief Oct tree based storage container
 */
class PointPool
{
public:

  explicit PointPool( double norm=1 ) : tree_( norm ) {}

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
  { return tree_.NewPoint( x, cut_edge, cut_side, tolerance ); }

  Point* GetPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
  { return tree_.GetPoint( x, cut_edge, cut_side, tolerance ); }

  unsigned size() const { return tree_.size(); }

  /*!
  \brief Collects the edges that are fully or partially contained in the edgebox
   */
  void CollectEdges( const BoundingBox & edgebox, plain_edge_set & edges )
  { tree_.CollectEdges( edgebox, edges ); }

  /*!
  \brief Collects the sides that are fully or partially contained in the sidebox
   */
  void CollectSides( const BoundingBox & sidebox, plain_side_set & sides )
  { tree_.CollectSides( sidebox, sides ); }

  /*!
  \brief Collects the elements that are fully or partially contained in the sidebox
   */
  void CollectElements( const BoundingBox & sidebox, plain_element_set & elements )
  { tree_.CollectElements( sidebox, elements ); }

  void ResetOutsidePoints()
  { tree_.ResetOutsidePoints(); }

  void Print( std::ostream & stream )
  { tree_.Print( 0, stream ); }

private:

  OctTreeNode tree_;
};

}
}

#endif
