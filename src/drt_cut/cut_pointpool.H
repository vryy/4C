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

// octtree to store points
//

class OctTreeNode
{
public:

  OctTreeNode( double norm ) : norm_( norm ) {}

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance );

  Point* GetPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance );

  unsigned size() const { return points_.size(); }

  bool IsLeaf() const { return nodes_[0]==Teuchos::null; }

  void AddPoint( const double * x, Teuchos::RCP<Point> p );

  void CollectEdges( const BoundingBox & edgebox, std::set<Edge*> & edges );

  void CollectSides( const BoundingBox & sidebox, std::set<Side*> & sides );

  void CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements );

  void ResetOutsidePoints();

  void Print( int level, std::ostream & stream );

private:

  Teuchos::RCP<Point> CreatePoint( unsigned newid, const double * x, Edge * cut_edge, Side * cut_side );

  OctTreeNode* Leaf( const double * x );

  void Split( int level );

  Teuchos::RCP<OctTreeNode> nodes_[8];

  std::set<Teuchos::RCP<Point>, PointPidLess> points_;

  BoundingBox bb_;

  LINALG::Matrix<3,1> splitpoint_;

  double norm_;
};

class PointPool
{
public:

  explicit PointPool( double norm=1 ) : tree_( norm ) {}

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
  { return tree_.NewPoint( x, cut_edge, cut_side, tolerance ); }

  Point* GetPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
  { return tree_.GetPoint( x, cut_edge, cut_side, tolerance ); }

  unsigned size() const { return tree_.size(); }

  void CollectEdges( const BoundingBox & edgebox, std::set<Edge*> & edges )
  { tree_.CollectEdges( edgebox, edges ); }

  void CollectSides( const BoundingBox & sidebox, std::set<Side*> & sides )
  { tree_.CollectSides( sidebox, sides ); }

  void CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements )
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
