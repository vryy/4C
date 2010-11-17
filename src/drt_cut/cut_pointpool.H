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

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side );

  unsigned size() const { return points_.size(); }

  bool IsLeaf() const { return nodes_[0]==Teuchos::null; }

  void AddPoint( const double * x, Teuchos::RCP<Point> p );

  void CollectSides( const BoundingBox & sidebox, std::set<Side*> & sides );

  void CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements );

  void ResetOutsidePoints();

private:

  Teuchos::RCP<Point> GetPoint( const double * x, Edge * cut_edge, Side * cut_side );

  Teuchos::RCP<Point> CreatePoint( unsigned newid, const double * x, Edge * cut_edge, Side * cut_side );

  OctTreeNode* Leaf( const double * x );

  void Split( int level );

  Teuchos::RCP<OctTreeNode> nodes_[8];

  std::set<Teuchos::RCP<Point>, PointPidLess> points_;

  BoundingBox bb_;

  LINALG::Matrix<3,1> splitpoint_;
};

class PointPool
{
public:

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side )
  { return tree_.NewPoint( x, cut_edge, cut_side ); }

  unsigned size() const { return tree_.size(); }

  void CollectSides( const BoundingBox & sidebox, std::set<Side*> & sides )
  { tree_.CollectSides( sidebox, sides ); }

  void CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements )
  { tree_.CollectElements( sidebox, elements ); }

  void ResetOutsidePoints()
  { tree_.ResetOutsidePoints(); }

private:

  OctTreeNode tree_;
};

}
}

#endif
