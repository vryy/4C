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

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side, bool nodalpoint );

  unsigned size() const { return points_.size(); }

  bool IsLeaf() const { return nodes_[0]==Teuchos::null; }

  void AddPoint( const double * x, Teuchos::RCP<Point> p );

  bool WithinBB( const Epetra_SerialDenseMatrix & xyz ) { return bb_.Within( xyz ); }

  void CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements );

private:

  Teuchos::RCP<Point> GetPoint( const double * x, Edge * cut_edge, Side * cut_side, bool nodalpoint );

  Teuchos::RCP<Point> CreatePoint( unsigned newid, const double * x, Edge * cut_edge, Side * cut_side, bool nodalpoint );

  OctTreeNode* Leaf( const double * x );

  void Split();

  Teuchos::RCP<OctTreeNode> nodes_[8];

  std::set<Teuchos::RCP<Point>, PointPidLess> points_;

  BoundingBox bb_;

  LINALG::Matrix<3,1> splitpoint_;
};

class PointPool
{
public:

  Point* NewPoint( const double * x, Edge * cut_edge, Side * cut_side, bool nodalpoint )
  { return tree_.NewPoint( x, cut_edge, cut_side, nodalpoint ); }

  unsigned size() const { return tree_.size(); }

  bool WithinBB( const Epetra_SerialDenseMatrix & xyz ) { return tree_.WithinBB( xyz ); }

  void CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements )
  { tree_.CollectElements( sidebox, elements ); }

private:

  OctTreeNode tree_;
};

}
}

#endif
