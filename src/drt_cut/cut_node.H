#ifndef CUT_NODE_H
#define CUT_NODE_H

#include <map>
#include <set>
#include <vector>

#include "cut_point.H"

namespace GEO
{
namespace CUT
{

class Node
{
public:

  Node( int nid, Point * point ) : nid_( nid ), point_( point ) {}

  int Id() const { return nid_; }

  void Register( Edge* edge ) { edges_.insert( edge ); point_->AddEdge( edge ); }

  void Register( Side* side ) { point_->AddSide( side ); }

  const std::set<Edge*> & Edges() { return edges_; }

  void Coordinates( double * x )
  { point_->Coordinates( x ); }

  void Print()
  {
    point_->Print();
  }

  void Plot( std::ofstream & f )
  {
    point_->Plot( f );
  }

  Point * point() { return point_; }

private:

  int nid_;
  Point* point_;

  std::set<Edge*> edges_;
};

}
}

#endif
