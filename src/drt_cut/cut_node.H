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
  class VolumeCell;

class Node
{
public:

  Node( int nid, Point * point, double lsv ) : nid_( nid ), point_( point ), lsv_( lsv ) {}

  int Id() const { return nid_; }

  void Register( Edge* edge ) { edges_.insert( edge ); point_->AddEdge( edge ); }

  void Register( Side* side ) { point_->AddSide( side ); }

  void Register( Element* element ) { point_->AddElement( element ); }

  const std::set<Edge*> & Edges() { return edges_; }

  void Coordinates( double * x )
  { point_->Coordinates( x ); }

  Point::PointPosition Position() const { return point_->Position(); }

  void Print()
  {
    point_->Print();
  }

  void Plot( std::ofstream & f )
  {
    point_->Plot( f );
  }

  Point * point() { return point_; }

  double LSV() const { return lsv_; }

  const std::set<Element*> & Elements() const { return point_->Elements(); }

  void FindDOFSets();

  int DofSetNumber( VolumeCell * cell );

  const std::vector<std::set<VolumeCell*> > & DofSets() const { return dofsets_; }

  int NumDofSets( bool include_inner );

private:

  void BuildDOFCellSets( Point * p,
                         const std::set<VolumeCell*> & cells,
                         const std::set<VolumeCell*> & nodal_cells,
                         std::set<VolumeCell*> & done );

  int nid_;
  Point* point_;
  double lsv_;

  std::set<Edge*> edges_;

  std::vector<std::set<VolumeCell*> > dofsets_;
};

}
}

#endif
