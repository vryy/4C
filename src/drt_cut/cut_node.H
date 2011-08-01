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

/// One node in a cut mesh
/*!
 * Nodes have a mesh-unique id, a levelset-value and a point. The point
 * contains the major part of the information.
 *
 * Furthermore, nodes have some idea how many dofsets are required. This is
 * optional information than can be used by the xfem code.
 */
class Node
{
public:

  Node( int nid, Point * point, double lsv ) : nid_( nid ), point_( point ), lsv_( lsv ) {}

  int Id() const { return nid_; }

  void Register( Edge* edge ) { edges_.insert( edge ); point_->AddEdge( edge ); }

  void Register( Side* side ) { point_->AddSide( side ); }

  void Register( Element* element ) { point_->AddElement( element ); }

  void RegisterCuts();

  const plain_edge_set & Edges() { return edges_; }

  void Coordinates( double * x )
  { point_->Coordinates( x ); }

  Point::PointPosition Position() const { return point_->Position(); }

  void Print( std::ostream & f=std::cout )
  {
    point_->Print();
  }

  void Plot( std::ostream & f )
  {
    point_->Plot( f );
  }

  Point * point() { return point_; }

  double LSV() const { return lsv_; }

  const plain_element_set & Elements() const { return point_->Elements(); }

  /// Find the number of dofsets required at this node.
  void FindDOFSets( bool include_inner );

  int DofSetNumber( VolumeCell * cell );

  const std::vector<plain_volumecell_set> & DofSets() const { return dofsets_; }

  int NumDofSets() const { return DofSets().size(); }

private:

  void BuildDOFCellSets( Point * p,
                         const plain_volumecell_set & cells,
                         const plain_volumecell_set & nodal_cells,
                         plain_volumecell_set & done );

  int nid_;
  Point* point_;
  double lsv_;

  plain_edge_set edges_;

  std::vector<plain_volumecell_set> dofsets_;
};

}
}

#endif
