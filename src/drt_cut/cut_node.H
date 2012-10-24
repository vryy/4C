#ifndef CUT_NODE_H
#define CUT_NODE_H

#include "cut_point.H"

#include "../drt_fluid/xfluid_defines.H"

namespace GEO
{
namespace CUT
{
  class VolumeCell;

  class Cmp
  {
  public:
    /// compare operator for plain_volumecell_sets
    bool operator()( plain_volumecell_set s1, plain_volumecell_set s2 );

    /// compare operator for sets of plain_volumecell_sets
    bool operator()( const std::set<plain_volumecell_set, Cmp>& set1, const std::set<plain_volumecell_set, Cmp>& set2);


  private:

    /// compare routine for sets of plain_volumecell_set
    bool Compare(const std::set<plain_volumecell_set, Cmp>& set1, const std::set<plain_volumecell_set, Cmp>& set2 );

    /// compare routine for plain_volumecell_sets
    bool Compare(const plain_volumecell_set& s1, const plain_volumecell_set& s2 );

    /// compare routine for volumecells based on point ids
    bool Compare( VolumeCell* vc1, VolumeCell* vc2 );


  };

/// One node in a cut mesh
/*!
 \brief Class to deal with nodes in the cut mesh. Nodes have a mesh-unique id, a levelset-value and a point. The point
 contains the major part of the information.

  Furthermore, nodes have some idea how many dofsets are required. This is
  optional information than can be used by the xfem code.
 */
class Node
{
public:

  Node( int nid, Point * point, double lsv ) : nid_( nid ), point_( point ), lsv_( lsv ), first_set_is_std_set_(true) {}

  int Id() const { return nid_; }

  void Register( Edge* edge ) { edges_.insert( edge ); point_->AddEdge( edge ); }

  void Register( Side* side ) { point_->AddSide( side ); }

  void Register( Element* element ) { point_->AddElement( element ); }

  void RegisterCuts();

  const plain_edge_set & Edges() { return edges_; }

  /*!
  \brief Get the coordinates of the node from its point
   */
  void Coordinates( double * x )
  { point_->Coordinates( x ); }

  /*!
  \brief Get the position of the node whether it is in fluid, structure or on the cut face
   */
  Point::PointPosition Position() const { return point_->Position(); }

  void Print( std::ostream & f=std::cout )
  {
    point_->Print();
  }

  void Plot( std::ostream & f )
  {
    point_->Plot( f );
  }

  /*!
  \brief Returns the point that defines the node
  */
  Point * point() { return point_; }

  /*!
  \brief Returns the levelset value at this node
   */
  double LSV() const { return lsv_; }

  const plain_element_set & Elements() const { return point_->Elements(); }

  /// Assign the vc_sets to the node if possible
  void AssignNodalCellSet( std::vector<plain_volumecell_set> & ele_vc_sets,
                                           std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets);

  /*!
  \brief Find the number of dofsets required at this node.
   */
  void FindDOFSets( bool include_inner );

  void FindDOFSetsNEW( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, std::vector<plain_volumecell_set> & cell_sets );

  int DofSetNumber( VolumeCell * cell );

  int DofSetNumberNEW( plain_volumecell_set & cells );

  const std::vector<plain_volumecell_set> & DofSets() const { return dofsets_; }

  const std::vector<std::set<plain_volumecell_set, Cmp> > & DofCellSets() const { return dof_cellsets_; }

  // This is only for test with new dofsetmanagement
#ifdef DOFSETS_NEW
  int NumDofSets() const { return DofCellSets().size(); }
#else
  int NumDofSets() const { return DofSets().size();}
#endif

  void SortDOFCellSets();

private:

  void BuildDOFCellSets( Point * p,
                         const std::vector<plain_volumecell_set> & cell_sets,
                         const plain_volumecell_set & cells,
                         const std::vector<plain_volumecell_set> & nodal_cells,
                         plain_volumecell_set & done,
                         bool isnodalcellset = false);

  void BuildDOFCellSets( Point * p,
                         const plain_volumecell_set & cells,
                         const plain_volumecell_set & nodal_cells,
                         plain_volumecell_set & done );



  int nid_;
  Point* point_;
  double lsv_;  //levelset value

  plain_edge_set edges_;

  std::vector<plain_volumecell_set> dofsets_;

  /// stores for each set of dofs, that has to be created, all its volume_sets which has been connected over all adjacent elements
  std::vector<std::set<plain_volumecell_set, Cmp> > dof_cellsets_;

  /// boolian if the first set in dof_cellsets_ is a std set or not. If it is, then the first set must not be changed during sorting
  /// elements without elementhandle assume the first set as std set!
  bool first_set_is_std_set_;
};


}
}

#endif
