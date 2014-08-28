#ifndef CUT_NODE_H
#define CUT_NODE_H

#include "cut_point.H"

#include "../drt_fluid_xfluid/xfluid_defines.H"

namespace GEO
{
namespace CUT
{
class VolumeCell;

/*!
  \brief internal comparator class used for sorting sets of plain_volumecell_sets and plain_volumecell_sets.
         This sorting is necessary to get the
         same order of sets, important for restarts, restarts of Newton schemes in monolithic solves but also for debugging
 */
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

  /// constructor
  Node( int nid, Point * point, double lsv ) : nid_( nid ), point_( point ), lsv_( lsv ), first_set_is_std_set_(false) { selfcutposition_ = Point::undecided; }

  /// get node's node Id
  int Id() const { return nid_; }

  /// register an edge adjacent to this node
  void Register( Edge* edge ) { edges_.insert( edge ); point_->AddEdge( edge ); }

  /// register a side adjacent to this node
  void Register( Side* side ) { point_->AddSide( side ); }

  /// register an element adjacent to this node
  void Register( Element* element ) { point_->AddElement( element ); }

  /// register cuts
  void RegisterCuts();


  /*========================================================================*/
  //! @name get routines for adjacent objects
  /*========================================================================*/

  /// get edges adjacent to this node
  const plain_edge_set & Edges() { return edges_; }

  /// Get the coordinates of the node from its point's information
  void Coordinates( double * x ) const { point_->Coordinates( x ); }

  /// Get the position of the node whether it is in fluid, structure or on the cut face
  Point::PointPosition Position() const { return point_->Position(); }

  /// Returns the point that defines the node
  Point * point() const { return point_; }

  /// Returns the levelset value at this node (if it is a levelset node)
  double LSV() const { return lsv_; }

  /// Returns sides that are connected at this node
  const plain_side_set & Sides () const { return point_->CutSides(); }

  /// get all elements adjacent to this node's point
  const plain_element_set & Elements() const { return point_->Elements(); }


  /*========================================================================*/
  //! @name print and plot routines
  /*========================================================================*/

  /// print node's or its point's information to the stream
  void Print( std::ostream & f=std::cout ){ point_->Print(); }

  /// plot node's or its point's information to the stream
  void Plot( std::ostream & f ){ point_->Plot( f, this->Id() ); };


  /*========================================================================*/
  //! @name DofSet-management
  /*========================================================================*/

  /// Assign the vc_sets to the node if possible
  void AssignNodalCellSet( std::vector<plain_volumecell_set> & ele_vc_sets,
      std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets);

  /// Find the dofsets required at this node. (old unused version)
  void FindDOFSets( bool include_inner );

  /// Find the dofsets required at this node.
  void FindDOFSetsNEW( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, std::vector<plain_volumecell_set> & cell_sets );

  /// get the dofset number of the Volumecell w.r.t this node (old unused version)
  int DofSetNumber( VolumeCell * cell );

  /// get the dofset number of the Volumecell w.r.t this node
  int DofSetNumberNEW( plain_volumecell_set & cells );

  /// return the sets of volumecells (old unsed version)
  const std::vector<plain_volumecell_set> & DofSets() const { return dofsets_; }

  /// return a vector of sets (corresponding to the dofsets) (connected via subelements) of Volumecells
  const std::vector<std::set<plain_volumecell_set, Cmp> > & DofCellSets() const { return dof_cellsets_; }

  /// sort all dofsets via xyz point coordinates (use compare functions in cut_node.H)
  void SortDOFCellSets();

  /// get the number of dofsets at this node
  int NumDofSets() const { return DofCellSets().size(); }

  /// old variant, unused: get the number of dofsets at this node
  // int NumDofSets() const { return DofSets().size();}

  /// Returns the selfcutposition of this node
  Point::PointPosition SelfCutPosition(){ return selfcutposition_; }

  /// Gives this node a selfcutposition and spread the positional information
  void SelfCutPosition( Point::PointPosition p );

  /// Erase the cutsideedge from this node because it is deleted in the selfcut
  void EraseCutSideEdge( Edge * cutsideedge ){ edges_.erase( cutsideedge ); }

  /// Return true is the given node has same position coordinates as this node
  bool isAtSameLocation ( const Node* nod1 ) const;



private:

  /// build sets of connected volumecells in a 1-ring around the node
  void BuildDOFCellSets( Point * p,
      const std::vector<plain_volumecell_set> & cell_sets,
      const plain_volumecell_set & cells,
      const std::vector<plain_volumecell_set> & nodal_cells,
      plain_volumecell_set & done,
      bool isnodalcellset = false);

  /// build sets of connected volumecells in a 1-ring around the node (old unused version)
  void BuildDOFCellSets( Point * p,
      const plain_volumecell_set & cells,
      const plain_volumecell_set & nodal_cells,
      plain_volumecell_set & done );


  /*========================================================================*/
  //! @name private class variables
  /*========================================================================*/

  int nid_;                ///< node's id (the same as in DRT::Discretization)
  Point* point_;           ///< pointer to the GEO::CUT::Point which defines the node
  double lsv_;             ///< levelset value if it is a levelset cut

  plain_edge_set edges_;   ///< all adjacent edges

  std::vector<plain_volumecell_set> dofsets_; ///< set of Volumecells (unused old version)

  /// stores for each node's dofset all its Volumecells (sets of Volumecells) which has been connected in a 1-ring of adjacent elements around the node
  std::vector<std::set<plain_volumecell_set, Cmp> > dof_cellsets_;

  /*!
  \brief boolian if the first set in dof_cellsets_ is a std set or not.
         If it is, then the first set must not be changed during sorting
         elements without elementhandle assume the first set as std set!
         Initialized to false and set to true in case that there is a standard dofset around this node
         REMARK: don't change the initialization to true!
   */
  bool first_set_is_std_set_;

  /// every cutsidenode knows its selfcutposition
  Point::PointPosition selfcutposition_;
};


}
}

#endif
