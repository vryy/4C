#ifdef STKADAPTIVE

#ifndef REFINE_HPP
#define REFINE_HPP

#include "stk_types.H"

namespace STK
{

  class Mesh;

/// Support data for mesh refinement algorithm
struct CreationData
{
  CreationData()
    : node_( NULL ),
      constraint_( NULL ),
      needconstraint_( false )
  {}

  // hanging node to be created
  stk::mesh::Entity* node_;
  stk::mesh::Entity* constraint_;

  bool needconstraint_;

  // participating processors
  // My rank is always included. If no communication is required only my rank.
  std::set<unsigned> procs_;

  EntitySet neighborElements_;

  bool new_hanging_node()
  {
    if ( node_==NULL )
    {
      if ( constraint_!=NULL )
        throw std::logic_error( "can not have constraint on new hanging node" );
      return true;
    }
    else
    {
      return false;
    }
  }

  unsigned owner_rank()
  {
    if ( procs_.size()==0 )
      throw std::runtime_error( "empty process list" );
    return *procs_.begin();
  }
};


/// parallel mesh refinement algorithm
/*!
  The Refine() method is called with a list of elements (actually entity keys)
  that are to be refined. All selected elements will be refined. The mesh
  state is preserved by greedy refinement. Neighboring elements are included
  in the refinement set if required.

  \author u.kue
 */
class RefineSet
{
public:

  explicit RefineSet(Mesh* mesh);

  void Refine(const std::vector<stk::mesh::EntityKey>& elements);

  void CreateHangingNode(const stk::mesh::Entity* e,
                         const std::vector<const stk::mesh::Entity*>& nodes);

  int AssignNodalIds(const std::vector<int>* nid);

  stk::mesh::EntityId NodeId(const stk::mesh::Entity* n1,
                             const stk::mesh::Entity* n2) const;

  stk::mesh::EntityId NodeId(const stk::mesh::Entity* n1,
                             const stk::mesh::Entity* n2,
                             const stk::mesh::Entity* n3,
                             const stk::mesh::Entity* n4) const;

  stk::mesh::EntityId NodeId(const stk::mesh::Entity* n1,
                             const stk::mesh::Entity* n2,
                             const stk::mesh::Entity* n3,
                             const stk::mesh::Entity* n4,
                             const stk::mesh::Entity* n5,
                             const stk::mesh::Entity* n6,
                             const stk::mesh::Entity* n7,
                             const stk::mesh::Entity* n8) const;

private:

  void RefineElement( stk::mesh::BulkData & bulk,
                      stk::mesh::Entity * e,
                      EntityVector::iterator & elements_iter,
                      EntitySet & unshared_nodes );

  void RefineTransitive( stk::mesh::Entity* e,
                         EntitySet& refine_elements,
                         EntitySet& remote_refine );

  void SyncronizeNodeCreation();

  void CommunicateDeactivation( const EntitySet & refine_elements );

  void CountNewEntities( const EntitySet & elements,
                         std::vector<std::size_t> & requests );
  void ConnectNewEntities( EntitySet & elements,
                           std::vector<std::size_t> & requests,
                           EntityVector & entities );

  bool FaceNeedsRefinement( stk::mesh::Entity * e );

  Mesh* mesh_;

  /// face node ids to hanging node and processor rank information
  std::map<EntityIdSet,CreationData> creation_;
};

}

#endif
#endif
