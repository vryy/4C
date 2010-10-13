#ifdef STKADAPTIVE

#ifndef STK_UNREFINE_H
#define STK_UNREFINE_H

#include "stk_types.H"
#include "stk_mesh.H"

namespace STK
{

/// parallel mesh unrefinement algorithm
/*!

  The Unrefine() method is called with a list of elements (actually entity
  keys) that are to be removed and replaced by their parent elements. The
  selected elements will be unrefined only if the mesh remains in a valid
  state. E.g. there is only one hanging node allowed on each element side.

  \author u.kue
 */
class UnrefineSet
{
public:
  explicit UnrefineSet(Mesh* mesh);

  ~UnrefineSet();

  void Unrefine(const std::vector<stk::mesh::EntityKey>& elements);

private:

  struct ParentInfo
  {
    EntitySet children;
    EntitySet nodes;
    EntitySet hanging;
  };

  struct HangingNodeInfo
  {
    HangingNodeInfo()
      : needconstraint( false ), constraint( NULL ) {}

    bool needconstraint;
    stk::mesh::Entity * constraint;
    EntitySet parents;
  };

  EntitySet remove_;
  EntitySet keep_;

  std::map<stk::mesh::Entity*, ParentInfo> parents_;
  std::map<stk::mesh::Entity*, HangingNodeInfo> hanging_;

  void FindUnrefineChildren(const std::vector<stk::mesh::EntityKey> & elements,
                            std::map<stk::mesh::Entity*,std::set<unsigned> > & unrefine_elements);

  void SendUnrefineChildren(std::map<stk::mesh::Entity*,std::set<unsigned> > & unrefine_elements,
                            std::map<stk::mesh::EntityKey,EntityKeySet> & transfer_children);

  void SelectParents( std::map<stk::mesh::EntityKey,EntityKeySet> & transfer_children );

  void FindCompleteParents();

  void FindConstraints();

  void RemoveImpossibleUnrefinement();

  bool UnrefinementPossible( stk::mesh::Entity * parent, ParentInfo & pi, EntitySet & reject );

  void TestSurfaceUnrefinement( stk::mesh::Entity * parent, ParentInfo & pi, EntitySet & reject );

  void KeepParent( stk::mesh::Entity * parent );

  void CreateConstraints();

  void DoUnrefine();

//   void CommunicateChildRemoval( EntitySet & remove_children );

//   void CommunicateParentActivation( EntitySet & activate_parent );

//   void RemoveChildren( std::vector<EntityKeySet> & child_remove_list );
//   void ActiveParent( EntitySet & activate_list );
//   void RemoveUnusedNodes( EntitySet & nodes );

//   void UpdateConstraint( stk::mesh::BulkData & bulk,
//                          stk::mesh::Part & active,
//                          const EntitySet& facenodes,
//                          stk::mesh::EntityId hangingNodeId );

//   void CollectUnrefine( stk::mesh::Entity * e, EntitySet & es );

  bool NeedConstraint( stk::mesh::Entity * hn );

  Mesh* mesh_;
};
}

#endif
#endif
