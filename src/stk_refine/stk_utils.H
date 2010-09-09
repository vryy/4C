#ifdef STKADAPTIVE

#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <set>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include "stk_mesh.H"
#include "stk_types.H"

namespace STK
{

void FindNeighborElements(const std::vector<const stk::mesh::Entity*>& nodes,
                          EntitySet & neighbors);
stk::mesh::Entity* FindConstraint( const std::vector<const stk::mesh::Entity*>& nodes );
void FindNeighbors( const std::vector<const stk::mesh::Entity*>& nodes,
                    EntitySet & neighbors,
                    unsigned entity_rank );
void FindNeighbors( stk::mesh::Entity * e,
                    EntitySet & neighbors,
                    unsigned entity_rank );

void FindCommonEntities( stk::mesh::Entity* e, EntitySet & entities, stk::mesh::EntityRank rank );
void FindCommonEntities( const EntitySet& nodes, EntitySet & entities, stk::mesh::EntityRank rank );

/// Find all processor ranks all nodes of this element share
void FindCommonNodeSharing(stk::mesh::Entity* element, std::set<unsigned> & common_sharing);

/// Find all processor ranks all nodes share
void FindCommonNodeSharing(const std::vector<const stk::mesh::Entity*> & nodes,
                           std::set<unsigned> & common_sharing);

/// Find inactive parent element that is connected to all nodes of its children.
stk::mesh::Entity* InactiveParent( stk::mesh::Part & active,
                                   const EntitySet & nodes,
                                   stk::mesh::EntityRank parent_rank );

void FindCommonParts( const std::vector<stk::mesh::Entity*> & nodes,
                      stk::mesh::PartVector & parts );

bool IsChildListComplete( stk::mesh::Entity* e,
                          const EntitySet& children );

bool is_child( stk::mesh::Entity * child, stk::mesh::Entity * parent );

stk::mesh::Entity & AddNode( Mesh & mesh, int nid, double* x );

void MiddlePoint( Mesh & mesh, stk::mesh::Entity* nid1, stk::mesh::Entity* nid2, double* x );

}

#endif
#endif
