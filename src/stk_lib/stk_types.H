#ifdef STKADAPTIVE

#ifndef TYPES_HPP
#define TYPES_HPP

#include <set>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

typedef std::vector<stk::mesh::Entity*> EntityVector;

typedef std::set<stk::mesh::Entity*,   stk::mesh::EntityLess> EntitySet;
typedef std::set<stk::mesh::EntityKey>                        EntityKeySet;
typedef std::set<stk::mesh::EntityId>                         EntityIdSet;
typedef std::set<stk::mesh::EntityProc , stk::mesh::EntityLess> EntityProcSet;

//typedef stk::PairIter< std::vector<stk::mesh::Entity *>::iterator > PairIterEntity ;

class RefineSet;
class Mesh;

#endif
#endif
