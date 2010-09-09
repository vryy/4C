#ifdef STKADAPTIVE

#include "stk_comm.H"

namespace STK
{

void DistributeEntityPackStrategy::communicate_sharing_partners()
{
  // do the owner to sharing partner communication
  Distribution<NewSharingStrategy, const EntitySet>
    comm( bulk(), newly_shared );
  comm();
}


void AuraNodesToSharedNodes( stk::mesh::BulkData & bulk , const EntitySet & unshared )
{
  Distribution<AuraToSharedStrategy, const EntitySet>
    comm( bulk, unshared );
  comm();
}


void FindAuraNodes( stk::mesh::BulkData & bulk, stk::mesh::Entity * e, EntitySet & aura )
{
  stk::mesh::Part & globally_shared = bulk.mesh_meta_data().globally_shared_part();

  for ( stk::mesh::PairIterRelation iter = e->relations( stk::mesh::Node );
        not iter.empty();
        ++iter )
  {
    stk::mesh::Entity * n = iter->entity();
    if ( not stk::mesh::has_superset( n->bucket(), globally_shared ) )
    {
      aura.insert( n );
    }
  }
}

}

#endif
