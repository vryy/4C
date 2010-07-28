#ifdef STKADAPTIVE

#include "stk_field_comm.H"

void STK::CommunicateFieldGhosting( stk::mesh::BulkData & bulk, const std::vector<stk::mesh::FieldBase*> & fields )
{
  const std::vector<stk::mesh::Entity*> & entity_comm = bulk.entity_comm();

  AuraDistribution<FieldDistribution, const std::vector<stk::mesh::Entity*> >
    dist( bulk, entity_comm );
  dist.set_fields( fields );
  dist();
}

#endif
