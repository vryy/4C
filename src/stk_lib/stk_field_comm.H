#ifdef STKADAPTIVE

#ifndef STK_FIELD_COMM_H
#define STK_FIELD_COMM_H

#include "stk_comm.H"

namespace STK
{

class FieldDistribution : public DefaultPackStrategy
{
public:
  explicit FieldDistribution( stk::mesh::BulkData & bulk )
    : DefaultPackStrategy( bulk )
  {
  }

  void set_fields( const std::vector<stk::mesh::FieldBase*> & fields )
  {
    m_fields.reserve( fields.size() );
    m_fields.assign( fields.begin(), fields.end() );
  }

protected:

  void reserve_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
    for ( std::vector<stk::mesh::FieldBase*>::const_iterator k=m_fields.begin();
          k!=m_fields.end();
          ++k )
    {
      const stk::mesh::FieldBase & f = **k;
      unsigned field_data_size = stk::mesh::field_data_size( f, e->bucket() ) / sizeof( double );
      send_buffer.skip<double>( field_data_size );
    }
  }

  void pack_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
    for ( std::vector<stk::mesh::FieldBase*>::const_iterator k=m_fields.begin();
          k!=m_fields.end();
          ++k )
    {
      const stk::mesh::FieldBase & f = **k;
      unsigned field_data_size = stk::mesh::field_data_size( f, e->bucket() ) / sizeof( double );
      double * data = reinterpret_cast<double*>( stk::mesh::field_data( f , *e ) );
      send_buffer.pack<double>( data, field_data_size );
    }
  }

  void unpack_details( const std::vector<stk::mesh::Entity*> & entities,
                       stk::CommBuffer & recv_buffer,
                       stk::mesh::Entity* e )
  {
    for ( std::vector<stk::mesh::FieldBase*>::const_iterator k=m_fields.begin();
          k!=m_fields.end();
          ++k )
    {
      const stk::mesh::FieldBase & f = **k;
      unsigned field_data_size = stk::mesh::field_data_size( f, e->bucket() ) / sizeof( double );
      double * data = reinterpret_cast<double*>( stk::mesh::field_data( f , *e ) );
      recv_buffer.unpack<double>( data, field_data_size );
    }
  }

private:
  std::vector<stk::mesh::FieldBase*> m_fields;
};


  void CommunicateFieldGhosting( stk::mesh::BulkData & bulk, const std::vector<stk::mesh::FieldBase*> & fields );
}

#endif
#endif
