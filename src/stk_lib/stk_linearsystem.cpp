
#ifdef STKADAPTIVE

#include "stk_linearsystem.H"

namespace stk {
namespace linsys {

void add_connectivities(stk::linsys::LinearSystemInterface& ls,
                        stk::mesh::EntityRank entity_rank,
                        stk::mesh::EntityRank connected_entity_rank,
                        const std::vector<stk::mesh::FieldBase*>& fields,
                        const stk::mesh::Selector& selector,
                        const stk::mesh::BulkData& mesh_bulk)
{
  const std::vector<mesh::Bucket*>& mesh_buckets = mesh_bulk.buckets(entity_rank);
  std::vector<mesh::Bucket*> part_buckets;
  stk::mesh::get_buckets(selector, mesh_buckets, part_buckets);

  if (part_buckets.empty()) return;

  stk::mesh::Entity& first_entity = *(part_buckets[0]->begin());
  stk::mesh::PairIterRelation rel = first_entity.relations(connected_entity_rank);
  int num_connected = rel.second - rel.first;

  DofMapper& dof_mapper = ls.get_DofMapper();

  for ( unsigned j=0; j<fields.size(); ++j )
  {
    dof_mapper.add_dof_mappings(mesh_bulk, selector, connected_entity_rank, *fields[j]);
  }

  // assume dense coupling of fields

  std::vector<int> numFieldsPerID;
  std::vector<int> fieldIDs;

  numFieldsPerID.reserve( num_connected );
  fieldIDs.reserve( num_connected*fields.size() );

  for ( int i=0; i<num_connected; ++i )
  {
    numFieldsPerID.push_back( fields.size() );
    for ( unsigned j=0; j<fields.size(); ++j )
    {
      fieldIDs.push_back( dof_mapper.get_field_id( *fields[j] ) );
    }
  }

  //int field_id = dof_mapper.get_field_id(field);

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();

  int pattern_id = matgraph->definePattern(num_connected, connected_entity_rank, &numFieldsPerID[0], &fieldIDs[0]);

  int num_entities = 0;
  for(size_t i=0; i<part_buckets.size(); ++i) {
    num_entities += part_buckets[i]->size();
  }

  int block_id = matgraph->initConnectivityBlock(num_entities, pattern_id);

  std::vector<int> connected_ids(num_connected);

  for(size_t i=0; i<part_buckets.size(); ++i) {
    stk::mesh::Bucket::iterator
      b_iter = part_buckets[i]->begin(),
      b_end  = part_buckets[i]->end();
    for(; b_iter != b_end; ++b_iter) {
      stk::mesh::Entity& entity = *b_iter;
      rel = entity.relations(connected_entity_rank);
      for(int j=0; rel.first != rel.second; ++rel.first, ++j) {
        connected_ids[j] = rel.first->entity()->identifier();
      }
      int conn_id = entity.identifier();
      matgraph->initConnectivity(block_id, conn_id, &connected_ids[0]);
    }
  }
}

}
}

STK::LinearSystem::LinearSystem( Mesh & mesh, stk::mesh::FieldBase * field )
  : m_mesh( mesh ),
    m_factory( new Factory_Trilinos( mesh.parallel() ) ),
    m_ls( mesh.parallel(), m_factory )
{
  m_fields.push_back( field );
  Setup();
}


STK::LinearSystem::LinearSystem( Mesh & mesh, stk::mesh::FieldBase * field1, stk::mesh::FieldBase * field2 )
  : m_mesh( mesh ),
    m_factory( new Factory_Trilinos( mesh.parallel() ) ),
    m_ls( mesh.parallel(), m_factory )
{
  m_fields.push_back( field1 );
  m_fields.push_back( field2 );
  Setup();
}


void STK::LinearSystem::Setup()
{
  stk::mesh::Selector sel = m_mesh.OwnedPart() & m_mesh.ActivePart();

  // Create the graph in advance.
  // do this for all fields and all combinations of element type parts

  stk::linsys::add_connectivities( m_ls, stk::mesh::Element, stk::mesh::Node,
                                   m_fields, sel, m_mesh.BulkData() );

  m_ls.synchronize_mappings_and_structure();
  m_ls.create_fei_LinearSystem();

  m_dof_mapper = & m_ls.get_DofMapper();
  m_fei_ls     = m_ls.get_fei_LinearSystem();
  m_fei_matrix = m_fei_ls->getMatrix();
  m_fei_rhs    = m_fei_ls->getRHS();
  m_fei_x      = m_fei_ls->getSolutionVector();
}


void STK::LinearSystem::Commit()
{
  //m_fei_matrix->putScalar( 1. );

  m_ls.finalize_assembly();

  m_fei_matrix->writeToFile("matrix.dump");
  //m_fei_rhs->writeToFile("rhs.dump");
}

int STK::LinearSystem::Solve( int & status, Teuchos::ParameterList & params )
{
  return m_ls.solve( status, params );
}

void STK::LinearSystem::Done()
{
  stk::linsys::copy_vector_to_mesh( *m_fei_x, *m_dof_mapper, m_mesh.BulkData() );
}

void STK::LinearSystem::Assemble( stk::mesh::Bucket & bucket,
                             Intrepid::FieldContainer<double> & localStiffMatrix,
                             stk::mesh::Selector & sel,
                             stk::mesh::Part & active,
                             stk::mesh::Part & dbc )
{
  int k=0;
  for (stk::mesh::Bucket::iterator j=bucket.begin();
       j!=bucket.end();
       ++j, ++k)
  {
    int row = 0;
    for (stk::mesh::PairIterRelation r=j->relations(stk::mesh::Node);
         not r.empty();
         ++r, ++row)
    {
      stk::mesh::Entity * rn = r->entity();
      if ( sel( rn->bucket() ) )
      {
        for ( unsigned rf = 0; rf < m_fields.size(); ++rf )
        {
          stk::mesh::FieldBase & rfield = *m_fields[rf];

          // This is a scalar or vector field. Take the number of
          // entries. Assume double values. Ignore any stride.
          int r_field_data_size = stk::mesh::field_data_size( rfield, rn->bucket() ) / sizeof( double );
          for ( int r=0; r<r_field_data_size; ++r )
          {
            int rowIndex = m_dof_mapper->get_global_index( rn->entity_rank(),
                                                           rn->identifier(),
                                                           rfield,
                                                           r );
            if ( not has_superset( rn->bucket(), dbc ) )
            {
              int col = 0;
              for (stk::mesh::PairIterRelation c=j->relations(stk::mesh::Node);
                   not c.empty();
                   ++c, ++col)
              {
                stk::mesh::Entity * cn = c->entity();
                if ( has_superset( cn->bucket(), active ) )
                {
                  for ( unsigned cf = 0; cf < m_fields.size(); ++cf )
                  {
                    stk::mesh::FieldBase & cfield = *m_fields[cf];

                    // This is a scalar or vector field. Take the number of
                    // entries. Assume double values. Ignore any stride.
                    int c_field_data_size = stk::mesh::field_data_size( cfield, cn->bucket() ) / sizeof( double );
                    for ( int c=0; c<c_field_data_size; ++c )
                    {
                      int colIndex = m_dof_mapper->get_global_index( cn->entity_rank(),
                                                                     cn->identifier(),
                                                                     cfield,
                                                                     c );
                      double* val = &localStiffMatrix(k,row,col);
                      if ( m_fei_matrix->sumIn(1,&rowIndex,1,&colIndex,&val)!=0)
                      {
                        throw std::runtime_error("matrix assemble failed");
                      }
                    }
                  }
                }
              }
            }
            else
            {
              double val = 1;
              double * pval = &val;
              if ( m_fei_matrix->sumIn( 1, &rowIndex, 1, &rowIndex, &pval )!=0 )
              {
                throw std::runtime_error( "matrix assemble failed" );
              }
            }
          }
        }
      }
    }
  }
}

void STK::LinearSystem::Assemble( stk::mesh::Bucket & bucket,
                             Intrepid::FieldContainer<double> & localRHS,
                             stk::mesh::Selector & sel,
                             stk::mesh::Part & dbc )
{
  int k=0;
  for ( stk::mesh::Bucket::iterator j=bucket.begin();
        j!=bucket.end();
        ++j, ++k )
  {
    int row = 0;
    for ( stk::mesh::PairIterRelation r=j->relations(stk::mesh::Node);
          not r.empty();
          ++r, ++row )
    {
      stk::mesh::Entity * rn = r->entity();
      if ( sel( rn->bucket() ) )
      {
        for ( unsigned rf = 0; rf < m_fields.size(); ++rf )
        {
          stk::mesh::FieldBase & rfield = *m_fields[rf];

          // This is a scalar or vector field. Take the number of
          // entries. Assume double values. Ignore any stride.
          int r_field_data_size = stk::mesh::field_data_size( rfield, rn->bucket() ) / sizeof( double );
          for ( int r=0; r<r_field_data_size; ++r )
          {
            int rowIndex = m_dof_mapper->get_global_index( rn->entity_rank(),
                                                           rn->identifier(),
                                                           rfield,
                                                           r );
            if ( not has_superset( rn->bucket(), dbc ) )
            {
              double val = -localRHS( k, row );
              m_fei_rhs->sumIn( 1, &rowIndex, &val );
            }
            else
            {
#if 0
              double val = 0;
              m_fei_rhs->sumIn( 1, &rowIndex, &val );
#endif
            }
          }
        }
      }
    }
  }
}

#endif
