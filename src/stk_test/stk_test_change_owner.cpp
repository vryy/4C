
#include <iostream>
#include <fstream>
#include <iomanip>

#include <mpi.h>
#include <stk_util/parallel/Parallel.hpp>

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

#include "stk_test_utils.H"

void test_change_owner( bool & success )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const int p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_rank_names() );

  stk::mesh::VectorField * coordinates_field =
    & stk::mesh::declare_vector_field_on_all_nodes( meta_data, "coordinates", 3 );

  stk::mesh::Part * block_quad = & meta_data.declare_part("quad" , stk::mesh::Element );

  stk::mesh::set_cell_topology< shards::Quadrilateral<4> >( *block_quad );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, pm, 1000 );
  bulk_data.modification_begin();

  int rows = 2;
  int cols = 2;

  int numnode = 0;
  int numele  = 0;
  if ( p_rank==0 )
  {
    numnode = (rows+1)*(cols+1);
    numele  = rows*cols;
  }

  std::vector<std::size_t> requests( stk::mesh::EntityRankEnd, 0 );
  requests[stk::mesh::Node] = numnode;
  requests[stk::mesh::Element] = numele;
  std::vector<stk::mesh::Entity *> requested_entities;
  bulk_data.generate_new_entities( requests, requested_entities );

  if ( p_rank==0 )
  {
    for (int i=0; i<rows+1; ++i)
    {
      for (int j=0; j<cols+1; ++j)
      {
        stk::mesh::Entity * n = requested_entities[ i+j*(rows+1) ];
        double * const coord = stk::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*i;
        coord[1] = .1*j;
        coord[2] = 0;
      }
    }

    std::vector<stk::mesh::Part*> add_parts    = std::vector<stk::mesh::Part*>();
    std::vector<stk::mesh::Part*> remove_parts = std::vector<stk::mesh::Part*>();
    add_parts.push_back( block_quad );

    for (int i=0; i<rows; ++i)
    {
      for (int j=0; j<cols; ++j)
      {
        stk::mesh::Entity * e = requested_entities[ numnode + i*cols+j ];
        bulk_data.change_entity_parts( *e, add_parts, remove_parts );
        bulk_data.declare_relation( *e , *requested_entities[ i+j*(rows+1) ] , 0 );
        bulk_data.declare_relation( *e , *requested_entities[ i+j*(rows+1)+1 ] , 1 );
        bulk_data.declare_relation( *e , *requested_entities[ i+(j+1)*(rows+1)+1 ] , 2 );
        bulk_data.declare_relation( *e , *requested_entities[ i+(j+1)*(rows+1) ] , 3 );
      }
    }
  }

  bulk_data.modification_end();

  if ( p_size==2 )
  {
    std::vector<stk::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 2 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 5 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 6 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 8 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 9 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 4 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    ep.clear();

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 1 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 4 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 7 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 1 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 2 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();
  }
}


STK_TEST_MAIN( test_change_owner )

