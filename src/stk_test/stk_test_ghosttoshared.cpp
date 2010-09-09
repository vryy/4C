
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

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "stk_test_utils.H"


void test_ghosttoshared( bool & success )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const int p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_rank_names() );

  stk::mesh::VectorField * coordinates_field =
    & stk::mesh::declare_vector_field_on_all_nodes( meta_data, "coordinates", 3 );

  stk::mesh::Part & owned_part = meta_data.locally_owned_part();
  stk::mesh::Part & quad_part  = meta_data.declare_part( "quad" , stk::mesh::Element );

  stk::mesh::set_cell_topology< shards::Quadrilateral<4> >( quad_part );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, pm, 100 );
  bulk_data.modification_begin();

  unsigned nx = 3;
  unsigned ny = 3;

  if ( p_rank==0 )
  {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::declare_element( bulk_data , quad_part , elem , nodes );
      }
    }

    for ( unsigned iy = 0 ; iy < ny+1 ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx+1 ; ++ix ) {
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity * n = bulk_data.get_entity( stk::mesh::Node, nid );
        double * const coord = stk::mesh::field_data( *coordinates_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  if ( p_rank==0 )
    std::cout << "original mesh:\n";
  Show( bulk_data, nx, ny );

  if ( p_size>1 )
  {
    std::vector<stk::mesh::EntityProc> ep;

    if ( p_rank==0 )
    {
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node,  3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node,  4 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node,  7 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node,  8 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 11 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 12 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 15 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 16 ), 1 ) );

      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 3 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 6 ), 1 ) );
      ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 9 ), 1 ) );
    }

    bulk_data.modification_begin();
    bulk_data.change_entity_owner( ep );
    bulk_data.modification_end();

    // output to debug

    Dump( bulk_data, "result" );
    if ( p_rank==0 )
      std::cout << "distributed mesh:\n";
    Show( bulk_data, nx, ny );

    bulk_data.modification_begin();

    if ( p_rank==1 )
    {
      // create constraint

      stk::mesh::Entity * n10 = bulk_data.get_entity( stk::mesh::Node, 10 );
      stk::mesh::Entity * n11 = bulk_data.get_entity( stk::mesh::Node, 11 );
      stk::mesh::Entity * n12 = bulk_data.get_entity( stk::mesh::Node, 12 );

      stk::mesh::PartVector add;
      add.push_back( &owned_part );
      stk::mesh::Entity & c = bulk_data.declare_entity( stk::mesh::Constraint, 1, add );
      bulk_data.declare_relation( c , *n10 , 0 );
      bulk_data.declare_relation( c , *n11 , 1 );
      bulk_data.declare_relation( c , *n12 , 2 );
    }

    bulk_data.modification_end();

    Dump( bulk_data, "shared" );

    stk::mesh::Entity * n10 = bulk_data.get_entity( stk::mesh::Node, 10 );

    if ( p_rank==0 or p_rank==1 )
    {
      if ( not stk::mesh::in_shared( *n10 ) )
        throw std::runtime_error( "NODE[10] not shared" );
    }

    bulk_data.modification_begin();

    if ( p_rank==1 )
    {
      // destroy constraint

      stk::mesh::Entity * c1 = bulk_data.get_entity( stk::mesh::Constraint, 1 );

      if ( not bulk_data.destroy_entity( c1 ) )
      {
        throw std::runtime_error( "failed to destroy constraint" );
      }
    }

    bulk_data.modification_end();

    Dump( bulk_data, "ghosted" );

    if ( p_rank==0 or p_rank==1 )
    {
      if ( stk::mesh::in_shared( *n10 ) )
        throw std::runtime_error( "NODE[10] shared" );
    }
  }
}


STK_TEST_MAIN( test_ghosttoshared )
