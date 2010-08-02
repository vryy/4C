
#include <mpi.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <cstring>
#include <cstdlib>

#include <stk_util/parallel/Parallel.hpp>

#include "../stk_lib/stk_io.H"
#include "../stk_lib/stk_mesh.H"
#include "../stk_lib/stk_refine.H"
#include "../stk_lib/stk_unrefine.H"
#include "../stk_lib/stk_iterator.H"

#include "../drt_stk/stk_algebra.H"

void execute( STK::Mesh & mesh,
              bool refine,
              const std::vector<stk::mesh::EntityKey> & eids )
{
  if ( refine )
  {
    std::cout << "Refine: " << eids.size() << "\n";
    STK::RefineSet rs( & mesh );
    mesh.Modify();
    rs.Refine( eids );
    mesh.Commit();
  }
  else
  {
    std::cout << "Unrefine: " << eids.size() << "\n";
    STK::UnrefineSet rs( & mesh );
    mesh.Modify();
    rs.Unrefine( eids );
    mesh.Commit();
  }
}


void load( STK::MetaMesh & meta, STK::Mesh & mesh, const std::string & filename )
{
  std::vector< const stk::mesh::FieldBase * > fields;
  //fields.push_back( & mesh.Coordinates() );

  meta.FileSchema().assign_indices( mesh.BulkData() );

  phdmesh::exodus::FileInput input( meta.FileSchema(), mesh.BulkData(), filename, fields );

  mesh.Commit();
}


void save( STK::MetaMesh & meta, STK::Mesh & mesh, const std::string & filename,
           std::vector< const stk::mesh::FieldBase * > fields = std::vector< const stk::mesh::FieldBase * >() )
{
  meta.FileSchema().assign_indices( mesh.BulkData() );

  int flags[ stk::mesh::EntityRankEnd ] = { 1 , 0 , 0 , 1 , 0 , 0 };

  phdmesh::exodus::FileOutput output( meta.FileSchema() , mesh.BulkData() ,
                                      mesh.ActivePart() ,
                                      filename , "Copy" ,
                                      false , fields , flags );
  output.write( 0, 0 );
}


void save( STK::MetaMesh & meta, STK::Mesh & mesh, int counter,
           std::vector< const stk::mesh::FieldBase * > fields )
{
  std::stringstream str;
  str << "result-" << counter << ".exo";
  save( meta, mesh, str.str(), fields );
}


void iterate( STK::Mesh & mesh )
{
  stk::mesh::BulkData & bulk = mesh.BulkData();
  stk::mesh::Part & active = mesh.ActivePart();

  for ( STK::ElementIterator i( bulk, active );
        not i.done();
        ++i )
  {
    stk::mesh::Entity & e = i.element();

    for ( STK::SideIterator j( i ); not j.done(); ++j )
    {
    }
  }
}


void mark( stk::mesh::BulkData & bulk,
           std::vector<stk::mesh::EntityKey> & eids,
           stk::mesh::ScalarField & constrained,
           stk::mesh::ScalarField & hanging,
           stk::mesh::ScalarField & error )
{
  algebra::PutScalar( bulk, constrained, 0.0 );
  algebra::PutScalar( bulk, hanging, 0.0 );

  const std::vector<stk::mesh::Bucket*> & constraints = bulk.buckets( stk::mesh::Constraint );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=constraints.begin();
        i!=constraints.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;

    for ( stk::mesh::Bucket::iterator j=bucket.begin();
          j!=bucket.end();
          ++j )
    {
      stk::mesh::Entity & c = *j;

      stk::mesh::PairIterRelation rel = c.relations( stk::mesh::Node );

      double & hnv = * stk::mesh::field_data( hanging, *rel->entity() );
      hnv += 1;

      for ( ; not rel.empty(); ++rel )
      {
        stk::mesh::Entity & n = *rel->entity();
        double & v = * stk::mesh::field_data( constrained, n );
        v += 1;
      }
    }
  }

  algebra::PutScalar( bulk, error, 0.0 );

  for ( std::vector<stk::mesh::EntityKey>::iterator i=eids.begin();
        i!=eids.end();
        ++i )
  {
    stk::mesh::EntityKey key = *i;
    stk::mesh::Entity * e = bulk.get_entity( key );
    if ( e != NULL )
    {
      double * err = stk::mesh::field_data( error , *e );
      if ( err!=NULL )
      {
        *err = 0.5;
      }
    }
  }
}


int main( int argc, char** argv )
{
  stk::ParallelMachine pm = stk::parallel_machine_init( &argc , &argv );

  if ( argc > 1 )
  {
    std::string filename = argv[1];

    STK::MetaMesh meta;

    stk::mesh::Part & active = *meta.MetaData().get_part( "active" );

    stk::mesh::ScalarField & error =
      stk::mesh::put_field( meta.MetaData().declare_field<stk::mesh::ScalarField>( "error" ),
                            stk::mesh::Element,
                            active );

    stk::mesh::ScalarField & constrained = declare_scalar_field_on_all_nodes( meta.MetaData(), "constrained" );
    stk::mesh::ScalarField & hanging     = declare_scalar_field_on_all_nodes( meta.MetaData(), "hanging" );

    meta.Commit();

    std::vector< const stk::mesh::FieldBase * > fields;
    fields.push_back( &constrained );
    fields.push_back( &hanging );
    fields.push_back( &error );

    STK::Mesh mesh( meta, pm );
    mesh.Modify();

    load( meta, mesh, filename );

    mesh.Modify();

    std::vector<stk::mesh::Part *> empty_parts;
    std::vector<stk::mesh::Part *> add_parts;
    add_parts.push_back( & mesh.ActivePart() );

    for ( int entity_type=0; entity_type<stk::mesh::EntityRankEnd; ++entity_type )
    {
      std::vector<stk::mesh::Entity*> store;
      const std::vector<stk::mesh::Bucket*> & entities = mesh.BulkData().buckets( entity_type );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator i=entities.begin();
            i!=entities.end();
            ++i)
      {
        stk::mesh::Bucket & bucket = **i;
        for ( stk::mesh::Bucket::iterator j=bucket.begin();
              j!=bucket.end();
              ++j )
        {
          stk::mesh::Entity & e = *j;
          store.push_back( &e );
        }
      }

      for ( std::vector<stk::mesh::Entity*>::iterator j=store.begin();
            j!=store.end();
            ++j )
      {
        stk::mesh::Entity & e = **j;
        mesh.BulkData().change_entity_parts( e, add_parts, empty_parts );
      }
    }

    mesh.Commit();

    std::ifstream log( "refine.log.save" );

    std::vector<stk::mesh::EntityKey> eids;
    bool refine = true;

    int counter = 0;
    while ( log )
    {
      std::string word;
      log >> word;

      if ( word=="Refine:" )
      {
        if ( eids.size()>0 )
        {
          execute( mesh, refine, eids );
          mark( mesh.BulkData(), eids, constrained, hanging, error );
          eids.clear();
          save( meta, mesh, counter, fields );
          iterate( mesh );
        }
        counter += 1;
        refine = true;
      }
      else if ( word=="Unrefine:" )
      {
        if ( eids.size()>0 )
        {
          execute( mesh, refine, eids );
          mark( mesh.BulkData(), eids, constrained, hanging, error );
          eids.clear();
          save( meta, mesh, counter, fields );
          iterate( mesh );
        }
        counter += 1;
        refine = false;
      }
      else if ( word=="" )
      {
      }
      else
      {
        if ( word.length()<10 )
          throw std::runtime_error( "word too short" );
        if ( word.substr( 0, 8 )!="ELEMENT[" )
          throw std::runtime_error( "not an element" );
        std::string number = word.substr( 8, word.length()-1 );
        int n = std::atoi( number.c_str() );
        stk::mesh::EntityKey key( stk::mesh::Element, n );
        eids.push_back( key );
      }
    }

    if ( eids.size()>0 )
    {
      execute( mesh, refine, eids );
      mark( mesh.BulkData(), eids, constrained, hanging, error );
      eids.clear();
    }

    // debug

    stk::mesh::Entity & e = * mesh.BulkData().get_entity( stk::mesh::Element, 1667 );
    double * err = stk::mesh::field_data( error , e );
    if ( err!=NULL )
    {
      *err = 1;
    }
    else
    {
      std::cout << "element 1667 not found\n";
    }

    fields.push_back( &error );

    save( meta, mesh, "success.exo", fields );
    iterate( mesh );
  }

  stk::parallel_machine_finalize();
  return 0;
}
