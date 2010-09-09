#ifdef STKADAPTIVE

#include <cmath>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

// linear algebra is needed

namespace algebra
{

template <class field_type>
double Norm2( const stk::mesh::BulkData & bulk,
              field_type & dest )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type data_type;

  const stk::mesh::DataTraits & dest_traits = dest.data_traits();

  if ( not dest_traits.is_floating_point )
    throw std::runtime_error( "floating point field expected" );

  const std::vector<stk::mesh::Bucket*>& all_buckets = bulk.buckets( stk::mesh::Node );

  double norm = 0;

  // Is there a way to ask each bucket family just once?

  for ( std::size_t i=0; i<all_buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *all_buckets[i];
    int dest_field_data_size = stk::mesh::field_data_size( dest, bucket );
    if ( dest_field_data_size == 0 ) continue;

    std::size_t size = bucket.size() * dest_field_data_size/sizeof( data_type );

    data_type * dest_data = reinterpret_cast<data_type*>( bucket.field_data_location( dest ) );

    // inner loops
    // maybe use kokkos or blas here

    for ( std::size_t j=0; j<size; ++j )
    {
      norm += dest_data[j]*dest_data[j];
    }
  }

  stk::all_reduce( bulk.parallel(), stk::ReduceSum<1>( &norm ) );

  return std::sqrt( norm );
}

template <class field_type>
void PutScalar( const stk::mesh::BulkData & bulk,
                field_type & dest,
                typename stk::mesh::FieldTraits< field_type >::data_type scalar )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type data_type;

  const stk::mesh::DataTraits & dest_traits = dest.data_traits();

  if ( not dest_traits.is_floating_point )
    throw std::runtime_error( "floating point field expected" );

  for ( std::vector<stk::mesh::FieldBase::Restriction>::const_iterator
          j  = dest.restrictions().begin() ;
          j != dest.restrictions().end() ; ++j )
  {
    const std::vector<stk::mesh::Bucket*>& all_buckets = bulk.buckets( entity_rank( j->key ) );

    // Is there a way to ask each bucket family just once?

    for ( std::size_t i=0; i<all_buckets.size(); ++i )
    {
      stk::mesh::Bucket & bucket = *all_buckets[i];
      int dest_field_data_size = stk::mesh::field_data_size( dest, bucket );
      if ( dest_field_data_size == 0 ) continue;

      std::size_t size = bucket.size() * dest_field_data_size/sizeof( data_type );

      data_type * dest_data = reinterpret_cast<data_type*>( bucket.field_data_location( dest ) );

      // inner loops
      // maybe use kokkos or blas here

      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] = scalar;
      }
    }
  }
}


template <class field_type>
void PutScalar( const stk::mesh::BulkData & bulk,
                const stk::mesh::Selector & select,
                field_type & dest,
                typename stk::mesh::FieldTraits< field_type >::data_type scalar )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type data_type;

  const stk::mesh::DataTraits & dest_traits = dest.data_traits();

  if ( not dest_traits.is_floating_point )
    throw std::runtime_error( "floating point field expected" );

  for ( std::vector<stk::mesh::FieldBase::Restriction>::const_iterator
          j  = dest.restrictions().begin() ;
        j != dest.restrictions().end() ; ++j )
  {
    const std::vector<stk::mesh::Bucket*>& all_buckets = bulk.buckets( entity_rank( j->key ) );

    // Is there a way to ask each bucket family just once?

    for ( std::size_t i=0; i<all_buckets.size(); ++i )
    {
      stk::mesh::Bucket & bucket = *all_buckets[i];

      if ( select( bucket ) )
      {
        int dest_field_data_size = stk::mesh::field_data_size( dest, bucket );
        if ( dest_field_data_size == 0 ) continue;

        std::size_t size = bucket.size() * dest_field_data_size/sizeof( data_type );

        data_type * dest_data = reinterpret_cast<data_type*>( bucket.field_data_location( dest ) );

        // inner loops
        // maybe use kokkos or blas here

        for ( std::size_t j=0; j<size; ++j )
        {
          dest_data[j] = scalar;
        }
      }
    }
  }
}


template <class field_type>
void Update( const stk::mesh::BulkData & bulk,
             //const stk::mesh::Selector& selector,
             field_type & dest,
             typename stk::mesh::FieldTraits< field_type >::data_type destscale,
             field_type & src,
             typename stk::mesh::FieldTraits< field_type >::data_type srcscale )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type data_type;

  if ( dest.mesh_meta_data_ordinal() == src.mesh_meta_data_ordinal() )
  {
    throw std::runtime_error( "cannot update self" );
  }

  const stk::mesh::DataTraits & dest_traits = dest.data_traits();
  const stk::mesh::DataTraits &  src_traits =  src.data_traits();

  if ( not dest_traits.is_floating_point or not src_traits.is_floating_point )
    throw std::runtime_error( "floating point fields expected" );

  const std::vector<stk::mesh::Bucket*>& all_buckets = bulk.buckets( stk::mesh::Node );
  //std::vector<stk::mesh::Bucket*> buckets;
  //stk::mesh::get_buckets( selector, all_buckets, buckets );

  // Is there a way to ask each bucket family just once?

  for ( std::size_t i=0; i<all_buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *all_buckets[i];
    int dest_field_data_size = stk::mesh::field_data_size( dest, bucket );
    int  src_field_data_size = stk::mesh::field_data_size(  src, bucket );

    if ( dest_field_data_size == 0 and src_field_data_size == 0 ) continue;
    if ( dest_field_data_size != src_field_data_size )
    {
      throw std::runtime_error( "fields do not match" );
    }

    std::size_t size = bucket.size() * dest_field_data_size/sizeof( data_type );

    data_type * dest_data = reinterpret_cast<data_type*>( bucket.field_data_location( dest ) );
    data_type *  src_data = reinterpret_cast<data_type*>( bucket.field_data_location(  src ) );

    // inner loops
    // maybe use kokkos or blas here

    if ( destscale==0 )
    {
      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] = srcscale * src_data[j];
      }
    }
    else
    {
      if ( destscale!=1 )
      {
        for ( std::size_t j=0; j<size; ++j )
        {
          dest_data[j] *= destscale;
        }
      }
      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] += srcscale * src_data[j];
      }
    }
  }
}


template <class field_type>
void Update( const stk::mesh::BulkData & bulk,
             field_type & dest,
             typename stk::mesh::FieldTraits< field_type >::data_type destscale,
             field_type & src1,
             typename stk::mesh::FieldTraits< field_type >::data_type srcscale1,
             field_type & src2,
             typename stk::mesh::FieldTraits< field_type >::data_type srcscale2 )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type data_type;

  if ( dest.mesh_meta_data_ordinal() == src1.mesh_meta_data_ordinal() or
       dest.mesh_meta_data_ordinal() == src2.mesh_meta_data_ordinal() )
  {
    throw std::runtime_error( "cannot update self" );
  }

  const stk::mesh::DataTraits & dest_traits = dest.data_traits();
  const stk::mesh::DataTraits & src1_traits = src1.data_traits();
  const stk::mesh::DataTraits & src2_traits = src2.data_traits();

  if ( not dest_traits.is_floating_point or
       not src1_traits.is_floating_point or
       not src2_traits.is_floating_point )
    throw std::runtime_error( "floating point fields expected" );

  const std::vector<stk::mesh::Bucket*>& all_buckets = bulk.buckets( stk::mesh::Node );

  // Is there a way to ask each bucket family just once?

  for ( std::size_t i=0; i<all_buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *all_buckets[i];
    int dest_field_data_size = stk::mesh::field_data_size( dest, bucket );
    int src1_field_data_size = stk::mesh::field_data_size( src1, bucket );
    int src2_field_data_size = stk::mesh::field_data_size( src2, bucket );

    if ( dest_field_data_size == 0 and
         src1_field_data_size == 0 and
         src2_field_data_size == 0 )
      continue;
    if ( dest_field_data_size != src1_field_data_size or
         dest_field_data_size != src2_field_data_size )
    {
      throw std::runtime_error( "fields do not match" );
    }

    std::size_t size = bucket.size() * dest_field_data_size/sizeof( data_type );

    data_type * dest_data = reinterpret_cast<data_type*>( bucket.field_data_location( dest ) );
    data_type * src1_data = reinterpret_cast<data_type*>( bucket.field_data_location( src1 ) );
    data_type * src2_data = reinterpret_cast<data_type*>( bucket.field_data_location( src2 ) );

    // inner loops
    // maybe use kokkos or blas here

    if ( destscale==0 )
    {
      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] = srcscale1 * src1_data[j] + srcscale2 * src2_data[j];
      }
    }
    else
    {
      if ( destscale!=1 )
      {
        for ( std::size_t j=0; j<size; ++j )
        {
          dest_data[j] *= destscale;
        }
      }
      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] += srcscale1 * src1_data[j] + srcscale2 * src2_data[j];
      }
    }
  }
}


template <class field_type>
void Update( const stk::mesh::BulkData & bulk,
             field_type & dest,
             typename stk::mesh::FieldTraits< field_type >::data_type destscale,
             field_type & src1,
             typename stk::mesh::FieldTraits< field_type >::data_type srcscale1,
             field_type & src2,
             typename stk::mesh::FieldTraits< field_type >::data_type srcscale2,
             field_type & src3,
             typename stk::mesh::FieldTraits< field_type >::data_type srcscale3 )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type data_type;

  if ( dest.mesh_meta_data_ordinal() == src1.mesh_meta_data_ordinal() or
       dest.mesh_meta_data_ordinal() == src2.mesh_meta_data_ordinal() or
       dest.mesh_meta_data_ordinal() == src3.mesh_meta_data_ordinal() )
  {
    throw std::runtime_error( "cannot update self" );
  }

  const stk::mesh::DataTraits & dest_traits = dest.data_traits();
  const stk::mesh::DataTraits & src1_traits = src1.data_traits();
  const stk::mesh::DataTraits & src2_traits = src2.data_traits();
  const stk::mesh::DataTraits & src3_traits = src3.data_traits();

  if ( not dest_traits.is_floating_point or
       not src1_traits.is_floating_point or
       not src2_traits.is_floating_point or
       not src3_traits.is_floating_point )
    throw std::runtime_error( "floating point fields expected" );

  const std::vector<stk::mesh::Bucket*>& all_buckets = bulk.buckets( stk::mesh::Node );

  // Is there a way to ask each bucket family just once?

  for ( std::size_t i=0; i<all_buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *all_buckets[i];
    int dest_field_data_size = stk::mesh::field_data_size( dest, bucket );
    int src1_field_data_size = stk::mesh::field_data_size( src1, bucket );
    int src2_field_data_size = stk::mesh::field_data_size( src2, bucket );
    int src3_field_data_size = stk::mesh::field_data_size( src3, bucket );

    if ( dest_field_data_size == 0 and
         src1_field_data_size == 0 and
         src2_field_data_size == 0 and
         src3_field_data_size == 0 )
      continue;
    if ( dest_field_data_size != src1_field_data_size or
         dest_field_data_size != src2_field_data_size or
         dest_field_data_size != src3_field_data_size )
    {
      throw std::runtime_error( "fields do not match" );
    }

    std::size_t size = bucket.size() * dest_field_data_size/sizeof( data_type );

    data_type * dest_data = reinterpret_cast<data_type*>( bucket.field_data_location( dest ) );
    data_type * src1_data = reinterpret_cast<data_type*>( bucket.field_data_location( src1 ) );
    data_type * src2_data = reinterpret_cast<data_type*>( bucket.field_data_location( src2 ) );
    data_type * src3_data = reinterpret_cast<data_type*>( bucket.field_data_location( src3 ) );

    // inner loops
    // maybe use kokkos or blas here

    if ( destscale==0 )
    {
      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] = srcscale1 * src1_data[j] + srcscale2 * src2_data[j] + srcscale3 * src3_data[j];
      }
    }
    else
    {
      if ( destscale!=1 )
      {
        for ( std::size_t j=0; j<size; ++j )
        {
          dest_data[j] *= destscale;
        }
      }
      for ( std::size_t j=0; j<size; ++j )
      {
        dest_data[j] += srcscale1 * src1_data[j] + srcscale2 * src2_data[j] + srcscale3 * src3_data[j];
      }
    }
  }
}


}

#endif
