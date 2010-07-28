#ifdef STKADAPTIVE

#ifndef COMM_HPP
#define COMM_HPP

#include "stk_types.H"
#include <boost/bind.hpp>

namespace STK
{

/// abstract communication template
template <class PackStrategy, class Container>
class Distribution : public PackStrategy
{
public:

  Distribution( stk::mesh::BulkData & bulk, Container & c )
    : PackStrategy( bulk ), m_c( c )
  {
  }

  bool operator()()
  {
    return this->communicate( m_c.begin(), m_c.end() );
  }

protected:

  Container & container() { return m_c; }

private:

  template <class iterator>
  bool communicate( const iterator & begin, const iterator & end )
  {
    stk::ParallelMachine pm = this->bulk().parallel();
    unsigned p_size = this->bulk().parallel_size();

    //int failed = 0;
    std::string msg;

    bool local = false;

    stk::CommAll all( pm );

#if 0
    try
    {
#endif
      for ( iterator i=begin; i!=end; ++i )
      {
        if ( this->matches( *i ) )
        {
          local = true;
          this->reserve( all, *i );
        }
      }
#if 0
    }
    catch ( std::exception &x )
    {
      failed = 1;
      msg = x.what();
    }

    stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );
    if ( failed )
    {
      throw std::runtime_error( msg );
    }
#endif

    bool global = all.allocate_buffers( p_size/4, false, local );
    if ( not global )
      return false;

#if 0
    try
    {
#endif
      for ( iterator i=begin; i!=end; ++i )
      {
        if ( this->matches( *i ) )
        {
          this->pack( all, *i );
        }
      }
#if 0
    }
    catch ( std::exception &x )
    {
      failed = 1;
      msg = x.what();
    }

    stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );
    if ( failed )
    {
      throw std::runtime_error( msg );
    }
#endif

    all.communicate();

#if 0
    try
    {
#endif
      for ( unsigned p=0; p<p_size; ++p )
      {
        stk::CommBuffer & recv_buffer = all.recv_buffer( p );
        while ( recv_buffer.remaining() )
        {
          this->unpack( this->container(), p, recv_buffer );
        }
      }
#if 0
    }
    catch ( std::exception &x )
    {
      failed = 1;
      msg = x.what();
    }

    stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );
    if ( failed )
    {
      throw std::runtime_error( msg );
    }
#endif

    this->done();

    return true;
  }

  Container & m_c;
};


//template <class C, class T>
class DefaultPackStrategy
{
public:

  explicit DefaultPackStrategy( stk::mesh::BulkData & bulk )
    : m_bulk( bulk )
  {
  }

protected:

  stk::mesh::BulkData & bulk() { return m_bulk; }

#if 0
  bool matches( T & t ) { return true; }

  void reserve( stk::CommAll & all, T & t ) { pack( all, t ); }

  void pack( stk::CommAll & all, T & t ) {}

  void unpack( C & container, unsigned p, stk::CommBuffer & recv_buffer ) {}
#endif

  void done() {}

private:
  stk::mesh::BulkData & m_bulk;
};


/// PackStrategy base class that is friend of stk::mesh::BulkData and stk::mesh::Entity
///
/// Communications that need to change internal STK data need to use a
/// PackStrategy derived from this class.
class DistributeEntityPackStrategy
{
public:

  explicit DistributeEntityPackStrategy( stk::mesh::BulkData & bulk )
    : m_bulk( bulk )
  {
  }

protected:

  stk::mesh::BulkData & bulk() { return m_bulk; }

  void set_add_part( stk::mesh::Part & part )
  {
    add_parts.push_back( &part );
  }

  void set_remove_part( stk::mesh::Part & part )
  {
    remove_parts.push_back( &part );
  }

  void change_entity_parts( stk::mesh::Entity * e )
  {
  }

  void make_ghost( stk::mesh::Entity * e, int ghosting, unsigned proc )
  {
  }

  void make_shared_plain( stk::mesh::Entity * e, unsigned proc )
  {
    make_ghost( e, 0, proc );
  }

  void make_shared( stk::mesh::Entity * e, unsigned proc )
  {
    make_shared_plain( e, proc );
    collect_new_sharing( e );
  }

  stk::mesh::Entity * add_new_ghost( stk::mesh::EntityKey key, int ghosting, unsigned owner )
  {
    stk::mesh::Entity & e = bulk().declare_entity( key.rank(), key.id(), add_parts );
    return &e;
  }

  stk::mesh::Entity * add_new_shared( stk::mesh::EntityKey key, unsigned owner )
  {
    return add_new_ghost( key, 0, owner );
  }

  void aura_to_shared( stk::mesh::Entity * e, unsigned proc )
  {
  }

  void done()
  {
  }

private:

  void communicate_sharing_partners();

  void collect_new_sharing( stk::mesh::Entity * e )
  {
    // If there are more than two processors, the entity might already be
    // shared. If I own the entity I need to introduce the other sharing
    // processors to each other.
    if ( bulk().parallel_size() > 2 and e->owner_rank()==bulk().parallel_rank() )
    {
      stk::mesh::PairIterEntityComm sharing = e->sharing();
      if ( sharing.size()>1 )
      {
        newly_shared.insert( e );
      }
    }
  }

  stk::mesh::BulkData & m_bulk;

  std::size_t entity_comm_size;

//   std::vector<stk::parallel::DistributedIndex::KeyType> add_new_keys;
//   std::vector<stk::parallel::DistributedIndex::KeyType> remove_existing_keys;

  std::vector<stk::mesh::Part*> add_parts;
  std::vector<stk::mesh::Part*> remove_parts;

  EntitySet newly_shared;
};


/// The owner of a shared entity introduces its sharing partners to each other.
///
/// This is used if a new process joins an existing sharing.
class NewSharingStrategy : public DistributeEntityPackStrategy
{
public:

  explicit NewSharingStrategy( stk::mesh::BulkData & bulk )
    : DistributeEntityPackStrategy( bulk )
  {
  }

protected:

  bool matches( stk::mesh::Entity * n )
  {
    return true;
  }

  void reserve( stk::CommAll & all, stk::mesh::Entity * n )
  {
    stk::mesh::PairIterEntityComm sharing = n->sharing();
    for ( stk::mesh::PairIterEntityComm::iterator i=sharing.begin();
          i!=sharing.end();
          ++i )
    {
      unsigned proc = i->proc;
      stk::CommBuffer & send_buffer = all.send_buffer( proc );
      for ( stk::mesh::PairIterEntityComm::iterator j=sharing.begin();
            j!=sharing.end();
            ++j )
      {
        if ( proc != j->proc )
        {
          send_buffer.skip<stk::mesh::EntityKey>( 1 )
            .skip<unsigned>( 1 );
        }
      }
    }
  }

  void pack( stk::CommAll & all, stk::mesh::Entity * n )
  {
    stk::mesh::PairIterEntityComm sharing = n->sharing();
    for ( stk::mesh::PairIterEntityComm::iterator i=sharing.begin();
          i!=sharing.end();
          ++i )
    {
      unsigned proc = i->proc;
      stk::CommBuffer & send_buffer = all.send_buffer( proc );
      for ( stk::mesh::PairIterEntityComm::iterator j=sharing.begin();
            j!=sharing.end();
            ++j )
      {
        if ( proc != j->proc )
        {
          send_buffer.pack<stk::mesh::EntityKey>( n->key() )
            .pack<unsigned>( j->proc );
        }
      }
    }
  }

  void unpack( const EntitySet & unshared,
               unsigned p,
               stk::CommBuffer & recv_buffer )
  {
    stk::mesh::EntityKey key;
    recv_buffer.unpack<stk::mesh::EntityKey>( key );

    stk::mesh::Entity * e = bulk().get_entity( key );
    if ( e==NULL )
    {
      throw std::runtime_error("receive entity unavailable");
    }

    unsigned proc;
    recv_buffer.unpack<unsigned>( proc );

    // create connection to other shared processor
    this->make_shared_plain( e, proc );
  }
};


/// Simple PackStrategy to distribute entities into the ghosting
template <class PackStrategy>
class GhostDistributionStrategy : public PackStrategy
{
public:

  explicit GhostDistributionStrategy( stk::mesh::BulkData & bulk )
    : PackStrategy( bulk )
  {
  }

protected:

  bool matches( stk::mesh::Entity * e )
  {
    return not e->comm().empty();
  }

  void reserve( stk::CommAll & all, stk::mesh::Entity * e )
  {
//     const std::vector<stk::mesh::Ghosting*> & ghosting = this->bulk().ghostings();

//     for ( std::vector<stk::mesh::Ghosting*>::const_iterator
//             ig = ghosting.begin() ; ig != ghosting.end() ; ++ig )
//     {
//       for ( stk::mesh::PairIterEntityComm comm = e->comm( **ig ) ;
//             not comm.empty() ;
//             ++comm )
//       {
    for ( stk::mesh::PairIterEntityComm comm=e->comm(); not comm.empty(); ++comm )
    {
      stk::CommBuffer & send_buffer = all.send_buffer( comm->proc );
      send_buffer.skip<stk::mesh::EntityKey>( 1 );

      this->reserve_details( send_buffer, e );
    }
  }

  void pack( stk::CommAll & all, stk::mesh::Entity * e )
  {
//     const std::vector<stk::mesh::Ghosting*> & ghosting = this->bulk().ghostings();

//     for ( std::vector<stk::mesh::Ghosting*>::const_iterator
//             ig = ghosting.begin() ; ig != ghosting.end() ; ++ig )
//     {
//       for ( stk::mesh::PairIterEntityComm comm = e->comm( **ig ) ;
//             not comm.empty() ;
//             ++comm )
//       {
    for ( stk::mesh::PairIterEntityComm comm=e->comm(); not comm.empty(); ++comm )
    {
      stk::CommBuffer & send_buffer = all.send_buffer( comm->proc );
      send_buffer.pack<stk::mesh::EntityKey>( e->key() );

      this->pack_details( send_buffer, e );
    }
  }

  template <class C>
  void unpack( C & entities, unsigned p, stk::CommBuffer & recv_buffer )
  {
    stk::mesh::EntityKey key;
    recv_buffer.unpack<stk::mesh::EntityKey>( key );

    stk::mesh::Entity * e = this->bulk().get_entity( key );
    if ( e==NULL )
    {
      throw std::runtime_error("receive entity unavailable");
    }

    this->unpack_details( entities, recv_buffer, e );
  }
};


/// Simple PackStrategy to distribute entities from their owners into the ghosting
template <class PackStrategy>
class AuraDistributionStrategy : public GhostDistributionStrategy<PackStrategy>
{
public:

  explicit AuraDistributionStrategy( stk::mesh::BulkData & bulk )
    : GhostDistributionStrategy<PackStrategy>( bulk )
  {
  }

protected:

  bool matches( stk::mesh::Entity * e )
  {
    return ( GhostDistributionStrategy<PackStrategy>::matches( e ) and
             e->owner_rank()==this->bulk().parallel_rank() );
  }
};


/// Abbreviation of aura communication template
template <class PackStrategy, class Container>
class AuraDistribution : public Distribution<AuraDistributionStrategy<PackStrategy>, Container>
{
public:

  AuraDistribution( stk::mesh::BulkData & bulk, Container & c )
    : Distribution<AuraDistributionStrategy<PackStrategy>, Container>( bulk, c )
  {
  }

};


/// simple collection of aura elements
class EntityDistribution : public DefaultPackStrategy
{
public:

  explicit EntityDistribution( stk::mesh::BulkData & bulk )
    : DefaultPackStrategy( bulk )
  {
  }

protected:

  void reserve_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
  }

  void pack_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
  }

  void unpack_details( EntitySet & entities,
                       stk::CommBuffer & recv_buffer,
                       stk::mesh::Entity* e )
  {
    entities.insert( e );
  }
};


/// communication of changed relations into the aura
///
/// Relation changes are not tracked by STK if the entities already exist.
///
/// Inherits DistributeEntityPackStrategy to allow private STK modifications
/// in derived classes.
class AuraCommNodeRelationsStrategy : public DistributeEntityPackStrategy
{
public:

  explicit AuraCommNodeRelationsStrategy( stk::mesh::BulkData & bulk )
    : DistributeEntityPackStrategy( bulk )
  {
  }

protected:

  void reserve_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
    // Send all relations. We filter on the receiving processor.
    stk::mesh::PairIterRelation ir = e->relations( stk::mesh::Node );
    send_buffer.skip<unsigned>( 1 ).skip<stk::mesh::EntityKey>( ir.size() );
  }

  void pack_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
    // Send all relations. We filter on the receiving processor.
    stk::mesh::PairIterRelation ir = e->relations( stk::mesh::Node );
    send_buffer.pack<unsigned>( ir.size() );
    for ( ; not ir.empty(); ++ir )
    {
      send_buffer.pack<stk::mesh::EntityKey>( ir->entity()->key() );
    }
  }

  template <class container>
  void unpack_details( container & entities,
                       stk::CommBuffer & recv_buffer,
                       stk::mesh::Entity* e )
  {
    stk::mesh::PairIterRelation ir = e->relations( stk::mesh::Node );
    EntitySet remove;
    std::transform( ir.begin(), ir.end(),
                    std::inserter( remove, remove.begin() ),
                    boost::bind( &stk::mesh::Relation::entity, _1 ) );
    unsigned size;
    recv_buffer.unpack<unsigned>( size );
    for ( unsigned i=0; i<size; ++i )
    {
      stk::mesh::EntityKey key;
      recv_buffer.unpack<stk::mesh::EntityKey>( key );
      stk::mesh::Entity * n = bulk().get_entity( key );

      // If we do not know the node we ignore it. It will be created and
      // connected by BulkData::modification_end().
      if ( n!=NULL )
      {
        // If the relation already exists, fine. Otherwise create it.
        bulk().declare_relation( *e , *n , i );
        remove.erase( n );
      }
    }
    for ( EntitySet::iterator i=remove.begin(); i!=remove.end(); ++i )
    {
      stk::mesh::Entity * n = *i;
      bulk().destroy_relation( *e, *n );
      stk::mesh::PairIterRelation rel = n->relations();
      if ( rel.empty() )
      {
        if ( not bulk().destroy_entity( n ) )
        {
          throw std::runtime_error( "failed to destroy node" );
        }
      }
    }
  }
};


class AuraToSharedStrategy : public DistributeEntityPackStrategy
{
public:

  explicit AuraToSharedStrategy( stk::mesh::BulkData & bulk )
    : DistributeEntityPackStrategy( bulk )
  {
  }

protected:

  bool matches( stk::mesh::Entity * n )
  {
    return true;
  }

  void reserve( stk::CommAll & all, stk::mesh::Entity * n )
  {
    stk::CommBuffer & send_buffer = all.send_buffer( n->owner_rank() );
    send_buffer.skip<stk::mesh::EntityKey>( 1 );
  }

  void pack( stk::CommAll & all, stk::mesh::Entity * n )
  {
    stk::CommBuffer & send_buffer = all.send_buffer( n->owner_rank() );
    send_buffer.pack<stk::mesh::EntityKey>( n->key() );

    this->aura_to_shared( n, n->owner_rank() );
  }

  void unpack( const EntitySet & unshared,
               unsigned p,
               stk::CommBuffer & recv_buffer )
  {
    stk::mesh::EntityKey key;
    recv_buffer.unpack<stk::mesh::EntityKey>( key );

    stk::mesh::Entity * n = bulk().get_entity( key );
    if ( n==NULL )
    {
      throw std::runtime_error("receive node unavailable");
    }

    this->aura_to_shared( n, p );
  }

};


void AuraNodesToSharedNodes( stk::mesh::BulkData & bulk , const EntitySet & unshared );

void FindAuraNodes( stk::mesh::BulkData & bulk, stk::mesh::Entity * e, EntitySet & aura );


/// Simple PackStrategy to distribute EntityProc collections
template <class PackStrategy>
class EntityProcDistributionStrategy : public PackStrategy
{
public:

  explicit EntityProcDistributionStrategy( stk::mesh::BulkData & bulk )
    : PackStrategy( bulk )
  {
  }

protected:

  bool matches( const stk::mesh::EntityProc & ep )
  {
    return true;
  }

  void reserve( stk::CommAll & all, const stk::mesh::EntityProc & ep )
  {
    stk::mesh::Entity * e = ep.first;
    unsigned proc = ep.second;
    stk::CommBuffer & send_buffer = all.send_buffer( proc );
    send_buffer.skip<stk::mesh::EntityKey>( 1 );

    this->reserve_details( send_buffer, e );
  }

  void pack( stk::CommAll & all, const stk::mesh::EntityProc & ep )
  {
    stk::mesh::Entity * e = ep.first;
    unsigned proc = ep.second;
    stk::CommBuffer & send_buffer = all.send_buffer( proc );
    send_buffer.pack<stk::mesh::EntityKey>( e->key() );

    this->pack_details( send_buffer, e );
  }

  template <class C>
  void unpack( C & entities, unsigned p, stk::CommBuffer & recv_buffer )
  {
    stk::mesh::EntityKey key;
    recv_buffer.unpack<stk::mesh::EntityKey>( key );

    stk::mesh::Entity * e = this->bulk().get_entity( key );
    if ( e==NULL )
    {
      throw std::runtime_error("receive entity unavailable");
    }

    this->unpack_details( entities, recv_buffer, e );
  }
};


/// Abbreviation of EntityProc communication template
template <class PackStrategy, class Container>
class EntityProcDistribution : public Distribution<EntityProcDistributionStrategy<PackStrategy>, Container>
{
public:

  EntityProcDistribution( stk::mesh::BulkData & bulk, Container & c )
    : Distribution<EntityProcDistributionStrategy<PackStrategy>, Container>( bulk, c )
  {
  }

};


/// communication of entity removal
///
/// Inherits DistributeEntityPackStrategy to allow private STK modifications
/// in derived classes.
class EntityRemovalStrategy : public DistributeEntityPackStrategy
{
public:

  explicit EntityRemovalStrategy( stk::mesh::BulkData & bulk )
    : DistributeEntityPackStrategy( bulk )
  {
  }

protected:

  void reserve_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
    // no additional load
  }

  void pack_details( stk::CommBuffer & send_buffer, stk::mesh::Entity* e )
  {
    // no additional load
  }

  template <class container>
  void unpack_details( container & entities,
                       stk::CommBuffer & recv_buffer,
                       stk::mesh::Entity* e )
  {
    if ( e->bucket().capacity() > 0 )
    {
      if ( not bulk().destroy_entity( e ) )
      {
        throw std::runtime_error( "failed to destroy entity" );
      }
    }
  }
};

}

#endif
#endif
