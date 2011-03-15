#ifdef STKADAPTIVE

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>

#include <boost/bind.hpp>

#include "stk_refine.H"
#include "stk_mesh.H"
#include "stk_topology.H"
#include "stk_utils.H"

#include "stk_comm.H"

namespace STK
{


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefineSet::RefineSet(Mesh* mesh)
  : mesh_(mesh)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::Refine(const std::vector<stk::mesh::EntityKey>& elements)
{
  // Each processor refines its own elements and communicates the refined
  // elements to its ghosting afterwards.

  stk::ParallelMachine pm = mesh_->parallel();
  unsigned p_size = mesh_->parallel_size();

  stk::mesh::BulkData & bulk = mesh_->BulkData();

  //stk::mesh::Part & owned  = mesh_->OwnedPart();
  //stk::mesh::Part & used   = mesh_->UsedPart();
  //stk::mesh::Part & active = mesh_->ActivePart();

  EntitySet refine_elements;
  EntitySet remote_refine;

  // Collect all local elements along with neighbor elements that need to be
  // refined, too.

  for (std::vector<stk::mesh::EntityKey>::const_iterator i=elements.begin();
       i!=elements.end();
       ++i)
  {
    stk::mesh::Entity * e = mesh_->BulkData().get_entity( *i );
    if ( e!=NULL )
    {
      RefineTransitive(e, refine_elements, remote_refine);
    }
  }

  // Each element refinement potentially spreads to the element's neighbors.
  // This needs to be communicated.

  for ( ;; )
  {
    stk::CommAll all( pm );

    for ( EntitySet::iterator i=remote_refine.begin();
          i!=remote_refine.end();
          ++i )
    {
      stk::mesh::Entity * e = *i;
      stk::CommBuffer & send_buffer = all.send_buffer( e->owner_rank() );
      send_buffer.pack<stk::mesh::EntityKey>( e->key() );
    }

    bool global = all.allocate_buffers( p_size/4, false, remote_refine.size()!=0 );
    if ( not global )
      break;

    for ( EntitySet::iterator i=remote_refine.begin();
          i!=remote_refine.end();
          ++i )
    {
      stk::mesh::Entity * e = *i;
      stk::CommBuffer & send_buffer = all.send_buffer( e->owner_rank() );
      send_buffer.pack<stk::mesh::EntityKey>( e->key() );
    }

    all.communicate();

    remote_refine.clear();

    for ( unsigned p=0; p<p_size; ++p )
    {
      stk::CommBuffer & recv_buffer = all.recv_buffer(p);
      while (recv_buffer.remaining())
      {
        stk::mesh::EntityKey key;
        recv_buffer.unpack<stk::mesh::EntityKey>(key);
        stk::mesh::Entity * e = bulk.get_entity( key );
        if ( e==NULL )
        {
          throw std::runtime_error("receive element unavailable");
        }

        RefineTransitive(e, refine_elements, remote_refine);
      }
    }
  }

  // find the nodes that must be created in order to refine

  for ( EntitySet::iterator i=refine_elements.begin();
        i!=refine_elements.end();
        ++i )
  {
    stk::mesh::Entity * e = *i;

    // Find nodes for all elements, including edge/face elements. The parallel
    // distribution might have some edge/face elements without their main
    // elements on one process.

    const TopologyType<void>* top = TopologyType<void>::GetTopologyType( e );
    top->FindChildNodes( *this, e );
  }

  // Each shared node knows the other processors it belongs to. But the new
  // hanging node might belong to just a subset of those. If only some of the
  // elements at these node are refined.

  SyncronizeNodeCreation();

  // communicate refinement to ghosts

  {
    AuraDistribution<EntityDistribution, EntitySet>
      comm( mesh_->BulkData(), refine_elements );
    comm();
  }

  // count new entities

  std::vector<std::size_t> requests( mesh_->MetaData().entity_rank_count(), 0 );
  CountNewEntities( refine_elements, requests );

  EntityVector requested_entities;
  mesh_->BulkData().generate_new_entities( requests, requested_entities );

  ConnectNewEntities( refine_elements, requests, requested_entities );

  // Communicate element deactivation to aura
  // At this point we communicate all relations between entities that already
  // exist in the aura. STK will add all relations to newly created entities
  // later on.

  CommunicateDeactivation( refine_elements );

  // interpolate field values to newly created nodes

  const stk::mesh::FieldVector & fields = mesh_->MetaData().get_fields();

  for ( std::map<EntityIdSet,CreationData>::iterator side=creation_.begin();
        side!=creation_.end();
        ++side )
  {
    const EntityIdSet & sideset = side->first;
    CreationData & data = side->second;

    if ( data.node_==NULL )
      throw std::runtime_error( "expect newly created node" );

    EntitySet sidenodes;

    for ( EntityIdSet::const_iterator j=sideset.begin();
          j!=sideset.end();
          ++j )
    {
      stk::mesh::Entity * n = bulk.get_entity( stk::mesh::Node, *j );
      if ( n==NULL )
        throw std::runtime_error( "node not found" );
      sidenodes.insert( n );
    }

    std::vector<double*> values;
    values.reserve( sideset.size() );

    for ( stk::mesh::FieldVector::const_iterator i=fields.begin();
          i!=fields.end();
          ++i )
    {
      stk::mesh::FieldBase & f = **i;
      int field_data_size = stk::mesh::field_data_size( f, data.node_->bucket() ) / sizeof( double );
      if ( field_data_size > 0 )
      {
        for ( EntitySet::const_iterator j=sidenodes.begin();
              j!=sidenodes.end();
              ++j )
        {
          stk::mesh::Entity * n = *j;

          int size = stk::mesh::field_data_size( f, n->bucket() ) / sizeof( double );
          if ( size == 0 or field_data_size != size )
          {
            // no field here!
            values.clear();
            break;
          }

          values.push_back( reinterpret_cast<double*>( stk::mesh::field_data( f , *n ) ) );
        }

        if ( values.size()>0 )
        {
          double fact = 1./values.size();
          double * v = reinterpret_cast<double*>( stk::mesh::field_data( f , *data.node_ ) );
          std::fill( v, v+field_data_size, 0 );
          for ( std::vector<double*>::iterator j=values.begin();
                j!=values.end();
                ++j )
          {
            const double * src = *j;
            for ( int k=0; k<field_data_size; ++k )
            {
              v[k] += fact*src[k];
            }
          }
        }

        values.clear();
      }
    }
  }

  //mesh_->Dump( "TestOneByOne-failed" );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::RefineTransitive( stk::mesh::Entity* e,
                                  EntitySet& refine_elements,
                                  EntitySet& remote_refine )
{
  stk::mesh::Part & owned  = mesh_->OwnedPart();
  stk::mesh::Part & active = mesh_->ActivePart();

  if ( has_superset( e->bucket(), active ) )
  {
    if ( has_superset( e->bucket(), owned ) )
    {
      std::pair<EntitySet::iterator, bool> result = refine_elements.insert( e );
      if ( result.second )
      {
        EntitySet nodes;
        stk::mesh::PairIterRelation iter = e->relations( stk::mesh::Node );
        std::transform( iter.begin(), iter.end(),
                        std::inserter( nodes, nodes.begin() ),
                        boost::bind( &stk::mesh::Relation::entity, _1 ) );
        stk::mesh::Entity * p = InactiveParent( active, nodes, e->entity_rank() );
        if ( p!=NULL )
        {
          // All active neighbors of the (passive) parent need to be refined, too.
          EntitySet neighbors;
          FindNeighbors( p, neighbors, stk::mesh::Element );
          for ( EntitySet::iterator i=neighbors.begin();
                i!=neighbors.end();
                ++i )
          {
            RefineTransitive( *i, refine_elements, remote_refine );
          }
        }

        // add edge/face elements that might have to be refined
        //
        // Rekursive! Thus I require face elements if edge elements are
        // needed.

        const CellTopologyData * celltopology = get_cell_topology( e->bucket() );
        int entity_rank = celltopology->dimension - 1;
        if ( entity_rank > stk::mesh::Node )
        {
          EntitySet faces;
          FindNeighbors( e, faces, entity_rank );
          for ( EntitySet::iterator i=faces.begin();
                i!=faces.end();
                ++i )
          {
            RefineTransitive( *i, refine_elements, remote_refine );
          }
        }
      }
    }
    else
    {
      remote_refine.insert( e );
    }
  }
}


class CommDeactivationStrategy : public AuraCommNodeRelationsStrategy
{
public:

  explicit CommDeactivationStrategy( stk::mesh::BulkData & bulk )
    : AuraCommNodeRelationsStrategy( bulk )
  {
  }

  void set_active_part( stk::mesh::Part & active )
  {
    this->set_remove_part( active );
  }

protected:

  template <class container>
  void unpack_details( container & entities,
                       stk::CommBuffer & recv_buffer,
                       stk::mesh::Entity* e )
  {
    // remove active part from element
    this->change_entity_parts( e );

    AuraCommNodeRelationsStrategy::unpack_details<container>( entities, recv_buffer, e );
  }
};


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::CommunicateDeactivation( const EntitySet & refine_elements )
{
  AuraDistribution<CommDeactivationStrategy, const EntitySet >
    comm( mesh_->BulkData(), refine_elements );
  comm.set_active_part( mesh_->ActivePart() );
  comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::CreateHangingNode(const stk::mesh::Entity* e,
                                  const std::vector<const stk::mesh::Entity*>& nodes)
{
  // See if the node set is already known. If it is, we are done.

  EntityIdSet nodeset;
  std::transform( nodes.begin(), nodes.end(),
                  std::inserter( nodeset, nodeset.begin() ),
                  boost::bind( &stk::mesh::Entity::identifier, _1 ) );

  CreationData & cd = creation_[nodeset];

  // if there are no neighbor elements the entry is new and needs to be setup
  if ( cd.neighborElements_.size()==0 )
  {
    // Get all elements connected to these nodes.
    FindNeighborElements( nodes, cd.neighborElements_ );
    if ( cd.neighborElements_.size()==0 )
    {
      throw std::logic_error( "there have to be neighbor elements" );
    }

    // find existing constraint attached to all nodes

    cd.constraint_ = FindConstraint( nodes );

    if ( cd.constraint_!=NULL )
    {
      stk::mesh::PairIterRelation ir = cd.constraint_->relations( stk::mesh::Node );
      if ( ir.empty() )
      {
        throw std::logic_error( "hanging node constraint without nodes" );
      }

      // The first node is the hang node. Always.
      cd.node_ = ir->entity();

      // If we have found a hanging node there is no need to create a new
      // one. It is already known on all participating processes thanks to the
      // aura. No communication!
      cd.procs_.clear();

      // add myself to be consistent.
      cd.procs_.insert( mesh_->parallel_rank() );

      // The constraint needs to be tested if it is still needed. That is if
      // all the constraint's elements are passive we have to remove it later
      // on.
    }
    else
    {
      // if the hanging node is connected to sharing nodes it will be a shared
      // node (or an aura node) and we need to communicate it

      FindCommonNodeSharing( nodes, cd.procs_ );

      // sharing includes other processors only, include myself.
      cd.procs_.insert( mesh_->parallel_rank() );

      // A new constraint is needed. At least it seems so. If all its elements
      // happen to become passive, we will not create it.
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::SyncronizeNodeCreation()
{
  // At this point we know some nodesets that need to be communicated among
  // different processors since the nodes or the adjoined elements belong to
  // different processors.
  //
  // We assume a valid ghosting, thus each processor that needs to create a
  // hanging node at the boundary of its domain knows all other processors
  // that share this boundary.
  //
  // If one of the other processors already created the hanging node, the
  // problem is solved. The ghosting has already communicated the hanging
  // node.
  //
  // If the hanging node is actually new, it must be assigned to one processor
  // only. However, it might be marked for creation on a subset of the
  // participating processors only. So it is required to find all processors
  // that actually need this node.

  unsigned p_rank = mesh_->parallel_rank();

  stk::CommAll all( mesh_->parallel() );

  for ( std::map<EntityIdSet,CreationData>::iterator j=creation_.begin();
        j!=creation_.end();
        ++j )
  {
    for (std::set<unsigned>::iterator i=j->second.procs_.begin(); i!=j->second.procs_.end(); ++i)
    {
      if ( *i!=p_rank )
      {
        stk::CommBuffer & send_buffer = all.send_buffer( *i );
        send_buffer.pack<unsigned>(j->first.size());
        std::for_each(j->first.begin(),j->first.end(),
                      boost::bind(&stk::CommBuffer::pack<stk::mesh::EntityId>,&send_buffer,_1));
      }
    }
  }

  all.allocate_buffers( mesh_->parallel_size()/4 );

  for ( std::map<EntityIdSet,CreationData>::iterator j=creation_.begin();
        j!=creation_.end();
        ++j )
  {
    for (std::set<unsigned>::iterator i=j->second.procs_.begin(); i!=j->second.procs_.end(); ++i)
    {
      if ( *i!=p_rank )
      {
        stk::CommBuffer & send_buffer = all.send_buffer( *i );
        send_buffer.pack<unsigned>(j->first.size());
        std::for_each(j->first.begin(),j->first.end(),
                      boost::bind(&stk::CommBuffer::pack<stk::mesh::EntityId>,&send_buffer,_1));
      }
    }
  }

  all.communicate();

  // Die Antworten zum Ausdünnen verwenden.

  std::map<EntityIdSet,std::set<unsigned> > newoffproc;

  unsigned p_size = mesh_->parallel_size();
  for ( unsigned p=0; p<p_size; ++p )
  {
    stk::CommBuffer & recv_buffer = all.recv_buffer( p );
    while ( recv_buffer.remaining() )
    {
      unsigned size;
      recv_buffer.unpack<unsigned>(size);
      EntityIdSet nodeset;
      for (unsigned i=0; i<size; ++i)
      {
        stk::mesh::EntityId id;
        recv_buffer.unpack<stk::mesh::EntityId>(id);
        nodeset.insert(id);
      }

      // Only mark this node set if I already know about it and the other side
      // tells about it as well.

      std::map<EntityIdSet,CreationData>::iterator ci = creation_.find( nodeset );
      if ( ci!=creation_.end() )
      {
        std::set<unsigned> & s = ci->second.procs_;
        if ( s.find(p) != s.end() )
        {
          newoffproc[nodeset].insert( p );
        }
      }
    }
  }

  // insert myself again
  for ( std::map<EntityIdSet,std::set<unsigned> >::iterator i=newoffproc.begin();
        i!=newoffproc.end();
        ++i )
  {
    i->second.insert( p_rank );
  }

  // keep the smaller maps

  for ( std::map<EntityIdSet,CreationData>::iterator j=creation_.begin();
        j!=creation_.end();
        ++j )
  {
    std::map<EntityIdSet,std::set<unsigned> >::iterator i=newoffproc.find( j->first );
    if ( i!=newoffproc.end() )
    {
      std::swap( i->second, j->second.procs_ );
    }
    else
    {
      // I did not receive anything about this one. So it is just me.
      j->second.procs_.clear();
      j->second.procs_.insert( p_rank );
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
stk::mesh::EntityId RefineSet::NodeId(const stk::mesh::Entity* n1,
                                      const stk::mesh::Entity* n2) const
{
  EntityIdSet nodeset;
  nodeset.insert(n1->identifier());
  nodeset.insert(n2->identifier());

  std::map<EntityIdSet,CreationData>::const_iterator i = creation_.find(nodeset);
  if ( i!=creation_.end() and i->second.node_!=NULL )
  {
    return i->second.node_->key().id();
  }
  else
  {
    return 0;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
stk::mesh::EntityId RefineSet::NodeId(const stk::mesh::Entity* n1,
                                      const stk::mesh::Entity* n2,
                                      const stk::mesh::Entity* n3,
                                      const stk::mesh::Entity* n4) const
{
  EntityIdSet nodeset;
  nodeset.insert(n1->identifier());
  nodeset.insert(n2->identifier());
  nodeset.insert(n3->identifier());
  nodeset.insert(n4->identifier());

  std::map<EntityIdSet,CreationData>::const_iterator i = creation_.find(nodeset);
  if ( i!=creation_.end() and i->second.node_!=NULL )
  {
    return i->second.node_->key().id();
  }
  else
  {
    return 0;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
stk::mesh::EntityId RefineSet::NodeId(const stk::mesh::Entity* n1,
                                      const stk::mesh::Entity* n2,
                                      const stk::mesh::Entity* n3,
                                      const stk::mesh::Entity* n4,
                                      const stk::mesh::Entity* n5,
                                      const stk::mesh::Entity* n6,
                                      const stk::mesh::Entity* n7,
                                      const stk::mesh::Entity* n8) const
{
  EntityIdSet nodeset;
  nodeset.insert(n1->identifier());
  nodeset.insert(n2->identifier());
  nodeset.insert(n3->identifier());
  nodeset.insert(n4->identifier());
  nodeset.insert(n5->identifier());
  nodeset.insert(n6->identifier());
  nodeset.insert(n7->identifier());
  nodeset.insert(n8->identifier());

  std::map<EntityIdSet,CreationData>::const_iterator i = creation_.find(nodeset);
  if ( i!=creation_.end() and i->second.node_!=NULL )
  {
    return i->second.node_->key().id();
  }
  else
  {
    return 0;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::CountNewEntities( const EntitySet & elements,
                                  std::vector<std::size_t> & requests )
{
  stk::mesh::BulkData & bulk = mesh_->BulkData();
  unsigned p_rank = mesh_->parallel_rank();

  // all information is available on all processors

  for (std::map<EntityIdSet,CreationData>::iterator i=creation_.begin();
       i!=creation_.end();
       ++i)
  {
    const EntityIdSet & facenodes = i->first;
    CreationData & cd = i->second;

    stk::mesh::Part & active = mesh_->ActivePart();

    cd.needconstraint_ = false;

    // Check if any of the elements connected to the parent set of nodes is
    // still active. If true we need a constraint to cover this set of nodes.
    //
    // Need to check all constraints here.
    for ( EntitySet::iterator j=cd.neighborElements_.begin();
          j!=cd.neighborElements_.end();
          ++j )
    {
      stk::mesh::Entity* e = *j;
      if ( stk::mesh::has_superset( e->bucket(), active ) and
           elements.find( e )==elements.end() )
      {
        cd.needconstraint_ = true;
        break;
      }
    }

    // If one of the side nodes happens to be a hanging node (and all side
    // nodes share a constraint) the new node will be a hanging node, since
    // the other side will be refined, too.

    for ( EntityIdSet::const_iterator j=facenodes.begin();
          not cd.needconstraint_ and j!=facenodes.end();
          ++j )
    {
      stk::mesh::Entity * n = bulk.get_entity( stk::mesh::Node, *j );
      stk::mesh::PairIterRelation rel = n->relations( stk::mesh::Constraint );
      for ( ; not rel.empty() ; ++rel )
      {
        stk::mesh::Entity & c = * rel->entity();
        stk::mesh::PairIterRelation nrel = c.relations( stk::mesh::Node );
        if ( nrel[0].entity()==n )
        {
          EntityIdSet constraint_nodes;
          for ( ; not nrel.empty(); ++nrel )
          {
            constraint_nodes.insert( nrel->entity()->key().id() );
          }

          bool found = true;
          for ( EntityIdSet::const_iterator k=facenodes.begin();
                k!=facenodes.end();
                ++k )
          {
            if ( constraint_nodes.count( *k )==0 )
            {
              found = false;
              break;
            }
          }

          if ( found )
          {
            cd.needconstraint_ = true;
            break;
          }
        }
      }
    }

    if ( cd.new_hanging_node() )
    {
      if ( cd.owner_rank()==p_rank )
      {
        requests[stk::mesh::Node] += 1;
        if ( cd.needconstraint_ )
        {
          requests[stk::mesh::Constraint] += 1;
        }
      }
    }
  }

  for ( EntitySet::const_iterator i=elements.begin();
        i!=elements.end();
        ++i )
  {
    stk::mesh::Entity* e = *i;
    if ( e->owner_rank()==p_rank )
    {
      if ( e->entity_rank()!=stk::mesh::Element )
      {
        //stk::mesh::print_entity_key( std::cout , mesh_->MetaData() , e->key() );
        if ( FaceNeedsRefinement( e ) )
        {
          const TopologyType<void> * top = TopologyType<void>::GetTopologyType( e );
          requests[ e->entity_rank() ] += top->NumChildElements();
        }
      }
      else
      {
        const TopologyType<void> * top = TopologyType<void>::GetTopologyType( e );
        requests[ stk::mesh::Element ] += top->NumChildElements();
      }
    }
  }
}


class NodesPackStrategy : public DistributeEntityPackStrategy
{
public:

  explicit NodesPackStrategy( stk::mesh::BulkData & bulk )
    : DistributeEntityPackStrategy( bulk )
  {
    p_rank = bulk.parallel_rank();

    // parts to set on newly created entities
    this->set_add_part( bulk.mesh_meta_data().universal_part() );
    this->set_add_part( bulk.mesh_meta_data().globally_shared_part() );
  }

  void set_active_part( stk::mesh::Part & active )
  {
    this->set_add_part( active );
  }

protected:

  bool matches( std::pair<const EntityIdSet,CreationData> & i )
  {
    CreationData & cd = i.second;
    return cd.owner_rank()==p_rank and cd.procs_.size()>1;
  }

  void reserve( stk::CommAll & all, std::pair<const EntityIdSet,CreationData> & i )
  {
    const EntityIdSet & face = i.first;
    CreationData & cd = i.second;

    std::set<unsigned>::iterator j = cd.procs_.begin();
    for ( ++j; j!=cd.procs_.end(); ++j )
    {
      stk::CommBuffer & send_buffer = all.send_buffer( *j );

      // node set
      send_buffer
        .skip<unsigned>( 1 )
        .skip<stk::mesh::EntityId>( face.size() )
        .skip<stk::mesh::EntityId>( 2 )
        ;
    }
  }

  void pack( stk::CommAll & all, std::pair<const EntityIdSet,CreationData> & i )
  {
    const EntityIdSet & face = i.first;
    CreationData & cd = i.second;

    stk::mesh::Entity * n = cd.node_;
    stk::mesh::Entity * c = NULL;

    if ( cd.needconstraint_ )
    {
      c = cd.constraint_;
    }

    std::set<unsigned>::iterator j = cd.procs_.begin();
    for ( ++j; j!=cd.procs_.end(); ++j )
    {
      unsigned proc = *j;
      stk::CommBuffer & send_buffer = all.send_buffer( proc );

      // node set
      send_buffer.pack<unsigned>( face.size() );
      std::for_each( face.begin(), face.end(),
                     boost::bind( &stk::CommBuffer::pack<stk::mesh::EntityId>,&send_buffer, _1 ) );

      // hanging node
      send_buffer.pack<stk::mesh::EntityId>( n->key().id() );

      // make this node a shared node, needs calls to add_new_shared() on the
      // receiving processor
      make_shared( n, proc );

      if ( cd.needconstraint_ )
      {
        send_buffer.pack<stk::mesh::EntityId>( c->key().id() );

        // make the constraint a ghost
        make_ghost( c, 1, proc );
      }
      else
      {
        send_buffer.pack<stk::mesh::EntityId>( 0 );
      }
    }
  }

  void unpack( std::map<EntityIdSet,CreationData> & creation,
               unsigned p,
               stk::CommBuffer & recv_buffer )
  {
    unsigned size;
    recv_buffer.unpack<unsigned>(size);
    EntityIdSet nodeset;
    for (unsigned i=0; i<size; ++i)
    {
      stk::mesh::EntityId id;
      recv_buffer.unpack<stk::mesh::EntityId>(id);
      nodeset.insert(id);
    }

    CreationData & cd = creation[nodeset];

    {
      stk::mesh::EntityId id;
      recv_buffer.unpack<stk::mesh::EntityId>(id);

      stk::mesh::EntityKey key( stk::mesh::Node, id );
      cd.node_ = add_new_shared( key, p );
    }

    {
      stk::mesh::EntityId id;
      recv_buffer.unpack<stk::mesh::EntityId>(id);

      if ( id>0 )
      {
        stk::mesh::EntityKey key( stk::mesh::Constraint, id );
        cd.constraint_ = add_new_ghost( key, 1, p );
        cd.needconstraint_ = true;
      }
      else
      {
        cd.constraint_ = NULL;
        cd.needconstraint_ = false;
      }
    }
  }

private:
  unsigned p_rank;
};


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::ConnectNewEntities( EntitySet & elements,
                                    std::vector<std::size_t> & requests,
                                    EntityVector & entities )
{
  unsigned p_rank = mesh_->parallel_rank();

  stk::mesh::BulkData & bulk = mesh_->BulkData();

  stk::mesh::Part & active = mesh_->ActivePart();

  std::vector<stk::mesh::Part*> empty;
  std::vector<stk::mesh::Part*> active_parts;
  active_parts.push_back( &active );

  // collect new entity iterators

  std::vector<EntityVector::iterator> new_entities;
  new_entities.reserve( requests.size() );
  EntityVector::iterator pos=entities.begin();
  for ( std::vector<std::size_t>::iterator i=requests.begin();
        i!=requests.end();
        ++i )
  {
    new_entities.push_back( pos );
    pos += *i;
  }

  // do nodes and constraints

  EntityVector::iterator & new_nodes       = new_entities[stk::mesh::Node      ];
  EntityVector::iterator & new_constraints = new_entities[stk::mesh::Constraint];

  for (std::map<EntityIdSet,CreationData>::iterator i=creation_.begin();
       i!=creation_.end();
       ++i)
  {
    //const EntityIdSet & face = i->first;
    CreationData & cd = i->second;

    if ( cd.new_hanging_node() )
    {
      if ( cd.owner_rank()==p_rank )
      {
        cd.node_ = *new_nodes;
        ++new_nodes;

        // add active part to new node
        bulk.change_entity_parts( *cd.node_, active_parts, empty );

        if ( cd.needconstraint_ )
        {
          cd.constraint_ = *new_constraints;
          ++new_constraints;
        }
      }
    }
  }

  {
    Distribution<NodesPackStrategy, std::map<EntityIdSet,CreationData> >
      comm( bulk, creation_ );
    comm.set_active_part( active );
    comm();
  }

  EntityProcSet removed_constraints;
  EntitySet created_constraints;
  EntitySet unshared_nodes;

  for (std::map<EntityIdSet,CreationData>::iterator i=creation_.begin();
       i!=creation_.end();
       ++i)
  {
    CreationData & cd = i->second;

    if ( cd.needconstraint_ )
    {
      // if a constraint is needed is has to be there already.
      stk::mesh::Entity * constraint = cd.constraint_;
      if ( constraint==NULL )
        throw std::logic_error( "constraint is needed, should be there at this point" );

      stk::mesh::PairIterRelation rel = constraint->relations( stk::mesh::Node );

      // If there are no relations it is a new constraint and needs to be filled.
      if ( rel.size()==0 and constraint->owner_rank()==p_rank )
      {
        int count = 0;
        bulk.declare_relation( *constraint, *cd.node_, count );

        for ( EntityIdSet::iterator j=i->first.begin();
              j!=i->first.end();
              ++j )
        {
          count += 1;
          stk::mesh::Entity * n = bulk.get_entity( stk::mesh::Node, *j );
          bulk.declare_relation( *constraint, *n, count );
        }

        created_constraints.insert( constraint );

        // The constraint might span aura nodes. Those need to become shared
        // nodes.

        FindAuraNodes( bulk, constraint, unshared_nodes );
      }
    }
    else
    {
      stk::mesh::Entity * constraint = cd.constraint_;
      //if ( constraint!=NULL and constraint->owner_rank()==p_rank )
      if ( constraint!=NULL )
      {
        if ( constraint->owner_rank()==p_rank )
        {
          if ( not bulk.destroy_entity( constraint ) )
          {
            throw std::runtime_error( "failed to destroy constraint" );
          }
        }
        else
        {
          stk::mesh::EntityProc tmp( constraint, constraint->owner_rank() );
          removed_constraints.insert( tmp );
        }
      }
    }
  }

  {
    // communicate new constraint connections to aura
    // this might be unnecessary
    AuraDistribution<AuraCommNodeRelationsStrategy, EntitySet>
      comm( bulk, created_constraints );
    comm();
  }

  {
    // communicate constraint removal to constraint owner
    EntityProcDistribution<EntityRemovalStrategy, EntityProcSet>
      comm( bulk, removed_constraints );
    comm();
  }

  // do elements
  // go backwards to make sure the main elements are refined before the
  // edge/face elements

  for ( EntitySet::iterator i=elements.end(); i!=elements.begin(); )
  {
    --i;
    stk::mesh::Entity* e = *i;
    if ( e->owner_rank()==p_rank )
    {
      if ( e->entity_rank()!=stk::mesh::Element )
      {
        if ( FaceNeedsRefinement( e ) )
        {
          EntityVector::iterator & elements_iter = new_entities[ e->entity_rank() ];
          EntityVector::iterator elements_begin = elements_iter;

          RefineElement( bulk, e, elements_iter, unshared_nodes );

          EntitySet children;
          std::copy( elements_begin, elements_iter, std::inserter( children, children.begin() ) );

          // find main element to new edge/face element and build connection

          for ( int rank=e->entity_rank()+1 ; rank<=stk::mesh::Element; ++rank )
          {
            stk::mesh::PairIterRelation rel = e->relations( rank );
            if ( not rel.empty() )
            {
              if ( rel.size()>1 )
              {
                throw std::logic_error( "can not have more than one main element" );
              }

              stk::mesh::Entity* main_parent = rel[0].entity();

              for ( EntitySet::iterator i=children.begin(); i!=children.end(); ++i )
              {
                stk::mesh::Entity * child = *i;
#if 0
                stk::mesh::Entity * main_parent = NULL;

                // Take the first parent element with a matching owner. We
                // demand matching owner at the moment.
                for ( stk::mesh::PairIterRelation::iterator ir=rel.begin();
                      ir!=rel.end();
                      ++ir )
                {
                  stk::mesh::Entity * main = ir->entity();
                  if ( main->owner_rank()==child->owner_rank() )
                  {
                    main_parent = main;
                    break;
                  }
                }

                if ( main_parent==NULL )
                {
                  throw std::runtime_error( "did not find a matching parent" );
                }
#endif

                // child teilt alle Knoten mit mindestens einem aktiven Element

                EntitySet es;
                FindCommonEntities( child, es, rank );

                for ( EntitySet::iterator j=es.begin(); j!=es.end(); ++j )
                {
                  stk::mesh::Entity * me = *j;
                  if ( stk::mesh::has_superset( me->bucket(), active ) and
                       me->owner_rank()==child->owner_rank() )
                  {
                    if ( is_child( me, main_parent ) )
                    {
                      // find side of child at me and connect
                      std::vector<stk::mesh::Entity*> node_entities;
                      stk::mesh::PairIterRelation child_nodes  = child->relations( stk::mesh::Node );
                      node_entities.reserve( child_nodes.size() );
                      std::transform( child_nodes.begin(), child_nodes.end(),
                                      std::back_inserter( node_entities ),
                                      boost::bind( &stk::mesh::Relation::entity, _1 ) );
                      //unsigned local_side_num = stk::mesh::element_local_side_id( *me, get_cell_topology( *child ), node_entities );
                      int local_side_num = stk::mesh::get_entity_subcell_id( *me,
                                                                             child->entity_rank(),
                                                                             get_cell_topology( *child ),
                                                                             node_entities );
                      stk::mesh::declare_element_side( *me, *child, local_side_num, NULL );
                    }
                  }
                }
              }

              break;
            }
          }
        }
        else
        {
          // So we have a surface element that is not refined. It must not be
          // communicated later on either. Remove it now.
          elements.erase( i++ );
        }
      }
      else
      {
        RefineElement( bulk, e, new_entities[ stk::mesh::Element ], unshared_nodes );
      }
    }
  }

  //AuraNodesToSharedNodes( bulk, unshared_nodes );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void RefineSet::RefineElement( stk::mesh::BulkData & bulk,
                               stk::mesh::Entity * e,
                               EntityVector::iterator & elements_iter,
                               EntitySet & unshared_nodes )
{
  const TopologyType<void> * top = TopologyType<void>::GetTopologyType( e );

  stk::mesh::PairIterRelation iter = e->relations( stk::mesh::Node );
  EntityVector nodes( iter.size() );
  std::transform(iter.begin(),iter.end(),
                 &nodes[0],
                 boost::bind(&stk::mesh::Relation::entity, _1));

  // get (or create) the new nodes that might be hanging
  EntityVector hanging;
  top->HangingNodes( this, mesh_, nodes, hanging );

  // create new active elements
  top->CreateElements( mesh_, e->entity_rank(), e->bucket().superset_part_ordinals(), nodes, hanging, elements_iter );

  // Add hanging nodes to parent element. So we have changing
  // relations. This is an implicit but feasible solution.

  unsigned count = nodes.size();
  for ( EntityVector::iterator j=hanging.begin();
        j!=hanging.end();
        ++j )
  {
    bulk.declare_relation( *e, **j, count );
    count += 1;
  }

  // At this point we could need some interpolation of element values.

  // deactivate new parent element

  std::vector<stk::mesh::Part*> empty;
  std::vector<stk::mesh::Part*> active_parts( 1, & mesh_->ActivePart() );

  bulk.change_entity_parts( *e, empty, active_parts );

  // Some hanging nodes might already be known as ghosts, however we are
  // about to use those nodes as shared nodes. This needs to be
  // communicated.

  FindAuraNodes( bulk, e, unshared_nodes );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool RefineSet::FaceNeedsRefinement( stk::mesh::Entity * e )
{
  stk::mesh::PairIterRelation nodes = e->relations( stk::mesh::Node );
  EntityIdSet nodeset;
  std::transform( nodes.begin(), nodes.end(),
                  std::inserter( nodeset, nodeset.begin() ),
                  boost::bind( &stk::mesh::Entity::identifier,
                               boost::bind( &stk::mesh::Relation::entity, _1 ) ) );
  std::map<EntityIdSet,CreationData>::iterator j = creation_.find( nodeset );
  if ( j==creation_.end() )
  {
    throw std::runtime_error( "face element with illegal node set" );
  }

  CreationData & cd = j->second;
  return not cd.needconstraint_;
}

}

#endif
