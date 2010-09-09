#ifdef STKADAPTIVE

#include <boost/bind.hpp>

#include "stk_unrefine.H"
#include "stk_topology.H"
#include "stk_utils.H"
#include "stk_comm.H"


/*

  Gibt es eine Analogie zum Refinement? Wie kann ich die Informationen zum
  Unrefinement zusammentragen und dann das Unrefinement ausführen?

  Die lokalen Elemente werden zum Unrefinement markiert.

  Die zugehörigen Elternelemente liegen auf anderen Prozessoren. Aber das
  Ghosting ist vorhanden.

  Die Kinder müssen zu den Eltern. Wenn alle Kinder eines Elternelements zum
  Unrefinement vorgesehen sind, ist dieses Element drin.

  Jeder hanging Knoten der alle aktiven Elemente verliert wird entfernt. Das
  zugehörige Constraint auch, wenn es bereits existiert.

  Jeder hanging Knoten der frei wird bekommt ein Constraint. Sollte ein
  Oberflächenelement anliegen, wird das entfernt.

  Wird eine der Bedingungen nicht erfüllt, kann das entsprechende Element
  nicht unrefined werden. Damit ändern sich auch die Bedingungen für die
  anliegenden Constraints und die Oberflächenelemente.

  ->

  Das Set der zu entfernenden Entities anlegen. Alle Entities herausnehmen,
  die nicht entfernt werden können. Das parallel abgleichen. Danach löschen
  und die parent-Elemente aktiv machen.

 */


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::UnrefineSet::UnrefineSet(STK::Mesh* mesh)
  : mesh_(mesh)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::UnrefineSet::~UnrefineSet()
{
  remove_.clear();
  keep_.clear();

  parents_.clear();
  hanging_.clear();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::Unrefine(const std::vector<stk::mesh::EntityKey>& elements)
{
  std::map<stk::mesh::Entity*,std::set<unsigned> > unrefine_elements;

  //
  // Alle eigenen aktiven Elemente suchen, die entfernt werden sollen.

  FindUnrefineChildren(elements,unrefine_elements);

  //
  // Jedes unrefine-Kind zu allen sharing-Knoten kommunizieren

  std::map<stk::mesh::EntityKey,EntityKeySet> transfer_children;

  SendUnrefineChildren(unrefine_elements,transfer_children);

  //
  // Die eigenen Eltern-Elemente suchen, deren Kinder vollständig refined
  // werden

  SelectParents( transfer_children );

  // find all parent elements those child list is complete

  FindCompleteParents();

  // find constraints

  FindConstraints();

  // remove all illegal unrefinements

  RemoveImpossibleUnrefinement();

  // Now the unrefinement information is there. Action.

  CreateConstraints();

  DoUnrefine();

  //mesh_->Dump( "TestHexMeshOneByOne-failed" );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::FindUnrefineChildren(const std::vector<stk::mesh::EntityKey> & elements,
                                            std::map<stk::mesh::Entity*,std::set<unsigned> > & unrefine_elements)
{
  stk::mesh::Part & owned  = mesh_->OwnedPart();
  stk::mesh::Part & active = mesh_->ActivePart();

  for (std::vector<stk::mesh::EntityKey>::const_iterator i=elements.begin();
       i!=elements.end();
       ++i)
  {
    stk::mesh::Entity * e = mesh_->BulkData().get_entity( *i );

    if ( e!=NULL and
         has_superset( e->bucket(), owned ) and
         has_superset( e->bucket(), active ) )
    {
      unrefine_elements[e];

      // find any attached surface element and add it to the list of remove
      // elements

      stk::mesh::PairIterRelation rel = e->relations( stk::mesh::Node );
      EntitySet nodes;
      std::transform( rel.begin(), rel.end(),
                      std::inserter( nodes, nodes.begin() ),
                      boost::bind( &stk::mesh::Relation::entity, _1 ) );

      EntitySet subcells;

      const CellTopologyData * celltopology = stk::mesh::get_cell_topology( *e );

      // find any possible subcell

      for ( EntitySet::iterator i=nodes.begin(); i!=nodes.end(); ++i )
      {
        stk::mesh::Entity * n = *i;
        for ( unsigned subcell_rank = celltopology->dimension - 1;
              subcell_rank > stk::mesh::Node;
              --subcell_rank )
        {
          stk::mesh::PairIterRelation rel = n->relations( subcell_rank );
          std::transform( rel.begin(), rel.end(),
                          std::inserter( subcells, subcells.begin() ),
                          boost::bind( &stk::mesh::Relation::entity, _1 ) );
        }
      }

      // ignore all subcells that have nodes other than my ones

      for ( EntitySet::iterator i=subcells.begin(); i!=subcells.end(); ++i )
      {
        stk::mesh::Entity * sc = *i;

        // subcell needs to be active but not owned. There is no guarantee we
        // are going to remove it on its own process.

        if ( has_superset( sc->bucket(), active ) )
        {
          bool included = true;

          for ( stk::mesh::PairIterRelation rel = sc->relations( stk::mesh::Node );
                not rel.empty();
                ++rel )
          {
            stk::mesh::Entity * n = rel->entity();
            if ( nodes.find( n )==nodes.end() )
            {
              included = false;
              break;
            }
          }

          if ( included )
          {
            unrefine_elements[ sc ];
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::SendUnrefineChildren(std::map<stk::mesh::Entity*,std::set<unsigned> > & unrefine_elements,
                                            std::map<stk::mesh::EntityKey,EntityKeySet> & transfer_children)
{
  // loop all elements that are to be removed
  for (std::map<stk::mesh::Entity*,std::set<unsigned> >::iterator i=unrefine_elements.begin();
       i!=unrefine_elements.end();
       ++i)
  {
    stk::mesh::Entity* e = i->first;

    // find the processors all nodes of this element have in common
    FindCommonNodeSharing( e, i->second );

    if ( transfer_children.find( e->key() )!=transfer_children.end() )
    {
      throw std::runtime_error("doubled child element key");
    }

    // collect the EntityKeys of the nodes

    EntityKeySet & nodes = transfer_children[e->key()];
    stk::mesh::PairIterRelation ir = e->relations( stk::mesh::Node );

    std::transform(ir.begin(),ir.end(),
                   std::inserter(nodes,nodes.begin()),
                   boost::bind( &stk::mesh::Entity::key ,
                                boost::bind( &stk::mesh::Relation::entity , _1 ) ) );
  }

  // send the nodal EntityKeys of all unrefine elements to all processors the
  // element lives on

  stk::CommAll all(mesh_->parallel());

  for (std::map<stk::mesh::Entity*,std::set<unsigned> >::iterator i=unrefine_elements.begin();
       i!=unrefine_elements.end();
       ++i)
  {
    stk::mesh::Entity* e = i->first;
    std::set<unsigned> & procs = i->second;
    EntityKeySet & nodes = transfer_children[e->key()];
    for (std::set<unsigned>::iterator ip=procs.begin();
         ip!=procs.end();
         ++ip)
    {
      stk::CommBuffer & send_buffer = all.send_buffer( *ip );
      send_buffer
        .skip<stk::mesh::EntityKey>( 1 )
        .skip<unsigned>( 1 )
        .skip<stk::mesh::EntityKey>( nodes.size() );
    }
  }

  all.allocate_buffers(mesh_->parallel_size()/4);

  for (std::map<stk::mesh::Entity*,std::set<unsigned> >::iterator i=unrefine_elements.begin();
       i!=unrefine_elements.end();
       ++i)
  {
    stk::mesh::Entity* e = i->first;
    std::set<unsigned> & procs = i->second;
    EntityKeySet & nodes = transfer_children[e->key()];
    for (std::set<unsigned>::iterator ip=procs.begin();
         ip!=procs.end();
         ++ip)
    {
      stk::CommBuffer & send_buffer = all.send_buffer( *ip );
      send_buffer.pack<stk::mesh::EntityKey>( e->key() );
      send_buffer.pack<unsigned>( nodes.size() );
      std::for_each(nodes.begin(),nodes.end(),
                    boost::bind(&stk::CommBuffer::pack<stk::mesh::EntityKey>,&send_buffer,_1));
    }
  }

  all.communicate();

  for (unsigned p=0; p<mesh_->parallel_size(); ++p)
  {
    stk::CommBuffer & recv_buffer = all.recv_buffer(p);
    while (recv_buffer.remaining())
    {
      stk::mesh::EntityKey key;
      recv_buffer.unpack<stk::mesh::EntityKey>(key);

      unsigned size;
      recv_buffer.unpack<unsigned>(size);

      EntityKeySet & nodes = transfer_children[key];
      for ( unsigned i=0; i<size; ++i )
      {
        recv_buffer.unpack<stk::mesh::EntityKey>(key);
        nodes.insert(key);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::SelectParents( std::map<stk::mesh::EntityKey,EntityKeySet> & transfer_children )
{
  //stk::mesh::Part & owned  = mesh_->OwnedPart();
  stk::mesh::Part & active = mesh_->ActivePart();
  stk::mesh::BulkData & bulk = mesh_->BulkData();
  for ( std::map<stk::mesh::EntityKey,EntityKeySet>::iterator i=transfer_children.begin();
        i!=transfer_children.end();
        ++i )
  {
    EntitySet nodes;

    stk::mesh::Entity * child = bulk.get_entity( i->first );
    if ( child==NULL )
      continue;

    std::transform( i->second.begin(),
                    i->second.end(),
                    std::inserter( nodes , nodes.begin() ),
                    boost::bind( &stk::mesh::BulkData::get_entity , &mesh_->BulkData() , _1 ) );
    if ( nodes.find( NULL ) != nodes.end() )
    {
      throw std::runtime_error("receive node unavailable");
    }

    // look for the parent element
    stk::mesh::Entity* parent = InactiveParent( active, nodes, child->entity_rank() );
    if ( parent!=NULL )
    {
      ParentInfo & pi = parents_[parent];

      // remember child key
      pi.children.insert( child );
    }
  }

  // remember child nodes and hanging nodes
  for ( std::map<stk::mesh::Entity*, ParentInfo>::iterator i=parents_.begin();
        i!=parents_.end();
        ++i )
  {
    stk::mesh::Entity * parent = i->first;
    ParentInfo & pi = i->second;

    stk::mesh::PairIterRelation rel = parent->relations( stk::mesh::Node );
    const CellTopologyData * topology = stk::mesh::get_cell_topology( *parent );

    for ( unsigned i=0; i<topology->node_count; ++i )
    {
      pi.nodes.insert( rel[i].entity() );
    }

    for ( unsigned i=topology->node_count; i<rel.size(); ++i )
    {
      pi.hanging.insert( rel[i].entity() );
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::FindCompleteParents()
{
  for ( std::map<stk::mesh::Entity*,ParentInfo>::iterator i=parents_.begin();
        i!=parents_.end(); )
  {
    stk::mesh::Entity* parent = i->first;
    ParentInfo & pi = i->second;
    if ( not IsChildListComplete( parent, pi.children ) )
    {
      parents_.erase( i++ );
    }
    else
    {
      // remember element modifications
      std::copy( pi.children.begin(), pi.children.end(),
                 std::inserter( remove_, remove_.begin() ) );
      ++i;
    }
  }

  // now communicate remove_ to all participating processes

  {
    // Communicate any Entity to its ghosting, regardless of the owner. This
    // ensures all child elements will know about their removal.
    Distribution<GhostDistributionStrategy<EntityDistribution>, EntitySet>
      comm( mesh_->BulkData(), remove_ );
    comm();
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::FindConstraints()
{
  for ( std::map<stk::mesh::Entity*,ParentInfo>::iterator i=parents_.begin();
        i!=parents_.end();
        ++i )
  {
    stk::mesh::Entity* parent = i->first;
    ParentInfo & pi = i->second;
    EntitySet & hanging = pi.hanging;
    for ( EntitySet::iterator j=hanging.begin(); j!=hanging.end(); ++j )
    {
      stk::mesh::Entity * hn = *j;

      // see if a constraint is needed for this hanging node

      HangingNodeInfo & hni = hanging_[ hn ];

      // The node will either become (stay) a hanging node or it is not used
      // any more and will go away.

      hni.needconstraint = NeedConstraint( hn );

      // remember which parent uses this hanging node
      hni.parents.insert( parent );

      // see if there is a constraint

      for ( stk::mesh::PairIterRelation rel = hn->relations( stk::mesh::Constraint );
            not rel.empty();
            ++rel )
      {
        // TODO: Actually there could be a lot of constraints. A distinct part
        // for hanging node constraints is needed.

        stk::mesh::Entity * c = rel->entity();
        stk::mesh::PairIterRelation nodes = c->relations( stk::mesh::Node );
        if ( nodes.size()>0 and nodes[0].entity()==hn )
        {
          if ( hni.constraint!=NULL and hni.constraint!=c )
          {
            throw std::runtime_error( "multiple constraints on hanging node" );
          }
          hni.constraint = c;
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STK::UnrefineSet::NeedConstraint( stk::mesh::Entity * hn )
{
  // A node that is know to lose at least one refined element will become (or
  // stay) a hanging node if there is at least one active element that stays.

  stk::mesh::Part & active = mesh_->ActivePart();
  for ( stk::mesh::PairIterRelation rel = hn->relations( stk::mesh::Element );
        not rel.empty();
        ++rel )
  {
    stk::mesh::Entity * e = rel->entity();
    if ( has_superset( e->bucket(), active ) and
         remove_.find( e )==remove_.end() )
    {
      return true;
    }
  }

  // If the node is part of a hanging node face and not the hanging node
  // itself, it will stay.

  for ( stk::mesh::PairIterRelation rel = hn->relations( stk::mesh::Constraint );
        not rel.empty();
        ++rel )
  {
    stk::mesh::Entity * c = rel->entity();
    stk::mesh::PairIterRelation nodes = c->relations( stk::mesh::Node );
    if ( nodes[0].entity() != hn )
    {
      return true;
    }
  }

  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::RemoveImpossibleUnrefinement()
{
  stk::mesh::BulkData & bulk = mesh_->BulkData();

  stk::ParallelMachine pm = bulk.parallel();
  unsigned p_size = bulk.parallel_size();

  for ( ;; )
  {
    EntitySet reject;

    for ( std::map<stk::mesh::Entity*,ParentInfo>::iterator i=parents_.begin();
          i!=parents_.end();
      )
    {
      stk::mesh::Entity * parent = i->first;
      ParentInfo & pi = i->second;

      ++i;

      // Test if this parent can be unrefined and remove it from the
      // unrefinement list if not.
      UnrefinementPossible( parent, pi, reject );
    }

    stk::CommAll all( pm );

    for ( EntitySet::iterator i=reject.begin(); i!=reject.end(); ++i )
    {
      stk::mesh::Entity * parent = *i;
      std::vector<unsigned> procs;
      stk::mesh::comm_procs( *parent , procs );
      for ( std::vector<unsigned>::iterator j=procs.begin();
            j!=procs.end();
            ++j )
      {
        stk::CommBuffer & send_buffer = all.send_buffer( *j );
        send_buffer.pack<stk::mesh::EntityKey>( parent->key() );
      }
    }

    // if there is nothing to change, we are done

    bool global = all.allocate_buffers( p_size/4, false, reject.size()!=0 );
    if ( not global )
      break;

    for ( EntitySet::iterator i=reject.begin(); i!=reject.end(); ++i )
    {
      stk::mesh::Entity * parent = *i;
      std::vector<unsigned> procs;
      stk::mesh::comm_procs( *parent , procs );
      for ( std::vector<unsigned>::iterator j=procs.begin();
            j!=procs.end();
            ++j )
      {
        stk::CommBuffer & send_buffer = all.send_buffer( *j );
        send_buffer.pack<stk::mesh::EntityKey>( parent->key() );
      }
    }

    all.communicate();

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

        // remove parent from list if it is still there
        KeepParent( e );
      }
    }
  }

  // now communicate keep_ to all participating processes and update the
  // remove_ list

  {
    // Communicate any Entity to its ghosting, regardless of the owner. This
    // ensures all child elements will know about their removal.
    Distribution<GhostDistributionStrategy<EntityDistribution>, EntitySet>
      comm( mesh_->BulkData(), keep_ );
    comm();
  }

  for ( EntitySet::iterator i=keep_.begin(); i!=keep_.end(); ++i )
  {
    stk::mesh::Entity * e = *i;
    remove_.erase( e );
  }

  keep_.clear();

  //remove_.erase( keep_.begin(), keep_.end() );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STK::UnrefineSet::UnrefinementPossible( stk::mesh::Entity * parent,
                                             ParentInfo & pi,
                                             EntitySet & reject )
{
  // if any of our children has one surface with a hanging node constraint
  // that is not going to disappear, the unrefinement is not possible

  for ( EntitySet::iterator j=pi.children.begin(); j!=pi.children.end(); ++j )
  {
    stk::mesh::Entity * child = *j;

    stk::mesh::PairIterRelation nodes = child->relations( stk::mesh::Node );
    const CellTopologyData * celltopology = stk::mesh::get_cell_topology( *child );

    for ( unsigned subcell_rank = celltopology->dimension - 1;
          subcell_rank > stk::mesh::Node;
          --subcell_rank )
    {
      // TODO: Maybe we do not need to loop the surfaces. If we just look at
      // all nodes of the parent at once? Constraints on a parent side are
      // fine, on the other hand.

      // loop the subcells
      for ( unsigned nitr = 0; nitr < celltopology->subcell_count[subcell_rank]; ++nitr )
      {
        // For the potentially common subcell, get it's nodes and num_nodes
        const unsigned* sub_nodes = celltopology->subcell[subcell_rank][nitr].node;
        unsigned num_nodes = celltopology->subcell[subcell_rank][nitr].topology->node_count;

        EntitySet node_entities;
        for ( unsigned itr = 0; itr < num_nodes; ++itr )
        {
          stk::mesh::Entity * n = nodes[sub_nodes[itr]].entity();
          node_entities.insert( n );
        }

        // TODO: respect hanging part

        EntitySet constraints;
        FindCommonEntities( node_entities, constraints, stk::mesh::Constraint );

        for ( EntitySet::iterator k=constraints.begin(); k!=constraints.end(); ++k )
        {
          stk::mesh::Entity * c = *k;
          stk::mesh::PairIterRelation rel = c->relations( stk::mesh::Node );

          // If the hanging node is not included we have a hanging node on
          // this childs surface
          // If the hanging node is marked for removal, the constraint does
          // not matter.
          if ( node_entities.size()==rel.size()-1 and
               node_entities.find( rel[0].entity() )==node_entities.end() and
               hanging_.find( rel[0].entity() )==hanging_.end() )
          {
            // This child can not be removed. So its parent, all children of
            // its parent and any connected surface element must stay.

            reject.insert( parent );
            KeepParent( parent );

            // now test this parents' surface elements

            TestSurfaceUnrefinement( parent, pi, reject );

            return false;
          }
        }
      }
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::TestSurfaceUnrefinement( stk::mesh::Entity * parent,
                                                ParentInfo & pi,
                                                EntitySet & reject )
{
  // We have a parent element that was marked for unrefinement but has been
  // rejected since there are hanging nodes on its children. Now see if there
  // are any surfaces that are not unrefined now either.

  stk::mesh::Part & active = mesh_->ActivePart();

  stk::mesh::PairIterRelation rel = parent->relations( stk::mesh::Node );
  EntitySet nodes;
  std::transform( rel.begin(), rel.end(),
                  std::inserter( nodes, nodes.begin() ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );

  EntitySet subcells;

  const CellTopologyData * celltopology = stk::mesh::get_cell_topology( *parent );

  // find any possible subcell

  for ( EntitySet::iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    stk::mesh::Entity * n = *i;
    for ( unsigned subcell_rank = celltopology->dimension - 1;
          subcell_rank > stk::mesh::Node;
          --subcell_rank )
    {
      stk::mesh::PairIterRelation rel = n->relations( subcell_rank );
      std::transform( rel.begin(), rel.end(),
                      std::inserter( subcells, subcells.begin() ),
                      boost::bind( &stk::mesh::Relation::entity, _1 ) );
    }
  }

  // Of all the surfaces connected to this element see if any one was marked
  // for unrefinement solely due to the unrefinement of this parent element.

  for ( EntitySet::iterator i=subcells.begin(); i!=subcells.end(); ++i )
  {
    stk::mesh::Entity * sc = *i;

    // ignore all subcells that are active or have nodes other than my ones

    bool included = not has_superset( sc->bucket(), active );

    for ( stk::mesh::PairIterRelation rel = sc->relations( stk::mesh::Node );
          included and not rel.empty();
          ++rel )
    {
      stk::mesh::Entity * n = rel->entity();
      if ( nodes.find( n )==nodes.end() )
      {
        included = false;
        break;
      }
    }

    if ( included )
    {
      // Now we have a surface element of the given parent element that is not
      // active. Lets see if it is connected to any of the other unrefine
      // parents. If not, it is not going to be unrefined.

      EntitySet elements;
      FindCommonEntities( sc, elements, stk::mesh::Element );

      included = true;
      for ( EntitySet::iterator i=elements.begin(); i!=elements.end(); ++i )
      {
        stk::mesh::Entity * e = *i;
        if ( e!=parent and remove_.find( e )!=remove_.end() )
        {
          included = false;
        }
      }

      if ( included )
      {
        // Ok. So we need to stop this subcell from being unrefined.

        reject.insert( sc );
        KeepParent( sc );
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::KeepParent( stk::mesh::Entity * parent )
{
  // It is fine to call this more than once.
  std::map<stk::mesh::Entity*, ParentInfo>::iterator i = parents_.find( parent );
  if ( i==parents_.end() )
  {
    return;
  }

  ParentInfo & pi = i->second;

  // remove children from remove list
  for ( EntitySet::iterator j=pi.children.begin(); j!=pi.children.end(); ++j )
  {
    stk::mesh::Entity * child = *j;
    //remove_.erase( child );
    keep_.insert( child );
  }

  // reevaluate element constraints

  EntitySet & hanging = pi.hanging;
  for ( EntitySet::iterator j=hanging.begin(); j!=hanging.end(); ++j )
  {
    stk::mesh::Entity * hn = *j;

    std::map<stk::mesh::Entity*, HangingNodeInfo>::iterator k=hanging_.find( hn );
    if ( k==hanging_.end() )
    {
      throw std::logic_error( "hanging node not found" );
    }

    HangingNodeInfo & hni = k->second;
    //hni.needconstraint = NeedConstraint( hn );
    hni.needconstraint = true;
    hni.parents.erase( parent );

    // If there are no more parent elements to this hanging node, drop the
    // information. The node is going to stay normal. (?)
    if ( hni.parents.size()==0 )
    {
      hanging_.erase( k );
    }
  }

  // Surface elements of any kind have already been dealt with.

  // This parent is gone now
  parents_.erase( i );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::CreateConstraints()
{
  stk::mesh::BulkData & bulk = mesh_->BulkData();
  stk::ParallelMachine pm = bulk.parallel();
  unsigned p_size = bulk.parallel_size();
  unsigned p_rank = bulk.parallel_rank();

  // find which processor creates the hanging node constraint

  stk::CommAll all( pm );

  for ( std::map<stk::mesh::Entity*, HangingNodeInfo>::iterator i=hanging_.begin();
        i!=hanging_.end();
        ++i )
  {
    stk::mesh::Entity * hn = i->first;
    HangingNodeInfo & hni = i->second;

    if ( hni.needconstraint and hni.constraint==NULL )
    {
      std::vector<unsigned> procs;
      stk::mesh::comm_procs( *hn , procs );
      for ( std::vector<unsigned>::iterator j=procs.begin();
            j!=procs.end();
            ++j )
      {
        stk::CommBuffer & send_buffer = all.send_buffer( *j );
        send_buffer.pack<stk::mesh::EntityKey>( hn->key() );
      }
    }
  }

  all.allocate_buffers( p_size/4 );

  for ( std::map<stk::mesh::Entity*, HangingNodeInfo>::iterator i=hanging_.begin();
        i!=hanging_.end();
        ++i )
  {
    stk::mesh::Entity * hn = i->first;
    HangingNodeInfo & hni = i->second;

    if ( hni.needconstraint and hni.constraint==NULL )
    {
      std::vector<unsigned> procs;
      stk::mesh::comm_procs( *hn , procs );
      for ( std::vector<unsigned>::iterator j=procs.begin();
            j!=procs.end();
            ++j )
      {
        stk::CommBuffer & send_buffer = all.send_buffer( *j );
        send_buffer.pack<stk::mesh::EntityKey>( hn->key() );
      }
    }
  }

  all.communicate();

  std::map<stk::mesh::Entity*,unsigned> hangingnodes;

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

      // the last processor that claims the constraint gets it
      hangingnodes[e] = p;
    }
  }

  std::vector<std::size_t> requests( mesh_->MetaData().entity_rank_count(), 0 );

  for ( std::map<stk::mesh::Entity*, HangingNodeInfo>::iterator i=hanging_.begin();
        i!=hanging_.end();
        ++i )
  {
    stk::mesh::Entity * hn = i->first;
    HangingNodeInfo & hni = i->second;

    if ( hni.needconstraint and hni.constraint==NULL )
    {
      std::map<stk::mesh::Entity*,unsigned>::iterator j = hangingnodes.find( hn );
      if ( j==hangingnodes.end() or j->second <= p_rank )
      {
        requests[stk::mesh::Constraint] += 1;
      }
    }
  }

  EntityVector requested_entities;
  bulk.generate_new_entities( requests, requested_entities );

  EntityVector::iterator ci = requested_entities.begin();

  for ( std::map<stk::mesh::Entity*, HangingNodeInfo>::iterator i=hanging_.begin();
        i!=hanging_.end();
        ++i )
  {
    stk::mesh::Entity * hn = i->first;
    HangingNodeInfo & hni = i->second;

    if ( hni.needconstraint and hni.constraint==NULL )
    {
      std::map<stk::mesh::Entity*,unsigned>::iterator j = hangingnodes.find( hn );
      if ( j==hangingnodes.end() or j->second <= p_rank )
      {
        hni.constraint = *ci;
        ++ci;

        // Do not create the node relations now. The hanging node is obvious but
        // the other nodes are best dealt with later.
      }
    }
  }
}


class CommActivationStrategy : public STK::AuraCommNodeRelationsStrategy
{
public:

  explicit CommActivationStrategy( stk::mesh::BulkData & bulk )
    : AuraCommNodeRelationsStrategy( bulk )
  {
  }

  void set_active_part( stk::mesh::Part & active )
  {
    this->set_add_part( active );
  }

protected:

  template <class container>
  void unpack_details( container & entities,
                       stk::CommBuffer & recv_buffer,
                       stk::mesh::Entity* e )
  {
    // add active part to element
    this->change_entity_parts( e );

    AuraCommNodeRelationsStrategy::unpack_details<container>( entities, recv_buffer, e );
  }
};


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::UnrefineSet::DoUnrefine()
{
  stk::mesh::Part & active = mesh_->ActivePart();

  std::vector<stk::mesh::Part*> add_parts;
  std::vector<stk::mesh::Part*> remove_parts;

  add_parts.push_back( &active );

  stk::mesh::BulkData & bulk = mesh_->BulkData();
  //unsigned p_size = bulk.parallel_size();
  unsigned p_rank = bulk.parallel_rank();

  EntitySet activate;
  EntitySet unshared_nodes;
  //EntityProcSet removed_children;

  // remove child elements first

  for ( EntitySet::iterator j=remove_.begin(); j!=remove_.end(); ++j )
  {
    stk::mesh::Entity * child = *j;

    // need to cut any surface <-> main element connections

    for ( int rank=child->entity_rank()+1 ; rank<=stk::mesh::Element; ++rank )
    {
      stk::mesh::PairIterRelation rel = child->relations( rank );
      if ( not rel.empty() )
      {
        if ( rel.size()>1 )
        {
          throw std::logic_error( "can not have more than one main element" );
        }

        stk::mesh::Entity* main = rel[0].entity();
        bulk.destroy_relation( *main, *child );
      }
    }

#if 0
    // need to communicate element removal to aura

    if ( child->owner_rank() == p_rank )
    {
      std::vector<unsigned> procs;
      comm_procs( *child , procs );
      for ( std::vector<unsigned>::iterator i=procs.begin(); i!=procs.end(); ++i )
      {
        stk::mesh::EntityProc tmp( child, *i );
        removed_children.insert( tmp );
      }
    }
#endif

    // and destroy entity locally

    if ( not bulk.destroy_entity( child ) )
    {
      throw std::runtime_error( "failed to destroy child element" );
    }
  }

  for ( std::map<stk::mesh::Entity*,ParentInfo>::iterator i=parents_.begin();
        i!=parents_.end();
        ++i )
  {
    stk::mesh::Entity * parent = i->first;
    ParentInfo & pi = i->second;

    const CellTopologyData * celltopology = stk::mesh::get_cell_topology( *parent );

    // relations to nodes (including hanging ones)
    stk::mesh::PairIterRelation nodes = parent->relations( stk::mesh::Node );

    // decide on hanging node constraints
    unsigned pos = 0;
    for ( stk::mesh::PairIterRelation::iterator k=nodes.begin()+celltopology->node_count;
          k!=nodes.end();
          ++k, ++pos )
    {
      stk::mesh::Entity * hn = k->entity();
      HangingNodeInfo & hni = hanging_[hn];

      if ( hni.needconstraint and hni.constraint != NULL )
      {
        stk::mesh::PairIterRelation rel = hni.constraint->relations( stk::mesh::Node );
        if ( rel.empty() )
        {
          // Here we have a new constraint. Fill its node relations. This is
          // find the nodes that make up the parent element side.

          unsigned subcell_rank;
          unsigned subcell_pos;

          if ( pos < celltopology->side_count )
          {
            // face or edge for 3d or 2d elements
            subcell_rank = celltopology->dimension - 1;
            subcell_pos  = pos;
          }
          else if ( celltopology->dimension == stk::mesh::Element and
                    pos < celltopology->side_count+celltopology->edge_count )
          {
            // we have a hanging node on an edge here
            subcell_rank = stk::mesh::Edge;
            subcell_pos  = pos - celltopology->side_count;
          }
          else
          {
            // We have a hanging node on the element that has no
            // constraint. At least not based on this element.
            continue;
          }

          const unsigned* sub_nodes = celltopology->subcell[subcell_rank][subcell_pos].node;
          unsigned num_nodes        = celltopology->subcell[subcell_rank][subcell_pos].topology->node_count;

          // fill EntitySet to get the right node order

          EntitySet node_entities;
          for ( unsigned itr = 0; itr < num_nodes; ++itr )
          {
            stk::mesh::Entity * n = nodes[sub_nodes[itr]].entity();
            node_entities.insert( n );
          }

          int count = 0;
          bulk.declare_relation( *hni.constraint, *hn, count );

          for ( EntitySet::iterator i=node_entities.begin(); i!=node_entities.end(); ++i )
          {
            count += 1;
            bulk.declare_relation( *hni.constraint, **i, count );
          }

          // The constraint might span aura nodes. Those need to become shared
          // nodes.

          FindAuraNodes( bulk, hni.constraint, unshared_nodes );
        }
      }
      else
      {
        if ( hni.constraint!=NULL )
        {
          if ( not bulk.destroy_entity( hni.constraint ) )
          {
            throw std::runtime_error( "failed to destroy constraint" );
          }
        }
      }
    }

    // change parent on owner only (?!) Needs to be communicated!
    if ( parent->owner_rank()==p_rank )
    {
      bulk.change_entity_parts( *parent, add_parts, remove_parts );
      activate.insert( parent );
    }

    // remove hanging node relations
    for ( EntitySet::iterator k=pi.hanging.end(); k!=pi.hanging.begin(); )
    {
      --k;
      bulk.destroy_relation( *parent , **k );
    }
  }

  {
    // communicate element activation
    // send relation changes along with part changes
    // nodes without relations will be destroyed
    AuraDistribution<CommActivationStrategy, const EntitySet>
    comm( bulk, activate );
    comm.set_active_part( active );
    comm();
  }

#if 0
  {
    // communicate element removal to owner
    EntityProcDistribution<EntityRemovalStrategy, EntityProcSet>
      comm( bulk, removed_children );
    comm();
  }
#endif

  //AuraNodesToSharedNodes( bulk, unshared_nodes );

  // finally destroy unused nodes

  for ( std::map<stk::mesh::Entity*,ParentInfo>::iterator i=parents_.begin();
        i!=parents_.end();
        ++i )
  {
    //stk::mesh::Entity * parent = i->first;
    ParentInfo & pi = i->second;

    for ( EntitySet::iterator k=pi.hanging.begin(); k!=pi.hanging.end(); ++k )
    {
      stk::mesh::Entity * hn = *k;
      stk::mesh::PairIterRelation rel = hn->relations();

      // destroy each node just once
      if ( rel.empty() and hn->bucket().capacity() > 0 )
      {
        if ( not bulk.destroy_entity( hn ) )
        {
          throw std::runtime_error( "failed to destroy node" );
        }
      }
    }
  }
}

#endif
