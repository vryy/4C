#ifdef STKADAPTIVE

#include "stk_iterator.H"
#include "../stk_refine/stk_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::ElementIterator::ElementIterator( stk::mesh::BulkData & bulk,
                                       const stk::mesh::Selector & selector )
  : bulk_( bulk ),
    selector_( selector ),
    elements_( bulk.buckets( stk::mesh::Element ) ),
    i_( elements_.begin() ),
    bucket_( NULL ),
    celltopology_( NULL )
{
  next_bucket();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::ElementIterator::next_bucket()
{
  while ( i_ != elements_.end() )
  {
    bucket_ = *i_;
    if ( selector_( *bucket_ ) )
    {
      celltopology_ = get_cell_topology( *bucket_ );

      j_ = bucket_->begin();
      return;
    }
    ++i_;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::SideIterator::SideIterator( ElementIterator & ei )
  : ei_( ei ),
    nodes_( ei_.element().relations( stk::mesh::Node ) ),
    side_dim_( ei_.topology().dimension - 1 ),
    side_num_( 0 )
{
  new_side();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STK::SideIterator::done()
{
  return side_num_ == ei_.topology().subcell_count[side_dim_];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::SideIterator::next()
{
  ++side_num_;
  if ( side_num_ < ei_.topology().subcell_count[side_dim_] )
  {
    new_side();
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::SideIterator::new_side()
{
  side_nodes_.clear();
  side_elements_.clear();
  other_side_num_.clear();
  other_topology_.clear();
  hanging_side_ = false;

  EntitySet side_elements;

  side_ = & ei_.topology().subcell[side_dim_][side_num_];
  side_topology_ = side_->topology;

  for ( unsigned j=0; j<side_topology_->node_count; ++j )
  {
    stk::mesh::Entity * n = nodes_[side_->node[j]].entity();
    side_nodes_.insert( n );
  }

  find_side_elements( side_nodes_, side_elements );

  if ( side_elements.count( &ei_.element() )==0 )
    throw std::runtime_error( "element not referenced by side nodes" );

  side_elements.erase( &ei_.element() );

  if ( side_elements.size()>1 )
    throw std::runtime_error( "element side with more that two elements" );

  if ( side_elements.size()==1 )
  {
    find_other_element_side( **side_elements.begin() );
  }
  else
  {
    // test for hanging node
    EntitySet constraints;
    FindCommonEntities( side_nodes_, constraints, stk::mesh::Constraint );
    if ( constraints.size()>0 )
    {
      if ( constraints.size()>1 )
      {
        for ( EntitySet::iterator i=constraints.begin();
              i!=constraints.end();
              ++i )
        {
          stk::mesh::Entity & c = **i;
          std::cout << "CONSTRAINT[" << c.key().id() << "]: ";
          stk::mesh::PairIterRelation rel = c.relations( stk::mesh::Node );
          for ( ; not rel.empty(); ++rel )
          {
            std::cout << rel->entity()->key().id() << " ";
          }
          std::cout << "\n";
        }
        throw std::runtime_error( "cannot have more than one constraint at side nodes" );
      }

      hanging_side_ = true;

      stk::mesh::Entity & constraint = **constraints.begin();
      stk::mesh::PairIterRelation constraint_nodes = constraint.relations( stk::mesh::Node );
      stk::mesh::Entity & n = *constraint_nodes[0].entity();
      if ( side_nodes_.count( &n )>0 )
      {
        // We own the hanging node here. Look for the other side
        // element. It has to be there.

        side_nodes_.erase( &n );
        for ( ++constraint_nodes; not constraint_nodes.empty(); ++constraint_nodes )
        {
          stk::mesh::Entity & n = *constraint_nodes->entity();
          side_nodes_.insert( &n );
        }

        side_elements.clear();

        find_side_elements( side_nodes_, side_elements );

        if ( side_elements.size()!=1 )
          throw std::runtime_error( "element side with hanging node failed to find other element" );

        find_other_element_side( **side_elements.begin() );

        if ( side_elements_.size()!=1 )
          throw std::runtime_error( "expect to find one element" );
      }
      else
      {
        // the hanging node belongs to the other side

        // Add the hanging node and look for all elements covered by
        // this set now.
        //
        // Look for all elements that share some nodes from side_nodes
        // with e along a side of the same dimension.
        side_nodes_.insert( &n );
        side_elements.clear();

        find_other_hanging_sides();

#if 0
        // If we are at a processor boundary, we might find just one side
        // element here. In this case we are inside an aura element anyway.
        if ( side_elements_.size() < 2 )
          throw std::runtime_error( "expect to find more than one element" );
#endif
      }
    }
    else
    {
      // boundary side, no neighbour
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::SideIterator::find_side_elements( const EntitySet & side_nodes,
                                            EntitySet & side_elements )
{
  EntitySet::iterator nodes = side_nodes.begin();
  stk::mesh::Entity * n = *nodes;
  for ( stk::mesh::PairIterRelation elements = n->relations( stk::mesh::Element );
        not elements.empty();
        ++elements )
  {
    stk::mesh::Entity & ele = *elements->entity();
    if ( ei_.selector()( ele.bucket() ) )
    {
      side_elements.insert( &ele );
    }
  }

  for ( ++nodes; nodes!=side_nodes.end(); ++nodes )
  {
    EntitySet this_side_elements;
    stk::mesh::Entity * n = *nodes;

    for ( stk::mesh::PairIterRelation elements = n->relations( stk::mesh::Element );
          not elements.empty();
          ++elements )
    {
      stk::mesh::Entity & ele = *elements->entity();
      if ( ei_.selector()( ele.bucket() ) )
      {
        if ( side_elements.count( &ele ) > 0 )
        {
          this_side_elements.insert( &ele );
        }
      }
    }

    std::swap( side_elements, this_side_elements );
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::SideIterator::find_other_element_side( stk::mesh::Entity & oe )
{
  // two elements that face each other
  stk::mesh::PairIterRelation onodes = oe.relations( stk::mesh::Node );

  const CellTopologyData * ocelltopology = get_cell_topology( oe.bucket() );
  unsigned dim = ocelltopology->dimension - 1;

  if ( ei_.topology().dimension - 1 != dim )
  {
    throw std::runtime_error( "expect neighbour elements of same dimension" );
  }

  const CellTopologyData_Subcell & side = ei_.topology().subcell[dim][side_num_];

  for ( unsigned oi=0; oi<ocelltopology->subcell_count[dim]; ++oi )
  {
    const CellTopologyData_Subcell & oside = ocelltopology->subcell[dim][oi] ;
    const CellTopologyData * oside_topology = oside.topology;

    if ( side.topology->node_count == oside_topology->node_count )
    {
      if ( contains_side( oside, onodes, side_nodes_ ) )
      {
        // at this point we have both elements and their respective
        // side numbers

        side_elements_.push_back( &oe );
        other_side_num_.push_back( oi );
        other_topology_.push_back( ocelltopology );

        return;
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::SideIterator::find_other_hanging_sides()
{
  EntitySet done_elements;

  // the main element is not considered any longer
  done_elements.insert( &ei_.element() );

  for ( EntitySet::iterator nodes = side_nodes_.begin(); nodes!=side_nodes_.end(); ++nodes )
  {
    stk::mesh::Entity * n = *nodes;
    for ( stk::mesh::PairIterRelation elements = n->relations( stk::mesh::Element );
          not elements.empty();
          ++elements )
    {
      stk::mesh::Entity & ele = *elements->entity();
      if ( ei_.selector()( ele.bucket() ) and
           done_elements.count( &ele )==0 )
      {
        done_elements.insert( &ele );

        stk::mesh::PairIterRelation nodes = ele.relations( stk::mesh::Node );

        const CellTopologyData * celltopology = get_cell_topology( ele.bucket() );
        unsigned dim = celltopology->dimension;

        if ( ei_.topology().dimension == dim )
        {
          for ( unsigned i=0; i<celltopology->subcell_count[dim-1]; ++i )
          {
            const CellTopologyData_Subcell & side = celltopology->subcell[dim-1][i] ;

            if ( contains_side( side, nodes, side_nodes_ ) )
            {
              side_elements_.push_back( &ele );
              other_side_num_.push_back( i );
              other_topology_.push_back( celltopology );
              break;
            }
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STK::SideIterator::contains_side( const CellTopologyData_Subcell & side,
                                       const stk::mesh::PairIterRelation & nodes,
                                       const EntitySet & side_nodes )
{
  const CellTopologyData * side_topology = side.topology;

  for ( unsigned j=0; j<side_topology->node_count; ++j )
  {
    stk::mesh::Entity * n = nodes[side.node[j]].entity();
    if ( side_nodes.count( n )==0 )
    {
      return false;
    }
  }
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::PrintMesh( std::ostream & stream, stk::mesh::BulkData & bulk, stk::mesh::Part & active )
{
  for ( ElementIterator i( bulk, active );
        not i.done();
        ++i )
  {
    stk::mesh::Entity & e = i.element();
    stream << "e" << e.key().id() << ": { ";

    for ( SideIterator j( i ); not j.done(); ++j )
    {
      unsigned side_number = j.side_number();
      stream << "s" << side_number << ":";

      if ( j.is_boundary_side() )
      {
        stream << "b";
      }
      else if ( j.is_inner_side() )
      {
        stk::mesh::Entity & oe = j.other_element();
        unsigned other_side_number = j.other_side_number();
        stream << "e" << oe.key().id() << "-s" << other_side_number;
      }
      else if ( j.is_small_hanging_side() )
      {
        stk::mesh::Entity & oe = j.other_element();
        unsigned other_side_number = j.other_side_number();
        stream << "h-e" << oe.key().id() << "-s" << other_side_number;
      }
      else if ( j.is_big_hanging_side() )
      {
        const EntityVector & side_elements = j.side_elements();
        const std::vector<unsigned> & other_side_numbers = j.other_side_numbers();
        if ( side_elements.size() != other_side_numbers.size() )
          throw std::runtime_error( "size mismatch" );
        stream << "( ";
        for ( unsigned k=0; k<side_elements.size(); ++k )
        {
          stk::mesh::Entity & oe = *side_elements[k];
          stream << "h-e" << oe.key().id() << "-s" << other_side_numbers[k] << " ";
        }
        stream << ")";
      }
      else
      {
        throw std::runtime_error( "side type undecided" );
      }
      stream << " ";
    }
    stream << "}\n";
  }
}

#endif
