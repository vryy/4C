#ifdef STKADAPTIVE

#include <algorithm>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include <boost/functional.hpp>
#include <boost/bind.hpp>

#include "stk_utils.H"
#include "stk_mesh.H"
#include "stk_topology.H"
#include "stk_comm.H"

namespace STK
{

void FindNeighborElements(const std::vector<const stk::mesh::Entity*>& nodes,
                          EntitySet & neighbors)
{
  FindNeighbors( nodes, neighbors, stk::mesh::Element );
}


stk::mesh::Entity* FindConstraint( const std::vector<const stk::mesh::Entity*>& nodes )
{
  stk::mesh::Entity* constraint = NULL;

  EntitySet neighbors;
  FindNeighbors( nodes, neighbors, stk::mesh::Constraint );

  // remove constraints that cover more that the given nodes

  for ( EntitySet::iterator i=neighbors.begin(); i!=neighbors.end(); )
  {
    stk::mesh::Entity * c = *i;
    stk::mesh::PairIterRelation rel = c->relations( stk::mesh::Node );
    if ( rel.size() > nodes.size()+1 )
    {
      neighbors.erase( i++ );
    }
    else
    {
      ++i;
    }
  }

  if ( neighbors.size()>1 )
    throw std::runtime_error("more that one constraint on node set");

  if ( neighbors.size()==1 )
  {
    constraint = *neighbors.begin();

    // If the set of nodes contains the constraint's hanging node, this is not
    // the constraint we are looking for.
    stk::mesh::PairIterRelation rel = constraint->relations( stk::mesh::Node );
    if ( not rel.empty() )
    {
      stk::mesh::Entity* hanging = rel[0].entity();
      if ( std::find( nodes.begin(), nodes.end(), hanging )!=nodes.end() )
        return NULL;
    }
  }
  return constraint;
}


void FindNeighbors( const std::vector<const stk::mesh::Entity*>& nodes,
                    EntitySet & neighbors,
                    unsigned entity_rank )
{
  if ( nodes.size()==0 )
    return;

  std::vector<const stk::mesh::Entity*>::const_iterator iter=nodes.begin();

  stk::mesh::PairIterRelation ir = (*iter)->relations( entity_rank );
  std::transform(ir.begin(),ir.end(),
                 std::inserter(neighbors,neighbors.begin()),
                 boost::bind(&stk::mesh::Relation::entity, _1));

  for ( ++iter ; iter!=nodes.end() and neighbors.size()>0 ; ++iter )
  {
    // get all inactive elements connected to node

    EntitySet nset;

    ir = (*iter)->relations( entity_rank );
    std::transform(ir.begin(),ir.end(),
                   std::inserter(nset,nset.begin()),
                   boost::bind(&stk::mesh::Relation::entity, _1));

    // get all inactive elements connected to all nodes from the set

    EntitySet intersection;

    std::set_intersection( nset.begin(), nset.end(),
                           neighbors.begin(), neighbors.end(),
                           std::inserter( intersection, intersection.begin() ),
                           EntitySet::key_compare() );

    std::swap(neighbors,intersection);
  }
}


void FindNeighbors( stk::mesh::Entity * e,
                    EntitySet & neighbors,
                    unsigned entity_rank )
{
  const CellTopologyData * celltopology = get_cell_topology( e->bucket() );
  stk::mesh::PairIterRelation nodes = e->relations( stk::mesh::Node );

  unsigned subcell_rank = celltopology->dimension - 1;
  for ( unsigned nitr = 0; nitr < celltopology->subcell_count[subcell_rank]; ++nitr )
  {
    // For the potentially common subcell, get it's nodes and num_nodes
    const unsigned* sub_nodes = celltopology->subcell[subcell_rank][nitr].node;
    unsigned num_nodes = celltopology->subcell[subcell_rank][nitr].topology->node_count;

    std::vector<const stk::mesh::Entity*> node_entities;
    node_entities.reserve( num_nodes );
    for ( unsigned itr = 0; itr < num_nodes; ++itr )
    {
      stk::mesh::Entity * n = nodes[sub_nodes[itr]].entity();
      node_entities.push_back( n );
    }

    EntitySet n;
    FindNeighbors( node_entities, n, entity_rank );

    std::copy( n.begin(), n.end(), std::inserter( neighbors, neighbors.begin() ) );
  }

  if ( entity_rank==e->entity_rank() )
    neighbors.erase( e );
}


void FindCommonEntities( stk::mesh::Entity* e, EntitySet & entities, stk::mesh::EntityRank rank )
{
  stk::mesh::PairIterRelation nodes = e->relations( stk::mesh::Node );

  stk::mesh::Entity * n = nodes->entity();
  stk::mesh::PairIterRelation rel = n->relations( rank );

  std::transform( rel.begin(), rel.end(),
                  std::inserter( entities, entities.begin() ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );

  for ( ++nodes; not nodes.empty() and entities.size()>0; ++nodes )
  {
    n = nodes->entity();
    rel = n->relations( rank );

    EntitySet es;
    std::transform( rel.begin(), rel.end(),
                    std::inserter( es, es.begin() ),
                    boost::bind( &stk::mesh::Relation::entity, _1 ) );

    EntitySet intersection;
    std::set_intersection( es.begin(), es.end(),
                           entities.begin(), entities.end(),
                           std::inserter( intersection, intersection.begin() ),
                           EntitySet::key_compare() );

    std::swap( entities, intersection );
  }

  if ( rank==e->entity_rank() )
    entities.erase( e );
}


void FindCommonEntities( const EntitySet & nodes, EntitySet & entities, stk::mesh::EntityRank rank )
{
  if ( nodes.size()==0 )
    return;

  EntitySet::const_iterator i=nodes.begin();

  stk::mesh::Entity * n = *i;
  stk::mesh::PairIterRelation rel = n->relations( rank );

  std::transform( rel.begin(), rel.end(),
                  std::inserter( entities, entities.begin() ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );

  for ( ++i; i!=nodes.end() and entities.size()>0; ++i )
  {
    n = *i;
    rel = n->relations( rank );

    EntitySet es;
    std::transform( rel.begin(), rel.end(),
                    std::inserter( es, es.begin() ),
                    boost::bind( &stk::mesh::Relation::entity, _1) );

    EntitySet intersection;
    std::set_intersection( es.begin(), es.end(),
                           entities.begin(), entities.end(),
                           std::inserter( intersection, intersection.begin() ),
                           EntitySet::key_compare() );

    std::swap( entities, intersection );
  }
}


void FindCommonNodeSharing(stk::mesh::Entity* element,
                           std::set<unsigned> & common_sharing)
{
  common_sharing.clear();

  stk::mesh::PairIterRelation ir = element->relations( stk::mesh::Node );
  if ( ir.empty() )
    return;

  // Include ghosted node as well. This way the connected elements will be
  // communicated to the aura.

  std::vector<unsigned> procs;
  stk::mesh::comm_procs( *ir->entity() , procs );
  std::copy( procs.begin(),
             procs.end(),
             std::inserter( common_sharing, common_sharing.begin() ) );
                  //boost::bind(&stk::mesh::PairIterEntityComm::value_type::proc,_1 ) );

  for ( ++ir ; not ir.empty() and common_sharing.size()>0 ; ++ir )
  {
    std::set<unsigned> mysharing;
    stk::mesh::comm_procs( *ir->entity() , procs );
    std::copy( procs.begin(),
               procs.end(),
               std::inserter( mysharing, mysharing.begin() ) );
                    //boost::bind(&stk::mesh::PairIterEntityComm::value_type::proc,_1 ) );

    std::set<unsigned> intersection;
    std::set_intersection(mysharing.begin(),mysharing.end(),
                          common_sharing.begin(),common_sharing.end(),
                          std::inserter(intersection,intersection.begin()));

    std::swap(common_sharing,intersection);
  }
}


void FindCommonNodeSharing(const std::vector<const stk::mesh::Entity*> & nodes,
                           std::set<unsigned> & common_sharing)
{
  common_sharing.clear();

  std::vector<const stk::mesh::Entity*>::const_iterator i=nodes.begin();
  if ( i==nodes.end() )
    return;

  const stk::mesh::PairIterEntityComm & sharing = ( *i )->sharing();
  std::transform(sharing.begin(),
                 sharing.end(),
                 std::inserter(common_sharing,common_sharing.begin()),
                 boost::bind(&stk::mesh::PairIterEntityComm::value_type::proc,_1));

  for ( ++i; i!=nodes.end(); ++i )
  {
    std::set<unsigned> mysharing;
    const stk::mesh::PairIterEntityComm & sharing = ( *i )->sharing();
    std::transform(sharing.begin(),
                   sharing.end(),
                   std::inserter(mysharing,mysharing.begin()),
                   boost::bind(&stk::mesh::PairIterEntityComm::value_type::proc,_1));

    std::set<unsigned> intersection;
    std::set_intersection(mysharing.begin(),mysharing.end(),
                          common_sharing.begin(),common_sharing.end(),
                          std::inserter(intersection,intersection.begin()));

    std::swap(common_sharing,intersection);
  }
}


stk::mesh::Entity* InactiveParent( stk::mesh::Part & active,
                                   const EntitySet & nodes,
                                   stk::mesh::EntityRank parent_rank )
{
  if ( nodes.size()==0 )
  {
    return NULL;
  }

  stk::mesh::Selector activesel = active;
  EntitySet::iterator i=nodes.begin();

  stk::mesh::PairIterRelation rn = (*i)->relations( parent_rank );
  if ( rn.empty() )
  {
    return NULL;
  }

  EntitySet elements;
  std::transform(rn.begin(),rn.end(),
                 std::inserter(elements,elements.begin()),
                 boost::bind(&stk::mesh::Relation::entity, _1));

  for ( ++i; i!=nodes.end(); ++i )
  {
    EntitySet new_elements;
    stk::mesh::Entity * n = *i;

    rn = n->relations( parent_rank );
    std::transform(rn.begin(),rn.end(),
                   std::inserter(new_elements,new_elements.begin()),
                   boost::bind(&stk::mesh::Relation::entity, _1));

    EntitySet intersection;
    std::set_intersection( new_elements.begin(), new_elements.end(),
                           elements.begin(), elements.end(),
                           std::inserter( intersection, intersection.begin() ),
                           EntitySet::key_compare() );

    std::swap( intersection, elements );

    // remove active elements since we are looking for the inactive parent
    for ( EntitySet::iterator j=elements.begin();
          j!=elements.end(); )
    {
      if ( activesel( ( *j )->bucket() ) )
      {
        elements.erase( j++ );
      }
      else
      {
        ++j;
      }
    }

    if ( elements.size()==0 )
    {
      return NULL;
    }
  }

  if ( elements.size()!=1 )
  {
    throw std::runtime_error("only one parent element");
  }
  return *elements.begin();
}


void FindCommonParts( const std::vector<stk::mesh::Entity*> & nodes,
                      stk::mesh::PartVector & parts )
{
  if ( nodes.size()==0 )
    return;

  std::vector<stk::mesh::Entity*>::const_iterator i=nodes.begin();
  std::set<unsigned> partids;
  std::pair<const unsigned *, const unsigned *> p = ( *i )->bucket().superset_part_ordinals();
  std::copy( p.first, p.second, std::inserter( partids, partids.begin() ) );

  for ( ++i ; i!=nodes.end() ; ++i )
  {
    std::set<unsigned> nodepartids;
    p = ( *i )->bucket().superset_part_ordinals();
    std::copy( p.first, p.second, std::inserter( nodepartids, nodepartids.begin() ) );

    std::set<unsigned> intersection;
    std::set_intersection( partids.begin(), partids.end(),
                           nodepartids.begin(), nodepartids.end(),
                           std::inserter( intersection, intersection.begin() ) );

    std::swap( intersection, partids );
  }

  // remove universal, ownes and uses
  for ( unsigned j=0; j<3; ++j )
  {
    std::set<unsigned>::iterator pi=partids.find( j );
    if ( pi!=partids.end() )
    {
      partids.erase( pi );
    }
  }

  const stk::mesh::MetaData & meta_data = nodes[0]->bucket().mesh().mesh_meta_data();
  parts.clear();
  for ( std::set<unsigned>::iterator pi=partids.begin();
        pi!=partids.end();
        ++pi )
  {
    stk::mesh::Part * part = & meta_data.get_part( *pi );
    if ( part->primary_entity_rank()==stk::mesh::Node or
         part->primary_entity_rank()==std::numeric_limits<unsigned>::max() )
    {
      parts.push_back( part );
    }
  }
}


bool IsChildListComplete( stk::mesh::Entity* e,
                          const EntitySet& children )
{
  EntitySet parentnodes;
  stk::mesh::PairIterRelation ip = e->relations( stk::mesh::Node );

  std::transform(ip.begin(),ip.end(),
                 std::inserter(parentnodes,parentnodes.begin()),
                 boost::bind( &stk::mesh::Relation::entity, _1 ) );

  // If parent and children use the same set of nodes, we have indeed all
  // children of the given parent. This requires that each children introduces
  // at least one new node to its parent's set of nodes.

  EntitySet nodes;

  for ( EntitySet::const_iterator i=children.begin(); i!=children.end(); ++i )
  {
    stk::mesh::Entity & e = **i;
    stk::mesh::PairIterRelation rel = e.relations( stk::mesh::Node );

    std::transform( rel.begin(), rel.end(),
                    std::inserter( nodes, nodes.begin() ),
                    boost::bind( &stk::mesh::Relation::entity, _1 ) );
  }

#if 0
  std::copy( childnodes.begin(), childnodes.end(),
             std::inserter( nodes, nodes.begin() ) );
  std::copy( hangingnodes.begin(), hangingnodes.end(),
             std::inserter( nodes, nodes.begin() ) );
#endif

  return parentnodes==nodes;
}


bool is_child( stk::mesh::Entity * child, stk::mesh::Entity * parent )
{
  // if the cild nodes are a subset of the parent nodes the child is actually
  // a child

  stk::mesh::PairIterRelation child_nodes  = child->relations( stk::mesh::Node );
  stk::mesh::PairIterRelation parent_nodes = parent->relations( stk::mesh::Node );

  EntitySet es;
  std::transform( parent_nodes.begin(), parent_nodes.end(),
                  std::inserter( es, es.begin() ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );
  std::transform( child_nodes.begin(), child_nodes.end(),
                  std::inserter( es, es.begin() ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );

  return es.size()==parent_nodes.size();
}


stk::mesh::Entity & AddNode( Mesh & mesh, int nid, double* x )
{
  stk::mesh::Entity * n = mesh.BulkData().get_entity( stk::mesh::Node, nid );
  if ( n!=NULL )
  {
#if 1
    double * const coord = stk::mesh::field_data( mesh.Coordinates() , *n );

    // this is a hack since now entities are always created beforehand
    if ( coord[0]==0 and coord[1]==0 and coord[2]==0 )
    {
      coord[0] = x[0];
      coord[1] = x[1];
      coord[2] = x[2];
    }
    else
    {
      double dx = coord[0] - x[0];
      double dy = coord[1] - x[1];
      double dz = coord[2] - x[2];

      if ( dx*dx+dy*dy+dz*dz > 1e-6 )
      {
        std::ostringstream str;
        str << "create node("
            << nid << "," << x[0] << "," << x[1] << "," << x[2]
            << ") failed because node("
            << nid << "," << coord[0] << "," << coord[1] << "," << coord[2]
            << ") exists";
        throw std::runtime_error( str.str() );
      }
    }
#endif
    return *n;
  }
  else
  {
    stk::mesh::PartVector empty ;
    stk::mesh::Entity & node = mesh.BulkData().declare_entity( stk::mesh::Node, nid, empty );
#if 1
    double * const coord = stk::mesh::field_data( mesh.Coordinates() , node );

    coord[0] = x[0];
    coord[1] = x[1];
    coord[2] = x[2];
#endif
    return node;
  }
}


void MiddlePoint( Mesh & mesh, stk::mesh::Entity* nid1, stk::mesh::Entity* nid2, double* x )
{
  stk::mesh::VectorField & coordinates = mesh.Coordinates();

  double * const c1 = stk::mesh::field_data( coordinates , *nid1 );
  double * const c2 = stk::mesh::field_data( coordinates , *nid2 );

  x[0] = .5*(c1[0] + c2[0]);
  x[1] = .5*(c1[1] + c2[1]);
  x[2] = .5*(c1[2] + c2[2]);
}


}

#endif
