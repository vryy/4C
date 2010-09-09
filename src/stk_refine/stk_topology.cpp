#ifdef STKADAPTIVE

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>

#include <boost/bind.hpp>

#include "stk_topology.H"
#include "stk_refine.H"
#include "stk_mesh.H"
#include "stk_utils.H"

namespace STK
{

const TopologyType<void>* TopologyType<void>::GetTopologyType(const stk::mesh::Entity* e)
{
  const CellTopologyData * top = get_cell_topology( e->bucket() );
  if (top==NULL)
  {
    std::ostringstream msg ;
    msg << "Entity " << e->entity_rank() << "," << e->identifier()
        << " has no topology.";
    throw std::runtime_error( msg.str() );
  }

  switch ( top->key )
  {
  case shards::Line<2>::key :
  {
    static TopologyType<shards::Line<2> > tt;
    return &tt;
  }

  case shards::Quadrilateral<4>::key :
  {
    static TopologyType<shards::Quadrilateral<4> > tt;
    return &tt;
  }

  case shards::Hexahedron<8>::key :
  {
    static TopologyType<shards::Hexahedron<8> > tt;
    return &tt;
  }

  case shards::Node::key :

  case shards::Line<3>::key :

  case shards::ShellLine<2>::key :

  case shards::ShellLine<3>::key :

  case shards::Triangle<3>::key :

  case shards::Triangle<6>::key :

  case shards::ShellTriangle<3>::key :

  case shards::ShellTriangle<6>::key :

  case shards::Quadrilateral<8>::key :

  case shards::Quadrilateral<9>::key :

  case shards::ShellQuadrilateral<4>::key :

  case shards::ShellQuadrilateral<8>::key :

  case shards::ShellQuadrilateral<9>::key :

  case shards::Tetrahedron< 4>::key :

  case shards::Tetrahedron<10>::key :

  case shards::Pyramid< 5>::key :

  case shards::Pyramid<13>::key :

  case shards::Pyramid<14>::key :

  case shards::Wedge< 6>::key :

  case shards::Wedge<15>::key :

  case shards::Wedge<18>::key :

  case shards::Hexahedron<20>::key :

  case shards::Hexahedron<27>::key :
    break;
  }

  throw std::runtime_error( "entity topology not supported" );
}


void TopologyType<shards::Line<2> >::FindChildNodes(RefineSet& rs,
                                                    const stk::mesh::Entity* e) const
{
  stk::mesh::PairIterRelation nodes = e->relations( stk::mesh::Node );
  std::vector<const stk::mesh::Entity*> n;
  n.reserve( 2 );
  std::transform( nodes.begin(), nodes.end(),
                  std::back_inserter( n ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );
  rs.CreateHangingNode( e, n );
}


void TopologyType<shards::Line<2> >::HangingNodes( const RefineSet* rs,
                                                   Mesh* mesh,
                                                   const std::vector<stk::mesh::Entity*> & nodes,
                                                   std::vector<stk::mesh::Entity*> & hanging ) const
{
  hanging.resize( 1 );
  double x[3];

  MiddlePoint( *mesh, nodes[0], nodes[1], x );
  hanging[0] = & AddNode( *mesh, rs->NodeId( nodes[0],
                                             nodes[1] ), x );

  if ( hanging[0]->owner_rank() == mesh->BulkData().parallel_rank() )
  {
    stk::mesh::PartVector add_parts;
    stk::mesh::PartVector remove_parts;

    FindCommonParts( nodes, add_parts );
    mesh->BulkData().change_entity_parts( *hanging[0], add_parts, remove_parts );
  }
}


void TopologyType<shards::Line<2> >::CreateElements( Mesh* mesh,
                                                     stk::mesh::EntityRank rank,
                                                     const std::pair<const unsigned *, const unsigned *> & part_ord,
                                                     std::vector<stk::mesh::Entity*> & nodes,
                                                     std::vector<stk::mesh::Entity*> & hanging,
                                                     std::vector<stk::mesh::Entity *>::iterator & elements_iter ) const
{
  stk::mesh::BulkData & bulk = mesh->BulkData();
  stk::mesh::MetaData & meta = mesh->MetaData();

//   stk::mesh::Part & active = mesh->ActivePart();
//   stk::mesh::Part & line2  = mesh->LinePart();

  std::vector<stk::mesh::Part*> add_parts;
  std::vector<stk::mesh::Part*> remove_parts;

  add_parts.reserve( part_ord.second - part_ord.first );

//   add_parts.push_back( &active );
//   add_parts.push_back( &line2 );

  for ( const unsigned * i=part_ord.first; i!=part_ord.second; ++i )
  {
    // ignore universal, uses and owns
    if ( *i > 2 )
    {
      stk::mesh::Part * p = & meta.get_part( *i );
      unsigned primary_rank = p->primary_entity_rank();
      if ( rank==primary_rank or primary_rank==std::numeric_limits<unsigned>::max() )
        add_parts.push_back( p );
    }
  }

  // declare elements

  stk::mesh::Entity * e;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *nodes[0]   , 0 );
  bulk.declare_relation( *e , *hanging[0] , 1 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *hanging[0] , 0 );
  bulk.declare_relation( *e , *nodes[1]   , 1 );

  ++elements_iter;
}


// =================================================================

void TopologyType<shards::Quadrilateral<4> >::FindChildNodes(RefineSet& rs,
                                                             const stk::mesh::Entity* e) const
{
  // We know the element that needs to be refined and are looking for the new
  // nodes. These can be hanging nodes that need to be created as well as
  // existing nodes that might even cease to be hanging nodes.

  stk::mesh::PairIterRelation nodes = e->relations( stk::mesh::Node );

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();

  std::vector<const stk::mesh::Entity*> n( 2 );
  for ( unsigned i=0; i<top_data->side_count; ++i )
  {
    const CellTopologyData_Subcell & side = top_data->side[i] ;

    n[0] = nodes[side.node[0]].entity();
    n[1] = nodes[side.node[1]].entity();

    rs.CreateHangingNode( e, n );
  }

  // middle node

  n.clear();
  std::transform( nodes.begin(), nodes.end(),
                  std::back_inserter( n ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );
  rs.CreateHangingNode( e, n );
}


void TopologyType<shards::Quadrilateral<4> >::HangingNodes(
  const RefineSet* rs,
  Mesh* mesh,
  const std::vector<stk::mesh::Entity*> & nodes,
  std::vector<stk::mesh::Entity*> & hanging ) const
{
  hanging.resize(5);
  double x[3];

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();

  stk::mesh::PartVector remove_parts;

  std::vector<stk::mesh::Entity*> n( 2 );

  // implicit definition of hanging node topology

  for ( unsigned i=0; i<top_data->side_count; ++i )
  {
    const CellTopologyData_Subcell & side = top_data->side[i] ;

    n[0] = nodes[side.node[0]];
    n[1] = nodes[side.node[1]];

    MiddlePoint( *mesh, n[0], n[1], x );
    hanging[i] = & AddNode( *mesh, rs->NodeId( n[0], n[1] ), x );

    if ( hanging[i]->owner_rank() == mesh->BulkData().parallel_rank() )
    {
      stk::mesh::PartVector add_parts;
      FindCommonParts( n, add_parts );
      mesh->BulkData().change_entity_parts( *hanging[i], add_parts, remove_parts );
    }
  }

  MiddlePoint( *mesh, hanging[0],hanging[2],x);
  hanging[4] = & AddNode( *mesh, rs->NodeId(nodes[0],
                                            nodes[1],
                                            nodes[2],
                                            nodes[3]), x );
  if ( hanging[4]->owner_rank() == mesh->BulkData().parallel_rank() )
  {
    stk::mesh::PartVector add_parts;
    FindCommonParts( nodes, add_parts );
    mesh->BulkData().change_entity_parts( *hanging[4], add_parts, remove_parts );
  }
}


void TopologyType<shards::Quadrilateral<4> >::CreateElements(
  Mesh* mesh,
  stk::mesh::EntityRank rank,
  const std::pair<const unsigned *, const unsigned *> & part_ord,
  std::vector<stk::mesh::Entity*> & nodes,
  std::vector<stk::mesh::Entity*> & hanging,
  std::vector<stk::mesh::Entity*>::iterator & elements_iter ) const
{
  stk::mesh::BulkData & bulk = mesh->BulkData();
  stk::mesh::MetaData & meta = mesh->MetaData();

//   stk::mesh::Part & active = mesh->ActivePart();
//   stk::mesh::Part & quad4  = mesh->QuadPart();

  std::vector<stk::mesh::Part*> add_parts;
  std::vector<stk::mesh::Part*> remove_parts;

  add_parts.reserve( part_ord.second - part_ord.first );

//   add_parts.push_back( &active );
//   add_parts.push_back( &quad4 );

  for ( const unsigned * i=part_ord.first; i!=part_ord.second; ++i )
  {
    // ignore universal, uses and owns
    if ( *i > 2 )
    {
      stk::mesh::Part * p = & meta.get_part( *i );
      unsigned primary_rank = p->primary_entity_rank();
      if ( rank==primary_rank or primary_rank==std::numeric_limits<unsigned>::max() )
        add_parts.push_back( p );
    }
  }

  // declare elements

  stk::mesh::Entity * e;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *nodes[0]   , 0 );
  bulk.declare_relation( *e , *hanging[0] , 1 );
  bulk.declare_relation( *e , *hanging[4] , 2 );
  bulk.declare_relation( *e , *hanging[3] , 3 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *hanging[0] , 0 );
  bulk.declare_relation( *e , *nodes[1]   , 1 );
  bulk.declare_relation( *e , *hanging[1] , 2 );
  bulk.declare_relation( *e , *hanging[4] , 3 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *hanging[4] , 0 );
  bulk.declare_relation( *e , *hanging[1] , 1 );
  bulk.declare_relation( *e , *nodes[2]   , 2 );
  bulk.declare_relation( *e , *hanging[2] , 3 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *hanging[3] , 0 );
  bulk.declare_relation( *e , *hanging[4] , 1 );
  bulk.declare_relation( *e , *hanging[2] , 2 );
  bulk.declare_relation( *e , *nodes[3]   , 3 );

  ++elements_iter;
}


// =================================================================

void TopologyType<shards::Hexahedron<8> >::FindChildNodes(RefineSet& rs,
                                                          const stk::mesh::Entity* e) const
{
  // We know the element that needs to be refined and are looking for the new
  // nodes. These can be hanging nodes that need to be created as well as
  // existing nodes that might even cease to be hanging nodes.

  stk::mesh::PairIterRelation nodes = e->relations( stk::mesh::Node );

  // Order is important here!
  // This method can be generalized.

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  // sides

  std::vector<const stk::mesh::Entity*> n( 4 );
  for ( unsigned i=0; i<top_data->side_count; ++i )
  {
    const CellTopologyData_Subcell & side = top_data->side[i] ;

    n[0] = nodes[side.node[0]].entity();
    n[1] = nodes[side.node[1]].entity();
    n[2] = nodes[side.node[2]].entity();
    n[3] = nodes[side.node[3]].entity();

    rs.CreateHangingNode( e, n );
  }

  // edges

  n.resize( 2 );
  for ( unsigned i=0; i<top_data->edge_count; ++i )
  {
    const CellTopologyData_Subcell & edge = top_data->edge[i] ;

    n[0] = nodes[edge.node[0]].entity();
    n[1] = nodes[edge.node[1]].entity();

    rs.CreateHangingNode( e, n );
  }

  // middle node

  n.clear();
  n.reserve( 8 );
  std::transform( nodes.begin(), nodes.end(),
                  std::back_inserter( n ),
                  boost::bind( &stk::mesh::Relation::entity, _1 ) );
  rs.CreateHangingNode( e, n );
}


void TopologyType<shards::Hexahedron<8> >::HangingNodes(
  const RefineSet* rs,
  Mesh* mesh,
  const std::vector<stk::mesh::Entity*> & nodes,
  std::vector<stk::mesh::Entity*> & hanging ) const
{
  hanging.resize(19);
  double x[3];

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  stk::mesh::PartVector remove_parts;

  std::vector<stk::mesh::Entity*> n( 4 );

  // implicit definition of hanging node topology

  unsigned sc = top_data->side_count;
  unsigned ec = top_data->edge_count;

  for ( unsigned i=0; i<sc; ++i )
  {
    const CellTopologyData_Subcell & side = top_data->side[i] ;

    n[0] = nodes[side.node[0]];
    n[1] = nodes[side.node[1]];
    n[2] = nodes[side.node[2]];
    n[3] = nodes[side.node[3]];

    MiddlePoint( *mesh, n[0], n[2], x );
    hanging[i] = & AddNode( *mesh, rs->NodeId( n[0], n[1], n[2], n[3] ), x );

    if ( hanging[i]->owner_rank() == mesh->BulkData().parallel_rank() )
    {
      stk::mesh::PartVector add_parts;
      FindCommonParts( n, add_parts );
      mesh->BulkData().change_entity_parts( *hanging[i], add_parts, remove_parts );
    }
  }

  n.resize( 2 );

  for ( unsigned i=0; i<ec; ++i )
  {
    const CellTopologyData_Subcell & edge = top_data->edge[i] ;

    n[0] = nodes[edge.node[0]];
    n[1] = nodes[edge.node[1]];

    MiddlePoint( *mesh, n[0], n[1], x );
    hanging[i+sc] = & AddNode( *mesh, rs->NodeId( n[0], n[1] ), x );

    if ( hanging[i+sc]->owner_rank() == mesh->BulkData().parallel_rank() )
    {
      stk::mesh::PartVector add_parts;
      FindCommonParts( n, add_parts );
      mesh->BulkData().change_entity_parts( *hanging[i+sc], add_parts, remove_parts );
    }
  }

  MiddlePoint( *mesh, hanging[0],hanging[2],x );
  hanging[sc+ec] = & AddNode( *mesh, rs->NodeId(nodes[0],
                                                nodes[1],
                                                nodes[2],
                                                nodes[3],
                                                nodes[4],
                                                nodes[5],
                                                nodes[6],
                                                nodes[7]), x );
  if ( hanging[sc+ec]->owner_rank() == mesh->BulkData().parallel_rank() )
  {
    stk::mesh::PartVector add_parts;
    FindCommonParts( nodes, add_parts );
    mesh->BulkData().change_entity_parts( *hanging[sc+ec], add_parts, remove_parts );
  }
}


void TopologyType<shards::Hexahedron<8> >::CreateElements(
  Mesh* mesh,
  stk::mesh::EntityRank rank,
  const std::pair<const unsigned *, const unsigned *> & part_ord,
  std::vector<stk::mesh::Entity*> & nodes,
  std::vector<stk::mesh::Entity*> & hanging,
  std::vector<stk::mesh::Entity*>::iterator & elements_iter ) const
{
  stk::mesh::BulkData & bulk = mesh->BulkData();
  stk::mesh::MetaData & meta = mesh->MetaData();

//   stk::mesh::Part & active = mesh->ActivePart();
//   stk::mesh::Part & hex8   = mesh->HexPart();

  std::vector<stk::mesh::Part*> add_parts;
  std::vector<stk::mesh::Part*> remove_parts;

  add_parts.reserve( part_ord.second - part_ord.first );

//   add_parts.push_back( &active );
//   add_parts.push_back( &hex8 );

  for ( const unsigned * i=part_ord.first; i!=part_ord.second; ++i )
  {
    // ignore universal, uses and owns
    if ( *i > 2 )
    {
      stk::mesh::Part * p = & meta.get_part( *i );
      unsigned primary_rank = p->primary_entity_rank();
      if ( rank==primary_rank or primary_rank==std::numeric_limits<unsigned>::max() )
        add_parts.push_back( p );
    }
  }

  stk::mesh::Entity ** side = &hanging[0];
  stk::mesh::Entity ** edge = &hanging[6];
  stk::mesh::Entity * center = hanging[6+12];

  std::vector<stk::mesh::Entity*> n;
  n.reserve( 27 );
  std::copy( nodes.begin(), nodes.end(), std::back_inserter( n ) );

  n.push_back( edge[ 0] );
  n.push_back( edge[ 1] );
  n.push_back( edge[ 2] );
  n.push_back( edge[ 3] );

  n.push_back( edge[ 8] );
  n.push_back( edge[ 9] );
  n.push_back( edge[10] );
  n.push_back( edge[11] );

  n.push_back( edge[4] );
  n.push_back( edge[5] );
  n.push_back( edge[6] );
  n.push_back( edge[7] );

  n.push_back( center );

  n.push_back( side[4] );
  n.push_back( side[5] );
  n.push_back( side[3] );
  n.push_back( side[1] );
  n.push_back( side[0] );
  n.push_back( side[2] );

  // declare elements

  stk::mesh::Entity * e;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[ 0] , 0 );
  bulk.declare_relation( *e , *n[ 8] , 1 );
  bulk.declare_relation( *e , *n[21] , 2 );
  bulk.declare_relation( *e , *n[11] , 3 );

  bulk.declare_relation( *e , *n[12] , 4 );
  bulk.declare_relation( *e , *n[25] , 5 );
  bulk.declare_relation( *e , *n[20] , 6 );
  bulk.declare_relation( *e , *n[23] , 7 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[ 8] , 0 );
  bulk.declare_relation( *e , *n[ 1] , 1 );
  bulk.declare_relation( *e , *n[ 9] , 2 );
  bulk.declare_relation( *e , *n[21] , 3 );

  bulk.declare_relation( *e , *n[25] , 4 );
  bulk.declare_relation( *e , *n[13] , 5 );
  bulk.declare_relation( *e , *n[24] , 6 );
  bulk.declare_relation( *e , *n[20] , 7 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[21] , 0 );
  bulk.declare_relation( *e , *n[ 9] , 1 );
  bulk.declare_relation( *e , *n[ 2] , 2 );
  bulk.declare_relation( *e , *n[10] , 3 );

  bulk.declare_relation( *e , *n[20] , 4 );
  bulk.declare_relation( *e , *n[24] , 5 );
  bulk.declare_relation( *e , *n[14] , 6 );
  bulk.declare_relation( *e , *n[26] , 7 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[11] , 0 );
  bulk.declare_relation( *e , *n[21] , 1 );
  bulk.declare_relation( *e , *n[10] , 2 );
  bulk.declare_relation( *e , *n[ 3] , 3 );

  bulk.declare_relation( *e , *n[23] , 4 );
  bulk.declare_relation( *e , *n[20] , 5 );
  bulk.declare_relation( *e , *n[26] , 6 );
  bulk.declare_relation( *e , *n[15] , 7 );

  ++elements_iter;

  // -----------------------------------------------------------------

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[12] , 0 );
  bulk.declare_relation( *e , *n[25] , 1 );
  bulk.declare_relation( *e , *n[20] , 2 );
  bulk.declare_relation( *e , *n[23] , 3 );

  bulk.declare_relation( *e , *n[ 4] , 4 );
  bulk.declare_relation( *e , *n[16] , 5 );
  bulk.declare_relation( *e , *n[22] , 6 );
  bulk.declare_relation( *e , *n[19] , 7 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[25] , 0 );
  bulk.declare_relation( *e , *n[13] , 1 );
  bulk.declare_relation( *e , *n[24] , 2 );
  bulk.declare_relation( *e , *n[20] , 3 );

  bulk.declare_relation( *e , *n[16] , 4 );
  bulk.declare_relation( *e , *n[ 5] , 5 );
  bulk.declare_relation( *e , *n[17] , 6 );
  bulk.declare_relation( *e , *n[22] , 7 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[20] , 0 );
  bulk.declare_relation( *e , *n[24] , 1 );
  bulk.declare_relation( *e , *n[14] , 2 );
  bulk.declare_relation( *e , *n[26] , 3 );

  bulk.declare_relation( *e , *n[22] , 4 );
  bulk.declare_relation( *e , *n[17] , 5 );
  bulk.declare_relation( *e , *n[ 6] , 6 );
  bulk.declare_relation( *e , *n[18] , 7 );

  ++elements_iter;

  e = *elements_iter;
  bulk.change_entity_parts( *e, add_parts, remove_parts );
  bulk.declare_relation( *e , *n[23] , 0 );
  bulk.declare_relation( *e , *n[20] , 1 );
  bulk.declare_relation( *e , *n[26] , 2 );
  bulk.declare_relation( *e , *n[15] , 3 );

  bulk.declare_relation( *e , *n[19] , 4 );
  bulk.declare_relation( *e , *n[22] , 5 );
  bulk.declare_relation( *e , *n[18] , 6 );
  bulk.declare_relation( *e , *n[ 7] , 7 );

  ++elements_iter;
}

}

#endif
