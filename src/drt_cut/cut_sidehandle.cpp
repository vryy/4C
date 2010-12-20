
#include "cut_sidehandle.H"
#include "cut_mesh.H"
#include "cut_node.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


GEO::CUT::Tri6SideHandle::Tri6SideHandle( Mesh & mesh, int sid, const std::vector<int> & nodes )
{
  subsides_.reserve( 4 );

  nodes_.reserve( 4 );
  for ( int i=0; i<4; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], NULL );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Triangle<3> >();

  std::vector<int> nids( 3 );

  nids[0] = nodes[0];
  nids[1] = nodes[3];
  nids[2] = nodes[5];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[3];
  nids[1] = nodes[1];
  nids[2] = nodes[4];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[3];
  nids[1] = nodes[4];
  nids[2] = nodes[5];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[5];
  nids[1] = nodes[4];
  nids[2] = nodes[2];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );
}

GEO::CUT::Quad8SideHandle::Quad8SideHandle( Mesh & mesh, int sid, const std::vector<int> & nodes )
{
  subsides_.reserve( 4 );

  LINALG::Matrix<3,8> xyze;
  nodes_.reserve( 8 );
  for ( int i=0; i<8; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], NULL );
    nodes_.push_back( n );
    n->Coordinates( &xyze( 0, i ) );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();

  // create middle node

  LINALG::Matrix<8,1> funct;
  DRT::UTILS::shape_function_2D( funct, 0, 0, DRT::Element::quad8 );

  LINALG::Matrix<3,1> xyz;
  xyz.Multiply( xyze, funct );

  std::set<int> node_nids;
  std::copy( nodes.begin(), nodes.end(), std::inserter( node_nids, node_nids.begin() ) );
  Node* middle = mesh.GetNode( node_nids, xyz.A() );
  int middle_id = middle->Id();

  std::vector<int> nids( 4 );

  nids[0] = nodes[0];
  nids[1] = nodes[4];
  nids[2] = middle_id;
  nids[3] = nodes[7];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[4];
  nids[1] = nodes[1];
  nids[2] = nodes[5];
  nids[3] = middle_id;
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = middle_id;
  nids[1] = nodes[5];
  nids[2] = nodes[2];
  nids[3] = nodes[6];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[7];
  nids[1] = middle_id;
  nids[2] = nodes[6];
  nids[3] = nodes[3];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );
}

GEO::CUT::Quad9SideHandle::Quad9SideHandle( Mesh & mesh, int sid, const std::vector<int> & nodes )
{
  subsides_.reserve( 4 );

  nodes_.reserve( 9 );
  for ( int i=0; i<9; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], NULL );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();

  std::vector<int> nids( 4 );

  nids[0] = nodes[0];
  nids[1] = nodes[4];
  nids[2] = nodes[8];
  nids[3] = nodes[7];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[4];
  nids[1] = nodes[1];
  nids[2] = nodes[5];
  nids[3] = nodes[8];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[8];
  nids[1] = nodes[5];
  nids[2] = nodes[2];
  nids[3] = nodes[6];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );

  nids[0] = nodes[7];
  nids[1] = nodes[8];
  nids[2] = nodes[6];
  nids[3] = nodes[3];
  subsides_.push_back( mesh.GetSide( sid, nids, top_data ) );
}
