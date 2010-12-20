
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "cut_elementhandle.H"
#include "cut_mesh.H"
#include "cut_node.H"
#include "cut_position.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

bool GEO::CUT::QuadraticElementHandle::IsCut()
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    if ( e->IsCut() )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::QuadraticElementHandle::GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->GetIntegrationCells( cells );
  }
}

void GEO::CUT::QuadraticElementHandle::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->GetBoundaryCells( bcells );
  }
}


GEO::CUT::Hex20ElementHandle::Hex20ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nodes )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  // create middle nodes

  LINALG::Matrix<3,1> xyz;

  std::set<int> node_nids;

  LINALG::Matrix<3,8> side_xyze;
  LINALG::Matrix<8,1> side_funct;
  std::vector<Node*> side_nodes( 8 );

  std::vector<Node*> center_nodes( 6 );

  for ( int localsideid = 0; localsideid < 6; ++localsideid )
  {
    node_nids.clear();
    for ( int i=0; i<8; ++i )
    {
      int localnodeid = DRT::UTILS::eleNodeNumbering_hex27_surfaces[localsideid][i];
      Node * n = mesh.GetNode( nodes[localnodeid], NULL );
      side_nodes[i] = n;
      node_nids.insert( nodes[localnodeid] );
      n->Coordinates( &side_xyze( 0, i ) );
    }

    DRT::UTILS::shape_function_2D( side_funct, 0, 0, DRT::Element::quad8 );
    xyz.Multiply( side_xyze, side_funct );

    center_nodes[localsideid] = mesh.GetNode( node_nids, xyz.A() );
  }

  Node* node20 = center_nodes[0];
  int node20_id = node20->Id();

  Node* node21 = center_nodes[1];
  int node21_id = node21->Id();

  Node* node22 = center_nodes[2];
  int node22_id = node22->Id();

  Node* node23 = center_nodes[3];
  int node23_id = node23->Id();

  Node* node24 = center_nodes[4];
  int node24_id = node24->Id();

  Node* node25 = center_nodes[5];
  int node25_id = node25->Id();

  LINALG::Matrix<3,20> xyze;
  nodes_.reserve( 20 );
  for ( int i=0; i<20; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], NULL );
    nodes_.push_back( n );
    n->Coordinates( &xyze( 0,i ) );
  }

  LINALG::Matrix<20,1> funct;
  DRT::UTILS::shape_function_3D( funct, 0, 0, 0, DRT::Element::hex20 );

  xyz.Multiply( xyze, funct );
  node_nids.clear();
  std::copy( nodes.begin(), nodes.end(), std::inserter( node_nids, node_nids.begin() ) );
  Node* node26 = mesh.GetNode( node_nids, xyz.A() );
  int node26_id = node26->Id();


  std::vector<int> nids( 8 );

  nids[0] = nodes[ 0];
  nids[1] = nodes[ 8];
  nids[2] = node20_id;
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = node21_id;
  nids[6] = node26_id;
  nids[7] = node24_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[ 8];
  nids[1] = nodes[ 1];
  nids[2] = nodes[ 9];
  nids[3] = node20_id;
  nids[4] = node21_id;
  nids[5] = nodes[13];
  nids[6] = node22_id;
  nids[7] = node26_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node20_id;
  nids[1] = nodes[ 9];
  nids[2] = nodes[ 2];
  nids[3] = nodes[10];
  nids[4] = node26_id;
  nids[5] = node22_id;
  nids[6] = nodes[14];
  nids[7] = node23_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[11];
  nids[1] = node20_id;
  nids[2] = nodes[10];
  nids[3] = nodes[ 3];
  nids[4] = node24_id;
  nids[5] = node26_id;
  nids[6] = node23_id;
  nids[7] = nodes[15];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = node21_id;
  nids[2] = node26_id;
  nids[3] = node24_id;
  nids[4] = nodes[ 4];
  nids[5] = nodes[16];
  nids[6] = node25_id;
  nids[7] = nodes[19];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node21_id;
  nids[1] = nodes[13];
  nids[2] = node22_id;
  nids[3] = node26_id;
  nids[4] = nodes[16];
  nids[5] = nodes[ 5];
  nids[6] = nodes[17];
  nids[7] = node25_id;
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node26_id;
  nids[1] = node22_id;
  nids[2] = nodes[14];
  nids[3] = node23_id;
  nids[4] = node25_id;
  nids[5] = nodes[17];
  nids[6] = nodes[ 6];
  nids[7] = nodes[18];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = node24_id;
  nids[1] = node26_id;
  nids[2] = node23_id;
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = node25_id;
  nids[6] = nodes[18];
  nids[7] = nodes[ 7];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
}

GEO::CUT::Hex27ElementHandle::Hex27ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nodes )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  nodes_.reserve( 27 );
  for ( int i=0; i<27; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], NULL );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  std::vector<int> nids( 8 );

  nids[0] = nodes[ 0];
  nids[1] = nodes[ 8];
  nids[2] = nodes[20];
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = nodes[21];
  nids[6] = nodes[26];
  nids[7] = nodes[24];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[ 8];
  nids[1] = nodes[ 1];
  nids[2] = nodes[ 9];
  nids[3] = nodes[20];
  nids[4] = nodes[21];
  nids[5] = nodes[13];
  nids[6] = nodes[22];
  nids[7] = nodes[26];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[20];
  nids[1] = nodes[ 9];
  nids[2] = nodes[ 2];
  nids[3] = nodes[10];
  nids[4] = nodes[26];
  nids[5] = nodes[22];
  nids[6] = nodes[14];
  nids[7] = nodes[23];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[11];
  nids[1] = nodes[20];
  nids[2] = nodes[10];
  nids[3] = nodes[ 3];
  nids[4] = nodes[24];
  nids[5] = nodes[26];
  nids[6] = nodes[23];
  nids[7] = nodes[15];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = nodes[21];
  nids[2] = nodes[26];
  nids[3] = nodes[24];
  nids[4] = nodes[ 4];
  nids[5] = nodes[16];
  nids[6] = nodes[25];
  nids[7] = nodes[19];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[21];
  nids[1] = nodes[13];
  nids[2] = nodes[22];
  nids[3] = nodes[26];
  nids[4] = nodes[16];
  nids[5] = nodes[ 5];
  nids[6] = nodes[17];
  nids[7] = nodes[25];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[26];
  nids[1] = nodes[22];
  nids[2] = nodes[14];
  nids[3] = nodes[23];
  nids[4] = nodes[25];
  nids[5] = nodes[17];
  nids[6] = nodes[ 6];
  nids[7] = nodes[18];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

  nids[0] = nodes[24];
  nids[1] = nodes[26];
  nids[2] = nodes[23];
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = nodes[25];
  nids[6] = nodes[18];
  nids[7] = nodes[ 7];
  subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
}

GEO::CUT::Tet10ElementHandle::Tet10ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  nodes_.reserve( 10 );
  for ( int i=0; i<10; ++i )
  {
    Node * n = mesh.GetNode( nids[i], NULL );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  std::vector<int> subnids( 4 );

  subnids[0] = nids[ 0];
  subnids[1] = nids[ 4];
  subnids[2] = nids[ 6];
  subnids[3] = nids[ 7];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 1];
  subnids[2] = nids[ 5];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 6];
  subnids[1] = nids[ 5];
  subnids[2] = nids[ 2];
  subnids[3] = nids[ 9];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 7];
  subnids[1] = nids[ 8];
  subnids[2] = nids[ 9];
  subnids[3] = nids[ 3];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  /////////////////////////////////////////////////////////////////

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 5];
  subnids[2] = nids[ 6];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 6];
  subnids[1] = nids[ 9];
  subnids[2] = nids[ 7];
  subnids[3] = nids[ 8];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 8];
  subnids[2] = nids[ 7];
  subnids[3] = nids[ 6];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );

  subnids[0] = nids[ 9];
  subnids[1] = nids[ 8];
  subnids[2] = nids[ 5];
  subnids[3] = nids[ 6];
  subelements_.push_back( mesh.GetElement( -1, subnids, *top_data ) );
}

void GEO::CUT::Hex20ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  LINALG::Matrix<3, 20> xyze;

  for ( int i=0; i<20; ++i )
  {
    Node * n = nodes_[i];
    n->Coordinates( &xyze( 0, i ) );
  }

  Position<DRT::Element::hex20> pos( xyze, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::Hex27ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  LINALG::Matrix<3, 27> xyze;

  for ( int i=0; i<27; ++i )
  {
    Node * n = nodes_[i];
    n->Coordinates( &xyze( 0, i ) );
  }

  Position<DRT::Element::hex27> pos( xyze, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::Tet10ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  LINALG::Matrix<3, 10> xyze;

  for ( int i=0; i<10; ++i )
  {
    Node * n = nodes_[i];
    n->Coordinates( &xyze( 0, i ) );
  }

  Position<DRT::Element::tet10> pos( xyze, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}

