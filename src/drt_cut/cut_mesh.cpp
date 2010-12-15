
#include <fstream>
#include <iostream>
#include <iterator>

#include "cut_mesh.H"

#include "cut_pointpool.H"
#include "cut_point.H"
#include "cut_point_impl.H"
#include "cut_line.H"
#include "cut_node.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_element.H"
#include "cut_facet.H"
#include "cut_levelsetside.H"
#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_integrationcell.H"

GEO::CUT::Mesh::Mesh( double norm, Teuchos::RCP<PointPool> pp, bool cutmesh )
  : setup_( true ),
    norm_( norm ),
    pp_( pp ),
    cutmesh_( cutmesh )
{
  if ( pp_ == Teuchos::null )
  {
    pp_ = Teuchos::rcp( new PointPool( norm ) );
  }
}

#ifdef QHULL

void GEO::CUT::Mesh::AddTetgen( const tetgenio & out )
{
  std::map<int, int> nodeidmap;

  // need to copy surface markers via side id
  for ( int i=0; i<out.numberoftrifaces; ++i )
  {
    std::vector<int> nids;
    nids.reserve( 3 );
    for ( int j=0; j<3; ++j )
    {
      int nid;
      int pointidx = out.trifacelist[i*3+j] * 3;
      std::map<int, int>::iterator ni = nodeidmap.find( pointidx );
      if ( ni==nodeidmap.end() )
      {
        nid = nodes_.size();
        if ( nodes_.find( nid )!=nodes_.end() )
        {
          throw std::runtime_error( "node ids not consecutive" );
        }
        LINALG::Matrix<3,1> x( &out.pointlist[pointidx] );
        x.Scale( 1./TETGENPOINTSCALE );
        //Node* n =
        GetNode( nid, x.A() );
        nodeidmap[pointidx] = nid;
      }
      else
      {
        nid = ni->second;
      }
      nids.push_back( nid );
    }
    //Side * side =
    CreateTri3( out.trifacemarkerlist[i], nids );
  }

  const int numTetNodes = 4;

  for ( int i=0; i<out.numberoftetrahedra; ++i )
  {
    std::vector<int> nids;
    nids.reserve( numTetNodes );
    for ( int j=0; j<numTetNodes; ++j )
    {
      int nid;
      int pointidx = out.tetrahedronlist[i*out.numberofcorners+j] * 3;
      std::map<int, int>::iterator ni = nodeidmap.find( pointidx );
      if ( ni==nodeidmap.end() )
      {
        nid = nodes_.size();
        if ( nodes_.find( nid )!=nodes_.end() )
        {
          throw std::runtime_error( "node ids not consecutive" );
        }
        LINALG::Matrix<3,1> x( &out.pointlist[pointidx] );
        x.Scale( 1./TETGENPOINTSCALE );
        //Node* n =
        GetNode( nid, x.A() );
        nodeidmap[pointidx] = nid;
      }
      else
      {
        nid = ni->second;
      }
      nids.push_back( nid );
    }
    int eid = elements_.size();
    if ( elements_.find( eid )!=elements_.end() )
    {
      throw std::runtime_error( "element ids not consecutive" );
    }
    //Element * e =
    CreateTet4( eid, nids );
  }
}

void GEO::CUT::Mesh::ExtractTetgen( tetgenio & out )
{
  const int dim = 3;

  out.numberofpoints = nodes_.size();
  out.pointlist = new double[out.numberofpoints * dim];
  out.pointmarkerlist = new int[out.numberofpoints];

  out.numberoftrifaces = sides_.size();
  out.trifacemarkerlist = new int[out.numberoftrifaces];
  out.trifacelist = new int[out.numberoftrifaces * 3];

  out.numberoftetrahedra = elements_.size();
  out.tetrahedronlist = new int[out.numberoftetrahedra * 4];

  std::map<int, int> nodeidmap;

  int count = 0;
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin();
        i!=nodes_.end();
        ++i )
  {
    Node * n = &*i->second;
    n->Coordinates( &out.pointlist[count*dim] );
    for ( int j=0; j<dim; ++j )
    {
      out.pointlist[count*dim+j] *= TETGENPOINTSCALE;
    }
    out.pointmarkerlist[count] = n->point()->Position();
    nodeidmap[n->Id()] = count;
    count += 1;
  }

  count = 0;
  for ( std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side * s = &*i->second;
    const std::vector<Node*> & side_nodes = s->Nodes();

    if ( side_nodes.size()!=3 )
    {
      throw std::runtime_error( "side not a tri3" );
    }

    for ( int j=0; j<3; ++j )
    {
      out.trifacelist[count*3+j] = nodeidmap[side_nodes[j]->Id()];
    }

    int sid = s->Id();
    if ( sid < 0 )
    {
      for ( int j=0; j<3; ++j )
      {
        sid = std::min( sid, out.pointmarkerlist[out.trifacelist[count*3+j]] );
      }
    }
    out.trifacemarkerlist[count] = sid;
    count += 1;
  }

  count = 0;
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element * e = &*i->second;
    const std::vector<Node*> & element_nodes = e->Nodes();

    if ( element_nodes.size()!=4 )
    {
      throw std::runtime_error( "element not a tet4" );
    }

    for ( int j=0; j<4; ++j )
    {
      out.tetrahedronlist[count*4+j] = nodeidmap[element_nodes[j]->Id()];
    }
    count += 1;
  }
}

#endif

void GEO::CUT::Mesh::FillComplete()
{
  setup_ = false;

  // make a local copy of me list of elements since the main list might change
  // inside the fill loop.

  std::map<std::set<int>, Teuchos::RCP<Side> > sides( sides_ );

  for ( std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i=sides.begin();
        i!=sides.end();
        ++i )
  {
    Side & s = *i->second;
    s.FillComplete( *this );
  }

  std::map<int, Teuchos::RCP<Element> > elements( elements_ );

  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements.begin();
        i!=elements.end();
        ++i )
  {
    Element & e = *i->second;
    e.FillComplete( *this );
  }
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateElement( int eid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
  switch ( distype )
  {
  case DRT::Element::hex8:
    return CreateHex8( eid, nids );
  case DRT::Element::hex20:
    return CreateHex20( eid, nids );
  case DRT::Element::hex27:
    return CreateHex27( eid, nids );
  case DRT::Element::tet4:
    return CreateTet4( eid, nids );
  case DRT::Element::tet10:
    return CreateTet10( eid, nids );
  case DRT::Element::pyramid5:
    return CreatePyramid5( eid, nids );
  case DRT::Element::wedge6:
    return CreateWedge6( eid, nids );
  case DRT::Element::wedge15:
  default:
    throw std::runtime_error( "unsupported distype" );
  }
  return NULL;
}

GEO::CUT::Side * GEO::CUT::Mesh::CreateSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype )
{
  switch ( distype )
  {
  case DRT::Element::quad4:
    return CreateQuad4( sid, nids );
  case DRT::Element::quad8:
    return CreateQuad8( sid, nids );
  case DRT::Element::quad9:
    return CreateQuad9( sid, nids );
  case DRT::Element::tri3:
    return CreateTri3( sid, nids );
  case DRT::Element::tri6:
    return CreateTri6( sid, nids );
  default:
    throw std::runtime_error( "unsupported distype" );
  }
  return NULL;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateTet4( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateTet10( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<10> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreatePyramid5( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Pyramid<5> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateWedge6( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Wedge<6> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateHex8( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateHex20( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<20> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}

GEO::CUT::Element * GEO::CUT::Mesh::CreateHex27( int eid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<27> >();

  Element * e = GetElement( eid, nids, *top_data );

  return e;
}


GEO::CUT::Side * GEO::CUT::Mesh::CreateTri3( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Triangle<3> >();
  return GetSide( sid, nids, top_data );
}

GEO::CUT::Side * GEO::CUT::Mesh::CreateTri6( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Triangle<6> >();
  return GetSide( sid, nids, top_data );
}

GEO::CUT::Side * GEO::CUT::Mesh::CreateQuad4( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  return GetSide( sid, nids, top_data );
}

GEO::CUT::Side * GEO::CUT::Mesh::CreateQuad8( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<8> >();
  return GetSide( sid, nids, top_data );
}

GEO::CUT::Side * GEO::CUT::Mesh::CreateQuad9( int sid, const std::vector<int> & nids )
{
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Quadrilateral<9> >();
  return GetSide( sid, nids, top_data );
}

GEO::CUT::Node* GEO::CUT::Mesh::GetNode( int nid, const double * xyz, double lsv )
{
  std::map<int, Teuchos::RCP<Node> >::iterator i = nodes_.find( nid );
  if ( i != nodes_.end() )
  {
    return &*i->second;
  }
  if ( xyz==NULL )
    throw std::runtime_error( "cannot create node without coordinates" );

  Point * p = NewPoint( xyz, NULL, NULL );
  Node * n = new Node( nid, p, lsv );
  nodes_[nid] = Teuchos::rcp( n );
  return n;
}

GEO::CUT::Node* GEO::CUT::Mesh::GetNode( const std::set<int> & nids, const double * xyz )
{
  std::map<std::set<int>, Node*>::iterator i=shadow_nodes_.find( nids );
  if ( i!=shadow_nodes_.end() )
  {
    return &*i->second;
  }
  int nid = - shadow_nodes_.size() - 1;
  if ( nodes_.find( nid )!=nodes_.end() )
  {
    throw std::runtime_error( "shadow node already exists" );
  }
  Node * n = GetNode( nid, xyz );
  shadow_nodes_[nids] = n;
  return n;
}

GEO::CUT::Edge* GEO::CUT::Mesh::GetEdge( Node* begin, Node* end )
{
  std::set<int> nids;
  nids.insert( begin->Id() );
  nids.insert( end->Id() );

  std::vector<Node*> nodes( 2 );
  nodes[0] = begin;
  nodes[1] = end;

  return GetEdge( nids, nodes, *shards::getCellTopologyData< shards::Line<2> >() );
}

GEO::CUT::Edge* GEO::CUT::Mesh::GetEdge( const std::set<int> & nids, const std::vector<Node*> & nodes, const CellTopologyData & edge_topology )
{
  std::map<std::set<int>, Teuchos::RCP<Edge> >::iterator i = edges_.find( nids );
  if ( i != edges_.end() )
  {
    return &*i->second;
  }
  Edge * e = NULL;
  switch ( edge_topology.key )
  {
  case shards::Line<2>::key :
    e = new ConcreteEdge<DRT::Element::line2>( nodes );
    break;
  case shards::Line<3>::key :
    e = new ConcreteEdge<DRT::Element::line3>( nodes );
    break;
  default:
    throw std::runtime_error( "unsupported edge topology" );
  }
  edges_[nids] = Teuchos::rcp( e );
  return e;
}

const std::vector<GEO::CUT::Side*> & GEO::CUT::Mesh::GetSides( int sid )
{
  std::map<int, std::vector<Side*> >::iterator i=cut_sides_.find( sid );
  if ( i!=cut_sides_.end() )
  {
    return i->second;
  }
  throw std::runtime_error( "no side with given id" );
}

GEO::CUT::Side* GEO::CUT::Mesh::GetSide( const std::set<int> & nids )
{
  std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i = sides_.find( nids );
  if ( i != sides_.end() )
  {
    return &*i->second;
  }
  return NULL;
}

GEO::CUT::Side* GEO::CUT::Mesh::GetSide( int sid, const std::vector<int> & nids, const CellTopologyData * top_data )
{
  unsigned nc = top_data->node_count;
  unsigned ec = top_data->edge_count;

  std::vector<Node*> nodes;
  std::vector<Edge*> edges;

  nodes.reserve( nc );
  edges.reserve( ec );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( nids[i], NULL ) );
  }

  for ( unsigned i=0; i<ec; ++i )
  {
    const CellTopologyData_Subcell & edge = top_data->edge[i] ;
    const CellTopologyData & edge_topology = *edge.topology;

    std::vector<Node*> edge_nodes;
    std::set<int> edge_nids;
    edge_nodes.reserve( edge_topology.node_count );
    for ( unsigned j=0; j<edge_topology.node_count; ++j )
    {
      edge_nids.insert( nids[edge.node[j]] );
      edge_nodes.push_back( nodes[edge.node[j]] );
    }
    edges.push_back( GetEdge( edge_nids, edge_nodes, edge_topology ) );
  }

  std::set<int> nidset;
  std::copy( nids.begin(), nids.end(),
             std::inserter( nidset, nidset.begin() ) );
  Side * s = GetSide( sid, nidset, nodes, edges, *top_data );

  return s;
}

GEO::CUT::Side* GEO::CUT::Mesh::GetSide( int sid,
                                         const std::set<int> & nids,
                                         const std::vector<Node*> & nodes,
                                         const std::vector<Edge*> & edges,
                                         const CellTopologyData & side_topology )
{
  std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i = sides_.find( nids );
  if ( i != sides_.end() )
  {
    return &*i->second;
  }
  //Facet * f = new Facet;
  //facets_.push_back( Teuchos::rcp( f ) );
  Side * s = NULL;
  switch ( side_topology.key )
  {
  case shards::Triangle<3>::key :
    s = new ConcreteSide<DRT::Element::tri3>( sid, nodes, edges );
    break;
  case shards::Triangle<6>::key :
    s = new ConcreteSide<DRT::Element::tri6>( sid, nodes, edges );
    break;
  case shards::Quadrilateral<4>::key :
    s = new ConcreteSide<DRT::Element::quad4>( sid, nodes, edges );
    break;
  case shards::Quadrilateral<8>::key :
    s = new ConcreteSide<DRT::Element::quad8>( sid, nodes, edges );
    break;
  case shards::Quadrilateral<9>::key :
    s = new ConcreteSide<DRT::Element::quad9>( sid, nodes, edges );
    break;
  default:
    throw std::runtime_error( "unsupported side topology" );
  }
  sides_[nids] = Teuchos::rcp( s );
  if ( sid > -1 )
  {
    cut_sides_[sid].push_back( s );
  }
  return s;
}

GEO::CUT::Element* GEO::CUT::Mesh::GetElement( int eid )
{
  std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find( eid );
  if ( ie != elements_.end() )
  {
    return &*ie->second;
  }
  return NULL;
}

GEO::CUT::Element* GEO::CUT::Mesh::GetElement( int eid,
                                               const std::vector<int> & nids,
                                               const CellTopologyData & top_data )
{
  static int shards_to_baci[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 26, 20, 25, 24, 22, 21, 23, -1
  };

  std::map<int, Teuchos::RCP<Element> >::iterator ie = elements_.find( eid );
  if ( ie != elements_.end() )
  {
    return &*ie->second;
  }

//   std::cout << "GetElement: ";
//   std::copy( nids.begin(), nids.end(), std::ostream_iterator<int>( std::cout, " " ) );
//   std::cout << "\n";

  unsigned sc = top_data.side_count;
//   unsigned ec = top_data.edge_count;
  unsigned nc = top_data.node_count;

  std::vector<Node*> nodes;
//   std::vector<Edge*> edges;
  std::vector<Side*> sides;

  nodes.reserve( nc );
//   edges.reserve( ec );
  sides.reserve( sc );

  for ( unsigned i=0; i<nc; ++i )
  {
    nodes.push_back( GetNode( nids[i], NULL ) );
  }

//   for ( unsigned i=0; i<ec; ++i )
//   {
//     const CellTopologyData_Subcell & edge = top_data->edge[i] ;
//     const CellTopologyData & edge_topology = *edge.topology;

//     std::vector<Node*> edge_nodes;
//     std::set<int> edge_nids;
//     edge_nodes.reserve( edge_topology.node_count );
//     for ( unsigned j=0; j<edge_topology.node_count; ++j )
//     {
//       edge_nids.add( nids[edge.node[j]] );
//       edge_nodes.push_back( nodes[edge.node[j]] );
//     }
//     edges.push_back( GetEdge( edge_nids, edge_nodes ) );
//   }

  for ( unsigned i=0; i<sc; ++i )
  {
    const CellTopologyData_Subcell & side = top_data.side[i] ;
    const CellTopologyData & side_topology = *side.topology;

    std::vector<Node*> side_nodes;
    std::vector<int> side_nids;
    side_nodes.reserve( side_topology.node_count );
    side_nids .reserve( side_topology.node_count );
    for ( unsigned j=0; j<side_topology.node_count; ++j )
    {
      int nid = shards_to_baci[side.node[j]];
      side_nids .push_back( nids [nid] );
      side_nodes.push_back( nodes[nid] );
    }

    std::vector<Edge*> side_edges;
    side_edges.reserve( side_topology.edge_count );
    for ( unsigned j=0; j<side_topology.edge_count; ++j )
    {
      const CellTopologyData_Subcell & edge = side_topology.edge[j] ;
      const CellTopologyData & edge_topology = *edge.topology;

      std::vector<Node*> edge_nodes;
      std::set<int> edge_nids;
      edge_nodes.reserve( edge_topology.node_count );
      for ( unsigned j=0; j<edge_topology.node_count; ++j )
      {
        edge_nids.insert( side_nids[edge.node[j]] );
        edge_nodes.push_back( side_nodes[edge.node[j]] );
      }

      side_edges.push_back( GetEdge( edge_nids, edge_nodes, edge_topology ) );
    }

    std::set<int> side_nidset;
    std::copy( side_nids.begin(), side_nids.end(),
               std::inserter( side_nidset, side_nidset.begin() ) );
    sides.push_back( GetSide( -1, side_nidset, side_nodes, side_edges, side_topology ) );
  }

  Element * e = NULL;
  switch ( top_data.key )
  {
  case shards::Tetrahedron<4>::key :
    e = new ConcreteElement<DRT::Element::tet4>( eid, sides, nodes );
    break;

  case shards::Hexahedron<8>::key :
    e = new ConcreteElement<DRT::Element::hex8>( eid, sides, nodes );
    break;

  case shards::Pyramid<5>::key :
    e = new ConcreteElement<DRT::Element::pyramid5>( eid, sides, nodes );
    break;

  case shards::Wedge<6>::key :
    e = new ConcreteElement<DRT::Element::wedge6>( eid, sides, nodes );
    break;

  case shards::Hexahedron<20>::key :
    e = new ConcreteElement<DRT::Element::hex20>( eid, sides, nodes );
    break;

  case shards::Hexahedron<27>::key :
    e = new ConcreteElement<DRT::Element::hex27>( eid, sides, nodes );
    break;

  case shards::Tetrahedron<10>::key :
    e = new ConcreteElement<DRT::Element::tet10>( eid, sides, nodes );
    break;

  case shards::Pyramid<13>::key :
  case shards::Pyramid<14>::key :
  case shards::Wedge<15>::key :
  case shards::Wedge<18>::key :
  default:
    throw std::runtime_error( "unsupported element topology" );
  }
  if ( eid > -1 )
  {
    elements_[eid] = Teuchos::rcp( e );
  }
  else
  {
    shadow_elements_.push_back( Teuchos::rcp( e ) );
  }
  return e;
}

GEO::CUT::Point* GEO::CUT::Mesh::NewPoint( const double * x, Edge * cut_edge, Side * cut_side )
{
  bb_.AddPoint( x );
  //Point* p = pp_->NewPoint( x, cut_edge, cut_side, setup_ ? SETUPNODECATCHTOL : MINIMALTOL );
  Point* p = pp_->NewPoint( x, cut_edge, cut_side, MINIMALTOL );
#if 0
  std::cout << "Mesh::NewPoint: ";
  p->Print();
  std::cout << "\n";
#endif
  return p;
}

GEO::CUT::Line* GEO::CUT::Mesh::NewLine( Point* p1, Point* p2, Side * cut_side1, Side * cut_side2, Element * cut_element )
{
  if ( p1==p2 )
    throw std::runtime_error( "no line between same point" );

  // If there is a line between those points already, return it. Otherwise
  // create a new one.
  Line * line = p1->CommonLine( p2 );
  if ( line==NULL )
  {
    line = new Line( p1, p2, cut_side1, cut_side2, cut_element );
    lines_.push_back( Teuchos::rcp( line ) );
#if 0
    std::cout << "Mesh::NewLine: ";
    p1->Print();
    std::cout << "--";
    p2->Print();
    std::cout << "\n";
#endif
  }
  else
  {
    if ( cut_side1 )
      line->AddSide( cut_side1 );
    if ( cut_side2 )
      line->AddSide( cut_side2 );
    if ( cut_element )
      line->AddElement( cut_element );
  }
  return line;
}

GEO::CUT::Facet* GEO::CUT::Mesh::NewFacet( const std::vector<Point*> & points, Side * side, bool cutsurface )
{
  Facet* f = new Facet( *this, points, side, cutsurface );
  facets_.push_back( Teuchos::rcp( f ) );
#if 0
  std::cout << "Mesh::NewFacet: ";
  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    p->Print();
    std::cout << " ";
  }
  std::cout << "\n";
#endif
  return f;
}

GEO::CUT::VolumeCell* GEO::CUT::Mesh::NewVolumeCell( const std::set<Facet*> & facets,
                                                     const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
                                                     LinearElement * element )
{
  VolumeCell * c = new VolumeCell( facets, volume_lines, element );
  cells_.push_back( Teuchos::rcp( c ) );
  return c;
}

GEO::CUT::Tri3BoundaryCell* GEO::CUT::Mesh::NewTri3Cell( const Epetra_SerialDenseMatrix & xyz, Facet * facet )
{
  Tri3BoundaryCell * bc = new Tri3BoundaryCell( xyz, facet );
  boundarycells_.push_back( Teuchos::rcp( bc ) );
  return bc;
}

GEO::CUT::Quad4BoundaryCell* GEO::CUT::Mesh::NewQuad4Cell( const Epetra_SerialDenseMatrix & xyz, Facet * facet )
{
  Quad4BoundaryCell * bc = new Quad4BoundaryCell( xyz, facet );
  boundarycells_.push_back( Teuchos::rcp( bc ) );
  return bc;
}

GEO::CUT::Hex8IntegrationCell* GEO::CUT::Mesh::NewHex8Cell( Point::PointPosition position,
                                                            const std::vector<Point*> & points,
                                                            VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }
  Hex8IntegrationCell * c = new Hex8IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );
  return c;
}

GEO::CUT::Tet4IntegrationCell* GEO::CUT::Mesh::NewTet4Cell( Point::PointPosition position,
                                                            const std::vector<Point*> & points,
                                                            VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }
  Tet4IntegrationCell * c = new Tet4IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );
  return c;
}

GEO::CUT::Tet4IntegrationCell* GEO::CUT::Mesh::NewTet4Cell( Point::PointPosition position,
                                                            const Epetra_SerialDenseMatrix & xyz,
                                                            VolumeCell * cell )
{
  std::vector<Point*> points;   // empty list of points
  Tet4IntegrationCell * c = new Tet4IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );
  return c;
}

GEO::CUT::Wedge6IntegrationCell* GEO::CUT::Mesh::NewWedge6Cell( Point::PointPosition position,
                                                                const std::vector<Point*> & points,
                                                                VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }
  Wedge6IntegrationCell * c = new Wedge6IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );
  return c;
}

GEO::CUT::Pyramid5IntegrationCell* GEO::CUT::Mesh::NewPyramid5Cell( Point::PointPosition position,
                                                                    const std::vector<Point*> & points,
                                                                    VolumeCell * cell )
{
  Epetra_SerialDenseMatrix xyz( 3, points.size() );
  for ( unsigned i=0; i<points.size(); ++i )
  {
    points[i]->Coordinates( &xyz( 0, i ) );
  }
  Pyramid5IntegrationCell * c = new Pyramid5IntegrationCell( position, xyz, points, cell );
  integrationcells_.push_back( Teuchos::rcp( c ) );
  return c;
}


void GEO::CUT::Mesh::SelfCut()
{
  std::set<Facet*> facets;
  for ( std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side & side = *i->second;
    LinearSide * ls1 = dynamic_cast<LinearSide*>( &side );
    if ( ls1!=NULL )
    {
      BoundingBox sidebox( side );
      std::set<Side*> sides;
      pp_->CollectSides( sidebox, sides );
      sides.erase( ls1 );
      for ( std::set<Side*>::iterator i=sides.begin(); i!=sides.end(); )
      {
        Side * s = *i;
        if ( ls1->HaveCommonEdge( *s ) )
        {
          sides.erase( i++ );
        }
        else
        {
          ++i;
        }
      }
      for ( std::set<Side*>::iterator i=sides.begin(); i!=sides.end(); ++i )
      {
        Side * s = *i;
        LinearSide * ls2 = dynamic_cast<LinearSide*>( s );
        if ( ls2!=NULL )
        {
          ls1->FindCutPoints( *this, NULL, *ls2 );
          ls2->FindCutPoints( *this, NULL, *ls1 );
        }
      }
      bool cut = false;
      for ( std::set<Side*>::iterator i=sides.begin(); i!=sides.end(); ++i )
      {
        Side * s = *i;
        LinearSide * ls2 = dynamic_cast<LinearSide*>( s );
        if ( ls2!=NULL )
        {
          bool normal_cut  = ls1->FindCutLines( *this, NULL, *ls2 );
          bool reverse_cut = ls2->FindCutLines( *this, NULL, *ls1 );
          if ( normal_cut or reverse_cut )
            cut = true;
        }
      }
      if ( cut )
      {
        for ( std::set<Side*>::iterator i=sides.begin(); i!=sides.end(); ++i )
        {
          Side * s = *i;
          LinearSide * ls2 = dynamic_cast<LinearSide*>( s );
          if ( ls2!=NULL )
          {
            SideSideCutFilter filter( ls1, ls2 );
            ls1->MakeOwnedSideFacets( *this, filter, facets );
          }
        }
        ls1->MakeSideCutFacets( *this, NULL, facets );
      }
    }
  }
  for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    f->CreateLinearElements( *this );
  }
}

void GEO::CUT::Mesh::Cut( Mesh & mesh, std::set<Element*> & elements_done )
{
  std::set<Element*> my_elements_done;

  for ( std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side & side = *i->second;
    LinearSide * ls = dynamic_cast<LinearSide*>( &side );
    if ( ls!=NULL )
    {
      mesh.Cut( *ls, elements_done, my_elements_done );
    }
  }

  std::copy( my_elements_done.begin(),
             my_elements_done.end(),
             std::inserter( elements_done, elements_done.begin() ) );
}

void GEO::CUT::Mesh::Cut( LinearSide & side, const std::set<Element*> & done, std::set<Element*> & elements_done )
{
  BoundingBox sidebox( side );
  std::set<Element*> elements;
  pp_->CollectElements( sidebox, elements );

  for ( std::set<Element*>::iterator i=elements.begin(); i!=elements.end(); ++i )
  {
    Element * e = *i;
    if ( done.count( e )==0 )
    {
      if ( e->Cut( *this, side ) )
      {
        elements_done.insert( e );
      }
    }
  }
}

void GEO::CUT::Mesh::Cut( LevelSetSide & side )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    e.Cut( *this, side );
  }
}

void GEO::CUT::Mesh::MakeFacets()
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    e.MakeFacets( *this );
  }
}

void GEO::CUT::Mesh::MakeVolumeCells()
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    e.MakeVolumeCells( *this );
  }
}

void GEO::CUT::Mesh::FindNodePositions()
{
  // On multiple cuts former outside positions can become inside
  // positions. Thus reset all outside positions.
  pp_->ResetOutsidePoints();

  // get nodal positions from elements
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    e.FindNodePositions();
  }

  // If there are any undecided facets left, those should be outside. This can
  // only happen in test cases where there is no edge cut to determine a
  // genuine facet position.
  for ( std::list<Teuchos::RCP<Facet> >::iterator i=facets_.begin();
        i!=facets_.end();
        ++i )
  {
    Facet * f = &**i;
    if ( f->Position()==Point::undecided )
    {
      f->Position( Point::outside );
    }
  }
}

void GEO::CUT::Mesh::FindLSNodePositions()
{
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin();
        i!=nodes_.end();
        ++i )
  {
    Node * n = &*i->second;
    Point * p = n->point();
    if ( p->Position()==Point::undecided )
    {
      double lsv = n->LSV();
      if ( lsv > 0 )
      {
        p->Position( Point::outside );
      }
      else if ( lsv < 0 )
      {
        p->Position( Point::inside );
      }
      else
      {
        throw std::runtime_error( "undecided nodal point on levelset surface" );
      }
    }
  }
}

void GEO::CUT::Mesh::FindNodalDOFSets()
{
  for ( std::map<int, Teuchos::RCP<Node> >::iterator i=nodes_.begin();
        i!=nodes_.end();
        ++i )
  {
    Node * n = &*i->second;
    n->FindDOFSets();
  }

  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * cell = &**i;
    cell->ConnectNodalDOFSets();
  }
}

void GEO::CUT::Mesh::CreateIntegrationCells()
{
  for ( std::list<Teuchos::RCP<VolumeCell> >::iterator i=cells_.begin();
        i!=cells_.end();
        ++i )
  {
    VolumeCell * cell = &**i;
    cell->CreateIntegrationCells( *this );
  }
}

void GEO::CUT::Mesh::Status()
{
  std::cout << "#points:    " << pp_->size() << "\n"
            << "#lines:     " << lines_.size() << "\n"
            << "#facets:    " << facets_.size() << "\n"
            << "#nodes:     " << nodes_.size() << "\n"
            << "#edges:     " << edges_.size() << "\n"
            << "#sides:     " << sides_.size() << "\n"
            << "#elements:  " << elements_.size() << "\n"
            << "#snodes:    " << shadow_nodes_.size() << "\n"
            << "#selements: " << shadow_elements_.size() << "\n"
            << "\n";

//   std::cout << "GetElement: ";
//   std::copy( nids.begin(), nids.end(), std::ostream_iterator<int>( std::cout, " " ) );
//   std::cout << "\n";

#if 0
  std::cout << "edges: {\n";
  for ( std::map<std::set<int>, Teuchos::RCP<Edge> >::iterator i=edges_.begin();
        i!=edges_.end();
        ++i )
  {
    const std::set<int> & s = i->first;
    std::cout << "  ";
    std::copy( s.begin(), s.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
  }
  std::cout << "}\n";
#endif

  std::string name;

  if ( cutmesh_ )
  {
    name = "cutmesh";
  }
  else
  {
    name = "mesh";
  }

  std::string linefile = name + "_line.plot";
  std::ofstream lf( linefile.c_str() );

  for ( std::list<Teuchos::RCP<Line > >::iterator i=lines_.begin();
        i!=lines_.end();
        ++i )
  {
    ( *i )->Plot( lf );
  }

  std::string edgefile = name + "_edge.plot";
  std::ofstream ef( edgefile.c_str() );

  for ( std::map<std::set<int>, Teuchos::RCP<Edge> >::iterator i=edges_.begin();
        i!=edges_.end();
        ++i )
  {
    i->second->Plot( ef );
  }
}

void GEO::CUT::Mesh::PrintFacets()
{
//   for ( std::list<Teuchos::RCP<Point> >::iterator i=points_.begin();
//         i!=points_.end();
//         ++i )
//   {
//     Point & p = **i;
//     p.Print();
//     std::cout << "\n";
//   }

  for ( std::list<Teuchos::RCP<Facet> >::iterator i=facets_.begin();
        i!=facets_.end();
        ++i )
  {
    Facet & f = **i;
    f.Print();
  }
}

void GEO::CUT::Mesh::DumpGmsh( std::string name )
{
  std::ofstream file( name.c_str() );
  file << "View \"" << name << "\" {\n";
  if ( elements_.size() > 0 )
  {
    for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
          i!=elements_.end();
          ++i )
    {
      Element & e = *i->second;
      //e.DumpGmsh( file );
      LinearElement * le = dynamic_cast<LinearElement*>( &e );
      if ( le!=NULL )
      {
        const std::vector<Node*> & nodes = e.Nodes();
        char elementtype;
        switch ( nodes.size() )
        {
        case 8:
          elementtype = 'H';
          break;
        case 4:
          elementtype = 'S';
          break;
        case 6:
          elementtype = 'I';
          break;
        default:
          throw std::runtime_error( "unknown element type" );
        }
        DumpGmsh( file, nodes, elementtype );
      }
    }
    for ( std::list<Teuchos::RCP<Element> >::iterator i=shadow_elements_.begin();
          i!=shadow_elements_.end();
          ++i )
    {
      Element & e = **i;
      const std::vector<Node*> & nodes = e.Nodes();
      char elementtype;
      switch ( nodes.size() )
      {
      case 8:
        elementtype = 'H';
        break;
      case 4:
        elementtype = 'S';
        break;
      case 6:
        elementtype = 'I';
        break;
      default:
        throw std::runtime_error( "unknown element type" );
      }
      DumpGmsh( file, nodes, elementtype );
    }
  }
  else
  {
    for ( std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i=sides_.begin();
          i!=sides_.end();
          ++i )
    {
      Side & s = *i->second;
      LinearSide * ls = dynamic_cast<LinearSide*>( &s );
      if ( ls!=NULL )
      {
        const std::vector<Node*> & nodes = ls->Nodes();
        char elementtype;
        switch ( nodes.size() )
        {
        case 3:
          elementtype = 'T';
          break;
        case 4:
          elementtype = 'Q';
          break;
        default:
          throw std::runtime_error( "unknown element type" );
        }
        DumpGmsh( file, nodes, elementtype );
      }
    }
  }
  file << "};\n";
}

void GEO::CUT::Mesh::DumpGmsh( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype )
{
  file << "S" << elementtype
       << "(";
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    double x[3];
    n->Coordinates( x );
    if ( i!=nodes.begin() )
      file << ",";
    file << x[0] << "," << x[1] << "," << x[2];
  }
  file << "){";
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    //Point * p = n->point();
    if ( i!=nodes.begin() )
      file << ",";
    //file << p->Position();
    file << n->DofSets().size();
  }
  file << "};\n";
}

void GEO::CUT::Mesh::DumpGmshIntegrationcells( std::string name )
{
  std::ofstream file( name.c_str() );
  file << "View \"" << name << "\" {\n";

  for ( std::list<Teuchos::RCP<IntegrationCell> >::iterator i=integrationcells_.begin();
        i!=integrationcells_.end();
        ++i )
  {
    IntegrationCell * ic = &**i;
    ic->DumpGmsh( file );
  }

  file << "};\n";
}

bool GEO::CUT::Mesh::WithinBB( const Epetra_SerialDenseMatrix & xyz )
{
  return bb_.Within( norm_, xyz );
}

bool GEO::CUT::Mesh::WithinBB( Element & element )
{
  return bb_.Within( norm_, element );
}
