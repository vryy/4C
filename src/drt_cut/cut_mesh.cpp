
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

GEO::CUT::Mesh::Mesh( Teuchos::RCP<PointPool> pp, bool cutmesh )
  : pp_( pp ),
    cutmesh_( cutmesh )
{
  if ( pp_ == Teuchos::null )
  {
    pp_ = Teuchos::rcp( new PointPool );
  }
}

void GEO::CUT::Mesh::FillComplete()
{
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
  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Triangle<6> >();
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

GEO::CUT::Node* GEO::CUT::Mesh::GetNode( int nid, const double * xyz )
{
  std::map<int, Teuchos::RCP<Node> >::iterator i = nodes_.find( nid );
  if ( i != nodes_.end() )
  {
    return &*i->second;
  }
  if ( xyz==NULL )
    throw std::runtime_error( "cannot create node without coordinates" );

  Point * p = NewPoint( xyz, NULL, NULL );
  Node * n = new Node( nid, p );
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
  return pp_->NewPoint( x, cut_edge, cut_side );
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
//     std::cout << "NewLine: ";
//     p1->Print();
//     std::cout << "--";
//     p2->Print();
//     std::cout << "\n";
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
  return f;
}

void GEO::CUT::Mesh::Cut( Mesh & mesh )
{
//   int count = 0;
//   std::cout << "\n";
  for ( std::map<std::set<int>, Teuchos::RCP<Side> >::iterator i=sides_.begin();
        i!=sides_.end();
        ++i )
  {
    Side & side = *i->second;
    LinearSide * ls = dynamic_cast<LinearSide*>( &side );
    if ( ls!=NULL )
    {
      mesh.Cut( *ls );

//       std::cout << "." << std::flush;
//       count += 1;
//       if ( count % 64 == 0 )
//         std::cout << "\n";
//       else if ( count % 8 == 0 )
//         std::cout << " " << std::flush;
    }
  }
//   std::cout << "\n";
}

void GEO::CUT::Mesh::Cut( LinearSide & side )
{
  BoundingBox sidebox( side );
  std::set<Element*> elements;
  pp_->CollectElements( sidebox, elements );

  for ( std::set<Element*>::iterator i=elements.begin(); i!=elements.end(); ++i )
  {
    Element * e = *i;
    e->Cut( *this, side );
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
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    e.FindNodePositions();
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
    Point * p = n->point();
    if ( i!=nodes.begin() )
      file << ",";
    file << p->Position();
  }
  file << "};\n";
}

void GEO::CUT::Mesh::GenerateTetgen( CellGenerator * generator )
{
  for ( std::map<int, Teuchos::RCP<Element> >::iterator i=elements_.begin();
        i!=elements_.end();
        ++i )
  {
    Element & e = *i->second;
    if ( e.IsCut() )
    {
      e.GenerateTetgen( *this, generator );
    }
  }
}

bool GEO::CUT::Mesh::WithinBB( const Epetra_SerialDenseMatrix & xyz )
{
  return bb_.Within( xyz );
}
