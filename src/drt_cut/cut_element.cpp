
//#include "../drt_geometry/intersection_templates.H"

#include "cut_tetgen.H"
#include "cut_intersection.H"
#include "cut_position.H"
#include "cut_facet.H"

#include <string>
#include <stack>

#include "cut_element.H"

void GEO::CUT::Element::FillComplete( Mesh & mesh )
{
  for ( std::vector<Side*>::iterator i=sides_.begin(); i!=sides_.end(); ++i )
  {
    Side & myside = **i;
    myside.FillComplete( mesh );
  }
}

bool GEO::CUT::LinearElement::Cut( Mesh & mesh, LinearSide & side )
{
  bool cut = false;

  // find nodal points inside the element
  const std::vector<Node*> side_nodes = side.Nodes();
  for ( std::vector<Node*>::const_iterator i=side_nodes.begin(); i!=side_nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();

    if ( not p->IsCut( this ) )
    {
      if ( PointInside( p ) )
      {
        p->AddElement( this );
        cut = true;
      }
    }
    else
    {
      cut = true;
    }
  }

  const std::vector<Side*> & sides = Sides();

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    LinearSide * s = dynamic_cast<LinearSide*>( *i );
    FindCutPoints( mesh, *s, side );
  }

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    LinearSide * s = dynamic_cast<LinearSide*>( *i );
    if ( FindCutLines( mesh, *s, side ) )
    {
      cut = true;
    }
  }

  // find lines inside the element
  const std::vector<Edge*> & side_edges = side.Edges();
  for ( std::vector<Edge*>::const_iterator i=side_edges.begin(); i!=side_edges.end(); ++i )
  {
    Edge * e = *i;
    std::vector<Point*> line;
    e->CutPointsInside( this, line );
    std::vector<Point*>::iterator i = line.begin();
    if ( i!=line.end() )
    {
      Point * bp = *i;
      for ( ++i; i!=line.end(); ++i )
      {
        Point * ep = *i;
        mesh.NewLine( bp, ep, &side, NULL, this );
        bp = ep;
      }
    }
  }

  if ( cut )
  {
    side.CreateLineSegment( mesh, this );
    cut_faces_.insert( &side );
    return true;
  }
  else
  {
    return false;
  }
}

bool GEO::CUT::LinearElement::FindCutPoints( Mesh & mesh, LinearSide & side, LinearSide & other )
{
  bool cut = side.FindCutPoints( mesh, this, other );
  bool reverse_cut = other.FindCutPoints( mesh, this, side );
  return cut or reverse_cut;
}

bool GEO::CUT::LinearElement::FindCutLines( Mesh & mesh, LinearSide & side, LinearSide & other )
{
  bool cut = side.FindCutLines( mesh, this, other );
  bool reverse_cut = other.FindCutLines( mesh, this, side );
  return cut or reverse_cut;
}

void GEO::CUT::LinearElement::MakeFacets( Mesh & mesh )
{
  if ( facets_.size()==0 )
  {
    const std::vector<Side*> & sides = Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side & side = **i;
      SideElementCutFilter filter( &side, this );
      side.MakeOwnedSideFacets( mesh, filter, facets_ );
    }
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side & side = **i;
      side.MakeSideCutFacets( mesh, this, facets_ );
    }
    for ( std::set<Side*>::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side & cut_side = **i;
      cut_side.MakeInternalFacets( mesh, this, facets_ );
    }
  }
}

void GEO::CUT::LinearElement::FindNodePositions()
{
  LINALG::Matrix<3,1> xyz;
  LINALG::Matrix<3,1> rst;

  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();
    Point::PointPosition pos = p->Position();
    if ( pos==Point::undecided )
    {
      bool done = false;
      const std::set<Facet*> & facets = p->Facets();
      for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        for ( std::set<Side*>::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
        {
          Side * s = *i;

          // Only take a side that belongs to one of this points facets and
          // shares a cut edge with this point. If there are multiple cut
          // sides within the element (facets), only the close one will always
          // give the right direction.
          if ( f->IsCutSide( s ) and p->CommonCutEdge( s )!=NULL )
          {
            if ( p->IsCut( s ) )
            {
              p->Position( Point::oncutsurface );
            }
            else
            {
              p->Coordinates( xyz.A() );
              s->LocalCoordinates( xyz, rst );
              double d = rst( 2, 0 );
              if ( fabs( d ) > MINIMALTOL )
              {
                if ( d > 0 )
                {
                  p->Position( Point::outside );
                }
                else
                {
                  p->Position( Point::inside );
                }
              }
              else
              {
                // within the cut plane but not cut by the side
                break;
              }
            }
            done = true;
            break;
          }
        }
        if ( done )
          break;
      }
      if ( p->Position()==Point::undecided )
      {
        // Still undecided! No facets with cut side attached! Will be set in a
        // minute.
      }
    }
    else if ( pos==Point::outside or pos==Point::inside )
    {
      // The nodal position is already known. Set it to my facets. If the
      // facets are already set, this will not have much effect anyway. But on
      // multiple cuts we avoid unset facets this way.
      const std::set<Facet*> & facets = p->Facets();
      for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        f->Position( pos );
      }
    }
  }
}

bool GEO::CUT::LinearElement::IsCut()
{
  if ( cut_faces_.size()>0 )
  {
    return true;
  }
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side & side = **i;
    if ( side.IsCut() )
    {
      return true;
    }
  }
  return false;
}

void GEO::CUT::LinearElement::GenerateTetgen( Mesh & mesh, CellGenerator * generator )
{
  GenerateTetgen( mesh, this, generator );
}

void GEO::CUT::LinearElement::GenerateTetgen( Mesh & mesh, Element * parent, CellGenerator * generator )
{
#ifdef QHULL
  const int dim = 3;
  tetgenio in;

  std::set<Point*, PointPidLess> points;
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet & f = **i;
    f.GetPoints( mesh, points );
  }

  // allocate pointlist
  in.numberofpoints = points.size();
  in.pointlist = new double[in.numberofpoints * dim];

  int pos = 0;
  for ( std::set<Point*, PointPidLess>::iterator i=points.begin();
        i!=points.end();
        ++i )
  {
    Point & p = **i;
    p.Coordinates( & in.pointlist[pos*dim] );
    for ( int j=0; j<dim; ++j )
    {
      in.pointlist[pos*dim+j] *= TETGENPOINTSCALE;
    }
    pos += 1;
  }

  in.pointmarkerlist = new int[in.numberofpoints];
  std::fill( in.pointmarkerlist, in.pointmarkerlist+in.numberofpoints, 0 );

  in.numberoffacets = 0;
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet & facet = **i;
    in.numberoffacets += facet.NumTetgenFacets( mesh );
  }
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  std::fill( in.facetmarkerlist, in.facetmarkerlist+in.numberoffacets, 2 );

  std::vector<Point*> pointlist;
  pointlist.reserve( points.size() );
  std::copy( points.begin(), points.end(), std::back_inserter( pointlist ) );

  pos = 0;
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet & facet = **i;
    int numtets = facet.NumTetgenFacets( mesh );
    for ( int j=0; j<numtets; ++j )
    {
      tetgenio::facet & f = in.facetlist[pos];

      facet.GenerateTetgen( mesh, this, f, j, in.facetmarkerlist[pos], pointlist );

      pos += 1;
    }
  }

  if ( generator==NULL )
  {
    // debug
    const char * name = "tet";
    in.save_nodes( const_cast<char*>( name ) );
    in.save_poly( const_cast<char*>( name ) );
  }

  char switches[] = "pQ";    // pQ o2 Y R nn
  tetgenio out;

  try
  {
    tetrahedralize( switches, &in, &out );
  }
  catch ( int & err )
  {
    if ( generator!=NULL )
    {
      const char * name = "tet";
      in.save_nodes( const_cast<char*>( name ) );
      in.save_poly( const_cast<char*>( name ) );
    }

    throw;
  }

  if ( generator!=NULL )
  {
    generator->Generate( parent, out );
  }
  else
  {
    // debug
    // Output mesh
    const char * name = "tetout";
    out.save_nodes( const_cast<char*>( name ) );
    out.save_elements( const_cast<char*>( name ) );
    out.save_faces( const_cast<char*>( name ) );
  }
#endif
}

bool GEO::CUT::LinearElement::OnSide( const std::vector<Point*> & facet_points )
{
  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Point*>::const_iterator i=facet_points.begin();
        i!=facet_points.end();
        ++i )
  {
    Point * p = *i;
    if ( not p->NodalPoint( nodes ) )
    {
      return false;
    }
  }

  std::set<Point*, PointPidLess> points;
  std::copy( facet_points.begin(), facet_points.end(),
             std::inserter( points, points.begin() ) );

  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side & side = **i;
    if ( side.OnSide( points ) )
    {
      return true;
    }
  }

  return false;
}

void GEO::CUT::LinearElement::GetCutPoints( std::set<Point*> & cut_points )
{
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side * side = *i;
    LinearSide * ls = dynamic_cast<LinearSide*>( side );
    if ( ls==NULL )
    {
      throw std::runtime_error( "linear element needs linear sides" );
    }

    for ( std::set<Side*>::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side * other = *i;
      LinearSide * ls_other = dynamic_cast<LinearSide*>( other );
      if ( ls_other==NULL )
      {
        throw std::runtime_error( "linear element needs linear side cuts" );
      }
      ls->GetCutPoints( this, *ls_other, cut_points );
    }
  }
}

void GEO::CUT::LinearElement::MakeVolumeCells( Mesh & mesh )
{
  std::map<std::pair<Point*, Point*>, std::set<Facet*> > lines;
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetLines( lines );
  }

  // Alle Facets einsammeln, die sich zu zweit eine Linie teilen. Jede Seite
  // nur ein Facet. Das sollte im Element eindeutig sein. Damit sollten die
  // Volumina erstellt werden können. Löcher sind noch ein extra Thema.

  std::set<Facet*> facets_done;

  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( facets_done.count( f )==0 and OwnedSide( f->ParentSide() ) )
    {
      std::stack<Facet*> new_facets;
      std::set<Facet*> collected_facets;
      std::set<Side*> sides_done;

      new_facets.push( f );

      while ( not new_facets.empty() )
      {
        Facet * f = new_facets.top();
        new_facets.pop();

        collected_facets.insert( f );
        sides_done.insert( f->ParentSide() );
        std::map<std::pair<Point*, Point*>, std::set<Facet*> > facet_lines;
        f->GetLines( facet_lines );

        for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=facet_lines.begin();
              i!=facet_lines.end();
              ++i )
        {
          const std::pair<Point*, Point*> & line = i->first;
          std::set<Facet*> & facets = lines[line];

          Facet * found_facet = NULL;
          for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
          {
            Facet * f = *i;
            if ( collected_facets.count( f )==0 and
                 ( not OwnedSide( f->ParentSide() ) or
                   sides_done.count( f->ParentSide() )==0 ) )
            {
              if ( found_facet==NULL )
              {
                found_facet = f;
              }
              else
              {
                // undecided. Ignore all matches.
                found_facet = NULL;
                break;
              }
            }
          }
          if ( found_facet!=NULL )
          {
            new_facets.push( found_facet );
          }
        }
      }

      // test for open lines in collected_facets

      std::map<std::pair<Point*, Point*>, std::set<Facet*> > volume_lines;
      for ( std::set<Facet*>::iterator i=collected_facets.begin();
            i!=collected_facets.end();
            ++i )
      {
        Facet * f = *i;
        f->GetLines( volume_lines );
      }

      for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=volume_lines.begin();
            i!=volume_lines.end();
            ++i )
      {
        std::set<Facet*> & facets = i->second;
        if ( facets.size()!=2 )
        {
          throw std::runtime_error( "not properly closed line in volume cell" );
        }
      }

      // Create new cell and remember done stuff!

      std::copy( collected_facets.begin(),
                 collected_facets.end(),
                 std::inserter( facets_done, facets_done.begin() ) );

      cells_.insert( mesh.NewVolumeCell( collected_facets, volume_lines, this ) );
    }
  }
}

bool GEO::CUT::QuadraticElement::Cut( Mesh & mesh, LinearSide & side )
{
  bool cut = false;
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    if ( e->Cut( mesh, side ) )
    {
      cut = true;
    }
  }
  return cut;
}

void GEO::CUT::QuadraticElement::MakeFacets( Mesh & mesh )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->MakeFacets( mesh );
  }
}

void GEO::CUT::QuadraticElement::FindNodePositions()
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->FindNodePositions();
  }
}

bool GEO::CUT::QuadraticElement::IsCut()
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    if ( e->IsCut() )
      return true;
  }
  return false;
}

void GEO::CUT::QuadraticElement::GenerateTetgen( Mesh & mesh, CellGenerator * generator )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
    //for ( std::vector<Element*>::reverse_iterator i=subelements_.rbegin(); i!=subelements_.rend(); ++i )
  {
    LinearElement * e = dynamic_cast<LinearElement*>( *i );
    e->GenerateTetgen( mesh, this, generator );
  }
}

bool GEO::CUT::QuadraticElement::OnSide( const std::vector<Point*> & facet_points )
{
  throw std::runtime_error( "not supposed to end up here" );
}


void GEO::CUT::ConcreteElement<DRT::Element::tet10>::FillComplete( Mesh & mesh )
{
  if ( subelements_.size()==0 )
  {
    subelements_.reserve( 8 );

    Element::FillComplete( mesh );

    const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();
    const std::vector<Node*> & nodes = Nodes();

    std::vector<int> nids( 4 );

    nids[0] = nodes[ 0]->Id();
    nids[1] = nodes[ 4]->Id();
    nids[2] = nodes[ 6]->Id();
    nids[3] = nodes[ 7]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 4]->Id();
    nids[1] = nodes[ 1]->Id();
    nids[2] = nodes[ 5]->Id();
    nids[3] = nodes[ 8]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 6]->Id();
    nids[1] = nodes[ 5]->Id();
    nids[2] = nodes[ 2]->Id();
    nids[3] = nodes[ 9]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 7]->Id();
    nids[1] = nodes[ 8]->Id();
    nids[2] = nodes[ 9]->Id();
    nids[3] = nodes[ 3]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    /////////////////////////////////////////////////////////////////

    nids[0] = nodes[ 4]->Id();
    nids[1] = nodes[ 5]->Id();
    nids[2] = nodes[ 6]->Id();
    nids[3] = nodes[ 8]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 6]->Id();
    nids[1] = nodes[ 9]->Id();
    nids[2] = nodes[ 7]->Id();
    nids[3] = nodes[ 8]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 4]->Id();
    nids[1] = nodes[ 8]->Id();
    nids[2] = nodes[ 7]->Id();
    nids[3] = nodes[ 6]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 9]->Id();
    nids[1] = nodes[ 8]->Id();
    nids[2] = nodes[ 5]->Id();
    nids[3] = nodes[ 6]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
  }
}

void GEO::CUT::ConcreteElement<DRT::Element::hex20>::FillComplete( Mesh & mesh )
{
  if ( subelements_.size()==0 )
  {
    subelements_.reserve( 8 );

    Element::FillComplete( mesh );

    const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();
    const std::vector<Node*> & nodes = Nodes();

    // create middle nodes

    std::set<int> node_nids;

    node_nids.clear();
    node_nids.insert( nodes[ 0]->Id() );
    node_nids.insert( nodes[ 8]->Id() );
    node_nids.insert( nodes[ 1]->Id() );
    node_nids.insert( nodes[ 9]->Id() );
    node_nids.insert( nodes[ 2]->Id() );
    node_nids.insert( nodes[10]->Id() );
    node_nids.insert( nodes[ 3]->Id() );
    node_nids.insert( nodes[11]->Id() );
    Node* node20 = mesh.GetNode( node_nids, NULL );
    int node20_id = node20->Id();

    node_nids.clear();
    node_nids.insert( nodes[ 0]->Id() );
    node_nids.insert( nodes[ 8]->Id() );
    node_nids.insert( nodes[ 1]->Id() );
    node_nids.insert( nodes[13]->Id() );
    node_nids.insert( nodes[ 5]->Id() );
    node_nids.insert( nodes[16]->Id() );
    node_nids.insert( nodes[ 4]->Id() );
    node_nids.insert( nodes[12]->Id() );
    Node* node21 = mesh.GetNode( node_nids, NULL );
    int node21_id = node21->Id();

    node_nids.clear();
    node_nids.insert( nodes[ 1]->Id() );
    node_nids.insert( nodes[ 9]->Id() );
    node_nids.insert( nodes[ 2]->Id() );
    node_nids.insert( nodes[14]->Id() );
    node_nids.insert( nodes[ 6]->Id() );
    node_nids.insert( nodes[17]->Id() );
    node_nids.insert( nodes[ 5]->Id() );
    node_nids.insert( nodes[13]->Id() );
    Node* node22 = mesh.GetNode( node_nids, NULL );
    int node22_id = node22->Id();

    node_nids.clear();
    node_nids.insert( nodes[ 2]->Id() );
    node_nids.insert( nodes[10]->Id() );
    node_nids.insert( nodes[ 3]->Id() );
    node_nids.insert( nodes[15]->Id() );
    node_nids.insert( nodes[ 7]->Id() );
    node_nids.insert( nodes[18]->Id() );
    node_nids.insert( nodes[ 6]->Id() );
    node_nids.insert( nodes[14]->Id() );
    Node* node23 = mesh.GetNode( node_nids, NULL );
    int node23_id = node23->Id();

    node_nids.clear();
    node_nids.insert( nodes[ 3]->Id() );
    node_nids.insert( nodes[11]->Id() );
    node_nids.insert( nodes[ 0]->Id() );
    node_nids.insert( nodes[12]->Id() );
    node_nids.insert( nodes[ 4]->Id() );
    node_nids.insert( nodes[19]->Id() );
    node_nids.insert( nodes[ 7]->Id() );
    node_nids.insert( nodes[15]->Id() );
    Node* node24 = mesh.GetNode( node_nids, NULL );
    int node24_id = node24->Id();

    node_nids.clear();
    node_nids.insert( nodes[ 4]->Id() );
    node_nids.insert( nodes[16]->Id() );
    node_nids.insert( nodes[ 5]->Id() );
    node_nids.insert( nodes[17]->Id() );
    node_nids.insert( nodes[ 6]->Id() );
    node_nids.insert( nodes[18]->Id() );
    node_nids.insert( nodes[ 7]->Id() );
    node_nids.insert( nodes[19]->Id() );
    Node* node25 = mesh.GetNode( node_nids, NULL );
    int node25_id = node25->Id();

    LINALG::Matrix<3,20> xyze;
    Coordinates( xyze );

    LINALG::Matrix<20,1> funct;
    DRT::UTILS::shape_function_3D( funct, 0, 0, 0, DRT::Element::hex20 );

    LINALG::Matrix<3,1> xyz;
    xyz.Multiply( xyze, funct );
    node_nids.clear();
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      node_nids.insert( n->Id() );
    }
    Node* node26 = mesh.GetNode( node_nids, xyz.A() );
    int node26_id = node26->Id();


    std::vector<int> nids( 8 );

    nids[0] = nodes[ 0]->Id();
    nids[1] = nodes[ 8]->Id();
    nids[2] = node20_id;
    nids[3] = nodes[11]->Id();
    nids[4] = nodes[12]->Id();
    nids[5] = node21_id;
    nids[6] = node26_id;
    nids[7] = node24_id;
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 8]->Id();
    nids[1] = nodes[ 1]->Id();
    nids[2] = nodes[ 9]->Id();
    nids[3] = node20_id;
    nids[4] = node21_id;
    nids[5] = nodes[13]->Id();
    nids[6] = node22_id;
    nids[7] = node26_id;
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = node20_id;
    nids[1] = nodes[ 9]->Id();
    nids[2] = nodes[ 2]->Id();
    nids[3] = nodes[10]->Id();
    nids[4] = node26_id;
    nids[5] = node22_id;
    nids[6] = nodes[14]->Id();
    nids[7] = node23_id;
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[11]->Id();
    nids[1] = node20_id;
    nids[2] = nodes[10]->Id();
    nids[3] = nodes[ 3]->Id();
    nids[4] = node24_id;
    nids[5] = node26_id;
    nids[6] = node23_id;
    nids[7] = nodes[15]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    /////////////////////////////////////////////////////////////////

    nids[0] = nodes[12]->Id();
    nids[1] = node21_id;
    nids[2] = node26_id;
    nids[3] = node24_id;
    nids[4] = nodes[ 4]->Id();
    nids[5] = nodes[16]->Id();
    nids[6] = node25_id;
    nids[7] = nodes[19]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = node21_id;
    nids[1] = nodes[13]->Id();
    nids[2] = node22_id;
    nids[3] = node26_id;
    nids[4] = nodes[16]->Id();
    nids[5] = nodes[ 5]->Id();
    nids[6] = nodes[17]->Id();
    nids[7] = node25_id;
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = node26_id;
    nids[1] = node22_id;
    nids[2] = nodes[14]->Id();
    nids[3] = node23_id;
    nids[4] = node25_id;
    nids[5] = nodes[17]->Id();
    nids[6] = nodes[ 6]->Id();
    nids[7] = nodes[18]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = node24_id;
    nids[1] = node26_id;
    nids[2] = node23_id;
    nids[3] = nodes[15]->Id();
    nids[4] = nodes[19]->Id();
    nids[5] = node25_id;
    nids[6] = nodes[18]->Id();
    nids[7] = nodes[ 7]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
  }
}

void GEO::CUT::ConcreteElement<DRT::Element::hex27>::FillComplete( Mesh & mesh )
{
  if ( subelements_.size()==0 )
  {
    subelements_.reserve( 8 );

    Element::FillComplete( mesh );

    const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();
    const std::vector<Node*> & nodes = Nodes();

    std::vector<int> nids( 8 );

    nids[0] = nodes[ 0]->Id();
    nids[1] = nodes[ 8]->Id();
    nids[2] = nodes[20]->Id();
    nids[3] = nodes[11]->Id();
    nids[4] = nodes[12]->Id();
    nids[5] = nodes[21]->Id();
    nids[6] = nodes[26]->Id();
    nids[7] = nodes[24]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[ 8]->Id();
    nids[1] = nodes[ 1]->Id();
    nids[2] = nodes[ 9]->Id();
    nids[3] = nodes[20]->Id();
    nids[4] = nodes[21]->Id();
    nids[5] = nodes[13]->Id();
    nids[6] = nodes[22]->Id();
    nids[7] = nodes[26]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[20]->Id();
    nids[1] = nodes[ 9]->Id();
    nids[2] = nodes[ 2]->Id();
    nids[3] = nodes[10]->Id();
    nids[4] = nodes[26]->Id();
    nids[5] = nodes[22]->Id();
    nids[6] = nodes[14]->Id();
    nids[7] = nodes[23]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[11]->Id();
    nids[1] = nodes[20]->Id();
    nids[2] = nodes[10]->Id();
    nids[3] = nodes[ 3]->Id();
    nids[4] = nodes[24]->Id();
    nids[5] = nodes[26]->Id();
    nids[6] = nodes[23]->Id();
    nids[7] = nodes[15]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    /////////////////////////////////////////////////////////////////

    nids[0] = nodes[12]->Id();
    nids[1] = nodes[21]->Id();
    nids[2] = nodes[26]->Id();
    nids[3] = nodes[24]->Id();
    nids[4] = nodes[ 4]->Id();
    nids[5] = nodes[16]->Id();
    nids[6] = nodes[25]->Id();
    nids[7] = nodes[19]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[21]->Id();
    nids[1] = nodes[13]->Id();
    nids[2] = nodes[22]->Id();
    nids[3] = nodes[26]->Id();
    nids[4] = nodes[16]->Id();
    nids[5] = nodes[ 5]->Id();
    nids[6] = nodes[17]->Id();
    nids[7] = nodes[25]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[26]->Id();
    nids[1] = nodes[22]->Id();
    nids[2] = nodes[14]->Id();
    nids[3] = nodes[23]->Id();
    nids[4] = nodes[25]->Id();
    nids[5] = nodes[17]->Id();
    nids[6] = nodes[ 6]->Id();
    nids[7] = nodes[18]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );

    nids[0] = nodes[24]->Id();
    nids[1] = nodes[26]->Id();
    nids[2] = nodes[23]->Id();
    nids[3] = nodes[15]->Id();
    nids[4] = nodes[19]->Id();
    nids[5] = nodes[25]->Id();
    nids[6] = nodes[18]->Id();
    nids[7] = nodes[ 7]->Id();
    subelements_.push_back( mesh.GetElement( -1, nids, *top_data ) );
  }
}

#ifdef QHULL
void GEO::CUT::ConcreteElement<DRT::Element::tet4>::FillTetgen( tetgenio & out )
{
  const int dim = 3;

  out.numberofpoints = 4;
  out.pointlist = new double[out.numberofpoints * dim];
  out.pointmarkerlist = new int[out.numberofpoints];
  std::fill( out.pointmarkerlist, out.pointmarkerlist+out.numberofpoints, 0 );

  out.numberoftrifaces = 4;
  out.trifacemarkerlist = new int[out.numberoftrifaces];
  out.trifacelist = new int[out.numberoftrifaces * dim];

  out.numberoftetrahedra = 1;
  out.tetrahedronlist = new int[out.numberoftetrahedra * 4];
  //out.tetrahedronmarkerlist = new int[out.numberoftetrahedra];
  //std::fill( out.tetrahedronmarkerlist, out.tetrahedronmarkerlist+out.numberoftetrahedra, 0 );

  const std::vector<Node*> & nodes = Nodes();
  for ( int i=0; i<4; ++i )
  {
    Node * n = nodes[i];
    n->Coordinates( &out.pointlist[i*dim] );
    out.pointmarkerlist[i] = n->point()->Position();
    out.tetrahedronlist[i] = i;
  }

  const std::vector<Side*> & sides = Sides();
  for ( int i=0; i<4; ++i )
  {
    Side * s = sides[i];
    const std::vector<Node*> & side_nodes = s->Nodes();

    for ( int j=0; j<3; ++j )
    {
      out.trifacelist[i*dim+j] = std::find( nodes.begin(), nodes.end(), side_nodes[j] ) - nodes.begin();
    }

    int sid = s->Id();
    if ( sid < 0 )
    {
      for ( int j=0; j<3; ++j )
      {
        sid = std::min( sid, out.pointmarkerlist[out.trifacelist[i*dim+j]] );
      }
    }
    out.trifacemarkerlist[i] = sid;
  }
}
#endif

bool GEO::CUT::ConcreteElement<DRT::Element::tet4>::PointInside( Point* p )
{
  Position<DRT::Element::tet4> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::hex8>::PointInside( Point* p )
{
  Position<DRT::Element::hex8> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::wedge6>::PointInside( Point* p )
{
  Position<DRT::Element::wedge6> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::pyramid5>::PointInside( Point* p )
{
  Position<DRT::Element::pyramid5> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::hex20>::PointInside( Point* p )
{
  Position<DRT::Element::hex20> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::hex27>::PointInside( Point* p )
{
  Position<DRT::Element::hex27> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::tet10>::PointInside( Point* p )
{
  Position<DRT::Element::tet10> pos( *this, *p );
  return pos.Compute();
}


void GEO::CUT::ConcreteElement<DRT::Element::tet4>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::tet4> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::hex8>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex8> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::wedge6>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::wedge6> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::pyramid5>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::pyramid5> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::hex20>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex20> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::hex27>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex27> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::tet10>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::tet10> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

