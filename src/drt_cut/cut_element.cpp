
//#include "../drt_geometry/intersection_templates.H"

#include "cut_intersection.H"
#include "cut_position.H"
#include "cut_facet.H"
#include "cut_volumecell.H"
#include "cut_tetmesh.H"
#include "cut_element.H"
#include "cut_options.H"
#include "cut_integrationcellcreator.H"
#include "cut_volumecellgenerator.H"
#include "cut_facetgraph.H"

#include <string>
#include <stack>

bool GEO::CUT::Element::Cut( Mesh & mesh, Side & side, int recursion )
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
    Side * s = *i;
    if ( FindCutPoints( mesh, *s, side, recursion ) )
    {
      cut = true;
    }
  }

  if ( cut )
  {
    cut_faces_.insert( &side );
    return true;
  }
  else
  {
    return false;
  }
}

void GEO::CUT::Element::MakeCutLines( Mesh & mesh, Creator & creator )
{
  for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
  {
    Side & side = **i;

    bool cut = false;

    const std::vector<Side*> & sides = Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side * s = *i;
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
      mesh.NewLinesBetween( line, &side, NULL, this );
    }

    if ( cut )
    {
      // create any remaining cut lines
      //side.CreateMissingLines( creator, this );
    }
  }
}

bool GEO::CUT::Element::FindCutPoints( Mesh & mesh, Side & side, Side & other, int recursion )
{
  bool cut = side.FindCutPoints( mesh, this, other, recursion );
  bool reverse_cut = other.FindCutPoints( mesh, this, side, recursion );
  return cut or reverse_cut;
}

bool GEO::CUT::Element::FindCutLines( Mesh & mesh, Side & side, Side & other )
{
  bool cut = side.FindCutLines( mesh, this, other );
  bool reverse_cut = other.FindCutLines( mesh, this, side );
  return cut or reverse_cut;
}

void GEO::CUT::Element::MakeFacets( Mesh & mesh )
{
  if ( facets_.size()==0 )
  {
    const std::vector<Side*> & sides = Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side & side = **i;
      side.MakeOwnedSideFacets( mesh, this, facets_ );
    }
    for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side & cut_side = **i;
      cut_side.MakeInternalFacets( mesh, this, facets_ );
    }
  }
}

void GEO::CUT::Element::FindNodePositions()
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
      const plain_facet_set & facets = p->Facets();
      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        for ( plain_side_set::const_iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
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
      const plain_facet_set & facets = p->Facets();
      for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        f->Position( pos );
      }
    }
  }
}

bool GEO::CUT::Element::IsCut()
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

bool GEO::CUT::Element::OnSide( Facet * f )
{
  if ( not f->HasHoles() )
  {
    return OnSide( f->Points() );
  }
  return false;
}

bool GEO::CUT::Element::OnSide( const std::vector<Point*> & facet_points )
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

  PointSet points;
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


void GEO::CUT::Element::GetIntegrationCells( plain_integrationcell_set & cells )
{
  for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = *i;
    vc->GetIntegrationCells( cells );
  }
}

void GEO::CUT::Element::GetBoundaryCells( plain_boundarycell_set & bcells )
{
  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( cut_faces_.count( f->ParentSide() )!= 0 )
    {
      f->GetBoundaryCells( bcells );
    }
  }
}

void GEO::CUT::Element::GetCutPoints( PointSet & cut_points )
{
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side * side = *i;

    for ( plain_side_set::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side * other = *i;
      side->GetCutPoints( this, *other, cut_points );
    }
  }
}

void GEO::CUT::Element::CreateIntegrationCells( Mesh & mesh, int count, bool levelset )
{
  if ( not active_ )
    return;

#ifdef DEBUGCUTLIBRARY
  {
    int volume_count = 0;
    for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); ++i )
    {
      VolumeCell * vc = *i;

      std::stringstream str;
      str << "volume-" << count << "-" << volume_count << ".plot";
      std::ofstream file( str.str().c_str() );
      vc->Print( file );
      volume_count += 1;
    }
  }
#endif

  if ( cells_.size()==1 )
  {
    VolumeCell * vc = *cells_.begin();
    if ( IntegrationCellCreator::CreateCell( mesh, Shape(), vc ) )
    {
      return;
    }
  }

  if ( mesh.CreateOptions().SimpleShapes() )
  {
    if ( IntegrationCellCreator::CreateCells( mesh, this, cells_ ) )
    {
      return;
    }
  }

  PointSet cut_points;

  // There are never holes in a cut facet. Furthermore, cut facets are
  // always convex, as all elements and sides are convex. Thus, we are free
  // to triangulate all cut facets. This needs to be done, so repeated cuts
  // work in the right way.

  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->OnCutSide() and f->HasHoles() )
      throw std::runtime_error( "no holes in cut facet possible" );
    //f->GetAllPoints( mesh, cut_points, f->OnCutSide() );
    f->GetAllPoints( mesh, cut_points, levelset and f->OnCutSide() );
  }

  std::vector<Point*> points;
  points.reserve( cut_points.size() );
  points.assign( cut_points.begin(), cut_points.end() );

  // sort points that go into qhull to obtain the same result independent of
  // pointer values (compiler flags, code structure, memory usage, ...)
  std::sort( points.begin(), points.end(), PointPidLess() );

  TetMesh tetmesh( points, facets_, false );
  tetmesh.CreateElementTets( mesh, this, cells_, cut_faces_, count, levelset );
}

void GEO::CUT::Element::RemoveEmptyVolumeCells()
{
  for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); )
  {
    VolumeCell * vc = *i;
    if ( vc->Empty() )
    {
      vc->Disconnect();
      set_erase( cells_, i );
    }
    else
    {
      ++i;
    }
  }
}

void GEO::CUT::Element::MakeVolumeCells( Mesh & mesh )
{
#if 0
#ifdef DEBUGCUTLIBRARY
  DumpFacets();
#endif
#endif

#if 0
  VolumeCellGenerator vcg( sides_, facets_ );
  vcg.CreateVolumeCells( mesh, this, cells_ );
#else
  FacetGraph fg( sides_, facets_ );
  fg.CreateVolumeCells( mesh, this, cells_ );
#endif
}

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

int GEO::CUT::Element::NumGaussPoints( DRT::Element::DiscretizationType shape )
{
  int numgp = 0;
  for ( plain_volumecell_set::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = *i;
    numgp += vc->NumGaussPoints( shape );
  }
  return numgp;
}

void GEO::CUT::Element::DebugDump()
{
  std::cout << "Problem in element " << Id() << " of shape " << Shape() << ":\n";
  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    //std::cout << n->LSV();
    n->Plot( std::cout );
  }
  std::cout << "\n";
  const plain_side_set & cutsides = CutSides();
  for ( plain_side_set::const_iterator i=cutsides.begin(); i!=cutsides.end(); ++i )
  {
    Side * s = *i;
    //s->Print();
    const std::vector<Node*> & side_nodes = s->Nodes();
    for ( std::vector<Node*>::const_iterator i=side_nodes.begin(); i!=side_nodes.end(); ++i )
    {
      Node * n = *i;
      n->Plot( std::cout );
    }
    std::cout << "\n";
  }
}

void GEO::CUT::Element::GnuplotDump()
{
  std::stringstream str;
  str << "element" << Id() << ".plot";
  std::ofstream file( str.str().c_str() );

  plain_edge_set all_edges;

  const std::vector<Side*> & sides = Sides();
  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    const std::vector<Edge*> & edges = s->Edges();

    std::copy( edges.begin(), edges.end(), std::inserter( all_edges, all_edges.begin() ) );
  }

  for ( plain_edge_set::iterator i=all_edges.begin(); i!=all_edges.end(); ++i )
  {
    Edge * e = *i;
    e->BeginNode()->point()->Plot( file );
    e->EndNode()  ->point()->Plot( file );
    file << "\n\n";
  }
}

void GEO::CUT::Element::DumpFacets()
{
  std::stringstream str;
  str << "facets" << Id() << ".plot";
  std::string name = str.str();

  std::cout << "write '" << name << "'\n";
  std::ofstream file( name.c_str() );

  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Print( file );
  }
}

