
#include "cut_volumecell.H"
#include "cut_boundarycell.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_tetmesh.H"


GEO::CUT::VolumeCell::VolumeCell( const std::set<Facet*> & facets,
                                  const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
                                  Element * element )
  : element_( element ),
    facets_( facets )
{
  for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Register( this );
  }
}

void GEO::CUT::VolumeCell::Neighbors( Point * p,
                                      const std::set<VolumeCell*> & cells,
                                      const std::set<VolumeCell*> & done,
                                      std::set<VolumeCell*> & connected,
                                      std::set<Element*> & elements )
{
  if ( done.count( this )==0 )
  {
    // this volume is included
    connected.insert( this );
    elements.insert( element_ );

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      if ( p==NULL or f->Contains( p ) )
      {
        f->Neighbors( p, cells, done, connected, elements );
      }
    }

    if ( p!=NULL )
    {
      for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
      {
        Facet * f = *i;
        if ( not f->Contains( p ) )
        {
          f->Neighbors( p, cells, done, connected, elements );
        }
      }
    }
  }
}

void GEO::CUT::VolumeCell::GetAllPoints( Mesh & mesh, std::set<Point*> & cut_points )
{
  if ( points_.size()==0 )
  {
    for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      f->GetAllPoints( mesh, cut_points );
    }
    points_.reserve( cut_points.size() );
    std::copy( cut_points.begin(), cut_points.end(), std::back_inserter( points_ ) );
  }
  else
  {
    std::copy( points_.begin(), points_.end(), std::inserter( cut_points, cut_points.begin() ) );
  }
}

void GEO::CUT::VolumeCell::CreateIntegrationCells( Mesh & mesh )
{
  if ( element_->IsCut() )
  {
    IntegrationCell * ic;

    if ( ( ic = Tet4IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
    {
      integrationcells_.insert( ic );
    }
#if 0
    else if ( ( ic = Hex8IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
    {
      integrationcells_.insert( ic );
    }
    else if ( ( ic = Wedge6IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
    {
      integrationcells_.insert( ic );
    }
    else if ( ( ic = Pyramid5IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
    {
      integrationcells_.insert( ic );
    }
#endif
    else
    {
      std::set<Point*> cut_points;
      GetAllPoints( mesh, cut_points );

      std::vector<Point*> points;
      points.reserve( cut_points.size() );
      points.assign( cut_points.begin(), cut_points.end() );

      // sort points that go into qhull to obtain the same result independent of
      // pointer values (compiler flags, code structure, memory usage, ...)
      std::sort( points.begin(), points.end(), PointPidLess() );

#ifdef DEBUGCUTLIBRARY
      {
        std::ofstream file( "volume.plot" );
        Print( file );
      }
#endif

      TetMesh tetmesh( points, facets_ );

      const std::vector<std::vector<int> > & tets = tetmesh.Tets();
      const std::map<Facet*, std::vector<Point*> > & sides_xyz = tetmesh.SidesXYZ();

      Point::PointPosition position = Position();

      for ( std::vector<std::vector<int> >::const_iterator i=tets.begin();
            i!=tets.end();
            ++i )
      {
        const std::vector<int> & t = *i;
        if ( t.size()==4 )
        {
          std::vector<Point*> tet( 4 );
          for ( int i=0; i<4; ++i )
          {
            tet[i] = points[t[i]];
          }
          IntegrationCell * ic = Tet4IntegrationCell::CreateCell( mesh, this, position, tet );
          integrationcells_.insert( ic );
        }
      }

      for ( std::map<Facet*, std::vector<Point*> >::const_iterator i=sides_xyz.begin();
            i!=sides_xyz.end();
            ++i )
      {
        Facet * f = i->first;
        const std::vector<Point*> & points = i->second;

        std::size_t length = points.size();
        if ( length % 3 != 0 )
          throw std::runtime_error( "expect list of triangles" );

        length /= 3;
        std::vector<Point*> p( 3 );
        for ( std::size_t i=0; i<length; ++i )
        {
          std::copy( &points[3*i], &points[3*( i+1 )], &p[0] );
          Tri3BoundaryCell::CreateCell( mesh, this, f, p );
        }
      }
    }
  }
  else
  {
    IntegrationCell * ic;
    switch ( element_->Shape() )
    {
    case DRT::Element::tet4:
      ic = Tet4IntegrationCell::CreateCell( mesh, this, facets_ );
      break;
    case DRT::Element::hex8:
      ic = Hex8IntegrationCell::CreateCell( mesh, this, facets_ );
      break;
    case DRT::Element::wedge6:
      ic = Wedge6IntegrationCell::CreateCell( mesh, this, facets_ );
      break;
    case DRT::Element::pyramid5:
      ic = Pyramid5IntegrationCell::CreateCell( mesh, this, facets_ );
      break;
    default:
      throw std::runtime_error( "unsupported element shape" );
    }

    if ( ic!=NULL )
    {
      integrationcells_.insert( ic );
    }
    else
    {
      throw std::runtime_error( "failed to recognize uncut element shape" );
    }
  }
}

void GEO::CUT::VolumeCell::GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells )
{
  std::copy( integrationcells_.begin(), integrationcells_.end(), std::inserter( cells, cells.begin() ) );
}

void GEO::CUT::VolumeCell::GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells )
{
  for ( std::set<BoundaryCell*>::iterator i=bcells_.begin(); i!=bcells_.end(); ++i )
  {
    BoundaryCell * bc = *i;
    Facet * f = bc->GetFacet();
    int sid = f->SideId();
    if ( sid > -1 )
    {
      bcells[sid].push_back( bc );
    }
  }
}

void GEO::CUT::VolumeCell::ConnectNodalDOFSets( bool include_inner )
{
  if ( not include_inner and Position()!=Point::outside )
    return;

  const std::vector<Node*> & nodes = element_->Nodes();
  nodaldofset_.reserve( nodes.size() );

  for ( std::vector<Node*>::const_iterator i=nodes.begin();
        i!=nodes.end();
        ++i )
  {
    Node * n = *i;
    nodaldofset_.push_back( n->DofSetNumber( this ) );
  }
}

GEO::CUT::Point::PointPosition GEO::CUT::VolumeCell::Position()
{
  GEO::CUT::Point::PointPosition position = GEO::CUT::Point::undecided;
  for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    GEO::CUT::Point::PointPosition fp = f->Position();
    switch ( fp )
    {
    case GEO::CUT::Point::undecided:
      throw std::runtime_error( "undecided facet position" );
    case GEO::CUT::Point::oncutsurface:
      break;
    case GEO::CUT::Point::inside:
    case GEO::CUT::Point::outside:
      if ( position!=GEO::CUT::Point::undecided and position!=fp )
      {
        throw std::runtime_error( "mixed facet set" );
      }
      position = fp;
    }
  }
  return position;
}

void GEO::CUT::VolumeCell::Print( std::ostream & stream )
{
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Print( stream );
  }
}

bool GEO::CUT::VolumeCell::Contains( Point * p )
{
  if ( points_.size() > 0 )
  {
    return std::binary_search( points_.begin(), points_.end(), p );
  }
  else
  {
    for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Contains( p ) )
      {
        return true;
      }
    }
  }
  return false;
}

void GEO::CUT::VolumeCell::NewTri3Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x )
{
  f->NewTri3Cell( mesh, this, x, bcells_ );
}

// void GEO::CUT::VolumeCell::NewTri3Cells( Mesh & mesh, Facet * f, const std::vector<Epetra_SerialDenseMatrix> & xyz )
// {
//   f->NewTri3Cells( mesh, this, xyz, bcells_ );
// }

void GEO::CUT::VolumeCell::NewQuad4Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x )
{
  f->NewQuad4Cell( mesh, this, x, bcells_ );
}

double GEO::CUT::VolumeCell::Volume()
{
  double volume = 0;
  for ( std::set<IntegrationCell*>::iterator i=integrationcells_.begin(); i!=integrationcells_.end(); ++i )
  {
    IntegrationCell * ic = *i;
    volume += ic->Volume();
  }
  return volume;
}
