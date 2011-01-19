
#include "cut_volumecell.H"
#include "cut_boundarycell.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"

GEO::CUT::VolumeCell::VolumeCell( const std::set<Facet*> & facets,
                                  const std::map<std::pair<Point*, Point*>, std::set<Facet*> > & volume_lines,
                                  Element * element )
  : element_( element ),
    facets_( facets )
{
/*

Das Ziel ist es einfache Volumen zu bekommen, die mit bekannten
Elementansätzen integriert werden können.

- Die Volumenzellen in Elementkoordinaten formulieren
- Die Facets jeder Schnittebene werden zu einem Facetset zusammengefaßt.
- Die Facets der Schnittebenen werden, wenn sie nicht tri3 oder quad4 sind in
  Dreiecke zerlegt.

Offene Fragen:

- In welche Richtung sollten die Ebenen angesehen werden? Wie kann man das
  Schneiden der Schnittfläche durch die Ebenen vermeiden?
- Was passiert mit (nahezu) achsparallelen Schnittflächen?
*/

//   for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::const_iterator i=volume_lines.begin();
//         i!=volume_lines.end();
//         ++i )
//   {
//     const std::pair<Point*, Point*> & ps = i->first;
//   }

  for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Register( this );
  }

  std::set<Point*> cut_points;
  GetAllPoints( cut_points );
}

// bool GEO::CUT::VolumeCell::Contains( Point * p )
// {
//   for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
//   {
//     Facet * f = *i;
//     if ( f->Contains( p ) )
//     {
//       return true;
//     }
//   }
//   return false;
// }

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

void GEO::CUT::VolumeCell::GetAllPoints( std::set<Point*> & cut_points )
{
  if ( points_.size()==0 )
  {
    for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      f->GetAllPoints( cut_points );
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
#if 0
  IntegrationCell * ic;

  if ( ( ic = Hex8IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
  {
    integrationcells_.insert( ic );
  }
  else if ( ( ic = Tet4IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
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
  else
#endif
  {
    Tet4IntegrationCell::CreateCells( mesh, element_, this, facets_, integrationcells_ );
  }
}

void GEO::CUT::VolumeCell::GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells )
{
  std::copy( integrationcells_.begin(), integrationcells_.end(), std::inserter( cells, cells.begin() ) );
}

void GEO::CUT::VolumeCell::SetIntegrationCells( Mesh & mesh,
                                                const std::vector<std::vector<std::vector<Point*> >::iterator> & tets,
                                                const std::map<std::vector<std::vector<Point*> >::iterator, std::vector<std::vector<Point*> > > & cut_sides )
{
  IntegrationCell * ic;

#if 0
  if ( ( ic = Hex8IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
  {
    integrationcells_.insert( ic );
  }
  else if ( ( ic = Tet4IntegrationCell::CreateCell( mesh, this, facets_ ) )!=NULL )
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
  else
#endif
  {
    Point::PointPosition position = Position();

    std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > sides_xyz;

    for ( std::vector<std::vector<std::vector<Point*> >::iterator>::const_iterator i=tets.begin();
          i!=tets.end();
          ++i )
    {
      std::vector<std::vector<Point*> >::iterator it = *i;
      std::vector<Point*> & t = *it;

      ic = Tet4IntegrationCell::CreateCell( mesh, this, position, t );
      integrationcells_.insert( ic );

      std::map<std::vector<std::vector<Point*> >::iterator, std::vector<std::vector<Point*> > >::const_iterator j=cut_sides.find( it );
      if ( j!=cut_sides.end() )
      {
        const std::vector<std::vector<Point*> > & sides = j->second;
        for ( std::vector<std::vector<Point*> >::const_iterator i=sides.begin(); i!=sides.end(); ++i )
        {
          const std::vector<Point*> & side = *i;
          Tri3BoundaryCell::CollectCoordinates( side, sides_xyz );
        }
      }
    }

    BoundaryCell::CreateCells( mesh, this, sides_xyz );
  }
}

void GEO::CUT::VolumeCell::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  std::copy( bcells_.begin(), bcells_.end(), std::inserter( bcells, bcells.begin() ) );
}

void GEO::CUT::VolumeCell::ConnectNodalDOFSets()
{
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

bool GEO::CUT::VolumeCell::OwnsFacet( Facet * f )
{
  return std::find( facets_.begin(), facets_.end(), f )!=facets_.end();
}

bool GEO::CUT::VolumeCell::OwnsTet( const std::vector<Point*> & tet )
{
  for ( std::vector<Point*>::const_iterator j=tet.begin(); j!=tet.end(); ++j )
  {
    Point * p = *j;
    if ( not Contains( p ) )
    {
      return false;
    }
  }

  // So we know all the nodes. But am I concave? Otherwise this tet might be
  // just outside of me.

#if 0
  for ( int i=0; i<4; ++i )
  {
    std::vector<Point*> side( tet );
    side.erase( side.begin() + i );

    bool all_oncutsurface = true;
    for ( std::vector<Point*>::iterator i=side.begin(); i!=side.end(); ++i )
    {
      Point * p = *i;
      if ( p->Position() != Point::oncutsurface )
      {
        all_oncutsurface = false;
        break;
      }
    }

    if ( all_oncutsurface )
      continue;

    bool found = false;
    for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Contains( side ) )
      {
        found = true;
        break;
      }
    }
    if ( not found )
    {
      return false;
    }
  }
#endif

  return true;
}

bool GEO::CUT::VolumeCell::Contains( Point * p )
{
  return std::binary_search( points_.begin(), points_.end(), p );
}

void GEO::CUT::VolumeCell::NewTri3Cell( Mesh & mesh, Facet * f, const Epetra_SerialDenseMatrix & x )
{
  f->NewTri3Cell( mesh, this, x, bcells_ );
}

void GEO::CUT::VolumeCell::NewTri3Cells( Mesh & mesh, Facet * f, const std::vector<Epetra_SerialDenseMatrix> & xyz )
{
  f->NewTri3Cells( mesh, this, xyz, bcells_ );
}

void GEO::CUT::VolumeCell::NewQuad4Cell( Mesh & mesh, Facet * f, const Epetra_SerialDenseMatrix & x )
{
  f->NewQuad4Cell( mesh, this, x, bcells_ );
}
