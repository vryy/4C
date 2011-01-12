
#include "cut_volumecell.H"
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
}

bool GEO::CUT::VolumeCell::Contains( Point * p )
{
  for ( std::set<Facet*>::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->Contains( p ) )
    {
      return true;
    }
  }
  return false;
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

void GEO::CUT::VolumeCell::CreateIntegrationCells( Mesh & mesh )
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
    Tet4IntegrationCell::CreateCells( mesh, element_, this, facets_, integrationcells_ );
  }
}

void GEO::CUT::VolumeCell::GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells )
{
  std::copy( integrationcells_.begin(), integrationcells_.end(), std::inserter( cells, cells.begin() ) );
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
