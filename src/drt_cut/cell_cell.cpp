
#include "../linalg/linalg_fixedsizematrix.H"

#include "cut_facet.H"
#include "cut_point.H"
#include "cut_element.H"

#include "cell_cell.H"


GEO::CELL::Cell::Cell( GEO::CUT::LinearElement * element )
  : element_( element )
{

}

void GEO::CELL::Cell::AddVolume( const std::set<GEO::CUT::Facet*> & collected_facets,
                                 const std::map<std::pair<GEO::CUT::Point*, GEO::CUT::Point*>, std::set<GEO::CUT::Facet*> > & volume_lines )
{
  for ( std::map<std::pair<GEO::CUT::Point*, GEO::CUT::Point*>, std::set<GEO::CUT::Facet*> >::const_iterator i=volume_lines.begin();
        i!=volume_lines.end();
        ++i )
  {
    const std::pair<GEO::CUT::Point*, GEO::CUT::Point*> & ps = i->first;
    Point* p1 = NewPoint( ps.first );
    Point* p2 = NewPoint( ps.second );
    //Line * l =
    NewLine( p1, p2 );
  }

  Volume * v = NewVolume();

  std::map<int, FacetSet*> fs;

  for ( std::set<GEO::CUT::Facet*>::const_iterator i=collected_facets.begin();
        i!=collected_facets.end();
        ++i )
  {
    GEO::CUT::Facet * f = *i;
    std::set<Facet *> nf;
    NewFacets( f, nf );
    int sid = f->PositionSideId();
    if ( fs.count( sid )==0 )
    {
      fs[sid] = v->NewFacetSet( sid );
    }
    fs[sid]->AddFacets( nf );
  }
}

GEO::CELL::Point* GEO::CELL::Cell::NewPoint( const GEO::CUT::Point* p )
{
  std::map<const GEO::CUT::Point*, Teuchos::RCP<Point> >::iterator i = points_.find( p );
  if ( i==points_.end() )
  {
    LINALG::Matrix<3,1> xyz( p->X() );
    LINALG::Matrix<3,1> rst;
    element_->LocalCoordinates( xyz, rst );
    Point * np = new Point( rst.A() );
    points_[p] = Teuchos::rcp( np );
    return np;
  }
  else
  {
    return &*i->second;
  }
}

GEO::CELL::Line* GEO::CELL::Cell::NewLine( Point* p1, Point* p2 )
{
  Line * l = p1->FindLine( p2 );
  if ( l==NULL )
  {
    l = new Line( p1, p2 );
    lines_.push_back( Teuchos::rcp( l ) );
  }
  return l;
}

GEO::CELL::Facet* GEO::CELL::Cell::NewFacet( GEO::CUT::Facet * f, const std::vector<Point*> & ps )
{
  unsigned length = ps.size();
  std::vector<Line*> ls;
  ls.reserve( length );

  for ( unsigned i=0; i<length; ++i )
  {
    unsigned j = ( i+1 ) % length;
    ls.push_back( NewLine( ps[i], ps[j] ) );
  }

  Facet * nf = new Facet( ps, ls );
  facets_[f].push_back( Teuchos::rcp( nf ) );
  return nf;
}

void GEO::CELL::Cell::NewFacets( GEO::CUT::Facet * f, std::set<Facet *> & nf )
{
  std::map<const GEO::CUT::Facet*, std::vector<Teuchos::RCP<Facet> > >::iterator i = facets_.find( f );
  if ( i!=facets_.end() )
  {
    std::vector<Teuchos::RCP<Facet> > & fs = i->second;
    for ( std::vector<Teuchos::RCP<Facet> >::iterator i=fs.begin(); i!=fs.end(); ++i )
    {
      Facet * mf = &**i;
      nf.insert( mf );
    }
  }
  else
  {
    f->CreateFacets( this, nf );
  }
}

GEO::CELL::Volume* GEO::CELL::Cell::NewVolume()
{
  Volume * v = new Volume;
  volumes_.push_back( Teuchos::rcp( v ) );
  return v;
}
