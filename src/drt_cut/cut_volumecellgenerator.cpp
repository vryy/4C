
#include "cut_volumecellgenerator.H"

#include "cut_facet.H"
#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_element.H"

GEO::CUT::VolumeCellGenerator::VolumeCellGenerator( const std::vector<Side*> & sides,
                                                    const plain_facet_set & facets )
{
  std::map<std::pair<Point*, Point*>, plain_facet_set> lines;
  for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    f->GetLines( lines );

    // connect side facet to all hole lines

    if ( f->HasHoles() )
    {
      const plain_facet_set & holes = f->Holes();

      std::map<std::pair<Point*, Point*>, plain_facet_set> hole_lines;
      for ( plain_facet_set::const_iterator i=holes.begin(); i!=holes.end(); ++i )
      {
        Facet * h = *i;
        h->GetLines( hole_lines );
      }
      for ( std::map<std::pair<Point*, Point*>, plain_facet_set>::iterator i=hole_lines.begin();
            i!=hole_lines.end();
            ++i )
      {
        const std::pair<Point*, Point*> & l = i->first;
        lines[l].insert( f );
      }
    }
  }

  plain_facet_set all_facets = facets;
  for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->HasHoles() )
    {
      const plain_facet_set & holes = f->Holes();
      std::copy( holes.begin(), holes.end(), std::inserter( all_facets, all_facets.begin() ) );
    }
  }

  // fix for very rare case

  for ( std::map<std::pair<Point*, Point*>, plain_facet_set>::iterator li=lines.begin(); li!=lines.end(); )
  {
    plain_facet_set & fs = li->second;
    if ( fs.size() < 2 )
    {
      for ( plain_facet_set::iterator i=fs.begin(); i!=fs.end(); ++i )
      {
        Facet * f = *i;
        all_facets.erase( f );
      }
      lines.erase( li++ );
    }
    else
    {
      ++li;
    }
  }

  //

  // extract all points

  for ( plain_facet_set::const_iterator i=all_facets.begin(); i!=all_facets.end(); ++i )
  {
    Facet * f = *i;
    const std::vector<Point*> & points = f->Points();
    for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
    {
      Point * p = *i;
      points_[p].p_ = p;
    }
  }

  // extract lines and facets

  lines_.reserve( lines.size() ); // make sure we do not re-alloc!

  for ( std::map<std::pair<Point*, Point*>, plain_facet_set>::iterator i=lines.begin();
        i!=lines.end();
        ++i )
  {
    Point * p1 = i->first.first;
    Point * p2 = i->first.second;
    lines_.push_back( GenLine() );
    GenLine & line = lines_.back();
    line.SetPoints( &points_[p1], &points_[p2] );

    plain_facet_set & linefacets = i->second;
    for ( plain_facet_set::iterator i=linefacets.begin(); i!=linefacets.end(); ++i )
    {
      Facet * f = *i;

      // do not care if it exists already
      GenFacet & gf = facets_[f];
      gf.facet_ = f;
      gf.openvolumes_ = f->OnCutSide() ? 2 : 1;
      gf.AddLine( &line );
    }
  }

  // correct open volumes counter for cut facets at element sides (touching)

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    const std::vector<Facet*> & side_facets = s->Facets();
    for ( std::vector<Facet*>::const_iterator i=side_facets.begin(); i!=side_facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->OnCutSide() )
      {
        facets_[f].openvolumes_ = 1;
      }
    }
  }

#if 0
#ifdef DEBUGCUTLIBRARY
  for ( std::vector<GenLine>::iterator i=lines_.begin(); i!=lines_.end(); ++i )
  {
    GenLine & gl = *i;
    std::cout << gl.Facets().size() << " ";
  }
  std::cout << "\n";
#endif
#endif
}

void GEO::CUT::VolumeCellGenerator::CreateVolumeCells( Mesh & mesh,
                                                       Element * element,
                                                       plain_volumecell_set & cells )
{
  while ( not Done() )
  {
    bool found = false;
    for ( std::map<Point*, GenPoint>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      GenPoint & gp = i->second;
      if ( gp.TestCorner( volumes_ ) )
      {
        found = true;
        break;
      }
    }
    if ( not found )
    {
      throw std::runtime_error( "no new corner found" );
    }
  }

  for ( std::list<GenVolume>::iterator i=volumes_.begin(); i!=volumes_.end(); ++i )
  {
    GenVolume & gv = *i;

    plain_facet_set facets;
    std::map<std::pair<Point*, Point*>, plain_facet_set> volume_lines;

    const GenFacetSet & volume_facets = gv.Facets();
    for ( GenFacetSet::const_iterator i=volume_facets.begin();
          i!=volume_facets.end();
          ++i )
    {
      GenFacet * gf = *i;
      Facet * f = gf->facet_;
      facets.insert( f );
      f->GetLines( volume_lines );
    }

    cells.insert( mesh.NewVolumeCell( facets, volume_lines, element ) );
  }
}

bool GEO::CUT::VolumeCellGenerator::Done()
{
  for ( std::map<Facet*, GenFacet>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    GenFacet & gf = i->second;
    if ( not gf.Done() )
    {
      return false;
    }
  }
  return true;
}

bool GEO::CUT::VolumeCellGenerator::GenPoint::TestCorner( std::list<GenVolume> & volumes )
{
  if ( lines_.size()==3 )
  {
    std::vector<GenLine*> lines;
    lines.reserve( 3 );
    lines.assign( lines_.begin(), lines_.end() );

    GenFacet * f1 = lines[0]->CommonFacet( lines[1] );
    if ( f1!=NULL and f1->openvolumes_==1 )
    {
      GenFacet * f2 = lines[0]->CommonFacet( lines[2] );
      if ( f2!=NULL and f2->openvolumes_==1 )
      {
        GenFacet * f3 = lines[1]->CommonFacet( lines[2] );
        if ( f3!=NULL and f3->openvolumes_==1 )
        {
          volumes.push_back( GenVolume() );
          GenVolume & gv = volumes.back();
          gv.AddFacet( f1 );
          gv.AddFacet( f2 );
          gv.AddFacet( f3 );
          gv.Finalize();
          return true;
        }
      }
    }
  }
  return false;
}

void GEO::CUT::VolumeCellGenerator::GenLine::CommonFacets( GenLine * other, GenFacetSet & facets )
{
  std::set_intersection( facets_.begin(), facets_.end(),
                         other->facets_.begin(), other->facets_.end(),
                         std::inserter( facets, facets.begin() ) );
}

GEO::CUT::VolumeCellGenerator::GenFacet *
GEO::CUT::VolumeCellGenerator::GenLine::CommonFacet( GenLine * other )
{
  GenFacetSet facets;
  CommonFacets( other, facets );
  if ( facets.size()==1 )
  {
    return *facets.begin();
  }
  return NULL;
}

void GEO::CUT::VolumeCellGenerator::GenVolume::Finalize()
{
  GenVolumeFinalizer finalizer( *this );
  for ( GenFacetSet::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    GenFacet * gf = *i;
    finalizer.AddFacet( gf );
  }
  finalizer.CloseVolume();
}

void GEO::CUT::VolumeCellGenerator::GenVolumeFinalizer::AddFacet( GenFacet * f )
{
  volume_.AddFacet( f );
  const GenLineSet & lines = f->Lines();
  for ( GenLineSet::const_iterator i=lines.begin(); i!=lines.end(); ++i )
  {
    GenLine * gl = *i;
    lines_[gl].insert( f );
    if ( lines_[gl].size() > 2 )
    {
      throw std::runtime_error( "too many facets at line" );
    }
    points_[gl->P1()].insert( gl );
    points_[gl->P2()].insert( gl );
    if ( openlines_.count( gl ) > 0 )
    {
      openlines_.erase( gl );
      openpoints_[gl->P1()].erase( gl );
      if ( openpoints_[gl->P1()].size()==0 )
        openpoints_.erase( gl->P1() );
      openpoints_[gl->P2()].erase( gl );
      if ( openpoints_[gl->P2()].size()==0 )
        openpoints_.erase( gl->P2() );
    }
    else
    {
      openlines_.insert( gl );
      openpoints_[gl->P1()].insert( gl );
      openpoints_[gl->P2()].insert( gl );
    }
  }
  if ( openpoints_.size() != openlines_.size() )
  {
    std::stringstream str;
    str << openpoints_.size() << " open points vs. "
        << openlines_.size() << " open lines";
    throw std::runtime_error( str.str() );
  }
  for ( std::map<GenPoint*, GenLineSet>::iterator i=openpoints_.begin();
        i!=openpoints_.end();
        ++i )
  {
    GenLineSet & gls = i->second;
    if ( gls.size() != 2 )
    {
      throw std::runtime_error( "need closed cycle of open lines" );
    }
  }
}

void GEO::CUT::VolumeCellGenerator::GenVolumeFinalizer::CloseVolume()
{
  int failed_counter = 1000;
  while ( openlines_.size() > 0 )
  {
    bool found = false;
    for ( GenLineSet::iterator i=openlines_.begin(); i!=openlines_.end(); ++i )
    {
      GenLine * gl = *i;
      GenPoint * gp = gl->P1();
      GenFacet * gf = CommonFacet( gp );
      if ( gf == NULL )
      {
        gp = gl->P2();
        gf = CommonFacet( gp );
      }
      if ( gf != NULL )
      {
        AddFacet( gf );
        found = true;
        break;
      }
    }
    if ( not found )
    {
      throw std::runtime_error( "no unique line in volume cell creation" );
    }
    failed_counter -= 1;
    if ( failed_counter==0 )
    {
      throw std::runtime_error( "too many iterations" );
    }
  }
}

GEO::CUT::VolumeCellGenerator::GenFacet *
GEO::CUT::VolumeCellGenerator::GenVolumeFinalizer::CommonFacet( GenPoint * p )
{
  GenLineSet & gls = openpoints_[p];

  std::vector<GenLine*> lines;
  lines.reserve( 2 );
  lines.assign( gls.begin(), gls.end() );

  GenFacetSet facets;
  lines[0]->CommonFacets( lines[1], facets );
  for ( GenFacetSet::iterator i=facets.begin(); i!=facets.end(); )
  {
    GenFacet * f = *i;
    if ( volume_.HasFacet( f ) )
    {
      set_erase( facets, i );
    }
    else
    {
      ++i;
    }
  }
  if ( facets.size()==1 )
  {
    return *facets.begin();
  }
  return NULL;
}

