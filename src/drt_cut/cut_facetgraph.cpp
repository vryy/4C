
#include <iterator>

#include "cut_facetgraph.H"
#include "cut_facet.H"
#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_element.H"

GEO::CUT::FacetGraph::FacetGraph( const std::vector<Side*> & sides, const std::set<Facet*> & facets )
  : graph_( facets.size() )
{

  std::map<std::pair<Point*, Point*>, std::set<Facet*> > lines;
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    f->GetLines( lines );

    // connect side facet to all hole lines

    if ( f->HasHoles() )
    {
      const std::set<Facet*> & holes = f->Holes();

      std::map<std::pair<Point*, Point*>, std::set<Facet*> > hole_lines;
      for ( std::set<Facet*>::const_iterator i=holes.begin(); i!=holes.end(); ++i )
      {
        Facet * h = *i;
        h->GetLines( hole_lines );
      }
      for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=hole_lines.begin();
            i!=hole_lines.end();
            ++i )
      {
        const std::pair<Point*, Point*> & l = i->first;
        lines[l].insert( f );
      }
    }
  }

  std::set<Facet*> all_facets = facets;
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->HasHoles() )
    {
      const std::set<Facet*> & holes = f->Holes();
      std::copy( holes.begin(), holes.end(), std::inserter( all_facets, all_facets.begin() ) );
    }
  }

  // fix for very rare case

  for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator li=lines.begin(); li!=lines.end(); )
  {
    std::set<Facet*> & fs = li->second;
    if ( fs.size() < 2 )
    {
      for ( std::set<Facet*>::iterator i=fs.begin(); i!=fs.end(); ++i )
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

  graph_.SetSplit( all_facets.size() );

  all_facets_.reserve( all_facets.size() );
  all_facets_.assign( all_facets.begin(), all_facets.end() );

  all_lines_.reserve( lines.size() );

  for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=lines.begin(); i!=lines.end(); ++i )
  {
    const std::pair<Point*, Point*> & l = i->first;
    const std::set<Facet*> & fs = i->second;

    int current = all_facets_.size() + all_lines_.size();
    all_lines_.push_back( l );

    for ( std::set<Facet*>::const_iterator i=fs.begin(); i!=fs.end(); ++i )
    {
      Facet * f = *i;
      graph_.Add( FacetId( f ), current );
    }
  }

  graph_.FixSingleLines();

  COLOREDGRAPH::Graph cycle( all_facets_.size() );
  //COLOREDGRAPH::Graph holes( all_facets_.size() );

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    const std::vector<Facet*> & side_facets = s->Facets();

    for ( std::vector<Facet*>::const_iterator i=side_facets.begin(); i!=side_facets.end(); ++i )
    {
      Facet * f = *i;

      if ( all_facets.count( f ) > 0 )
      {
        int p1 = FacetId( f );
        std::set<int> & row = graph_[p1];
        for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int p2 = *i;
          cycle.Add( p1, p2 );
        }
      }

      if ( f->HasHoles() )
      {
        const std::set<Facet*> & hs = f->Holes();
        for ( std::set<Facet*>::iterator i=hs.begin(); i!=hs.end(); ++i )
        {
          Facet * h = *i;

          if ( all_facets.count( h ) > 0 )
          {
            int p1 = FacetId( h );
            std::set<int> & row = graph_[p1];
            for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
            {
              int p2 = *i;
              cycle.Add( p1, p2 );
            }
          }
        }
      }
    }
  }

//   std::cout << "cycle:\n";
//   cycle.Print();

  std::set<int> free;
  graph_.GetAll( free );

  for ( COLOREDGRAPH::Graph::const_iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    int p = i->first;
    free.erase( p );
  }

//   std::cout << "free: ";
//   std::copy( free.begin(), free.end(), std::ostream_iterator<int>( std::cout, " " ) );
//   std::cout << "\n";

  COLOREDGRAPH::Graph used( cycle );
  cycle_list_.AddPoints( graph_, used, cycle, free );
}

void GEO::CUT::FacetGraph::CreateVolumeCells( Mesh & mesh, Element * element, std::set<VolumeCell*> & cells )
{
  std::vector<std::set<Facet*> > volumes;
  volumes.reserve( cycle_list_.size() );

  std::vector<int> counter( all_facets_.size(), 0 );

  for ( COLOREDGRAPH::CycleList::iterator i=cycle_list_.begin(); i!=cycle_list_.end(); ++i )
  {
    COLOREDGRAPH::Graph & g = *i;

    volumes.push_back( std::set<Facet*>() );
    std::set<Facet*> & collected_facets = volumes.back();

    for ( COLOREDGRAPH::Graph::const_iterator i=g.begin(); i!=g.end(); ++i )
    {
      int p = i->first;
      if ( p >= g.Split() )
      {
        break;
      }
      collected_facets.insert( all_facets_[p] );
      counter[p] += 1;
    }
  }

  for ( std::vector<int>::iterator i=counter.begin(); i!=counter.end(); ++i )
  {
    int p = *i;
    if ( p > 2 )
    {
      // assume this is a degenerated cell we do not need to build
      return;
    }
  }

  for ( std::vector<std::set<Facet*> >::iterator i=volumes.begin(); i!=volumes.end(); ++i )
  {
    std::set<Facet*> & collected_facets = *i;

    std::map<std::pair<Point*, Point*>, std::set<Facet*> > volume_lines;
    for ( std::set<Facet*>::iterator i=collected_facets.begin();
          i!=collected_facets.end();
          ++i )
    {
      Facet * f = *i;
      f->GetLines( volume_lines );
    }

    cells.insert( mesh.NewVolumeCell( collected_facets, volume_lines, element ) );
#ifdef DEBUGCUTLIBRARY
    all_collected_facets_.push_back( collected_facets );
#endif
  }

#if 0
  for ( std::vector<int>::iterator i=debug_counter.begin(); i!=debug_counter.end(); ++i )
  {
    int p = *i;
    if ( p > 2 )
    {
#ifdef DEBUGCUTLIBRARY
      PrintAllCollected();
#endif
      throw std::runtime_error( "failed to split graph properly" );
    }
  }
#endif
}

#ifdef DEBUGCUTLIBRARY

bool GEO::CUT::FacetGraph::InCollectedFacets( const std::set<Facet*> & collected_facets )
{
  return std::find( all_collected_facets_.begin(), all_collected_facets_.end(), collected_facets )!=all_collected_facets_.end();
}

void GEO::CUT::FacetGraph::PrintAllCollected()
{
  int count = 0;
  for ( std::vector<std::set<Facet*> >::iterator i=all_collected_facets_.begin(); i!=all_collected_facets_.end(); ++i )
  {
    std::set<Facet*> & v = *i;

    count += 1;

    std::stringstream str;
    str << "gv-" << count << ".plot";
    std::ofstream file( str.str().c_str() );

    for ( std::set<Facet*>::iterator i=v.begin(); i!=v.end(); ++i )
    {
      Facet * f = *i;
      f->Print( file );
    }
  }
}

#endif
