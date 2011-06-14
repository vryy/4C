
#include <iterator>

#include "cut_facetgraph.H"
#include "cut_facet.H"
#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_element.H"

GEO::CUT::FacetGraph::FacetGraph( const std::vector<Side*> & sides, const plain_facet_set & facets )
  : graph_( facets.size() )
{

  std::map<std::pair<Point*, Point*>, plain_facet_set > lines;
  for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    f->GetLines( lines );

    // connect side facet to all hole lines

    if ( f->HasHoles() )
    {
      const plain_facet_set & holes = f->Holes();

      std::map<std::pair<Point*, Point*>, plain_facet_set > hole_lines;
      for ( plain_facet_set::const_iterator i=holes.begin(); i!=holes.end(); ++i )
      {
        Facet * h = *i;
        h->GetLines( hole_lines );
      }
      for ( std::map<std::pair<Point*, Point*>, plain_facet_set >::iterator i=hole_lines.begin();
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

  for ( std::map<std::pair<Point*, Point*>, plain_facet_set >::iterator li=lines.begin(); li!=lines.end(); )
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

  graph_.SetSplit( all_facets.size() );

  all_facets_.reserve( all_facets.size() );
  all_facets_.assign( all_facets.begin(), all_facets.end() );

  all_lines_.reserve( lines.size() );

  for ( std::map<std::pair<Point*, Point*>, plain_facet_set >::iterator i=lines.begin(); i!=lines.end(); ++i )
  {
    const std::pair<Point*, Point*> & l = i->first;
    const plain_facet_set & fs = i->second;

    int current = all_facets_.size() + all_lines_.size();
    all_lines_.push_back( l );

    for ( plain_facet_set::const_iterator i=fs.begin(); i!=fs.end(); ++i )
    {
      Facet * f = *i;
      if ( all_facets.count( f ) > 0 )
      {
        graph_.Add( FacetId( f ), current );
      }
    }
  }

#if 0
  graph_.FixSingleLines();
#else
  graph_.TestClosed();
#endif

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
        plain_int_set & row = graph_[p1];
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int p2 = *i;
          cycle.Add( p1, p2 );
        }
      }

      if ( f->HasHoles() )
      {
        const plain_facet_set & hs = f->Holes();
        for ( plain_facet_set::iterator i=hs.begin(); i!=hs.end(); ++i )
        {
          Facet * h = *i;

          if ( all_facets.count( h ) > 0 )
          {
            int p1 = FacetId( h );
            plain_int_set & row = graph_[p1];
            for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
            {
              int p2 = *i;
              cycle.Add( p1, p2 );
            }
          }
        }
      }
    }
  }

#if 0
#ifdef DEBUGCUTLIBRARY
  std::cout << "graph:\n";
  graph_.Print();

  std::cout << "cycle:\n";
  cycle.Print();
#endif
#endif

  plain_int_set free;
  graph_.GetAll( free );

  for ( COLOREDGRAPH::Graph::const_iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    int p = i->first;
    free.erase( p );
  }

#if 0
#ifdef DEBUGCUTLIBRARY
  std::cout << "free: ";
  std::copy( free.begin(), free.end(), std::ostream_iterator<int>( std::cout, " " ) );
  std::cout << "\n";

  for ( std::vector<Facet*>::iterator i=all_facets_.begin(); i!=all_facets_.end(); ++i )
  {
    Facet * f = *i;
    std::stringstream str;
    str << "facet-" << std::distance( all_facets_.begin(), i ) << ".plot";
    std::cout << str.str() << "\n";
    std::ofstream file( str.str().c_str() );
    f->Print( file );
  }
#endif
#endif

  COLOREDGRAPH::Graph used( cycle );
  cycle_list_.AddPoints( graph_, used, cycle, free );
}

void GEO::CUT::FacetGraph::CreateVolumeCells( Mesh & mesh, Element * element, plain_volumecell_set & cells )
{
  std::vector<plain_facet_set> volumes;
  volumes.reserve( cycle_list_.size() );

  std::vector<int> counter( all_facets_.size(), 0 );

  for ( COLOREDGRAPH::CycleList::iterator i=cycle_list_.begin(); i!=cycle_list_.end(); ++i )
  {
    COLOREDGRAPH::Graph & g = *i;

#ifdef DEBUGCUTLIBRARY
    g.TestSplit();
#endif

    volumes.push_back( plain_facet_set() );
    plain_facet_set & collected_facets = volumes.back();

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

//   std::cout << "generated volumes " << volumes.size() << " vs " << cycle_list_.size() << "\n";

//   if ( cycle_list_.size()==1 )
//   {
//     graph_.Print();
//   }

  for ( std::vector<plain_facet_set>::iterator i=volumes.begin(); i!=volumes.end(); ++i )
  {
    plain_facet_set & collected_facets = *i;

    std::map<std::pair<Point*, Point*>, plain_facet_set> volume_lines;
    for ( plain_facet_set::iterator i=collected_facets.begin();
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

bool GEO::CUT::FacetGraph::InCollectedFacets( const plain_facet_set & collected_facets )
{
  return std::find( all_collected_facets_.begin(), all_collected_facets_.end(), collected_facets )!=all_collected_facets_.end();
}

void GEO::CUT::FacetGraph::PrintAllCollected()
{
  int count = 0;
  for ( std::vector<plain_facet_set>::iterator i=all_collected_facets_.begin(); i!=all_collected_facets_.end(); ++i )
  {
    plain_facet_set & v = *i;

    count += 1;

    std::stringstream str;
    str << "gv-" << count << ".plot";
    std::ofstream file( str.str().c_str() );

    for ( plain_facet_set::iterator i=v.begin(); i!=v.end(); ++i )
    {
      Facet * f = *i;
      f->Print( file );
    }
  }
}

#endif
