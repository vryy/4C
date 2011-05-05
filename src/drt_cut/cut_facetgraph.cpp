
#include <iterator>

#include "cut_facetgraph.H"
#include "cut_facet.H"
#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_element.H"

GEO::CUT::FacetGraph::FacetGraph( const std::vector<Side*> & sides, const std::set<Facet*> & facets )
  : graph_( facets.size() )
{
  all_facets_.reserve( facets.size() );
  all_facets_.assign( facets.begin(), facets.end() );

  std::map<std::pair<Point*, Point*>, int> all_lines;

  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;

    std::map<std::pair<Point*, Point*>, std::set<Facet*> > lines;
    f->GetLines( lines );

    for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=lines.begin(); i!=lines.end(); ++i )
    {
      const std::pair<Point*, Point*> & l = i->first;
      std::map<std::pair<Point*, Point*>, int>::iterator j = all_lines.find( l );
      if ( j==all_lines.end() )
      {
        int current = all_facets_.size() + all_lines.size();
        all_lines[l] = current;
      }

      graph_.Add( FacetId( f ), all_lines[l] );
    }
  }

  all_lines_.resize( all_lines.size() );
  for ( std::map<std::pair<Point*, Point*>, int>::iterator i=all_lines.begin(); i!=all_lines.end(); ++i )
  {
    all_lines_[i->second - all_facets_.size()] = i->first;
  }

  COLOREDGRAPH::Graph cycle( facets.size() );

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    const std::vector<Facet*> & side_facets = s->Facets();

    for ( std::vector<Facet*>::const_iterator i=side_facets.begin(); i!=side_facets.end(); ++i )
    {
      Facet * f = *i;

      std::map<std::pair<Point*, Point*>, std::set<Facet*> > lines;
      f->GetLines( lines );

      for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=lines.begin(); i!=lines.end(); ++i )
      {
        const std::pair<Point*, Point*> & l = i->first;
        cycle.Add( FacetId( f ), all_lines[l] );
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
  for ( COLOREDGRAPH::CycleList::iterator i=cycle_list_.begin(); i!=cycle_list_.end(); ++i )
  {
    COLOREDGRAPH::Graph & g = *i;

    std::set<Facet*> collected_facets;

    for ( COLOREDGRAPH::Graph::const_iterator i=g.begin(); i!=g.end(); ++i )
    {
      int p = i->first;
      if ( p >= g.Split() )
      {
        break;
      }
      collected_facets.insert( all_facets_[p] );
    }

    std::map<std::pair<Point*, Point*>, std::set<Facet*> > volume_lines;
    for ( std::set<Facet*>::iterator i=collected_facets.begin();
          i!=collected_facets.end();
          ++i )
    {
      Facet * f = *i;
      f->GetLines( volume_lines );
    }

    cells.insert( mesh.NewVolumeCell( collected_facets, volume_lines, element ) );
  }
}

