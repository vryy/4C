
#include <set>
#include <vector>
#include <iostream>
#include <iterator>

#include "../../src/drt_cut/cut_coloredgraph.H"
#include "../../src/drt_cut/cut_graph.H"
#include "cut_test_utils.H"

void test_graph()
{
  GEO::CUT::GRAPH::Graph g;

  g.Add( 1, 14 );
  g.Add( 1, 20 );
  g.Add( 1, 36 );

  g.Add( 14, 19 );
  g.Add( 14, 28 );
  g.Add( 14, 35 );

  g.Add( 19, 20 );
  g.Add( 19, 35 );
  g.Add( 19, 37 );

  g.Add( 20, 27 );
  g.Add( 20, 36 );
  g.Add( 20, 37 );

  g.Add( 26, 27 );
  g.Add( 26, 28 );

  g.Add( 35, 37 );

  g.Add( 36, 37 );

  g.Print();

  std::vector<int> cycle;
  cycle.push_back( 27 );
  cycle.push_back( 20 );
  cycle.push_back( 1 );
  cycle.push_back( 36 );
  cycle.push_back( 37 );
  cycle.push_back( 35 );
  cycle.push_back( 14 );
  cycle.push_back( 28 );
  cycle.push_back( 26 );

  std::set<int> free;
  g.GetAll( free );

  for ( std::vector<int>::iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    int p = *i;
    free.erase( p );
  }

  GEO::CUT::GRAPH::Graph used;
  used.AddCycle( cycle );

  GEO::CUT::GRAPH::CycleList facet_cycles;

  facet_cycles.AddPoints( g, used, cycle, free );

  for ( GEO::CUT::GRAPH::CycleListIterator i=facet_cycles.begin();
        i!=facet_cycles.end();
        ++i )
  {
    const std::vector<int> & c = *i;
    std::copy( c.begin(), c.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
  }
}

void test_graph2()
{
  GEO::CUT::GRAPH::Graph g;

  g.Add( 0, 10 );
  g.Add( 0, 11 );
  g.Add( 0, 23 );

  g.Add( 9, 10 );
  g.Add( 9, 11 );
  g.Add( 9, 23 );

  g.Add( 10, 20 );
  g.Add( 10, 23 );

  g.Add( 11, 21 );
  g.Add( 11, 23 );

  g.Add( 20, 22 );

  g.Add( 21, 22 );

  g.Print();

  std::vector<int> cycle;
  cycle.push_back( 10 );
  cycle.push_back(  0 );
  cycle.push_back( 11 );
  cycle.push_back( 21 );
  cycle.push_back( 22 );
  cycle.push_back( 20 );

  std::set<int> free;
  g.GetAll( free );

  for ( std::vector<int>::iterator i=cycle.begin(); i!=cycle.end(); ++i )
  {
    int p = *i;
    free.erase( p );
  }

  GEO::CUT::GRAPH::Graph used;
  used.AddCycle( cycle );

  GEO::CUT::GRAPH::CycleList facet_cycles;

  facet_cycles.AddPoints( g, used, cycle, free );

  for ( GEO::CUT::GRAPH::CycleListIterator i=facet_cycles.begin();
        i!=facet_cycles.end();
        ++i )
  {
    const std::vector<int> & c = *i;
    std::copy( c.begin(), c.end(), std::ostream_iterator<int>( std::cout, " " ) );
    std::cout << "\n";
  }
}

void test_colored_graph()
{
  GEO::CUT::COLOREDGRAPH::Graph g( 7 );

  g.Add( 0, 7 );
  g.Add( 0, 9 );
  g.Add( 0, 12 );
  g.Add( 0, 15 );

  g.Add( 1, 11 );
  g.Add( 1, 12 );
  g.Add( 1, 15 );

  g.Add( 2, 10 );
  g.Add( 2, 13 );
  g.Add( 2, 14 );
  g.Add( 2, 16 );

  g.Add( 3, 11 );
  g.Add( 3, 13 );
  g.Add( 3, 16 );

  g.Add( 4, 7 );
  g.Add( 4, 8 );
  g.Add( 4, 10 );

  g.Add( 5, 15 );
  g.Add( 5, 16 );

  g.Add( 6, 8 );
  g.Add( 6, 9 );
  g.Add( 6, 14 );

  g.FixSingleLines();
  g.Print();
  g.TestClosed();

  GEO::CUT::COLOREDGRAPH::Graph c( 7 );

  c.Add( 0, 7 );
  c.Add( 0, 9 );
  c.Add( 0, 12 );
  c.Add( 0, 15 );

  c.Add( 1, 11 );
  c.Add( 1, 12 );
  c.Add( 1, 15 );

  c.Add( 2, 10 );
  c.Add( 2, 13 );
  c.Add( 2, 14 );
  c.Add( 2, 16 );

  c.Add( 3, 11 );
  c.Add( 3, 13 );
  c.Add( 3, 16 );

  c.Add( 4, 7 );
  c.Add( 4, 8 );
  c.Add( 4, 10 );

  c.Add( 6, 8 );
  c.Add( 6, 9 );
  c.Add( 6, 14 );

  c.TestClosed();

  std::set<int> free;
  g.GetAll( free );

  for ( GEO::CUT::COLOREDGRAPH::Graph::const_iterator i=c.begin(); i!=c.end(); ++i )
  {
    int p = i->first;
    free.erase( p );
  }

  GEO::CUT::COLOREDGRAPH::Graph used( c );

  GEO::CUT::COLOREDGRAPH::CycleList cycle_list;
  cycle_list.AddPoints( g, used, c, free );

  for ( GEO::CUT::COLOREDGRAPH::CycleList::iterator i=cycle_list.begin(); i!=cycle_list.end(); ++i )
  {
    GEO::CUT::COLOREDGRAPH::Graph & g = *i;

    g.Print();
    g.TestClosed();
    g.TestSplit();
  }
}
