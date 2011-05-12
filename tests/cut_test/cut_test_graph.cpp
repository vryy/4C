
#include <set>

#include "../../src/drt_cut/cut_coloredgraph.H"
#include "cut_test_utils.H"

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
