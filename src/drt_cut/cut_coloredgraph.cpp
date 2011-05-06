
#include <iostream>
#include <stdexcept>
#include <stack>

#include "cut_coloredgraph.H"

bool GEO::CUT::COLOREDGRAPH::ForkFinder::operator()( const std::pair<const int, std::set<int> > & point )
{
  if ( point.first < graph_.Split() )
    return false;

  std::set<int> & row = graph_[point.first];
  if ( row.size() > 2 )
  {
    std::set<int> & used_row = used_[point.first];
    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      if ( used_row.count( p ) == 0 )
      {
        if ( free_.count( p ) > 0 )
          return true;
        if ( cycle_.count( p ) > 0 )
          return true;
      }
    }
    return false;
  }
  return false;
}

void GEO::CUT::COLOREDGRAPH::Graph::Add( int row, int col )
{
  if ( row >= color_split_ and col >= color_split_ )
    throw std::runtime_error( "two lines connected" );
  if ( row < color_split_ and col < color_split_ )
    throw std::runtime_error( "two facets connected" );
  graph_[row].insert( col );
  graph_[col].insert( row );
}

void GEO::CUT::COLOREDGRAPH::Graph::Add( int p, const std::set<int> & row )
{
  for ( std::set<int>::const_iterator i=row.begin(); i!=row.end(); ++i )
  {
    Add( p, *i );
  }
}

int GEO::CUT::COLOREDGRAPH::Graph::FindNext( Graph & used, int point, Graph & cycle, const std::set<int> & free )
{
  std::set<int> & row      = graph_[point];
  std::set<int> & used_row = used[point];
  for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
  {
    int p = *i;
    if ( used_row.count( p ) == 0 )
    {
      if ( free.count( p ) > 0 )
        return p;
      if ( cycle.count( p ) > 0 )
        return p;
    }
  }
  return -1;
}

void GEO::CUT::COLOREDGRAPH::Graph::GetAll( std::set<int> & all )
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    all.insert( p );
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::FixSingleLines()
{
  for ( ;; )
  {
    std::map<int, std::set<int> >::iterator j = std::find_if( graph_.begin(), graph_.end(), SingeLineFinder() );
    if ( j==graph_.end() )
    {
      return;
    }

    int p1 = j->first;
    std::set<int> & row = j->second;
    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p2 = *i;
      std::set<int> & row2 = graph_[p2];
      row2.erase( p1 );
      if ( row2.size()==0 )
        graph_.erase( p2 );
    }
    graph_.erase( j );
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestClosed()
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    std::set<int> & row = i->second;
    if ( row.size() < 2 )
    {
      throw std::runtime_error( "open point in colored graph" );
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::Print()
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    std::set<int> & row = i->second;
    std::cout << p << ": ";
    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      std::cout << p << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void GEO::CUT::COLOREDGRAPH::Cycle::Split( Graph & graph, Graph & used, CycleList & cycles, std::set<int> & free )
{
  Graph::const_iterator start = std::find_if( cycle_.begin(), cycle_.end(), ForkFinder( graph, used, cycle_, free ) );
  if ( start!=cycle_.end() )
  {
    Graph connection( graph.Split() );
    std::set<int> open;

    int p1 = start->first;
    while ( p1 > -1 )
    {
      int p2 = graph.FindNext( used, p1, cycle_, free );
      if ( p2 > -1 )
      {
        std::set<int> & row = graph[p2];
        for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int p3 = *i;

          if ( used.count( p3 ) == 0 )
            open.insert( p3 );
          else
            open.erase( p3 );

          used.Add( p3, p2 );
          connection.Add( p3, p2 );
        }

        if ( open.size() > 0 )
        {
          p1 = *open.begin();
        }
        else
        {
          // split graph:

          Graph c1( graph.Split() );
          Graph c2( graph.Split() );

          cycle_.Split( connection, c1, c2 );

          for ( Graph::const_iterator i=connection.begin(); i!=connection.end(); ++i )
          {
            int p = i->first;

            std::set<int> & row = graph[p];
            for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
            {
              int p = *i;
              free.erase( p );
            }
            free.erase( p );
          }

          cycles.AddPoints( graph, used, c1, free );
          cycles.AddPoints( graph, used, c2, free );

          active_ = false;
          p1 = -1;
        }
      }
      else
      {
        throw std::runtime_error( "failed to find next point" );
      }
    }
  }
  else
  {
#ifdef DEBUGCUTLIBRARY
    cycle_.TestClosed();
#endif
    active_ = true;
  }
}

void GEO::CUT::COLOREDGRAPH::Cycle::Print()
{
  std::cout << "Cycle: " << active_ << "\n";
  cycle_.Print();
  std::cout << "\n";
}

void GEO::CUT::COLOREDGRAPH::Graph::Split( Graph & connection, Graph & c1, Graph & c2 )
{
  // find split trace

  std::set<int> split_trace;
  for ( Graph::const_iterator i=connection.begin(); i!=connection.end(); ++i )
  {
    int p = i->first;
    const std::set<int> & row = i->second;
    if ( p >= connection.Split() )
    {
      if ( row.size()==1 )
      {
        split_trace.insert( p );
      }
    }
  }

  // find lhs and rhs starting from split trace

  int p = *split_trace.begin();
  std::set<int> & row = at( p );

  if ( row.size()!=2 )
    throw std::runtime_error( "expect two facets at line" );

  std::set<int>::iterator i = row.begin();
  Fill( split_trace, connection, *i, c1 );
  ++i;
  Fill( split_trace, connection, *i, c2 );

  for ( std::set<int>::iterator i=split_trace.begin(); i!=split_trace.end(); ++i )
  {
    int p = *i;
    unsigned c1rowlen = c1[p].size();
    unsigned c2rowlen = c2[p].size();

    Graph * fine = NULL;
    Graph * open = NULL;

    if ( c1rowlen>1 and c2rowlen>1 )
    {
      // fine
    }
    else if ( c1rowlen==1 and c2rowlen>1 )
    {
      fine = &c2;
      open = &c1;
    }
    else if ( c1rowlen>1 and c2rowlen==1 )
    {
      fine = &c1;
      open = &c2;
    }
    else
    {
      throw std::runtime_error( "open line after graph split" );
    }

    if ( open!=NULL )
    {
      std::set<int> & row = at( p );
      if ( row.size()!=2 )
        throw std::runtime_error( "expect two facets at line" );
      std::set<int>::iterator i = row.begin();
      int f1 = *i;
      ++i;
      int f2 = *i;
      if ( fine->count( f1 ) > 0 )
      {
        open->Add( f2, at( f2 ) );
      }
      else if ( fine->count( f2 ) > 0 )
      {
        open->Add( f1, at( f1 ) );
      }
      else
      {
        throw std::runtime_error( "confused" );
      }
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::Fill( const std::set<int> & split_trace, Graph & connection, int seed, Graph & c )
{
  std::set<int> done( split_trace );
  std::stack<int> stack;

  stack.push( seed );

  while ( not stack.empty() )
  {
    int f = stack.top();
    stack.pop();

    done.insert( f );

    std::set<int> & row = at( f );
    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      c.Add( f, p );
      if ( done.count( p )==0 )
      {
        std::set<int> & row = at( p );
        for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int f = *i;
          if ( done.count( f )==0 )
          {
            stack.push( f );
          }
        }
      }
    }
  }

  for ( Graph::const_iterator i=connection.begin(); i!=connection.end(); ++i )
  {
    int p = i->first;
    const std::set<int> & row = i->second;
    c.Add( p, row );
  }
}

void GEO::CUT::COLOREDGRAPH::CycleList::AddPoints( Graph & graph, Graph & used, Graph & cycle, std::set<int> & free )
{
  cycles_.push_back( Cycle( graph.Split() ) );
  Cycle & c = cycles_.back();
  c.Assign( cycle );

  c.Split( graph, used, *this, free );
}

void GEO::CUT::COLOREDGRAPH::CycleList::Print()
{
  for ( std::list<Cycle>::iterator i=cycles_.begin(); i!=cycles_.end(); ++i )
  {
    Cycle & c = *i;
    c.Print();
  }
}
