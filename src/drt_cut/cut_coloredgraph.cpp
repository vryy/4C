
#include <iostream>
#include <stdexcept>
#include <stack>

#include "cut_coloredgraph.H"

bool GEO::CUT::COLOREDGRAPH::ForkFinder::operator()( const std::pair<const int, plain_int_set > & point )
{
  if ( point.first < graph_.Split() )
    return false;

  plain_int_set & row = graph_[point.first];
  if ( row.size() > 2 )
  {
    plain_int_set & used_row = used_[point.first];
    for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
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

void GEO::CUT::COLOREDGRAPH::Graph::Add( int p, const plain_int_set & row )
{
  for ( plain_int_set::const_iterator i=row.begin(); i!=row.end(); ++i )
  {
    Add( p, *i );
  }
}

int GEO::CUT::COLOREDGRAPH::Graph::FindNext( Graph & used, int point, Graph & cycle, const plain_int_set & free )
{
  plain_int_set & row      = graph_[point];
  plain_int_set & used_row = used[point];
  for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
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

void GEO::CUT::COLOREDGRAPH::Graph::GetAll( plain_int_set & all )
{
  for ( std::map<int, plain_int_set >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    all.insert( p );
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::FixSingleLines()
{
  for ( ;; )
  {
    std::map<int, plain_int_set >::iterator j = std::find_if( graph_.begin(), graph_.end(), SingeLineFinder( color_split_ ) );
    if ( j==graph_.end() )
    {
      return;
    }

    int p1 = j->first;
    plain_int_set & row = j->second;
    for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p2 = *i;
      plain_int_set & row2 = graph_[p2];
      row2.erase( p1 );
      if ( row2.size()==0 )
        graph_.erase( p2 );
    }
    graph_.erase( j );
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestClosed()
{
  for ( std::map<int, plain_int_set >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    plain_int_set & row = i->second;
    if ( row.size() < 2 )
    {
      throw std::runtime_error( "open point in colored graph" );
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestSplit()
{
  for ( std::map<int, plain_int_set >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    if ( p >= color_split_ )
    {
      plain_int_set & row = i->second;
      if ( row.size() > 2 )
      {
        throw std::runtime_error( "colored graph not properly split" );
      }
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestFacets()
{
  for ( std::map<int, plain_int_set >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    if ( p < color_split_ )
    {
      plain_int_set & row = i->second;
      if ( row.size() < 3 )
      {
        throw std::runtime_error( "facets need at least three lines" );
      }
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::Print()
{
  for ( std::map<int, plain_int_set >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    plain_int_set & row = i->second;
    std::cout << p << ": ";
    for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      std::cout << p << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void GEO::CUT::COLOREDGRAPH::Graph::FindFreeFacets( Graph & graph, Graph & used, plain_int_set & free )
{
  int free_facet = *free.begin();
  if ( free_facet >= Split() )
  {
    throw std::runtime_error( "no free facet but free lines" );
  }

  std::stack<int> stack;
  stack.push( free_facet );

  while ( not stack.empty() )
  {
    int facet = stack.top();
    stack.pop();

    plain_int_set & row = graph[facet];
    for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int line = *i;

      if ( used.count( line ) == 0 and free.count( line ) > 0 )
      {
        plain_int_set & row = graph[line];
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int f = *i;
          if ( f != facet and used.count( f ) == 0 and free.count( f ) > 0 )
          {
            stack.push( f );

            // We can never have more than one new facet per line. This will
            // not work in very extreme cases, when many internal lines have
            // multiple choices.
            break;
          }
        }
      }

      used.Add( line, facet );
      Add( line, facet );
    }
  }

  for ( Graph::const_iterator i=begin(); i!=end(); ++i )
  {
    int p = i->first;

    const plain_int_set & row = i->second;
    for ( plain_int_set::const_iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      free.erase( p );
    }
    free.erase( p );
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::FindSplitTrace( std::vector<int> & split_trace )
{
  for ( Graph::const_iterator i=begin(); i!=end(); ++i )
  {
    int p = i->first;
    const plain_int_set & row = i->second;
    if ( p >= Split() )
    {
      if ( row.size()==1 )
      {
        split_trace.push_back( p );
      }
    }
  }
}

bool GEO::CUT::COLOREDGRAPH::Graph::ContainsTrace( const std::vector<int> & split_trace )
{
  for ( std::vector<int>::const_iterator i=split_trace.begin(); i!=split_trace.end(); ++i )
  {
    int p = *i;
    if ( count( p )==0 )
      return false;
  }
  return true;
}

void GEO::CUT::COLOREDGRAPH::Cycle::Print()
{
  std::cout << "Cycle:\n";
  cycle_.Print();
  std::cout << "\n";
}

void GEO::CUT::COLOREDGRAPH::Graph::Split( Graph & connection, const std::vector<int> & split_trace, Graph & c1, Graph & c2 )
{
  // find lhs and rhs starting from split trace

  int p = split_trace.front();
  plain_int_set & row = at( p );

  if ( row.size()!=2 )
    throw std::runtime_error( "expect two facets at line" );

  plain_int_set::iterator i = row.begin();
  Fill( split_trace, connection, *i, c1 );
  ++i;
  Fill( split_trace, connection, *i, c2 );

  for ( std::vector<int>::const_iterator i=split_trace.begin(); i!=split_trace.end(); ++i )
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
      plain_int_set & row = at( p );
      if ( row.size()!=2 )
        throw std::runtime_error( "expect two facets at line" );
      plain_int_set::iterator i = row.begin();
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

void GEO::CUT::COLOREDGRAPH::Graph::Fill( const std::vector<int> & split_trace, Graph & connection, int seed, Graph & c )
{
  plain_int_set done;
  done.insert( split_trace.begin(), split_trace.end() );
  std::stack<int> stack;

  stack.push( seed );

  while ( not stack.empty() )
  {
    int f = stack.top();
    stack.pop();

    done.insert( f );

    plain_int_set & row = at( f );
    for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      c.Add( f, p );
      if ( done.count( p )==0 )
      {
        plain_int_set & row = at( p );
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
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
    const plain_int_set & row = i->second;
    c.Add( p, row );
  }
}

void GEO::CUT::COLOREDGRAPH::CycleList::AddPoints( Graph & graph, Graph & used, Graph & cycle, plain_int_set & free )
{
  PushBack( cycle );

  while ( free.size() > 0 )
  {
    Graph connection( graph.Split() );
    connection.FindFreeFacets( graph, used, free );

    std::vector<int> split_trace;
    connection.FindSplitTrace( split_trace );

    bool found = false;
    for ( std::list<Cycle>::iterator i=cycles_.begin(); i!=cycles_.end(); ++i )
    {
      Cycle & c = *i;

      // There can be just one graph that contains the trace. This prohibits self-cuts.
      if ( c.ContainsTrace( split_trace ) )
      {
        Graph c1( graph.Split() );
        Graph c2( graph.Split() );

        c.Split( connection, split_trace, c1, c2 );

        cycles_.erase( i );

        PushBack( c1 );
        PushBack( c2 );

        found = true;
        break;
      }
    }
    if ( not found )
    {
      throw std::runtime_error( "did not find volume that contains split facets" );
    }
  }
}

void GEO::CUT::COLOREDGRAPH::CycleList::PushBack( Graph & g )
{
  cycles_.push_back( Cycle( g.Split() ) );
  Cycle & c = cycles_.back();
  c.Assign( g );
}

void GEO::CUT::COLOREDGRAPH::CycleList::Print()
{
  for ( std::list<Cycle>::iterator i=cycles_.begin(); i!=cycles_.end(); ++i )
  {
    Cycle & c = *i;
    c.Print();
  }
}
