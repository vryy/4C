
#include <iostream>
#include <stdexcept>
#include <stack>
#include <queue>

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

namespace GEO
{
  namespace CUT
  {
    namespace COLOREDGRAPH
    {

      bool IsFree( Graph & used, plain_int_set & free, int i )
      {
        return used.count( i ) == 0 and free.count( i ) > 0;
      }

      int FindFirstFreeFacet( Graph & graph, Graph & used, plain_int_set & free )
      {
        for ( plain_int_set::iterator i=free.begin(); i!=free.end(); ++i )
        {
          int facet = *i;
          if ( facet >= graph.Split() )
          {
            throw std::runtime_error( "no free facet but free lines" );
          }
          plain_int_set & row = graph[facet];
          for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
          {
            int line = *i;
            if ( not IsFree( used, free, line ) )
            {
              return facet;
            }
          }
        }
        throw std::runtime_error( "empty free set" );
      }

      bool VisitFacetBFS( Graph & graph,
                          Graph & used,
                          plain_int_set & free,
                          int facet,
                          std::vector<int> & visited,
                          int & num_split_lines )
      {
        std::queue<int> facets;
        facets.push( facet );

        while ( not facets.empty() )
        {
          facet = facets.front();
          facets.pop();

          plain_int_set & row = graph[facet];

          // mark facet and lines
          visited[facet] += 1;
          for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
          {
            int line = *i;
            if ( not IsFree( used, free, line ) )
              num_split_lines += 1;
            visited[line] += 1;
          }

          // test for success: only non-free lines are open
          if ( num_split_lines==0 )
            return true;

          // try neighbouring facets
          for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
          {
            int line = *i;
            if ( visited[line] < 2 and IsFree( used, free, line ) )
            {
              plain_int_set & row = graph[line];
              for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
              {
                int f = *i;
                if ( visited[f] == 0 )
                {
                  facets.push( f );
                }
              }
            }
          }

          // unmark facet and lines
          for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
          {
            int line = *i;
            if ( not IsFree( used, free, line ) )
              num_split_lines -= 1;
            visited[line] -= 1;
          }
          visited[facet] -= 1;
        }

        return false;
      }

      bool IsValidFacet( plain_int_set & row, std::vector<int> & visited )
      {
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int line = *i;
          if ( visited[line] >= 2 )
          {
            return false;
          }
        }
        return true;
      }

      void MarkFacet( Graph & used,
                      plain_int_set & free,
                      int facet,
                      plain_int_set & row,
                      std::vector<int> & visited,
                      int & num_split_lines )
      {
        // mark facet and lines
        visited[facet] += 1;
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int line = *i;
          if ( IsFree( used, free, line ) )
          {
            if ( visited[line]==0 )
              num_split_lines += 1;
            else if ( visited[line]==1 )
              num_split_lines -= 1;
          }
          visited[line] += 1;
        }
      }

      void UnMarkFacet( Graph & used,
                        plain_int_set & free,
                        int facet,
                        plain_int_set & row,
                        std::vector<int> & visited,
                        int & num_split_lines )
      {
        // unmark facet and lines
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int line = *i;
          if ( IsFree( used, free, line ) )
          {
            if ( visited[line]==1 )
              num_split_lines -= 1;
            else if ( visited[line]==2 )
              num_split_lines += 1;
          }
          visited[line] -= 1;
        }
        visited[facet] -= 1;
      }

#if 1
      bool VisitFacetDFS( Graph & graph,
                          Graph & used,
                          plain_int_set & free,
                          int facet,
                          std::vector<int> & visited,
                          int & num_split_lines )
      {
        std::vector<int> facet_stack;
        facet_stack.push_back( facet );

        std::vector<int> facet_color( graph.size(), 0 );

        while ( not facet_stack.empty() )
        {
          facet = facet_stack.back();
          facet_stack.pop_back();

          if ( facet_color[facet] == 0 ) // white
          {
            plain_int_set & row = graph[facet];

            if ( IsValidFacet( row, visited ) )
            {
              MarkFacet( used, free, facet, row, visited, num_split_lines );

              // test for success: only non-free lines are open
              if ( num_split_lines==0 )
                return true;

              facet_color[facet] = 1;
              facet_stack.push_back( facet );

              // try neighbouring facets
              for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
              {
                int line = *i;
                if ( visited[line] < 2 and IsFree( used, free, line ) )
                {
                  plain_int_set & row = graph[line];
                  for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
                  {
                    int f = *i;
                    if ( facet_color[f] == 0 )
                    {
                      facet_stack.push_back( f );
                    }
                  }
                }
              }
            }
          }
          else if ( facet_color[facet] == 1 ) // gray
          {
            // search for any facet within the stack that can be added
            int pos = facet_stack.size()-1;
            for ( ; pos >=0; --pos )
            {
              int f = facet_stack[pos];
              if ( facet_color[f] == 0 and IsValidFacet( graph[f], visited ) )
              {
                facet_stack.push_back( facet );
                facet_stack.push_back( f );
                break;
              }
            }

            // clear current facet if there is nothing more to add
            if ( pos < 0 )
            {
              plain_int_set & row = graph[facet];
              UnMarkFacet( used, free, facet, row, visited, num_split_lines );

              if ( std::find( facet_stack.begin(), facet_stack.end(), facet )==facet_stack.end() )
                facet_color[facet] = 0;
              else
                facet_color[facet] = 2;
            }
          }
          else if ( facet_color[facet] == 2 ) // black
          {
            // if a black facet is poped for the last time (it is not any more
            // on the stack), we make it available again.
            if ( std::find( facet_stack.begin(), facet_stack.end(), facet )==facet_stack.end() )
            {
              facet_color[facet] = 0;
            }
          }
        }
        return false;
      }
#else
      bool VisitFacetDFS( Graph & graph,
                          Graph & used,
                          plain_int_set & free,
                          int facet,
                          std::vector<int> & visited,
                          int & num_split_lines )
      {
        plain_int_set & row = graph[facet];

        // mark facet and lines
        visited[facet] += 1;
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int line = *i;
          if ( not IsFree( used, free, line ) )
            num_split_lines += 1;
          visited[line] += 1;
        }

        // test for success: only non-free lines are open
#if 1
        if ( num_split_lines==0 )
          return true;
#else
        bool success = true;
        for ( std::vector<int>::iterator i=visited.begin() + graph.Split();
              i!=visited.end();
              ++i )
        {
          if ( *i == 1 )
          {
            int line = i - visited.begin();
            if ( IsFree( used, free, line ) )
            {
              success = false;
              break;
            }
          }
          else if ( *i > 2 )
          {
            throw std::runtime_error( "illegal facet combination" );
          }
        }
        if ( success )
        {
          return true;
        }
#endif

        // try neighbouring facets
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int line = *i;
          if ( visited[line] < 2 and IsFree( used, free, line ) )
          {
            plain_int_set & row = graph[line];
            for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
            {
              int f = *i;
              if ( visited[f] == 0 )
              {
                if ( VisitFacetDFS( graph, used, free, f, visited ) )
                  return true;
              }
            }
          }
        }

        // unmark facet and lines
        for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int line = *i;
          if ( not IsFree( used, free, line ) )
            num_split_lines -= 1;
          visited[line] -= 1;
        }
        visited[facet] -= 1;

        return false;
      }
#endif

    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::FindFreeFacets( Graph & graph,
                                                    Graph & used,
                                                    plain_int_set & free,
                                                    std::vector<int> & split_trace )
{
  int free_facet = FindFirstFreeFacet( graph, used, free );

  std::vector<int> visited( graph.size(), 0 );

  int num_split_lines = 0;
  if ( not VisitFacetDFS( graph, used, free, free_facet, visited, num_split_lines ) )
  {
    throw std::runtime_error( "failed to find volume split" );
  }

  for ( std::vector<int>::iterator i=visited.begin() + graph.Split();
        i!=visited.end();
        ++i )
  {
    if ( *i == 1 )
    {
      int line = i - visited.begin();
      if ( not IsFree( used, free, line ) )
      {
        split_trace.push_back( line );
      }
    }
  }
  if ( split_trace.size()==0 )
  {
    throw std::runtime_error( "no split trace" );
  }
  for ( std::vector<int>::iterator i=visited.begin();
        i!=visited.begin() + graph.Split();
        ++i )
  {
    if ( *i == 1 )
    {
      int facet = i - visited.begin();
      plain_int_set & row = graph[facet];
      for ( plain_int_set::iterator i=row.begin(); i!=row.end(); ++i )
      {
        int line = *i;
        used.Add( line, facet );
        Add( line, facet );
        free.erase( line );
      }
      free.erase( facet );
    }
    else if ( *i > 1 )
    {
      throw std::runtime_error( "same facet twice" );
    }
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
    std::vector<int> split_trace;
    connection.FindFreeFacets( graph, used, free, split_trace );

#if 0
#ifdef DEBUGCUTLIBRARY
    std::cout << "connection:\n";
    connection.Print();
#endif
#endif

    //connection.FindSplitTrace( split_trace );

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
