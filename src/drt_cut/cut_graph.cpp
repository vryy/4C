
#include <iostream>
#include <stdexcept>
#include "cut_graph.H"

bool GEO::CUT::GRAPH::ForkFinder::operator()( int point )
{
  std::set<int> & row = graph_[point];
  if ( row.size() > 2 )
  {
    std::set<int> & used_row = used_[point];
    for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
    {
      int p = *i;
      if ( used_row.count( p ) == 0 )
      {
        if ( free_.count( p ) > 0 )
        {
          std::vector<int> history( 1, point );
          if ( BackToCycle( history, p ) )
            return true;
        }
        else if ( std::find( cycle_.begin(), cycle_.end(), p )!=cycle_.end() )
          return true;
      }
    }
    return false;
  }
  return false;
}

bool GEO::CUT::GRAPH::ForkFinder::BackToCycle( std::vector<int> & history, int p )
{
  std::set<int> & row = graph_[p];
  std::set<int> & used_row = used_[p];
  for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
  {
    int p1 = *i;
    if ( used_row.count( p1 ) == 0 )
    {
      std::vector<int>::iterator h = std::find( history.begin(), history.end(), p1 );
      if ( h == history.end() )
      {
        if ( free_.count( p1 ) > 0 )
        {
          history.push_back( p );
          if ( BackToCycle( history, p1 ) )
          {
            return true;
          }
          history.resize( history.size()-1 );
        }
        else if ( std::find( cycle_.begin(), cycle_.end(), p1 )!=cycle_.end() )
          return true;
      }
      else if ( h!=history.end()-1 )
      {
        return true;
      }
    }
  }
  return false;
}

void GEO::CUT::GRAPH::Graph::Add( int row, int col )
{
  graph_[row].insert( col );
  graph_[col].insert( row );
}

void GEO::CUT::GRAPH::Graph::Erase( int p )
{
  std::set<int> & row = graph_[p];
  for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
  {
    int p2 = *i;
    std::set<int> & row2 = graph_[p2];
    row2.erase( p );
    if ( row2.size()==0 )
    {
      graph_.erase( p2 );
    }
  }
  graph_.erase( p );
}

void GEO::CUT::GRAPH::Graph::AddCycle( const std::vector<int> & cycle )
{
  unsigned size = cycle.size();
  for ( unsigned i=0; i<size; ++i )
  {
    int p1 = cycle[i];
    int p2 = cycle[( i+1 ) % size];

    Add( p1, p2 );
  }
}

int GEO::CUT::GRAPH::Graph::FindNext( Graph & used, int point, const std::vector<int> & cycle, const std::set<int> & free )
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
      if ( std::find( cycle.begin(), cycle.end(), p )!=cycle.end() )
        return p;
    }
  }
  return -1;
}

void GEO::CUT::GRAPH::Graph::GetAll( std::set<int> & all )
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    int p = i->first;
    all.insert( p );
  }
}

void GEO::CUT::GRAPH::Graph::FixSinglePoints()
{
  for ( ;; )
  {
    bool found = false;
    for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
    {
      int p = i->first;
      std::set<int> & row = i->second;
      if ( row.size() < 2 )
      {
        found = true;
        for ( std::set<int>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          int p2 = *i;
          std::set<int> & row2 = graph_[p2];
          row2.erase( p );
          if ( row2.size()==0 )
            graph_.erase( p2 );
        }
        graph_.erase( p );
        break;
      }
    }
    if ( not found )
    {
      return;
    }
  }
}

void GEO::CUT::GRAPH::Graph::TestClosed()
{
  for ( std::map<int, std::set<int> >::iterator i=graph_.begin(); i!=graph_.end(); ++i )
  {
    std::set<int> & row = i->second;
    if ( row.size() < 2 )
    {
      throw std::runtime_error( "open point in graph" );
    }
  }
}

void GEO::CUT::GRAPH::Graph::Print()
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

void GEO::CUT::GRAPH::Cycle::Split( Graph & graph, Graph & used, CycleList & cycles, std::set<int> & free )
{
  std::vector<int>::iterator start = std::find_if( cycle_.begin(), cycle_.end(), ForkFinder( graph, used, cycle_, free ) );
  if ( start!=cycle_.end() )
  {
    std::vector<int> connection;

    int p1 = *start;
    while ( p1 > -1 )
    {
      int p2 = graph.FindNext( used, p1, cycle_, free );
      if ( p2 > -1 )
      {
        used.Add( p1, p2 );

        std::vector<int>::iterator end = std::find( cycle_.begin(), cycle_.end(), p2 );
        if ( end!=cycle_.end() )
        {
          // split here.

          if ( end < start )
            throw std::runtime_error( "iterator confusion" );

          std::vector<int> c1;
          std::vector<int> c2;

          std::copy( cycle_.begin(), start, std::back_inserter( c1 ) );
          c1.push_back( *start );
          std::copy( connection.begin(), connection.end(), std::back_inserter( c1 ) );
          std::copy( end, cycle_.end(), std::back_inserter( c1 ) );

          std::copy( start, end, std::back_inserter( c2 ) );
          c2.push_back( *end );
          std::copy( connection.rbegin(), connection.rend(), std::back_inserter( c2 ) );

          for ( std::vector<int>::iterator i=connection.begin(); i!=connection.end(); ++i )
          {
            int p = *i;
            free.erase( p );
          }

          cycles.AddPoints( graph, used, c1, free );
          cycles.AddPoints( graph, used, c2, free );

          active_ = false;

          p1 = -1;
        }
        else if ( ( end = std::find( connection.begin(), connection.end(), p2 ) )!=connection.end() )
        {
          // loop

          std::vector<int> c1;
          std::vector<int> c2;

          std::copy( cycle_.begin(), start, std::back_inserter( c1 ) );
          c1.push_back( *start );
          std::copy( connection.begin(), connection.end(), std::back_inserter( c1 ) );
          for ( std::vector<int>::iterator i=end; i!=connection.begin(); --i )
          {
            c1.push_back( *i );
          }
          c1.push_back( *connection.begin() );
          std::copy( start, cycle_.end(), std::back_inserter( c1 ) );

          std::copy( end, connection.end(), std::back_inserter( c2 ) );

          for ( std::vector<int>::iterator i=connection.begin(); i!=connection.end(); ++i )
          {
            int p = *i;
            free.erase( p );
          }

          cycles.AddPoints( graph, used, c1, free );
          cycles.AddPoints( graph, used, c2, free );

          active_ = false;

          p1 = -1;
        }
        else
        {
          connection.push_back( p2 );
          p1 = p2;
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
    active_ = true;
  }
}

void GEO::CUT::GRAPH::Cycle::Print()
{
  std::cout << active_;
  for ( std::vector<int>::iterator i=cycle_.begin(); i!=cycle_.end(); ++i )
  {
    int p = *i;
    std::cout << " " << p;
  }
  std::cout << "\n";
}

void GEO::CUT::GRAPH::CycleList::AddPoints( Graph & graph, Graph & used, std::vector<int> & cycle, std::set<int> & free )
{
  cycles_.push_back( Cycle() );
  Cycle & c = cycles_.back();
  c.Assign( cycle );

  c.Split( graph, used, *this, free );
}

void GEO::CUT::GRAPH::CycleList::AddFreePoints( Graph & graph, Graph & used, std::set<int> & free )
{
  while ( free.size() > 0 )
  {
    std::vector<int> cycle;

    int p1 = *free.begin();
    free.erase( p1 );
    cycle.push_back( p1 );
    while ( p1 != -1 )
    {
      int p2 = graph.FindNext( used, p1, cycle, free );
      if ( p2 > -1 )
      {
        used.Add( p1, p2 );
        std::vector<int>::iterator end = std::find( cycle.begin(), cycle.end(), p2 );
        if ( end!=cycle.end() )
        {
          std::vector<int> c2;
          if ( end==cycle.begin() )
          {
            std::swap( cycle, c2 );
          }
          else
          {
            for ( std::vector<int>::iterator i=cycle.begin(); i!=end; ++i )
            {
              int p3 = *i;
              free.insert( p3 );
              used.Erase( p3 );
            }
            c2.reserve( cycle.end()-end );
            std::copy( end, cycle.end(), std::back_inserter( c2 ) );
          }

          if ( c2.size() < 3 )
            throw std::runtime_error( "degenerated inner cycle" );

          AddPoints( graph, used, c2, free );

          p1 = -1;
        }
        else
        {
          free.erase( p2 );
          cycle.push_back( p2 );
          p1 = p2;
        }
      }
      else
      {
        throw std::runtime_error( "failed to find next point in inner cycle" );
      }
    }
  }
}

void GEO::CUT::GRAPH::CycleList::Print()
{
  for ( std::list<Cycle>::iterator i=cycles_.begin(); i!=cycles_.end(); ++i )
  {
    Cycle & c = *i;
    c.Print();
  }
}

unsigned GEO::CUT::GRAPH::CycleList::ActiveCount()
{
  unsigned count = 0;
  for ( CycleListIterator i=begin(); i!=end(); ++i )
  {
    count += 1;
  }
  return count;
}

