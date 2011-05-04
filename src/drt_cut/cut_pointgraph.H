#ifndef CUT_POINTGRAPH_H
#define CUT_POINTGRAPH_H

#include <map>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

namespace GEO
{
namespace CUT
{
  class Point;
  class Line;
  class Side;

class PointGraph
{

  struct Graph;
  struct FacetCycleList;
  struct FacetCycleListIterator;

  struct ForkFinder
  {
    ForkFinder( Graph & graph, Graph & used, const std::vector<Point*> & cycle, const std::set<Point*> & free )
      : graph_( graph ),
        used_( used ),
        cycle_( cycle ),
        free_( free )
    {
    }

    bool operator()( Point * point )
    {
      std::set<Point*> & row = graph_[point];
      if ( row.size() > 2 )
      {
        std::set<Point*> & used_row = used_[point];
        for ( std::set<Point*>::iterator i=row.begin(); i!=row.end(); ++i )
        {
          Point * p = *i;
          if ( used_row.count( p ) == 0 )
          {
            if ( free_.count( p ) > 0 )
              return true;
            if ( std::find( cycle_.begin(), cycle_.end(), p )!=cycle_.end() )
              return true;
          }
        }
        return false;
      }
      return false;
    }

    Graph & graph_;
    Graph & used_;
    const std::vector<Point*> & cycle_;
    const std::set<Point*> & free_;
  };

  struct FacetCycle
  {
    FacetCycle()
      : active_( true )
    {
    }

    void Assign( const std::vector<Point*> & cycle )
    {
      cycle_.reserve( cycle.size() );
      cycle_.assign( cycle.begin(), cycle.end() );
    }

    void Split( Graph & graph, Graph & used, FacetCycleList & facet_cycles, std::set<Point*> & free );

    void Print();

    bool active_;

    std::vector<Point*> cycle_;
  };

  struct Graph
  {
    void Add( Point * p1, Point * p2 );

    void AddAll( Side * side, Point * p1, Point * p2 );

    void AddCycle( const std::vector<Point*> & cycle );

    Point * FindNext( Graph & used, Point * p, const std::vector<Point*> & cycle, const std::set<Point*> & free );

    void IsValid();

    void Print();

    std::set<Point*> & operator[]( Point * p ) { return graph_[p]; }

    unsigned count( Point * p ) { return graph_.count( p ); }

    std::map<Point*, std::set<Point*> > graph_;
  };

  struct FacetCycleListIterator
  {
    FacetCycleListIterator( std::list<FacetCycle> facet_cycles, std::list<FacetCycle>::iterator i )
      : facet_cycles_( facet_cycles ),
        i_( i )
    {
      NextActive();
    }

    void NextActive()
    {
      while ( i_ != facet_cycles_.end() )
      {
        if ( i_->active_ )
          break;
        ++i_;
      }
    }

    FacetCycleListIterator & operator++()
    {
      ++i_;
      NextActive();
      return *this;
    }

    const std::vector<Point*> & operator*()
    {
      return i_->cycle_;
    }

    bool operator!=( const FacetCycleListIterator & other )
    {
      return i_ != other.i_;
    }

    std::list<FacetCycle> & facet_cycles_;
    std::list<FacetCycle>::iterator i_;
  };

  struct FacetCycleList
  {
    void AddFacetPoints( Graph & graph, Graph & used, std::vector<Point*> & cycle, std::set<Point*> & free );

    void Print();

    FacetCycleListIterator begin() { return FacetCycleListIterator( facet_cycles_, facet_cycles_.begin() ); }

    FacetCycleListIterator end()   { return FacetCycleListIterator( facet_cycles_, facet_cycles_.end() ); }

    std::list<FacetCycle> facet_cycles_;
  };

public:

  typedef FacetCycleListIterator iterator;

  PointGraph( Side * side );

  iterator begin() { return facet_cycles_.begin(); }

  iterator end()   { return facet_cycles_.end(); }

  void Print()
  {
    graph_.Print();
    facet_cycles_.Print();
  }

private:

  void AddFacetPoints( std::vector<Point*> & cycle, std::set<Point*> & free );

  Side * side_;
  Graph graph_;
  std::vector<std::vector<Point*> > facet_points_;

  FacetCycleList facet_cycles_;
};

}
}


#endif
