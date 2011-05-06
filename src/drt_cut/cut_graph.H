#ifndef CUT_GRAPH_H
#define CUT_GRAPH_H

#include <map>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

namespace GEO
{
  namespace CUT
  {
    namespace GRAPH
    {

      class Graph;
      class CycleList;
      class CycleListIterator;

      class ForkFinder
      {
      public:
        ForkFinder( Graph & graph, Graph & used, const std::vector<int> & cycle, const std::set<int> & free )
          : graph_( graph ),
            used_( used ),
            cycle_( cycle ),
            free_( free )
        {
        }

        bool operator()( int point );

      private:

        bool BackToCycle( std::vector<int> & history, int p );

        Graph & graph_;
        Graph & used_;
        const std::vector<int> & cycle_;
        const std::set<int> & free_;
      };

      class Cycle
      {
        friend class CycleListIterator;
      public:

        Cycle()
          : active_( true )
        {
        }

        void Assign( const std::vector<int> & cycle )
        {
          cycle_.reserve( cycle.size() );
          cycle_.assign( cycle.begin(), cycle.end() );
        }

        void Split( Graph & graph, Graph & used, CycleList & cycles, std::set<int> & free );

        void Print();

        bool Active() const { return active_; }

      private:

        bool active_;

        std::vector<int> cycle_;
      };

      class Graph
      {
      public:

        typedef std::map<int, std::set<int> >::const_iterator const_iterator;

        void Add( int row, int col );

        void Erase( int p );

        void AddCycle( const std::vector<int> & cycle );

        int FindNext( Graph & used, int p, const std::vector<int> & cycle, const std::set<int> & free );

        void GetAll( std::set<int> & all );

        void FixSinglePoints();

        void TestClosed();

        void Print();

        std::set<int> & at( int p ) { return graph_[p]; }

        std::set<int> & operator[]( int p ) { return graph_[p]; }

        unsigned count( int p ) { return graph_.count( p ); }

        std::map<int, std::set<int> >::const_iterator begin() const { return graph_.begin(); }

        std::map<int, std::set<int> >::const_iterator end()   const { return graph_.end(); }

        unsigned size() const { return graph_.size(); }

      private:

        std::map<int, std::set<int> > graph_;
      };

      class CycleListIterator
      {
      public:

        CycleListIterator( std::list<Cycle> & cycles, std::list<Cycle>::iterator i )
          : cycles_( cycles ),
            i_( i )
        {
          NextActive();
        }

        void NextActive()
        {
          while ( i_ != cycles_.end() )
          {
            if ( i_->Active() )
              break;
            ++i_;
          }
        }

        CycleListIterator & operator++()
        {
          ++i_;
          NextActive();
          return *this;
        }

        const std::vector<int> & operator*()
        {
          return i_->cycle_;
        }

        bool operator!=( const CycleListIterator & other )
        {
          return i_ != other.i_;
        }

      private:
        std::list<Cycle> & cycles_;
        std::list<Cycle>::iterator i_;
      };

      class CycleList
      {
      public:

        void AddPoints( Graph & graph, Graph & used, std::vector<int> & cycle, std::set<int> & free );

        void AddFreePoints( Graph & graph, Graph & used, std::set<int> & free );

        void Print();

        CycleListIterator begin() { return CycleListIterator( cycles_, cycles_.begin() ); }

        CycleListIterator end()   { return CycleListIterator( cycles_, cycles_.end() ); }

        unsigned size() { return cycles_.size(); }

        unsigned ActiveCount();

      private:
        std::list<Cycle> cycles_;
      };

    }
  }
}

#endif
