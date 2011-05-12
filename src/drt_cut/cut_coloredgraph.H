#ifndef CUT_COLOREDGRAPH_H
#define CUT_COLOREDGRAPH_H

#include <map>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

namespace GEO
{
  namespace CUT
  {
    namespace COLOREDGRAPH
    {

      class Graph;
      class CycleList;
      class CycleListIterator;

      class ForkFinder
      {
      public:
        ForkFinder( Graph & graph, Graph & used, Graph & cycle, const std::set<int> & free )
          : graph_( graph ),
            used_( used ),
            cycle_( cycle ),
            free_( free )
        {
        }

        bool operator()( const std::pair<const int, std::set<int> > & point );

        Graph & graph_;
        Graph & used_;
        Graph & cycle_;
        const std::set<int> & free_;
      };

      class SingeLineFinder
      {
      public:

        SingeLineFinder( int color_split )
          : color_split_( color_split )
        {}

        bool operator()( const std::pair<const int, std::set<int> > & i )
        {
          int p = i.first;
          if ( p < color_split_ )
          {
            return i.second.size() < 3;
          }
          else
          {
            return i.second.size() < 2;
          }
        }

      private:
        int color_split_;
      };

      class Graph
      {
      public:

        typedef std::map<int, std::set<int> >::const_iterator const_iterator;

        explicit Graph( int color_split )
          : color_split_( color_split )
        {}

        void Add( int row, int col );

        void Add( int p, const std::set<int> & row );

        int FindNext( Graph & used, int p, Graph & cycle, const std::set<int> & free );

        void GetAll( std::set<int> & all );

        void FixSingleLines();

        void TestClosed();

        void TestSplit();

        void TestFacets();

        void Print();

        std::set<int> & at( int p ) { return graph_[p]; }

        std::set<int> & operator[]( int p ) { return graph_[p]; }

        unsigned count( int p ) { return graph_.count( p ); }

        std::map<int, std::set<int> >::const_iterator begin() const { return graph_.begin(); }

        std::map<int, std::set<int> >::const_iterator end()   const { return graph_.end(); }

        void Swap( Graph & other )
        {
          std::swap( graph_, other.graph_ );
          std::swap( color_split_, other.color_split_ );
        }

        int Split() { return color_split_; }

        void SetSplit( int color_split ) { color_split_ = color_split; }

//         void Split( int color_split ) { color_split_ = color_split; }

        void Split( Graph & connection, Graph & c1, Graph & c2 );

      private:

        void Fill( const std::set<int> & split_trace, Graph & connection, int seed, Graph & c );

        std::map<int, std::set<int> > graph_;
        int color_split_;
      };

      class Cycle
      {
        friend class CycleListIterator;
      public:

        explicit Cycle( int color_split )
          : active_( true ),
            cycle_( color_split )
        {
        }

        void Assign( Graph & cycle )
        {
          cycle_.Swap( cycle );
        }

        void Split( Graph & graph, Graph & used, CycleList & cycles, std::set<int> & free );

        void Print();

        bool Active() const { return active_; }

      private:

        bool active_;

        Graph cycle_;
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

        Graph & operator*()
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
        typedef CycleListIterator iterator;

        void AddPoints( Graph & graph, Graph & used, Graph & cycle, std::set<int> & free );

        void Print();

        unsigned size() const { return cycles_.size(); }

        CycleListIterator begin() { return CycleListIterator( cycles_, cycles_.begin() ); }

        CycleListIterator end()   { return CycleListIterator( cycles_, cycles_.end() ); }

      private:
        std::list<Cycle> cycles_;
      };

    }
  }
}

#endif
