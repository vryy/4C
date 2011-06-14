#ifndef CUT_COLOREDGRAPH_H
#define CUT_COLOREDGRAPH_H

#include <map>
#include <set>
#include <vector>
#include <list>
#include <algorithm>

#include "cut_utils.H"

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
        ForkFinder( Graph & graph, Graph & used, Graph & cycle, const plain_int_set & free )
          : graph_( graph ),
            used_( used ),
            cycle_( cycle ),
            free_( free )
        {
        }

        bool operator()( const std::pair<const int, plain_int_set > & point );

        Graph & graph_;
        Graph & used_;
        Graph & cycle_;
        const plain_int_set & free_;
      };

      class SingeLineFinder
      {
      public:

        SingeLineFinder( int color_split )
          : color_split_( color_split )
        {}

        bool operator()( const std::pair<const int, plain_int_set > & i )
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

        typedef std::map<int, plain_int_set >::const_iterator const_iterator;

        explicit Graph( int color_split )
          : color_split_( color_split )
        {}

        void Add( int row, int col );

        void Add( int p, const plain_int_set & row );

        int FindNext( Graph & used, int p, Graph & cycle, const plain_int_set & free );

        void FindFreeFacets( Graph & graph, Graph & used, plain_int_set & free );

        void FindSplitTrace( std::vector<int> & split_trace );

        void GetAll( plain_int_set & all );

        void FixSingleLines();

        void TestClosed();

        void TestSplit();

        void TestFacets();

        void Print();

        plain_int_set & at( int p ) { return graph_[p]; }

        plain_int_set & operator[]( int p ) { return graph_[p]; }

        unsigned count( int p ) { return graph_.count( p ); }

        std::map<int, plain_int_set >::const_iterator begin() const { return graph_.begin(); }

        std::map<int, plain_int_set >::const_iterator end()   const { return graph_.end(); }

        void Swap( Graph & other )
        {
          std::swap( graph_, other.graph_ );
          std::swap( color_split_, other.color_split_ );
        }

        int Split() { return color_split_; }

        void SetSplit( int color_split ) { color_split_ = color_split; }

//         void Split( int color_split ) { color_split_ = color_split; }

        void Split( Graph & connection, const std::vector<int> & split_trace, Graph & c1, Graph & c2 );

        bool ContainsTrace( const std::vector<int> & split_trace );

      private:

        void Fill( const std::vector<int> & split_trace, Graph & connection, int seed, Graph & c );

        std::map<int, plain_int_set > graph_;
        int color_split_;
      };

      class Cycle
      {
        friend class CycleListIterator;
      public:

        explicit Cycle( int color_split )
          : cycle_( color_split )
        {
        }

        void Assign( Graph & cycle )
        {
          cycle_.Swap( cycle );
        }

        void Print();

        void Split( Graph & connection, const std::vector<int> & split_trace, Graph & c1, Graph & c2 )
        { cycle_.Split( connection, split_trace, c1, c2 ); }

        bool ContainsTrace( const std::vector<int> & split_trace )
        { return cycle_.ContainsTrace( split_trace ); }

      private:

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

        void AddPoints( Graph & graph, Graph & used, Graph & cycle, plain_int_set & free );

        void Print();

        unsigned size() const { return cycles_.size(); }

        CycleListIterator begin() { return CycleListIterator( cycles_, cycles_.begin() ); }

        CycleListIterator end()   { return CycleListIterator( cycles_, cycles_.end() ); }

      private:

        void PushBack( Graph & g );

        std::list<Cycle> cycles_;
      };

    }
  }
}

#endif
