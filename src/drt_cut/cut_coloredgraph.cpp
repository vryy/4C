/*---------------------------------------------------------------------*/
/*!
\file cut_coloredgraph.cpp

\brief colored graph to create volumecells from facets
\warning do we still use this graph?

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include <fstream>
#include <stack>
#include <queue>
#include <algorithm>
#include <stdexcept>

#include "cut_coloredgraph.H"

bool GEO::CUT::COLOREDGRAPH::ForkFinder::operator()(
    const std::pair<const int, plain_int_set>& point)
{
  if (point.first < graph_.Split()) return false;

  plain_int_set& row = graph_[point.first];
  if (row.size() > 2)
  {
    plain_int_set& used_row = used_[point.first];
    for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
    {
      int p = *i;
      if (used_row.count(p) == 0)
      {
        if (free_.count(p) > 0) return true;
        if (cycle_.count(p) > 0) return true;
      }
    }
    return false;
  }
  return false;
}

void GEO::CUT::COLOREDGRAPH::Graph::Add(int row, int col)
{
  if (row >= color_split_ and col >= color_split_) run_time_error("two lines connected");
  if (row < color_split_ and col < color_split_) run_time_error("two facets connected");
  graph_[row].insert(col);
  graph_[col].insert(row);
}

void GEO::CUT::COLOREDGRAPH::Graph::Add(int p, const plain_int_set& row)
{
  for (plain_int_set::const_iterator i = row.begin(); i != row.end(); ++i)
  {
    Add(p, *i);
  }
}

int GEO::CUT::COLOREDGRAPH::Graph::FindNext(
    Graph& used, int point, Graph& cycle, const plain_int_set& free)
{
  plain_int_set& row = graph_[point];
  plain_int_set& used_row = used[point];
  for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
  {
    int p = *i;
    if (used_row.count(p) == 0)
    {
      if (free.count(p) > 0) return p;
      if (cycle.count(p) > 0) return p;
    }
  }
  return -1;
}

void GEO::CUT::COLOREDGRAPH::Graph::GetAll(plain_int_set& all)
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int p = i->first;
    all.insert(p);
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::FixSingleLines()
{
  for (;;)
  {
    std::map<int, plain_int_set>::iterator j =
        std::find_if(graph_.begin(), graph_.end(), SingeLineFinder(color_split_));
    if (j == graph_.end())
    {
      return;
    }

    int p1 = j->first;
    plain_int_set& row = j->second;
    for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
    {
      int p2 = *i;
      plain_int_set& row2 = graph_[p2];
      row2.erase(p1);
      if (row2.size() == 0) graph_.erase(p2);
    }
    graph_.erase(j);
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestClosed()
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    plain_int_set& row = i->second;
    if (row.size() < 2)
    {
      // -------------------------------------------------------------------    (sudhakar 08/14)
      // one of the possible reasons this may occur is the following:
      // A background element is cut by two cut sides, that are not connected
      // This is not a multiple cut situation ==> this works with cut algorithm
      // See the 2D example here
      //
      //        ++++++++++++++++              ++++++++++++++
      //        +              +              +            +
      //        +              +              +            +
      //  o------------------------o          +            +
      //        +              +          o-------o   o-----------o
      //        +              +              +            +
      //  o------------------------o          +            +
      //        +              +              +            +
      //        ++++++++++++++++              ++++++++++++++
      //     (multiple cut --> okay)     (open-point in colored graph)
      //
      // In such situations, geometrically two separate volumecells can't be formed
      // Check your cut_mesh from cut_mesh*.pos
      // -------------------------------------------------------------------
      dserror("open point in colored graph ( facet-id = %d )", i->first);
      run_time_error("open point in colored graph");
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestFacets()
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int p = i->first;
    if (p < color_split_)
    {
      plain_int_set& row = i->second;
      if (row.size() < 3)
      {
        run_time_error("facets need at least three lines");
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::COLOREDGRAPH::Graph::Print() const
{
  for (std::map<int, plain_int_set>::const_iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int p = i->first;
    const plain_int_set& row = i->second;
    std::cout << p << ": ";
    for (plain_int_set::const_iterator j = row.begin(); j != row.end(); ++j)
    {
      int pp = *j;
      std::cout << pp << " ";
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
      bool IsFree(Graph& used, plain_int_set& free, int i)
      {
        return used.count(i) == 0 and free.count(i) > 0;
      }

      int FindFirstFreeFacet(Graph& graph, Graph& used, plain_int_set& free)
      {
        for (plain_int_set::iterator i = free.begin(); i != free.end(); ++i)
        {
          int facet = *i;
          if (facet >= graph.Split())
          {
            throw std::runtime_error("no free facet but free lines");
          }
          plain_int_set& row = graph[facet];
          for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
          {
            int line = *i;
            if (not IsFree(used, free, line))
            {
              return facet;
            }
          }
        }
        throw std::runtime_error("empty free set");
      }

      bool IsValidFacet(plain_int_set& row, const std::vector<int>& visited)
      {
        for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
        {
          int line = *i;
          if (visited[line] >= 2)
          {
            return false;
          }
        }
        return true;
      }

      void MarkFacet(Graph& used, plain_int_set& free, int facet, plain_int_set& row,
          std::vector<int>& visited, int& num_split_lines)
      {
        // mark facet and lines
        visited[facet] += 1;
        if (visited[facet] > 1) throw std::runtime_error("facet visited more than once");
        for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
        {
          int line = *i;
          if (IsFree(used, free, line))
          {
            if (visited[line] == 0)
              num_split_lines += 1;
            else if (visited[line] == 1)
              num_split_lines -= 1;
          }
          visited[line] += 1;
          if (visited[line] > 2) throw std::runtime_error("too many facets at line");
        }
      }

      void UnMarkFacet(Graph& used, plain_int_set& free, int facet, plain_int_set& row,
          std::vector<int>& visited, int& num_split_lines)
      {
        // unmark facet and lines
        for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
        {
          int line = *i;
          if (IsFree(used, free, line))
          {
            if (visited[line] == 1)
              num_split_lines -= 1;
            else if (visited[line] == 2)
              num_split_lines += 1;
          }
          visited[line] -= 1;
          if (visited[line] < 0) throw std::runtime_error("too few facets at line");
        }
        visited[facet] -= 1;
        if (visited[facet] < 0) run_time_error("facet left more than once");
      }

#if 0
      bool IsOneCycle( const std::vector<std::pair<Point*, Point*> > & all_lines,
                       std::vector<int> & split_trace,
                       int split_color )
      {
        // See if those split lines form one loop. Just one.

        std::map<Point*, std::vector<int> > points;
        for ( std::vector<int>::iterator i=split_trace.begin(); i!=split_trace.end(); ++i )
        {
          int line = *i;
          int pos = line - split_color;
          if ( pos < 0 or pos >= static_cast<int>( all_lines.size() ) )
          {
            run_time_error( "line index error" );
          }
          points[all_lines[pos].first ].push_back( pos );
          points[all_lines[pos].second].push_back( pos );
        }

        Point * first = points.begin()->first;
        int l = points[first][0];

        unsigned count = 0;
        for ( Point * next = ( all_lines[l].first == first ) ? all_lines[l].second : all_lines[l].first;
              next != first;
              next = ( all_lines[l].first == next ) ? all_lines[l].second : all_lines[l].first )
        {
          if ( points[next].size() != 2 )
          {
            //throw std::runtime_error( "not a cycle" );
            return false;
          }
          l = ( points[next][0] == l ) ? points[next][1] : points[next][0];
          count += 1;
          if ( count >= split_trace.size() )
          {
            run_time_error( "overflow." );
          }
        }

        return count == split_trace.size()-1;
      }
#endif

      bool VisitFacetDFS(Graph& graph, Graph& used, plain_int_set& free, int facet,
          const std::vector<std::pair<Point*, Point*>>& all_lines, std::vector<int>& visited,
          std::vector<int>& split_trace)
      {
        std::vector<int> facet_stack;
        facet_stack.push_back(facet);

        int num_split_lines = 0;
        std::vector<int> facet_color(graph.size(), 0);

        while (not facet_stack.empty())
        {
          facet = facet_stack.back();
          facet_stack.pop_back();

          if (facet_color[facet] == 0)  // white
          {
            plain_int_set& row = graph[facet];

            if (IsValidFacet(row, visited))
            {
              MarkFacet(used, free, facet, row, visited, num_split_lines);

              facet_color[facet] = 1;
              facet_stack.push_back(facet);

              // test for success: only non-free lines are open
              if (num_split_lines == 0)
              {
                // build split_trace
                for (std::vector<int>::iterator i = visited.begin() + graph.Split();
                     i != visited.end(); ++i)
                {
                  if (*i == 1)
                  {
                    int line = i - visited.begin();
                    if (not IsFree(used, free, line))
                    {
                      split_trace.push_back(line);
                    }
                  }
                }
                if (split_trace.size() == 0)
                {
                  run_time_error("no split trace");
                }

                // ok. all points connected.
                // if ( IsOneCycle( all_lines, split_trace, graph.Split() ) )
                return true;
              }

              // try neighbouring facets
              for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
              {
                int line = *i;
                if (visited[line] < 2 and IsFree(used, free, line))
                {
                  plain_int_set& row = graph[line];
                  for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
                  {
                    int f = *i;
                    if (facet_color[f] == 0 and IsFree(used, free, f))
                    {
                      facet_stack.push_back(f);
                    }
                  }
                }
              }
            }
          }
          else if (facet_color[facet] == 1)  // gray
          {
            // search for any facet within the stack that can be added
            int pos = facet_stack.size() - 1;
            for (; pos >= 0; --pos)
            {
              int f = facet_stack[pos];
              if (facet_color[f] == 0 and IsValidFacet(graph[f], visited))
              {
                facet_stack.push_back(facet);
                facet_stack.push_back(f);
                break;
              }
            }

            // clear current facet if there is nothing more to add
            if (pos < 0)
            {
              plain_int_set& row = graph[facet];
              UnMarkFacet(used, free, facet, row, visited, num_split_lines);

              if (std::find(facet_stack.begin(), facet_stack.end(), facet) == facet_stack.end())
                facet_color[facet] = 0;
              else
                facet_color[facet] = 2;
            }
          }
          else if (facet_color[facet] == 2)  // black
          {
            // if a black facet is poped for the last time (it is not any more
            // on the stack), we make it available again.
            if (std::find(facet_stack.begin(), facet_stack.end(), facet) == facet_stack.end())
            {
              facet_color[facet] = 0;
            }
          }
        }
        return false;
      }

    }  // namespace COLOREDGRAPH
  }    // namespace CUT
}  // namespace GEO

void GEO::CUT::COLOREDGRAPH::Graph::FindFreeFacets(Graph& graph, Graph& used, plain_int_set& free,
    const std::vector<std::pair<Point*, Point*>>& all_lines, std::vector<int>& split_trace)
{
  int free_facet = FindFirstFreeFacet(graph, used, free);

  std::vector<int> visited(graph.size(), 0);

  if (not VisitFacetDFS(graph, used, free, free_facet, all_lines, visited, split_trace))
  {
    throw std::runtime_error("failed to find volume split");
  }

  for (std::vector<int>::iterator i = visited.begin(); i != visited.begin() + graph.Split(); ++i)
  {
    if (*i == 1)
    {
      int facet = i - visited.begin();
      plain_int_set& row = graph[facet];
      for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
      {
        int line = *i;
        used.Add(line, facet);
        Add(line, facet);
        free.erase(line);
      }
      free.erase(facet);
    }
    else if (*i > 1)
    {
      run_time_error("same facet twice");
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::FindSplitTrace(std::vector<int>& split_trace)
{
  for (Graph::const_iterator i = begin(); i != end(); ++i)
  {
    int p = i->first;
    const plain_int_set& row = i->second;
    if (p >= Split())
    {
      if (row.size() == 1)
      {
        split_trace.push_back(p);
      }
    }
  }
}

bool GEO::CUT::COLOREDGRAPH::Graph::ContainsTrace(const std::vector<int>& split_trace)
{
  for (std::vector<int>::const_iterator i = split_trace.begin(); i != split_trace.end(); ++i)
  {
    int p = *i;
    if (count(p) == 0) return false;
  }
  return true;
}

void GEO::CUT::COLOREDGRAPH::Cycle::Print() const
{
  std::cout << "Cycle:\n";
  cycle_.Print();
  std::cout << "\n";
}

void GEO::CUT::COLOREDGRAPH::Graph::Split(Graph& used, plain_int_set& free, Graph& connection,
    const std::vector<int>& split_trace, Graph& c1, Graph& c2)
{
  // find lhs and rhs starting from split trace

  plain_int_set* row = NULL;
  for (std::vector<int>::const_iterator i = split_trace.begin(); i != split_trace.end(); ++i)
  {
    int p = *i;
    row = &at(p);
    if (row->size() == 2)
    {
      break;
    }
  }

  if (row->size() != 2)
  {
    // This might happen and it might be a valid split. How to deal with it?
    run_time_error("expect two facets at line");
  }

  plain_int_set::iterator i = row->begin();
  Fill(split_trace, used, free, connection, *i, c1);
  ++i;
  Fill(split_trace, used, free, connection, *i, c2);

  for (std::vector<int>::const_iterator i = split_trace.begin(); i != split_trace.end(); ++i)
  {
    int p = *i;
    unsigned c1rowlen = c1[p].size();
    unsigned c2rowlen = c2[p].size();

    Graph* fine = NULL;
    Graph* open = NULL;

    if (c1rowlen > 1 and c2rowlen > 1)
    {
      // fine
    }
    else if (c1rowlen == 1 and c2rowlen > 1)
    {
      fine = &c2;
      open = &c1;
    }
    else if (c1rowlen > 1 and c2rowlen == 1)
    {
      fine = &c1;
      open = &c2;
    }
    else
    {
      run_time_error("open line after graph split");
    }

    if (open != NULL)
    {
      plain_int_set& row = at(p);
      if (row.size() != 2) run_time_error("expect two facets at line");
      plain_int_set::iterator i = row.begin();
      int f1 = *i;
      ++i;
      int f2 = *i;
      if (fine->count(f1) > 0)
      {
        open->Add(f2, at(f2));
      }
      else if (fine->count(f2) > 0)
      {
        open->Add(f1, at(f1));
      }
      else
      {
        run_time_error("confused");
      }
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::Fill(const std::vector<int>& split_trace, Graph& used,
    plain_int_set& free, Graph& connection, int seed, Graph& c)
{
  plain_int_set done;
  done.insert(split_trace.begin(), split_trace.end());
  std::stack<int> stack;

  stack.push(seed);

  while (not stack.empty())
  {
    int f = stack.top();
    stack.pop();

    done.insert(f);

    plain_int_set& row = at(f);
    for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
    {
      int p = *i;
      c.Add(f, p);
      if (done.count(p) == 0)  // and IsFree( used, free, p ) )
      {
        plain_int_set& row = at(p);
        for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
        {
          int f = *i;
          if (done.count(f) == 0)
          {
            stack.push(f);
          }
        }
      }
    }
  }

  for (Graph::const_iterator i = connection.begin(); i != connection.end(); ++i)
  {
    int p = i->first;
    const plain_int_set& row = i->second;
    c.Add(p, row);
  }
}

void GEO::CUT::COLOREDGRAPH::CycleList::AddPoints(Graph& graph, Graph& used, Graph& cycle,
    plain_int_set& free, const std::vector<std::pair<Point*, Point*>>& all_lines)
{
  PushBack(cycle);

  while (free.size() > 0)
  {
    Graph connection(graph.Split());

    // find connection graph and trace lines
    std::vector<int> split_trace;
    connection.FindFreeFacets(graph, used, free, all_lines, split_trace);

    // There might be multiple matches. Only one of those is the one we are
    // looking for.
    std::vector<std::list<Cycle>::iterator> matching;
    for (std::list<Cycle>::iterator i = cycles_.begin(); i != cycles_.end(); ++i)
    {
      Cycle& c = *i;
      if (c.ContainsTrace(split_trace))
      {
        matching.push_back(i);
      }
    }

    bool found = false;

    for (std::vector<std::list<Cycle>::iterator>::iterator ilist = matching.begin();
         ilist != matching.end(); ++ilist)
    {
      Cycle& c = **ilist;

      Graph c1(graph.Split());
      Graph c2(graph.Split());

#ifdef DEBUGCUTLIBRARY
      c().DumpGraph("cycle.py");
      connection.DumpGraph("connection.py");
#endif

      c.Split(used, free, connection, split_trace, c1, c2);

#ifdef DEBUGCUTLIBRARY
      c1.DumpGraph("cycle1.py");
      c2.DumpGraph("cycle2.py");
#endif

      if (c1 == c2)
      {
        if (matching.size() == 1)
        {
          PushBack(c1);

          run_time_error("bad luck");

          cycles_.erase(*ilist);
          found = true;
          break;
        }
      }
      else
      {
        // sanity test
        for (Graph::const_iterator i = c().begin(); i != c().end(); ++i)
        {
          int f = i->first;
          if (f >= c().Split()) break;
          if (connection.count(f) == 0 and c1.count(f) > 0 and c2.count(f) > 0)
          {
            run_time_error("not a valid split");
          }
        }

        PushBack(c1);
        PushBack(c2);

        cycles_.erase(*ilist);
        found = true;
        break;
      }
    }
    if (not found)
    {
      throw std::runtime_error("did not find volume that contains split facets");
    }
  }
}

void GEO::CUT::COLOREDGRAPH::CycleList::PushBack(Graph& g)
{
  cycles_.push_back(Cycle(g.Split()));
  Cycle& c = cycles_.back();
  c.Assign(g);
}

void GEO::CUT::COLOREDGRAPH::CycleList::Print() const
{
  for (std::list<Cycle>::const_iterator i = cycles_.begin(); i != cycles_.end(); ++i)
  {
    const Cycle& c = *i;
    c.Print();
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::TestSplit()
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int p = i->first;
    if (p >= color_split_)
    {
      plain_int_set& row = i->second;
      if (row.size() > 2)
      {
#ifdef DEBUGCUTLIBRARY
        DumpGraph("failedgraph.py");
#endif
        run_time_error("colored graph not properly split");
      }
    }
  }
}

void GEO::CUT::COLOREDGRAPH::Graph::DumpGraph(const std::string& name)
{
  std::ofstream file(name.c_str());
  file << "color_split = " << Split() << "\n";
  file << "graph = [";
  for (const_iterator i = begin(); i != end(); ++i)
  {
    file << i->first << ",";
  }
  file << "]\n";
  file << "data = {\n";
  for (const_iterator i = begin(); i != end(); ++i)
  {
    file << "    " << i->first << ": [";
    std::copy(i->second.begin(), i->second.end(), std::ostream_iterator<int>(file, ","));
    file << "],\n";
  }
  file << "}\n";
}
