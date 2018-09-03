
#include "graph_test.H"

class SpanningTreeCreator : public boost::default_bfs_visitor
{
 public:
  SpanningTreeCreator(Graph& st) : st_(st) {}

  void tree_edge(edge_t e, const Graph& g)
  {
    vertex_t u = boost::source(e, g);
    vertex_t v = boost::target(e, g);

    if (u > v)
    {
      std::swap(u, v);
    }

    // std::cout << "tree_edge: " << e << "\n";

    // assume same vertex ids
    boost::add_edge(u, v, st_);
  }

 private:
  Graph& st_;
};


void create_spanning_tree(Graph& g, Graph& st)
{
  name_map_t g_name_map = boost::get(boost::vertex_name, g);
  name_map_t st_name_map = boost::get(boost::vertex_name, st);

  color_map_t color_map = boost::get(boost::vertex_color, g);

  vertex_iterator vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
  {
    color_map[*vi] = color_t::white();
    st_name_map[*vi] = g_name_map[*vi];
  }

  for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
  {
    vertex_t src = *vi;
    if (color_map[src] == color_t::white())
    {
      boost::breadth_first_visit(
          g, src, boost::color_map(color_map).visitor(SpanningTreeCreator(st)));
    }
  }

  //   boost::print_graph( st, st_name_map );
}

cycle_t* intersect(std::set<cycle_t*>& u, std::set<cycle_t*>& v)
{
  std::set<cycle_t*> intersection;
  std::set_intersection(
      u.begin(), u.end(), v.begin(), v.end(), std::inserter(intersection, intersection.begin()));

#ifdef DEBUG
  std::cout << "intersect:\n";

  for (std::set<cycle_t*>::iterator i = intersection.begin(); i != intersection.end(); ++i)
  {
    cycle_t* c = *i;
    print_cycle(c);
  }
#endif

  if (intersection.size() == 1)
  {
    return *intersection.begin();
  }
  else if (intersection.size() == 0)
  {
    return NULL;
  }

  throw std::runtime_error("more than one cycle in common");
}


class SpanningTreePathFinder : public boost::default_dfs_visitor
{
 public:
  SpanningTreePathFinder(vertex_t v, std::vector<vertex_t>& path) : done_(false), v_(v), path_(path)
  {
  }

  void discover_vertex(vertex_t u, const Graph& g)
  {
    if (not done_)
    {
      if (u == v_) done_ = true;
      path_.push_back(u);
    }
  }

  bool operator()(vertex_t u, const Graph& g) { return done_; }

  void finish_vertex(vertex_t u, const Graph& g)
  {
    if (not done_)
    {
      if (path_.back() != u) throw std::runtime_error("confused");
      path_.pop_back();
    }
  }

 private:
  bool done_;
  vertex_t v_;
  std::vector<vertex_t>& path_;
};

cycle_t* find_path(vertex_t u, vertex_t v, Graph& st)
{
  color_map_t color_map = boost::get(boost::vertex_color, st);

  vertex_iterator vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(st); vi != vi_end; ++vi)
  {
    color_map[*vi] = color_t::white();
  }

  std::vector<vertex_t>* path = new std::vector<vertex_t>();
  SpanningTreePathFinder vis(v, *path);
  boost::depth_first_visit(st, u, vis, color_map, vis);

#ifdef DEBUG
  std::cout << "find_path between " << u << " and " << v << ":\n";
  print_cycle(path);
#endif

  return path;
}

void add_edges(const cycle_t& path, std::set<std::pair<vertex_t, vertex_t>>& edges)
{
  for (unsigned i = 0; i != path.size(); ++i)
  {
    unsigned j = (i + 1) % path.size();
    vertex_t n1 = path.at(i);
    vertex_t n2 = path.at(j);
    if (n1 > n2) std::swap(n1, n2);
    std::pair<vertex_t, vertex_t> e = std::make_pair(n1, n2);
    if (edges.count(e) == 0)
    {
      edges.insert(e);
    }
    else
    {
      edges.erase(e);
    }
  }
}

cycle_t* substract(cycle_t* path, const std::set<cycle_t*>& partial)
{
  std::set<std::pair<vertex_t, vertex_t>> edges;

  add_edges(*path, edges);

  for (std::set<cycle_t*>::const_iterator i = partial.begin(); i != partial.end(); ++i)
  {
    cycle_t* c = *i;
    add_edges(*c, edges);
  }

  std::set<vertex_t> vertices;
  for (std::set<std::pair<vertex_t, vertex_t>>::iterator i = edges.begin(); i != edges.end(); ++i)
  {
    vertices.insert(i->first);
    vertices.insert(i->second);
  }

  cycle_t* newpath = new cycle_t;

  for (cycle_t::iterator i = path->begin(); i != path->end(); ++i)
  {
    vertex_t u = *i;
    if (vertices.count(u) > 0)
    {
      newpath->push_back(u);
    }
  }

#ifdef DEBUG
  std::cout << "substract:\n";

  for (std::set<std::pair<vertex_t, vertex_t>>::iterator i = edges.begin(); i != edges.end(); ++i)
  {
    std::cout << "(" << i->first << "," << i->second << ") ";
  }
  std::cout << "\n";

  print_cycle(path);
  print_cycle(newpath);
  for (std::set<cycle_t*>::const_iterator i = partial.begin(); i != partial.end(); ++i)
  {
    cycle_t* c = *i;
    print_cycle(c);
  }
#endif

  delete path;
  return newpath;
}

bool issuperset(const std::set<vertex_t>& pathset, const cycle_t& c)
{
  for (cycle_t::const_iterator i = c.begin(); i != c.end(); ++i)
  {
    vertex_t n = *i;
    if (pathset.count(n) == 0)
    {
      return false;
    }
  }
  return true;
}

void split(const cycle_t* cycle_to_split, vertex_t u, vertex_t v, cycle_t& cycle1, cycle_t& cycle2)
{
  cycle_t::const_iterator p1 = std::find(cycle_to_split->begin(), cycle_to_split->end(), u);
  cycle_t::const_iterator p2 = std::find(cycle_to_split->begin(), cycle_to_split->end(), v);

  if (p1 == cycle_to_split->end() or p2 == cycle_to_split->end())
    throw std::runtime_error("edge vertex not in cycle");

  if (p1 > p2)
  {
    std::swap(p1, p2);
    std::swap(u, v);
  }

  cycle1.reserve(p1 - cycle_to_split->begin() + 1 + cycle_to_split->end() - p2);
  cycle2.reserve(p2 - p1 + 1);

  std::copy(cycle_to_split->begin(), p1, std::back_inserter(cycle1));
  cycle1.push_back(u);
  std::copy(p2, cycle_to_split->end(), std::back_inserter(cycle1));

  std::copy(p1, p2, std::back_inserter(cycle2));
  cycle2.push_back(v);

#ifdef DEBUG
  std::cout << "split:\n";
  print_cycle(cycle_to_split);
  print_cycle(&cycle1);
  print_cycle(&cycle2);
#endif
}

void find_cycles(
    Graph& g, std::map<int, std::vector<double>>& coords, std::set<cycle_t*>& base_cycles)
{
  Graph st(boost::num_vertices(g));
  std::set<std::pair<vertex_t, vertex_t>> st_edges;

  create_spanning_tree(g, st);

  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(st); ei != ei_end; ++ei)
  {
    edge_t e = *ei;
    vertex_t u = boost::source(e, st);
    vertex_t v = boost::target(e, st);

    if (u > v)
    {
      std::swap(u, v);
    }

    st_edges.insert(std::make_pair(u, v));
  }

  gnuplot_graph("st.plot", st, coords);

  //   std::copy( non_tree_edges.begin(), non_tree_edges.end(), std::ostream_iterator<edge_t>(
  //   std::cout, " " ) ); std::cout << "\n";

  std::vector<std::set<cycle_t*>> cycles(boost::num_vertices(g));

  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    edge_t e = *ei;
    vertex_t u = boost::source(e, g);
    vertex_t v = boost::target(e, g);

    if (u > v)
    {
      std::swap(u, v);
    }

    if (st_edges.count(std::make_pair(u, v)) > 0) continue;

#ifdef DEBUG
    std::cout << "edge: " << u << "," << v << "\n";
#endif

    cycle_t* cycle_to_split = intersect(cycles[u], cycles[v]);
    if (cycle_to_split == NULL)
    {
      cycle_t* path = find_path(u, v, st);

      std::set<cycle_t*> known;
      for (cycle_t::iterator i = path->begin(); i != path->end(); ++i)
      {
        vertex_t n = *i;
        std::set<cycle_t*>& c = cycles[n];
        std::copy(c.begin(), c.end(), std::inserter(known, known.begin()));
      }

      std::set<vertex_t> pathset;
      std::copy(path->begin(), path->end(), std::inserter(pathset, pathset.begin()));

      std::set<cycle_t*> partial;
      for (std::set<cycle_t*>::iterator i = known.begin(); i != known.end(); ++i)
      {
        cycle_t* c = *i;
        if (issuperset(pathset, *c))
        {
          partial.insert(c);
        }
      }

      if (partial.size() > 0)
      {
        path = substract(path, partial);
      }

      for (cycle_t::iterator i = path->begin(); i != path->end(); ++i)
      {
        vertex_t n = *i;
        std::set<cycle_t*>& c = cycles[n];
        c.insert(path);
      }
    }
    else
    {
      cycle_t* cycle1 = new cycle_t;
      cycle_t* cycle2 = new cycle_t;
      split(cycle_to_split, u, v, *cycle1, *cycle2);
      for (cycle_t::iterator i = cycle_to_split->begin(); i != cycle_to_split->end(); ++i)
      {
        vertex_t n = *i;
        std::set<cycle_t*>& c = cycles[n];
        c.erase(cycle_to_split);
        if (std::find(cycle1->begin(), cycle1->end(), n) != cycle1->end())
        {
          c.insert(cycle1);
        }
        if (std::find(cycle2->begin(), cycle2->end(), n) != cycle2->end())
        {
          c.insert(cycle2);
        }
      }
      delete cycle_to_split;
    }
  }

  for (std::vector<std::set<cycle_t*>>::iterator i = cycles.begin(); i != cycles.end(); ++i)
  {
    std::set<cycle_t*>& cs = *i;
    for (std::set<cycle_t*>::iterator i = cs.begin(); i != cs.end(); ++i)
    {
      cycle_t* c = *i;
      base_cycles.insert(c);
    }
    cs.clear();
  }
  cycles.clear();
}

void cleanup(Graph& g, std::map<int, std::vector<double>>& coords, std::set<cycle_t*>& base_cycles)
{
#ifdef DEBUG
  std::cout << "result:\n";

  int count = 0;
#endif

  for (std::set<cycle_t*>::iterator i = base_cycles.begin(); i != base_cycles.end(); ++i)
  {
    cycle_t* c = *i;

#ifdef DEBUG
    print_cycle(c);

    std::stringstream str;
    str << "cycle" << count << ".plot";
    gnuplot_cycle(str.str().c_str(), g, c, coords);
    count += 1;
#endif

    delete c;
  }
}
