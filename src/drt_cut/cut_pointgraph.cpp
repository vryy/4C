/*---------------------------------------------------------------------*/
/*! \file
\brief PointGraph, Graph Algorithm to create Facets from lines and edges

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include <iostream>
#include <iterator>

#include <cmath>

#include "cut_pointgraph.H"
#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_output.H"

#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/graphviz.hpp>

#include "cut_pointgraph_simple.H"


#define DEBUG_POINTGRAPH false
//#define CLN_CALC_OUTSIDE_KERNEL
#ifdef CLN_CALC_OUTSIDE_KERNEL
#include "cut_clnwrapper.H"
#endif

/*----------------------------------------------------------------------------*
 * Constructor for the selfcut                                     wirtz 05/13
 *----------------------------------------------------------------------------*/
GEO::CUT::IMPL::PointGraph::PointGraph(Side *side) : graph_(CreateGraph(side->Dim()))
{
  Cycle cycle;
  FillGraph(side, cycle);
  if (GetGraph().HasSinglePoints(element_side))  // if any edge in graph has single point
  {
    GetGraph().FixSinglePoints(cycle);  // delete single point edges
  }
  GetGraph().FindCycles(side, cycle);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::IMPL::PointGraph::PointGraph(
    Mesh &mesh, Element *element, Side *side, Location location, Strategy strategy)
    : graph_(CreateGraph(element->Dim()))
{
  // here we create the facets...
  Cycle cycle;
  FillGraph(element, side, cycle, strategy);
#if 1
#ifdef DEBUGCUTLIBRARY
  {
    std::ofstream f("all_points0.plot");
    GetGraph().PlotAllPoints(f);
  }
  {
    std::ofstream f("graph0.txt");
    GetGraph().Print(f);
  }
  {
    std::ofstream f("cycle0.txt");
    f << cycle;
  }
#endif
#endif

  // if any edge in graph has single point
  if (GetGraph().HasSinglePoints(location))
  {
    // NOTE: levelset method does not have complicated check for the single point. Feel free to
    // exend it
    if (side->IsLevelSetSide() or GetGraph().SimplifyConnections(element, side) or
        (GetGraph().HasTouchingEdge(element, side)))
    {
#if DEBUG_POINTGRAPH
      std::cout << "WARNING: Deleting single point in the pointgraph" << std::endl;
#endif
      GetGraph().FixSinglePoints(cycle);  // delete single point edges
    }
    else
    {
      std::ofstream f("graph0.txt");
      GetGraph().Print(f);
      dserror("Pointgraph has single point.This shouldn't happen or we should understand why!");
    }
    //    GetGraph().TestClosed();
  }

#ifdef DEBUGCUTLIBRARY
  {
    std::ofstream f("all_points.plot");
    GetGraph().PlotAllPoints(f);
  }
  {
    std::ofstream f("graph.txt");
    GetGraph().Print(f);
  }
  {
    std::ofstream f("cycle.txt");
    f << cycle;
  }
#endif

  try
  {
    GetGraph().FindCycles(element, side, cycle, location, strategy);
  }
  catch (std::runtime_error &err)
  {
    std::ofstream file("failed_pointgraph.pos");
    GEO::CUT::OUTPUT::GmshSideDump(file, side, std::string("Side"));

    // add cut lines to graph
    const std::vector<Line *> &cut_lines = side->CutLines();

    for (std::vector<Line *>::const_iterator i = cut_lines.begin(); i != cut_lines.end(); ++i)
    {
      int line_index = i - cut_lines.begin();
      std::stringstream section_name;
      section_name << "Cut_lines" << line_index;
      Line *l = *i;
      GEO::CUT::OUTPUT::GmshNewSection(file, section_name.str());
      GEO::CUT::OUTPUT::GmshLineDump(file, l, false, NULL);
      GEO::CUT::OUTPUT::GmshEndSection(file, false);
      // output distance between points of the line
      file << "// Distance between points of the line is"
           << GEO::CUT::DistanceBetweenPoints(l->BeginPoint(), l->EndPoint()) << std::endl;
    }
    file.close();
    throw err;
  }
#if 0
  cycle.Print();
  std::cout << "Main-Cycles\n";
  for ( std::vector<Cycle>::const_iterator cit = GetGraph().main_cycles_.begin();
        cit != GetGraph().main_cycles_.end(); ++cit )
    cit->Print();
#endif
}

/*-------------------------------------------------------------------------------------*
 * Graph is filled wihl all edges of the selfcut: uncutted edges, selfcutedges
 * and new splitted edges; but no the cutted edges                          wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::FillGraph(Side *side, Cycle &cycle)
{
  const std::vector<Node *> &nodes = side->Nodes();
  const std::vector<Edge *> &edges = side->Edges();
  int end_pos = 0;

  // loop over all edges of the parent side
  for (std::vector<Edge *>::const_iterator i = edges.begin(); i != edges.end(); ++i)
  {
    Edge *e = *i;

    // get start and end node numbers corresponding to this edge
    int begin_pos = end_pos;
    end_pos = (end_pos + 1) % nodes.size();
    std::vector<Point *> edge_points;

    // get all points on this edge including start and end points
    // points are already sorted
    e->CutPoint(nodes[begin_pos], nodes[end_pos], edge_points);

    // when edge of a side has "n" cut points, the edge itself is split into (n+1) edges
    // store all (n+1) edges to graph
    for (unsigned i = 1; i < edge_points.size(); ++i)  // no of edges = no of points-1
    {
      Point *p1 = edge_points[i - 1];
      Point *p2 = edge_points[i];
      GetGraph().AddEdge(p1, p2);
    }
    for (std::vector<Point *>::iterator i = edge_points.begin() + 1; i != edge_points.end(); ++i)
    {
      Point *p = *i;
      cycle.push_back(p);
    }
  }
  const plain_edge_set &selfcutedges = side->SelfCutEdges();
  for (plain_edge_set::const_iterator i = selfcutedges.begin(); i != selfcutedges.end(); ++i)
  {
    Edge *selfcutedge = *i;
    GetGraph().AddEdge(selfcutedge->BeginNode()->point(), selfcutedge->EndNode()->point());
  }
}

/*----------------------------------------------------------------------------*
 * Get all edges created on this side after cut, store cycle of points on this
 * side to create facet. Also add cut lines to the graph
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::FillGraph(
    Element *element, Side *side, Cycle &cycle, Strategy strategy)
{
  const std::vector<Node *> &nodes = side->Nodes();
  const std::vector<Edge *> &edges = side->Edges();
  int end_pos = 0;

#if DEBUG_POINTGRAPH
  std::cout << "Filling graph" << std::endl;
#endif

  // loop over all edges of the parent side
  for (std::vector<Edge *>::const_iterator i = edges.begin(); i != edges.end(); ++i)
  {
    Edge *e = *i;

#if DEBUG_POINTGRAPH
    int index = i - edges.begin();
    std::cout << "Processing edge with index " << index << " and Id=" << e->Id() << std::endl;
#endif

    // get start and end node numbers corresponding to this edge
    int begin_pos = end_pos;
    end_pos = (end_pos + 1) % nodes.size();

    std::vector<Point *> edge_points;

    // get all points on this edge including start and end points
    // points are already sorted
    e->CutPoint(nodes[begin_pos], nodes[end_pos], edge_points);

#if DEBUG_POINTGRAPH
    std::cout << "Number of points on the current edge is " << edge_points.size() << std::endl;
#endif

    // when edge of a side has "n" cut points, the edge itself is split into (n+1) edges
    // store all (n+1) edges to graph
    for (unsigned i = 1; i < edge_points.size(); ++i)  // number of edges = number of points-1
    {
      Point *p1 = edge_points[i - 1];
      Point *p2 = edge_points[i];
#if DEBUG_POINTGRAPH
      std::cout << "Adding line betweeen points with ids " << p1->Id() << " and " << p2->Id()
                << std::endl;
#endif
      GetGraph().AddEdge(p1, p2);
    }

    BuildCycle(edge_points, cycle);
  }

  AddCutLinesToGraph(element, side, strategy);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::BuildCycle(
    const std::vector<Point *> &edge_points, Cycle &cycle) const
{
  for (std::vector<Point *>::const_iterator i = edge_points.begin() + 1; i != edge_points.end();
       ++i)
  {
    Point *p = *i;
    cycle.push_back(p);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::AddCutLinesToGraph(Element *element, Side *side, Strategy strategy)
{
  const std::vector<Line *> &cut_lines = side->CutLines();

  // add cut lines to graph
#if DEBUG_POINTGRAPH
  std::cout << "Adding cut lines to the graph " << std::endl;
#endif
  for (std::vector<Line *>::const_iterator i = cut_lines.begin(); i != cut_lines.end(); ++i)
  {
    Line *l = *i;

    // GetGraph().AddEdge(l->BeginPoint(), l->EndPoint());
    bool element_cut = l->IsCut(element);
    if (strategy == all_lines or element_cut) GetGraph().AddEdge(l->BeginPoint(), l->EndPoint());
#if DEBUG_POINTGRAPH
    l->BeginPoint()->Print();
    l->EndPoint()->Print();
#endif


#ifdef DEBUGCUTLIBRARY
    if (element_cut)
    {
      Point *p1 = l->BeginPoint();
      Point *p2 = l->EndPoint();
      if (not p1->IsCut(element) or not p2->IsCut(element))
      {
        std::stringstream str;
        str << "line between " << (*p1) << " and " << (*p2)
            << " is cut by element, but point cuts are: " << p1->IsCut(element) << " and "
            << p2->IsCut(element);
        throw std::runtime_error(str.str());
      }
    }
#endif
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::AddEdge(int row, int col)
{
  graph_[row].insert(col);
  graph_[col].insert(row);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::AddEdge(Point *p1, Point *p2)
{
  all_points_[p1->Id()] = p1;
  all_points_[p2->Id()] = p2;

  AddEdge(p1->Id(), p2->Id());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::Print(std::ostream &stream)
{
  stream << "--- PointGraph::Graph ---\n";
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int p = i->first;
    plain_int_set &row = i->second;
    stream << p << ": ";
    for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
    {
      int p = *i;
      stream << p << " ";
    }
    stream << "\n";
  }
  stream << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::PlotAllPoints(std::ostream &stream)
{
  for (std::map<int, Point *>::iterator i = all_points_.begin(); i != all_points_.end(); ++i)
  {
    i->second->Plot(stream);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::PlotPoints(Element *element)
{
  for (std::map<int, Point *>::iterator i = all_points_.begin(); i != all_points_.end(); ++i)
  {
    Point *p = i->second;
    std::cout << p->Id() << "(" << p->IsCut(element) << ") ";
  }
  std::cout << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IMPL::FindCycles(graph_t &g, GEO::CUT::Cycle &cycle,
    std::map<vertex_t, LINALG::Matrix<3, 1>> &local,
    std::vector<GEO::CUT::Cycle> &cycles) /* non-member function */
{
  name_map_t name_map = boost::get(boost::vertex_name, g);

  // Initialize the interior edge index
  edge_index_map_t e_index = boost::get(boost::edge_index, g);
  boost::graph_traits<graph_t>::edges_size_type edge_count = 0;
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  // updating property map of edges with indexes
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    boost::put(e_index, *ei, edge_count++);

  typedef std::vector<edge_t> vec_t;
  std::vector<vec_t> embedding(boost::num_vertices(g));


  // Use geometry to build embedding. The only safe way to do it.

#ifdef CLN_CALC_OUTSIDE_KERNEL
  typedef ClnWrapper floatType;
  // NOTE: Cln can be used, if one get problem with double arc and there is no other way to fix it
  // However, if running cln with as custom memory manager, this should be changed ( mostly to free
  // objects in a container, similarly as in cut_kernel )
#else
  typedef double floatType;
#endif

  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
  {
    const LINALG::Matrix<3, 1> &pos = local[*vi];
#if DEBUG_POINTGRAPH
    std::cout << "First coordinate before substraction " << std::setprecision(16) << pos
              << std::endl;
#endif

#if DEBUG_POINTGRAPH
    std::cout << "First point is " << name_map[*vi]->Id() << std::endl;
#endif
    std::map<floatType, vertex_t> arcs;
    adjacency_iterator ai, ai_end;
    for (boost::tie(ai, ai_end) = boost::adjacent_vertices(*vi, g); ai != ai_end; ++ai)
    {
      LINALG::Matrix<3, 1> d = local[*ai];

#if DEBUG_POINTGRAPH
      std::cout << "Adjacent point is " << name_map[*ai]->Id() << std::endl;
      std::cout << "Second coordinate before substrication " << std::setprecision(16) << d
                << std::endl;
#endif
      d.Update(-1, pos, 1);

#ifdef CLN_CALC_OUTSIDE_KERNEL
      //                      Order of arguments is  __x___  ___y___
      floatType arc = cln::atan(cln::cl_float(d(0), cln::float_format(ClnWrapper::precision_)),
          cln::cl_float(d(1), cln::float_format(ClnWrapper::precision_)));
#else
      //                    Order of arguments is  __y___  ___x___
      floatType arc = std::atan2(d(1), d(0));
#endif

#if DEBUG_POINTGRAPH
      std::cout << "Arc is equal to " << arc << std::endl;
#endif

      std::map<floatType, vertex_t>::iterator j = arcs.find(arc);

      if (j != arcs.end())
      {
        // this can occure once when more than one nodes of the background element
        // has same coordinates (sudhakar)
        // check input file for two nodes (in same domain) having same coordinates
        std::stringstream err_msg;
        Point *first = name_map[*vi];
        Point *second = name_map[*ai];
        Point *previous = name_map[j->second];

        std::ofstream file("double_arc.pos");

        GEO::CUT::OUTPUT::GmshNewSection(file, "NewLine");
        GEO::CUT::OUTPUT::GmshLineDump(file, first, second, first->Id(), second->Id(), false, NULL);
        GEO::CUT::OUTPUT::GmshEndSection(file);
        GEO::CUT::OUTPUT::GmshNewSection(file, "OldLine");
        GEO::CUT::OUTPUT::GmshLineDump(
            file, first, previous, first->Id(), previous->Id(), false, NULL);
        GEO::CUT::OUTPUT::GmshEndSection(file, true);

        err_msg
            << "Numerical error: double arc when trying to create arc with between points with Id="
            << first->Id() << " and  " << second->Id()
            << "!. Arc of the same length exists between Ids " << first->Id() << " and   "
            << previous->Id() << std::endl;
        file.close();

        run_time_error(err_msg.str());
      }

      arcs[arc] = *ai;
    }

    vec_t &em = embedding[*vi];

// NOTE: We want to have embedding with clockwise ordering of edges. Otherwise it will produce an
// error
#if DEBUG_POINTGRAPH
    std::cout << "For vertex " << name_map[*vi]->Id() << " planar graph is edges with indexes: ";
#endif

    for (std::map<floatType, vertex_t>::iterator i = arcs.begin(); i != arcs.end(); ++i)
    {
      out_edge_iterator oi, oi_end;
      for (boost::tie(oi, oi_end) = boost::out_edges(*vi, g); oi != oi_end; ++oi)
      {
        edge_t e = *oi;
        if (boost::target(e, g) == i->second)
        {
#if DEBUG_POINTGRAPH
          std::cout << e_index[e] << " ; ";
#endif
          em.push_back(e);
          break;
        }
      }
    }
#if DEBUG_POINTGRAPH
    std::cout << "\n";
#endif
  }
#ifdef CLN_CALC_OUTSIDE_KERNEL
  // ClnWrapper::precision_ = 30;
#endif

#if 0
  // this check actually modifies the graph!  if the input graph is not planar, it create kuratowski subgraph
  if ( not boost::boyer_myrvold_planarity_test( boost::boyer_myrvold_params::graph = g,
                                                boost::boyer_myrvold_params::embedding = &embedding[0] ) )
  // we can better use
  if (not boost::boyer_myrvold_planarity_test( g ) )
  {
    throw std::runtime_error( "input graph is not planar" );
  }
#endif
  //#endif

#if 0
#ifdef DEBUGCUTLIBRARY
  std::cout << "embedding:\n";
  for ( std::vector<vec_t>::iterator i=embedding.begin(); i!=embedding.end(); ++i )
  {
    vec_t & em = *i;
    std::copy( em.begin(), em.end(), std::ostream_iterator<edge_t>( std::cout, " " ) );
    std::cout << "\n";
  }
#endif
#endif

  face_visitor vis(name_map, cycles);
  boost::planar_face_traversal(g, &embedding[0], vis);


//#ifdef DEBUGCUTLIBRARY
#if DEBUG_POINTGRAPH
  for (std::vector<Cycle>::iterator i = cycles.begin(); i != cycles.end(); ++i)
  {
    Cycle &c = *i;
    c.TestUnique();
  }
#endif
  //#endif


  // boost face traversal will produce two cycles, in case if there  is one planar face ( surface
  // with no cut lines ), hence we need to remove redundant one and the other will serve as a facet
  // for us in case of normal configuration (more then one planar face),  we need to remove all the
  // instances of the full cycle produced by boost face traversal (which should be one), since the
  // set of small cycles would be the one creating facets

  bool save_first = cycles.size() == 2;

  int erase_count = 0;
  for (std::vector<Cycle>::iterator i = cycles.begin(); i != cycles.end();)
  {
    Cycle &c = *i;
    if (cycle.Equals(c))
    {
      if (save_first and erase_count == 0)
      {
        ++i;
      }
      else
      {
        cycles.erase(i);
      }
      erase_count += 1;
    }
    else
    {
      ++i;
    }
  }

  if (erase_count > (save_first ? 2 : 1))
  {
    throw std::runtime_error("more than one back facet");
  }

#if DEBUG_POINTGRAPH
  if (erase_count == 0)
  {
    std::cout << "ERASED 0 cycles ( no main cycle in the pointgraph)" << std::endl;
    std::cout << "Number of cycles is" << cycles.size() << std::endl;
    int counter = 0;
    for (std::vector<Cycle>::iterator i = cycles.begin(); i != cycles.end(); ++i, ++counter)
    {
      Cycle &c = *i;
      std::stringstream s;
      s << "Cycle_" << counter << ".pos";
      std::ofstream file(s.str());
      c.GmshDump(file);
      file.close();
    }
    std::ofstream file("main_cycle.pos");
    cycle.GmshDump(file);
    file.close();
  }
#endif

  return erase_count != 0;
}

/*-------------------------------------------------------------------------------------*
 * Creates maincycles (outer polygons) and holecycles (inner polygons = holes)
 * of the selfcut graph                                                     wirtz 05/13
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::FindCycles(Side *side, Cycle &cycle)
{
  graph_t g;

  // create boost graph

  name_map_t name_map = boost::get(boost::vertex_name, g);
  edge_index_map_t edge_index_map = boost::get(boost::edge_index, g);

  std::map<int, vertex_t> vertex_map;

  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int n = i->first;

    Point *p = GetPoint(n);
    vertex_t u = add_vertex(g);
    name_map[u] = p;
    vertex_map[n] = u;
  }

  int counter = 0;

  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int u = i->first;

    //    Point * p1 = GetPoint( u );

    plain_int_set &row = i->second;

    for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
    {
      int v = *i;
      //      Point * p2 = GetPoint( v );

      if (u < v)
      {
        edge_t e;
        bool inserted;
        boost::tie(e, inserted) = boost::add_edge(vertex_map[u], vertex_map[v], g);
        if (inserted)
        {
          edge_index_map[e] = counter;
          counter += 1;
        }
      }
    }
  }

  // All vertices are connected. If there is no cycle, done.
  if (boost::num_vertices(g) > boost::num_edges(g))
  {
    return;
  }


  // Use geometry to find the right embedding and find the cycles.
  // find local coordinates

  std::map<vertex_t, LINALG::Matrix<3, 1>> local;

  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
  {
    // prepare vars
    Point *p = name_map[*vi];
    LINALG::Matrix<3, 1> xyz(p->X());
    LINALG::Matrix<3, 1> tmpmat;

    // get coords
    side->LocalCoordinates(xyz, tmpmat);

    // add to map
    std::pair<vertex_t, LINALG::Matrix<3, 1>> tmppair(*vi, tmpmat);
    local.insert(tmppair);
  }

  // find unconnected components (main facet(s) and holes)

  std::vector<int> component(boost::num_vertices(g));

  int num_comp = boost::connected_components(
      g, boost::make_iterator_property_map(component.begin(), boost::get(boost::vertex_index, g)));

  // find cycles on each component

  if (num_comp == 1)
  {
    GEO::CUT::IMPL::FindCycles(g, cycle, local, main_cycles_);
  }
  else if (num_comp > 1)
  {
    for (int i = 0; i < num_comp; ++i)
    {
      typedef boost::filtered_graph<graph_t, edge_filter> filtered_graph_t;
      edge_filter filter(g, component, i);
      filtered_graph_t fg(g, filter);

      std::vector<Cycle> filtered_cycles;

      graph_t cg;
      boost::copy_graph(fg, cg);

      bool main_cycle = GEO::CUT::IMPL::FindCycles(cg, cycle, local, filtered_cycles);

      if (main_cycle)
      {
        if (main_cycles_.size() != 0)
        {
          run_time_error("one set of main cycles only");
        }
        std::swap(main_cycles_, filtered_cycles);
      }
      else
      {
        hole_cycles_.push_back(std::vector<Cycle>());
        std::swap(hole_cycles_.back(), filtered_cycles);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::FindCycles(
    Element *element, Side *side, Cycle &cycle, Location location, Strategy strategy)
{
  graph_t g;

  // create boost graph

  name_map_t name_map = boost::get(boost::vertex_name, g);
  edge_index_map_t edge_index_map = boost::get(boost::edge_index, g);

  std::map<int, vertex_t> vertex_map;

  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int n = i->first;

    Point *p = GetPoint(n);

    // this check if it is  not levelset add all points and remove later, otherwise use this
    // strategy
    if (!(strategy == own_lines) or ((location == element_side) or (p->IsCut(element))))
    {
      vertex_t u = add_vertex(g);
      name_map[u] = p;
      vertex_map[n] = u;
    }
  }

  int counter = 0;
#if DEBUG_POINTGRAPH
  std::cout << "\n";
#endif
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    int u = i->first;

    Point *p1 = GetPoint(u);

    if (!(strategy == own_lines) or ((location == element_side) or (p1->IsCut(element))))
    {
      plain_int_set &row = i->second;

      for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
      {
        int v = *i;
        Point *p2 = GetPoint(v);

        if (!(strategy == own_lines) or ((location == element_side) or (p2->IsCut(element))))
        {
          if (u < v)
          {
            edge_t e;
            bool inserted;
            boost::tie(e, inserted) = boost::add_edge(vertex_map[u], vertex_map[v], g);
            if (inserted)
            {
#if DEBUG_POINTGRAPH
              std::cout << "Inserting edge with edge_index " << counter << " between points "
                        << p1->Id() << " and " << p2->Id() << std::endl;
#endif
              edge_index_map[e] = counter;
              counter += 1;
            }
          }
        }
      }
    }
  }

  // All vertices are connected. If there is no cycle, done.
  if (boost::num_vertices(g) > boost::num_edges(g))
  {
#if 0
#ifdef DEBUGCUTLIBRARY
    if ( boost::num_vertices( g ) > 2 )
    {
      std::cout << "failed graph: num_vertices=" << boost::num_vertices( g )
                << "   num_edges=" << boost::num_edges( g )
                << "\n"
                << cycle << "\n";
      boost::print_graph( g, boost::get( boost::vertex_name, g ) );
    }
#endif
#endif
    return;
  }

  if (strategy == own_lines)
  {
    // If just the lines owned by the element are here, use a "simpler"
    // algorithm that does not depend on geometry. This is required for
    // levelset cut sides than do not posses geometrical information.

    plain_cycle_set base_cycles;
    find_cycles(g, base_cycles);

    main_cycles_.reserve(base_cycles.size());

    for (plain_cycle_set::iterator i = base_cycles.begin(); i != base_cycles.end(); ++i)
    {
      cycle_t *c = *i;

      main_cycles_.push_back(Cycle());
      Cycle &pc = main_cycles_.back();
      pc.reserve(c->size());

      for (cycle_t::iterator i = c->begin(); i != c->end(); ++i)
      {
        vertex_t u = *i;
        pc.push_back(name_map[u]);
      }

      delete c;
    }
  }
  else
  {
    // Use geometry to find the right embedding and find the cycles.

    // find local coordinates

    std::map<vertex_t, LINALG::Matrix<3, 1>> local;

    vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
    {
      // prepare vars
      Point *p = name_map[*vi];
      const LINALG::Matrix<3, 1> xyz(p->X(), true);
      LINALG::Matrix<3, 1> rst;

      // get local coordinates from the side element
      side->LocalCoordinates(xyz, rst);
#if DEBUG_POINTGRAPH
      std::cout << "For point" << p->Id() << std::endl;
      std::cout << "Local coordinate on the side are" << rst << std::endl;
#endif

      std::pair<vertex_t, LINALG::Matrix<3, 1>> tmppair(*vi, rst);
      local.insert(tmppair);
    }

    // find unconnected components (main facet(s) and holes)

    std::vector<int> component(boost::num_vertices(g));

    int num_comp = boost::connected_components(g,
        boost::make_iterator_property_map(component.begin(), boost::get(boost::vertex_index, g)));

    // find cycles on each component

    if (num_comp == 1)
    {
      bool main_cycle = GEO::CUT::IMPL::FindCycles(g, cycle, local, main_cycles_);
      if (location == element_side and not main_cycle)
      {
        GnuplotDumpCycles("cycles", main_cycles_);
        boost::print_graph(g, boost::get(boost::vertex_name, g));

#if 0
        // Output the graph in DOT format
        boost::dynamic_properties dp;
        dp.property( "label", boost::get( boost::vertex_index, g ) );
        std::ofstream out( "side-graph.dot" );
        boost::write_graphviz( out, g, dp, std::string(), boost::get( boost::vertex_index, g ) );
#endif

        run_time_error("cycle needs to contain side edges");
      }
    }
    else if (num_comp > 1)
    {
      for (int i = 0; i < num_comp; ++i)
      {
        typedef boost::filtered_graph<graph_t, edge_filter> filtered_graph_t;
        edge_filter filter(g, component, i);
        filtered_graph_t fg(g, filter);

        std::vector<Cycle> filtered_cycles;

        graph_t cg;
        boost::copy_graph(fg, cg);
        bool main_cycle = GEO::CUT::IMPL::FindCycles(cg, cycle, local, filtered_cycles);

        if (main_cycle)
        {
          if (main_cycles_.size() != 0)
          {
            run_time_error("one set of main cycles only");
          }
          std::swap(main_cycles_, filtered_cycles);
        }
        else
        {
          hole_cycles_.push_back(std::vector<Cycle>());
          std::swap(hole_cycles_.back(), filtered_cycles);
        }
      }

      if (location == element_side and main_cycles_.size() == 0)
      {
        run_time_error("cycle needs to contain side edges");
      }
    }
    else
    {
      if (location == element_side) run_time_error("empty graph discovered");
    }
  }


  // filtering hole cycles and maincycles to include only internal cut_facets when creating facets
  // on the cut side
  if ((location == cut_side) and (!(strategy == own_lines)))
  {
    std::vector<Cycle> erased_cycles;
    for (std::vector<Cycle>::iterator it = main_cycles_.begin(); it != main_cycles_.end();)
    {
      const std::vector<Point *> cycle_points = (*it)();
      bool to_erase = false;
      for (std::vector<Point *>::const_iterator ip = cycle_points.begin(); ip != cycle_points.end();
           ++ip)
      {
        if (not(*ip)->IsCut(element))
        {
          erased_cycles.push_back(*it);
          to_erase = true;
          break;
        }
      }
      if (to_erase)
        it = main_cycles_.erase(it);
      else
        ++it;
    }
  }
}

/*---------------------------------------------------------------------------------*
 * In graph, if any edge has a single point, it will be deleted
 *---------------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::FixSinglePoints(Cycle &cycle)
{
  for (;;)
  {
    bool found = false;
    for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
    {
      int p = i->first;
      plain_int_set &row = i->second;
      if (row.size() < 2)
      {
        found = true;
        for (plain_int_set::iterator i = row.begin(); i != row.end(); ++i)
        {
          int p2 = *i;
          plain_int_set &row2 = graph_[p2];
          row2.erase(p);
          if (row2.size() == 0) graph_.erase(p2);
        }
        graph_.erase(p);

        // There are degenerated cases. A very sharp triangle with one and the
        // same cut point on two edges close to the sharp node. In this case
        // the node will be dropped. The cycle will contain the cut point
        // twice. This needs to be fixed.

        cycle.DropPoint(GetPoint(p));

        break;
      }
    }
    if (not found)
    {
      return;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IMPL::PointGraph::Graph::HasSinglePoints(
    GEO::CUT::IMPL::PointGraph::Location location)
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    plain_int_set &row = i->second;
    if (row.size() < 2)
    {
      return true;
    }
  }
  return false;
}


// Check if this side has single point in the pointgraph, because other
// side was touched by the "tip" at this point
bool GEO::CUT::IMPL::PointGraph::Graph::HasTouchingEdge(Element *element, Side *side)
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    plain_int_set &row = i->second;
    if (row.size() < 2)
    {
      LINALG::Matrix<3, 1> cut_pointxyz;
      // if there is  point in the poingraph, that have no less than neighbors
      Point *cut_point = all_points_[i->first];
      cut_point->Coordinates(cut_pointxyz.A());

      for (plain_edge_set::const_iterator e = cut_point->CutEdges().begin();
           e != cut_point->CutEdges().end(); ++e)
      {
        LINALG::Matrix<3, 1> edge_vector;
        Edge *ed = *e;

        // get the vector from opposite node of the edge to cutting point
        if (cut_point->NodalPoint(ed->Nodes()))
        {
          if (ed->Nodes()[0]->point() == cut_point)
          {
            (ed->Nodes()[1]->point())->Coordinates(edge_vector.A());
          }

          else if (ed->Nodes()[1]->point() == cut_point)
          {
            (ed->Nodes()[0]->point())->Coordinates(edge_vector.A());
          }
          else
            dserror("This should not happen or not implemented");

          edge_vector.Update(-1.0, cut_pointxyz, 1.0);
          for (plain_side_set::const_iterator s = ed->Sides().begin(); s != ed->Sides().end(); ++s)
          {
            // getting side normal with respect to resp(0,0) by default local coordiantes
            LINALG::Matrix<2, 1> resp;
            LINALG::Matrix<3, 1> norm_vec;
            Side *sd = *s;
            sd->Normal(resp, norm_vec);

            for (plain_element_set::const_iterator el = sd->Elements().begin();
                 el != sd->Elements().end(); ++el)
            {
              // getting element center for this element
              LINALG::Matrix<3, 1> element_center;
              Element *elmnt = *el;
              if (elmnt->Shape() != DRT::Element::hex8)
              {
                std::cout << "==| WARNING: Element Type != hex8 not supported by check "
                             "Graph::HasTouchingEdge! |==\n"
                          << "==| WARNING: Therefore we skip this test, please implement if you "
                             "use another element type! |=="
                          << std::endl;
                continue;
              }
              elmnt->ElementCenter(element_center);
              // getting vector pointing outward the element
              LINALG::Matrix<3, 1> out_vec;
              out_vec.Update(1.0, cut_pointxyz, -1.0, element_center);
              // if normal is pointing inwards, reverse it to point outwards
              if (out_vec.Dot(norm_vec) < 0) norm_vec.Scale(-1.0);
              if (norm_vec.Dot(edge_vector) < 0)
              {
                dserror("Single point problem, one elements is going inside another");
              }
            }
          }
        }
        else
        {
#ifdef DEBUG_POINTGRAPH
          std::ofstream file("touchign_element_detectiong_failed.pos");

          GEO::CUT::OUTPUT::GmshNewSection(file, "Element");
          GEO::CUT::OUTPUT::GmshElementDump(file, element, false);
          GEO::CUT::OUTPUT::GmshEndSection(file, false);

          GEO::CUT::OUTPUT::GmshSideDump(file, side, std::string("Side"));
          GEO::CUT::OUTPUT::GmshEdgeDump(file, ed, std::string("EdgeContainingPoint"));
          GEO::CUT::OUTPUT::GmshPointDump(
              file, cut_point, cut_point->Id(), std::string("CutPoint"), false, NULL);


          if (row.size() == 1)
          {
            Point *next_point = all_points_[row[0]];
            file << "//Next point has Id" << next_point->Id() << std::endl;
            cut_point->DumpConnectivityInfo();
            next_point->DumpConnectivityInfo();
          }
          else
            file << "//This point is not connected to anything\n";

          /// dump all the points of the pointgraph
          std::ofstream file_pgraph("pointgraph_dump.pos");
          for (std::map<int, Point *>::iterator it = all_points_.begin(); it != all_points_.end();
               ++it)
          {
            std::stringstream point_section_name;
            point_section_name << "Point" << (it->second)->Id();
            GEO::CUT::OUTPUT::GmshNewSection(file_pgraph, point_section_name.str());
            GEO::CUT::OUTPUT::GmshPointDump(
                file_pgraph, (it->second), (it->second)->Id(), false, NULL);
            GEO::CUT::OUTPUT::GmshEndSection(file_pgraph, false);
            (it->second)->DumpConnectivityInfo();
          }

          file_pgraph.close();
          file.close();
#endif

          std::stringstream err_msg;
          err_msg << "The single cut point in pointgraph(Id=" << cut_point->Id() << ")"
                  << " is not a nodal point of any of the edges connected to it (Not Touching)\n\
            This can for instance happen if your cut surface is not closed, so check your geometry first!\n";
          throw std::runtime_error(err_msg.str());
        }
      }
    }
  }
  return true;
}

bool GEO::CUT::IMPL::PointGraph::Graph::SimplifyConnections(Element *element, Side *side)
{
  for (std::map<int, plain_int_set>::iterator i = graph_.begin(); i != graph_.end(); ++i)
  {
    plain_int_set &row = i->second;
    if (row.size() < 2)
    {
      if (row.size() == 1)
      {
        Point *single = all_points_[i->first];
        Point *other = all_points_[row.front()];
        if (single->NodalPoint(side->Nodes()))
        {
          // get touching edges of nodal point
          const std::vector<Edge *> side_edges = side->Edges();
          std::vector<Edge *> point_side_edges;
          for (std::vector<Edge *>::const_iterator it = side_edges.begin(); it != side_edges.end();
               ++it)
          {
            if (single->NodalPoint((*it)->Nodes())) point_side_edges.push_back(*it);
          }
          std::vector<Edge *>::iterator it = point_side_edges.begin();
          // we are fine if single point touches all touching edges (on this side) of the cut point
          for (; it != point_side_edges.end(); ++it)
          {
            if (not other->IsCut(*it)) break;
          }

          if (it != point_side_edges.end())
            return false;
          else
            return true;
        }
        else
          return false;
      }
      else
        dserror("Point in pointgraph is not connected to anything. Look into it!");
    }
  }
  return false;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::PointGraph::Graph::GnuplotDumpCycles(
    const std::string &filename, const std::vector<Cycle> &cycles)
{
  int counter = 0;
  for (std::vector<Cycle>::const_iterator i = cycles.begin(); i != cycles.end(); ++i)
  {
    const Cycle &points = *i;

    std::stringstream str;
    str << filename << counter << ".plot";
    std::cout << str.str() << "\n";
    std::ofstream file(str.str().c_str());
    points.GnuplotDump(file);

    counter += 1;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Point *GEO::CUT::IMPL::PointGraph::Graph::GetPoint(int i)
{
  std::map<int, Point *>::iterator j = all_points_.find(i);
  if (j != all_points_.end()) return j->second;
  return NULL;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::IMPL::PointGraph *GEO::CUT::IMPL::PointGraph::Create(Mesh &mesh, Element *element,
    Side *side, PointGraph::Location location, PointGraph::Strategy strategy)
{
  PointGraph *pg = NULL;
  const unsigned dim = element->Dim();
  switch (dim)
  {
    case 1:
      pg = new SimplePointGraph_1D(mesh, element, side, location, strategy);
      break;
    case 2:
      pg = new SimplePointGraph_2D(mesh, element, side, location, strategy);
      break;
    case 3:
      pg = new PointGraph(mesh, element, side, location, strategy);
      break;
    default:
      dserror("Unsupported element dimension! ( dim = %d )", dim);
      break;
  }
  return pg;
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::IMPL::PointGraph::Graph> GEO::CUT::IMPL::PointGraph::CreateGraph(
    unsigned dim)
{
  switch (dim)
  {
    case 1:
      return Teuchos::rcp(new SimplePointGraph_1D::Graph());
    case 2:
      return Teuchos::rcp(new SimplePointGraph_2D::Graph());
    case 3:
      return Teuchos::rcp(new PointGraph::Graph());
    default:
      dserror("Unsupported element dimension!");
      exit(EXIT_FAILURE);
  }
}
