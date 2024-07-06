/*----------------------------------------------------------------------*/
/*! \file

\brief class representing a geometrical edge

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_cut_facet.hpp"
#include "4C_cut_intersection.hpp"
#include "4C_cut_levelsetside.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_point_impl.hpp"
#include "4C_cut_position.hpp"
#include "4C_global_data.hpp"

#include <stack>
#include <string>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Edge> Core::Geo::Cut::Edge::create(
    Core::FE::CellType edgetype, const std::vector<Node*>& nodes)
{
  EdgeFactory factory;
  return factory.create_edge(edgetype, nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Edge> Core::Geo::Cut::Edge::create(
    unsigned shardskey, const std::vector<Node*>& nodes)
{
  return Edge::create(Core::Elements::ShardsKeyToDisType(shardskey), nodes);
}

bool Core::Geo::Cut::Edge::find_cut_points_level_set(
    Mesh& mesh, Element* element, Side& side, Side& other)
{
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->is_cut(&other))  // cut points may already obtained when some other side was considered
    {
      p->add_element(element);
    }
  }

  // test for the cut of edge and side first time in this mesh

  PointSet cut_points;
  other.cut(mesh, *this, cut_points);
  for (PointSet::iterator i = cut_points.begin(); i != cut_points.end(); ++i)
  {
    Point* p = *i;

    p->add_edge(this);
    p->add_element(element);

    // These adds are implicitly done, but for documentation do them all explicitly.

    p->add_side(&side);
    p->add_side(&other);
    add_point(p);
  }

  return cut_points.size() > 0;
}

bool Core::Geo::Cut::Edge::find_cut_points(Mesh& mesh, Element* element, Side& side, Side& other)
{
  // dispatch function call between side = LevelSetSide and  normal side
  return other.find_cut_points_dispatch(mesh, element, side, *this);
}


/*-----------------------------------------------------------------------------*
 *  Find points at which this edge which is in "side" cuts the "other"
 *-----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Edge::find_cut_points_mesh_cut(
    Mesh& mesh, Element* element, Side& side, Side& other, PointSet* cutpoints)
{
  bool cut = false;
  bool parallel_cut = false;
  PointSet cut_points;
  std::vector<Point*> point_stack;

  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    // cut points may already obtained when some other side was considered
    if (p->is_cut(&other) and p->is_cut(&other, this) and p->is_cut(this))
    {
      point_stack.push_back(p);
    }
  }
  // everything is handled

  // this should not be possible
  if (point_stack.size() > 2)
  {
    Core::Geo::Cut::Output::DebugDump_MoreThanTwoIntersectionPoints(this, &other, point_stack);
    FOUR_C_THROW("Line x Side has more than 2 intersection points.Namely %u", point_stack.size());
  }

  // if 2 points found - we are done
  else if (point_stack.size() == 2)
  {
    point_stack[0]->add_element(element);
    point_stack[1]->add_element(element);
    cut = true;
  }
  // might be possible to have another intersection point (if edge is parallel to the size)
  else if (point_stack.size() == 1)
  {
    Point* p = point_stack[0];
    if ((begin_node()->point() != p) and (end_node()->point() != p))
    {
      // try to compute assuming parallel cut ( no real intersection )
      parallel_cut = other.just_parallel_cut(mesh, *this, cut_points);
      // if edge is not parallel, nothing to do here, one point already intersected and no more can
      // be found for normal cases, otherwise we need to add some extra connection
      if (not parallel_cut)
      {
        cut = true;
        p->add_element(element);
      }
    }
    // else we need to check for paralelism
    else
    {
      /* \   /
       |  \ /
       |   x
       |  /  \  ,
     --x-x    \
       |       \   ,
       |___ _ __\ */
      // one point can be from another edge that is connected to it and be on side. However there is
      // possibility of another not yet found intersection point, if this edge is parallel to the
      // side
      int id_skip = (begin_node()->point() == p) ? 0 : 1;
      parallel_cut = other.just_parallel_cut(mesh, *this, cut_points, id_skip);
      if (not parallel_cut) cut = true;
    }
  }
  else
  {
    // not yet cut at all, so do full cut
    other.cut(mesh, *this, cut_points);
  }

  if ((not parallel_cut) and cut)
  {
    if (cutpoints)
    {
      *cutpoints = cut_points;
      std::copy(
          point_stack.begin(), point_stack.end(), std::inserter(*cutpoints, cutpoints->end()));
    }
    return true;
  }

  for (PointSet::iterator i = cut_points.begin(); i != cut_points.end(); ++i)
  {
    Point* p = *i;

    p->add_edge(this);
    p->add_element(element);

    // These adds are implicitly done, but for documentation do them all explicitly.

    p->add_side(&side);
    p->add_side(&other);

    p->add_pair(&other, this);

    add_point(p);
  }

  if (cutpoints)
  {
    *cutpoints = cut_points;
    std::copy(point_stack.begin(), point_stack.end(), std::inserter(*cutpoints, cutpoints->end()));
  }
  return cut_points.size() > 0;
}

/*----------------------------------------------------------------------------*
 * Cut points falling on this edge that are common to the two given sides are
 * extracted
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::get_cut_points(Element* element, Side& side, Side& other, PointSet& cuts)
{
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->is_cut(&other) and p->is_cut(element) and p->is_cut(&other, this))
    {
      cuts.insert(p);
    }
  }
}

/*----------------------------------------------------------------------------*
 * Cut points falling on this edge that are common to the given edge are
 * extracted
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::get_cut_points(Edge* other, PointSet& cuts)
{
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->is_cut(other))
    {
      cuts.insert(p);
    }
  }
}

/// Find common points (excluding cut_points points) between two edges
void Core::Geo::Cut::Edge::common_nodal_points(Edge* edge, std::vector<Point*>& common)
{
  const std::vector<Node*>& other_nodes = edge->nodes();
  const std::vector<Node*>& my_nodes = nodes();
  for (std::vector<Node*>::const_iterator it = other_nodes.begin(); it != other_nodes.end(); ++it)
  {
    for (std::vector<Node*>::const_iterator jt = my_nodes.begin(); jt != my_nodes.end(); ++jt)
    {
      if ((*it) == (*jt)) common.push_back((*it)->point());
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::add_point(Point* cut_point)
{
  // make sure the position of the point on this edge is known
  cut_point->t(this);

  cut_points_.insert(cut_point);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::cut_point(
    Node* edge_start, Node* edge_end, std::vector<Point*>& edge_points)
{
  Point* bp = begin_node()->point();
  Point* ep = end_node()->point();
  if (bp == edge_start->point() and ep == edge_end->point())
  {
    edge_points.clear();
    edge_points.assign(cut_points_.begin(), cut_points_.end());
  }
  else if (ep == edge_start->point() and bp == edge_end->point())
  {
    edge_points.clear();
    edge_points.assign(cut_points_.rbegin(), cut_points_.rend());
  }
  else
  {
    FOUR_C_THROW("not a valid edge");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::cut_points(Side* side, PointSet& cut_points)
{
  Impl::SideCutFilter filter(side);
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    plain_line_set cut_lines;
    p->cut_lines(filter, cut_lines);
    if (cut_lines.size() > 0)
    {
      cut_points.insert(p);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::cut_points_between(Point* begin, Point* end, std::vector<Point*>& line)
{
  //   PointPositionLess::iterator bi = cut_points_.find( begin );
  //   PointPositionLess::iterator ei = cut_points_.find( end );
  PointPositionSet::iterator bi =
      std::lower_bound(cut_points_.begin(), cut_points_.end(), begin, PointPositionLess(this));
  PointPositionSet::iterator ei =
      std::lower_bound(cut_points_.begin(), cut_points_.end(), end, PointPositionLess(this));

  if (*bi != begin) bi = cut_points_.end();

  if (*ei != end) ei = cut_points_.end();

  double bt = begin->t(this);
  double et = end->t(this);

  if (bt < et)
  {
    if (bi != cut_points_.end())
    {
      ++bi;
    }
    // std::copy( bi, ei, std::back_inserter( line ) );
    line.insert(line.end(), bi, ei);
  }
  else if (bt > et)
  {
    if (ei != cut_points_.end())
    {
      ++ei;
    }
    // std::copy( ei, bi, std::back_inserter( line ) );
    line.insert(line.end(), ei, bi);
  }
  else
  {
    if (begin != end)
    {
      FOUR_C_THROW("different points at the same place");
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::cut_points_including(Point* begin, Point* end, std::vector<Point*>& line)
{
  PointPositionSet::iterator bi =
      std::lower_bound(cut_points_.begin(), cut_points_.end(), begin, PointPositionLess(this));
  PointPositionSet::iterator ei =
      std::lower_bound(cut_points_.begin(), cut_points_.end(), end, PointPositionLess(this));

  double bt = begin->t(this);
  double et = end->t(this);

  PointPositionSet::iterator bi_or = bi;
  PointPositionSet::iterator et_or = ei;

  if (bt < et)
  {
    ++ei;
    // std::copy( bi, ei, std::back_inserter( line ) );
    line.insert(line.end(), bi, ei);
  }
  else if (bt > et)
  {
    ++bi;
    // std::copy( ei, bi, std::inserter( line, line.begin() ) );
    line.insert(line.begin(), ei, bi);
  }
  else
  {
    if (begin != end)
    {
      FOUR_C_THROW("different points at the same place");
    }
  }

  // to detect if cut point local coordinates is outside of edges endpoints
  if (*bi_or != begin) FOUR_C_THROW("beginning point not on edge");
  if (*et_or != end) FOUR_C_THROW("end point not on edge");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::cut_points_inside(Element* element, std::vector<Point*>& line)
{
  Point* first = nullptr;
  Point* last = nullptr;
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->is_cut(element))
    {
      if (first == nullptr)
      {
        first = last = p;
      }
      else
      {
        last = p;
      }
    }
  }
  if (first != nullptr and first != last)
  {
    cut_points_including(first, last, line);
    if (line.size() > 2)
      for (std::vector<Point*>::iterator it = line.begin(); it != line.end(); /* */)
      {
        if (not(*it)->is_cut(element))
        {
#if EXTENDED_CUT_DEBUG_OUTPUT
          std::cout << "NOTICE: Got non-belonging point" << (*it)->Id()
                    << " for the line\n"
                       "From the element"
                    << element->Id() << "\n Erasing it" << std::endl;
#endif

          it = line.erase(it);
        }
        else
          ++it;
      }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Edge::is_cut(Side* side)
{
  // cutpoints contains end-points and internal cut-points
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->is_cut(side))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::Point* Core::Geo::Cut::Edge::node_in_element(Element* element, Point* other)
{
  Point* p = begin_node()->point();
  if (p != other and p->is_cut(element))
  {
    return p;
  }
  p = end_node()->point();
  if (p != other and p->is_cut(element))
  {
    return p;
  }
  return nullptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::rectify_cut_numerics()
{
  if (cut_points_.size() > 2)
  {
    // Rectify numerical problems that occasionally occur. There might be middle
    // points that do not know their position on an element side.
    //
    // This assumes linear sides. There might be a problem with quad4
    // sides. Those are actually not supported.

    std::map<Side*, PointPositionSet::iterator> sidecuts;

    for (PointPositionSet::iterator pi = cut_points_.begin(); pi != cut_points_.end(); ++pi)
    {
      Point* p = *pi;
      const plain_side_set& cutsides = p->cut_sides();
      for (plain_side_set::const_iterator i = cutsides.begin(); i != cutsides.end(); ++i)
      {
        Side* s = *i;
        std::map<Side*, PointPositionSet::iterator>::iterator j = sidecuts.find(s);
        if (j != sidecuts.end())
        {
          // if ( std::distance( j->second, pi ) > 1 )
          {
            PointPositionSet::iterator next = j->second;
            next++;

            for (PointPositionSet::const_iterator i = next; i != pi; ++i)
            {
              Point* p = *i;
              p->add_side(s);
            }
          }
          j->second = pi;
        }
        else
        {
          sidecuts[s] = pi;
        }
      }
    }
  }
}

/*------------------------------------------------------------------------*
 *  Gives this edge a selfcutposition and spreads the positional
 *  information                                                 wirtz 05/13
 *------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::self_cut_position(Point::PointPosition pos)
{
  FOUR_C_ASSERT(IsCutPositionUnchanged(selfcutposition_, pos),
      "Are you sure that you want to change the edge-position from inside to outside or vice "
      "versa?");

  if (selfcutposition_ == Point::undecided)
  {
    if (selfcutposition_ != pos)
    {
      selfcutposition_ = pos;
      if (pos == Point::outside or pos == Point::inside)
      {
        for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
        {
          Node* n = *i;
          if (n->self_cut_position() == Point::undecided)
          {
            n->self_cut_position(pos);
          }
        }
        for (plain_side_set::iterator i = sides_.begin(); i != sides_.end(); ++i)
        {
          Side* s = *i;
          if (s->id() >= 0) s->get_self_cut_position(pos);
        }
      }
    }
  }
}

/*------------------------------------------------------------------------*
 *  Changes the selfcutposition of this edge and spreads the positional
 *  information                                                 wirtz 07/16
 *------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::change_self_cut_position(Point::PointPosition pos)
{
  if (selfcutposition_ != pos)
  {
    selfcutposition_ = pos;
    for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
    {
      Node* n = *i;
      n->change_self_cut_position(pos);
    }
    for (plain_side_set::iterator i = sides_.begin(); i != sides_.end(); ++i)
    {
      Side* s = *i;
      s->change_self_cut_position(pos);
    }
  }
}

/*------------------------------------------------------------------------*
 *  Replaces the node "nod" of the edge with given node "replwith"
 *                                                              sudhakar 09/13
 *------------------------------------------------------------------------*/
void Core::Geo::Cut::Edge::replace_node(Node* nod, Node* replwith)
{
  for (unsigned i = 0; i < nodes_.size(); i++)
  {
    Node* orig = nodes_[i];

    if (orig->id() == nod->id())
    {
      nodes_[i] = replwith;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned prob_dim, Core::FE::CellType edge_type, unsigned dim_edge,
    unsigned num_nodes_edge>
bool Core::Geo::Cut::ConcreteEdge<prob_dim, edge_type, dim_edge, num_nodes_edge>::cut(
    Mesh& mesh, Side& side, PointSet& cuts)
{
  Teuchos::RCP<Core::Geo::Cut::IntersectionBase> inter_ptr = intersection_ptr(side.shape());

  inter_ptr->init(&mesh, this, &side, false, false, true);
  return inter_ptr->intersect(cuts);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned prob_dim, Core::FE::CellType edge_type, unsigned dim_edge,
    unsigned num_nodes_edge>
bool Core::Geo::Cut::ConcreteEdge<prob_dim, edge_type, dim_edge, num_nodes_edge>::just_parallel_cut(
    Mesh& mesh, Side& side, PointSet& cuts, int skip_id)
{
  Teuchos::RCP<Core::Geo::Cut::IntersectionBase> inter_ptr = intersection_ptr(side.shape());

  inter_ptr->init(&mesh, this, &side, false, false, false);
  return (inter_ptr->handle_parallel_intersection(cuts, skip_id) > 0);
}

template <unsigned prob_dim, Core::FE::CellType edge_type, unsigned dim_edge,
    unsigned num_nodes_edge>
bool Core::Geo::Cut::ConcreteEdge<prob_dim, edge_type, dim_edge,
    num_nodes_edge>::handle_parallel_cut(Edge* other, Side* side, PointSet* cut_points,
    Core::Geo::Cut::CutFloatType floattype)
{
  PointSet parallel_cuts;
  if (num_nodes_edge != 2 or other->num_nodes() != 2)
    FOUR_C_THROW("This method works only for 2 line2 intersections!");

  std::vector<Node*> touch_nodes;
  get_touching_points(other->nodes(), touch_nodes, floattype);

  unsigned int nodes_e1 = touch_nodes.size();
  other->get_touching_points(nodes(), touch_nodes, floattype);

  // nodes that belong to e1
  for (unsigned int i = 0; i < nodes_e1; ++i)
  {
    Point* cut_point = Point::insert_cut(this, side, touch_nodes[i]);
    std::stringstream msg;
#if CUT_CREATION_INFO
    msg << "//Added because " << cut_point->Id() << " lies on the edge" << this->Id() << std::endl;
    msg << "//Added during parallel edges checking" << std::endl;
#endif
    cut_point->add_edge_intersection(this, other, side, this, msg.str());
    if (cut_points) cut_points->insert(cut_point);

      // if a nodal point inserted on the other edge, we also want to add all sides of that node
#if CUT_CREATION_INFO
    msg << "//Added because of intersection of the edges of the sides of this nodal point with the "
           "edge "
        << this->Id() << std::endl;
#endif
    const plain_edge_set& edges = touch_nodes[i]->edges();
    for (plain_edge_set::const_iterator it = edges.begin(); it != edges.end(); ++it)
    {
      Edge* e = *it;
      cut_point->add_edge_intersection(this, e, side, this, msg.str());
    }
  }

  // nodes that belong to e2
  for (unsigned int i = nodes_e1; i < touch_nodes.size(); ++i)
  {
    Point* cut_point = Point::insert_cut(other, side, touch_nodes[i]);
    std::stringstream msg;
#if CUT_CREATION_INFO
    msg << "//Added because " << cut_point->Id() << " lies on the edge" << other->Id() << std::endl;
    msg << "//Added during parallel edges checking" << std::endl;
#endif
    cut_point->add_edge_intersection(this, other, side, this, msg.str());
    if (cut_points) cut_points->insert(cut_point);
      // if a nodal point inserted on the other edge, we also want to add all sides of that node
#if CUT_CREATION_INFO
    msg << "//Added because of intersection of the edges of the sides of this nodal point with the "
           "edge "
        << other->Id() << std::endl;
#endif
    const plain_edge_set& edges = touch_nodes[i]->edges();
    for (plain_edge_set::const_iterator it = edges.begin(); it != edges.end(); ++it)
    {
      Edge* e = *it;
      cut_point->add_edge_intersection(this, e, side, this, msg.str());
    }
  }

  switch (touch_nodes.size())
  {
    case 1:
    case 2:
      return true;
    case 0:
      return false;
    default:
      FOUR_C_THROW("Something went wrong!");
  }
  return false;
}

template <unsigned prob_dim, Core::FE::CellType edge_type, unsigned dim_edge,
    unsigned num_nodes_edge>
void Core::Geo::Cut::ConcreteEdge<prob_dim, edge_type, dim_edge,
    num_nodes_edge>::get_touching_points(const std::vector<Node*>& nodes,
    std::vector<Node*>& touch_nodes, Core::Geo::Cut::CutFloatType floattype)
{
  bool signeddistance = false;
  double distance = 0;
  Core::LinAlg::Matrix<prob_dim, 1> xsi;
  Core::LinAlg::Matrix<prob_dim, num_nodes_edge> xyze_edge;
  coordinates(xyze_edge.data());

  const std::vector<Node*> edge_nodes = Edge::nodes();

  for (std::vector<Node*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
  {
    xsi = 0.0;
    distance = 0;

    Node* node = *it;
    Point* p = node->point();
    Core::LinAlg::Matrix<prob_dim, 1> p_coord;
    p->coordinates(p_coord.data());

    bool conv = false;
    Kernel::PointOnSurfaceLoc loc;
    bool SurfaceWithinLimits0 = false;
    switch (floattype)
    {
      case Core::Geo::Cut::floattype_double:
      {
        Kernel::ComputeDistance<prob_dim, edge_type, false> cd(xsi);
        conv = cd(xyze_edge, p_coord, distance, signeddistance);
        loc = cd.get_side_location();
        if (conv && loc.within_side() && loc.on_side())
          SurfaceWithinLimits0 = cd.surface_within_limits(0.0);
        break;
      }
      case Core::Geo::Cut::floattype_cln:
      {
        Kernel::ComputeDistance<prob_dim, edge_type, true> cd(xsi);
        conv = cd(xyze_edge, p_coord, distance, signeddistance);
        loc = cd.get_side_location();
        if (conv && loc.within_side() && loc.on_side())
          SurfaceWithinLimits0 = cd.surface_within_limits(0.0);
        break;
      }
      default:
      {
        FOUR_C_THROW("Unexpected floattype for Kernel::compute_distance!");
      }
    }

    if (!conv) FOUR_C_THROW("Newton did not converge for parallelism detection of lines!");

    if (SurfaceWithinLimits0) touch_nodes.push_back(node);

    // Extra check of the connection to the nodal points of the edge
    for (std::vector<Node*>::const_iterator jt = edge_nodes.begin(); jt != edge_nodes.end(); ++jt)
    {
      Node* edge_node = *jt;
      Point* edge_point = edge_node->point();
      if (Core::Geo::Cut::DistanceBetweenPoints(p, edge_point) < SIDE_DETECTION_TOLERANCE)
      {
        p->dump_connectivity_info();
        edge_point->dump_connectivity_info();
        std::stringstream err_msg;
        err_msg << "Distance between points is " << std::setprecision(15)
                << Core::Geo::Cut::DistanceBetweenPoints(p, edge_point)
                << " This two points should have been merged!";
        FOUR_C_THROW(err_msg.str());
      }
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned prob_dim, Core::FE::CellType edge_type, unsigned dim_edge,
    unsigned num_nodes_edge>
bool Core::Geo::Cut::ConcreteEdge<prob_dim, edge_type, dim_edge, num_nodes_edge>::compute_cut(
    Mesh* mesh, Edge* other, Side* side, PointSet* cut_points, double& tolerance)
{
#if CUT_CREATION_INFO
  static unsigned long int edge_inter_id = 0;
  edge_inter_id++;
#endif

  // nodal edge coordinates
  Core::LinAlg::Matrix<prob_dim, num_nodes_edge> xyze_this;
  coordinates(xyze_this.data());

  // nodal pseudo side coordinates
  /* Note: We insert the this edge as side. In this way we can use the
   * standard intersection object. */
  Core::LinAlg::SerialDenseMatrix xyze_other;
  xyze_other.shape(other->n_prob_dim(), other->num_nodes());
  other->coordinates(xyze_other);

  Teuchos::RCP<Core::Geo::Cut::IntersectionBase> inter_ptr = intersection_ptr(other->shape());
  // other is line element, this is surface
  inter_ptr->init(xyze_other, xyze_this, false, false, false, &(mesh->get_options()));

  bool edges_parallel =
      handle_parallel_cut(other, side, cut_points, mesh->get_options().geom_distance_floattype());

  if (not edges_parallel)
  {
    // perform the actual edge - edge intersection
    const enum IntersectionStatus istatus =
        inter_ptr->compute_edge_side_intersection(tolerance, true, nullptr);

    enum IntersectionStatus retstatus = istatus;

    switch (istatus)
    {
      case intersect_single_cut_point:
      {
        double* x_ptr = nullptr;

        x_ptr = inter_ptr->final_point_edge_edge();

        const double& pos = inter_ptr->local_coordinates()[0];

        Point* p = nullptr;
        p = Point::new_point(*mesh, x_ptr, pos, this, side, tolerance);

        p->add_edge(this);
        // NOTE: Might cause errors in the future, since local coordinate of the point on other edge
        // is based on the merged point, not original one
        p->add_edge(other);
        // might remove this one, as there is already one in cut_edge.cpp::find_cut_points cln
        p->add_pair(side, this);

        std::stringstream msg;
        std::stringstream msg_extra;

#if CUT_CREATION_INFO

        msg << "// Added because of REAL edge-edge intersection of " << this->Id() << " and "
            << other->Id() << std::endl;
        if (side) msg << "// Other side Id is " << side->Id();
        msg << "// intersection id is " << edge_inter_id;
        msg_extra << "// Added because of REAL edge-edge intersection of " << this->Id() << " and "
                  << other->Id() << std::endl;
        msg_extra << " // And this follows";
        if (side) msg_extra << "// Other side Id is " << side->Id();
        msg_extra << "// intersection id is " << edge_inter_id;
        msg << "// Point before merging " << x_ptr[0] << " " << x_ptr[1] << " " << x_ptr[2];
        msg_extra << "// Point before merging " << x_ptr[0] << " " << x_ptr[1] << " " << x_ptr[2];
        p->AddCreationInfo(std::make_pair(side, this), msg.str());
#endif

        if (cut_points) cut_points->insert(p);

        p->add_edge_intersection(this, other, side, this, msg_extra.str());

        break;
      }
      case intersect_multiple_cut_points:
      {
        // happens when one edge lies on another

        std::vector<Core::LinAlg::Matrix<prob_dim, 1>> xyz_cuts;
        std::vector<Core::LinAlg::Matrix<dim_edge, 1>> r_cuts;
        inter_ptr->final_points(xyz_cuts);
        inter_ptr->local_side_coordinates(r_cuts);

        FOUR_C_ASSERT(xyz_cuts.size() == r_cuts.size(), "Size mismatch!");

        Point* p = nullptr;
        for (unsigned i = 0; i < xyz_cuts.size(); ++i)
        {
          p = Point::new_point(*mesh, xyz_cuts[i].data(), r_cuts[i](0), this, side, tolerance);

          p->add_edge(this);
          p->add_edge(other);
          // might remove this one, as there is already one in cut_edge.cpp::find_cut_points cln
          p->add_pair(side, this);
          p->add_edge_intersection(this, other, side, this, "MULTIPLE CUT POINTS");

          if (cut_points) cut_points->insert(p);

#if CUT_CREATION_INFO
          p->AddCreationInfo(std::make_pair(side, this),
              "// Added because of edge-edge intersection. But multiple cut points");
#endif
        }

        break;
      }
      case intersect_no_cut_point:
      case intersect_newton_failed:
      {
        // if newton fails for edge-edge intersetion, it means least-square failed => no
        // intersection point
        retstatus = intersect_no_cut_point;
        /* do nothing */
        break;
      }
      default:
        FOUR_C_THROW("Unsupported intersection status! ( status = %d )", istatus);
        exit(EXIT_FAILURE);
    }
    return (static_cast<bool>(retstatus));
  }
  else
    return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned prob_dim, Core::FE::CellType edge_type, unsigned dim_edge,
    unsigned num_nodes_edge>
Teuchos::RCP<Core::Geo::Cut::IntersectionBase>
Core::Geo::Cut::ConcreteEdge<prob_dim, edge_type, dim_edge, num_nodes_edge>::intersection_ptr(
    const Core::FE::CellType& sidetype) const
{
  return Core::Geo::Cut::IntersectionBase::create(edge_type, sidetype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::Cut::Edge> Core::Geo::Cut::EdgeFactory::create_edge(
    const Core::FE::CellType& edgetype, const std::vector<Node*>& nodes) const
{
  Teuchos::RCP<Edge> cedge_ptr = Teuchos::null;
  const int probdim = Global::Problem::instance()->n_dim();
  switch (edgetype)
  {
    case Core::FE::CellType::line2:
    {
      cedge_ptr = Teuchos::rcp(create_concrete_edge<Core::FE::CellType::line2>(nodes, probdim));
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported edge type! ( %d | %s )", edgetype,
          Core::FE::CellTypeToString(edgetype).c_str());
      break;
    }
  }
  return cedge_ptr;
}

template class Core::Geo::Cut::ConcreteEdge<2, Core::FE::CellType::line2>;
template class Core::Geo::Cut::ConcreteEdge<3, Core::FE::CellType::line2>;

FOUR_C_NAMESPACE_CLOSE
