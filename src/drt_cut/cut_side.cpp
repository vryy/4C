/*----------------------------------------------------------------------*/
/*! \file
 \brief class representing a geometrical side

\level 2

\maintainer Martin Kronbichler
 *----------------------------------------------------------------------*/

#include "cut_clnwrapper.H"
#include "cut_side.H"
#include "cut_position.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"
#include "cut_pointgraph.H"
#include "cut_levelsetside.H"

#include "../drt_lib/drt_globalproblem.H"


#include <string>
#include <stack>
#include <list>

//#define DEBUG_PARALLEL_CUT_SURFACE


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side::Side(int sid, const std::vector<Node*>& nodes, const std::vector<Edge*>& edges)
    : sid_(sid), nodes_(nodes), edges_(edges)
{
  for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end(); ++i)
  {
    Edge* e = *i;
    e->Register(this);
  }
  selfcutposition_ = Point::undecided;

  if (sid > -1)
  {
    boundingvolume_ = Teuchos::rcp(BoundingBox::Create(*this));
  }
}

/*----------------------------------------------------------------------------*
      Returns the edge of this side with given begin and end points
 *----------------------------------------------------------------------------*/
GEO::CUT::Edge* GEO::CUT::Side::FindEdge(Point* begin, Point* end)
{
  for (std::vector<Edge*>::iterator i = edges_.begin(); i != edges_.end(); ++i)
  {
    Edge* e = *i;
    if (e->Matches(begin, end))
    {
      return e;
    }
  }
  return NULL;
}

bool GEO::CUT::Side::FindCutPointsDispatch(Mesh& mesh, Element* element, Side& side, Edge& e)
{
  return e.FindCutPointsMeshCut(mesh, element, side, *this);
}

/*----------------------------------------------------------------------------*
 * Calculate the points at which the other side intersects with this considered
 * side
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindCutPoints(Mesh& mesh, Element* element, Side& other)
{
  bool cut = false;
  const std::vector<Edge*>& edges = Edges();
  for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end(); ++i)
  {
    Edge* e = *i;
    if (e->FindCutPoints(mesh, element, *this, other))
    {
      cut = true;
    }
  }
  return cut;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindCutLines(Mesh& mesh, Element* element, Side& other)
{
  /* check whether cut lines are already created (still need to create lines new
   * in case we create AmbiguousCutLines!) */
  for (std::vector<Line*>::iterator i = cut_lines_.begin(); i != cut_lines_.end(); ++i)
  {
    Line* l = *i;
    if (l->IsCut(this, &other))
    {
      if (not l->IsCut(element))
      {
        dserror(
            "Line (%d, %d) is cut by both sides but not by the element, check this "
            "situation as it is not expected!",
            l->BeginPoint()->Id(), l->EndPoint()->Id());
        // l->AddElement( element );
      }
      other.AddLine(l);
    }
  }

  // creating cut lines between the given sides for the first time
  PointSet cuts;
  GetCutPoints(element, other, cuts);

  switch (cuts.size())
  {
    // ------------------------------------------------------------------------
    // no point --> no line!
    // ------------------------------------------------------------------------
    case 0:
    {
      return false;
    }
    // ------------------------------------------------------------------------
    // just touching point, there shouldn't be any line!
    // ------------------------------------------------------------------------
    /* This is also the default case for element.Dim() == 1. But in this case
     * the construction of cut_lines_ seems unnecessary, since the CutPoints
     * should be sufficient.
     *
     * (see also IMPL::PointGraph1D::AddCutLinesToGraph)     hiermeier 11/16 */
    case 1:
    {
      if (IsTouched(other, cuts[0]))
        return false;
      else
        dserror("Single point between two sides is not a touching point!");
    }
      // ------------------------------------------------------------------------
      // The normal case. A straight cut.
      // ------------------------------------------------------------------------
    case 2:
    {
      std::vector<Point*> c;
      c.reserve(2);
      c.assign(cuts.begin(), cuts.end());
      mesh.NewLine(c[0], c[1], this, &other, element);
      return true;
    }
    default:
    {
      return other.FindAmbiguousCutLines(mesh, element, *this, cuts);
    }
  }

  dserror("How did you get here?");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::AllOnNodes(const PointSet& points)
{
  const std::vector<Node*>& nodes = Nodes();
  for (PointSet::const_iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    if (not p->NodalPoint(nodes))
    {
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::GetCutPoints(Element* element, Side& other, PointSet& cuts)
{
  // Get all Cut Points intersecting this side and an edge of the other side
  {
    const std::vector<Edge*>& edges = Edges();
    for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end(); ++i)
    {
      Edge* e = *i;
      e->GetCutPoints(element, *this, other, cuts);
    }
  }
  // Get all Cut Points intersecting the other side and an edge of this side
  {
    const std::vector<Edge*>& edges = other.Edges();
    for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end(); ++i)
    {
      Edge* e = *i;
      e->GetCutPoints(element, other, *this, cuts);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::AddPoint(Point* cut_point) { cut_points_.insert(cut_point); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::AddLine(Line* cut_line)
{
  if (std::find(cut_lines_.begin(), cut_lines_.end(), cut_line) == cut_lines_.end())
  {
    cut_lines_.push_back(cut_line);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Facet* GEO::CUT::Side::FindFacet(const std::vector<Point*>& facet_points)
{
  for (std::vector<Facet*>::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    if (f->Equals(facet_points))
    {
      return f;
    }
  }
  return NULL;
}



bool GEO::CUT::Side::IsTouched(Side& other, Point* p)
{
  // first we check if it is a nodal point of one of the surface
  if (not((p->NodalPoint(other.Nodes())) || (p->NodalPoint(Nodes()))))
  {
    // now we check if it was merged, and is no longer Nodal point, however it is still
    // topologically touchign point
    if ((IsTouchedAt(&other, p)) || (other.IsTouchedAt(this, p))) return true;

    // if not the previous case, check if the cut_point lie on the edges of both sides
    // for my side
    bool on_my_edge = false;
    for (unsigned ledge = 0; ledge < Edges().size(); ++ledge)
    {
      for (unsigned lpoint = 0; lpoint < Edges()[ledge]->CutPoints().size(); ++lpoint)
      {
        if (Edges()[ledge]->CutPoints()[lpoint]->Id() == p->Id())
        {
          on_my_edge = true;
        }
      }
    }

    // for other side
    bool on_other_edge = false;
    for (unsigned ledge = 0; ledge < other.Edges().size(); ++ledge)
    {
      for (unsigned lpoint = 0; lpoint < other.Edges()[ledge]->CutPoints().size(); ++lpoint)
      {
        if (other.Edges()[ledge]->CutPoints()[lpoint]->Id() == p->Id())
        {
          on_other_edge = true;
        }
      }
    }

    if (!on_my_edge || !on_other_edge)
    {
      std::ofstream file("touching_point_between_two_sides_is_not_on_edge.pos");
      GEO::CUT::OUTPUT::GmshSideDump(file, this, std::string("ThisSide"));
      GEO::CUT::OUTPUT::GmshSideDump(file, &other, std::string("OtherSide"));
      GEO::CUT::OUTPUT::GmshPointDump(file, p, p->Id(), std::string("CutPoint"), false, NULL);
      p->DumpConnectivityInfo();
      GEO::CUT::OUTPUT::GmshWriteSection(file, "Edge", p->CutEdges());
      file.close();
      dserror(
          "Touching point between two sides does not lie on edge. This case should be analyzed");
    }
  }
  return true;
}



bool GEO::CUT::Side::IsTouchedAt(Side* other, Point* p)
{
  // count of edges from other side touching at this point
  int count = 0;
  for (unsigned int ledge = 0; ledge < other->Edges().size(); ++ledge)
  {
    if (p->IsCut(this, other->Edges()[ledge])) count++;
  }

  switch (count)
  {
    case 0:
      return false;
    case 1:
      // if touched only by one edge it can be edge-edge intersection, so we need to check current
      // side on itersection with other side, however we test it later on
      return false;
    case 2:
      // this is touching point then
      return true;
    default:
      dserror("This case should be investigated");
      return false;
  }
}

/*----------------------------------------------------------------------------*/
// Find Cut Lines for two Cut Sides, which have more than two common cut points!
//(This happens if the cutsides are in the same plane )              ager 08/15
/*----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindTouchingCutLines(
    Mesh& mesh, Element* element, Side& side, const PointSet& cut)
{
  // More that two cut points show a touch.
  //
  // If all nodes are caught and nothing else, the cut surface has hit this
  // surface exactly. No need to cut anything. However, the surface might be
  // required for integration.
  {
    // find if this side is completely inside the other side
    const std::vector<Node*>& nodes = Nodes();
    if (cut.size() == nodes.size() and AllOnNodes(cut))
    {
      for (unsigned i = 0; i < nodes.size(); ++i)
      {
        unsigned j = (i + 1) % nodes.size();
        mesh.NewLine(nodes[i]->point(), nodes[j]->point(), this, &side, element);
      }
      return true;
    }
  }
  {
    // find if the other side is completely inside this side
    const std::vector<Node*>& nodes_o = side.Nodes();
    if (cut.size() == nodes_o.size() and side.AllOnNodes(cut))
    {
      for (unsigned i = 0; i < nodes_o.size(); ++i)
      {
        unsigned j = (i + 1) % nodes_o.size();
        mesh.NewLine(nodes_o[i]->point(), nodes_o[j]->point(), this, &side, element);
      }
      return true;
    }
  }
  return false;
}


// Try to find intersection of edges of "this" side  with side.
bool GEO::CUT::Side::FindParallelIntersection(
    Mesh& mesh, Element* element, Side& side, const PointSet& cut, point_line_set& new_lines)
{
  // number of lines that could be created
  unsigned int line_counter = 0;
  // queue of points that will be inserted as cut lines, if everything goes OK.
  // we want to create only when we know, that this is normal case, and not a special cases,
  // that needs to be handled separately
  std::vector<std::pair<Point*, Point*>> cut_lines_queue;
  PointSet c;

  // loop over the edges
  for (std::vector<Edge*>::const_iterator e_it = Edges().begin(); e_it != Edges().end(); ++e_it)
  {
    Edge* e = *e_it;
    c.clear();
    const PointPositionSet& edge_cut_points = (*e_it)->CutPoints();
    for (PointPositionSet::const_iterator p_it = edge_cut_points.begin();
         p_it != edge_cut_points.end(); ++p_it)
    {
      if (cut.find(*p_it) != cut.end()) c.insert(*p_it);  //
    }

    switch (c.size())
    {
      case 0:
        break;
      case 1:
      {
        if (not(IsTouchedAt(&side, c[0])))
        {
          // if this is not a nodal point we definitely throw an error
          if ((c[0] != e->BeginNode()->point()) and (c[0] != e->EndNode()->point()))
          {
            if (not CreateParallelCutSurface(mesh, element, side, cut))
            {
              // maybe that other side will create require number of lines;
              return false;
            }
            // else we should have generated all the parallel surface and we are done with this side
            else
            {
              return true;
            }
          }
        }
        break;
      }
      case 2:
      {
        cut_lines_queue.push_back(std::make_pair(c[0], c[1]));
        break;
      }
      default:
      {
        // this can happen when edge is parallel and nodal point of the side lies on the edge
        for (PointSet::iterator c_it = c.begin(); c_it != c.end(); ++c_it)
        {
          if ((not(*c_it)->NodalPoint(e->Nodes())) and (not(*c_it)->NodalPoint(side.Nodes())))
          {
            GEO::CUT::OUTPUT::DebugDump_ThreePointsOnEdge(this, &side, e, c[0], cut);
            dserror(
                "Uknown case of three intersection point lying on side's edge in side-side "
                "intersection, this case should be reported");
          }
        }

#if EXTENDED_CUT_DEBUG_OUTPUT
        std::cout << "NOTICE: Creating cut line on a edge with multiple ( " << c.size()
                  << " ) cut_points " << std::endl;
#endif
        // sorting points first
        PointPositionSet sorted_cut_points(c.begin(), c.end(), PointPositionLess(e));

        PointPositionSet::iterator next = sorted_cut_points.begin();
        for (PointPositionSet::iterator c_it = next; next != (--sorted_cut_points.end());
             c_it = next)
        {
          ++next;
          cut_lines_queue.push_back(std::make_pair(*c_it, *next));
        }

        break;
      }
    }
  }

  // If there were no anomalies, we add all point as cut lines.
  for (std::vector<std::pair<Point*, Point*>>::iterator it = cut_lines_queue.begin();
       it != cut_lines_queue.end(); ++it)
  {
    Point* first = it->first;
    Point* second = it->second;

    mesh.NewLine(first, second, this, &side, element);
    line_counter++;
    if (first->Id() > second->Id())
      new_lines.insert(std::make_pair(first, second));
    else
      new_lines.insert(std::make_pair(second, first));
  }

  // if number of lines equal to number of cut points the whole surface is parallel
  if (line_counter == cut.size())
  {
    return true;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*
 * Find Cut Lines for two Cut Sides specially based on a discretization,
 * which have more than two common cut points!
 *   (This happens if the cut sides are in the same plane or due to numerical
 *   tolerances!)                                                    ager 08/15
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::FindAmbiguousCutLines(
    Mesh& mesh, Element* element, Side& side, const PointSet& cut)
{
  /* More than two cut points shows a touch.
   *
   * (1) If all nodes are caught and nothing else, the cut surface has hit this
   *     surface exactly. No need to cut anything. However, the surface might be
   *     required for integration. */
  if (FindTouchingCutLines(mesh, element, side, cut))
  {
    return true;
  }

  /* (2) Trying to find parallel intersection between two surfaces */
  point_line_set new_lines;
  PointSet c;

  std::vector<bool> lies_inside(2);

  // try to find intersection of edges of this side with other side
  lies_inside[0] = this->FindParallelIntersection(mesh, element, side, cut, new_lines);

  // edges of other side with this side
  lies_inside[1] = side.FindParallelIntersection(mesh, element, *this, cut, new_lines);

  // either first side lies inside second, or second inside first, or they share parts of each other
  if (lies_inside[0] || lies_inside[1] || (new_lines.size() == cut.size()))
    return true;

  else
  {
#ifdef DEBUG_PARALLEL_CUT_SURFACE
    std::cout << "Trying to create parallel surface from both sides" << std::endl;
#endif

    // either we missed one side that is partially parallel ( all cut points are on their edges)
    // or it can happen (special case for 3 points) that cut points are split
    // between edges of two sides
    std::vector<Point*> cut_point_on_edges_1;
    std::vector<Point*> cut_point_on_edges_2;

    lies_inside[1] =
        side.CreateParallelCutSurface(mesh, element, *this, cut, &cut_point_on_edges_1);
    lies_inside[0] =
        this->CreateParallelCutSurface(mesh, element, side, cut, &cut_point_on_edges_2);

    // we create cut lines during one of CreateParallelCutSurface calls and we are done
    if (lies_inside[0] or lies_inside[1]) return true;

    // we try to create cut points from the combination
    else
    {
      PointSet collected_points(cut_point_on_edges_1.begin(), cut_point_on_edges_1.end());
      std::copy(cut_point_on_edges_2.begin(), cut_point_on_edges_2.end(),
          std::inserter(collected_points, collected_points.end()));
      // if we collected this all points with parallel cut surface function calls we handle this
      // case
      if (collected_points == cut)
      {
        if (cut.size() == 3)
        {
#ifdef DEBUG_PARALLEL_CUT_SURFACE
          std::cout << "NOTICE: Creating parallel cut surface" << std::endl;
#endif
          // note: here we connect three points in any direction since it is always closed cycles
          for (PointSet::const_iterator it = cut.begin(); it != (--cut.end()); ++it)
          {
            PointSet::const_iterator next = it;
            ++next;

            Point* p1 = *it;
            Point* p2 = *next;
            if (p1 == p2) dserror("Trying to create line between two points, which are the same!");
            mesh.NewLine(p1, p2, this, &side, element);
          }

          // connect last to first to close the cycle
          if (cut.front() == cut.back())
          {
            dserror("Trying to create line between two points, which are the same!");
          }
          else
            mesh.NewLine(cut.front(), cut.back(), this, &side, element);

          return true;
        }
        else
        {
          std::ofstream file("more_than_3_points_in_the_parallel_surface.pos");
          GEO::CUT::OUTPUT::GmshSideDump(file, this, std::string("ThisSide"));
          GEO::CUT::OUTPUT::GmshSideDump(file, &side, std::string("OtherSide"));
          for (PointSet::const_iterator it = cut.begin(); it != cut.end(); ++it)
          {
            GEO::CUT::OUTPUT::GmshPointDump(
                file, (*it), (*it)->Id(), std::string("CutPoint"), false, NULL);
          }
          file.close();
          dserror("This case is not handled yet, generate GMSH output and look into it!");
        }
      }
      GEO::CUT::OUTPUT::DebugDump_MultipleCutPointsSpecial(
          this, &side, cut, collected_points, new_lines);
      throw std::runtime_error("This case should be reported probably");
    }
  }
  // should not reach here
  return false;
}


// Create paralel cut surface based on intersection with the other side, and "cut" cut_points
bool GEO::CUT::Side::CreateParallelCutSurface(Mesh& mesh, Element* element, Side& other,
    const PointSet& cut, std::vector<Point*>* cut_point_for_lines_out)
{
  std::list<Edge*> edges_cycle;
  const std::vector<Edge*>& edges = Edges();

  Edge* first = edges[0];
  edges_cycle.push_back(first);
  Node* begin_node = first->BeginNode();
  Node* end_node = first->EndNode();

  // creating list of edges, so that every previous edges is connected with the next
  for (std::vector<Edge*>::const_iterator it = (++edges.begin()); it != edges.end(); ++it)
  {
    Edge* e = *it;

    // if items has a node, that is equal to the start of the list, we push it to the
    // front of list and update begin node
    if ((begin_node == e->BeginNode()) || (begin_node == e->EndNode()))
    {
      begin_node = (begin_node == e->BeginNode()) ? e->EndNode() : e->BeginNode();
      edges_cycle.push_front(e);
    }

    // otherwise we push it to the end of the list
    else if ((end_node == e->BeginNode()) || (end_node == e->EndNode()))
    {
      end_node = (end_node == e->EndNode()) ? e->BeginNode() : e->EndNode();
      edges_cycle.push_back(e);
    }
  }

  // (since PointPositionSet on the edge is sorted (-1,1) = (BeginPoint(), EndPoint() )
  // we add points into the cycles based on orientation of first edge begin_edge1->end_edge1
  // if end_edge1 == begin_edge2 we add point in same direction as we collected,
  // Simiarly done for all  check and possibly reverse direction of inserting for remaining edges
  std::vector<Point*> cut_points_for_lines;
  Edge* first_edge = *edges_cycle.begin();
  Point* first_edge_start = first_edge->BeginNode()->point();
  Edge* second_edge = *std::next(edges_cycle.begin());

  // establish proper beginning of the cycle
  Point* tail = first_edge_start->NodalPoint(second_edge->Nodes())
                    ? first_edge->EndNode()->point()
                    : first_edge->BeginNode()->point();
  for (std::list<Edge*>::iterator it = edges_cycle.begin(); it != edges_cycle.end(); ++it)
  {
    Edge* e = *it;
    std::vector<Point*> points;
    bool reverse = false;

    // if the end node of previous vector is end node of current vector, we reverse direction of
    // point insertion
    if (tail != e->BeginNode()->point()) reverse = true;

    const PointPositionSet& edge_points = e->CutPoints();
    for (PointPositionSet::const_iterator pt = edge_points.begin(); pt != edge_points.end(); ++pt)
    {
      Point* p = *pt;
      /* p != tail so end points, that are cut points, are added only once as heads, and no tails */
      if ((cut.find(p) != cut.end()) and (p != tail))
      {
        points.push_back(p);
      }
    }

    if (reverse)
      cut_points_for_lines.insert(cut_points_for_lines.end(), points.rbegin(), points.rend());
    else
      cut_points_for_lines.insert(cut_points_for_lines.end(), points.begin(), points.end());

    if (not reverse)
      tail = e->EndNode()->point();
    else
      tail = e->BeginNode()->point();
  }

  // if last tail is the cut point, when we  connect it to the front, which was the same point -->
  // we obtain duplicate which needs to be erased
  if ((cut.find(tail) != cut.end()) and (cut_points_for_lines.back() == tail) and
      (cut_points_for_lines.front() == tail))
    cut_points_for_lines.erase(--cut_points_for_lines.end());

  bool create_simplified_geometry = false;

#ifdef DEBUG_PARALLEL_CUT_SURFACE
  std::cout << "Points before erasing are: \n";
  for (std::vector<Point*>::iterator it = cut_points_for_lines.begin();
       it != cut_points_for_lines.end(); ++it)
  {
    std::cout << (*it)->Id() << std::endl;
  }

  std::cout << "Cut points are : \n";
  for (PointSet::const_iterator it = cut.begin(); it != cut.end(); ++it)
  {
    std::cout << (*it)->Id() << std::endl;
  }
#endif

  // if we collected more points than cut points, maybe some points
  // are collected twice because they are cut points of multiple edge
  if (cut_points_for_lines.size() > cut.size())
  {
    // map to find out which location id  on parallel surface each point has
    std::map<Point*, unsigned int> unique_points;
    typedef std::map<Point*, unsigned int> unique_points_map;

    // trying to find duplicate points and remove points between them
    for (std::vector<Point*>::iterator it = cut_points_for_lines.begin();
         it != cut_points_for_lines.end();
        /* */)
    {
      int counter = std::distance(cut_points_for_lines.begin(), it);
      Point* dublicate = *it;
      std::pair<unique_points_map::iterator, bool> inserted =
          unique_points.insert(std::make_pair(dublicate, counter));
      // we found dublicate
      if (not inserted.second)
      {
#ifdef DEBUG_PARALLEL_CUT_SURFACE
        std::cout << "Found dublicate!" << std::endl;
        std::cout << "Dublicate has id of" << dublicate->Id() << std::endl;
#endif

        bool erased = false;
        // get exisiting point
        const std::pair<Point*, unsigned int>& existing_el = *(inserted.first);
        bool current_erased = false;

        // try to find these duplicates by iterating over edges
        for (std::list<Edge*>::iterator e_it = edges_cycle.begin(); e_it != edges_cycle.end();
             ++e_it)
        {
          if (not current_erased)
          {
            // if the existing point was cut by this edge, try to find current one (duplicate), on
            // one of the next edges. then remove everything in between to avoid problems with
            // facets creation e.g. in the pointgraph
            if (existing_el.first->IsCut(*e_it))
            {
              // number of cut points in between
              int n_cuts_between;
              // next cut point between existing and duplicate
              int next_cut;

              std::list<Edge*>::iterator e_next = e_it;

              // normal ordering
              if (*e_it != edges_cycle.back())
              {
                n_cuts_between = (counter - existing_el.second - 1);
                next_cut = existing_el.second + 1;
                ++e_next;
              }
              // is cut by the last edge, and next edges is first edge of the cycle
              // Then number of points in betwen is distance to the end of cycle from duplicate id +
              // distance from the beginning to the existing point
              else
              {
                n_cuts_between = (existing_el.second + (cut_points_for_lines.size() - counter - 1));
                next_cut = (counter + 1) % cut_points_for_lines.size();
                e_next = edges_cycle.begin();
              }

              if (existing_el.first->IsCut(*e_next))
              {
                std::vector<Point*> common_points;
                // safety check that we have only nodal point in between
                (*e_next)->CommonNodalPoints(*e_it, common_points);
                if (common_points.size() != 1) dserror("This should not be possible!");

                Point* nodal = common_points[0];

                // NOTE: we only handle the case, when the point in between existing is a single
                // nodal point From what it seems we can also erase even if there is other cut point
                // in between, but only if it touches both edges. If it touches just one edge
                // something different shoudl be done...
                if ((n_cuts_between > 0) and
                    (!((n_cuts_between == 1) and (cut_points_for_lines[next_cut] == nodal))))
                  dserror("There is cut point that is close to cut_point that touches both edges");

                else
                {
                  // erase one of duplicates and everything in between from the vector ( leave one
                  // for connection to other lines)
                  if (next_cut > counter)
                  {
#ifdef DEBUG_PARALLEL_CUT_SURFACE
                    std::cout << "Erasing points from end to beginning" << std::endl;
#endif
                    cut_points_for_lines.erase(cut_points_for_lines.begin(),
                        cut_points_for_lines.begin() + existing_el.second);
                    cut_points_for_lines.erase(it, cut_points_for_lines.end());
                    // obviously we are in the end now, so we don't want to process anything else
                    it = cut_points_for_lines.end();
                  }
                  else
                  {
#ifdef DEBUG_PARALLEL_CUT_SURFACE
                    std::cout << "Erasing points from beginning to end" << std::endl;
#endif
                    // we want to be pointintg to the next point after the dublicate now
                    it = cut_points_for_lines.erase(
                        cut_points_for_lines.begin() + existing_el.second, it);
                    ++it;
                  }


#ifdef DEBUG_PARALLEL_CUT_SURFACE
                  std::cout << "New point list is " << std::endl;
                  for (std::vector<Point*>::iterator it = cut_points_for_lines.begin();
                       it != cut_points_for_lines.end(); ++it)
                  {
                    std::cout << (*it)->Id() << std::endl;
                  }
#endif

#ifdef DEBUG_PARALLEL_CUT_SURFACE
                  std::cout << "Current *it id is " << (*it)->Id() << std::endl;
#endif

                  erased = true;
                  current_erased = true;
                }
              }
            }
          }
        }
        if (not erased)
          dserror("Found duplicate cut_point on different edge, but nothing got erased!");
      }
      else
        ++it;
    }

    // if we erased duplicates and some points for simplification
    if (cut_points_for_lines.size() <= cut.size()) create_simplified_geometry = true;
  }

#ifdef DEBUG_PARALLEL_CUT_SURFACE
  std::cout << "Points after erasing are" << std::endl;
  for (std::vector<Point*>::iterator it = cut_points_for_lines.begin();
       it != cut_points_for_lines.end(); ++it)
  {
    // we found matching point
    std::cout << (*it)->Id() << std::endl;
  }
#endif

  if (cut_point_for_lines_out)
  {
    *cut_point_for_lines_out = cut_points_for_lines;
  }

  // if we captured all cut_points on the edges of the side create a parallel cut surface
  if (cut_points_for_lines.size() == cut.size() || create_simplified_geometry)
  {
    // there is different parallel cut surface already on one of the sides
    std::set<Point*> new_surface(cut_points_for_lines.begin(), cut_points_for_lines.end());
    if (HasMixedParallelCutSurface(new_surface) || other.HasMixedParallelCutSurface(new_surface))
    {
      SimplifyMixedParallelCutSurface(mesh, element, other, new_surface, cut_points_for_lines);
    }

    this->parallel_cut_surfaces_.insert(new_surface);
    other.parallel_cut_surfaces_.insert(new_surface);
#ifdef DEBUG_PARALLEL_CUT_SURFACE
    std::cout << "Creating parallel cut lines" << std::endl;
#endif

    for (std::vector<Point*>::iterator it = cut_points_for_lines.begin();
         it != (--cut_points_for_lines.end()); ++it)
    {
      std::vector<Point*>::iterator next = it;
      ++next;

      Point* p1 = *it;
      Point* p2 = *next;
      if (p1 == p2)
        dserror("Trying to create line between two points, which are the same! Point id is %d",
            p1->Id());
      mesh.NewLine(p1, p2, this, &other, element);
    }

    if (cut_points_for_lines[0] == cut_points_for_lines.back())
    {
      std::ofstream file("parallel_dump_first_equal_last.pos");
      // Dump cut points for lines and normal cut points
      for (std::vector<Point*>::iterator it = cut_points_for_lines.begin();
           it != cut_points_for_lines.end(); ++it)
      {
        std::stringstream pname;
        pname << "Point_for_line" << (*it)->Id();
        GEO::CUT::OUTPUT::GmshPointDump(file, *it, (*it)->Id(), pname.str(), false, NULL);
      }
      for (PointSet::const_iterator it = cut.begin(); it != cut.end(); ++it)
      {
        std::stringstream pname;
        pname << "CutPoint" << (*it)->Id();
        GEO::CUT::OUTPUT::GmshPointDump(file, *it, (*it)->Id(), pname.str(), false, NULL);
        (*it)->DumpConnectivityInfo();
      }
      GEO::CUT::OUTPUT::GmshSideDump(file, &other, std::string("OtherSide"));
      GEO::CUT::OUTPUT::GmshSideDump(file, this, std::string("ThisSide"));
      file.close();
      dserror("Trying to create line between two points, which are the same!");
    }
    mesh.NewLine(cut_points_for_lines[0], cut_points_for_lines.back(), this, &other, element);
    return true;
  }
  // else we still might consider different side
  else
  {
    return false;
  }
}

void GEO::CUT::Side::SimplifyMixedParallelCutSurface(Mesh& mesh, Element* element, Side& other,
    std::set<Point*>& new_surface, std::vector<Point*>& cut_points_for_lines)
{
  auto& existing_surfaces = HasMixedParallelCutSurface(new_surface) ? this->parallel_cut_surfaces_
                                                                    : other.parallel_cut_surfaces_;

  for (const std::set<Point*>& existing_surface : existing_surfaces)
  {
    std::vector<Point*> common_points;
    std::set_intersection(new_surface.begin(), new_surface.end(), existing_surface.begin(),
        existing_surface.end(), std::back_inserter(common_points));

    // one or zero common points should not cause problems
    if (common_points.size() > 1)
    {
      // might be that case with 3 common points will occur, however one would need to think
      // what exactly to merge then.
      if (common_points.size() == 2)
      {
        // we probably don't want to merge already merged-into or nodal points
        auto p_delete_it = std::find_if(common_points.begin(), common_points.end(),
            [](Point* p) { return (p->CutNode() == NULL && p->GetMergedPoints().size() == 0); });

        if (p_delete_it == common_points.end()) dserror("Cannot decide which point to merge");
        Point* p_delete = *p_delete_it;
        // pick up the next point in the list of common points
        Point* p_keep = common_points[(std::distance(common_points.begin(), p_delete_it) + 1) %
                                      common_points.size()];

        if (GEO::CUT::DistanceBetweenPoints(p_delete, p_keep) > 10.0 * MERGING_TOLERANCE)
          dserror("Trying to merge points, that are too far");

        std::cout << "WARNING: Perfoming simplification of the parallel surface geometry by "
                     "merging point "
                  << p_delete->Id() << " into " << p_keep->Id() << std::endl;

        p_delete->Replace(p_keep);

        cut_points_for_lines.erase(
            std::remove(cut_points_for_lines.begin(), cut_points_for_lines.end(), p_delete),
            cut_points_for_lines.end());
        new_surface.erase(p_delete);
        // NOTE: At this point we modify existing surfaces so we certainly don't want want to
        // continue, since original surfaces do not exist anymore
        break;
      }
      else
      {
        // the surface is already created, we are just here for the second call (maybe on another
        // side)
        if (std::set<Point*>(common_points.begin(), common_points.end()) != new_surface)
        {
          dserror(
              "There are more than 2 common points in the parallel cut surfaces. Dont know how to "
              "simplify!");
        }
      }
    }
  }
}

void GEO::CUT::Side::GetBoundaryCells(plain_boundarycell_set& bcells)
{
  for (std::vector<Facet*>::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->GetBoundaryCells(bcells);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeOwnedSideFacets(Mesh& mesh, Element* element, plain_facet_set& facets)
{
  // --- facet already created from another element ---------------------------
  // fill the output variable and return
  if (facets_.size() != 0)
  {
    std::copy(facets_.begin(), facets_.end(), std::inserter(facets, facets.begin()));
    return;
  }

  // create the pointgraph for the element_side
  Teuchos::RCP<IMPL::PointGraph> pg = Teuchos::rcp(IMPL::PointGraph::Create(
      mesh, element, this, IMPL::PointGraph::element_side, IMPL::PointGraph::all_lines));

  for (IMPL::PointGraph::facet_iterator i = pg->fbegin(); i != pg->fend(); ++i)
  {
    const Cycle& points = *i;

    Facet* f = mesh.NewFacet(points(), this, IsCutSide());
    if (f == NULL) run_time_error("failed to create owned facet");
    facets_.push_back(f);
  }

  for (IMPL::PointGraph::hole_iterator i = pg->hbegin(); i != pg->hend(); ++i)
  {
    const std::vector<Cycle>& hole = *i;

    // If we have a hole and multiple cuts we have to test which facet the
    // hole belongs to. Not supported now.
    unsigned int facetid = 0;
    if (facets_.size() != 1)
    {
      for (std::vector<Facet*>::iterator i = facets_.begin(); i != facets_.end(); ++i)
      {
        Facet* facet = *i;
        if (HoleOfFacet(*facet, hole))
        {
          break;
        }
        facetid++;
      }
      if (facetid == facets_.size())
      {
        run_time_error("failed to find the facet of the hole");
      }
    }

    for (std::vector<Cycle>::const_iterator i = hole.begin(); i != hole.end(); ++i)
    {
      const Cycle& points = *i;

      Facet* h = mesh.NewFacet(points(), this, false);
      facets_[facetid]->AddHole(h);
    }
  }

  std::copy(facets_.begin(), facets_.end(), std::inserter(facets, facets.begin()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeInternalFacets(Mesh& mesh, Element* element, plain_facet_set& facets)
{
  // create the pointgraph for the cut_side
  Teuchos::RCP<IMPL::PointGraph> pg = Teuchos::rcp(IMPL::PointGraph::Create(
      mesh, element, this, IMPL::PointGraph::cut_side, IMPL::PointGraph::all_lines));

  for (IMPL::PointGraph::facet_iterator i = pg->fbegin(); i != pg->fend(); ++i)
  {
    const Cycle& points = *i;
    MakeInternalFacets(mesh, element, points, facets);
  }
  for (IMPL::PointGraph::hole_iterator i = pg->hbegin(); i != pg->hend(); ++i)
  {
    const std::vector<Cycle>& hole = *i;
    for (std::vector<Cycle>::const_iterator i = hole.begin(); i != hole.end(); ++i)
    {
      const Cycle& points = *i;
      MakeInternalFacets(mesh, element, points, facets);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::MakeInternalFacets(
    Mesh& mesh, Element* element, const Cycle& points, plain_facet_set& facets)
{
  /* ignore cycles with all points on one and the same edge
   * ( ToDo seems meaningful for 3-D, check it for the other cases ) */
  // ignore cycles with points outside the current element
  if ((not points.IsValid() and element->Dim() == 3) or not points.IsCut(element))
  {
    // dserror("Shouldn't reach this case!");
    return;
  }

  Side* s = NULL;

  plain_side_set sides(element->Sides().begin(), element->Sides().end());
  points.Intersection(sides);


  if (sides.size() > 1)
  {
#if 0
    std::stringstream str;
    str << "can touch exactly one element side: "
      << points
      << "found sides:\n";
    std::copy( sides.begin(), sides.end(), std::ostream_iterator<Side*>( str, "\n" ) );


    int counter = 0;
    std::ofstream file("multiple_touched_sides.pos");
    for (plain_side_set::iterator it = sides.begin(); it != sides.end(); ++it, ++counter)
    {
      file << "View \"ElementSide" << counter << "\" {\n";
      GEO::CUT::OUTPUT::GmshSideDump(file, *it, NULL);
      file << "};\n";
    }
    file << "View \"ThisSide" << "\" {\n";
    GEO::CUT::OUTPUT::GmshSideDump(file, this, NULL);
    file << "};\n";
    file << "//Cycle dump\n";

    points.GmshDump( file );

    file.close();

    dserror( str.str() );
#endif
    {
      Facet* f = NULL;
      for (plain_side_set::iterator it = sides.begin(); it != sides.end(); ++it)
      {
        f = (*it)->FindFacet(points());
        if (f != NULL) break;
      }

      if (f != NULL)
        f->ExchangeSide(this, true);
      else
        f = mesh.NewFacet(points(), this, true);

      facets.insert(f);
      facets_.push_back(f);

#if EXTENDED_CUT_DEBUG_OUTPUT
      std::cout << "NOTICE:"
                << "cut_side facet touching multiple (" << sides.size() << ")"
                << "element sides" << std::endl;
#endif

      return;
    }

    // run_time_error( str.str() );
  }
  else if (sides.size() == 1)
  {
    s = *sides.begin();
  }
  // else s==NULL :-)

  if (s != NULL)
  {
    Facet* f = s->FindFacet(points());
    if (f != NULL)
    {
      f->ExchangeSide(this, true);
      facets.insert(f);
      facets_.push_back(f);
    }
    /* this case means that side pointers show, that the cut_facet lies on
     * the side of the element, but we cannot find this facet??? */
    else
    {
      // we throw an error at this point since this is not allowed -- hiermeier 11/16
      {
        int counter = 0;
        std::ofstream file("points_are_not_facet..pos");
        for (plain_side_set::iterator it = sides.begin(); it != sides.end(); ++it, ++counter)
        {
          std::stringstream side_name;
          side_name << "ElementSide" << counter;
          GEO::CUT::OUTPUT::GmshSideDump(file, *it, side_name.str());
        }
        GEO::CUT::OUTPUT::GmshSideDump(file, this, std::string("ThisSide"));
        file << "//Cycle dump\n";

        const std::vector<Facet*> side_facets = s->Facets();
        counter = 0;
        for (std::vector<Facet*>::const_iterator it = side_facets.begin(); it != side_facets.end();
             ++it, ++counter)
        {
          file << "View \"ThisCycle" << counter << "\" {\n";
          CUT::OUTPUT::GmshFacetDump(file, *it, "lines", true, false, NULL);
          file << "};\n";
        }

        points.GmshDump(file);

        file.close();


        std::cout << "\nOn cut side " << this;
        std::cout << "\nFound element side " << s;
        std::cout << "\nIs cut = " << (s->IsCut() ? "TRUE" : "FALSE");

        std::cout << "\n --- Mesh " << &mesh << " ---\n";
        for (std::map<int, Teuchos::RCP<Node>>::const_iterator cit = mesh.Nodes().begin();
             cit != mesh.Nodes().end(); ++cit)
          std::cout << "Point" << cit->second->point() << "\n";

        std::cout << "\n --- Element " << element << " ---\n";
        for (std::vector<Point*>::const_iterator cit = element->Points().begin();
             cit != element->Points().end(); ++cit)
          std::cout << "Point " << (*cit) << "\n";

        std::cout << "\n --- Side " << s << " ---\n";
        s->Print();

        std::cout << "\n\n --- PointCycle ---" << std::endl;
        for (std::vector<Point*>::const_iterator cit = points().begin(); cit != points().end();
             ++cit)
          std::cout << "Point " << (*cit) << "\n";


        std::cout << "\n\n -- Distance between points is --- " << std::endl;
        for (unsigned int i = 0; i < points().size(); ++i)
        {
          Point* first = points()[i];
          Point* next = points()[(i + 1) % points.size()];
          std::cout << "Between" << first->Id() << " and " << next->Id() << " is "
                    << GEO::CUT::DistanceBetweenPoints(first, next) << std::endl;
        }


        std::cout << "\n\n -- Distance between element points and cut_points --- " << std::endl;
        for (std::vector<Point*>::const_iterator eit = element->Points().begin();
             eit != element->Points().end(); ++eit)
        {
          Point* element_point = *eit;
          for (std::vector<Point*>::const_iterator cit = points().begin(); cit != points().end();
               ++cit)
          {
            Point* cycle_point = *cit;
            std::cout << "Between" << element_point->Id() << " and " << cycle_point->Id() << " is "
                      << GEO::CUT::DistanceBetweenPoints(element_point, cycle_point) << std::endl;
          }
        }

        std::cout << std::endl;

        run_time_error(
            "The point cycle was found on the side, but the facet was not found "
            "on the side. This shouldn't happen, since the facet has to be owned by the side "
            "in this case! Properly there is some hanging node in your cut mesh. Take a closer "
            "look.");
      }

      // multiple facets on one cut side within one element
      Facet* f = mesh.NewFacet(points(), this, true);
      facets.insert(f);
      facets_.push_back(f);
      //    dserror("must have matching facet on side");
    }
  }
  // if comon side is null, meaning cut_side does not lie on the element_side
  else
  {
    // insert new internal facet
    Facet* f = mesh.NewFacet(points(), this, true);
    facets.insert(f);
    facets_.push_back(f);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::OnSide(const PointSet& points)
{
  if (nodes_.size() == points.size())
  {
    for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
    {
      Node* n = *i;
      if (points.count(n->point()) == 0)
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::OnEdge(Point* point)
{
  for (std::vector<Edge*>::const_iterator i = edges_.begin(); i != edges_.end(); ++i)
  {
    Edge* e = *i;
    if (point->IsCut(e))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::OnEdge(Line* line)
{
  for (std::vector<Edge*>::const_iterator i = edges_.begin(); i != edges_.end(); ++i)
  {
    Edge* e = *i;
    if (line->OnEdge(e))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::AllPointsCommon(Side& side)
{
  const std::vector<Node*>& other_nodes = side.Nodes();
  for (std::vector<Node*>::const_iterator i = nodes_.begin(); i != nodes_.end(); ++i)
  {
    Point* p = (*i)->point();
    bool found = false;
    for (std::vector<Node*>::const_iterator j = other_nodes.begin(); j != other_nodes.end(); ++j)
    {
      Point* other_p = (*j)->point();
      if (other_p->Id() == p->Id())
      {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::HaveCommonNode(Side& side)
{
  const std::vector<Node*>& other_nodes = side.Nodes();
  for (std::vector<Node*>::const_iterator i = nodes_.begin(); i != nodes_.end(); ++i)
  {
    Node* e = *i;
    if (std::find(other_nodes.begin(), other_nodes.end(), e) != other_nodes.end())
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::HaveCommonEdge(Side& side)
{
  const std::vector<Edge*>& other_edges = side.Edges();
  for (std::vector<Edge*>::const_iterator i = edges_.begin(); i != edges_.end(); ++i)
  {
    Edge* e = *i;
    if (std::find(other_edges.begin(), other_edges.end(), e) != other_edges.end())
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Element* GEO::CUT::Side::CommonElement(Side* other)
{
  plain_element_set intersection;
  std::set_intersection(elements_.begin(), elements_.end(), other->elements_.begin(),
      other->elements_.end(), std::inserter(intersection, intersection.begin()));
  switch (intersection.size())
  {
    case 0:
      return NULL;
    case 1:
      return *intersection.begin();
    default:
      throw std::runtime_error("sides with more than one element in common");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::Print()
{
  std::cout << "[ ";
  for (std::vector<Edge*>::iterator i = edges_.begin(); i != edges_.end(); ++i)
  {
    Edge* e = *i;
    e->Print();
    std::cout << " ; ";
    if (i + 1 != edges_.end()) std::cout << "\n  ";
  }
  std::cout << " ]";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned GEO::CUT::Side::UncutFacetNumberPerSide() const
{
  if (elements_.size() < 1)
    dserror(
        "We need at least one registered element to be able to tell the "
        "default facet number of this side.");

  const unsigned ele_dim = elements_[0]->Dim();
  switch (ele_dim)
  {
    // 3-D elements: in the uncut case we have one facet per side
    case 3:
      return 1;
    /* 2-D elements: in the uncut case we have four facets per side ( lines
     *               surrounding the element / side ) */
    case 2:
      return 4;
    /* 1-D elements: in the uncut case we have two facets per side ( corner
     *               points of the element / side / edge ) */
    case 1:
      return 2;
    default:
      dserror("Unsupported parent element dimension! (ele->Dim() = %d )", ele_dim);
      exit(EXIT_FAILURE);
  }
  // can never be reached
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::IsCut()
{
  if (facets_.size() > UncutFacetNumberPerSide())
  {
    return true;
  }
  if (facets_[0]->OnCutSide())
  {
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, GEO::CUT::Side& s)
{
  stream << "side: {";
  const std::vector<GEO::CUT::Node*>& nodes = s.Nodes();
  for (std::vector<GEO::CUT::Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    GEO::CUT::Node* n = *i;
    n->point()->Print(stream);
    stream << ",";
  }
  stream << "}";
  return stream;
}

/*----------------------------------------------------------------------------*
 *  Gets the selfcutposition and spreads the positional information wirtz 05/13
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::GetSelfCutPosition(Point::PointPosition position)
{
  if (selfcutposition_ != position)
  {
    if (position == Point::oncutsurface) dserror("if (position == Point::oncutsurface)");

    if (Id() < 0)  // Propagate Selfcutposition only on Cutsides :)
      return;

    selfcutposition_ = position;

    for (std::vector<Edge*>::iterator i = edges_.begin(); i != edges_.end(); ++i)
    {
      Edge* e = *i;
      Point::PointPosition ep = e->SelfCutPosition();
      if (ep == Point::undecided)
      {
        e->SelfCutPosition(position);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *  Changes the selfcutposition of this cutside and spreads the positional
 *  information                                                    wirtz 07/16
 *----------------------------------------------------------------------------*/
void GEO::CUT::Side::ChangeSelfCutPosition(Point::PointPosition position)
{
  if (selfcutposition_ != position)
  {
    if (position == Point::oncutsurface) dserror("if (position == Point::oncutsurface)");
    selfcutposition_ = position;

    for (std::vector<Edge*>::iterator i = edges_.begin(); i != edges_.end(); ++i)
    {
      Edge* e = *i;
      e->ChangeSelfCutPosition(position);
    }
  }
}

/*----------------------------------------------------------------------------*
 * returns true if the hole is inside the facet                     wirtz 05/13
 *----------------------------------------------------------------------------*/
bool GEO::CUT::Side::HoleOfFacet(Facet& facet, const std::vector<Cycle>& hole)
{
  int intersectioncount = 0;
  bool intersectioninpoint = true;
  std::vector<Point*> facetpoints = facet.Points();
  int facetsize = facetpoints.size();
  std::vector<LINALG::Matrix<3, 1>> facetpointslocalcoord;
  facetpointslocalcoord.reserve(facetsize);
  for (std::vector<Point*>::iterator i = facetpoints.begin(); i != facetpoints.end(); ++i)
  {
    Point* facetpoint = *i;
    LINALG::Matrix<3, 1> pointcoord;
    facetpoint->Coordinates(pointcoord.A());
    LINALG::Matrix<3, 1> pointlocalcoord;
    LocalCoordinates(pointcoord, pointlocalcoord, false);
    facetpointslocalcoord.push_back(pointlocalcoord);
  }
  LINALG::Matrix<3, 1> holepointcoord;
  LINALG::Matrix<3, 1> holepointlocalcoord;
  hole[0]()[0]->Coordinates(holepointcoord.A());
  LocalCoordinates(holepointcoord, holepointlocalcoord, false);
  double epsilon = 0;
  while (intersectioninpoint)
  {
    intersectioninpoint = false;
    for (std::vector<LINALG::Matrix<3, 1>>::iterator i = facetpointslocalcoord.begin();
         i != facetpointslocalcoord.end(); ++i)
    {
      LINALG::Matrix<3, 1> facetpoint1 = *i;
      LINALG::Matrix<3, 1> facetpoint2;
      if (i + 1 != facetpointslocalcoord.end())
      {
        facetpoint2 = *(i + 1);
      }
      else
      {
        facetpoint2 = *(i + 1 - facetsize);
      }
      double A = facetpoint1(0) - facetpoint2(0);
      double B = facetpoint1(1) - facetpoint2(1);
      double C = facetpoint1(0) + facetpoint2(0) - 2 * holepointlocalcoord(0) - 2;
      double D = facetpoint1(1) + facetpoint2(1) - 2 * holepointlocalcoord(1) - epsilon;
      double N = 2 * B - epsilon * A;
      if (abs(N) > REFERENCETOL)
      {
        double eta = (B * C - A * D) / N;
        double xsi = (2 * D - epsilon * C) / N;
        if (eta < 1 and eta > -1 and xsi < 1 and xsi > -1)
        {
          intersectioncount++;
          double xlocalcoord = holepointlocalcoord(0) + 1 + eta;
          double ylocalcoord = (2 * holepointlocalcoord(1) + epsilon + epsilon * eta) / 2;
          if ((abs(xlocalcoord - facetpoint1(0)) < REFERENCETOL and
                  abs(ylocalcoord - facetpoint1(1)) < REFERENCETOL) or
              (abs(xlocalcoord - facetpoint2(0)) < REFERENCETOL and
                  abs(ylocalcoord - facetpoint2(1)) < REFERENCETOL))
          {
            intersectioninpoint = true;
            intersectioncount = 0;
            epsilon += REFERENCETOL;
            break;
          }
        }
      }
    }
  }
  if (intersectioncount % 2 == 0)
  {
    return false;
  }
  else
  {
    return true;
  }
}

/*--------------------------------------------------------------------------*
 * replace the Node "nod" with the new node "replwith"         sudhakar 10/13
 * Modify also the edge informations correspondingly
 *--------------------------------------------------------------------------*/
void GEO::CUT::Side::replaceNodes(Node* nod, Node* replwith)
{
  bool replaced = false;

  for (unsigned i = 0; i < nodes_.size(); i++)
  {
    Node* orig = nodes_[i];

    if (orig->Id() == nod->Id())
    {
      nodes_[i] = replwith;
      replaced = true;
    }
  }

  if (not replaced) return;

  // also modify the corresponding edge information
  for (std::vector<Edge*>::iterator it = edges_.begin(); it != edges_.end(); it++)
  {
    Edge* ed = *it;
    ed->replaceNode(nod, replwith);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType sidetype, unsigned numNodesSide,
    unsigned dim>
bool GEO::CUT::ConcreteSide<probdim, sidetype, numNodesSide, dim>::IsCloserSide(
    const LINALG::Matrix<probdim, 1>& startpoint_xyz, GEO::CUT::Side* other, bool& is_closer)
{
  /* shoot a ray starting from the startpoint through the midpoint of this side
   * and find an intersection point with the other side */
  LINALG::Matrix<probdim, 1> ray_point_xyz(true);
  // as second point on the ray we define the midpoint of this side


  /* choose a point inside the side such that the angle between the normal of
   * the other side and the vector between start-point and ray-point is next to
   * 0 or 180 degree to avoid orthogonal ray and normal vector therefore, for
   * almost parallel sides and a start-point next to the common point of the sides
   * the side-center might lead to a ray which is almost parallel to the other side */

  LINALG::Matrix<dim, numNodesSide> corner_coords_rst(true);
  this->LocalCornerCoordinates(corner_coords_rst.A());

  /* shrink/perturb the local coordinates around the center point with a given
   * tolerance to obtain points which are next to the corner points however slightly
   * inside */
  LINALG::Matrix<dim, 1> rst_center = DRT::UTILS::getLocalCenterPosition<dim>(this->Shape());

  //-----------------------------
  // get perturbed coordinates
  //-----------------------------
  LINALG::Matrix<dim, numNodesSide> inner_corner_coords_rst(true);

  // 1. transform such that coordinates center is located in the element center
  for (unsigned i = 0; i < numNodesSide; ++i)
  {
    for (unsigned j = 0; j < dim; ++j)
      inner_corner_coords_rst(j, i) = corner_coords_rst(j, i) - rst_center(j);
  }

  /* 2. shrink the side coordinates (hard-coded tolerance w.r.t. local coordinates of
   *    element, this should be fine) */
  const double TOL = 1e-003;
  const double scalefac = 1.0 - TOL;
  inner_corner_coords_rst.Scale(scalefac);

  // 3. transform the element back
  for (unsigned i = 0; i < numNodesSide; ++i)
  {
    for (unsigned j = 0; j < dim; ++j) inner_corner_coords_rst(j, i) += rst_center(j);
  }

  //-----------------------------
  /* choose the center point point or a perturbed inner corner point
   * such that the ray between startpoint and this point is as perpendicular
   * as possible to guarantee well-conditioned systems for finding ray-cut points */
  //-----------------------------

  LINALG::Matrix<probdim, 1> xyz(true);
  LINALG::Matrix<probdim, 1> ray_dir(true);

  // get normal of other side at its center
  LINALG::Matrix<probdim, 1> t1, t2, n(true);
  other->BasisAtCenter(t1, t2, n);


  // start with center point
  this->SideCenter(xyz);  // as second point on the ray we define the midpoint of this side
  ray_dir.Update(1.0, xyz, -1.0, startpoint_xyz, 0.0);
  ray_point_xyz.Update(1.0, xyz, 0.0);

  double cosine = ray_dir.Dot(n) / ray_dir.Norm2();  // n is normalized


  // loop corner nodes
  for (unsigned i = 0; i < numNodesSide; i++)
  {
    // get ray vector (endpoint-startpoint)
    LINALG::Matrix<dim, 1> rst_inner_corner(&inner_corner_coords_rst(0, i), true);
    PointAt(rst_inner_corner, xyz);
    ray_dir.Update(1.0, xyz, -1.0, startpoint_xyz, 0.0);

    double cosine_tmp = ray_dir.Dot(n) / ray_dir.Norm2();  // n is normalized

    /* maximize the absolute value of the cosine to choose the ray which
     * is as perpendicular to the side as possible */
    if (fabs(cosine_tmp) > fabs(cosine))
    {
      ray_point_xyz.Update(1.0, xyz, 0.0);
      cosine = cosine_tmp;
    }
  }
  //-----------------------------
  /* shoot the ray and find a cutpoint with the other side's plane or curved
   * surface space */
  //-----------------------------
  LINALG::Matrix<dim, 1> rs(true);
  double line_xi = 0.0;

  bool cut_found = other->RayCut(startpoint_xyz, ray_point_xyz, rs, line_xi);

  if (!cut_found)
    return false;
  else
  {
    // The main decision if the side lies closer to the start-point than the other side

    //    std::cout << "line_xi " << line_xi << std::endl;

    if (line_xi > 1.0 + REFERENCETOL)
    {
      // the first side is closer to the start point than the second side
      is_closer = true;
      return true;
    }
    else if (fabs(line_xi - 1.0) <= REFERENCETOL)
    {
      // the found intersection point on the other is the midpoint of the fist side
      // in that case both sides lie within one plane
      // this case is catched in SameNormal afterwards

      // std::cout << "line_xi " << line_xi << std::endl;
      // std::cout << "check if both sides lie in one plane " << std::endl;

      return false;
    }
    else if (line_xi < 1.0 - REFERENCETOL and line_xi > -1.0 + REFERENCETOL)
    {
      // the most safe check (accept the other side as the nearest side)
      is_closer = false;
      return true;
    }
    else if (fabs(line_xi + 1.0) <= REFERENCETOL)
    {
      /* the intersection point is the same as the start-point of the ray
       * the other side contains the start-point and the cut-point shared
       * with the original side in that case the line between the start-point
       * and the cut-point lies in the second side then the side is orthogonal
       * to the line this side should be removed in */
      std::cout << "line_xi " << line_xi << std::endl;
      std::cout << "start-point: " << startpoint_xyz << std::endl;
      std::cout << "side orthogonal ? " << std::endl;
      other->Print();

      throw std::runtime_error("IsCloserSide along the ray-tracing line failed! ");

      return false;
    }
    else if (line_xi < -1.0 - REFERENCETOL)
    {
      // both sides lead to the same result
      is_closer = true;  // false would be also okay
      return true;
    }
    else
    {
      // undermined range of local coordinates!

      std::cout << "line_xi " << line_xi << std::endl;
      std::cout << "cut point found, but the local line coordinates along the "
                   "ray-tracing line lies in undefined region"
                << std::endl;

      throw std::runtime_error("IsCloserSide along the ray-tracing line failed! ");
    }

    return false;
  }

  return false;  // return not successful
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType sidetype, unsigned numNodesSide,
    unsigned dim>
///  lies point with given coordinates within this side?
bool GEO::CUT::ConcreteSide<probdim, sidetype, numNodesSide, dim>::WithinSide(
    const LINALG::Matrix<probdim, 1>& xyz, LINALG::Matrix<dim, 1>& rs, double& dist)
{
  dserror("Do we use this function?");
  Teuchos::RCP<Position> pos = PositionFactory::BuildPosition<probdim, sidetype>(*this, xyz);
  bool success = pos->IsGivenPointWithinElement();
  if (not success)
  {
    throw std::runtime_error("ComputeDistance w.r.t side not successful");
    rs = 0;
    dist = 9999;  // set large value
    return false;
  }
  pos->LocalCoordinates(rs);
  dist = pos->Distance();


  return pos->WithinLimits(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType sidetype, unsigned numNodesSide,
    unsigned dim>
///  lies point with given coordinates within this side?
bool GEO::CUT::ConcreteSide<probdim, sidetype, numNodesSide, dim>::LocalCoordinates(
    const LINALG::Matrix<probdim, 1>& xyz, LINALG::Matrix<probdim, 1>& rsd, bool allow_dist,
    double tol)
{
  Teuchos::RCP<Position> pos = PositionFactory::BuildPosition<probdim, sidetype>(*this, xyz);
  bool success = pos->Compute(tol, allow_dist);
  LINALG::Matrix<dim, 1> rs(true);
  if (pos->Status() == Position::position_valid) pos->LocalCoordinates(rs);
  // copy the position
  std::copy(rs.A(), rs.A() + dim, &rsd(0));
  // copy the distance
  switch (dim)
  {
    /* Side dimension is by a factor of one smaller than the problem dimension,
     * i.e. surface element in 3-D or line element in 2-D. */
    case probdim - 1:
    {
      if (pos->Status() >= Position::position_distance_valid) rsd(dim) = pos->Distance();
      break;
    }
    /* Side dimension is by a factor of two smaller than the problem dimension,
     * i.e. line element in 3-D. */
    case probdim - 2:
    {
      dserror(
          "Unsupported dim / probdim combination. I think this can't happen "
          "for side elements. If I'm wrong, just ask me. -- hiermeier");
      exit(EXIT_FAILURE);
    }
    /* For a 1-D and a 2-D problem, the side dimension is equal to the problem
     * dimension */
    default:
    {
      /* do nothing, there is no distance in the 2-D non-embedded standard case */
      break;
    }
  }

  return success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType sidetype, unsigned numNodesSide,
    unsigned dim>
///  lies point with given coordinates within this side?
bool GEO::CUT::ConcreteSide<probdim, sidetype, numNodesSide, dim>::RayCut(
    const LINALG::Matrix<probdim, 1>& p1_xyz, const LINALG::Matrix<probdim, 1>& p2_xyz,
    LINALG::Matrix<dim, 1>& rs, double& line_xi)
{
  LINALG::Matrix<probdim, numNodesSide> xyze_surface(true);
  this->Coordinates(xyze_surface);

  LINALG::Matrix<probdim, 2> xyze_line(true);
  for (unsigned i = 0; i < probdim; ++i)
  {
    xyze_line(i, 0) = p1_xyz(i);
    xyze_line(i, 1) = p2_xyz(i);
  }
  /* The dim+1 entry corresponds to the parameter space coordinate of the
   * 1-D line element. */
  LINALG::Matrix<dim + 1, 1> xsi(true);

  // do not check for within-limits during the Newton-scheme, since the cut-point is
  // allowed to be not within the side and line
  bool checklimits = false;
  // do not use cln here
  GEO::CUT::KERNEL::ComputeIntersection<probdim, DRT::Element::line2, sidetype, false> ci(
      xsi, checklimits);
  // GEO::CUT::KERNEL::DebugComputeIntersection<probdim,DRT::Element::line2, sidetype,false>
  //      ci( xsi, checklimits );

  // successful line-side intersection
  if (ci(xyze_surface, xyze_line))
  {
    std::copy(xsi.A(), xsi.A() + dim, &rs(0));
    line_xi = xsi(dim);

    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::SideFactory::CreateSide(DRT::Element::DiscretizationType sidetype,
    int sid, const std::vector<Node*>& nodes, const std::vector<Edge*>& edges) const
{
  Side* s = NULL;
  const int probdim = DRT::Problem::Instance()->NDim();
  switch (sidetype)
  {
    case DRT::Element::line2:
    {
      s = CreateConcreteSide<DRT::Element::line2>(sid, nodes, edges, probdim);
      break;
    }
    case DRT::Element::tri3:
      s = CreateConcreteSide<DRT::Element::tri3>(sid, nodes, edges, probdim);
      break;
    case DRT::Element::quad4:
      s = CreateConcreteSide<DRT::Element::quad4>(sid, nodes, edges, probdim);
      break;
    default:
    {
      dserror(
          "Unsupported side type! ( %d | %s )", sidetype, DRT::DistypeToString(sidetype).c_str());
      break;
    }
  }
  return s;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Side::CreateLevelSetSide(const int& sid)
{
  Side* lvs_side_ptr = NULL;

  const int probdim = DRT::Problem::Instance()->NDim();
  switch (probdim)
  {
    case 2:
      lvs_side_ptr = new LevelSetSide<2>(sid);
      break;
    case 3:
      lvs_side_ptr = new LevelSetSide<3>(sid);
      break;
    default:
      dserror("Unsupported problem dimension! (probdim=%d)", probdim);
      break;
  }
  return lvs_side_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Side::Create(const DRT::Element::DiscretizationType& sidetype,
    const int& sid, const std::vector<Node*>& nodes, const std::vector<Edge*>& edges)
{
  SideFactory factory;
  return factory.CreateSide(sidetype, sid, nodes, edges);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Side* GEO::CUT::Side::Create(const unsigned& shardskey, const int& sid,
    const std::vector<Node*>& nodes, const std::vector<Edge*>& edges)
{
  return Create(DRT::ShardsKeyToDisType(shardskey), sid, nodes, edges);
}


template class GEO::CUT::ConcreteSide<2, DRT::Element::line2>;
template class GEO::CUT::ConcreteSide<3, DRT::Element::line2>;
template class GEO::CUT::ConcreteSide<3, DRT::Element::quad4>;
template class GEO::CUT::ConcreteSide<3, DRT::Element::tri3>;
