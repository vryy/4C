/*!-----------------------------------------------------------------------------------------------*
 \file cut_edge.cpp

 \brief class representing a geometrical edge

\level 2

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/

#include "cut_position.H"
#include "cut_intersection.H"
#include "cut_facet.H"
#include "cut_point_impl.H"
#include "cut_levelsetside.H"
#include "cut_options.H"

#include "../drt_lib/drt_globalproblem.H"

#include <string>
#include <stack>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Edge> GEO::CUT::Edge::Create(
    DRT::Element::DiscretizationType edgetype, const std::vector<Node*>& nodes)
{
  EdgeFactory factory;
  return factory.CreateEdge(edgetype, nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Edge> GEO::CUT::Edge::Create(
    unsigned shardskey, const std::vector<Node*>& nodes)
{
  return Edge::Create(DRT::ShardsKeyToDisType(shardskey), nodes);
}

bool GEO::CUT::Edge::FindCutPointsLevelSet(Mesh& mesh, Element* element, Side& side, Side& other)
{
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->IsCut(&other))  // cut points may already obtained when some other side was considered
    {
      p->AddElement(element);
    }
  }

  // test for the cut of edge and side first time in this mesh

  PointSet cut_points;
  other.Cut(mesh, *this, cut_points);
  for (PointSet::iterator i = cut_points.begin(); i != cut_points.end(); ++i)
  {
    Point* p = *i;

    p->AddEdge(this);
    p->AddElement(element);

    // These adds are implicitly done, but for documentation do them all explicitly.

    p->AddSide(&side);
    p->AddSide(&other);
    AddPoint(p);
  }

  return cut_points.size() > 0;
}

bool GEO::CUT::Edge::FindCutPoints(Mesh& mesh, Element* element, Side& side, Side& other)
{
  // dispatch function call between side = LevelSetSide and  normal side
  return other.FindCutPointsDispatch(mesh, element, side, *this);
}


/*-----------------------------------------------------------------------------*
 *  Find points at which this edge which is in "side" cuts the "other"
 *-----------------------------------------------------------------------------*/
bool GEO::CUT::Edge::FindCutPointsMeshCut(
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
    if (p->IsCut(&other) and p->IsCut(&other, this) and p->IsCut(this))
    {
      point_stack.push_back(p);
    }
  }
  // everything is handled

  // this should not be possible
  if (point_stack.size() > 2)
  {
    GEO::CUT::OUTPUT::DebugDump_MoreThanTwoIntersectionPoints(this, &other, point_stack);
    dserror("Line x Side has more than 2 intersection points.Namely %u", point_stack.size());
  }

  // if 2 points found - we are done
  else if (point_stack.size() == 2)
  {
    point_stack[0]->AddElement(element);
    point_stack[1]->AddElement(element);
    cut = true;
  }
  // might be possible to have another intersection point (if edge is parallel to the size)
  else if (point_stack.size() == 1)
  {
    Point* p = point_stack[0];
    if ((BeginNode()->point() != p) and (EndNode()->point() != p))
    {
      // try to compute assuming parallel cut ( no real intersection )
      parallel_cut = other.JustParallelCut(mesh, *this, cut_points);
      // if edge is not parallel, nothing to do here, one point already intersected and no more can
      // be found for normal cases, otherwise we need to add some extra connection
      if (not parallel_cut)
      {
        cut = true;
        p->AddElement(element);
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
      int id_skip = (BeginNode()->point() == p) ? 0 : 1;
      parallel_cut = other.JustParallelCut(mesh, *this, cut_points, id_skip);
      if (not parallel_cut) cut = true;
    }
  }
  else
  {
    // not yet cut at all, so do full cut
    other.Cut(mesh, *this, cut_points);
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

    p->AddEdge(this);
    p->AddElement(element);

    // These adds are implicitly done, but for documentation do them all explicitly.

    p->AddSide(&side);
    p->AddSide(&other);

    p->AddPair(&other, this);

    AddPoint(p);
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
void GEO::CUT::Edge::GetCutPoints(Element* element, Side& side, Side& other, PointSet& cuts)
{
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->IsCut(&other) and p->IsCut(element) and p->IsCut(&other, this))
    {
      cuts.insert(p);
    }
  }
}

/*----------------------------------------------------------------------------*
 * Cut points falling on this edge that are common to the given edge are
 * extracted
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::GetCutPoints(Edge* other, PointSet& cuts)
{
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->IsCut(other))
    {
      cuts.insert(p);
    }
  }
}

/// Find common points (excluding cut_points points) between two edges
void GEO::CUT::Edge::CommonNodalPoints(Edge* other, std::vector<Point*>& common)
{
  const std::vector<Node*>& other_nodes = other->Nodes();
  const std::vector<Node*>& my_nodes = Nodes();
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
void GEO::CUT::Edge::AddPoint(Point* cut_point)
{
  // make sure the position of the point on this edge is known
  cut_point->t(this);

#if 0
  std::vector<Point*>::iterator j = std::lower_bound( cut_points_.begin(), cut_points_.end(), cut_point, PointPositionLess( this ) );
  if ( j==cut_points_.end() or *j!=cut_point )
  {
    cut_points_.push_back( cut_point );
    std::sort( cut_points_.begin(), cut_points_.end(), PointPositionLess( this ) );
  }
#else
  cut_points_.insert(cut_point);
#endif

#if 0
#ifdef DEBUGCUTLIBRARY
  PointSet cp;
  std::copy( cut_points_.begin(), cut_points_.end(), std::inserter( cp, cp.begin() ) );
  if ( cut_points_.size() != cp.size() )
    throw std::runtime_error( "broken cut points" );
#endif
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPoint(Node* edge_start, Node* edge_end, std::vector<Point*>& edge_points)
{
  Point* bp = BeginNode()->point();
  Point* ep = EndNode()->point();
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
    throw std::runtime_error("not a valid edge");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPoints(Side* side, PointSet& cut_points)
{
  IMPL::SideCutFilter filter(side);
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    plain_line_set cut_lines;
    p->CutLines(filter, cut_lines);
    if (cut_lines.size() > 0)
    {
      cut_points.insert(p);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPointsBetween(Point* begin, Point* end, std::vector<Point*>& line)
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
      throw std::runtime_error("different points at the same place");
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPointsIncluding(Point* begin, Point* end, std::vector<Point*>& line)
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
      throw std::runtime_error("different points at the same place");
    }
  }

  // to detect if cut point local coordinates is outside of edges endpoints
  if (*bi_or != begin) throw std::runtime_error("beginning point not on edge");
  if (*et_or != end) throw std::runtime_error("end point not on edge");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::CutPointsInside(Element* element, std::vector<Point*>& line)
{
  Point* first = NULL;
  Point* last = NULL;
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->IsCut(element))
    {
      if (first == NULL)
      {
        first = last = p;
      }
      else
      {
        last = p;
      }
    }
  }
  if (first != NULL and first != last)
  {
    CutPointsIncluding(first, last, line);
    if (line.size() > 2)
      for (std::vector<Point*>::iterator it = line.begin(); it != line.end(); /* */)
      {
        if (not(*it)->IsCut(element))
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
bool GEO::CUT::Edge::IsCut(Side* side)
{
  // cutpoints contains end-points and internal cut-points
  for (PointPositionSet::iterator i = cut_points_.begin(); i != cut_points_.end(); ++i)
  {
    Point* p = *i;
    if (p->IsCut(side))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::Edge::NodeInElement(Element* element, Point* other)
{
  Point* p = BeginNode()->point();
  if (p != other and p->IsCut(element))
  {
    return p;
  }
  p = EndNode()->point();
  if (p != other and p->IsCut(element))
  {
    return p;
  }
  return NULL;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Edge::RectifyCutNumerics()
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
      const plain_side_set& cutsides = p->CutSides();
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
              p->AddSide(s);
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
void GEO::CUT::Edge::SelfCutPosition(Point::PointPosition pos)
{
#ifdef DEBUGCUTLIBRARY
  if ((selfcutposition_ == Point::inside and pos == Point::outside) or
      (selfcutposition_ == Point::outside and pos == Point::inside))
  {
    dserror(
        "Are you sure that you want to change the edge-position from inside to outside or vice "
        "versa?");
  }
#endif

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
          if (n->SelfCutPosition() == Point::undecided)
          {
            n->SelfCutPosition(pos);
          }
        }
        for (plain_side_set::iterator i = sides_.begin(); i != sides_.end(); ++i)
        {
          Side* s = *i;
          s->GetSelfCutPosition(pos);
        }
      }
    }
  }
}

/*------------------------------------------------------------------------*
 *  Changes the selfcutposition of this edge and spreads the positional
 *  information                                                 wirtz 07/16
 *------------------------------------------------------------------------*/
void GEO::CUT::Edge::ChangeSelfCutPosition(Point::PointPosition pos)
{
  if (selfcutposition_ != pos)
  {
    selfcutposition_ = pos;
    for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
    {
      Node* n = *i;
      n->ChangeSelfCutPosition(pos);
    }
    for (plain_side_set::iterator i = sides_.begin(); i != sides_.end(); ++i)
    {
      Side* s = *i;
      s->ChangeSelfCutPosition(pos);
    }
  }
}

/*------------------------------------------------------------------------*
 *  Replaces the node "nod" of the edge with given node "replwith"
 *                                                              sudhakar 09/13
 *------------------------------------------------------------------------*/
void GEO::CUT::Edge::replaceNode(Node* nod, Node* replwith)
{
  for (unsigned i = 0; i < nodes_.size(); i++)
  {
    Node* orig = nodes_[i];

    if (orig->Id() == nod->Id())
    {
      nodes_[i] = replwith;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probDim, DRT::Element::DiscretizationType edgeType, unsigned dimEdge,
    unsigned numNodesEdge>
bool GEO::CUT::ConcreteEdge<probDim, edgeType, dimEdge, numNodesEdge>::Cut(
    Mesh& mesh, Side& side, PointSet& cuts)
{
  Teuchos::RCP<GEO::CUT::IntersectionBase> inter_ptr = IntersectionPtr(side.Shape());

  inter_ptr->Init(&mesh, this, &side, false, false, true);
  return inter_ptr->Intersect(cuts);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probDim, DRT::Element::DiscretizationType edgeType, unsigned dimEdge,
    unsigned numNodesEdge>
bool GEO::CUT::ConcreteEdge<probDim, edgeType, dimEdge, numNodesEdge>::JustParallelCut(
    Mesh& mesh, Side& side, PointSet& cuts, int skip_id)
{
  Teuchos::RCP<GEO::CUT::IntersectionBase> inter_ptr = IntersectionPtr(side.Shape());

  inter_ptr->Init(&mesh, this, &side, false, false, false);
  return (inter_ptr->HandleParallelIntersection(cuts, skip_id) > 0);
}

template <unsigned probDim, DRT::Element::DiscretizationType edgeType, unsigned dimEdge,
    unsigned numNodesEdge>
bool GEO::CUT::ConcreteEdge<probDim, edgeType, dimEdge, numNodesEdge>::HandleParallelCut(
    Edge* other, Side* side, PointSet* cut_points, INPAR::CUT::CUT_Floattype floattype)
{
  PointSet parallel_cuts;
  if (numNodesEdge != 2 or other->NumNodes() != 2)
    dserror("This method works only for 2 line2 intersections!");

  std::vector<Node*> touch_nodes;
  GetTouchingPoints(other->Nodes(), touch_nodes, floattype);

  unsigned int nodes_e1 = touch_nodes.size();
  other->GetTouchingPoints(Nodes(), touch_nodes, floattype);

  // nodes that belong to e1
  for (unsigned int i = 0; i < nodes_e1; ++i)
  {
    Point* cut_point = Point::InsertCut(this, side, touch_nodes[i]);
    std::stringstream msg;
#if CUT_CREATION_INFO
    msg << "//Added because " << cut_point->Id() << " lies on the edge" << this->Id() << std::endl;
    msg << "//Added during parallel edges checking" << std::endl;
#endif
    cut_point->AddEdgeIntersection(this, other, side, this, msg.str());
    if (cut_points) cut_points->insert(cut_point);

      // if a nodal point inserted on the other edge, we also want to add all sides of that node
#if CUT_CREATION_INFO
    msg << "//Added because of intersection of the edges of the sides of this nodal point with the "
           "edge "
        << this->Id() << std::endl;
#endif
    const plain_edge_set& edges = touch_nodes[i]->Edges();
    for (plain_edge_set::const_iterator it = edges.begin(); it != edges.end(); ++it)
    {
      Edge* e = *it;
      cut_point->AddEdgeIntersection(this, e, side, this, msg.str());
    }
  }

  // nodes that belong to e2
  for (unsigned int i = nodes_e1; i < touch_nodes.size(); ++i)
  {
    Point* cut_point = Point::InsertCut(other, side, touch_nodes[i]);
    std::stringstream msg;
#if CUT_CREATION_INFO
    msg << "//Added because " << cut_point->Id() << " lies on the edge" << other->Id() << std::endl;
    msg << "//Added during parallel edges checking" << std::endl;
#endif
    cut_point->AddEdgeIntersection(this, other, side, this, msg.str());
    if (cut_points) cut_points->insert(cut_point);
      // if a nodal point inserted on the other edge, we also want to add all sides of that node
#if CUT_CREATION_INFO
    msg << "//Added because of intersection of the edges of the sides of this nodal point with the "
           "edge "
        << other->Id() << std::endl;
#endif
    const plain_edge_set& edges = touch_nodes[i]->Edges();
    for (plain_edge_set::const_iterator it = edges.begin(); it != edges.end(); ++it)
    {
      Edge* e = *it;
      cut_point->AddEdgeIntersection(this, e, side, this, msg.str());
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
      dserror("Something went wrong!");
  }
  return false;
}

template <unsigned probDim, DRT::Element::DiscretizationType edgeType, unsigned dimEdge,
    unsigned numNodesEdge>
void GEO::CUT::ConcreteEdge<probDim, edgeType, dimEdge, numNodesEdge>::GetTouchingPoints(
    const std::vector<Node*>& nodes, std::vector<Node*>& touch_nodes,
    INPAR::CUT::CUT_Floattype floattype)
{
  bool signeddistance = false;
  double distance = 0;
  LINALG::Matrix<probDim, 1> xsi;
  LINALG::Matrix<probDim, numNodesEdge> xyze_edge;
  Coordinates(xyze_edge.A());

  const std::vector<Node*> edge_nodes = Nodes();

  for (std::vector<Node*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
  {
    xsi = 0.0;
    distance = 0;

    Node* node = *it;
    Point* p = node->point();
    LINALG::Matrix<probDim, 1> p_coord;
    p->Coordinates(p_coord.A());

    bool conv = false;
    KERNEL::PointOnSurfaceLoc loc;
    bool SurfaceWithinLimits0 = false;
    switch (floattype)
    {
      case INPAR::CUT::floattype_double:
      {
        KERNEL::ComputeDistance<probDim, edgeType, false> cd(xsi);
        conv = cd(xyze_edge, p_coord, distance, signeddistance);
        loc = cd.GetSideLocation();
        if (conv && loc.WithinSide() && loc.OnSide())
          SurfaceWithinLimits0 = cd.SurfaceWithinLimits(0.0);
        break;
      }
      case INPAR::CUT::floattype_cln:
      {
        KERNEL::ComputeDistance<probDim, edgeType, true> cd(xsi);
        conv = cd(xyze_edge, p_coord, distance, signeddistance);
        loc = cd.GetSideLocation();
        if (conv && loc.WithinSide() && loc.OnSide())
          SurfaceWithinLimits0 = cd.SurfaceWithinLimits(0.0);
        break;
      }
      default:
      {
        throw std::runtime_error("Unexpected floattype for KERNEL::ComputeDistance!");
      }
    }

    if (!conv)
      throw std::runtime_error("Newton did not converge for parallelism detection of lines!");

    if (SurfaceWithinLimits0) touch_nodes.push_back(node);

    // Extra check of the connection to the nodal points of the edge
    for (std::vector<Node*>::const_iterator jt = edge_nodes.begin(); jt != edge_nodes.end(); ++jt)
    {
      Node* edge_node = *jt;
      Point* edge_point = edge_node->point();
      if (GEO::CUT::DistanceBetweenPoints(p, edge_point) < SIDE_DETECTION_TOLERANCE)
      {
        p->DumpConnectivityInfo();
        edge_point->DumpConnectivityInfo();
        std::stringstream err_msg;
        err_msg << "Distance between points is " << std::setprecision(15)
                << GEO::CUT::DistanceBetweenPoints(p, edge_point)
                << " This two points should have been merged!";
        dserror(err_msg.str());
      }
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probDim, DRT::Element::DiscretizationType edgeType, unsigned dimEdge,
    unsigned numNodesEdge>
bool GEO::CUT::ConcreteEdge<probDim, edgeType, dimEdge, numNodesEdge>::ComputeCut(
    Mesh* mesh, Edge* other, Side* side, PointSet* cut_points, double& tolerance)
{
#if CUT_CREATION_INFO
  static unsigned long int edge_inter_id = 0;
  edge_inter_id++;
#endif

  // nodal edge coordinates
  LINALG::Matrix<probDim, numNodesEdge> xyze_this;
  Coordinates(xyze_this.A());

  // nodal pseudo side coordinates
  /* Note: We insert the this edge as side. In this way we can use the
   * standard intersection object. */
  Epetra_SerialDenseMatrix xyze_other;
  xyze_other.Shape(other->ProbDim(), other->NumNodes());
  other->Coordinates(xyze_other);

  Teuchos::RCP<GEO::CUT::IntersectionBase> inter_ptr = IntersectionPtr(other->Shape());
  // other is line element, this is surface
  inter_ptr->Init(xyze_other, xyze_this, false, false, false, &(mesh->GetOptions()));

  bool edges_parallel =
      HandleParallelCut(other, side, cut_points, mesh->GetOptions().GeomDistance_Floattype());

  if (not edges_parallel)
  {
    // perform the actual edge - edge intersection
    const enum IntersectionStatus istatus =
        inter_ptr->ComputeEdgeSideIntersection(tolerance, true, NULL);

    enum IntersectionStatus retstatus = istatus;

    switch (istatus)
    {
      case intersect_single_cut_point:
      {
        double* x_ptr = NULL;

        x_ptr = inter_ptr->FinalPointEdgeEdge();

        const double& pos = inter_ptr->LocalCoordinates()[0];

        Point* p = NULL;
        p = Point::NewPoint(*mesh, x_ptr, pos, this, side, tolerance);

        p->AddEdge(this);
        // NOTE: Might cause errors in the future, since local coordinate of the point on other edge
        // is based on the merged point, not original one
        p->AddEdge(other);
        // might remove this one, as there is already one in cut_edge.cpp::FindCutPoints cln
        p->AddPair(side, this);

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

        p->AddEdgeIntersection(this, other, side, this, msg_extra.str());

        break;
      }
      case intersect_multiple_cut_points:
      {
        // happens when one edge lies on another

        std::vector<LINALG::Matrix<probDim, 1>> xyz_cuts;
        std::vector<LINALG::Matrix<dimEdge, 1>> r_cuts;
        inter_ptr->FinalPoints(xyz_cuts);
        inter_ptr->LocalSideCoordinates(r_cuts);

        dsassert(xyz_cuts.size() == r_cuts.size(), "Size mismatch!");

        Point* p = NULL;
        for (unsigned i = 0; i < xyz_cuts.size(); ++i)
        {
          p = Point::NewPoint(*mesh, xyz_cuts[i].A(), r_cuts[i](0), this, side, tolerance);

          p->AddEdge(this);
          p->AddEdge(other);
          // might remove this one, as there is already one in cut_edge.cpp::FindCutPoints cln
          p->AddPair(side, this);
          p->AddEdgeIntersection(this, other, side, this, "MULTIPLE CUT POINTS");

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
        dserror("Unsupported intersection status! ( status = %d )", istatus);
        exit(EXIT_FAILURE);
    }
    return (static_cast<bool>(retstatus));
  }
  else
    return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probDim, DRT::Element::DiscretizationType edgeType, unsigned dimEdge,
    unsigned numNodesEdge>
Teuchos::RCP<GEO::CUT::IntersectionBase>
GEO::CUT::ConcreteEdge<probDim, edgeType, dimEdge, numNodesEdge>::IntersectionPtr(
    const DRT::Element::DiscretizationType& sidetype) const
{
  return GEO::CUT::IntersectionBase::Create(edgeType, sidetype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Edge> GEO::CUT::EdgeFactory::CreateEdge(
    const DRT::Element::DiscretizationType& edgetype, const std::vector<Node*>& nodes) const
{
  Teuchos::RCP<Edge> cedge_ptr = Teuchos::null;
  const int probdim = DRT::Problem::Instance()->NDim();
  switch (edgetype)
  {
    case DRT::Element::line2:
    {
      cedge_ptr = Teuchos::rcp(CreateConcreteEdge<DRT::Element::line2>(nodes, probdim));
      break;
    }
    default:
    {
      dserror(
          "Unsupported edge type! ( %d | %s )", edgetype, DRT::DistypeToString(edgetype).c_str());
      break;
    }
  }
  return cedge_ptr;
}

template class GEO::CUT::ConcreteEdge<2, DRT::Element::line2>;
template class GEO::CUT::ConcreteEdge<3, DRT::Element::line2>;
