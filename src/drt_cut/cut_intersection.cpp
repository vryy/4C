/*----------------------------------------------------------------------*/
/*!
\file cut_intersection.cpp

\brief here the intersection of a (plane) surface with a line is performed

\level 2

\maintainer Christoph Ager
*/
/*----------------------------------------------------------------------*/

#include "cut_intersection.H"
#include "cut_boundingbox.H"
#include "cut_side.H"
#include "cut_position.H"
#include "cut_utils.H"

#include "../drt_lib/drt_globalproblem.H"

// whether to perform triangulation during intersection
#define TRIANGULATED_INTERSECTION true

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::IntersectionBase> GEO::CUT::IntersectionBase::Create(
    const DRT::Element::DiscretizationType& edge_type,
    const DRT::Element::DiscretizationType& side_type)
{
  const IntersectionFactory factory;
  return factory.CreateIntersection(edge_type, side_type);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
void GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::SetCoordinates()
{
  GetEdge().Coordinates(xyze_lineElement_);
  GetSide().Coordinates(xyze_surfaceElement_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::ComputeCut(Edge* sedge, Edge* eedge, Side* side, PointSet& ee_cut_points,
    double& tolerance)
{
  return eedge->ComputeCut(GetMeshPtr(), sedge, side, &ee_cut_points, tolerance);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
void GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::GetConnectivityInfo(const LINALG::Matrix<probdim, 1>& xreal,
    const std::vector<int>& touched_edges_ids, std::set<std::pair<Side*, Edge*>>& out)
{
  const std::vector<Edge*> side_edges = GetSide().Edges();
  for (std::vector<int>::const_iterator i = touched_edges_ids.begin(); i != touched_edges_ids.end();
       ++i)
  {
    // mapping indeces to real edges
    Edge* e = side_edges[*i];

    // if this edge cuts other edge in this point, this edge cuts all sides of the edge there
    const plain_side_set& touched_edge_sides = e->Sides();
    for (plain_side_set::const_iterator j = touched_edge_sides.begin();
         j != touched_edge_sides.end(); ++j)
    {
      Side* s = *j;
      out.insert(std::make_pair(s, &GetEdge()));
    }

    // also add connection of < touched edge x sides of the cut_edge >
    const plain_side_set& cut_edge_sides = GetEdge().Sides();
    for (plain_side_set::const_iterator j = cut_edge_sides.begin(); j != cut_edge_sides.end(); ++j)
    {
      Side* s = *j;
      out.insert(std::make_pair(s, e));
    }
  }
}

template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
void GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::AddConnectivityInfo(Point* p, const LINALG::Matrix<probdim, 1>& xreal,
    const std::vector<int>& touched_edges_ids,
    const std::set<std::pair<Side*, Edge*>>& touched_cut_pairs)
{
  // original intersection pair
  std::pair<Side*, Edge*> or_int_pair = std::make_pair(&GetSide(), &GetEdge());
  const std::vector<Edge*>& side_edges = GetSide().Edges();
  std::set<Edge*> touched_edges;

  for (std::vector<int>::const_iterator i = touched_edges_ids.begin(); i != touched_edges_ids.end();
       ++i)
  {
    // mapping indeces to real edges
    Edge* e = side_edges[*i];
    touched_edges.insert(e);

    // NOTE: This might cause some problems in the future, due to the way p->t is used
    if (probdim == 3)
    {
      LINALG::Matrix<3, 1> A;
      for (unsigned int cnt = 0; cnt < probdim; ++cnt) A(cnt, 0) = xreal(cnt, 0);
      p->t(e, A);
    }
    else
    {
      dserror("Not supported");
    }
    p->AddEdge(e);
  }

  for (std::set<std::pair<Side*, Edge*>>::const_iterator it = touched_cut_pairs.begin();
       it != touched_cut_pairs.end(); ++it)
  {
#if CUT_CREATION_INFO
    std::stringstream msg;
    if (touched_edges.find((*it).second) != touched_edges.end())
      msg << "Pair added because during intersection of side" << (or_int_pair.first)->Id()
          << " and edge " << (or_int_pair.second)->Id() << " edge " << (*it).second->Id()
          << " is a touching edge. Hence it cuts all sides of cut edge";
    else
      msg << "Pair added because during intersection of side" << (or_int_pair.first)->Id()
          << " and edge " << (or_int_pair.second)->Id() << " side " << (*it).first->Id()
          << " is a side of touching edge. Hence cut_edge cuts it";
    p->AddCreationInfo(*it, msg.str());
#endif
    //  p->AddEdge((*it).second);
    p->AddSide((*it).first);
    p->AddPair((*it).first, (*it).second, or_int_pair);
  }
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckBoundingBoxOverlap()
{
  Teuchos::RCP<BoundingBox> sbb = Teuchos::rcp(BoundingBox::Create(GetSide()));
  Teuchos::RCP<BoundingBox> ebb = Teuchos::rcp(BoundingBox::Create(GetEdge()));

  return CheckBoundingBoxOverlap(*sbb, *ebb);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckBoundingBoxOverlap(BoundingBox& ebb, BoundingBox& sbb) const
{
  return (not sbb.Within(1.0, ebb));
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckParallelism(std::vector<LINALG::Matrix<dimside, 1>>& side_rs_intersect,
    std::vector<LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance)
{
  switch (dimside)
  {
    case 1:
    {
      // NOTE: Edge-Edge parallelism is directly handled during ComputeEdgeIntersection
      return false;
    }
    case 2:
    {
      return CheckParallelismBetweenSideAndEdge(side_rs_intersect, edge_r_intersect, tolerance);
    }
    default:
    {
      dserror(
          "The given side element type is currently unsupported! \n"
          "( dim = %d | sideType = %s ",
          dimside, DRT::DistypeToString(sidetype).c_str());
      exit(EXIT_FAILURE);
    }
  }
  // this cannot be reached
  exit(EXIT_FAILURE);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckParallelismBetweenSideAndEdge(std::vector<LINALG::Matrix<dimside, 1>>&
                                                          side_rs_intersect,
    std::vector<LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance)
{
  side_rs_intersect.clear();
  edge_r_intersect.clear();
  tolerance = 0.0;

  //  if ( CheckAngleCriterionBetweenSideNormalAndEdge() )
  //    dserror("This shouldn't happen at this point! -- hiermeier");

  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckAngleCriterionBetweenSideNormalAndEdge()
{
  // calculate the normal at the side element center
  const LINALG::Matrix<dimside, 1> rst_side_center(
      DRT::UTILS::getLocalCenterPosition<dimside>(sidetype));

  LINALG::Matrix<probdim, numNodesSide> side_deriv1;
  LINALG::Matrix<probdim, probdim> xjm;
  LINALG::Matrix<probdim, 1> normal_center;
  GEO::CUT::EvalDerivsInParameterSpace<probdim, sidetype, double>(
      xyze_surfaceElement_, rst_side_center, side_deriv1, xjm, NULL, &normal_center, NULL, true);

  // calcualte the direction vector of the edge element
  LINALG::Matrix<probdim, 1> dedge(&xyze_lineElement_(0, 0), false);
  const LINALG::Matrix<probdim, 1> e_endpoint(&xyze_lineElement_(0, 1), true);
  dedge.Update(-1.0, e_endpoint, 1.0);
  const double e_nrm2 = dedge.Norm2();
  if (e_nrm2 == 0.0)
    dserror("The 1-st edge length is zero!");
  else
    dedge.Scale(1.0 / e_nrm2);

  // calculate the inner product and check the angle between the normal and
  // the edge

  const double inner_product = dedge.Dot(normal_center);
  /* If the angle between the normal and the edge is smaller than 89°
   * ( cos( 89° ) = 0.017472406... ), the vectors are definitely not parallel,
   * otherwise it's possible and we do another test. Note, that both vectors have
   * unit length! */
  {
    double a = std::acos(std::abs(inner_product)) * 180 / (2.0 * std::acos(0.0));
    std::cout << "angle between side normal and edge = " << a << "°" << std::endl;
  }
  return (std::abs(inner_product) <= 0.01747);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckCollinearity(std::vector<LINALG::Matrix<dimside, 1>>&
                                         side_rs_corner_intersect,
    std::vector<LINALG::Matrix<dimedge, 1>>& edge_r_corner_intersect, double& tolerance)
{
  dserror("Is this used?");
  if (numNodesEdge != 2 or numNodesSide != 2)
    dserror(
        "Two line2 elements are expected, but instead a %s (edge) and %s (side) "
        "element were given.",
        DRT::DistypeToString(edgetype).c_str(), DRT::DistypeToString(sidetype).c_str());

  side_rs_corner_intersect.clear();
  edge_r_corner_intersect.clear();


  // quick check first -- check angle between the two edges, if the check fails,
  // the two edges are definitely not parallel.
  if (not CheckAngleCriterionBetweenTwoEdges()) return false;

  /* If the lines seem to be parallel, we do the expensive calculation and
   * check if the end points of the edge are on the side edge ( and are collinear ). */
  bool is_collinear = true;

  std::vector<LINALG::Matrix<probdim, 1>> side_xyz_corner_intersect;
  side_xyz_corner_intersect.reserve(2);

  side_rs_corner_intersect.reserve(2);
  edge_r_corner_intersect.reserve(2);

  LINALG::TMatrix<LINALG::Matrix<2, 1>, 2, 1> e_corner_distance;
  for (unsigned i = 0; i < 2; ++i)
  {
    const LINALG::Matrix<probdim, 1> e_cornerpoint(&xyze_lineElement_(0, i), true);
    Teuchos::RCP<GEO::CUT::Position> pos =
        GEO::CUT::Position::Create(xyze_surfaceElement_, e_cornerpoint, sidetype);

    bool withinlimits = pos->Compute(true);

    pos->Distance(e_corner_distance(i));
    const double pos_tolerance = pos->NewtonTolerance();

    // keep the largest up-coming tolerance
    if (i == 0 or tolerance < pos_tolerance) tolerance = pos_tolerance;

    if (e_corner_distance(i).Norm2() < pos_tolerance)
    {
      // the edge end point seems to lie on the side edge, and the end-point is
      // also within the side limits
      if (withinlimits)
      {
        side_xyz_corner_intersect.push_back(e_cornerpoint);
        side_rs_corner_intersect.push_back(LINALG::Matrix<dimside, 1>(true));
        pos->LocalCoordinates(*(side_rs_corner_intersect.end() - 1));
      }
    }
    // if one distance is larger than the given tolerance, the edge cannot be
    // collinear
    else
      is_collinear = false;
  }
  // set number of cut points
  num_cut_points_ = side_rs_corner_intersect.size();

  // The two edges are collinear and we found no end points of the edge on
  // the side.
  // --> No more intersections possible for two edges besides the already
  // existing end-points of the side edge.
  if (is_collinear and side_rs_corner_intersect.size() == 0)
  {
    return true;
  }

  // End points of the edge on the side edge. In this case
  // we are done, as well, and it doesn't matter if the two objects
  // are parallel or not.
  if (num_cut_points_ > 0)
  {
    for (unsigned i = 0; i < num_cut_points_; ++i)
    {
      xsi_side_.Update(side_rs_corner_intersect[0]);
      if (not FindLocalCoordinateOfEdgeEndPoint(
              xsi_edge_(0), side_xyz_corner_intersect[i], tolerance))
        dserror("We couldn't find the correct edge end-point!");
      if (num_cut_points_ > 1)
        edge_r_corner_intersect.push_back(LINALG::Matrix<dimedge, 1>(xsi_edge_.A(), false));

      // safety
      if (not(LineWithinLimits() and SurfaceWithinLimits())) dserror("Something went wrong!");
    }

    return true;
  }

  if (not is_collinear)
  {
    double prod = 0.0;
    // check the signs of the distance values
    for (unsigned i = 0; i < (probdim - dimside); ++i)
    {
      prod = e_corner_distance(0)(i) * e_corner_distance(1)(i);
      // if the off-set direction is exactly in one of the normal directions, we
      // use the second criterion to detect numerical artifacts ( values close
      // to zero but different signs )
      if ((prod < 0.0) and std::abs(e_corner_distance(0)(i)) > tolerance and
          std::abs(e_corner_distance(1)(i)) > tolerance)
      {
        // there is a intersection possible, although the lines are almost
        // parallel
        return false;
      }
    }
    // all distance values have the same sign ( or are close to zero ) and the
    // edges are not collinear
    // --> no intersection possible
    return true;
  }

  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::CheckAngleCriterionBetweenTwoEdges()
{
  LINALG::Matrix<probdim, 1> deedge(&xyze_lineElement_(0, 0), false);
  LINALG::Matrix<probdim, 1> dsedge(&xyze_surfaceElement_(0, 0), false);

  const LINALG::Matrix<probdim, 1> e_endpoint(&xyze_lineElement_(0, 1), true);
  deedge.Update(-1.0, e_endpoint, 1.0);
  const double e_nrm2 = deedge.Norm2();
  if (e_nrm2 == 0.0) dserror("The 1-st edge length is zero!");

  const LINALG::Matrix<probdim, 1> s_endpoint(&xyze_surfaceElement_(0, 1), true);
  dsedge.Update(-1.0, s_endpoint, 1.0);
  const double s_nrm2 = dsedge.Norm2();
  if (s_nrm2 == 0.0) dserror("The 2-nd edge length is zero!");

  const double inner_product = deedge.Dot(dsedge);
  /* If the angle between the two lines is larger than 1° ( cos( 1° ) = 0.99866695... ),
   * the two edges are definitely not parallel, otherwise it's possible and
   * we take a closer look. */
  if (debug)
  {
    std::cout << e_nrm2 << ", " << s_nrm2 << ", " << inner_product << std::endl;
    double a = std::abs(inner_product) / (e_nrm2 * s_nrm2);
    a = std::acos((a > 1.0 ? 1.0 : a)) * 180 / (2.0 * std::acos(0.0));
    std::cout << "angle between edge and edge = " << a << "°" << std::endl;
  }
  return (std::abs(inner_product) >= e_nrm2 * s_nrm2 * 0.9998);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::FindLocalCoordinateOfEdgeEndPoint(double& pos,
    const LINALG::Matrix<probdim, 1>& xyz, const double& tolerance) const
{
  const double r_endpoints[2] = {-1.0, 1.0};
  LINALG::Matrix<probdim, 1> dist;
  LINALG::Matrix<numNodesEdge, 1> lineFunct;

  for (unsigned i = 0; i < 2; ++i)
  {
    pos = r_endpoints[i];

    DRT::UTILS::shape_function_1D(lineFunct, pos, edgetype);
    dist.MultiplyNN(xyze_lineElement_, lineFunct);

    dist.Update(1.0, xyz, -1.0);

    if (dist.Norm2() < tolerance) return true;
  }
  return false;
}

template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
GEO::CUT::ParallelIntersectionStatus GEO::CUT::Intersection<probdim, edgetype, sidetype, debug,
    dimedge, dimside, numNodesEdge, numNodesSide>::HandleParallelIntersection(PointSet& cuts,
    int skip_id, bool output)
{
  Node* id_to_node[2] = {GetEdge().BeginNode(), GetEdge().EndNode()};
  int id_start = 0, id_end = 2;
  int count_inside = 0;  // number of points inside the surface

  bool found_intersections = false;

  bool zeroarea = false;
  // local coordinates are inside element
  LINALG::TMatrix<bool, 2, 1> lineendpoint_within_surfacelimits(true);
  // point is really inside element (normal distance = 0)
  LINALG::TMatrix<bool, 2, 1> lineendpoint_in_surface(true);
  LINALG::Matrix<2, 1> lineendpoint_dist;
  LINALG::Matrix<2, 1> lineendpoint_tol;
  LINALG::TMatrix<bool, 2, 1> lineendpoint_conv(true);
  LINALG::Matrix<dimside + dimedge, 2> lineendpoint_xsi(true);
  std::vector<std::vector<int>> lineendpoint_touched_edges(2);
  std::vector<KERNEL::PointOnSurfaceLoc> lineendpoint_location_kernel(2);

  // if one intersection point was found already for this cut_side & edge, and we are searching for
  // another one
  if (skip_id != -1)
  {
    if (skip_id == 0)
      id_start = 1;
    else if (skip_id == 1)
      id_end = 1;
    else
      dserror(
          "Trying to skip id = %d of the edge nodes. You can only skip BeginNode(0) or "
          "EndNode(1)\n",
          skip_id);

    lineendpoint_conv(skip_id) = true;
    lineendpoint_location_kernel[skip_id] = KERNEL::PointOnSurfaceLoc(true, true);
    count_inside++;
    found_intersections = true;
  }

  for (int lp = id_start; lp < id_end; ++lp)
  {
    GEO::CUT::Point* actpoint = id_to_node[lp]->point();

    switch (sidetype)
    {
      // --------------------------------------------------------------------
      /* handle parallel case for TRI3 elements in the 3-dimensional space
       * handle parallel case for 1-D elements in the 2-dimensional space */
      // --------------------------------------------------------------------
      case DRT::Element::tri3:
      case DRT::Element::line2:
      {
        lineendpoint_conv(lp) =
            ComputeDistance(actpoint, lineendpoint_dist(lp), lineendpoint_tol(lp), zeroarea,
                lineendpoint_location_kernel[lp], lineendpoint_touched_edges[lp], true);
        lineendpoint_within_surfacelimits(lp) = SurfaceWithinLimits();
        // copy xsi_ in the lp column of lineendpoint_xsi
        std::copy(xsi_.A(), xsi_.A() + (dimside + dimedge), &lineendpoint_xsi(0, lp));
        break;
      }  // end: case (DRT::Element::tri3) and 1-D elements
         // --------------------------------------------------------------------
         // handle parallel case for QUAD4 elements
         // --------------------------------------------------------------------
      case DRT::Element::quad4:
      {
        /* problem here is that ComputeDistance might change the normal direction if
         * the projected point lies outside the element, that why we prefer to compute
         * the distance onto tri3! */
        LINALG::Matrix<2, 1> tri_dist;
        LINALG::Matrix<2, 1> tri_tol;
        LINALG::TMatrix<bool, 2, 1> tri_conv;

        std::vector<KERNEL::PointOnSurfaceLoc> tri_location_kernel(2);
        std::vector<std::vector<int>> tri_touched_edges(2);

        bool extended_tri_tolerance_loc_triangle_split[2] = {false, false};
        bool on_side[2] = {false, false};
        bool within_side[2] = {false, false};

        unsigned tri = 0;
        // try to calculate by splitting quad4 intro 2 triangles
        for (tri = 0; tri < 2; ++tri)
        {
          tri_conv(tri) = ComputeDistance(actpoint, tri_dist(tri), tri_tol(tri), zeroarea,
              tri_location_kernel[tri], tri_touched_edges[tri], true, tri,
              extended_tri_tolerance_loc_triangle_split[tri]);

          // we prefer to get on_side data from the side the point is within
          if (tri_location_kernel[tri].WithinSide())
          {
            within_side[tri] = true;
            on_side[tri] = tri_location_kernel[tri].OnSide();
          }

          if (tri_location_kernel[tri].OnSide()) on_side[tri] = true;
        }

        // if it inside quad4
        if (within_side[0] or within_side[1])
        {
          // triangle in which point was inside will be served as reference
          tri = within_side[0] ? 0 : 1;

          // safety check to make sure obtained edges are different (which they should be, since
          // diagonal one is not added)
          std::set<int> touched_edges(tri_touched_edges[0].begin(), tri_touched_edges[0].end());
          std::copy(tri_touched_edges[1].begin(), tri_touched_edges[1].end(),
              std::inserter(touched_edges, touched_edges.begin()));
          if (touched_edges.size() != (tri_touched_edges[1].size() + tri_touched_edges[0].size()))
            dserror("Got duplicate edges from triangle split");

          // set up all the topology data
          std::vector<int> tri_touched_edges_full(touched_edges.begin(), touched_edges.end());
          lineendpoint_touched_edges[lp] = tri_touched_edges_full;
          lineendpoint_location_kernel[lp] = tri_location_kernel[tri];
          lineendpoint_conv(lp) = tri_conv(tri);
          lineendpoint_dist(lp) = tri_dist(tri);
          lineendpoint_tol(lp) = tri_tol(tri);
          lineendpoint_within_surfacelimits(lp) = tri_location_kernel[tri].WithinSide();

          std::copy(xsi_.A(), xsi_.A() + (dimside + dimedge), &lineendpoint_xsi(0, lp));
        }
        // it is outside for sure
        else
        {
          // first check if it is not close to the common edge of two triangles and we missed it
          // if the point is parellel to both tri, we needle to handle extreme case of point being
          // close to shared edge of both tri (diagonal of the quad)
          if ((tri_location_kernel[0].OnSide()) and (tri_location_kernel[1].OnSide()))
          {
            if (extended_tri_tolerance_loc_triangle_split[0] and
                extended_tri_tolerance_loc_triangle_split[1])
            {
              GenerateGmshDump();
              dserror(
                  "Here we need to create point on diagonal line between two triangles in "
                  "triangulated quad4! This is not yet implemented");
            }
          }

          // basically we can take any triangle now as a refernce, since point is outside of both

          // reset to previous sucessfull value
          if (on_side[0])
            tri = 0;
          else
            tri = 1;

          lineendpoint_location_kernel[lp] = tri_location_kernel[tri];
          lineendpoint_touched_edges[lp] = tri_touched_edges[tri];
          lineendpoint_conv(lp) = tri_conv(tri);
          lineendpoint_dist(lp) = tri_dist(tri);
          lineendpoint_tol(lp) = tri_tol(tri);
        }

        break;
      }  // end: case (DRT::Element::quad4)
      default:
      {
        dserror(
            "Intersection Error: Other surfaces than TRI3 and QUAD4 are not supported "
            "to avoid huge problems!!!");
        break;
      }
    }
  }


  if (!(lineendpoint_location_kernel[0].OnSide()) and
      !(lineendpoint_location_kernel[1].OnSide()) and
      lineendpoint_dist(0) * lineendpoint_dist(1) > 0 and not zeroarea)
  {
    if (debug)
    {
      std::cout << "No Cut Point, because we detected following normal distance "
                   "at the line end points: ... lineendpoint_dist[0]: "
                << lineendpoint_dist(0) << " lineendpoint_dist[1]: " << lineendpoint_dist(1)
                << "(Codeline:" << __LINE__ << ")" << std::endl;
    }
    return intersection_not_possible;
  }


  if (lineendpoint_conv(0) and lineendpoint_conv(1))
  {
    // ------------------------------------------------------------------------
    // starting- and end-point of the edge are inside/on the side
    // ------------------------------------------------------------------------

    if (lineendpoint_location_kernel[0].WithinSide() and lineendpoint_location_kernel[0].OnSide())
      lineendpoint_in_surface(0) = true;


    if (lineendpoint_location_kernel[1].WithinSide() and lineendpoint_location_kernel[1].OnSide())
      lineendpoint_in_surface(1) = true;

    for (int lp = id_start; lp < id_end; ++lp)
    {
      if (lineendpoint_in_surface(lp))
      {
        LINALG::Matrix<probdim, 1> x;
        id_to_node[lp]->point()->Coordinates(x.A());
        InsertCut(id_to_node[lp], cuts);

#if CUT_CREATION_INFO
        std::stringstream msg;
        msg << "// Added because point" << id_to_node[lp]->point()->Id() << " of the edge "
            << GetEdge().Id() << " lies on the side " << GetSide().Id() << std::endl;
        msg << "// Point " << id_to_node[lp]->point()->Id() << " has "
            << lineendpoint_touched_edges[lp].size() << " touching edges of the side " << std::endl;
        msg << "// Computed from ComputeDistance" << std::endl;
        if (not lineendpoint_within_surfacelimits(lp))
          msg << "// Is NOT within surface limits" << std::endl;
        else
          msg << "// Within surface limits" << std::endl;
        id_to_node[lp]->point()->AddCreationInfo(std::make_pair(&GetSide(), &GetEdge()), msg.str());
#endif
        // Adding connection that are close to this nodal point
        std::set<std::pair<Side*, Edge*>> topological_cut_pairs;
        GetConnectivityInfo(x, lineendpoint_touched_edges[lp], topological_cut_pairs);
        AddConnectivityInfo(
            id_to_node[lp]->point(), x, lineendpoint_touched_edges[lp], topological_cut_pairs);

        ++count_inside;
        found_intersections = true;
      }
    }

    // reset for looping over edges later
    found_intersections = false;

    // if both points inside, no need to do anything
    if (count_inside == 2)
    {
      return intersection_found;
    }

    else
    {
      // ------------------------------------------------------------------------
      /* Try to find intersection of the current edge with all edges of the side   */
      // ------------------------------------------------------------------------
      Teuchos::RCP<BoundingBox> edge_bb = Teuchos::rcp(BoundingBox::Create(GetEdge()));
      const std::vector<Edge*>& side_edges = GetSide().Edges();
      for (std::vector<Edge*>::const_iterator i = side_edges.begin(); i != side_edges.end(); ++i)
      {
        Edge* e = *i;
        Teuchos::RCP<BoundingBox> side_edge_bb = Teuchos::rcp(BoundingBox::Create(*e));
        ;
        // skip this edge if it is too far or intersection has been done already
        if (not side_edge_bb->Within(1.0, *edge_bb)) continue;

        PointSet cut_points;
        GetEdge().GetCutPoints(e, cut_points);
        // cut point obtained e.g. from other neigboring edge (e.g. it is end point) or from the
        // other intersection

        std::stringstream msg;
#if CUT_CREATION_INFO
        msg << "// Edge-edge intersection was happenning in intersection between edge "
            << GetEdge().Id() << " and "
            << "edges of the side" << GetSide().Id() << "( edge " << e->Id() << ")" << std::endl;
#endif
        if (cut_points.size() > 0)
        {
          if (cut_points.size() == 1)
          {
#if CUT_CREATION_INFO
            msg << "// Edge-edge intersection is added because point" << (*cut_points.begin())->Id()
                << " was already cut by connecting edge" << std::endl;
#endif
            /* Then in order for proper cut we need to add all the necessary connectivity
             * infromation */
            /* NOTE: Here it might be possible to better transfer original intersection from the
             * another pair */
            (*cut_points.begin())
                ->AddEdgeIntersection(&GetEdge(), e, &GetSide(), &GetEdge(), msg.str());
            /* Nothing to be done. There are already cut points between these
             * edges. We cannot find new ones. */
            cuts.insert(*(cut_points.begin()));
          }
          else if (cut_points.size() == 2)
          {
            PointSet::iterator cp = cut_points.begin();
            cuts.insert(*cp);
#if CUT_CREATION_INFO
            msg << "Edge-edge intersection is added because point" << (*cp)->Id()
                << " was already cut by connecting edge. In fact two common points were found!"
                << std::endl;
#endif
            (*cp)->AddEdgeIntersection(&GetEdge(), e, &GetSide(), &GetEdge(), msg.str());
            cp++;
            cuts.insert(*cp);
#if CUT_CREATION_INFO
            msg << "Edge-edge intersection is added because point" << (*cp)->Id()
                << " was already cut by connecting edge. In fact two common points were found!"
                << std::endl;
#endif
            (*cp)->AddEdgeIntersection(&GetEdge(), e, &GetSide(), &GetEdge(), msg.str());
          }
          else
          {
            if (debug)
            {
              for (PointSet::iterator i = cut_points.begin(); i != cut_points.end(); ++i)
                (*i)->Print(std::cout);
              GetEdge().BeginNode()->point()->Print(std::cout);
              GetEdge().EndNode()->point()->Print(std::cout);
              e->BeginNode()->point()->Print(std::cout);
              e->EndNode()->point()->Print(std::cout);
            }
            dserror("Two Edges have more than two cutpoint, sounds strange!");
          }
          if (debug)
          {
            std::cout << "Cut points found, by intersection of edges!"
                      << " (Codeline:" << __LINE__ << ")" << std::endl;
          }
          found_intersections = true;
        }  // if ( cut_points.size() > 0 )
        // given edge has currently no cut points
        else
        {
          double tolerance = 0.0;
          // check if the given edge intersects with one of the side edges
          if (ComputeCut(e, GetEdgePtr(), GetSidePtr(), cuts, tolerance))
          {
            if (debug)
            {
              std::cout << "Cut points found, by intersection of edges (edges of side "
                           "and given edge)!"
                        << " (Codeline:" << __LINE__ << ")" << std::endl;
            }
            found_intersections = true;
          }
        }
      }  // end: loop over side edges

      if (found_intersections)
      {
        return intersection_found;
      }

      // determining return status based compute distance data, that we calculated earlier
      else
      {
        // if both nodes of the edge are on the side and there was no intersection till now, it
        // means that they lie outside and paralell to the side --> Newton method for the
        // intersection of cut_side x cut_edge  will not converge so we don't want to continue from
        // here
        /* both points outside the plane on the same side --> no intersection
         * Note: distances of both end-points have the same sign. */

        if ((lineendpoint_location_kernel[0].OnSide()) and
            (lineendpoint_location_kernel[1].OnSide()))
        {
          if ((not(lineendpoint_location_kernel[0].WithinSide())) and
              (not(lineendpoint_location_kernel[1].WithinSide())))
          {
            return intersection_not_possible;
          }
          // if at least one point is inside, there are two options:
          // 1. Both are inside, and this should have been definitely handled before ( count_inside
          // == 2 )
          // 2. One is inside and parallel to one tri3 split, one is outside and parallel to another
          // tri3
          //    only in such case no edge-edge intersection would hapen. The quad4 is need to be
          //    slightly distorted (e.g. even 1e-16) Otherwise there would be cut point from edge
          //    edge intersection
          else
          {
            std::pair<bool, bool> is_distorted = IsQuad4Distorted();
            if (not(is_distorted.first or is_distorted.second))
            {
              GenerateGmshDump();
              dserror("This case should have been handled before with edge-edge intersection!");
            }
            else
              return intersection_found;
          }
        }
        // If one point is inside and there is no intersection with edges , we are sure that nothing
        // else can be done. Otherwise there is still possibility for a cut point
        else if (count_inside == 1)
        {
          return intersection_found;
        }
        // if one point lies on the extended side, but there was no edge-edge intersection, there is
        // no other possible intersection
        else if ((lineendpoint_location_kernel[0].OnSide()) or
                 (lineendpoint_location_kernel[1].OnSide()))
        {
          return intersection_not_possible;
        }

        else
          return intersection_not_found;
      }
    }
  }

  else if (zeroarea)
  {
    if (debug)
    {
      std::cout << "No Cut Point, because surface has no area and edges don't have "
                   "any intersection! (Codeline:"
                << __LINE__ << ")" << std::endl;
    }
    return intersection_not_found;
  }

  return intersection_not_found;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
void GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::GenerateGmshDump()
{
  // just to create gmsh output of the failed intersection!!!
  try
  {
    // compute intersection with no cln, just to make sure we don't fail
    KERNEL::ComputeIntersection<probdim, edgetype, sidetype, false> ci(xsi_);
    ci(xyze_surfaceElement_, xyze_lineElement_);
    std::string filename(OUTPUT::GenerateGmshOutputFilename(".intersection_CUTFAIL_end.pos"));
    std::ofstream file(filename.c_str());
    ci.WritetoGmsh(file);
    file.close();
  }
  catch (std::runtime_error& e)
  {
    dserror("Cautch error in the cut_intersection:  \n%s . Current tolerance must be increased",
        e.what());
  }
}

template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::HandleSpecialCases()
{
  LINALG::Matrix<dimside + dimedge, 1> i_xsi;
  i_xsi = xsi_;

  /* point is outside the element and is not part of the interpolation
   * space (just for QUAD4 try triangulation) */
  if (sidetype == DRT::Element::quad4)
  {
    LINALG::Matrix<2, 1> tri_conv;
    KERNEL::PointOnSurfaceLoc location_status[2];
    for (unsigned tri = 0; tri < 2; ++tri)
    {
      tri_conv(tri) = ComputeEdgeTri3Intersection(tri, location_status[tri]);
      /* we expect the QUAD4 to converge in case that the projected point is
       * inside the QUAD4 or TRI3! */
      if (location_status[tri].WithinSide())
      {
        GenerateGmshDump();
        bool is_quad4_distorted = IsQuad4Distorted().first;
        if (is_quad4_distorted) std::cout << "NOTICE: Element is distorted" << std::endl;

        const std::vector<Node*>& side_nodes = GetSide().Nodes();
        for (std::vector<Node*>::const_iterator it = side_nodes.begin(); it != side_nodes.end();
             ++it)
        {
          (*it)->point()->DumpConnectivityInfo();
        }
        dserror(
            "ComputeEdgeSideIntersection for Quad4 didn't converge, "
            "but ComputeEdgeTri3Intersection for triangulation (id=%d) is inside the Element!",
            tri);
      }
    }
    // both converged and are outside!!! ( tested directly )
    if (tri_conv(0) and tri_conv(1))
    {
      if (debug)
      {
        std::cout << "No Cut Point found by intersection with triangulated quad4! "
                     "(Codeline:"
                  << __LINE__ << ")" << std::endl;
      }
      // point outside the interpolation space of QUAD4
      return false;
    }
  }
  // If we reach here, this is probably distorted  element
  dserror("Distorted element. To enable support for it switch TRIANGULATED_INTERSECTION to true");
  return false;
}


template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::TriangulatedIntersection(PointSet& cuts)
{
  double itol = 0.0;
  switch (sidetype)
  {
    case DRT::Element::tri3:
    case DRT::Element::line2:
    {
      IntersectionStatus conv;
      try
      {
        conv = ComputeEdgeSideIntersection(itol, true);
      }
      catch (std::runtime_error& e)
      {
        GenerateGmshDump();
        dserror("Cautch error in cut kernel. Current tolerance must be increased! Error is: \n %s ",
            e.what());
      }

      // if newton did not fail, we trust the result
      if (conv != intersect_newton_failed)
      {
        // if the point is inside ( or inside tolerance)
        if (conv == intersect_single_cut_point)
        {
          conv = intersect_single_cut_point;
          FinalPoint();

          Point* p = Point::NewPoint(
              GetMesh(), x_.A(), xsi_(dimside, 0), GetEdgePtr(), GetSidePtr(), itol);

          cuts.insert(p);
          return true;
        }
        // no intersection
        else
          return false;
      }
      else
        dserror(
            "Newton did not converge for edge-side(tri3) intersection and there is not cut point!");
      break;
    }
    case DRT::Element::quad4:
    {
      std::vector<LINALG::Matrix<3, 1>> final_points;
      std::vector<LINALG::Matrix<dimedge, 1>> edge_coords;
      LINALG::Matrix<dimside + dimedge, 1> i_xsi;
      bool close_to_shared_edge[2] = {false, false};
      i_xsi = xsi_;
      IntersectionStatus tri_status[2];
      // compute intersection for each of the triangles
      int cut_point_count = 0;
      bool first_close_to_shared = false;
      for (unsigned tri = 0; tri < 2; ++tri)
      {
        try
        {
          tri_status[tri] = ComputeEdgeTri3IntersectionQuad4Split(tri, close_to_shared_edge + tri);
        }
        catch (std::runtime_error& e)
        {
          GenerateGmshDump();
          dserror(
              "Cautch error in cut kernel. Current tolerance must be increased! Error is: \n %s ",
              e.what());
        }
        if (tri_status[tri] != intersect_newton_failed)
        {
          if (tri_status[tri] == intersect_single_cut_point)
          {
            istatus_ = intersect_single_cut_point;
            final_points.push_back(LINALG::Matrix<3, 1>(FinalPoint()));
            edge_coords.push_back(LINALG::Matrix<dimedge, 1>(xsi_edge_.A()));
            // we found cut point
            cut_point_count++;
          }
          // if have not obtained a definite point already and reasonably close to shared edge
          else if ((cut_point_count == 0) and (close_to_shared_edge[tri]))
          {
            // if first already satisfy condition and we are on tri3, and can safely create point
            if (first_close_to_shared)
            {
              if (tri == 1)
              {
                GenerateGmshDump();
                dserror(
                    "Here we should create point on diagonal line between two triangles in "
                    "triangulated quad4! This is not yet implemented!");
                // possible reference implementation
                istatus_ = intersect_single_cut_point;
                final_points.push_back(LINALG::Matrix<3, 1>(FinalPoint()));
                edge_coords.push_back(LINALG::Matrix<dimedge, 1>(xsi_edge_.A()));
                // we found cut point
                cut_point_count++;
              }
              else
                dserror("This should not happen");
            }
            // if tri equal to 1 we just let the check go to the next triangle
            if (tri == 0) first_close_to_shared = true;
          }
        }
      }

      switch (cut_point_count)
      {
        case 0:
        {
          // no intersection
          return false;
          break;
        }
        case 1:
        {
          // intersection
          Point* p = Point::NewPoint(GetMesh(), final_points[0].A(), edge_coords[0](0, 0),
              GetEdgePtr(), GetSidePtr(), itol);
          cuts.insert(p);
          // update status, in case of the other triangle returning different one
          istatus_ = intersect_single_cut_point;
          return true;
          break;
        }
        case 2:
        {
          // fisrt check if this point is not on the split between triangles (diagonal)
          LINALG::Matrix<3, 1> diff = final_points[0];
          diff.Update(-1, final_points[1], 1);
          // this case is unsupported
          if (diff.Norm2() > TOPOLOGICAL_TOLERANCE)
          {
            std::stringstream err_msg;
            err_msg << " Go two intersection points during intersection of triangulated quad4.  "
                       "This case is not supported! Global coordinates are \n";
            for (std::vector<LINALG::Matrix<3, 1>>::iterator it = final_points.begin();
                 it != final_points.end(); ++it)
            {
              err_msg << (*it);
            }
            err_msg << "Local coordinates on the edge are \n";
            for (typename std::vector<LINALG::Matrix<dimedge, 1>>::iterator it =
                     edge_coords.begin();
                 it != edge_coords.end(); ++it)
            {
              err_msg << (*it);
            }
            GenerateGmshDump();

            throw std::runtime_error(err_msg.str());
          }
          else
          {
            Point* p = Point::NewPoint(GetMesh(), final_points[0].A(), edge_coords[0](0, 0),
                GetEdgePtr(), GetSidePtr(), itol);
            cuts.insert(p);
            return true;
          }
          break;
        }
        default:
          dserror("This should not be possible!");
      }
      break;
    }
    default:
    {
      dserror(
          "Intersection Error: Other surfaces than TRI3 and QUAD4 are not supported "
          "to avoid huge problems!!!");
      break;
    }
  }
  return false;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::Intersect(PointSet& cuts)
{
  CheckInit();
  // ------------------------------------------------------------------------
  //  bounding box check
  // ------------------------------------------------------------------------
  bool debugbb = CheckBoundingBoxOverlap();
  if (debugbb)
  {
    if (debug)
    {
      std::cout << "No Cut Point, as not even the Bounding Boxes are "
                   "Overlapping! ( Codeline:"
                << __LINE__ << " )" << std::endl;
    }

    if (UseBoundingBox())
    {
      return false; /* not even the bounding boxes are overlapping, no need to
                     * search for an intersection point! */
    }
    else
      dserror("We should use bounding box now!");
  }

  // try to find out if this intersection can be handled just by ComputeDistance
  // and edge-edge intersections
  ParallelIntersectionStatus parallel_res = HandleParallelIntersection(cuts);

  if (parallel_res == intersection_found)
    return true;

  else
  {
    // we don't want to continue then, because  newton will fail
    if (parallel_res == intersection_not_possible) return false;
      // ------------------------------------------------------------------------
      // Try to calculate the intersection point directly with Newton
      // ------------------------------------------------------------------------
#if (TRIANGULATED_INTERSECTION)
    return TriangulatedIntersection(cuts);
#else

    double itol = 0.0;
    IntersectionStatus conv;
    try
    {
      conv = ComputeEdgeSideIntersection(itol, true);
    }
    catch (const std::runtime_error& e)
    {
      GenerateGmshDump();
      dserror("Cautch error in cut kernel. Current tolerance must be increased! Error is: \n %s ",
          e.what());
    }

    // if newton did not fail, we trust the result
    if (conv != intersect_newton_failed)
    {
      // if the point is inside ( or inside tolerance)
      if (conv == intersect_single_cut_point)
      {
        conv = intersect_single_cut_point;
        FinalPoint();

        Point* p =
            Point::NewPoint(GetMesh(), x_.A(), xsi_(dimside, 0), GetEdgePtr(), GetSidePtr(), itol);

#if CUT_CREATION_INFO
        std::stringstream msg;
        msg << "//Added because of the direct side " << GetSide().Id() << " and " << GetEdge().Id()
            << " intersection";
        p->AddCreationInfo(std::make_pair(&GetSide(), &GetEdge()), msg.str());
#endif

        cuts.insert(p);
        if (debug)
        {
          std::cout << "Cut points found by intersection, with local coords: xsi: " << xsi_(0, 0)
                    << " / eta: " << xsi_(1, 0) << " / alpha: " << xsi_(2, 0)
                    << " (Codeline:" << __LINE__ << ")" << std::endl;
        }
        if (debugbb)
        {
          std::cout
              << "\n ==== \n Bounding boxes said that there is no intersection point! \n ==== \n"
              << std::endl;
          GenerateGmshDump();
          dserror("Bounding Boxes said that there is no intersection point!");
        }
        return true;
      }
      else
      {
        if (debug)
        {
          std::cout << "No Cut Point found by intersection, with local coords: xsi: " << xsi_(0, 0)
                    << " / eta: " << xsi_(1, 0) << " / alpha: " << xsi_(2, 0)
                    << " and a tolerance of " << itol << " (Codeline:" << __LINE__ << ")"
                    << std::endl;
        }
        // there is no intersection between surface and line!!!
        return false;
      }
    }
    else
    {
      // check if point is not outside quad4 by splitting quad into 2 triangles and calculating
      // distances
      return HandleSpecialCases();
    }
#endif
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
bool GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside, numNodesEdge,
    numNodesSide>::RefinedBBOverlapCheck(int maxstep)
{
  if (sidetype != DRT::Element::tri3) dserror("RefinedBBOverlapCheck is made for distored tri3s!");

  std::vector<LINALG::Matrix<3, 1>> surfpoints;
  std::vector<LINALG::Matrix<3, 1>> linepoints;
  LINALG::Matrix<3, 1> actpoint;
  for (unsigned nid = 0; nid < numNodesSide; nid++)
  {
    (GetSide().Nodes()[nid])->Coordinates(actpoint.A());
    surfpoints.push_back(actpoint);
  }
  // find shortest edge on the surface!
  int smallestedge;
  {
    LINALG::Matrix<3, 1> d1, d2, d3;
    d1.Update(1.0, surfpoints[1], -1.0, surfpoints[2]);
    d2.Update(1.0, surfpoints[2], -1.0, surfpoints[0]);
    d3.Update(1.0, surfpoints[0], -1.0, surfpoints[1]);
    std::vector<double> lengths;
    lengths.push_back(d1.Norm2());
    lengths.push_back(d2.Norm2());
    lengths.push_back(d3.Norm2());
    if (lengths[0] < lengths[1] && lengths[0] < lengths[2])  // 0 is smallest
      smallestedge = 0;
    else if (lengths[1] < lengths[2] && lengths[1] < lengths[0])  // 1 is smallest
      smallestedge = 1;
    else  // 2 is smallest plus all the rest
      smallestedge = 2;
  }

  LINALG::Matrix<3, 1> v1(true);
  LINALG::Matrix<3, 1> v2(true);
  v1.Update(1.0, surfpoints[(smallestedge + 1) % 3], -1.0, surfpoints[smallestedge]);
  v2.Update(1.0, surfpoints[(smallestedge + 2) % 3], -1.0, surfpoints[smallestedge]);

  for (unsigned nid = 0; nid < numNodesEdge; nid++)
  {
    (GetEdge().Nodes()[nid])->Coordinates(actpoint.A());
    linepoints.push_back(actpoint);
  }
  LINALG::Matrix<3, 1> v3(true);
  v3.Update(1.0, linepoints[1], -1.0, linepoints[0]);

  bool overlap = false;  // we have an overlap of the refined bounding boxes?

  int pow_fac = 4;
  int max_steps = pow(pow_fac, maxstep - 1);
  double dmax_steps = 1.0 / max_steps;
  int act_steps = 0;
  uint act_boxidx = 0;
  std::vector<std::vector<int>> overlappingidx;
  std::vector<std::vector<int>> newoverlappingidx;
  std::vector<int> tmpoverlappingidx;
  tmpoverlappingidx.push_back(0);          // min idx with overlap for tri3
  tmpoverlappingidx.push_back(max_steps);  // max idx with overlap tri3
  tmpoverlappingidx.push_back(0);          // min idx with overlap for line2
  tmpoverlappingidx.push_back(max_steps);  // max idx with overlap for line2
  overlappingidx.push_back(tmpoverlappingidx);

  LINALG::Matrix<3, 1> p1, p2, p3, p4, lp1, lp2;

  for (int refinestep = 1; refinestep <= maxstep; refinestep++)
  {
    newoverlappingidx.clear();
    act_boxidx = 0;
    act_steps = pow(pow_fac, (maxstep - refinestep));

    overlap = false;
    for (int surfstep = 0; surfstep < max_steps; surfstep += act_steps)
    {
      tmpoverlappingidx.clear();
      while ((act_boxidx + 1) < overlappingidx.size() and
             (overlappingidx[act_boxidx + 1])[0] < surfstep)
        act_boxidx++;
      if (!(overlappingidx[act_boxidx][0] <= surfstep and
              overlappingidx[act_boxidx][1] >= surfstep + act_steps))
        continue;
      double alpha = surfstep * dmax_steps;
      double alphap = (surfstep + act_steps) * dmax_steps;
      p1.Update(1.0, surfpoints[smallestedge], alpha, v1);
      p2.Update(1.0, surfpoints[smallestedge], alphap, v1);
      p3.Update(1.0, surfpoints[smallestedge], alpha, v2);
      p4.Update(1.0, surfpoints[smallestedge], alphap, v2);
      Teuchos::RCP<BoundingBox> sbb = Teuchos::rcp(BoundingBox::Create());
      sbb->AddPoint(p1);
      sbb->AddPoint(p2);
      sbb->AddPoint(p3);
      sbb->AddPoint(p4);
      for (int linestep = 0; linestep < max_steps; linestep += act_steps)
      {
        if (!(overlappingidx[act_boxidx][2] <= linestep and
                overlappingidx[act_boxidx][3] >= linestep + act_steps))
          continue;
        double lalpha = linestep * dmax_steps;
        double lalphap = (linestep + act_steps) * dmax_steps;
        lp1.Update(1.0, linepoints[0], lalpha, v3);
        lp2.Update(1.0, linepoints[0], lalphap, v3);
        Teuchos::RCP<BoundingBox> ebb = Teuchos::rcp(BoundingBox::Create());
        ebb->AddPoint(lp1);
        ebb->AddPoint(lp2);

        if (sbb->Within(POSITIONTOL / BOXOVERLAP, *ebb))
        {
          overlap = true;
          if (!tmpoverlappingidx.size())
          {
            tmpoverlappingidx.push_back(surfstep);              // min idx with overlap for tri3
            tmpoverlappingidx.push_back(surfstep + act_steps);  // max idx with overlap for tri3
            tmpoverlappingidx.push_back(linestep);              // min idx with overlap for line2
            tmpoverlappingidx.push_back(linestep + act_steps);  // max idx with overlap for line2
          }
          else
          {
            if (tmpoverlappingidx[2] > linestep) tmpoverlappingidx[2] = linestep;
            if (tmpoverlappingidx[3] < linestep + act_steps)
              tmpoverlappingidx[3] = linestep + act_steps;
          }
        }
      }
      if (tmpoverlappingidx.size()) newoverlappingidx.push_back(tmpoverlappingidx);
    }
    overlappingidx = newoverlappingidx;

    std::cout << "RefinedBBOverlapCheck: Refinement Level " << refinestep << " there is " << overlap
              << " overlap!" << std::endl;
    if (!overlap) break;
  }
  return overlap;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::IntersectionBase> GEO::CUT::IntersectionFactory::CreateIntersection(
    DRT::Element::DiscretizationType edge_type, DRT::Element::DiscretizationType side_type) const
{
  const int probdim = DRT::Problem::Instance()->NDim();
  switch (edge_type)
  {
    case DRT::Element::line2:
      return Teuchos::rcp(CreateIntersection<DRT::Element::line2>(side_type, probdim));
    default:
      dserror(
          "Unsupported edgeType! If meaningful, add your edgeType here. \n"
          "Given edgeType = %s",
          DRT::DistypeToString(edge_type).c_str());
      break;
  }
  exit(EXIT_FAILURE);
}


template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
std::pair<bool, bool> GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge, dimside,
    numNodesEdge, numNodesSide>::IsQuad4Distorted()
{
  // we assume problem dimension equal to 3 here
  KERNEL::PointOnSurfaceLoc loc[2];
  KERNEL::PointOnSurfaceLoc loc_strict[2];

  double distance[2];

  for (int tri_id = 0; tri_id < 2; ++tri_id)
  {
    bool signeddistance = true;
    int other_point_index = (tri_id == 0) ? 3 : 1;
    LINALG::Matrix<3, 1> point(xyze_surfaceElement_.A() + other_point_index * 3, true);
    LINALG::Matrix<3, 1> xsi;
    LINALG::Matrix<3, 3> xyze_triElement;
    GetTriangle(xyze_triElement, tri_id);

    bool conv = false;
    switch (GetOptionsPtr()->GeomDistance_Floattype())
    {
      case INPAR::CUT::floattype_double:
      {
        KERNEL::ComputeDistance<3, DRT::Element::tri3, false> cd(xsi);
        conv = cd(xyze_triElement, point, distance[tri_id], signeddistance);
        if (conv) loc[tri_id] = cd.GetSideLocation();
        break;
      }
      case INPAR::CUT::floattype_cln:
      {
        KERNEL::ComputeDistance<3, DRT::Element::tri3, true> cd(xsi);
        conv = cd(xyze_triElement, point, distance[tri_id], signeddistance);
        if (conv) loc[tri_id] = cd.GetSideLocation();
        break;
      }
      default:
      {
        throw std::runtime_error("Unexpected floattype for KERNEL::ComputeDistance!");
      }
    }
    if (conv)
    {
      loc_strict[tri_id] = loc[tri_id];
    }
    else
    {
#if EXTENDED_CUT_DEBUG_OUTPUT
      std::cout << "NOTICE: Cut distance did not converge, when detecting distortion of the quad4"
                << std::endl;
#endif
    }
    if (std::abs(distance[tri_id]) >= 1e-30)
    {
      // it is outside of the plane of the side
      loc_strict[tri_id] = KERNEL::PointOnSurfaceLoc(loc[tri_id].WithinSide(), false);
    }
  }

  bool result = false;
  bool result_strict = false;

  for (int tri_id = 0; tri_id < 2; ++tri_id)
  {
    if (not(loc[tri_id].OnSide()))
    {
      result = true;
#if EXTENDED_CUT_DEBUG_OUTPUT

      int other_point_index = (tri_id == 0) ? 3 : 1;
      Point* out_of_plane = GetSidePtr()->Nodes()[other_point_index]->point();
      LINALG::Matrix<3, 1> coord;
      out_of_plane->Coordinates(coord.A());

      for (std::vector<Node*>::const_iterator it = GetSidePtr()->Nodes().begin();
           it != GetSidePtr()->Nodes().end(); ++it)
      {
        Point* p = (*it)->point();
        LINALG::Matrix<3, 1> distc;
        p->Coordinates(distc.A());
        distc.Update(-1, coord, 1);
      }

      std::cout << "Point of the quad4 with Id = " << out_of_plane->Id()
                << " does not lie on the plane created by other 3 points" << std::endl;
      std::cout << "Distance is between out of plane point and side, created by 3 other points is  "
                << distance[tri_id] << std::endl;
#endif
    }
    if (not(loc_strict[tri_id].OnSide())) result_strict = true;
  }
  return std::make_pair(result, result_strict);
}



template class GEO::CUT::Intersection<2, DRT::Element::line2, DRT::Element::line2>;
template class GEO::CUT::Intersection<3, DRT::Element::line2, DRT::Element::line2>;
// template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::quad4>;
template class GEO::CUT::Intersection<3, DRT::Element::line2, DRT::Element::quad4>;
// template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::quad8>;
template class GEO::CUT::Intersection<3, DRT::Element::line2, DRT::Element::quad8>;
// template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::quad9>;
template class GEO::CUT::Intersection<3, DRT::Element::line2, DRT::Element::quad9>;
// template class GEO::CUT::Intersection<2,DRT::Element::line2,DRT::Element::tri3>;
template class GEO::CUT::Intersection<3, DRT::Element::line2, DRT::Element::tri3>;
