/*----------------------------------------------------------------------------*/
/*! \file

\brief Create a simple point graph for 1-D and 2-D elements ( embedded in a
       higher dimensional space )


\date Nov 12, 2016

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "cut_pointgraph_simple.H"

#include "cut_side.H"
#include "cut_line.H"
#include "cut_facet.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::IMPL::SimplePointGraph_1D::SimplePointGraph_1D(Mesh &mesh, Element *element, Side *side,
    PointGraph::Location location, PointGraph::Strategy strategy)
    : PointGraph::PointGraph(1) /* explicit call of the protected base
                                   class constructor */
{
  if (element->Dim() != 1) dserror("This class is only meaningful for 1-D elements!");

  Cycle cycle;
  FillGraph(element, side, cycle, strategy);

  GetGraph().FindCycles(element, side, cycle, location, strategy);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_1D::BuildCycle(
    const std::vector<Point *> &edge_points, Cycle &cycle) const
{
  /* since a closed cycle is not possible for 1-D elements, we have to consider
   * also the starting edge point */
  for (std::vector<Point *>::const_iterator cit = edge_points.begin(); cit != edge_points.end();
       ++cit)
  {
    Point *p = *cit;
    cycle.push_back(p);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_1D::AddCutLinesToGraph(
    Element *element, Side *side, Strategy strategy, Cycle &cycle)
{
  // this is completely unneccesary for the element_side
  if (not side->IsCutSide()) return;

  AddCutPointsToCycle(element, side, cycle);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_1D::AddCutPointsToCycle(
    Element *element, Side *side, Cycle &cycle)
{
  std::vector<Side *> ele_sides = element->Sides();
  if (ele_sides.size() != 1) dserror("There should be only one side in the 1-D case!");

  Side &ele_side = **ele_sides.begin();
  PointSet cuts;

  // get the cut point(s) (probably only one)
  ele_side.GetCutPoints(element, *side, cuts);

  if (cuts.size() > 1)
    dserror("There are a multiple cut points for the 2 line cut, check this case.");

  // add the found cut points and actual in the case of a cut_side, no other points
  // were added till now, to the cycle.
  for (PointSet::const_iterator cit = cuts.begin(); cit != cuts.end(); ++cit) cycle.push_back(*cit);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_1D::Graph::FindCycles(
    Element *element, Side *side, Cycle &cycle, Location location, Strategy strategy)
{
  // each point defines a own main cycle in 1D
  Edge *edge = element->Sides()[0]->Edges()[0];
  std::map<double, Cycle> sorted_cycles;

  for (std::vector<Point *>::const_iterator cit = cycle().begin(); cit != cycle().end(); ++cit)
  {
    Point *p = *cit;
    double t = p->t(edge);

    Cycle pseudo_1d_cycle(std::vector<Point *>(1, p));
    sorted_cycles[t] = pseudo_1d_cycle;
  }

  main_cycles_.reserve(sorted_cycles.size());
  for (std::map<double, Cycle>::const_iterator cit = sorted_cycles.begin();
       cit != sorted_cycles.end(); ++cit)
    main_cycles_.push_back(cit->second);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::IMPL::SimplePointGraph_2D::SimplePointGraph_2D(Mesh &mesh, Element *element, Side *side,
    PointGraph::Location location, PointGraph::Strategy strategy)
    : PointGraph::PointGraph(mesh, element, side, location, strategy),
      graph_2d_(Teuchos::rcp_dynamic_cast<SimplePointGraph_2D::Graph>(GraphPtr(), true))
{
  if (element->Dim() != 2) dserror("This class is only meaningful for 2-D elements!");
  /* intentionally left blank, the work is done in the base class */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IMPL::SimplePointGraph_2D::Graph::HasSinglePoints(Location location)
{
  // it is allowed for the 2-D case
  if (location == cut_side) return false;

  return PointGraph::Graph::HasSinglePoints(location);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_2D::Graph::FindCycles(
    Element *element, Side *side, Cycle &cycle, Location location, Strategy strategy)
{
  if (location == element_side)
  {
    GEO::CUT::IMPL::PointGraph::Graph::FindCycles(element, side, cycle, location, strategy);

    if (correct_rotation_direction_) CorrectRotationDirection(element->Sides()[0], main_cycles_);

    SplitMainCyclesIntoLineCycles();
  }
  else
  {
    if (all_points_.size() > 2)
    {
      for (auto i = all_points_.begin(); i != all_points_.end(); ++i) i->second->Print();
      dserror(
          "A cut_side point cycle with more than 2 points is currently not "
          "supported in the 2-D case, since it cannot happen for level-set "
          "cut cases ( always straight cuts ). Check this scenario if you run "
          "into this and extend the functionality.");
    }

    unsigned lcount = 0;
    std::vector<Point *> line(2, NULL);
    for (std::map<int, Point *>::const_iterator cit = all_points_.begin(); cit != all_points_.end();
         ++cit)
    {
      line[lcount++] = cit->second;
    }
    main_cycles_ = std::vector<Cycle>(1, Cycle(line));
  }

#if 0
  std::cout << "New main-Cycles\n";
  for ( std::vector<Cycle>::const_iterator cit = main_cycles_.begin();
        cit != main_cycles_.end(); ++cit )
    cit->Print();
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_2D::Graph::SplitMainCyclesIntoLineCycles()
{
  std::vector<Cycle> line_main_cycles;
  for (std::vector<Cycle>::const_iterator ic = main_cycles_.begin(); ic != main_cycles_.end(); ++ic)
  {
    Cycle closed_cycle = *ic;
    const unsigned cycle_length = closed_cycle.size();

    std::vector<Point *> line(2, NULL);
    // split the connected cycle into line segments
    for (unsigned i = 0; i < cycle_length; ++i)
    {
      line[0] = closed_cycle()[i % cycle_length];
      line[1] = closed_cycle()[(i + 1) % cycle_length];
      line_main_cycles.push_back(Cycle(line));
    }
  }
  // store the previously generated surface main cycles
  surface_main_cycles_ = main_cycles_;

  // replace the main cycles by the splitted line main cycles
  main_cycles_ = line_main_cycles;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::IMPL::SimplePointGraph_2D::SimplePointGraph_2D()
    : GEO::CUT::IMPL::PointGraph(2),
      graph_2d_(Teuchos::rcp_dynamic_cast<SimplePointGraph_2D::Graph>(GraphPtr(), true))
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_2D::FindLineFacetCycles(
    const plain_facet_set &line_facets, Element *parent_element)
{
  Cycle cycle;
  FillGraphAndCycleWithLineFacets(line_facets, cycle);

  graph_2d_->SetCorrectRotationDirection(true);

  FindCycles(parent_element, cycle);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_2D::FillGraphAndCycleWithLineFacets(
    const plain_facet_set &line_facets, Cycle &cycle)
{
  plain_point_set point_set;
  for (plain_facet_set::const_iterator cit = line_facets.begin(); cit != line_facets.end(); ++cit)
  {
    Facet *f = *cit;

    if (not f->Equals(DRT::Element::line2)) dserror("This function works only for line facets!");

    const std::vector<Point *> &line_points = f->Points();
    GetGraph().AddEdge(line_points[0], line_points[1]);

    point_set.insert(line_points.begin(), line_points.end());
  }

  cycle = Cycle(std::vector<Point *>(point_set.begin(), point_set.end()));

#if 0
  GetGraph().Print();
  cycle.Print();
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_2D::FindCycles(Element *element, Cycle &cycle)
{
  Side *side = element->Sides()[0];
  GetGraph().FindCycles(element, side, cycle, PointGraph::element_side, PointGraph::all_lines);
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IMPL::SimplePointGraph_2D::CorrectRotationDirection(
    const Side *side, std::vector<Cycle> &cycles)
{
  LINALG::Matrix<2, 1> rs = DRT::UTILS::getLocalCenterPosition<2>(side->Shape());
  LINALG::Matrix<3, 1> normal_side(false);

  LINALG::SerialDenseMatrix xyze_side(3, side->NumNodes());
  side->Coordinates(xyze_side);

  // get the normal direction of the underlying side
  EvalNormalVectors<3>(side->Shape(), xyze_side, rs, normal_side);

  std::vector<Point *> invert_cycle;
  for (std::vector<Cycle>::iterator it = cycles.begin(); it != cycles.end(); ++it)
  {
    Cycle &cycle = *it;

    const unsigned cycle_size = cycle().size();
    if (cycle_size < 3) dserror("The cycle needs at least three points!");

    std::vector<Point *>::const_iterator it_begin = cycle().begin();
    std::vector<Point *>::const_iterator it_end = cycle().end() - 1;

    const LINALG::Matrix<3, 1> x1((*it_begin)->X(), true);
    const LINALG::Matrix<3, 1> x2((*(++it_begin))->X(), true);
    const LINALG::Matrix<3, 1> x3((*it_end)->X(), true);

    LINALG::Matrix<3, 1> x12;
    LINALG::Matrix<3, 1> x13;

    x12.Update(1.0, x2, -1.0, x1);
    x13.Update(1.0, x3, -1.0, x1);

    LINALG::Matrix<3, 1> normal_cycle;
    normal_cycle.CrossProduct(x12, x13);

    // not really necessary, but we do it anyway ...
    double nrm2 = normal_cycle.Norm2();
    normal_cycle.Scale(1.0 / nrm2);

    // check orientation and skip the remaining part, if the orientation is okay
    if (normal_cycle.Dot(normal_side) > 0.0) continue;

    invert_cycle.reserve(cycle_size);

    it_begin = cycle().end() - 1;
    it_end = cycle().begin() - 1;
    for (std::vector<Point *>::const_iterator cit = it_begin; cit != it_end; --cit)
    {
      invert_cycle.push_back(*cit);
    }

    // replace the old the cycle by the inverted one
    cycle = Cycle(invert_cycle);
    invert_cycle.clear();
  }
}
