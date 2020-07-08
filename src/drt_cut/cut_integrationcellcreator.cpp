/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells for the Tessellation routine

\level 3


*----------------------------------------------------------------------*/

#include "cut_integrationcellcreator.H"
#include "cut_position.H"
#include "cut_integrationcell.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_side.H"
#include "cut_options.H"

#include "cut_pointgraph_simple.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateCells(
    Mesh& mesh, Element* element, const plain_volumecell_set& cells)
{
  IntegrationCellCreator creator;

  for (plain_volumecell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    VolumeCell* cell = *i;
    bool found = false;
    switch (element->Dim())
    {
      case 1:
        found = (creator.CreateLine2Cell(mesh, cell, cell->Facets()) or
                 creator.CreatePoint1Cell(mesh, cell, cell->Facets()));
        if (not found)
          dserror(
              "No 1-D cell could be generated! Seems impossible to happen in 1-D! "
              "-- hiermeier 01/17");
        break;
      case 2:
      {
        found = (creator.Create2DCell<DRT::Element::tri3>(mesh, cell, cell->Facets()) or
                 creator.Create2DCell<DRT::Element::quad4>(mesh, cell, cell->Facets()));
        if (not found)
          dserror(
              "No 2-D cell could be generated and tessellation is currently "
              "unsupported! Thus the given cut case is not yet supported! -- hiermeier 01/17");
        break;
      }
      case 3:
      {
        found = (creator.CreateTet4Cell(mesh, cell, cell->Facets()) or
                 creator.CreateHex8Cell(mesh, cell, cell->Facets()) or
                 creator.CreateWedge6Cell(mesh, cell, cell->Facets()) or
                 creator.CreateSpecialCases(mesh, cell, cell->Facets()));
        break;
      }
      default:
        dserror("Wrong element dimension! ( element dim = %d )", element->Dim());
        exit(EXIT_FAILURE);
    }

    // pyramids are not save right now.
    // Pyramid5IntegrationCell::CreateCell( mesh, this, facets_, creator ) );
    // creator.CreatePyramid5Cell( mesh, cell, cell->Facets() )

    if (not found)
    {
      /* return false in case that not for all volumecells simple-shaped integration
       * cells could be created */
      return false;
    }
  }
  creator.Execute(mesh);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateCell(
    Mesh& mesh, DRT::Element::DiscretizationType shape, VolumeCell* cell)
{
  IntegrationCellCreator creator;

  bool success;
  switch (shape)
  {
    case DRT::Element::tet4:
      success = creator.CreateTet4Cell(mesh, cell, cell->Facets());
      break;
    case DRT::Element::hex8:
      success = creator.CreateHex8Cell(mesh, cell, cell->Facets());
      break;
    case DRT::Element::wedge6:
      success = creator.CreateWedge6Cell(mesh, cell, cell->Facets());
      break;
    case DRT::Element::pyramid5:
      // success = creator.CreatePyramid5Cell( mesh, cell, cell->Facets() );
      success = false;
      break;
    case DRT::Element::line2:
      success = creator.CreateLine2Cell(mesh, cell, cell->Facets());
      break;
    case DRT::Element::tri3:
      success = creator.Create2DCell<DRT::Element::tri3>(mesh, cell, cell->Facets());
      break;
    case DRT::Element::quad4:
      success = creator.Create2DCell<DRT::Element::quad4>(mesh, cell, cell->Facets());
      break;
    default:
      dserror("unsupported element shape ( shape = %s )", DRT::DistypeToString(shape).c_str());
      exit(EXIT_FAILURE);
  }
  // if the create process was successful, we can finally create the integration cell
  if (success)
  {
    creator.Execute(mesh);
  }

  return success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IntegrationCellCreator::Execute(Mesh& mesh)
{
  for (std::map<VolumeCell*, volume>::iterator it = cells_.begin(); it != cells_.end(); ++it)
  {
    VolumeCell* vc = it->first;
    volume& cell = it->second;
    cell.Execute(mesh, vc);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreatePoint1Cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check the facet number
  if (facets.size() != 1) return false;

  // check the actual facet type
  Facet* f = *facets.begin();
  if (not f->Equals(DRT::Element::point1)) return false;

  // add the side for the boundary integration cell creation
  if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::point1, f->CornerPoints());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::IntegrationCellCreator::AddSide(INPAR::CUT::BoundaryCellPosition bcell_position,
    VolumeCell* vc, Facet* facet, DRT::Element::DiscretizationType shape,
    const std::vector<Point*>& side)
{
  switch (bcell_position)
  {
    case INPAR::CUT::bcells_on_cut_side:
    {
      if (facet->OnBoundaryCellSide()) AddSide(vc, facet, shape, side);
      break;
    }
    case INPAR::CUT::bcells_on_all_sides:
    {
      if (vc->ParentElement()->IsCut()) AddSide(vc, facet, shape, side);
      break;
    }
    case INPAR::CUT::bcells_none:
      /* do nothing */
      break;
    default:
    {
      dserror("Unknown boundary creation position type! ( enum = %d )", bcell_position);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateLine2Cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check the facet number
  if (facets.size() != 2) return false;

  std::vector<Point*> line_corner_points;
  line_corner_points.reserve(2);

  const enum INPAR::CUT::BoundaryCellPosition bcell_pos =
      mesh.CreateOptions().GenBoundaryCellPosition();

  for (plain_facet_set::const_iterator cit = facets.begin(); cit != facets.end(); ++cit)
  {
    // check the actual facet type
    Facet* f = *cit;
    if (not f->Equals(DRT::Element::point1)) return false;

    // add the side for the boundary integration cell creation
    AddSide(bcell_pos, cell, f, DRT::Element::point1, f->CornerPoints());

    // collect the facet points
    line_corner_points.push_back(*f->CornerPoints().begin());
  }

  // check the two points for uniqueness
  if (line_corner_points.size() == 2)
  {
    if (line_corner_points[0] == line_corner_points[1])
      dserror("The line is not well defined! ( same corner points )");
  }
  else
    dserror(
        "There are more than two line corner points. You shouldn't "
        "reach this point!");

  // add the actual cell
  Add(cell, DRT::Element::line2, line_corner_points);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
#if __cplusplus < 201103L
template <DRT::Element::DiscretizationType celltype>
bool GEO::CUT::IntegrationCellCreator::Create2DCell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  const DRT::Element::DiscretizationType facetype =
      DRT::UTILS::DisTypeToFaceShapeType<celltype>::shape;
  const unsigned numfaces = DRT::UTILS::DisTypeToNumFaces<celltype>::numFaces;
#else
template <DRT::Element::DiscretizationType celltype, DRT::Element::DiscretizationType facetype,
    unsigned numfaces>
bool GEO::CUT::IntegrationCellCreator::Create2DCell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
#endif
  // check the facet number
  if (facets.size() != numfaces) return false;

  const enum INPAR::CUT::BoundaryCellPosition bcell_pos =
      mesh.CreateOptions().GenBoundaryCellPosition();

  IMPL::SimplePointGraph_2D pg = IMPL::SimplePointGraph_2D();
  pg.FindLineFacetCycles(facets, cell->ParentElement());

  for (IMPL::PointGraph::facet_iterator it = pg.fbegin(); it != pg.fend(); ++it)
  {
    Cycle& line_cycle = *it;

    if (line_cycle().size() != 2) dserror("The line cycle has the wrong length!");

    // find corresponding facet
    Facet* f = FindFacet(facets, line_cycle());
    if (not f) dserror("Could not find the corresponding line facet!");

    AddSide(bcell_pos, cell, f, facetype, line_cycle());
  }

  if (pg.NumSurfaces() != 1) dserror("There shouldn't be more than one surface cycle!");

  const Cycle& vol_cycle = (*pg.sbegin());
  Add(cell, celltype, vol_cycle());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateTet4Cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check if the volumecell has 4 facets and each facet is tri3
  if (facets.size() != 4) return false;

  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    if (not f->Equals(DRT::Element::tri3))
    {
      return false;
    }
  }

  // So we have a tet4 here. Create it.
  plain_facet_set::const_iterator cit_facet = facets.begin();
  Facet* bot = *cit_facet;

  const std::vector<Point*>& bot_points = bot->CornerPoints();
  Point* top_point = NULL;

  for (++cit_facet; cit_facet != facets.end(); ++cit_facet)
  {
    Facet* f = *cit_facet;

    const std::vector<Point*>& side_points = f->CornerPoints();
    for (std::vector<Point*>::const_iterator cit_sp = side_points.begin();
         cit_sp != side_points.end(); ++cit_sp)
    {
      Point* p = *cit_sp;
      if (std::find(bot_points.begin(), bot_points.end(), p) == bot_points.end())
      {
        if (top_point == NULL)
        {
          top_point = p;
        }
        else if (top_point != p)
        {
          dserror(
              "Illegal tet4 cell. We found a side point of TET4 which "
              "belongs neither to the bottom facet nor is it the top node!");
        }
      }
    }
  }

  //     for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  //     {
  //       Facet * f = *i;
  //       f->NewTri3Cell( mesh );
  //     }

  LINALG::Matrix<3, 3> bot_xyze;
  LINALG::Matrix<3, 1> top_xyz;

  bot->CornerCoordinates(bot_xyze.A());
  top_point->Coordinates(top_xyz.A());

  Teuchos::RCP<GEO::CUT::Position> bot_distance =
      GEO::CUT::Position::Create(bot_xyze, top_xyz, DRT::Element::tri3);

  bot_distance->Compute(true);
  const bool invalid_pos = (bot_distance->Status() < Position::position_distance_valid);

  if (invalid_pos or bot_distance->Distance() >= 0)
  {
    std::vector<Point*> points;
    points.reserve(4);
    // std::copy( bot_points.begin(), bot_points.end(), std::back_inserter( points ) );
    points.assign(bot_points.begin(), bot_points.end());
    points.push_back(top_point);

    for (int i = 0; i < 4; ++i)
    {
      std::vector<Point*> side(3);
      for (int j = 0; j < 3; ++j)
      {
        side[j] = points[DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j]];
      }
      Facet* f = FindFacet(facets, side);
      if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::tri3, side);
    }

    /* We create no TET4 cell, if the position calculation failed or the cell
     * is planar. */
    if (not invalid_pos and bot_distance->Distance() > 0) Add(cell, DRT::Element::tet4, points);

    return true;
  }
  else if (bot_distance->Distance() < 0)
  {
    std::vector<Point*> points;
    points.reserve(4);
    // std::copy( bot_points.rbegin(), bot_points.rend(), std::back_inserter( points ) );
    points.assign(bot_points.rbegin(), bot_points.rend());
    points.push_back(top_point);

    for (int i = 0; i < 4; ++i)
    {
      std::vector<Point*> side(3);
      for (int j = 0; j < 3; ++j)
      {
        side[j] = points[DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j]];
      }
      Facet* f = FindFacet(facets, side);
      if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::tri3, side);
    }

    Add(cell, DRT::Element::tet4, points);
    return true;
  }
  else
  {
    dserror("This cannot happen!");
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateHex8Cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  if (facets.size() == 6)
  {
    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (not f->Equals(DRT::Element::quad4))
      {
        return false;
      }
    }

    // Must not be concave?!

    // So we have a hex8 here. Create it.
    //
    // Find two sides that do not share a common point, find positive
    // directions for both, and find the edges between those sides.

    plain_facet_set::const_iterator i = facets.begin();
    Facet* bot = *i;
    Facet* top = NULL;
    for (++i; i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (not f->Touches(bot))
      {
        top = f;
        break;
      }
    }
    if (top == NULL)
    {
      run_time_error("illegal hex8 cell");
    }

    const std::vector<Point*>& bot_points = bot->CornerPoints();
    const std::vector<Point*>& top_points = top->CornerPoints();

    std::vector<Point*> points(8, static_cast<Point*>(NULL));

    std::copy(bot_points.begin(), bot_points.end(), &points[0]);

    std::vector<Point*>::iterator end = points.begin() + 4;

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f != bot and f != top)
      {
        const std::vector<Point*>& side_points = f->CornerPoints();
        for (std::vector<Point*>::const_iterator i = side_points.begin(); i != side_points.end();
             ++i)
        {
          Point* p = *i;
          std::vector<Point*>::iterator pointpos1 = std::find(points.begin(), end, p);
          if (pointpos1 != end)
          {
            ++i;
            if (i == side_points.end())
            {
              run_time_error("illegal hex8 cell");
            }

            std::vector<Point*>::iterator pointpos2 = std::find(points.begin(), end, *i);
            if (pointpos2 == end)
            {
              i = side_points.end();
              pointpos2 = std::find(points.begin(), end, side_points[3]);
              if (pointpos2 == end)
              {
                run_time_error("illegal hex8 cell");
              }
              std::swap(pointpos1, pointpos2);
            }

            unsigned pos = i - side_points.begin() - 1;
            Point* top_point2 = side_points[(pos + 2) % 4];
            Point* top_point1 = side_points[(pos + 3) % 4];

            if (std::find(top_points.begin(), top_points.end(), top_point1) == top_points.end() or
                std::find(top_points.begin(), top_points.end(), top_point2) == top_points.end())
            {
              run_time_error("illegal hex8 cell");
            }

            pos = (pointpos1 - points.begin()) + 4;
            if (points[pos] == NULL)
            {
              points[pos] = top_point1;
            }
            else if (points[pos] != top_point1)
            {
              run_time_error("illegal hex8 cell");
            }

            pos = (pointpos2 - points.begin()) + 4;
            if (points[pos] == NULL)
            {
              points[pos] = top_point2;
            }
            else if (points[pos] != top_point2)
            {
              run_time_error("illegal hex8 cell");
            }

            break;
          }
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    LINALG::Matrix<3, 4> bot_xyze;
    LINALG::Matrix<3, 4> top_xyze;

    bot->CornerCoordinates(bot_xyze.A());
    top->CornerCoordinates(top_xyze.A());

    int distance_counter = 0;
    for (int i = 0; i < 4; ++i)
    {
      LINALG::Matrix<3, 1> top_xyz(&top_xyze(0, i), true);
      Teuchos::RCP<GEO::CUT::Position> bot_distance =
          GEO::CUT::Position::Create(bot_xyze, top_xyz, DRT::Element::quad4);

      bot_distance->Compute(true);

      // If the distance calculation is not possible, we return false.
      const bool invalid_pos = (bot_distance->Status() < Position::position_distance_valid);
      if (invalid_pos) return false;

      if (bot_distance->Distance() > 0)
      {
        distance_counter += 1;
      }
      if (bot_distance->Distance() < 0)
      {
        distance_counter -= 1;
      }
    }

    if (distance_counter == -4)
    {
      std::vector<Point*> rpoints(8);
      rpoints[0] = points[3];
      rpoints[1] = points[2];
      rpoints[2] = points[1];
      rpoints[3] = points[0];
      rpoints[4] = points[7];
      rpoints[5] = points[6];
      rpoints[6] = points[5];
      rpoints[7] = points[4];

      for (int i = 0; i < 6; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = rpoints[DRT::UTILS::eleNodeNumbering_hex27_surfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::quad4, side);
      }

      Add(cell, DRT::Element::hex8, rpoints);
      return true;
    }
    else if (distance_counter == 4)
    {
      for (int i = 0; i < 6; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_hex27_surfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::quad4, side);
      }

      Add(cell, DRT::Element::hex8, points);
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateWedge6Cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  if (facets.size() == 5)
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f->Equals(DRT::Element::tri3))
      {
        tris.push_back(f);
      }
      else if (f->Equals(DRT::Element::quad4))
      {
        quads.push_back(f);
      }
      else
      {
        return false;
      }
    }

    if (tris.size() != 2 or quads.size() != 3)
    {
      return false;
    }

    Facet* bot = tris[0];
    Facet* top = tris[1];

    const std::vector<Point*>& bot_points = bot->CornerPoints();
    const std::vector<Point*>& top_points = top->CornerPoints();

    std::vector<Point*> points(6, static_cast<Point*>(NULL));

    std::copy(bot_points.begin(), bot_points.end(), &points[0]);

    std::vector<Point*>::iterator end = points.begin() + 3;

    for (std::vector<Facet*>::const_iterator i = quads.begin(); i != quads.end(); ++i)
    {
      Facet* f = *i;
      const std::vector<Point*>& side_points = f->CornerPoints();
      for (std::vector<Point*>::const_iterator i = side_points.begin(); i != side_points.end(); ++i)
      {
        Point* p = *i;
        std::vector<Point*>::iterator pointpos1 = std::find(points.begin(), end, p);
        if (pointpos1 != end)
        {
          ++i;
          if (i == side_points.end())
          {
            run_time_error("illegal wedge6 cell");
          }

          std::vector<Point*>::iterator pointpos2 = std::find(points.begin(), end, *i);
          if (pointpos2 == end)
          {
            i = side_points.end();
            pointpos2 = std::find(points.begin(), end, side_points[3]);
            if (pointpos2 == end)
            {
              run_time_error("illegal wedge6 cell");
            }
            std::swap(pointpos1, pointpos2);
          }

          unsigned pos = i - side_points.begin() - 1;
          Point* top_point2 = side_points[(pos + 2) % 4];
          Point* top_point1 = side_points[(pos + 3) % 4];

          if (std::find(top_points.begin(), top_points.end(), top_point1) == top_points.end() or
              std::find(top_points.begin(), top_points.end(), top_point2) == top_points.end())
          {
            run_time_error("illegal wedge6 cell");
          }

          pos = (pointpos1 - points.begin()) + 3;
          if (points[pos] == NULL)
          {
            points[pos] = top_point1;
          }
          else if (points[pos] != top_point1)
          {
            run_time_error("illegal wedge6 cell");
          }

          pos = (pointpos2 - points.begin()) + 3;
          if (points[pos] == NULL)
          {
            points[pos] = top_point2;
          }
          else if (points[pos] != top_point2)
          {
            run_time_error("illegal wedge6 cell");
          }

          break;
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    LINALG::Matrix<3, 3> bot_xyze;
    LINALG::Matrix<3, 3> top_xyze;

    bot->CornerCoordinates(bot_xyze.A());
    top->CornerCoordinates(top_xyze.A());

    int distance_counter = 0;
    for (int i = 0; i < 3; ++i)
    {
      LINALG::Matrix<3, 1> top_xyz(&top_xyze(0, i), true);
      Teuchos::RCP<GEO::CUT::Position> bot_distance =
          GEO::CUT::Position::Create(bot_xyze, top_xyz, DRT::Element::tri3);

      bot_distance->Compute(true);

      // If the distance calculation is not possible, we return false.
      const bool invalid_pos = (bot_distance->Status() < Position::position_distance_valid);
      if (invalid_pos) return false;

      if (bot_distance->Distance() > 0)
      {
        distance_counter += 1;
      }
      if (bot_distance->Distance() < 0)
      {
        distance_counter -= 1;
      }
    }

    //     if ( distance_counter==-3 or
    //          distance_counter==3 )
    //     {
    //       for ( std::vector<Facet*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
    //       {
    //         Facet * f = *i;
    //         f->NewTri3Cell( mesh );
    //       }
    //       for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
    //       {
    //         Facet * f = *i;
    //         f->NewQuad4Cell( mesh );
    //       }
    //     }

    if (distance_counter == -3)
    {
      std::vector<Point*> rpoints(6);
      rpoints[0] = points[2];
      rpoints[1] = points[1];
      rpoints[2] = points[0];
      rpoints[3] = points[5];
      rpoints[4] = points[4];
      rpoints[5] = points[3];

      for (int i = 0; i < 2; ++i)
      {
        std::vector<Point*> side(3);
        for (int j = 0; j < 3; ++j)
        {
          side[j] = rpoints[DRT::UTILS::eleNodeNumbering_wedge18_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::tri3, side);
      }
      for (int i = 0; i < 3; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = rpoints[DRT::UTILS::eleNodeNumbering_wedge18_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::quad4, side);
      }

      Add(cell, DRT::Element::wedge6, rpoints);
      return true;
    }
    else if (distance_counter == 3)
    {
      for (int i = 0; i < 2; ++i)
      {
        std::vector<Point*> side(3);
        for (int j = 0; j < 3; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_wedge18_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::tri3, side);
      }
      for (int i = 0; i < 3; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_wedge18_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::quad4, side);
      }

      Add(cell, DRT::Element::wedge6, points);
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreatePyramid5Cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  if (facets.size() == 5)
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f->Equals(DRT::Element::tri3))
      {
        tris.push_back(f);
      }
      else if (f->Equals(DRT::Element::quad4))
      {
        quads.push_back(f);
      }
      else
      {
        return false;
      }
    }

    if (tris.size() != 4 or quads.size() != 1)
    {
      return false;
    }

    Facet* bot = quads[0];

    const std::vector<Point*>& bot_points = bot->CornerPoints();
    Point* top_point = NULL;

    for (std::vector<Facet*>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
      Facet* f = *i;

      const std::vector<Point*>& side_points = f->CornerPoints();
      for (std::vector<Point*>::const_iterator i = side_points.begin(); i != side_points.end(); ++i)
      {
        Point* p = *i;
        if (std::find(bot_points.begin(), bot_points.end(), p) == bot_points.end())
        {
          if (top_point == NULL)
          {
            top_point = p;
          }
          else if (top_point != p)
          {
            // Corner point confusion. This is actually not a pyramid5.
            //
            // throw std::runtime_error( "illegal pyramid5 cell" );
            return false;
          }
        }
      }
    }

    //     for ( std::vector<Facet*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
    //     {
    //       Facet * f = *i;
    //       f->NewTri3Cell( mesh );
    //     }
    //     for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
    //     {
    //       Facet * f = *i;
    //       f->NewQuad4Cell( mesh );
    //     }

    LINALG::Matrix<3, 4> bot_xyze;
    LINALG::Matrix<3, 1> top_xyze;

    bot->CornerCoordinates(bot_xyze.A());
    top_point->Coordinates(top_xyze.A());

    Teuchos::RCP<Position> bot_distance =
        GEO::CUT::Position::Create(bot_xyze, top_xyze, DRT::Element::quad4);
    bot_distance->Compute(true);

    // check the status of the position computation
    const bool invalid_pos = (bot_distance->Status() < Position::position_distance_valid);

    if (invalid_pos or bot_distance->Distance() >= 0)
    {
      std::vector<Point*> points;
      points.reserve(5);
      // std::copy( bot_points.begin(), bot_points.end(), std::back_inserter( points ) );
      points.assign(bot_points.begin(), bot_points.end());
      points.push_back(top_point);

      for (int i = 0; i < 4; ++i)
      {
        std::vector<Point*> side(3);
        for (int j = 0; j < 3; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::tri3, side);
      }
      for (int i = 0; i < 1; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::quad4, side);
      }

      /* We create no PYRAMID5 cell, if the position calculation failed or the cell
       * is planar. */
      if ((not invalid_pos) and (bot_distance->Distance() > 0))
        Add(cell, DRT::Element::pyramid5, points);
      return true;
    }
    else if (bot_distance->Distance() < 0)
    {
      std::vector<Point*> points;
      points.reserve(5);
      // std::copy( bot_points.rbegin(), bot_points.rend(), std::back_inserter( points ) );
      points.assign(bot_points.rbegin(), bot_points.rend());
      points.push_back(top_point);

      for (int i = 0; i < 4; ++i)
      {
        std::vector<Point*> side(3);
        for (int j = 0; j < 3; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::tri3, side);
      }
      for (int i = 0; i < 1; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->OnBoundaryCellSide()) AddSide(cell, f, DRT::Element::quad4, side);
      }

      // cell->NewPyramid5Cell( mesh, points );
      Add(cell, DRT::Element::pyramid5, points);
      return true;
    }
    else
    {
      dserror("This is not possible!");
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::CreateSpecialCases(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    if (f->HasHoles())
    {
      return false;
    }
  }

  Element* parent = cell->ParentElement();
  const std::vector<Side*>& sides = parent->Sides();

  switch (parent->Shape())
  {
    case DRT::Element::hex8:
    {
      // find how many element sides are touched by this volume cell and how
      // often those sides are touched.
      std::vector<int> touched(6, 0);
      for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
      {
        Facet* f = *i;
        if (not f->OnBoundaryCellSide())
        {
          Side* s = f->ParentSide();
          std::vector<Side*>::const_iterator pos = std::find(sides.begin(), sides.end(), s);
          if (pos != sides.end())
          {
            touched[pos - sides.begin()] += 1;
          }
        }
      }

      int touched_size = 0;
      for (std::vector<int>::iterator i = touched.begin(); i != touched.end(); ++i)
      {
        if (*i > 0)
        {
          touched_size += 1;
        }
        if (*i > 1)
        {
          return false;
        }
      }

      int uncutcount = 0;
      std::vector<int> cut;
      cut.reserve(6);
      const std::vector<Side*>& sides = parent->Sides();
      for (std::vector<Side*>::const_iterator i = sides.begin(); i != sides.end(); ++i)
      {
        Side* s = *i;
        cut.push_back(s->IsCut());
        if (not cut.back()) uncutcount += 1;
      }

      if (uncutcount == 2)
      {
        double r = 0.;
        int axis = -1;

        // We use shards face numbering. Be careful.
        if (not cut[4] and not cut[5])
        {
          if (touched_size != 5)
          {
            return false;
          }

          axis = 2;
          if (touched[5] > 0)
          {
            r = 1;
          }
          else
          {
            r = -1;
          }
          return Hex8HorizontalCut(mesh, parent, cell, facets, axis, r);
        }
        else if (not cut[0] and not cut[2])
        {
          if (touched_size != 5)
          {
            return false;
          }

          axis = 1;
          if (touched[2] > 0)
          {
            r = 1;
          }
          else
          {
            r = -1;
          }
          return Hex8HorizontalCut(mesh, parent, cell, facets, axis, r);
        }
        else if (not cut[1] and not cut[3])
        {
          if (touched_size != 5)
          {
            return false;
          }

          axis = 0;
          if (touched[1] > 0)
          {
            r = 1;
          }
          else
          {
            r = -1;
          }
          return Hex8HorizontalCut(mesh, parent, cell, facets, axis, r);
        }
      }
      return false;
    }
    default:
      break;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::IntegrationCellCreator::Hex8HorizontalCut(Mesh& mesh, Element* element,
    VolumeCell* cell, const plain_facet_set& facets, int axis, double r)
{
  //  Point::PointPosition position = cell->Position();

  PointSet cut_points;
  cell->GetAllPoints(mesh, cut_points);

  std::vector<Point*> points;
  points.reserve(cut_points.size());
  points.assign(cut_points.begin(), cut_points.end());

  // std::cout << "hex8 projection along axis " << axis << " to " << r << "\n";

  // find all inner points that need projecting

  PointSet inner;

  std::vector<Facet*> inner_facets;

  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    if (f->OnBoundaryCellSide())
    {
      inner_facets.push_back(f);
      const std::vector<Point*>& points = f->CornerPoints();
      std::copy(points.begin(), points.end(), std::inserter(inner, inner.begin()));
    }
  }

  std::vector<Point*> inner_points;
  std::vector<Point*> projected_points;

  inner_points.reserve(inner.size());
  // std::copy( inner.begin(), inner.end(), std::back_inserter( inner_points ) );
  inner_points.assign(inner.begin(), inner.end());

  // project along given axis to r

  LINALG::Matrix<3, 1> rst;

  std::vector<LINALG::Matrix<3, 1>> local_points;
  local_points.reserve(inner_points.size());
  projected_points.reserve(inner_points.size());

  for (std::vector<Point*>::iterator i = inner_points.begin(); i != inner_points.end(); ++i)
  {
    Point* p = *i;
    LINALG::Matrix<3, 1> xyz;
    p->Coordinates(xyz.A());
    element->LocalCoordinates(xyz, rst);

    local_points.push_back(rst);

    // projection
    rst(axis) = r;

    // create new points
    element->GlobalCoordinates(rst, xyz);
    projected_points.push_back(mesh.NewPoint(xyz.A(), NULL, NULL, 0.0));
    // change Benedikt: do not set the position for additionally created points
    /* REMARK:
     * the propagation of the inside/outside position to facets and volume cells
     * can destroy on-cut-surface points! At the moment we do not need the
     * inside/outside information of newly created points for tetrahedralization */
    //    projected_points.back()->Position( position );
  }

  // create integration cells

  for (std::vector<Facet*>::iterator i = inner_facets.begin(); i != inner_facets.end(); ++i)
  {
    Facet* f = *i;
    if (f->Equals(DRT::Element::tri3))
    {
      std::vector<Point*> points;
      points.reserve(6);

      const std::vector<Point*>& corner_points = f->CornerPoints();

      // find facet orientation in local coordinates
      double drs = 0;
      LINALG::Matrix<3, 3> xyze;
      LINALG::Matrix<3, 1> normal;
      LINALG::Matrix<2, 3> deriv;
      LINALG::Matrix<2, 2> metrictensor;

      double* x = xyze.A();

      int sidepos = -1;
      if (r > 0)
      {
        sidepos = 0;
        // std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
        points.assign(corner_points.begin(), corner_points.end());
      }

      for (std::vector<Point*>::const_iterator i = corner_points.begin(); i != corner_points.end();
           ++i)
      {
        Point* p = *i;

        std::vector<Point*>::iterator pos = std::find(inner_points.begin(), inner_points.end(), p);
        if (pos == inner_points.end())
        {
          run_time_error("inner point missing");
        }

        points.push_back(projected_points[pos - inner_points.begin()]);

        const LINALG::Matrix<3, 1>& rst = local_points[pos - inner_points.begin()];
        x = std::copy(rst.A(), rst.A() + 3, x);
      }

      if (r < 0)
      {
        sidepos = 1;
        // std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
        points.insert(points.end(), corner_points.begin(), corner_points.end());
      }

      DRT::UTILS::shape_function_2D_deriv1(deriv, 0., 0., DRT::Element::tri3);
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::tri3>(
          xyze, deriv, metrictensor, drs, &normal);

      if (normal(axis) < 0)
      {
        std::swap(points[1], points[2]);
        std::swap(points[4], points[5]);
      }

      std::vector<Point*> side(3);
      for (int j = 0; j < 3; ++j)
      {
        side[j] = points[DRT::UTILS::eleNodeNumbering_wedge18_trisurfaces[sidepos][j]];
      }
      // Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      AddSide(cell, f, DRT::Element::tri3, side);

      // cell->NewWedge6Cell( mesh, points );
      Add(cell, DRT::Element::wedge6, points);
    }
    else if (f->Equals(DRT::Element::quad4))
    {
      std::vector<Point*> points;
      points.reserve(8);

      const std::vector<Point*>& corner_points = f->CornerPoints();

      // find facet orientation in local coordinates
      double drs = 0;
      LINALG::Matrix<3, 4> xyze;
      LINALG::Matrix<3, 1> normal;
      LINALG::Matrix<2, 4> deriv;
      LINALG::Matrix<2, 2> metrictensor;

      double* x = xyze.A();

      int sidepos = -1;
      if (r > 0)
      {
        sidepos = 0;
        // std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
        points.assign(corner_points.begin(), corner_points.end());
      }

      for (std::vector<Point*>::const_iterator i = corner_points.begin(); i != corner_points.end();
           ++i)
      {
        Point* p = *i;

        std::vector<Point*>::iterator pos = std::find(inner_points.begin(), inner_points.end(), p);
        if (pos == inner_points.end())
        {
          run_time_error("inner point missing");
        }

        points.push_back(projected_points[pos - inner_points.begin()]);

        const LINALG::Matrix<3, 1>& rst = local_points[pos - inner_points.begin()];
        x = std::copy(rst.A(), rst.A() + 3, x);
      }

      if (r < 0)
      {
        sidepos = 5;
        // std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
        points.insert(points.end(), corner_points.begin(), corner_points.end());
      }

      DRT::UTILS::shape_function_2D_deriv1(deriv, 0., 0., DRT::Element::quad4);
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad4>(
          xyze, deriv, metrictensor, drs, &normal);

      if (normal(axis) < 0)
      {
        std::swap(points[1], points[3]);
        std::swap(points[5], points[7]);
      }

      std::vector<Point*> side(4);
      for (int j = 0; j < 4; ++j)
      {
        side[j] = points[DRT::UTILS::eleNodeNumbering_hex27_surfaces[sidepos][j]];
      }
      // Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      AddSide(cell, f, DRT::Element::quad4, side);

      // cell->NewHex8Cell( mesh, points );
      Add(cell, DRT::Element::hex8, points);
    }
    else
    {
#if 0
      PointSet points;
      const std::vector<Point*> & corner_points = f->CornerPoints();

      f->AllPoints( points );

      for ( std::vector<Point*>::const_iterator i=corner_points.begin();
            i!=corner_points.end();
            ++i )
      {
        Point * p = *i;

        std::vector<Point*>::iterator pos = std::find( inner_points.begin(), inner_points.end(), p );
        if ( pos==inner_points.end() )
        {
          throw std::runtime_error( "inner point missing" );
        }

        points.insert( projected_points[pos - inner_points.begin()] );
      }

      std::vector<Point*> cell_points;
      cell_points.reserve( points.size() );
      //std::copy( points.begin(), points.end(), std::back_inserter( cell_points ) );
      cell_points.assign( points.begin(), points.end() );

      // sort points that go into qhull to obtain the same result independent of
      // pointer values (compiler flags, code structure, memory usage, ...)
      std::sort( cell_points.begin(), cell_points.end(), PointPidLess() );

      plain_facet_set cell_facets;
      cell_facets.insert( f );
      cell->CreateTet4IntegrationCells( mesh, position, cell_points, cell_facets, true );
#endif
      return false;
    }
  }
  return true;
}

template bool GEO::CUT::IntegrationCellCreator::Create2DCell<DRT::Element::tri3>(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);
template bool GEO::CUT::IntegrationCellCreator::Create2DCell<DRT::Element::quad4>(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);
