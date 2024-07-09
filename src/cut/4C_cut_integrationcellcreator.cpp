/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells for the Tessellation routine

\level 3


*----------------------------------------------------------------------*/

#include "4C_cut_integrationcellcreator.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_integrationcell.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_pointgraph_simple.hpp"
#include "4C_cut_position.hpp"
#include "4C_cut_side.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_cells(
    Mesh& mesh, Element* element, const plain_volumecell_set& cells)
{
  IntegrationCellCreator creator;

  for (plain_volumecell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    VolumeCell* cell = *i;
    bool found = false;
    switch (element->n_dim())
    {
      case 1:
        found = (creator.create_line2_cell(mesh, cell, cell->facets()) or
                 creator.create_point1_cell(mesh, cell, cell->facets()));
        if (not found)
          FOUR_C_THROW(
              "No 1-D cell could be generated! Seems impossible to happen in 1-D! "
              "-- hiermeier 01/17");
        break;
      case 2:
      {
        found = (creator.create2_d_cell<Core::FE::CellType::tri3>(mesh, cell, cell->facets()) or
                 creator.create2_d_cell<Core::FE::CellType::quad4>(mesh, cell, cell->facets()));
        if (not found)
          FOUR_C_THROW(
              "No 2-D cell could be generated and tessellation is currently "
              "unsupported! Thus the given cut case is not yet supported! -- hiermeier 01/17");
        break;
      }
      case 3:
      {
        found = (creator.create_tet4_cell(mesh, cell, cell->facets()) or
                 creator.create_hex8_cell(mesh, cell, cell->facets()) or
                 creator.create_wedge6_cell(mesh, cell, cell->facets()) or
                 creator.create_special_cases(mesh, cell, cell->facets()));
        break;
      }
      default:
        FOUR_C_THROW("Wrong element dimension! ( element dim = %d )", element->n_dim());
        exit(EXIT_FAILURE);
    }

    // pyramids are not save right now.
    // Pyramid5IntegrationCell::CreateCell( mesh, this, facets_, creator ) );
    // creator.create_pyramid5_cell( mesh, cell, cell->Facets() )

    if (not found)
    {
      /* return false in case that not for all volumecells simple-shaped integration
       * cells could be created */
      return false;
    }
  }
  creator.execute(mesh);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_cell(
    Mesh& mesh, Core::FE::CellType shape, VolumeCell* cell)
{
  IntegrationCellCreator creator;

  bool success;
  switch (shape)
  {
    case Core::FE::CellType::tet4:
      success = creator.create_tet4_cell(mesh, cell, cell->facets());
      break;
    case Core::FE::CellType::hex8:
      success = creator.create_hex8_cell(mesh, cell, cell->facets());
      break;
    case Core::FE::CellType::wedge6:
      success = creator.create_wedge6_cell(mesh, cell, cell->facets());
      break;
    case Core::FE::CellType::pyramid5:
      // success = creator.create_pyramid5_cell( mesh, cell, cell->Facets() );
      success = false;
      break;
    case Core::FE::CellType::line2:
      success = creator.create_line2_cell(mesh, cell, cell->facets());
      break;
    case Core::FE::CellType::tri3:
      success = creator.create2_d_cell<Core::FE::CellType::tri3>(mesh, cell, cell->facets());
      break;
    case Core::FE::CellType::quad4:
      success = creator.create2_d_cell<Core::FE::CellType::quad4>(mesh, cell, cell->facets());
      break;
    default:
      FOUR_C_THROW(
          "unsupported element shape ( shape = %s )", Core::FE::CellTypeToString(shape).c_str());
      exit(EXIT_FAILURE);
  }
  // if the create process was successful, we can finally create the integration cell
  if (success)
  {
    creator.execute(mesh);
  }

  return success;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::IntegrationCellCreator::execute(Mesh& mesh)
{
  for (std::map<VolumeCell*, Volume>::iterator it = cells_.begin(); it != cells_.end(); ++it)
  {
    VolumeCell* vc = it->first;
    Volume& cell = it->second;
    cell.execute(mesh, vc);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_point1_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check the facet number
  if (facets.size() != 1) return false;

  // check the actual facet type
  Facet* f = *facets.begin();
  if (not f->equals(Core::FE::CellType::point1)) return false;

  // add the side for the boundary integration cell creation
  if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::point1, f->corner_points());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::IntegrationCellCreator::add_side(
    Core::Geo::Cut::BoundaryCellPosition bcell_position, VolumeCell* vc, Facet* facet,
    Core::FE::CellType shape, const std::vector<Point*>& side)
{
  switch (bcell_position)
  {
    case Core::Geo::Cut::bcells_on_cut_side:
    {
      if (facet->on_boundary_cell_side()) add_side(vc, facet, shape, side);
      break;
    }
    case Core::Geo::Cut::bcells_on_all_sides:
    {
      if (vc->parent_element()->is_cut()) add_side(vc, facet, shape, side);
      break;
    }
    case Core::Geo::Cut::bcells_none:
      /* do nothing */
      break;
    default:
    {
      FOUR_C_THROW("Unknown boundary creation position type! ( enum = %d )", bcell_position);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_line2_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check the facet number
  if (facets.size() != 2) return false;

  std::vector<Point*> line_corner_points;
  line_corner_points.reserve(2);

  const enum Core::Geo::Cut::BoundaryCellPosition bcell_pos =
      mesh.create_options().gen_boundary_cell_position();

  for (plain_facet_set::const_iterator cit = facets.begin(); cit != facets.end(); ++cit)
  {
    // check the actual facet type
    Facet* f = *cit;
    if (not f->equals(Core::FE::CellType::point1)) return false;

    // add the side for the boundary integration cell creation
    add_side(bcell_pos, cell, f, Core::FE::CellType::point1, f->corner_points());

    // collect the facet points
    line_corner_points.push_back(*f->corner_points().begin());
  }

  // check the two points for uniqueness
  if (line_corner_points.size() == 2)
  {
    if (line_corner_points[0] == line_corner_points[1])
      FOUR_C_THROW("The line is not well defined! ( same corner points )");
  }
  else
    FOUR_C_THROW(
        "There are more than two line corner points. You shouldn't "
        "reach this point!");

  // add the actual cell
  add(cell, Core::FE::CellType::line2, line_corner_points);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <Core::FE::CellType celltype, Core::FE::CellType facetype, unsigned numfaces>
bool Core::Geo::Cut::IntegrationCellCreator::create2_d_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check the facet number
  if (facets.size() != numfaces) return false;

  const enum Core::Geo::Cut::BoundaryCellPosition bcell_pos =
      mesh.create_options().gen_boundary_cell_position();

  Impl::SimplePointGraph2D pg = Impl::SimplePointGraph2D();
  pg.find_line_facet_cycles(facets, cell->parent_element());

  for (Impl::PointGraph::facet_iterator it = pg.fbegin(); it != pg.fend(); ++it)
  {
    Cycle& line_cycle = *it;

    if (line_cycle().size() != 2) FOUR_C_THROW("The line cycle has the wrong length!");

    // find corresponding facet
    Facet* f = FindFacet(facets, line_cycle());
    if (not f) FOUR_C_THROW("Could not find the corresponding line facet!");

    add_side(bcell_pos, cell, f, facetype, line_cycle());
  }

  if (pg.num_surfaces() != 1) FOUR_C_THROW("There shouldn't be more than one surface cycle!");

  const Cycle& vol_cycle = (*pg.sbegin());
  add(cell, celltype, vol_cycle());

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_tet4_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  // check if the volumecell has 4 facets and each facet is tri3
  if (facets.size() != 4) return false;

  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    if (not f->equals(Core::FE::CellType::tri3))
    {
      return false;
    }
  }

  // So we have a tet4 here. Create it.
  plain_facet_set::const_iterator cit_facet = facets.begin();
  Facet* bot = *cit_facet;

  const std::vector<Point*>& bot_points = bot->corner_points();
  Point* top_point = nullptr;

  for (++cit_facet; cit_facet != facets.end(); ++cit_facet)
  {
    Facet* f = *cit_facet;

    const std::vector<Point*>& side_points = f->corner_points();
    for (std::vector<Point*>::const_iterator cit_sp = side_points.begin();
         cit_sp != side_points.end(); ++cit_sp)
    {
      Point* p = *cit_sp;
      if (std::find(bot_points.begin(), bot_points.end(), p) == bot_points.end())
      {
        if (top_point == nullptr)
        {
          top_point = p;
        }
        else if (top_point != p)
        {
          FOUR_C_THROW(
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

  Core::LinAlg::Matrix<3, 3> bot_xyze;
  Core::LinAlg::Matrix<3, 1> top_xyz;

  bot->corner_coordinates(bot_xyze.data());
  top_point->coordinates(top_xyz.data());

  Teuchos::RCP<Core::Geo::Cut::Position> bot_distance =
      Core::Geo::Cut::Position::create(bot_xyze, top_xyz, Core::FE::CellType::tri3);

  bot_distance->compute(true);
  const bool invalid_pos = (bot_distance->status() < Position::position_distance_valid);

  if (invalid_pos or bot_distance->distance() >= 0)
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
        side[j] = points[Core::FE::eleNodeNumbering_tet10_surfaces[i][j]];
      }
      Facet* f = FindFacet(facets, side);
      if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::tri3, side);
    }

    /* We create no TET4 cell, if the position calculation failed or the cell
     * is planar. */
    if (not invalid_pos and bot_distance->distance() > 0)
      add(cell, Core::FE::CellType::tet4, points);

    return true;
  }
  else if (bot_distance->distance() < 0)
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
        side[j] = points[Core::FE::eleNodeNumbering_tet10_surfaces[i][j]];
      }
      Facet* f = FindFacet(facets, side);
      if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::tri3, side);
    }

    add(cell, Core::FE::CellType::tet4, points);
    return true;
  }
  else
  {
    FOUR_C_THROW("This cannot happen!");
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_hex8_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  if (facets.size() == 6)
  {
    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (not f->equals(Core::FE::CellType::quad4))
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
    Facet* top = nullptr;
    for (++i; i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (not f->touches(bot))
      {
        top = f;
        break;
      }
    }
    if (top == nullptr)
    {
      FOUR_C_THROW("illegal hex8 cell");
    }

    const std::vector<Point*>& bot_points = bot->corner_points();
    const std::vector<Point*>& top_points = top->corner_points();

    std::vector<Point*> points(8, static_cast<Point*>(nullptr));

    std::copy(bot_points.begin(), bot_points.end(), points.data());

    std::vector<Point*>::iterator end = points.begin() + 4;

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f != bot and f != top)
      {
        const std::vector<Point*>& side_points = f->corner_points();
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
              FOUR_C_THROW("illegal hex8 cell");
            }

            std::vector<Point*>::iterator pointpos2 = std::find(points.begin(), end, *i);
            if (pointpos2 == end)
            {
              i = side_points.end();
              pointpos2 = std::find(points.begin(), end, side_points[3]);
              if (pointpos2 == end)
              {
                FOUR_C_THROW("illegal hex8 cell");
              }
              std::swap(pointpos1, pointpos2);
            }

            unsigned pos = i - side_points.begin() - 1;
            Point* top_point2 = side_points[(pos + 2) % 4];
            Point* top_point1 = side_points[(pos + 3) % 4];

            if (std::find(top_points.begin(), top_points.end(), top_point1) == top_points.end() or
                std::find(top_points.begin(), top_points.end(), top_point2) == top_points.end())
            {
              FOUR_C_THROW("illegal hex8 cell");
            }

            pos = (pointpos1 - points.begin()) + 4;
            if (points[pos] == nullptr)
            {
              points[pos] = top_point1;
            }
            else if (points[pos] != top_point1)
            {
              FOUR_C_THROW("illegal hex8 cell");
            }

            pos = (pointpos2 - points.begin()) + 4;
            if (points[pos] == nullptr)
            {
              points[pos] = top_point2;
            }
            else if (points[pos] != top_point2)
            {
              FOUR_C_THROW("illegal hex8 cell");
            }

            break;
          }
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    Core::LinAlg::Matrix<3, 4> bot_xyze;
    Core::LinAlg::Matrix<3, 4> top_xyze;

    bot->corner_coordinates(bot_xyze.data());
    top->corner_coordinates(top_xyze.data());

    int distance_counter = 0;
    for (int i = 0; i < 4; ++i)
    {
      Core::LinAlg::Matrix<3, 1> top_xyz(&top_xyze(0, i), true);
      Teuchos::RCP<Core::Geo::Cut::Position> bot_distance =
          Core::Geo::Cut::Position::create(bot_xyze, top_xyz, Core::FE::CellType::quad4);

      bot_distance->compute(true);

      // If the distance calculation is not possible, we return false.
      const bool invalid_pos = (bot_distance->status() < Position::position_distance_valid);
      if (invalid_pos) return false;

      if (bot_distance->distance() > 0)
      {
        distance_counter += 1;
      }
      if (bot_distance->distance() < 0)
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
          side[j] = rpoints[Core::FE::eleNodeNumbering_hex27_surfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::quad4, side);
      }

      add(cell, Core::FE::CellType::hex8, rpoints);
      return true;
    }
    else if (distance_counter == 4)
    {
      for (int i = 0; i < 6; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[Core::FE::eleNodeNumbering_hex27_surfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::quad4, side);
      }

      add(cell, Core::FE::CellType::hex8, points);
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
bool Core::Geo::Cut::IntegrationCellCreator::create_wedge6_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  if (facets.size() == 5)
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f->equals(Core::FE::CellType::tri3))
      {
        tris.push_back(f);
      }
      else if (f->equals(Core::FE::CellType::quad4))
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

    const std::vector<Point*>& bot_points = bot->corner_points();
    const std::vector<Point*>& top_points = top->corner_points();

    std::vector<Point*> points(6, static_cast<Point*>(nullptr));

    std::copy(bot_points.begin(), bot_points.end(), points.data());

    std::vector<Point*>::iterator end = points.begin() + 3;

    for (std::vector<Facet*>::const_iterator i = quads.begin(); i != quads.end(); ++i)
    {
      Facet* f = *i;
      const std::vector<Point*>& side_points = f->corner_points();
      for (std::vector<Point*>::const_iterator i = side_points.begin(); i != side_points.end(); ++i)
      {
        Point* p = *i;
        std::vector<Point*>::iterator pointpos1 = std::find(points.begin(), end, p);
        if (pointpos1 != end)
        {
          ++i;
          if (i == side_points.end())
          {
            FOUR_C_THROW("illegal wedge6 cell");
          }

          std::vector<Point*>::iterator pointpos2 = std::find(points.begin(), end, *i);
          if (pointpos2 == end)
          {
            i = side_points.end();
            pointpos2 = std::find(points.begin(), end, side_points[3]);
            if (pointpos2 == end)
            {
              FOUR_C_THROW("illegal wedge6 cell");
            }
            std::swap(pointpos1, pointpos2);
          }

          unsigned pos = i - side_points.begin() - 1;
          Point* top_point2 = side_points[(pos + 2) % 4];
          Point* top_point1 = side_points[(pos + 3) % 4];

          if (std::find(top_points.begin(), top_points.end(), top_point1) == top_points.end() or
              std::find(top_points.begin(), top_points.end(), top_point2) == top_points.end())
          {
            FOUR_C_THROW("illegal wedge6 cell");
          }

          pos = (pointpos1 - points.begin()) + 3;
          if (points[pos] == nullptr)
          {
            points[pos] = top_point1;
          }
          else if (points[pos] != top_point1)
          {
            FOUR_C_THROW("illegal wedge6 cell");
          }

          pos = (pointpos2 - points.begin()) + 3;
          if (points[pos] == nullptr)
          {
            points[pos] = top_point2;
          }
          else if (points[pos] != top_point2)
          {
            FOUR_C_THROW("illegal wedge6 cell");
          }

          break;
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    Core::LinAlg::Matrix<3, 3> bot_xyze;
    Core::LinAlg::Matrix<3, 3> top_xyze;

    bot->corner_coordinates(bot_xyze.data());
    top->corner_coordinates(top_xyze.data());

    int distance_counter = 0;
    for (int i = 0; i < 3; ++i)
    {
      Core::LinAlg::Matrix<3, 1> top_xyz(&top_xyze(0, i), true);
      Teuchos::RCP<Core::Geo::Cut::Position> bot_distance =
          Core::Geo::Cut::Position::create(bot_xyze, top_xyz, Core::FE::CellType::tri3);

      bot_distance->compute(true);

      // If the distance calculation is not possible, we return false.
      const bool invalid_pos = (bot_distance->status() < Position::position_distance_valid);
      if (invalid_pos) return false;

      if (bot_distance->distance() > 0)
      {
        distance_counter += 1;
      }
      if (bot_distance->distance() < 0)
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
          side[j] = rpoints[Core::FE::eleNodeNumbering_wedge18_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::tri3, side);
      }
      for (int i = 0; i < 3; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = rpoints[Core::FE::eleNodeNumbering_wedge18_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::quad4, side);
      }

      add(cell, Core::FE::CellType::wedge6, rpoints);
      return true;
    }
    else if (distance_counter == 3)
    {
      for (int i = 0; i < 2; ++i)
      {
        std::vector<Point*> side(3);
        for (int j = 0; j < 3; ++j)
        {
          side[j] = points[Core::FE::eleNodeNumbering_wedge18_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::tri3, side);
      }
      for (int i = 0; i < 3; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[Core::FE::eleNodeNumbering_wedge18_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::quad4, side);
      }

      add(cell, Core::FE::CellType::wedge6, points);
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
bool Core::Geo::Cut::IntegrationCellCreator::create_pyramid5_cell(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  if (facets.size() == 5)
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* f = *i;
      if (f->equals(Core::FE::CellType::tri3))
      {
        tris.push_back(f);
      }
      else if (f->equals(Core::FE::CellType::quad4))
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

    const std::vector<Point*>& bot_points = bot->corner_points();
    Point* top_point = nullptr;

    for (std::vector<Facet*>::iterator i = tris.begin(); i != tris.end(); ++i)
    {
      Facet* f = *i;

      const std::vector<Point*>& side_points = f->corner_points();
      for (std::vector<Point*>::const_iterator i = side_points.begin(); i != side_points.end(); ++i)
      {
        Point* p = *i;
        if (std::find(bot_points.begin(), bot_points.end(), p) == bot_points.end())
        {
          if (top_point == nullptr)
          {
            top_point = p;
          }
          else if (top_point != p)
          {
            // Corner point confusion. This is actually not a pyramid5.
            //
            // FOUR_C_THROW( "illegal pyramid5 cell" );
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

    Core::LinAlg::Matrix<3, 4> bot_xyze;
    Core::LinAlg::Matrix<3, 1> top_xyze;

    bot->corner_coordinates(bot_xyze.data());
    top_point->coordinates(top_xyze.data());

    Teuchos::RCP<Position> bot_distance =
        Core::Geo::Cut::Position::create(bot_xyze, top_xyze, Core::FE::CellType::quad4);
    bot_distance->compute(true);

    // check the status of the position computation
    const bool invalid_pos = (bot_distance->status() < Position::position_distance_valid);

    if (invalid_pos or bot_distance->distance() >= 0)
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
          side[j] = points[Core::FE::eleNodeNumbering_pyramid5_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::tri3, side);
      }
      for (int i = 0; i < 1; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[Core::FE::eleNodeNumbering_pyramid5_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::quad4, side);
      }

      /* We create no PYRAMID5 cell, if the position calculation failed or the cell
       * is planar. */
      if ((not invalid_pos) and (bot_distance->distance() > 0))
        add(cell, Core::FE::CellType::pyramid5, points);
      return true;
    }
    else if (bot_distance->distance() < 0)
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
          side[j] = points[Core::FE::eleNodeNumbering_pyramid5_trisurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::tri3, side);
      }
      for (int i = 0; i < 1; ++i)
      {
        std::vector<Point*> side(4);
        for (int j = 0; j < 4; ++j)
        {
          side[j] = points[Core::FE::eleNodeNumbering_pyramid5_quadsurfaces[i][j]];
        }
        Facet* f = FindFacet(facets, side);
        if (f->on_boundary_cell_side()) add_side(cell, f, Core::FE::CellType::quad4, side);
      }

      // cell->NewPyramid5Cell( mesh, points );
      add(cell, Core::FE::CellType::pyramid5, points);
      return true;
    }
    else
    {
      FOUR_C_THROW("This is not possible!");
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCellCreator::create_special_cases(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets)
{
  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    Facet* f = *i;
    if (f->has_holes())
    {
      return false;
    }
  }

  Element* parent = cell->parent_element();
  const std::vector<Side*>& sides = parent->sides();

  switch (parent->shape())
  {
    case Core::FE::CellType::hex8:
    {
      // find how many element sides are touched by this volume cell and how
      // often those sides are touched.
      std::vector<int> touched(6, 0);
      for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
      {
        Facet* f = *i;
        if (not f->on_boundary_cell_side())
        {
          Side* s = f->parent_side();
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
      const std::vector<Side*>& sides = parent->sides();
      for (std::vector<Side*>::const_iterator i = sides.begin(); i != sides.end(); ++i)
      {
        Side* s = *i;
        cut.push_back(s->is_cut());
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
          return hex8_horizontal_cut(mesh, parent, cell, facets, axis, r);
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
          return hex8_horizontal_cut(mesh, parent, cell, facets, axis, r);
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
          return hex8_horizontal_cut(mesh, parent, cell, facets, axis, r);
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
bool Core::Geo::Cut::IntegrationCellCreator::hex8_horizontal_cut(Mesh& mesh, Element* element,
    VolumeCell* cell, const plain_facet_set& facets, int axis, double r)
{
  //  Point::PointPosition position = cell->Position();

  PointSet cut_points;
  cell->get_all_points(mesh, cut_points);

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
    if (f->on_boundary_cell_side())
    {
      inner_facets.push_back(f);
      const std::vector<Point*>& points = f->corner_points();
      std::copy(points.begin(), points.end(), std::inserter(inner, inner.begin()));
    }
  }

  std::vector<Point*> inner_points;
  std::vector<Point*> projected_points;

  inner_points.reserve(inner.size());
  // std::copy( inner.begin(), inner.end(), std::back_inserter( inner_points ) );
  inner_points.assign(inner.begin(), inner.end());

  // project along given axis to r

  Core::LinAlg::Matrix<3, 1> rst;

  std::vector<Core::LinAlg::Matrix<3, 1>> local_points;
  local_points.reserve(inner_points.size());
  projected_points.reserve(inner_points.size());

  for (std::vector<Point*>::iterator i = inner_points.begin(); i != inner_points.end(); ++i)
  {
    Point* p = *i;
    Core::LinAlg::Matrix<3, 1> xyz;
    p->coordinates(xyz.data());
    element->local_coordinates(xyz, rst);

    local_points.push_back(rst);

    // projection
    rst(axis) = r;

    // create new points
    element->global_coordinates(rst, xyz);
    projected_points.push_back(mesh.new_point(xyz.data(), nullptr, nullptr, 0.0));
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
    if (f->equals(Core::FE::CellType::tri3))
    {
      std::vector<Point*> points;
      points.reserve(6);

      const std::vector<Point*>& corner_points = f->corner_points();

      // find facet orientation in local coordinates
      double drs = 0;
      Core::LinAlg::Matrix<3, 3> xyze;
      Core::LinAlg::Matrix<3, 1> normal;
      Core::LinAlg::Matrix<2, 3> deriv;
      Core::LinAlg::Matrix<2, 2> metrictensor;

      double* x = xyze.data();

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
          FOUR_C_THROW("inner point missing");
        }

        points.push_back(projected_points[pos - inner_points.begin()]);

        const Core::LinAlg::Matrix<3, 1>& rst = local_points[pos - inner_points.begin()];
        x = std::copy(rst.data(), rst.data() + 3, x);
      }

      if (r < 0)
      {
        sidepos = 1;
        // std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
        points.insert(points.end(), corner_points.begin(), corner_points.end());
      }

      Core::FE::shape_function_2D_deriv1(deriv, 0., 0., Core::FE::CellType::tri3);
      Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::tri3>(
          xyze, deriv, metrictensor, drs, &normal);

      if (normal(axis) < 0)
      {
        std::swap(points[1], points[2]);
        std::swap(points[4], points[5]);
      }

      std::vector<Point*> side(3);
      for (int j = 0; j < 3; ++j)
      {
        side[j] = points[Core::FE::eleNodeNumbering_wedge18_trisurfaces[sidepos][j]];
      }
      // Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      add_side(cell, f, Core::FE::CellType::tri3, side);

      // cell->NewWedge6Cell( mesh, points );
      add(cell, Core::FE::CellType::wedge6, points);
    }
    else if (f->equals(Core::FE::CellType::quad4))
    {
      std::vector<Point*> points;
      points.reserve(8);

      const std::vector<Point*>& corner_points = f->corner_points();

      // find facet orientation in local coordinates
      double drs = 0;
      Core::LinAlg::Matrix<3, 4> xyze;
      Core::LinAlg::Matrix<3, 1> normal;
      Core::LinAlg::Matrix<2, 4> deriv;
      Core::LinAlg::Matrix<2, 2> metrictensor;

      double* x = xyze.data();

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
          FOUR_C_THROW("inner point missing");
        }

        points.push_back(projected_points[pos - inner_points.begin()]);

        const Core::LinAlg::Matrix<3, 1>& rst = local_points[pos - inner_points.begin()];
        x = std::copy(rst.data(), rst.data() + 3, x);
      }

      if (r < 0)
      {
        sidepos = 5;
        // std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
        points.insert(points.end(), corner_points.begin(), corner_points.end());
      }

      Core::FE::shape_function_2D_deriv1(deriv, 0., 0., Core::FE::CellType::quad4);
      Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::quad4>(
          xyze, deriv, metrictensor, drs, &normal);

      if (normal(axis) < 0)
      {
        std::swap(points[1], points[3]);
        std::swap(points[5], points[7]);
      }

      std::vector<Point*> side(4);
      for (int j = 0; j < 4; ++j)
      {
        side[j] = points[Core::FE::eleNodeNumbering_hex27_surfaces[sidepos][j]];
      }
      // Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      add_side(cell, f, Core::FE::CellType::quad4, side);

      // cell->NewHex8Cell( mesh, points );
      add(cell, Core::FE::CellType::hex8, points);
    }
    else
      return false;
  }
  return true;
}

template bool Core::Geo::Cut::IntegrationCellCreator::create2_d_cell<Core::FE::CellType::tri3>(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);
template bool Core::Geo::Cut::IntegrationCellCreator::create2_d_cell<Core::FE::CellType::quad4>(
    Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

FOUR_C_NAMESPACE_CLOSE
