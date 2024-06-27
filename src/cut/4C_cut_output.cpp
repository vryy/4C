/*---------------------------------------------------------------------*/
/*! \file

\brief Handles file writing of all cut related stuff (gmsh)

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_output.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_cycle.hpp"
#include "4C_cut_edge.hpp"
#include "4C_cut_element.hpp"
#include "4C_cut_facet.hpp"
#include "4C_cut_kernel.hpp"
#include "4C_cut_line.hpp"
#include "4C_cut_point.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_volumecell.hpp"

#include <iosfwd>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshFacetsOnly(
    const plain_facet_set& facets, Element* ele, const std::string& file_affix)
{
  // write details of volume cells
  std::stringstream str;
  str << ".facets" << (not file_affix.empty() ? "_" + file_affix : "") << "_CUTFAIL.pos";
  std::string filename(Core::Geo::Cut::Output::GenerateGmshOutputFilename(str.str()));

  std::cout << "\nFacets of element " << ele->Id() << " are written to " << filename << "\n";

  std::ofstream file(filename.c_str());
  int count = 0;
  for (plain_facet_set::const_iterator it = facets.begin(); it != facets.end(); ++it)
  {
    std::ostringstream sectionname;
    sectionname << "Facet_" << count++;
    Output::GmshNewSection(file, sectionname.str(), false);
    Output::GmshFacetDump(file, *it, "sides", true, false, ele);
    Output::GmshEndSection(file);
  }
  file.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshEdgesOnly(const plain_edge_set& edges)
{
  std::stringstream str;
  str << ".edges"
      << "_CUTFAIL.pos";
  std::string filename(Core::Geo::Cut::Output::GenerateGmshOutputFilename(str.str()));
  std::ofstream file(filename.c_str());
  Core::IO::cout << "\nEdges are written to " << filename << "\n";

  int count = 0;
  for (plain_edge_set::const_iterator it = edges.begin(); it != edges.end(); ++it)
  {
    std::ostringstream sectionname;
    sectionname << "Edge_" << count++;
    Output::GmshNewSection(file, sectionname.str(), false);
    Output::GmshEdgeDump(file, *it);
    Output::GmshEndSection(file);
  }
  file.close();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshVolumeCellsOnly(const plain_volumecell_set& vcells)
{
  // write details of volume cells
  Core::IO::cout << "\nVolumeCells are written to [...].volumecells_CUTFAIL.pos\n";

  std::stringstream str;
  str << ".volumecells"
      << "_CUTFAIL.pos";
  std::string filename(Core::Geo::Cut::Output::GenerateGmshOutputFilename(str.str()));
  std::ofstream file(filename.c_str());
  int count = 0;
  for (plain_volumecell_set::const_iterator its = vcells.begin(); its != vcells.end(); ++its)
  {
    std::ostringstream sectionname;
    sectionname << "VolumeCell_" << count++;
    Output::GmshNewSection(file, sectionname.str(), false);
    Output::GmshVolumecellDump(file, *its, "sides", true, false, (*its)->parent_element());
    Output::GmshEndSection(file);
  }
  file.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
char Core::Geo::Cut::Output::GmshElementType(Core::FE::CellType shape)
{
  switch (shape)
  {
    case Core::FE::CellType::point1:
    {
      return 'P';
    }
    case Core::FE::CellType::line2:
    {
      return 'L';
    }
    case Core::FE::CellType::tri3:
    {
      return 'T';
    }
    case Core::FE::CellType::quad4:
    {
      return 'Q';
    }
    case Core::FE::CellType::hex8:
    {
      return 'H';
    }
    case Core::FE::CellType::tet4:
    {
      return 'S';
    }
    case Core::FE::CellType::wedge6:
    {
      return 'I';
    }
    case Core::FE::CellType::pyramid5:
    {
      return 'P';
    }
    case Core::FE::CellType::dis_none:
    {
      // arbitrary cells are not yet supported
      return ' ';
    }
    default:
      FOUR_C_THROW(
          "Unsupported cell shape! ( shape = %s )", Core::FE::CellTypeToString(shape).c_str());
      exit(EXIT_FAILURE);
  }
  // impossible to reach this point
  exit(EXIT_FAILURE);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given element                                           sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshElementDump(std::ofstream& file, Element* ele, bool to_local)
{
  const std::vector<Node*>& nodes = ele->Nodes();
  char elementtype = GmshElementType(ele->Shape());
  GmshElementDump(file, nodes, elementtype, to_local, ele);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given element                                           sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshElementDump(std::ofstream& file,
    const std::vector<Core::Geo::Cut::Node*>& nodes, char elementtype, bool to_local, Element* ele)
{
  file << "S" << elementtype << "(";
  for (std::vector<Core::Geo::Cut::Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    if (i != nodes.begin()) file << ",";
    GmshWriteCoords(file, *i, to_local, ele);
  }
  file << "){";
  for (std::vector<Core::Geo::Cut::Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Core::Geo::Cut::Node* n = *i;
    Core::Geo::Cut::Point* p = n->point();
    if (i != nodes.begin()) file << ",";
#if CUT_DEVELOP
    file << p->Id();
#else
    file << p->Position();
#endif
  }
  file << "};\n";
}



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshCellDump(std::ofstream& file, Core::FE::CellType shape,
    const Core::LinAlg::SerialDenseMatrix& xyze, const Point::PointPosition* position,
    const int* value)
{
  char elementtype = GmshElementType(shape);

  file.precision(16);
  file << "S" << elementtype << "(";
  for (unsigned i = 0; i < static_cast<unsigned>(xyze.numCols()); ++i)
  {
    if (i > 0) file << ", ";
    for (unsigned j = 0; j < static_cast<unsigned>(xyze.numRows()); ++j)
    {
      if (j > 0) file << ",";
      file << xyze(j, i);
    }
  }
  file << "){";
  for (unsigned i = 0; i < static_cast<unsigned>(xyze.numCols()); ++i)
  {
    if (i > 0) file << ",";
    if (value)
      file << (*value);
    else if (position)
      file << (*position);
    // dummy value
    else
      file << 0.0;
  }
  file << "};\n";
}

void Core::Geo::Cut::Output::GmshSideDump(
    std::ofstream& file, const Side* s, const std::string& sname, bool to_local, Element* ele)
{
  GmshNewSection(file, sname, false);
  GmshSideDump(file, s, false, ele);
  file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given side                                       sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshSideDump(
    std::ofstream& file, const Side* s, bool to_local, Element* ele)
{
  const std::vector<Node*>& nodes = s->Nodes();
  char elementtype;
  switch (nodes.size())
  {
    case 0:
      return;  // I'm a Levelset Side - do nothing!
    case 2:
      elementtype = 'L';
      break;
    case 3:
      elementtype = 'T';
      break;
    case 4:
      elementtype = 'Q';
      break;
    default:
      std::stringstream str;
      str << "unknown element type in GmshSideDump for " << nodes.size() << " nodes!";
      FOUR_C_THROW(str.str());
      exit(EXIT_FAILURE);
  }
  GmshElementDump(file, nodes, elementtype, to_local, ele);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given side                                                ager 04/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshTriSideDump(
    std::ofstream& file, std::vector<Point*> points, bool to_local, Element* ele)
{
  char elementtype;
  switch (points.size())
  {
    case 3:
      elementtype = 'T';
      break;
    case 4:
      elementtype = 'Q';
      break;
    default:
    {
      std::stringstream str;
      str << "unknown element type in GmshTriSideDump for " << points.size() << " points!";
      FOUR_C_THROW(str.str());
      exit(EXIT_FAILURE);
    }
  }

  file << "S" << elementtype << "(";
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    if (i != 0) file << ",";
    GmshWriteCoords(file, points[i], to_local, ele);
  }
  file << "){";
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    Core::Geo::Cut::Point* p = points[i];
    if (i != 0) file << ",";
    file << p->Position();
  }
  file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given facet                                                ager 08/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshFacetDump(std::ofstream& file, Facet* facet,
    const std::string& visualizationtype, bool print_all, bool to_local, Element* ele)
{
  if (to_local)
  {
    if (ele == nullptr)
      FOUR_C_THROW("GmshWriteCoords: Didn't get a parent element for the Coordinate!");
  }

  if (visualizationtype == "sides")
  {
    if (facet->Points().size() == 1)
    {
      GmshFacetDump(file, facet, "points", print_all, to_local, ele);
      return;
    }

    if (facet->Points().size() == 2)
    {
      GmshFacetDump(file, facet, "lines", print_all, to_local, ele);
      return;
    }

    if (facet->IsTriangulated())
    {
      for (unsigned j = 0; j < facet->Triangulation().size(); ++j)
      {
        Core::Geo::Cut::Output::GmshTriSideDump(file, facet->Triangulation()[j], to_local, ele);
      }
    }
    else if (facet->IsFacetSplit())
    {
      for (unsigned j = 0; j < facet->GetSplitCells().size(); ++j)
      {
        Core::Geo::Cut::Output::GmshTriSideDump(file, facet->GetSplitCells()[j], to_local, ele);
      }
    }
    else if (facet->belongs_to_level_set_side() || facet->CornerPoints().size() == 3 ||
             facet->CornerPoints().size() == 4)
    {
      Core::Geo::Cut::Output::GmshTriSideDump(file, facet->CornerPoints(), to_local, ele);
    }
    else if (facet->CornerPoints().size() > 2 &&
             print_all)  // do mitpoint triangulation just for visualization reasons! (not usedfull
                         // if you want to check if a triangulation exists!)
    {
      std::vector<double> xmid(3, 0);
      for (unsigned int i = 0; i < facet->CornerPoints().size(); ++i)
      {
        for (unsigned int dim = 0; dim < 3; ++dim) xmid[dim] += facet->CornerPoints()[i]->X()[dim];
      }
      for (unsigned int dim = 0; dim < 3; ++dim)
        xmid[dim] = xmid[dim] / facet->CornerPoints().size();

      ConcretePoint<3> midpoint = ConcretePoint<3>(-1, xmid.data(), nullptr, nullptr, 0.0);

      std::vector<Point*> tri;
      for (unsigned int i = 0; i < facet->CornerPoints().size(); ++i)
      {
        tri.clear();
        tri.push_back(facet->CornerPoints()[i]);
        tri.push_back(facet->CornerPoints()[(i + 1) % facet->CornerPoints().size()]);
        tri.push_back(&midpoint);
        Core::Geo::Cut::Output::GmshTriSideDump(file, tri, to_local, ele);
      }
    }
  }
  else if (visualizationtype == "lines")
  {
    for (std::size_t pidx = 0; pidx < facet->Points().size(); ++pidx)
    {
      Core::Geo::Cut::Output::GmshLineDump(file, facet->Points()[pidx],
          facet->Points()[(pidx + 1) % facet->Points().size()], to_local, ele);
    }
  }
  else if (visualizationtype == "points")
  {
    for (std::size_t pidx = 0; pidx < facet->Points().size(); ++pidx)
      Core::Geo::Cut::Output::GmshPointDump(
          file, facet->Points()[pidx], facet->SideId(), to_local, ele);
  }
  else
    FOUR_C_THROW("GmshFacetDump: unknown visualizationtype!");
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given volumecell                                       ager 08/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshVolumecellDump(std::ofstream& file, VolumeCell* VC,
    const std::string& visualizationtype, bool print_all, bool to_local, Element* ele)
{
  for (plain_facet_set::const_iterator j = VC->Facets().begin(); j != VC->Facets().end(); j++)
    GmshFacetDump(file, *j, visualizationtype, print_all, to_local, ele);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given cycle                                                ager 08/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshCycleDump(std::ofstream& file, Cycle* cycle,
    const std::string& visualizationtype, bool to_local, Element* ele)
{
  if (visualizationtype == "points")
  {
    for (unsigned i = 0; i != (*cycle)().size(); ++i)
      GmshPointDump(file, (*cycle)()[i], to_local, ele);
  }
  else if (visualizationtype == "lines")
  {
    for (unsigned i = 0; i != (*cycle)().size(); ++i)
      GmshLineDump(file, (*cycle)()[i], (*cycle)()[(i + 1) % (*cycle)().size()], to_local, ele);
  }
  else
    FOUR_C_THROW("GmshFacetDump: unknown visualizationtype!");
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of element along with all its cut sides                sudhakar 03/14
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshCompleteCutElement(
    std::ofstream& file, Element* ele, bool to_local)
{
  // write details of background element
  GmshNewSection(file, "Element");
  GmshElementDump(file, ele, to_local);

  // write details of points
  GmshNewSection(file, "Points", true);
  const std::vector<Point*>& ps = ele->Points();
  for (std::vector<Point*>::const_iterator its = ps.begin(); its != ps.end(); its++)
    GmshPointDump(file, *its, to_local, ele);

  // write details of cut sides
  GmshNewSection(file, "Cut_Facets", true);
  const plain_facet_set& facets = ele->Facets();
  for (plain_facet_set::const_iterator its = facets.begin(); its != facets.end(); its++)
  {
    if ((*its)->ParentSide()->IsCutSide()) GmshFacetDump(file, *its, "sides", true, to_local, ele);
  }

  // write details of cut sides
  GmshNewSection(file, "Ele_Facets", true);
  for (plain_facet_set::const_iterator its = facets.begin(); its != facets.end(); its++)
  {
    if (!(*its)->ParentSide()->IsCutSide()) GmshFacetDump(file, *its, "sides", true, to_local, ele);
  }

  // write details of volumecells
  // GmshNewSection( file, "Volumecells",true);
  const plain_volumecell_set& vcs = ele->VolumeCells();
  for (plain_volumecell_set::const_iterator its = vcs.begin(); its != vcs.end(); its++)
  {
    GmshNewSection(file, "Volumecells", true);
    GmshVolumecellDump(file, *its, "sides", true, to_local, ele);
    GmshEndSection(file);
  }

  // write details of cut sides
  GmshNewSection(file, "Cut sides", false);
  const plain_side_set& cutsides = ele->CutSides();
  for (plain_side_set::const_iterator its = cutsides.begin(); its != cutsides.end(); its++)
    GmshSideDump(file, *its, to_local, ele);
  GmshEndSection(file);

  if (ele->HasLevelSetSide())
  {
    GmshNewSection(file, "LevelSetValues");
    GmshLevelSetValueDump(file, ele, true, to_local);  // true -> dumps LS values at nodes as well.

    GmshNewSection(file, "LevelSetGradient", true);
    GmshLevelSetGradientDump(file, ele, to_local);

    GmshNewSection(file, "LevelSetOrientation", true);
    GmshLevelSetOrientationDump(file, ele, to_local);

    GmshNewSection(file, "LevelSetZeroShape", true);
    GmshLevelSetValueZeroSurfaceDump(file, ele, to_local);
    GmshEndSection(file);
  }
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given line                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshLineDump(
    std::ofstream& file, Core::Geo::Cut::Line* line, bool to_local, Element* ele)
{
  GmshLineDump(file, line->BeginPoint(), line->EndPoint(), to_local, ele);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given line                                           ager 08/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshLineDump(std::ofstream& file, Core::Geo::Cut::Point* p1,
    Core::Geo::Cut::Point* p2, int idx1, int idx2, bool to_local, Element* ele)
{
  file << "SL (";
  GmshWriteCoords(file, p1, to_local, ele);
  file << ",";
  GmshWriteCoords(file, p2, to_local, ele);
  file << "){";
  file << idx1 << ",";
  file << idx2;
  file << "};\n";
}


void Core::Geo::Cut::Output::GmshCutPairDump(
    std::ofstream& file, Side* side, Edge* edge, int id, const std::string& suffix)
{
  std::stringstream side_section_name;
  side_section_name << "Side" << suffix << id;
  Core::Geo::Cut::Output::GmshNewSection(file, side_section_name.str());
  Core::Geo::Cut::Output::GmshSideDump(file, side, false, nullptr);
  Core::Geo::Cut::Output::GmshEndSection(file, false);

  std::stringstream edge_section_name;
  edge_section_name << "Edge" << suffix << id;
  Core::Geo::Cut::Output::GmshNewSection(file, edge_section_name.str());
  Core::Geo::Cut::Output::GmshEdgeDump(file, edge, false, nullptr);
  Core::Geo::Cut::Output::GmshEndSection(file, false);
}


void Core::Geo::Cut::Output::GmshCutPairDump(
    std::ofstream& file, const std::pair<Side*, Edge*>& pair, int id, const std::string& suffix)
{
  GmshCutPairDump(file, pair.first, pair.second, id, suffix);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given edge                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshEdgeDump(
    std::ofstream& file, Core::Geo::Cut::Edge* edge, bool to_local, Element* ele)
{
#if CUT_DEVELOP
  GmshLineDump(file, edge->BeginNode()->point(), edge->EndNode()->point(),
      edge->BeginNode()->point()->Id(), edge->EndNode()->point()->Id(), to_local, ele);
#else
  GmshLineDump(file, edge->BeginNode()->point(), edge->EndNode()->point(), edge->BeginNode()->Id(),
      edge->EndNode()->Id(), to_local, ele);
#endif
}

void Core::Geo::Cut::Output::GmshEdgeDump(std::ofstream& file, Core::Geo::Cut::Edge* edge,
    const std::string& ename, bool to_local, Element* ele)
{
  GmshNewSection(file, ename, false);
  GmshEdgeDump(file, edge, to_local, ele);
  file << "};\n";
}


/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given node                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshNodeDump(
    std::ofstream& file, Core::Geo::Cut::Node* node, bool to_local, Element* ele)
{
  GmshPointDump(file, node->point(), node->Id(), to_local, ele);
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given point                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshPointDump(
    std::ofstream& file, Core::Geo::Cut::Point* point, int idx, bool to_local, Element* ele)
{
  file << "SP (";
  GmshWriteCoords(file, point, to_local, ele);
  file << "){";
  file << idx;
  file << "};\n";
}



void Core::Geo::Cut::Output::GmshPointDump(std::ofstream& file, Core::Geo::Cut::Point* point,
    int idx, const std::string& pname, bool to_local, Element* ele)
{
  GmshNewSection(file, pname, false);
  GmshPointDump(file, point, idx, to_local, ele);
  file << "};\n";
}


/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given point                                           ager 04/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshPointDump(
    std::ofstream& file, Core::Geo::Cut::Point* point, bool to_local, Element* ele)
{
  GmshPointDump(file, point, point->Position(), to_local, ele);
}

/*--------------------------------------------------------------------------------------*
 * Write Level Set Gradient for given Element
 *
 * The gradients are written at the midpoint of the facets and if the facet is triangulated,
 * also in the midpoint of the triangles.
 *                                                                           winter 07/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshLevelSetGradientDump(
    std::ofstream& file, Element* ele, bool to_local)
{
  const plain_facet_set facets = ele->Facets();
  for (plain_facet_set::const_iterator j = facets.begin(); j != facets.end(); j++)
  {
    Facet* facet = *j;

    std::vector<double> normal_triag_midp;
    if (facet->OnCutSide())
    {
      Core::LinAlg::Matrix<3, 1> facet_triang_midpoint_coord(true);

      if (facet->IsTriangulated())
      {
        std::vector<std::vector<Point*>> facet_triang = facet->Triangulation();
        Point* facet_triang_midpoint = (facet_triang[0])[0];

        facet_triang_midpoint->Coordinates(&facet_triang_midpoint_coord(0, 0));
        normal_triag_midp = ele->get_level_set_gradient(facet_triang_midpoint_coord);

        for (std::vector<std::vector<Point*>>::iterator k = facet_triang.begin();
             k != facet_triang.end(); k++)
        {
          std::vector<Point*> facet_triang_tri = *k;

          Core::LinAlg::Matrix<3, 1> cur;
          Core::LinAlg::Matrix<3, 1> f_triang_tri_midp(true);
          for (std::vector<Point*>::iterator i = facet_triang_tri.begin();
               i != facet_triang_tri.end(); i++)
          {
            Point* p1 = *i;
            p1->Coordinates(cur.data());
            f_triang_tri_midp.update(1.0, cur, 1.0);
          }
          f_triang_tri_midp.scale(1.0 / facet_triang_tri.size());

          std::vector<double> normal = ele->get_level_set_gradient(f_triang_tri_midp);

          GmshVector(file, f_triang_tri_midp, normal, true, to_local, ele);
        }
      }
      else
      {
        Core::LinAlg::Matrix<3, 1> cur;
        std::vector<Point*> pts = facet->Points();
        for (std::vector<Point*>::iterator i = pts.begin(); i != pts.end(); i++)
        {
          Point* p1 = *i;
          p1->Coordinates(cur.data());
          facet_triang_midpoint_coord.update(1.0, cur, 1.0);
        }
        facet_triang_midpoint_coord.scale(1.0 / pts.size());
        normal_triag_midp = ele->get_level_set_gradient(facet_triang_midpoint_coord);
      }

      std::vector<double> normal = ele->get_level_set_gradient(facet_triang_midpoint_coord);
      GmshVector(file, facet_triang_midpoint_coord, normal, true, to_local, ele);

      // Write Corner-points of LS:
      std::vector<Point*> cornerpts = facet->CornerPoints();
      for (std::vector<Point*>::iterator i = cornerpts.begin(); i != cornerpts.end(); i++)
      {
        Core::LinAlg::Matrix<3, 1> cornercoord;
        Point* p1 = *i;
        p1->Coordinates(cornercoord.data());
        std::vector<double> normal = ele->get_level_set_gradient(cornercoord);

        GmshVector(file, cornercoord, normal, true, to_local, ele);
      }
    }
  }
}

/*--------------------------------------------------------------------------------------*
 * Write Level Set Values for given Element
 *
 * The LS-value written at the midpoint of the facets and if the facet is triangulated,
 * also in the midpoint of the triangles.
 * Values at the nodes are also written.
 *                                                                           winter 07/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshLevelSetValueDump(
    std::ofstream& file, Element* ele, bool dumpnodevalues, bool to_local)
{
  const plain_facet_set facets = ele->Facets();
  for (plain_facet_set::const_iterator j = facets.begin(); j != facets.end(); j++)
  {
    Facet* facet = *j;

    if (facet->OnCutSide())
    {
      Core::LinAlg::Matrix<3, 1> facet_triang_midpoint_coord(true);

      if (facet->IsTriangulated())
      {
        std::vector<std::vector<Point*>> facet_triang = facet->Triangulation();
        for (std::vector<std::vector<Point*>>::iterator k = facet_triang.begin();
             k != facet_triang.end(); k++)
        {
          std::vector<Point*> facet_triang_tri = *k;

          Core::LinAlg::Matrix<3, 1> cur;
          Core::LinAlg::Matrix<3, 1> f_triang_tri_midp(true);
          for (std::vector<Point*>::iterator i = facet_triang_tri.begin();
               i != facet_triang_tri.end(); i++)
          {
            Point* p1 = *i;
            p1->Coordinates(cur.data());
            f_triang_tri_midp.update(1.0, cur, 1.0);
          }
          f_triang_tri_midp.scale(1.0 / facet_triang_tri.size());

          double ls_value = ele->get_level_set_value(f_triang_tri_midp);
          GmshScalar(file, f_triang_tri_midp, ls_value, to_local, ele);
        }
        Point* facet_triang_midpoint = (facet_triang[0])[0];
        facet_triang_midpoint->Coordinates(&facet_triang_midpoint_coord(0, 0));
      }
      else
      {
        Core::LinAlg::Matrix<3, 1> cur;
        std::vector<Point*> pts = facet->Points();
        for (std::vector<Point*>::iterator i = pts.begin(); i != pts.end(); i++)
        {
          Point* p1 = *i;
          p1->Coordinates(cur.data());
          facet_triang_midpoint_coord.update(1.0, cur, 1.0);
        }
        facet_triang_midpoint_coord.scale(1.0 / pts.size());
      }

      double ls_value = ele->get_level_set_value(facet_triang_midpoint_coord);
      GmshScalar(file, facet_triang_midpoint_coord, ls_value, to_local, ele);
    }
  }

  if (dumpnodevalues)
  {
    std::vector<Node*> nodes = ele->Nodes();
    for (std::vector<Node*>::iterator j = nodes.begin(); j != nodes.end(); j++)
    {
      Node* node = *j;
      Core::LinAlg::Matrix<3, 1> node_coord(true);
      node->Coordinates(&node_coord(0, 0));

      GmshScalar(file, node_coord, node->LSV(), to_local, ele);
    }
  }
}

/*--------------------------------------------------------------------------------------*
 * Write GMSH output of given coord as point                                   ager 02/17
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshCoordDump(
    std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord, double idx, bool to_local, Element* ele)
{
  file << "SP (";
  GmshWriteCoords(file, coord, to_local, ele);
  file << "){";
  file << idx;
  file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * Write Level Set Values for given Element
 *
 * The LS-value written at the midpoint of the facets and if the facet is triangulated,
 * also in the midpoint of the triangles.
 * Values at the nodes are also written.
 *                                                                           winter 07/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshLevelSetValueZeroSurfaceDump(
    std::ofstream& file, Element* ele, bool to_local)
{
  std::vector<double> lsv_value(8);

  std::vector<Node*> nodes = ele->Nodes();
  int mm = 0;
  for (std::vector<Node*>::iterator j = nodes.begin(); j != nodes.end(); j++)
  {
    Node* node = *j;
    lsv_value[mm] = node->LSV();
    mm++;
  }

  double lsv_max = lsv_value[0];
  double lsv_min = lsv_value[0];

  for (unsigned l = 1; l < lsv_value.size(); l++)
  {
    if (lsv_max < lsv_value[l]) lsv_max = lsv_value[l];

    if (lsv_min > lsv_value[l]) lsv_min = lsv_value[l];
  }

  // localcoord [-1,-1,-1] x [1,1,1]
  int z_sp = 150;
  int y_sp = 150;
  int x_sp = 150;


  double fac = 5 * 1e-3;
  double tolerance = (lsv_max - lsv_min) * fac;  //(0.001;

  //  double* x(3);
  Core::LinAlg::Matrix<3, 1> coord;

  for (int i = 0; i < x_sp; i++)
  {
    // std::cout << "i: " << i << std::endl;
    coord(0, 0) = -1.0 + (2.0 / (double(x_sp) - 1)) * double(i);
    for (int j = 0; j < y_sp; j++)
    {
      // std::cout << "j: " << j << std::endl;
      coord(1, 0) = -1.0 + (2.0 / (double(y_sp) - 1)) * double(j);

      for (int k = 0; k < z_sp; k++)
      {
        // std::cout << "k: " << k << std::endl;
        coord(2, 0) = -1.0 + (2.0 / (double(z_sp) - 1)) * double(k);

        double ls_value = ele->get_level_set_value_at_local_coords(coord);
        if (fabs(ls_value) < tolerance)
        {
          Core::LinAlg::Matrix<3, 1> coord_global;
          ele->global_coordinates(coord, coord_global);
          GmshScalar(file, coord_global, ls_value, to_local, ele);
        }
      }
    }
  }
}


/*--------------------------------------------------------------------------------------*
 * Write Level Set Gradient Orientation of Boundary-Cell Normal and LevelSet
 *                                                                           winter 07/15
 *--------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshLevelSetOrientationDump(
    std::ofstream& file, Element* ele, bool to_local)
{
  const plain_volumecell_set volcells = ele->VolumeCells();
  for (plain_volumecell_set::const_iterator i = volcells.begin(); i != volcells.end(); i++)
  {
    VolumeCell* volcell = *i;

    if (volcell->Position() == Core::Geo::Cut::Point::inside) continue;

    //    const plain_facet_set facets = volcells->Facets();

    volcell->BoundaryCells();
    plain_boundarycell_set bc_cells = volcell->BoundaryCells();
    for (plain_boundarycell_set::iterator j = bc_cells.begin(); j != bc_cells.end(); ++j)
    {
      BoundaryCell* bc = *j;

      //      Facet *facet = *bc->GetFacet();
      Core::LinAlg::Matrix<3, 1> midpoint_bc;
      bc->element_center(midpoint_bc);

      Core::LinAlg::Matrix<3, 1> normal_bc;
      Core::LinAlg::Matrix<2, 1> xsi;
      bc->Normal(xsi, normal_bc);

      std::vector<std::vector<double>> coords_bc = bc->CoordinatesV();
      // const Core::LinAlg::SerialDenseMatrix ls_coordEp = bc->Coordinates();
      Core::LinAlg::Matrix<3, 1> ls_coord(true);
      ls_coord(0, 0) = coords_bc[1][0];
      ls_coord(1, 0) = coords_bc[1][1];
      ls_coord(2, 0) = coords_bc[1][2];

      std::vector<double> normal_ls = ele->get_level_set_gradient(ls_coord);

      double dotProduct = 0.0;
      for (unsigned d = 0; d < normal_ls.size(); ++d) dotProduct += normal_ls[d] * normal_bc(d, 0);

      double normalized_dotproduct = 0;
      if (dotProduct < 1e-15)
        normalized_dotproduct = 0;
      else
        normalized_dotproduct = dotProduct / fabs(dotProduct);

      GmshScalar(file, midpoint_bc, normalized_dotproduct, to_local, ele);
    }
  }
}


/*!
\brief Write Eqn of plane normal for facet (used for DirectDivergence).
 */
void Core::Geo::Cut::Output::GmshEqnPlaneNormalDump(
    std::ofstream& file, Element* ele, bool normalize, bool to_local)
{
  const plain_facet_set facets = ele->Facets();
  for (plain_facet_set::const_iterator j = facets.begin(); j != facets.end(); j++)
  {
    Facet* facet = *j;
    GmshEqnPlaneNormalDump(file, facet, normalize, to_local, ele);
  }
}

/*!
\brief Write Eqn of plane normal for all facets (used for DirectDivergence).
 */
void Core::Geo::Cut::Output::GmshEqnPlaneNormalDump(
    std::ofstream& file, Facet* facet, bool normalize, bool to_local, Element* ele)
{
  Core::LinAlg::Matrix<3, 1> facet_triang_midpoint_coord(true);
  std::vector<Point*> f_cornpts = facet->CornerPoints();
  std::vector<double> eqn_plane = GetEqOfPlane(f_cornpts);

  if (facet->IsTriangulated())
  {
    std::vector<std::vector<Point*>> facet_triang = facet->Triangulation();
    Point* facet_triang_midpoint = (facet_triang[0])[0];
    facet_triang_midpoint->Coordinates(&facet_triang_midpoint_coord(0, 0));

    for (std::vector<std::vector<Point*>>::iterator k = facet_triang.begin();
         k != facet_triang.end(); k++)
    {
      std::vector<Point*> facet_triang_tri = *k;

      Core::LinAlg::Matrix<3, 1> cur;
      Core::LinAlg::Matrix<3, 1> f_triang_tri_midp(true);
      for (std::vector<Point*>::iterator i = facet_triang_tri.begin(); i != facet_triang_tri.end();
           i++)
      {
        Point* p1 = *i;
        p1->Coordinates(cur.data());
        f_triang_tri_midp.update(1.0, cur, 1.0);
      }
      f_triang_tri_midp.scale(1.0 / facet_triang_tri.size());

      GmshVector(file, f_triang_tri_midp, GetEqOfPlane(facet_triang_tri), normalize, to_local, ele);
    }
  }
  else
  {
    Core::LinAlg::Matrix<3, 1> cur;
    std::vector<Point*> pts = facet->Points();
    for (std::vector<Point*>::iterator i = pts.begin(); i != pts.end(); i++)
    {
      Point* p1 = *i;
      p1->Coordinates(cur.data());
      facet_triang_midpoint_coord.update(1.0, cur, 1.0);
    }
    facet_triang_midpoint_coord.scale(1.0 / pts.size());
  }

  GmshVector(file, facet_triang_midpoint_coord, eqn_plane, normalize, to_local, ele);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshScalar(std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord,
    double scalar, bool to_local, Element* ele)  // use gmshpoint?
{
  file << "SP(";
  GmshWriteCoords(file, coord, to_local, ele);
  file << "){";
  file << scalar;
  file << "};\n";
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshVector(std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord,
    std::vector<double> vector, bool normalize, bool to_local, Element* ele)
{
  file << "VP(";
  GmshWriteCoords(file, coord, to_local, ele);
  file << "){";

  if (normalize)
  {
    double norm = 0.0;
    for (unsigned i = 0; i < vector.size(); ++i) norm += vector[i] * vector[i];

    norm = sqrt(norm);
    if (norm > 1e-15)
    {
      for (unsigned i = 0; i < vector.size(); ++i) vector[i] = vector[i] / norm;
    }
    else
    {
      std::fill(vector.begin(), vector.end(), 0.0);
    }
  }
  GmshWriteCoords(file, vector, to_local, ele);
  file << "};\n";
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshWriteCoords(
    std::ofstream& file, std::vector<double> coord, bool to_local, Element* ele)
{
  Core::LinAlg::Matrix<3, 1> xyz(true);

  if (coord.size() <= 3)
    std::copy(coord.begin(), coord.end(), xyz.data());
  else
    FOUR_C_THROW("The coord vector dimension is wrong! (coord.size() = %d)", coord.size());

  if (to_local)
  {
    if (ele == nullptr)
      FOUR_C_THROW("GmshWriteCoords: Didn't get a parent element for the Coordinate!");

    Core::LinAlg::Matrix<3, 1> rst(true);

    ele->local_coordinates(xyz, rst);
    GmshWriteCoords(file, rst, false, nullptr);  // rst are already local coords!
    return;
  }
  file << std::setprecision(15) << xyz(0) << "," << xyz(1) << "," << xyz(2);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshWriteCoords(
    std::ofstream& file, Core::LinAlg::Matrix<3, 1> coord, bool to_local, Element* ele)
{
  if (to_local)
  {
    if (ele == nullptr)
      FOUR_C_THROW("GmshWriteCoords: Didn't get a parent element for the Coordinate!");

    Core::LinAlg::Matrix<3, 1> xyz = coord;
    ele->local_coordinates(xyz, coord);
  }
  file << std::setprecision(15) << coord(0, 0) << "," << coord(1, 0) << "," << coord(2, 0);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshWriteCoords(
    std::ofstream& file, Node* node, bool to_local, Element* ele)
{
  GmshWriteCoords(file, node->point(), to_local, ele);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshWriteCoords(
    std::ofstream& file, Point* point, bool to_local, Element* ele)
{
  Core::LinAlg::Matrix<3, 1> coord;
  point->Coordinates(coord.data());

  GmshWriteCoords(file, coord, to_local, ele);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string Core::Geo::Cut::Output::GenerateGmshOutputFilename(const std::string& filename_tail)
{
  //  std::string filename = ::Global::Problem::Instance()->OutputControlFile()->file_name();
  std::string filename("xxx");
  filename.append(filename_tail);
  return filename;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshNewSection(
    std::ofstream& file, const std::string& section, bool first_endsection)
{
  if (first_endsection) file << "};\n";
  file << "View \"" << section << "\" {\n";
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshEndSection(std::ofstream& file, bool close_file)
{
  file << "};\n";
  if (close_file) file.close();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<double> Core::Geo::Cut::Output::GetEqOfPlane(std::vector<Point*> pts)
{
  int mm = 0;

  std::vector<std::vector<double>> corners(pts.size());

  for (std::vector<Point*>::iterator k = pts.begin(); k != pts.end(); k++)
  {
    Point* p1 = *k;
    Core::LinAlg::Matrix<3, 1> cur;
    p1->Coordinates(cur.data());

    std::vector<double> pt(3);

    pt[0] = cur(0, 0);
    pt[1] = cur(1, 0);
    pt[2] = cur(2, 0);

    corners[mm] = pt;
    mm++;
  }
  return Kernel::EqnPlaneOfPolygon(corners);
}


/*-------------------------------------------------------------------------------*
 * Write cuttest for this element!                                     ager 04/15
 *-------------------------------------------------------------------------------*/
void Core::Geo::Cut::Output::GmshElementCutTest(
    std::ofstream& file, Core::Geo::Cut::Element* ele, bool haslevelsetside)
{
  std::cout << "Write Cut Test for Element " << ele->Id() << " ... " << std::flush;

  // default precision for coordinates
  file << std::setprecision(20);
  // -- 1 -- header of cut_test -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  file << "// This test was automatically generated by Cut::Output::GmshElementCutTest(), "
       << "\n";
  file << "// as the cut crashed for this configuration!"
       << "\n";
  file << ""
       << "\n";
  file << "#include <iostream>"
       << "\n";
  file << "#include <map>"
       << "\n";
  file << "#include <string>"
       << "\n";
  file << "#include <vector>"
       << "\n";
  file << ""
       << "\n";
  file << "#include \"cut_test_utils.H\""
       << "\n";
  file << ""
       << "\n";
  file << "#include \"../../src/cut/cut_side.H\""
       << "\n";
  file << "#include \"../../src/cut/cut_meshintersection.H\""
       << "\n";
  file << "#include \"../../src/cut/cut_levelsetintersection.H\""
       << "\n";
  file << "#include \"../../src/cut/cut_combintersection.H\""
       << "\n";
  file << "#include \"../../src/cut/cut_tetmeshintersection.H\""
       << "\n";
  file << "#include \"../../src/cut/cut_options.H\""
       << "\n";
  file << "#include \"../../src/cut/cut_volumecell.H\""
       << "\n";
  file << ""
       << "\n";
  file << "#include \"../../src/fem_general/utils_local_connectivity_matrices.H\""
       << "\n";
  file << ""
       << "\n";
  file << "void test_generated_" << ele->Id() << "()"
       << "\n";
  file << "{"
       << "\n";
  if (not haslevelsetside)
    file << "  Core::Geo::Cut::MeshIntersection intersection;"
         << "\n";
  else
    file << "  Core::Geo::Cut::CombIntersection intersection(-1);"
         << "\n";
  file << "  intersection.GetOptions().Init_for_Cuttests();  // use full cln\n";
  file << "  std::vector<int> nids;"
       << "\n";
  file << ""
       << "\n";
  file << "  int sidecount = 0;"
       << "\n";
  file << "  std::vector<double> lsvs(" << ele->Nodes().size() << ");"
       << "\n";

  // --- get all neighboring elements and cutsides
  plain_element_set eles;
  plain_side_set csides;
  for (std::vector<Side*>::const_iterator sit = ele->Sides().begin(); sit != ele->Sides().end();
       ++sit)
  {
    for (plain_element_set::const_iterator eit = (*sit)->Elements().begin();
         eit != (*sit)->Elements().end(); ++eit)
    {
      eles.insert(*eit);
      for (plain_side_set::const_iterator csit = (*eit)->CutSides().begin();
           csit != (*eit)->CutSides().end(); ++csit)
        csides.insert(*csit);
    }
  }


  if (not haslevelsetside)
  {
    // -- 2 -- add sides -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    const plain_side_set& cutsides = csides;
    for (plain_side_set::const_iterator i = cutsides.begin(); i != cutsides.end(); ++i)
    {
      file << "  {"
           << "\n";
      file << "    Core::LinAlg::SerialDenseMatrix tri3_xyze( 3, 3 );"
           << "\n";
      file << ""
           << "\n";
      Side* s = *i;
      const std::vector<Node*>& side_nodes = s->Nodes();
      int nodelid = -1;
      file << "    nids.clear();"
           << "\n";
      for (std::vector<Node*>::const_iterator j = side_nodes.begin(); j != side_nodes.end(); ++j)
      {
        nodelid++;
        Node* n = *j;
        for (int dim = 0; dim < 3; ++dim)
        {
          file << "    tri3_xyze(" << dim << "," << nodelid << ") = " << n->point()->X()[dim] << ";"
               << "\n";
        }
        file << "    nids.push_back( " << n->Id() << " );"
             << "\n";
      }
      file << "    intersection.AddCutSide( ++sidecount, nids, tri3_xyze, "
              "Core::FE::CellType::tri3 );"
           << "\n";
      file << "  }"
           << "\n";
    }
  }
  else
  {
    file << "  ci.AddLevelSetSide(1);"
         << "\n";
    for (std::size_t i = 0; i < ele->Nodes().size(); ++i)
    {
      file << "  lsvs[" << i << "] = " << ele->Nodes()[i]->LSV() << ";"
           << "\n";
    }
  }

  // -- 3 -- add background element -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  const plain_element_set& elements = eles;
  for (plain_element_set::const_iterator eit = elements.begin(); eit != elements.end(); ++eit)
  {
    Element* aele = *eit;
    file << "  {"
         << "\n";
    file << "  LinAlg::SerialDenseMatrix hex" << aele->Nodes().size() << "_xyze( 3, "
         << aele->Nodes().size() << " );"
         << "\n";
    file << ""
         << "\n";
    file << "    nids.clear();"
         << "\n";
    for (std::size_t i = 0; i < aele->Nodes().size(); ++i)
    {
      for (unsigned dim = 0; dim < 3; ++dim)
      {
        file << "  hex8_xyze(" << dim << "," << i << ") = " << aele->Nodes()[i]->point()->X()[dim]
             << ";"
             << "\n";
      }
      file << "  nids.push_back( " << aele->Nodes()[i]->Id() << " );"
           << "\n";
    }
    file << ""
         << "\n";
    if (not haslevelsetside)
      file << "  intersection.add_element( " << aele->Id()
           << ", nids, hex8_xyze, Core::FE::CellType::hex8);"
           << "\n";
    else
      file << "  intersection.add_element( " << aele->Id()
           << ", nids, hex8_xyze, Core::FE::CellType::hex8, &lsvs[0], false );"
           << "\n";
    file << "  }"
         << "\n";
    file << ""
         << "\n";
  }
  file << "  intersection.CutTest_Cut( true,Inpar::Cut::VCellGaussPts_DirectDivergence, "
          "Inpar::Cut::BCellGaussPts_Tessellation );"
       << "\n";
  file << "  intersection.Cut_Finalize( true, Inpar::Cut::VCellGaussPts_DirectDivergence, "
          "Inpar::Cut::BCellGaussPts_Tessellation, false, true );"
       << "\n";
  file << ""
       << "\n";

  if (not haslevelsetside && 0)
  {
    // -- 4 -- compare integration methods -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    file << "  std::vector<double> tessVol,momFitVol,dirDivVol;"
         << "\n";
    file << ""
         << "\n";
    file << "  Core::Geo::Cut::Mesh mesh = intersection.NormalMesh();"
         << "\n";
    file << "  const std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell> > & other_cells = "
            "mesh.VolumeCells();"
         << "\n";
    file << "  for ( std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell> >::const_iterator "
            "i=other_cells.begin();"
         << "\n";
    file << "        i!=other_cells.end();"
         << "\n";
    file << "        ++i )"
         << "\n";
    file << "  {"
         << "\n";
    file << "    Core::Geo::Cut::VolumeCell * vc = &**i;"
         << "\n";
    file << "    tessVol.push_back(vc->Volume());"
         << "\n";
    file << "  for ( std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell> >::const_iterator "
            "i=other_cells.begin();"
         << "\n";
    file << "        i!=other_cells.end();"
         << "\n";
    file << "        ++i )"
         << "\n";
    file << "  {"
         << "\n";
    file << "    Core::Geo::Cut::VolumeCell * vc = &**i;"
         << "\n";
    file << "    "
            "vc->moment_fit_gauss_weights(vc->parent_element(),mesh,true,Inpar::Cut::BCellGaussPts_"
            "Tessellation);"
         << "\n";
    file << "    momFitVol.push_back(vc->Volume());"
         << "\n";
    file << "  }"
         << "\n";
    file << ""
         << "\n";
    file << "  for ( std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell> >::const_iterator "
            "i=other_cells.begin();"
         << "\n";
    file << "           i!=other_cells.end();"
         << "\n";
    file << "           ++i )"
         << "\n";
    file << "   {"
         << "\n";
    file << "     Core::Geo::Cut::VolumeCell * vc = &**i;"
         << "\n";
    file << "     "
            "vc->direct_divergence_gauss_rule(vc->parent_element(),mesh,true,Inpar::Cut::"
            "BCellGaussPts_"
            "Tessellation);"
         << "\n";
    file << "     dirDivVol.push_back(vc->Volume());"
         << "\n";
    file << "   }"
         << "\n";
    file << ""
         << "\n";
    file << "  std::cout<<\"the volumes predicted by\\n tessellation \\t MomentFitting \\t "
            "DirectDivergence\\n\";"
         << "\n";
    file << "  for(unsigned i=0;i<tessVol.size();i++)"
         << "\n";
    file << "  {"
         << "\n";
    file << "    std::cout<<tessVol[i]<<\"\\t\"<<momFitVol[i]<<\"\\t\"<<dirDivVol[i]<<\"\\n\";"
         << "\n";
    file << "    if( fabs(tessVol[i]-momFitVol[i])>1e-9 || fabs(dirDivVol[i]-momFitVol[i])>1e-9 )"
         << "\n";
    file << "    {"
         << "\n";
    file << "      mesh.DumpGmsh(\"Cuttest_Debug_Output.pos\");"
         << "\n";
    file << "      intersection.CutMesh().GetElement(1)->DebugDump();"
         << "\n";
    file << "      FOUR_C_THROW(\"volume predicted by either one of the method is wrong\");"
         << "\n";
    file << "      }"
         << "\n";
    file << "    }"
         << "\n";
  }
  file << "}"
       << "\n";
  std::cout << "done " << std::endl;
}

// Debug related functions

void Core::Geo::Cut::Output::DebugDump_ThreePointsOnEdge(
    Side* first, Side* second, Edge* e, Point* p, const PointSet& cut)
{
  std::ofstream file("three_points_on_the_edge.pos");
  Core::Geo::Cut::Output::GmshSideDump(file, first, std::string("FirstSide"));
  Core::Geo::Cut::Output::GmshSideDump(file, second, std::string("SecondSide"));
  file << "// Edge is " << e->BeginNode()->point() << "->" << e->EndNode()->point() << std::endl;
  Core::Geo::Cut::Output::GmshEdgeDump(file, e, std::string("NotTouchedEdge"));
  Core::Geo::Cut::Output::GmshPointDump(
      file, p, p->Id(), std::string("NotTouchedPoint"), false, nullptr);
  for (PointSet::const_iterator it = cut.begin(); it != cut.end(); ++it)
  {
    std::stringstream point_name;
    point_name << "CutPoint" << std::distance(cut.begin(), it);
    Core::Geo::Cut::Output::GmshPointDump(file, *it, (*it)->Id(), point_name.str(), false, nullptr);
    (*it)->dump_connectivity_info();
  }

  const PointPositionSet& edge_cut_points = e->CutPoints();
  for (PointPositionSet::const_iterator it = edge_cut_points.begin(); it != edge_cut_points.end();
       ++it)
  {
    std::stringstream point_name;
    point_name << "EdgeCutPoint" << std::distance(edge_cut_points.begin(), it);
    Core::Geo::Cut::Output::GmshPointDump(file, *it, (*it)->Id(), point_name.str(), false, nullptr);
  }
  file.close();
}

void Core::Geo::Cut::Output::DebugDump_MoreThanTwoIntersectionPoints(
    Edge* edge, Side* other, const std::vector<Point*>& point_stack)
{
  std::ofstream file("multiple_intersection_points.pos");
  // dump everything
  for (std::vector<Point*>::const_iterator it = point_stack.begin(); it != point_stack.end(); ++it)
  {
    std::stringstream section_name;
    section_name << "Point" << (*it)->Id();
    Core::Geo::Cut::Output::GmshNewSection(file, section_name.str());
    Core::Geo::Cut::Output::GmshPointDump(file, *it, (*it)->Id(), false, nullptr);
    Core::Geo::Cut::Output::GmshEndSection(file, false);
#if CUT_CREATION_INFO
    const std::pair<Side*, Edge*>& cu_pair = std::make_pair(other, edge);
    const std::pair<Side*, Edge*>& or_pair = (*it)->AddedFrom(cu_pair);
    Core::Geo::Cut::Output::GmshCutPairDump(file, or_pair, 0, std::string("added_from"));
    if (or_pair != cu_pair)
    {
      file << "// original added becase: " << (*it)->GetCreationInfo(or_pair) << "\n";
      file << "// common added because: " << (*it)->GetCreationInfo(cu_pair) << "\n";
    }
    else
      file << "// original added becase: " << (*it)->GetCreationInfo(cu_pair) << "\n";
#endif

    (*it)->dump_connectivity_info();
  }
  // Compute all possible differences between points

  for (std::vector<Point*>::const_iterator it = point_stack.begin(); it != point_stack.end(); ++it)
  {
    for (std::vector<Point*>::const_iterator jt = point_stack.begin(); jt != point_stack.end();
         ++jt)
    {
      if (*it != *jt)
      {
        std::cout << "Diference between " << (*it)->Id() << " and " << (*jt)->Id() << " is "
                  << Core::Geo::Cut::DistanceBetweenPoints(*jt, *it) << std::endl;
      }
    }
  }
  Core::Geo::Cut::Output::GmshCutPairDump(file, other, edge, 0, std::string("common"));
  file.close();
}

void Core::Geo::Cut::Output::DebugDump_MultipleCutPointsSpecial(Side* first, Side* second,
    const PointSet& cut, const PointSet& collected_points, const point_line_set& new_lines)
{
  std::ofstream file("special_case_multiple_cut_points.pos");
  Core::Geo::Cut::Output::GmshSideDump(file, first, std::string("ThisSide"));
  Core::Geo::Cut::Output::GmshSideDump(file, second, std::string("OtherSide"));
  for (PointSet::const_iterator it = cut.begin(); it != cut.end(); ++it)
  {
    std::stringstream pname;
    pname << "CutPoint" << std::distance(cut.begin(), it);
    Core::Geo::Cut::Output::GmshPointDump(file, *it, (*it)->Id(), pname.str(), false, nullptr);
    (*it)->dump_connectivity_info();
  }
  file.close();

  std::cout << "Collected " << collected_points.size() << " points " << std::endl;
  for (PointSet::const_iterator it = collected_points.begin(); it != collected_points.end(); ++it)
  {
    std::cout << (*it)->Id() << "  " << std::endl;
  }
  std::cout << "Totally there are  " << cut.size() << " points " << std::endl;
  for (PointSet::const_iterator it = cut.begin(); it != cut.end(); ++it)
  {
    std::cout << (*it)->Id() << "  " << std::endl;
  }
  std::cout << "Got " << new_lines.size() << " cut lines"
            << " for " << cut.size() << " cut points\n";
  std::cout << "Cut lines are: \n";
  for (point_line_set::const_iterator it = new_lines.begin(); it != new_lines.end(); ++it)
  {
    std::cout << it->first << "--" << it->second << std::endl;
  }
}

FOUR_C_NAMESPACE_CLOSE
