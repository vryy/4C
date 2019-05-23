/*---------------------------------------------------------------------*/
/*!

\brief Cut Volumecell (3d-Object described by 2d surfaces(facets))

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include "cut_point.H"
#include "cut_cycle.H"
#include "cut_volumecell.H"
#include "cut_boundarycell.H"
#include "cut_integrationcell.H"
#include "cut_tetmesh.H"
#include "cut_mesh.H"
#include "cut_options.H"
#include "cut_kernel.H"
#include "cut_triangulateFacet.H"
#include "volume_integration.H"
#include "boundarycell_integration.H"
#include "direct_divergence.H"
#include "cut_utils.H"

#include <Teuchos_TimeMonitor.hpp>


#include <algorithm>


int GEO::CUT::VolumeCell::hex8totet4[5][4] = {
    {0, 1, 3, 4}, {1, 2, 3, 6}, {4, 5, 1, 6}, {6, 7, 3, 4}, {1, 6, 3, 4}};

int GEO::CUT::VolumeCell::wedge6totet4[3][4] = {{0, 1, 2, 3}, {3, 4, 1, 5}, {1, 5, 2, 3}};


int GEO::CUT::VolumeCell::pyramid5totet4[2][4] = {{0, 1, 3, 4}, {1, 2, 3, 4}};


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::VolumeCell::VolumeCell(const plain_facet_set& facets,
    const std::map<std::pair<Point*, Point*>, plain_facet_set>& volume_lines, Element* element)
    : element_(element),
      position_(Point::undecided),
      facets_(facets),
      volume_(0.0),
      isNegligibleSmall_(false)
{
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->Register(this);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::VolumeCell::IsEqual(const plain_facet_set& vcell) const
{
  bool isequal = false;
  for (plain_facet_set::const_iterator ci = facets_.begin(); ci != facets_.end(); ++ci)
  {
    const Facet& f = **ci;
    for (plain_facet_set::const_iterator cii = vcell.begin(); cii != vcell.end(); ++cii)
    {
      const Facet& ff = **cii;
      const std::vector<Point*>& ffpoints = ff.Points();

      isequal = (f.Points().size() == ffpoints.size() and f.Contains(ffpoints));
      if (isequal) break;
    }

    if (not isequal) return false;
  }

  // if it reaches this point, the volume cells are identical
  return isequal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::Neighbors(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& done, plain_volumecell_set& connected, plain_element_set& elements)
{
  if (done.count(this) == 0)
  {
    // this volume is included
    connected.insert(this);
    elements.insert(element_);

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
    {
      Facet* f = *i;
      if (p == NULL or f->Contains(p))
      {
        f->Neighbors(p, cells, done, connected, elements);
      }
    }

    if (p != NULL)
    {
      for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
      {
        Facet* f = *i;
        if (not f->Contains(p))
        {
          f->Neighbors(p, cells, done, connected, elements);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
// without check for elements
void GEO::CUT::VolumeCell::Neighbors(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& done, plain_volumecell_set& connected)
{
  if (done.count(this) == 0)
  {
    // this volume is included
    connected.insert(this);

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
    {
      Facet* f = *i;

      if (p == NULL or f->Contains(p))
      {
        f->Neighbors(p, cells, done, connected);
      }
    }

    if (p != NULL)
    {
      for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
      {
        Facet* f = *i;
        if (not f->Contains(p))
        {
          f->Neighbors(p, cells, done, connected);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::GetAllPoints(Mesh& mesh, PointSet& cut_points)
{
  for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->GetAllPoints(mesh, cut_points);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::VolumeCell::Contains(Point* p)
{
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    if (f->Contains(p))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::VolumeCell::Contains(LINALG::Matrix<3, 1>& x)
{
  if (integrationcells_.size() == 0)
    dserror(
        "no integrationcells for volumecell stored, implement Contains check without "
        "integrationcells");

  for (GEO::CUT::plain_integrationcell_set::iterator it = integrationcells_.begin();
       it != integrationcells_.end(); it++)
  {
    GEO::CUT::IntegrationCell* intcell = *it;

    if (intcell->Contains(x)) return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::CreateTet4IntegrationCells(Mesh& mesh,
    const std::vector<std::vector<Point*>>& tets,
    const std::map<Facet*, std::vector<Point*>>& sides_xyz)
{
  for (std::vector<std::vector<Point*>>::const_iterator i = tets.begin(); i != tets.end(); ++i)
  {
    const std::vector<Point*>& tet = *i;
    if (tet.size() != 4)
    {
      throw std::runtime_error("tet expected");
    }
    NewTet4Cell(mesh, tet);
  }

  for (std::map<Facet*, std::vector<Point*>>::const_iterator i = sides_xyz.begin();
       i != sides_xyz.end(); ++i)
  {
    Facet* f = i->first;
    const std::vector<Point*>& points = i->second;

    std::size_t length = points.size();
    if (length % 3 != 0) throw std::runtime_error("expect list of triangles");

    length /= 3;
    std::vector<Point*> p(3);
    for (std::size_t i = 0; i < length; ++i)  // loop the list of triangles
    {
      std::copy(&points[3 * i], &points[3 * (i + 1)], &p[0]);
      // Tri3BoundaryCell::CreateCell( mesh, this, f, p );
      NewTri3Cell(mesh, f, p);  // create tri3 cell
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::GetIntegrationCells(plain_integrationcell_set& cells)
{
  std::copy(
      integrationcells_.begin(), integrationcells_.end(), std::inserter(cells, cells.begin()));
}

/// get a map of boundary cells for all cutting sides, key= side-Id, value= vector of boundary cells
/// note that the boundary cells of subsides with the same side id are stored now in one key
void GEO::CUT::VolumeCell::GetBoundaryCells(
    std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells)
{
  for (plain_boundarycell_set::iterator i = bcells_.begin(); i != bcells_.end(); ++i)
  {
    BoundaryCell* bc = *i;
    Facet* f = bc->GetFacet();
    if (f->OnCutSide())
    {
      int sid = f->SideId();  // f->OnCutSide => sid>-1
      // usually there are more facets with the same side id as the cutting sides have been
      // subdivided into subsides which produce own facets for each subside
      bcells[sid].push_back(bc);
    }
  }
}


/// get a map of boundary cells for all cutting sides, key= side-Id, value= vector of boundary cells
/// note that the boundary cells of subsides with the same side id are stored now in one key
void GEO::CUT::VolumeCell::GetBoundaryCellsToBeIntegrated(
    std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells)
{
  for (plain_boundarycell_set::iterator i = bcells_.begin(); i != bcells_.end(); ++i)
  {
    BoundaryCell* bc = *i;
    Facet* f = bc->GetFacet();
    // Get all bc's for cuts from only the outside vc's
    //  as to not integrate twice over the same surface
    if ((f->OnCutSide() and Position() == GEO::CUT::Point::outside))
    {
      int sid = f->SideId();  // f->OnCutSide => sid>-1
      // usually there are more facets with the same side id as the cutting sides have been
      // subdivided into subsides which produce own facets for each subside
      bcells[sid].push_back(bc);

      //@ Christoph: Here one could add marked sides for the cut-mesh!
    }
    else if (f->OnMarkedBackgroundSide())
    {
      // Loop over all marked actions and extract bc's for corresponding coupling object.
      for (std::map<GEO::CUT::MarkedActions, int>::iterator markit =
               f->ParentSide()->GetMarkedsidemap().begin();
           markit != f->ParentSide()->GetMarkedsidemap().end(); ++markit)
      {
        if (markit->first == GEO::CUT::mark_and_create_boundarycells)
          bcells[markit->second].push_back(bc);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 * SideId() of Facet (used for timeintegration)
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::CollectCutSides(plain_int_set& cutside_ids)
{
  for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    if (f->OnCutSide())
    {
      int sid = f->SideId();  // f->OnCutSide => sid>-1
      // usually there are more facets with the same side id as the cutting sides have been
      // subdivided into subsides which produce own facets for each subside
      cutside_ids.insert(sid);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::VolumeCell::GetParentElementId() const { return element_->GetParentId(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::ConnectNodalDOFSets(bool include_inner)
{
  //   if ( Empty() )
  //     return;
  if (not include_inner and Position() != Point::outside) return;

  const std::vector<Node*>& nodes = element_->Nodes();
  nodaldofset_.reserve(nodes.size());

  for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Node* n = *i;
    nodaldofset_.push_back(n->DofSetNumber(this));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::Position(Point::PointPosition position)
{
  if (position_ != position)
  {
    position_ = position;

    for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
    {
      Facet* f = *i;
      Point::PointPosition fp = f->Position();
      if (fp == Point::undecided)
      {
        f->Position(position);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::Print(std::ostream& stream) const
{
  stream << "\n==========================================\n";
  stream << "=== VolumeCell ( address: " << std::setw(10) << this << " ) ===\n";
  stream << "# VolumeCell: "
         << " pos: " << Point::PointPosition2String(position_) << " "
         << "#facets: " << facets_.size() << " "
         << "#intcells: " << integrationcells_.size() << " "
         << "#bcells: " << bcells_.size() << "\n";
  unsigned count = 0;
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    stream << "\n# Facet " << count++ << " of VolumeCell:\n";
    f->Print(stream);
  }

  count = 0;
  for (plain_boundarycell_set::const_iterator i = bcells_.begin(); i != bcells_.end(); ++i)
  {
    BoundaryCell* bcell = *i;
    stream << "\n# BoundaryCell " << count++ << " of VolumeCell:\n";
    bcell->Print(stream);
  }

  count = 0;
  for (plain_integrationcell_set::const_iterator i = integrationcells_.begin();
       i != integrationcells_.end(); ++i)
  {
    IntegrationCell* icell = *i;
    stream << "\n# IntegrationCell " << count++ << " of VolumeCell:\n";
    icell->Print(stream);
  }

  stream << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewBoundaryCell(
    Mesh& mesh, DRT::Element::DiscretizationType shape, Facet* f, const std::vector<Point*>& x)
{
  if (facets_.count(f) == 0)
  {
    run_time_error("facet does not belong to volume cell");
  }
  switch (shape)
  {
    case DRT::Element::point1:
      NewPoint1Cell(mesh, f, x);
      break;
    case DRT::Element::line2:
      NewLine2Cell(mesh, f, x);
      break;
    case DRT::Element::tri3:
      NewTri3Cell(mesh, f, x);
      break;
    case DRT::Element::quad4:
      NewQuad4Cell(mesh, f, x);
      break;
    default:
      dserror("Unsupported shape ( shape = %s )", DRT::DistypeToString(shape).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewPoint1Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->NewPoint1Cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewLine2Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->NewLine2Cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewTri3Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->NewTri3Cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewQuad4Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->NewQuad4Cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewArbitraryCell(Mesh& mesh, Facet* f, const std::vector<Point*>& x,
    const DRT::UTILS::GaussIntegration& gp, const LINALG::Matrix<3, 1>& normal)
{
  f->NewArbitraryCell(mesh, this, x, bcells_, gp, normal);
}

/*double GEO::CUT::VolumeCell::Volume()
{
  double volume = 0;
  for ( plain_integrationcell_set::iterator i=integrationcells_.begin(); i!=integrationcells_.end();
++i )
  {
    IntegrationCell * ic = *i;
    volume += ic->Volume();
  }
  return volume;
}*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int GEO::CUT::VolumeCell::NumGaussPoints(DRT::Element::DiscretizationType shape)
{
  int numgp = 0;

  for (plain_integrationcell_set::const_iterator i = integrationcells_.begin();
       i != integrationcells_.end(); ++i)
  {
    IntegrationCell* ic = *i;

    // Create (unmodified) gauss points for integration cell with requested
    // polynomial order. This is supposed to be fast, since there is a cache.
    DRT::UTILS::GaussIntegration gi(ic->Shape(), ic->CubatureDegree(shape));

    // we just need the number of points per cell
    numgp += gi.NumPoints();
  }

  return numgp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::Disconnect()
{
  for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->DisconnectVolume(this);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewIntegrationCell(
    Mesh& mesh, DRT::Element::DiscretizationType shape, const std::vector<Point*>& x)
{
  switch (shape)
  {
    // --- 1-D elements ---
    case DRT::Element::line2:
      NewLine2Cell(mesh, x);
      break;
    // --- 2-D elements ---
    case DRT::Element::tri3:
      NewTri3Cell(mesh, x);
      break;
    case DRT::Element::quad4:
      NewQuad4Cell(mesh, x);
      break;
    // --- 3-D elements ---
    case DRT::Element::hex8:
      NewHex8Cell(mesh, x);
      break;
    case DRT::Element::tet4:
      NewTet4Cell(mesh, x);
      break;
    case DRT::Element::wedge6:
      NewWedge6Cell(mesh, x);
      break;
    case DRT::Element::pyramid5:
      NewPyramid5Cell(mesh, x);
      break;
    default:
      dserror("Unsupported shape ( shape = %s )", DRT::DistypeToString(shape).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewLine2Cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  integrationcells_.insert(mesh.NewLine2Cell(position, points, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewTri3Cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  integrationcells_.insert(mesh.NewTri3Cell(position, points, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewQuad4Cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  integrationcells_.insert(mesh.NewQuad4Cell(position, points, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewHex8Cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  if (mesh.CreateOptions().GenHex8())
  {
    integrationcells_.insert(mesh.NewHex8Cell(position, points, this));
  }
  else
  {
    std::vector<Point*> tet4_points(4);
    for (int i = 0; i < 5; ++i)
    {
      SetTetPoints(hex8totet4[i], points, tet4_points);
      integrationcells_.insert(mesh.NewTet4Cell(position, tet4_points, this));
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::IntegrationCell* GEO::CUT::VolumeCell::NewTet4Cell(
    Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  IntegrationCell* ic = mesh.NewTet4Cell(position, points, this);
  integrationcells_.insert(ic);
  return ic;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewWedge6Cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  if (mesh.CreateOptions().GenWedge6())
  {
    integrationcells_.insert(mesh.NewWedge6Cell(position, points, this));
  }
  else
  {
    std::vector<Point*> tet4_points(4);
    for (int i = 0; i < 3; ++i)
    {
      SetTetPoints(wedge6totet4[i], points, tet4_points);
      integrationcells_.insert(mesh.NewTet4Cell(position, tet4_points, this));
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::NewPyramid5Cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = Position();
  if (mesh.CreateOptions().GenPyramid5())
  {
    integrationcells_.insert(mesh.NewPyramid5Cell(position, points, this));
  }
  else
  {
    std::vector<Point*> tet4_points(4);
    for (int i = 0; i < 2; ++i)
    {
      SetTetPoints(pyramid5totet4[i], points, tet4_points);
      integrationcells_.insert(mesh.NewTet4Cell(position, tet4_points, this));
    }
  }
}

/*--------------------------------------------------------------------*
 * Check wheter the point is inside, outside or on the boundary
 * of this volumecelll                                    sudhakar 07/12
 *--------------------------------------------------------------------*/
std::string GEO::CUT::VolumeCell::IsThisPointInside(Point* pt)
{
  LINALG::Matrix<3, 1> xglo;
  pt->Coordinates(xglo.A());
  std::string inside = IsThisPointInside(xglo);
  return inside;
}

/*-----------------------------------------------------------------------------------------------*
 * Check whether the point with this global coordinates is inside, outside or on the boundary
 * of this volumecell                                                               sudhakar 07/12
 *-----------------------------------------------------------------------------------------------*/
std::string GEO::CUT::VolumeCell::IsThisPointInside(LINALG::Matrix<3, 1>& xglo)
{
  LINALG::Matrix<3, 1> xloc;
  element_->LocalCoordinates(xglo, xloc);

  const GEO::CUT::Point::PointPosition posi = Position();
  if (posi == 0) dserror("undefined position for the volumecell");

  VolumeIntegration vc(this, element_, posi, 0);
  std::string inside = vc.IsPointInside(xloc);
  return inside;
}

/*
 * There is a flaw in this test. If there are points added through TetMeshIntersection on the facet
 * boundary the test will fail.
 *
 * In order to fix this test one has to check if the remaining lines are connected. I.e. are the
 * points added from the tesselation connected correctly?
 *
 */
void GEO::CUT::VolumeCell::TestSurface()
{
  if (Empty())
  {
    // This is an artificial cell with zero volume. It should not exist in the
    // first place.
    return;
  }

  // see if all lines are closed
  //
  // This finds all the degenerated cases that where dropped before. Thus the
  // test complains a lot.

  for (plain_facet_set::iterator j = facets_.begin(); j != facets_.end(); ++j)
  {
    Facet* f = *j;

    if (f->OnCutSide())
    {
      if (f->IsTriangulated())
      {
        //        std::cout << "f->Triangulation().size(): " << f->Triangulation().size();
        //
      }
      if (f->HasHoles())
      {
        //
      }

      point_line_set lines;

      const std::vector<Point*>& points = f->Points();
      Cycle cycle(points);
      cycle.Add(lines);

      point_line_set facetlines = lines;

      for (plain_boundarycell_set::iterator i = bcells_.begin(); i != bcells_.end(); ++i)
      {
        BoundaryCell* bc = *i;
        if (bc->GetFacet() == f)
        {
          //          std::cout << "Printing boundary cell: " << std::endl;
          //          bc->Print();
          //          std::cout << std::endl;
          const std::vector<Point*>& points = bc->Points();
          Cycle cycle(points);
          cycle.Add(lines);
        }
      }

      if (lines.size() != 0)
      {
        //        std::cout << "Problem Facet: " << std::endl;
        //        f->Print(std::cout);
        //        std::cout << std::endl;
        //        std::cout << "lines.size(): " << lines.size() << std::endl;
        //        for(unsigned k=0; k< lines.size(); k++)
        //        {
        //          std::cout << "#: " << k <<std::endl;
        //          lines[k].first->Print(std::cout);
        //          lines[k].second->Print(std::cout);
        //          std::cout << std::endl;
        //        }
        //
        // Q: What line in the facetlines is not connected?
        int numberoflines = 0;
        std::vector<int> facetlineindex;
        for (unsigned k = 0; k < lines.size(); k++)
        {
          for (unsigned l = 0; l < facetlines.size(); l++)
          {
#if (0)
            if (lines[k] == facetlines[l])
            {
              numberoflines++;
              facetlineindex.push_back(l);
            }
#else
            dserror(
                "not supported when using std::set in definition of point_line_set. Adapt this. "
                "Anyway this routine is not bugfree! This about before using this function!!!");
#endif
          }
        }
        std::cout << "numberoflines: " << numberoflines << std::endl;

        // test that no line is the same.
        for (unsigned k = 0; k < facetlineindex.size(); k++)
        {
          for (unsigned l = 0; l < facetlineindex.size(); l++)
            if (facetlineindex[k] == facetlineindex[l] and k != l)
              throw std::runtime_error("volume cut facets not closed!!");
        }

        //        //Find the connection.
        //        for(unsigned k=0; k< numberoflines; k++)
        //        {
        //          int index1 = facetlines[facetlineindex[k]].first->Id();
        //          int index2 = facetlines[facetlineindex[k]].second->Id();
        //
        //
        //
        //
        //        }


        // One would have to check if the facetline(s) found here is/are connected in lines. Would
        // probably be best implemented
        // with some sort of tree structure.

        throw std::runtime_error("volume cut facets not closed");
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
                Write the bounding lines of volumecell details for visualization
                Gausspoints of moment fitting are not included
*--------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DumpGmsh(std::ofstream& file)
{
  const plain_facet_set& facete = Facets();

  file << "View \"Volume Cell \" {\n";
  for (plain_facet_set::const_iterator j = facete.begin(); j != facete.end(); j++)
  {
    Facet* ref = *j;
#ifndef OUTPUT_GLOBAL_DIVERGENCE_CELLS
    std::vector<std::vector<double>> corners;
    ref->CornerPointsLocal(ParentElement(), corners);
#else
    const std::vector<std::vector<double>> corners = ref->CornerPointsGlobal(ParentElement());
#endif
    for (unsigned i = 0; i < corners.size(); i++)
    {
      const std::vector<double> coords1 = corners[i];
      const std::vector<double> coords2 = corners[(i + 1) % corners.size()];
      file << "SL(" << coords1[0] << "," << coords1[1] << "," << coords1[2] << "," << coords2[0]
           << "," << coords2[1] << "," << coords2[2] << ")"
           << "{0,0};\n";
    }
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,255} };\n";  // Changing color to red
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                     // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";                 // Disable legend
  file << "View[PostProcessing.NbViews-1].LineWidth = 3.0;";               // increase line width
}

/*---------------------------------------------------------------------------------------*
 * Write the geometry of the volumecell based on facet surfaces            sudhakar 07/15
 * Can be used to check if the geometry of vc is correct or not
 *---------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DumpGmshSolid(std::ofstream& file, Mesh& mesh)
{
  const plain_facet_set& facete = Facets();
  file << "View \"Volume Cell \" {\n";
  for (plain_facet_set::const_iterator j = facete.begin(); j != facete.end(); j++)
  {
    Facet* fac = *j;
    std::vector<Point*> corners = fac->CornerPoints();

    if (corners.size() == 3)
    {
      file << "ST(";
      for (unsigned ipt = 0; ipt < corners.size(); ipt++)
      {
        Point* pt = corners[ipt];
        const double* x = pt->X();
        file << x[0] << "," << x[1] << "," << x[2];
        if (ipt != (corners.size() - 1)) file << ",";
      }
      file << "){0.0,0.0,0.0};\n";
    }
    else
    {
      if (!fac->IsTriangulated()) fac->DoTriangulation(mesh, corners);
      const std::vector<std::vector<Point*>>& triangulation = fac->Triangulation();

      for (std::vector<std::vector<Point*>>::const_iterator j = triangulation.begin();
           j != triangulation.end(); ++j)
      {
        std::vector<Point*> tri = *j;

        if (tri.size() == 3)
          file << "ST(";
        else
          file << "SQ(";

        for (unsigned ipt = 0; ipt < tri.size(); ipt++)
        {
          Point* pt = tri[ipt];
          const double* x = pt->X();
          file << x[0] << "," << x[1] << "," << x[2];
          if (ipt != (tri.size() - 1)) file << ",";
        }
        if (tri.size() == 3)
          file << "){0.0,0.0,0.0};\n";
        else
          file << "){0.0,0.0,0.0,0.0};\n";
      }
    }
  }
  file << "};\n";
}

/*--------------------------------------------------------------------------------------------------------*
        write the boundaries of volumecell and the positions of Gauss points when using
        Moment fitting for visualization a separate file with "Mom_volcell" prefix is generated
        for every volumecell as the gausspoint distribution can be clearly seen
*---------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DumpGmshGaussPointsMomFit(
    const std::vector<std::vector<double>>& gauspts)
{
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "Mom_volcell" << sideno << ".pos";
  std::ofstream file(str.str().c_str());

  DumpGmsh(file);

  file << "View \" Gauss Points \" {\n";
  for (unsigned i = 0; i < gauspts.size(); i++)
  {
    file << "SP(" << gauspts[i][0] << "," << gauspts[i][1] << "," << gauspts[i][2] << ",1){0.0};\n";
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,0} };\n";  // Changing color to black
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                   // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";               // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 4.0;";             // fix point size
  file.close();
}

/*--------------------------------------------------------------------------------------------------------*
        write the boundaries of volumecell and the positions of Gauss points for visualization
        a separate file when using tessellation
*---------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DumpGmshGaussPointsTessellation()
{
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "MTes_volcell" << sideno << ".pos";
  std::ofstream file(str.str().c_str());

  DumpGmsh(file);

  Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
      Teuchos::rcp(new DRT::UTILS::GaussPointsComposite(0));

  const plain_integrationcell_set& cells = IntegrationCells();
  for (plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    GEO::CUT::IntegrationCell* ic = *i;
    switch (ic->Shape())
    {
      case DRT::Element::hex8:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>(ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::tet4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>(ic);
        gpc->Append(gp);
        break;
      }
      default:
      {
        dserror("Include this element here");
        break;
      }
    }
  }

  DRT::UTILS::GaussIntegration gpv(gpc);

  file << "View \" Gauss Points \" {\n";
  for (DRT::UTILS::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    const LINALG::Matrix<3, 1> eta(iquad.Point());
    file << "SP(" << eta(0, 0) << "," << eta(1, 0) << "," << eta(2, 0) << ",1){0.0};\n";
  }

  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,0} };\n";  // Changing color to black
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                   // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";               // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 4.0;";             // fix point size
  file.close();
}

/*----------------------------------------------------------------------------------------------------------*
 * Perform integration of a pre-defined function over this vc using gauss points          Sudhakar
 *01/13 generated from tessellation
 *----------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::integrateSpecificFunctionsTessellation()
{
  Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
      Teuchos::rcp(new DRT::UTILS::GaussPointsComposite(0));

  const plain_integrationcell_set& cells = IntegrationCells();
  for (plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    GEO::CUT::IntegrationCell* ic = *i;
    switch (ic->Shape())
    {
      case DRT::Element::hex8:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>(ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::tet4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>(ic);
        gpc->Append(gp);
        break;
      }
      default:
      {
        dserror("Include this element here");
        break;
      }
    }
  }

  DRT::UTILS::GaussIntegration gpv(gpc);

  double intVal = 0.0;
  for (DRT::UTILS::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    double weight = iquad.Weight();

    const LINALG::Matrix<3, 1> eta(iquad.Point());
    double xx = eta(0, 0);
    double yy = eta(1, 0);
    double zz = eta(2, 0);

    intVal +=
        (pow(xx, 6) + xx * pow(yy, 4) * zz + xx * xx * yy * yy * zz * zz + pow(zz, 6)) * weight;
  }
  std::cout << std::setprecision(20) << "TESSELLATION Integration = " << intVal << "\n";
}

template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::VolumeCell::CreateProjected(
    GEO::CUT::IntegrationCell* ic)
{
  const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<3, nen> xie;

  const std::vector<GEO::CUT::Point*>& cpoints = ic->Points();
  if (cpoints.size() != nen) throw std::runtime_error("non-matching number of points");

  for (unsigned i = 0; i < nen; ++i)
  {
    GEO::CUT::Point* p = cpoints[i];
    LINALG::Matrix<3, 1> xg, xi;
    p->Coordinates(xg.A());
    element_->LocalCoordinates(xg, xi);
    std::copy(xi.A(), xi.A() + 3, &xie(0, i));
  }

  Teuchos::RCP<DRT::UTILS::GaussPoints> gp = DRT::UTILS::GaussIntegration::CreateProjected<distype>(
      xie, ic->CubatureDegree(element_->Shape()));
  return gp;
}

/*------------------------------------------------------------------------------------------------------*
    convert the Gaussian points and weights into appropriate Gauss rule as per BACI implementation
*-------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::VolumeCell::GaussPointsFitting()
{
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp =
      Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints(0));

  for (unsigned i = 0; i < gausPts_.size(); i++)
  {
    LINALG::Matrix<3, 1> xe, xei;
    xe(0, 0) = gausPts_[i][0];
    xe(1, 0) = gausPts_[i][1];
    xe(2, 0) = gausPts_[i][2];

    cgp->Append(xe, weights_(i));
  }

  return cgp;
}

/*--------------------------------------------------------------------------------------------*
                 Generate boundary cells for the cut facets of the volumecell
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::GenerateBoundaryCells(Mesh& mesh,
    const GEO::CUT::Point::PointPosition posi, Element* elem, int BaseNos,
    INPAR::CUT::BCellGaussPts BCellgausstype)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::VolumeCell::GenerateBoundaryCells" );


  // TODO: we have to restructure the creation of boundary cells.
  // bcs should not be stored for a volumecell but for the respective facet, see comments in
  // f->GetBoundaryCells

  const plain_facet_set& facete = Facets();
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fac = *i;

    if (fac->OnBoundaryCellSide() == false)  // we need boundary cells only for the cut facets
      continue;

    // For LevelSetSides generate Boundary Cells in own way.
    if (fac->BelongsToLevelSetSide())
    {
      GenerateBoundaryCellsLevelSetSide(mesh, posi, elem, fac, BaseNos, BCellgausstype);
      continue;
    }

    //--------------------------------------------------------------------
    // Normal vector from parent side is used to identify whether normals
    // from facet is in appropriate direction or not
    //--------------------------------------------------------------------
    const Side* parside = fac->ParentSide();
    const std::vector<Node*>& par_nodes = parside->Nodes();

    std::vector<double> eqnpar(4), eqnfac(4);
    bool rever = false;

#if 1
    std::vector<Point*> pts_par(par_nodes.size());
    for (unsigned parnode = 0; parnode < par_nodes.size(); parnode++)
      pts_par[parnode] = par_nodes[parnode]->point();
    eqnpar = KERNEL::EqnPlaneOfPolygon(pts_par);

    std::vector<Point*> corners = fac->CornerPoints();
    std::vector<Point*> cornersTemp(corners);

    // when finding eqn of plane for the facet, inline points should not be taken
    CUT::KERNEL::DeleteInlinePts(cornersTemp);

    if (cornersTemp.size() != 0)
    {
      eqnfac = KERNEL::EqnPlaneOfPolygon(corners);
      rever = ToReverse(posi, eqnpar, eqnfac);
    }
#endif

#if 0
    std::vector<Point*> parpts(3);

    parpts[0] = par_nodes[0]->point();
    parpts[1] = par_nodes[1]->point();
    parpts[2] = par_nodes[2]->point();

    std::vector<double> eqnpar(4),eqnfac(4);
    // equation of plane denotes normal direction
    eqnpar = KERNEL::EqnPlane( parpts[0], parpts[1], parpts[2] );

    std::vector<Point*> corners = fac->CornerPoints();
    std::vector<Point*> cornersTemp (corners);

    // when finding eqn of plane for the facet, inline points should not be taken
   CUT::KERNEL::DeleteInlinePts( cornersTemp );

   bool rever = false;
   if( cornersTemp.size()!=0 )
   {
     eqnfac = KERNEL::EqnPlanePolygon( cornersTemp );
     rever = ToReverse( posi, eqnpar, eqnfac );
   }
#endif

    // For Marked sides the boundary-cells on outside vc's need to be the same as for inside.
    if (fac->OnMarkedBackgroundSide() and posi == GEO::CUT::Point::outside) rever = !rever;

    if (rever)  // normal from facet is in wrong direction
    {
      std::reverse(corners.begin(), corners.end());  // change ordering to correct this
      std::reverse(cornersTemp.begin(),
          cornersTemp.end());  // Todo: Ager ... I don't understand the sense of this cornersTemp???
    }

    // if no of corners are 3 or 4, just add them as boundary integrationcells directly
    if (corners.size() == 3)
    {
      double areaCell = GEO::CUT::KERNEL::getAreaTri(corners);
      if (areaCell <
          REF_AREA_BCELL)  // What is this Ref_AREA_CELL TOLERANCE?! Make it viable for GLOBAL!!!
        continue;
      NewTri3Cell(mesh, fac, corners);
    }
    else
    {
      if (BCellgausstype == INPAR::CUT::BCellGaussPts_Tessellation)  // generate boundarycell
                                                                     // gausspoints by triangulation
      {
#if 1  // Call split facets -- Remember: the facet can be  concave shaped
        if (!fac->IsTriangulated()) fac->DoTriangulation(mesh, corners);
        const std::vector<std::vector<Point*>>& triangulation = fac->Triangulation();
#else  // creates both tri and quad. less no of Gauss points - but deleted some points leads to
       // error

        if (cornersTemp.size() == 0) continue;
        // Use "cornersTemp" for triangulation - deleted some points results in error
        if (!fac->IsFacetSplit()) fac->SplitFacet(cornersTemp);

        const std::vector<std::vector<Point*>> triangulation = fac->GetSplitCells();
#endif

        for (std::vector<std::vector<Point*>>::const_iterator j = triangulation.begin();
             j != triangulation.end(); ++j)
        {
          std::vector<Point*> tri = *j;

          if (tri.size() == 3)
          {
            double areaCell = GEO::CUT::KERNEL::getAreaTri(tri);
            if (areaCell < REF_AREA_BCELL) continue;
            NewTri3Cell(mesh, fac, tri);
          }
          else if (tri.size() == 4)
          {
            double areaCell = GEO::CUT::KERNEL::getAreaConvexQuad(tri);
            if (areaCell < REF_AREA_BCELL) continue;
            NewQuad4Cell(mesh, fac, tri);
          }
          else
            dserror("Triangulation created neither tri3 or quad4");
        }
      }

      else if (BCellgausstype ==
               INPAR::CUT::BCellGaussPts_MomentFitting)  // generate boundarycell gausspoints by
                                                         // solving moment fitting equations
      {
        BoundarycellIntegration bcell_inte(elem, fac, posi, BaseNos);
        Bcellweights_ = bcell_inte.GenerateBoundaryCellIntegrationRule();
        BcellgausPts_ = bcell_inte.getBcellGaussPointLocation();

        // the boundarycell integration is carriedout in the local coord of the element
        // to project the coordinates of Gauss points, shape functions of element can be used
        //
        //                                            area of facet in global coordinates
        // but to transform the weight, the jacobian = -----------------------------------
        //                                            area of facet in local coordinates
        FacetIntegration bcellLocal(fac, elem, posi, true, false);
        bcellLocal.set_integ_number(1);
        double areaLocal = bcellLocal.integrate_facet();

        FacetIntegration bcellGlobal(fac, elem, posi, true, true);
        bcellGlobal.set_integ_number(1);
        double areaGlobal = bcellGlobal.integrate_facet();
        double jaco = areaGlobal / areaLocal;

        int numBcellpts = BcellgausPts_.size();
        Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp =
            Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints(numBcellpts));

        LINALG::Matrix<3, 1> xeLocal, xeGlobal;
        for (unsigned i = 0; i < BcellgausPts_.size(); i++)
        {
          xeLocal(0, 0) = BcellgausPts_[i][0];
          xeLocal(1, 0) = BcellgausPts_[i][1];
          xeLocal(2, 0) = BcellgausPts_[i][2];

          elem->GlobalCoordinates(xeLocal, xeGlobal);

          cgp->Append(xeGlobal, Bcellweights_(i) * jaco);
        }

        LINALG::Matrix<3, 1> normal;
        double normalFac;
        if (rever)
        {
          // std::reverse(corners.begin(),corners.end());
          normalFac = -1.0;
        }
        else
          normalFac = 1.0;

        normalFac =
            normalFac * sqrt(eqnfac[0] * eqnfac[0] + eqnfac[1] * eqnfac[1] + eqnfac[2] * eqnfac[2]);
        for (unsigned i = 0; i < 3; i++) normal(i, 0) = eqnfac[i] / normalFac;

        DRT::UTILS::GaussIntegration gi(cgp);
        NewArbitraryCell(mesh, fac, corners, gi, normal);
      }
    }
  }
}

/*--------------------------------------------------------------------------------------------*
                 Generate boundary cells for facet cut by level set side

               COMMENT:  Might need to rethink BC-creation, as it generates a lot of tris now.
                         Could probably be enough with quads some times?
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::GenerateBoundaryCellsLevelSetSide(Mesh& mesh,
    const GEO::CUT::Point::PointPosition posi, Element* elem, Facet* fac, int BaseNos,
    INPAR::CUT::BCellGaussPts BCellgausstype)
{
  if (not fac->BelongsToLevelSetSide())
    dserror("Why would you call BC-creation for LS-Side without a LS side?");

  if (BCellgausstype == INPAR::CUT::BCellGaussPts_MomentFitting)
    dserror("Not supported for BC-Cell creation for LevelSetSides.");

  // Is the facet split/triangulated and if it consists of 4 corners is it planar.
  //  Then decompose and create Boundary Cells from triangulation/split.
  //  Otherwise, create quad4/tri3 boundary cell.
  bool istriangulated = (fac->IsTriangulated() or fac->IsFacetSplit());
  bool notsimpleshape =
      !(fac->CornerPoints().size() == 4 and (fac->IsPlanar(mesh, fac->CornerPoints())));

  // UNCOMMENT THIS FOR ONLY TRIANGLES ON SURFACE!
  // DO FULL TRIANGULATION OF BC-SURFACE
  //---------------------------------------
  //  if(fac->CornerPoints().size()>3)
  //  {
  ////    std::cout << "fac->CornerPoints().size()>3?!" << std::endl;
  //    if(not fac->IsTriangulated())
  //    {
  ////      std::cout << "not fac->IsTriangulated()" << std::endl;
  //      fac->DoTriangulation(mesh,fac->Points());  //triangulate everything!
  ////      std::cout << "fac->IsTriangulated(): " << fac->IsTriangulated() << std::endl;
  //    }
  //  }
  //  notsimpleshape = true; //Create BC from triangulation solely
  //  istriangulated = true;
  //---------------------------------------

  //  if( (fac->IsTriangulated() or fac->IsFacetSplit()) and (fac->CornerPoints().size() == 4 and
  //  !(fac->IsPlanar(mesh, fac->CornerPoints())))  )
  if (istriangulated and notsimpleshape)
  {
    std::vector<std::vector<Point*>> facet_triang;
    if (fac->IsTriangulated())
      facet_triang = fac->Triangulation();
    else
      facet_triang = fac->GetSplitCells();

    for (std::vector<std::vector<Point*>>::const_iterator j = facet_triang.begin();
         j != facet_triang.end(); ++j)
    {
      std::vector<Point*> tri = *j;
      std::vector<Point*> tri_temp(tri);  // could be quad?

      //    // when finding eqn of plane for the facet, inline points should not be taken
      //    CUT::KERNEL::DeleteInlinePts( cornersTemp );

      std::vector<double> fac_tri_normal(4);
      fac_tri_normal = KERNEL::EqnPlaneOfPolygon(tri_temp);

      LINALG::Matrix<3, 1> ls_coord(tri_temp[1]->X());
      const std::vector<double> fac_ls_normal =
          elem->GetLevelSetGradient(ls_coord);  // fac->GetLevelSetFacetNormal(elem);
      double dotProduct = fac_tri_normal[0] * fac_ls_normal[0] +
                          fac_tri_normal[1] * fac_ls_normal[1] +
                          fac_tri_normal[2] * fac_ls_normal[2];
      if (posi == GEO::CUT::Point::outside)
      {
        if (dotProduct > 0.0)
        {
          std::reverse(tri_temp.begin(), tri_temp.end());
        }
      }
      else if (posi == GEO::CUT::Point::inside)
      {
        if (dotProduct < 0.0)  // ( < ) should be correct solution.
          std::reverse(tri_temp.begin(), tri_temp.end());
      }
      if (tri_temp.size() == 3)
      {
        double areaCell = GEO::CUT::KERNEL::getAreaTri(tri_temp);
        if (areaCell < REF_AREA_BCELL)
        {
          std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
          continue;
        }
        NewTri3Cell(mesh, fac, tri_temp);
      }
      else if (tri_temp.size() == 4)
      {
        double areaCell = GEO::CUT::KERNEL::getAreaConvexQuad(tri_temp);
        if (areaCell < REF_AREA_BCELL)
        {
          std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
          continue;
        }
        NewQuad4Cell(mesh, fac, tri_temp);
      }
      else
        dserror("Triangulation created neither tri3 or quad4");
    }
  }
  else
  {
    // For non-triangulated facets-> i.e. planar quad4 or tri3. (does quad4 exist here?)

    std::vector<Point*> tri = fac->CornerPoints();
    std::vector<Point*> tri_temp(tri);  // could be quad?


    //    // when finding eqn of plane for the facet, inline points should not be taken
    //    CUT::KERNEL::DeleteInlinePts( cornersTemp );

    // Fix normal direction
    std::vector<double> fac_tri_normal(4);
    fac_tri_normal = KERNEL::EqnPlaneOfPolygon(tri_temp);

    LINALG::Matrix<3, 1> ls_coord(tri_temp[1]->X());
    const std::vector<double> fac_ls_normal =
        elem->GetLevelSetGradient(ls_coord);  // fac->GetLevelSetFacetNormal(elem);
    double dotProduct = fac_tri_normal[0] * fac_ls_normal[0] +
                        fac_tri_normal[1] * fac_ls_normal[1] + fac_tri_normal[2] * fac_ls_normal[2];
    if (posi == GEO::CUT::Point::outside)
    {
      if (dotProduct > 0.0)
      {
        std::reverse(tri_temp.begin(), tri_temp.end());
      }
    }
    else if (posi == GEO::CUT::Point::inside)
    {
      if (dotProduct < 0.0)  // ( < ) should be correct solution.
        std::reverse(tri_temp.begin(), tri_temp.end());
    }

    // Add boundary cell
    if (tri_temp.size() == 3)
    {
      double areaCell = GEO::CUT::KERNEL::getAreaTri(tri_temp);
      if (areaCell < REF_AREA_BCELL)
      {
        std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
      }
      NewTri3Cell(mesh, fac, tri_temp);
    }
    else if (tri_temp.size() == 4)
    {
      double areaCell = GEO::CUT::KERNEL::getAreaConvexQuad(tri_temp);
      if (areaCell < REF_AREA_BCELL)
      {
        std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
      }
      NewQuad4Cell(mesh, fac, tri_temp);
    }
    else
      dserror("Triangulation created neither tri3 or quad4");
  }
}

/*--------------------------------------------------------------------------------------------------------*
    This is to check whether the corner points of the cut side facet is aligned to give outward
normal
*---------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::VolumeCell::ToReverse(const GEO::CUT::Point::PointPosition posi,
    const std::vector<double>& parEqn, const std::vector<double>& facetEqn)
{
  bool rever = false;

  // position is inside
  if (posi == GEO::CUT::Point::outside)  //-3 before...
  {
    if (fabs(parEqn[0]) > TOL_EQN_PLANE && parEqn[0] * facetEqn[0] > 0.0)
      rever = true;
    else if (fabs(parEqn[1]) > TOL_EQN_PLANE && parEqn[1] * facetEqn[1] > 0.0)
      rever = true;
    else if (fabs(parEqn[2]) > TOL_EQN_PLANE && parEqn[2] * facetEqn[2] > 0.0)
      rever = true;
    else
      rever = false;
  }

  // position is outside
  else if (posi == GEO::CUT::Point::inside)  //-2 before...
  {
    if (fabs(parEqn[0]) > TOL_EQN_PLANE && parEqn[0] * facetEqn[0] < 0.0)
      rever = true;
    else if (fabs(parEqn[1]) > TOL_EQN_PLANE && parEqn[1] * facetEqn[1] < 0.0)
      rever = true;
    else if (fabs(parEqn[2]) > TOL_EQN_PLANE && parEqn[2] * facetEqn[2] < 0.0)
      rever = true;
    else
      rever = false;
  }
  return rever;
}

/*------------------------------------------------------------------------------------------*
   When DirectDivergence method is used for gauss point generation, for every gauss point
   on the facet, an internal gauss rule is to be generated to find the modified integrand
*-------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::VolumeCell::GenerateInternalGaussRule(
    Teuchos::RCP<DRT::UTILS::GaussPoints>& gp)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::VolumeCell::GenerateInternalGaussRule" );


  DRT::UTILS::GaussIntegration grule(gp);

  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp =
      Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints(0));

  for (DRT::UTILS::GaussIntegration::iterator quadint = grule.begin(); quadint != grule.end();
       ++quadint)
  {
    const LINALG::Matrix<3, 1> etaFacet(
        quadint.Point());  // coordinates and weight of main gauss point
    LINALG::Matrix<3, 1> intpt(etaFacet);

    DRT::UTILS::GaussIntegration gi(
        DRT::Element::line2, (DIRECTDIV_GAUSSRULE - 1));  // internal gauss rule for interval (-1,1)

    // x-coordinate of main Gauss point is projected in the reference plane
    double xbegin =
        (RefEqnPlane_[3] - RefEqnPlane_[1] * etaFacet(1, 0) - RefEqnPlane_[2] * etaFacet(2, 0)) /
        RefEqnPlane_[0];

    double jac = fabs(xbegin - etaFacet(0, 0)) * 0.5;  // jacobian for 1D transformation rule

    // -----------------------------------------------------------------------------
    // project internal gauss point from interval (-1,1) to the actual interval
    // -----------------------------------------------------------------------------
    for (DRT::UTILS::GaussIntegration::iterator iqu = gi.begin(); iqu != gi.end(); ++iqu)
    {
      const LINALG::Matrix<1, 1> eta(iqu.Point());
      double weight = iqu.Weight();

      double xmid = 0.5 * (xbegin + etaFacet(0, 0));
      intpt(0, 0) = (xmid - xbegin) * eta(0, 0) + xmid;  // location of internal gauss points

      weight = weight * jac;  // weight of internal gauss points
      if (xbegin > etaFacet(0, 0)) weight = -1.0 * weight;

      weight =
          weight * quadint.Weight();  // multiply with weight of main Gauss points so that internal
                                      // and main pts can be combined into a single data structure
      cgp->Append(intpt, weight);
    }
  }
  return cgp;
}

/*------------------------------------------------------------------------------------------*
   Moment fitting equations are solved at each volume cell to construct integration rules
*-------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::MomentFitGaussWeights(
    Element* elem, Mesh& mesh, bool include_inner, INPAR::CUT::BCellGaussPts BCellgausstype)
{
#ifdef LOCAL
  // position is used to decide whether the ordering of points are in clockwise or not
  if (Position() == Point::undecided)
  {
    if (!SetPositionCutSideBased())
    {
      dserror("undefined position for the volumecell");
    }
  }

  // if the volumecell is inside and include_inner is false, no need to compute the Gaussian points
  // as this vc will never be computed in xfem algorithm
  if (Position() == Point::inside && include_inner == false) return;

  int BaseNos = 84;  // number of base functions to be used in the integration
  VolumeIntegration vc_inte(this, elem, Position(), BaseNos);

  weights_ = vc_inte.compute_weights();        // obtain the integration weight at all points
  gausPts_ = vc_inte.getGaussPointLocation();  // get the coordinates of all the Gauss points

  gp_ = GaussPointsFitting();  // convert the weight and the location to Gauss rule

  // generate boundary cells -- when using tessellation this is automatically done
  GenerateBoundaryCells(mesh, Position(), elem, BaseNos, BCellgausstype);

  // std::cout<<"MOMENT FITTING ::: Number of points = "<<weights_.Length()<<"\n";
#else

  // std::cout << "DirectDivergence Is Used!!!! MomFitting not functional in Global coordinates. "
  // << std::endl;

  std::cout << "MomFitting not functional in Global coordinates. NOTHING IS DONE!" << std::endl;

//  DirectDivergenceGaussRule( elem, mesh, include_inner, BCellgausstype );
#endif
}

/*---------------------------------------------------------------------------------------------------------------*
                     The facets that have non-zero x-component normal is triangulated. sudhakar
03/12 The gauss integration rules are generated by applying divergence theorem The reference facet
is identified which will be used to find the modified integral in fluid integration
*----------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DirectDivergenceGaussRule(
    Element* elem, Mesh& mesh, bool include_inner, INPAR::CUT::BCellGaussPts BCellgausstype)
{
  if (elem->Shape() != DRT::Element::hex8 && elem->Shape() != DRT::Element::hex20)
    throw std::runtime_error("DirectDivergenceGaussRule: Just hex8 and hex20 avaiable yet in DD!");

  if (BCellgausstype != INPAR::CUT::BCellGaussPts_Tessellation)
    dserror(
        "DirectDivergenceGaussRule: just INPAR::CUT::BCellGaussPts_Tessellation supported at the "
        "moment!");

  // position is used to decide whether the ordering of points are in clockwise or not
  if (Position() == Point::undecided)
  {
    if (!SetPositionCutSideBased())
    {
      dserror("undefined position for the volumecell");
    }
  }

  // if the volumecell is inside and includeinner is false, no need to compute the Gaussian points
  // as this vc will never be computed in xfem algorithm
  if (Position() == Point::inside and include_inner == false) return;

  // If the Volume Cell consists of less than 4 facets, it can't span a volume in 3D.
  if (Facets().size() < 4)
    dserror("If the Volume Cell consists of less than 4 facets, it can't span a volume in 3D?");
  // return;

  isNegligibleSmall_ = false;


  DirectDivergence dd(this, elem, Position(), mesh);

  RefEqnPlane_.reserve(4);  // it has to store a,b,c,d in ax+by+cz=d

  Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
      dd.VCIntegrationRule(RefEqnPlane_);  // compute main gauss points
  gp_ = GenerateInternalGaussRule(gp);  // compute internal gauss points for every main gauss point

#ifdef DEBUGCUTLIBRARY  // write volumecell, main and internal Gauss points
  {
    DRT::UTILS::GaussIntegration gpi(gp_);
    dd.DivengenceCellsGMSH(gpi, gp);
  }
#endif

  // compute volume of this cell
  // also check whether generated gauss rule predicts volume accurately
  // also check this vc can be eliminated due to its very small volume
  bool isNegVol = false;
  {
    DRT::UTILS::GaussIntegration gpi(gp_);
    dd.DebugVolume(gpi, isNegVol);

    // then this vol is extremely small that we erase the gauss points
    if (isNegVol)
    {
      gp_.reset();
      isNegligibleSmall_ = true;
    }
  }

  if (not isNegVol)
  {
#ifdef LOCAL

#else
    // we have generated the integration rule in global coordinates of the element
    // Now we map this rule to local coodinates since the weak form evaluation is done on local
    // coord
    ProjectGaussPointsToLocalCoodinates();
#endif

#if 0  // integrate a predefined function
    DRT::UTILS::GaussIntegration gpi(gp_);
    dd.IntegrateSpecificFuntions( gpi );
#endif
  }
  // generate boundary cells -- when using tessellation this is automatically done
  GenerateBoundaryCells(mesh, Position(), elem, 0, BCellgausstype);
}

/*----------------------------------------------------------------------------------------------------*
 * Project the integration rule generated on global coordinate system of sudhakar 05/15 the
 *background element to its local coordinates
 *----------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::ProjectGaussPointsToLocalCoodinates()
{
  if (element_->Shape() != DRT::Element::hex8)
    dserror("Currently Direct divergence in global coordinates works only for hex8 elements\n");

  DRT::UTILS::GaussIntegration intpoints(gp_);

  if (element_->isShadow() && (element_->getQuadShape() == DRT::Element::hex20 ||
                                  element_->getQuadShape() == DRT::Element::hex27))
  {
    switch (element_->getQuadShape())
    {
      case DRT::Element::hex20:
      {
        LINALG::Matrix<3, 20> xyze;
        element_->CoordinatesQuad(xyze.A());
        gp_ = DRT::UTILS::GaussIntegration::ProjectGaussPointsGlobalToLocal<DRT::Element::hex20>(
            xyze, intpoints, false);
        break;
      }
      case DRT::Element::hex27:
      {
        LINALG::Matrix<3, 27> xyze;
        element_->CoordinatesQuad(xyze.A());
        gp_ = DRT::UTILS::GaussIntegration::ProjectGaussPointsGlobalToLocal<DRT::Element::hex27>(
            xyze, intpoints, false);
        break;
      }
      default:
      {
        dserror(
            "Currently Direct divergence in global coordinates works only for hex8, hex20, hex27 "
            "elements\n");
        break;
      }
    }
  }
  else
  {
    LINALG::Matrix<3, 8> xyze;
    element_->Coordinates(xyze.A());
    gp_ = DRT::UTILS::GaussIntegration::ProjectGaussPointsGlobalToLocal<DRT::Element::hex8>(
        xyze, intpoints, false);
  }
}

/*-------------------------------------------------------------------------------------*
| Return Ids of all the points associated with this volumecell           shahmiri 06/12
*--------------------------------------------------------------------------------------*/
const std::set<int>& GEO::CUT::VolumeCell::VolumeCellPointIds()
{
  if (vcpoints_ids_.size() != 0)
  {
    return vcpoints_ids_;
  }
  else
  {
    const plain_facet_set& facete = Facets();

    // loop over facets
    for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
    {
      Facet* fe = *i;
      const std::vector<Point*>& corners = fe->CornerPoints();

      for (std::vector<Point*>::const_iterator c = corners.begin(); c != corners.end(); c++)
      {
        Point* pt = *c;
        vcpoints_ids_.insert(pt->Id());
      }
    }
  }

  if (vcpoints_ids_.size() == 0) dserror("The size of volumecell points is zero!!");

  return vcpoints_ids_;
}

/*-------------------------------------------------------------------------------------*
| Find Position of the Volumecell based on the orientation of the cut_sides   ager 08/15
*--------------------------------------------------------------------------------------*/
bool GEO::CUT::VolumeCell::SetPositionCutSideBased()
{
  if (Position() != Point::undecided)
    dserror(
        "Do not call FindPositionCutSideBased() if Position for volumecell is already set (%d)!",
        (int)Position());

  std::map<Facet*, bool> outsidenormal;

  const plain_facet_set& facets = Facets();

  // First find a facet based on a background side!
  Facet* f;
  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    f = *i;
    if (!f->OnCutSide())
    {
      //-----
      // STEP 1: Get centroid of the parent element
      //-----
      LINALG::Matrix<3, 1> elecen;
      ParentElement()->ElementCenter(elecen);

      //-----
      // STEP 2: Get geometric center point of the facet
      // For concave facets, this is not the actual center, but it does not matter
      //-----
      LINALG::Matrix<3, 1> facecen;
      for (std::vector<Point*>::const_iterator fit = f->CornerPoints().begin();
           fit != f->CornerPoints().end(); fit++)
      {
        for (unsigned dim = 0; dim < 3; dim++) facecen(dim, 0) += (*fit)->X()[dim];
      }

      for (unsigned dim = 0; dim < 3; dim++)
        facecen(dim, 0) = facecen(dim, 0) / f->CornerPoints().size();

      //-----
      // STEP 3: Construct a vector that points FROM element centre TO the facet centre
      // This reference vector points outside the background element! (for LOCAL this is a unit
      // normal vector)
      //-----
      LINALG::Matrix<3, 1> ref_vec;
      ref_vec.Update(1.0, facecen, -1.0, elecen);

      //-----
      // STEP 4: get point from parent side as the normal orientation should be calculated from the
      // side!!!
      // get unit normal vector
      std::vector<double> eqn_plane = KERNEL::EqnPlaneOfPolygon(f->CornerPoints());


      //-----
      // STEP 5: Take dot product with the normal of facet
      // If both are in the opposite direction, then the facet nodes are arranged clockwise
      //-----
      LINALG::Matrix<3, 1> norm_fac;
      for (unsigned dim = 0; dim < 3; dim++) norm_fac(dim, 0) = eqn_plane[dim];

      double dotProduct =
          ref_vec.Dot(norm_fac);  // > 0 ...Side Normal is pointing outside of the Element //< 0
                                  // Side Normal points inside the element!
      if (abs(dotProduct) < BASICTOL)
      {
        std::cout << "Reference Vector: (" << ref_vec(0, 0) << "," << ref_vec(1, 0) << ","
                  << ref_vec(2, 0) << ")" << std::endl;
        std::cout << "Facet Vector: (" << norm_fac(0, 0) << "," << norm_fac(1, 0) << ","
                  << norm_fac(2, 0) << ")" << std::endl;
        dserror("Check this really small dotProduct! %d", dotProduct);
      }

      outsidenormal[f] = (dotProduct > 0);
    }  // if not cutfacet
  }    // for facets

  int iter = 0;
  bool done = false;
  while (outsidenormal.size() != facets.size() && iter < 1000)
  {
    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* ff = *i;
      if (ff->OnCutSide() && outsidenormal.find(ff) == outsidenormal.end())  // cutside
      {
        for (std::map<Facet*, bool>::iterator on = outsidenormal.begin(); on != outsidenormal.end();
             ++on)
        {
          bool consistant_normal = false;
          if (ff->HaveConsistantNormal(on->first, consistant_normal))
          {
            if (consistant_normal)
              outsidenormal[ff] = outsidenormal[on->first];
            else
              outsidenormal[ff] = !outsidenormal[on->first];
            done = true;  // If we can a least identify the direction on one CutSide facet this is
                          // fine in principle
            break;
          }
        }
      }
    }
    iter++;
  }

  if (iter == 1000 && !done)
  {
    throw std::runtime_error(
        "SetPositionCutSideBased failed: too many iterations (theoretically a facet with many "
        "points could also lead to this)!");
    return false;
  }

  Point::PointPosition posi = Point::undecided;
  for (std::map<Facet*, bool>::iterator on = outsidenormal.begin(); on != outsidenormal.end(); ++on)
  {
    if (on->first->OnCutSide())  // cutside
    {
      double prod = 0.0;
      std::vector<double> eqn_plane_facet = KERNEL::EqnPlaneOfPolygon(on->first->CornerPoints());
      std::vector<Point*> spoints(on->first->ParentSide()->Nodes().size());

      for (unsigned int i = 0; i < on->first->ParentSide()->Nodes().size(); ++i)
        spoints[i] = on->first->ParentSide()->Nodes()[i]->point();

      // get unit normal vector
      std::vector<double> eqn_plane_side = KERNEL::EqnPlaneOfPolygon(spoints);

      for (unsigned int dim = 0; dim < 3; ++dim) prod += eqn_plane_facet[dim] * eqn_plane_side[dim];

      if ((on->second && prod > 0) || (!on->second && prod < 0))  // this means that the
      {
        if (posi != Point::undecided && posi != Point::inside)
          throw std::runtime_error(
              "SetPositionCutSideBased: posi != Point::undecided && posi != Point::inside (Are all "
              "you Cut Sides oriented correct?)");
        // dserror("SetPositionCutSideBased: posi != Point::undecided && posi != Point::inside (Are
        // all you Cut Sides oriented correct?)");
        posi = Point::inside;
      }
      else
      {
        if (posi != Point::undecided && posi != Point::outside)
          throw std::runtime_error(
              "SetPositionCutSideBased: posi != Point::undecided && posi != Point::outside (Are "
              "all you Cut Sides oriented correct?)");
        // dserror("SetPositionCutSideBased: posi != Point::undecided && posi != Point::outside (Are
        // all you Cut Sides oriented correct?)");
        posi = Point::outside;
      }
    }
  }

  if (posi != Point::undecided)
  {
    position_ = posi;
    return true;
  }
  else
  {
    return false;
  }
}
