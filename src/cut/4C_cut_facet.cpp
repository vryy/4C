/*---------------------------------------------------------------------*/
/*! \file

\brief cut facet (surface descripted by a cycle of points)

\level 3


*----------------------------------------------------------------------*/
#include "4C_cut_boundarycell.hpp"
#include "4C_cut_kernel.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_triangulateFacet.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_linalg_gauss.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::Facet::Facet(
    Mesh& mesh, const std::vector<Point*>& points, Side* side, bool cutsurface)
    : points_(points),
      parentside_(side),
      planar_(false),
      planar_known_(false),
      position_(cutsurface ? Point::oncutsurface : Point::undecided),
      is_planar_computed_(false)
{
  find_corner_points();

  if (cutsurface)
  {
    for (std::vector<Point*>::const_iterator i = points.begin(); i != points.end(); ++i)
    {
      Point* p = *i;
      p->Position(Point::oncutsurface);
    }
  }
  else
  {
    // On multiple cuts there are facets on element sides that belong to an
    // old cut surface. Thus if all nodes are on a cut surface, the facet is
    // as well.
    bool allonsurface = true;
    for (std::vector<Point*>::const_iterator i = points.begin(); i != points.end(); ++i)
    {
      Point* p = *i;
      if (p->Position() != Point::oncutsurface)
      {
        allonsurface = false;
        break;
      }
    }
    if (allonsurface)
    {
      // If my side has an id this facet is actually on a cut
      // surface. Otherwise it could be an inside or outside facet. The actual
      // decision does not matter much.
      if (OnCutSide())
      {
        position_ = Point::oncutsurface;
      }
      else
      {
        position_ = Point::undecided;
      }
    }
  }
  for (std::vector<Point*>::const_iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    p->Register(this);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::Register(VolumeCell* cell)
{
  cells_.insert(cell);
  if (cells_.size() > 2)
  {
    this->print();
    for (plain_volumecell_set::const_iterator ic = cells_.begin(); ic != cells_.end(); ++ic)
    {
      std::cout << "\n\nVolumeCell " << *ic << std::endl;
      (*ic)->print(std::cout);
    }

    // write details of volume cells
    Output::GmshVolumeCellsOnly(cells_);

    std::ostringstream ostr;
    ostr << "Too many volume cells at facet! ( num cells = " << cells_.size() << " )";
    FOUR_C_THROW(ostr.str().c_str());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::DisconnectVolume(VolumeCell* cell) { cells_.erase(cell); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::OnCutSide() const { return parentside_->IsCutSide(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::OnBoundaryCellSide() const { return parentside_->IsBoundaryCellSide(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::on_marked_background_side() const
{
  return parentside_->is_marked_background_side();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Facet::SideId() const { return parentside_->Id(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::Coordinates(double* x)
{
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    p->Coordinates(x);
    x += 3;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::CornerCoordinates(double* x)
{
  find_corner_points();
  for (std::vector<Point*>::const_iterator i = corner_points_.begin(); i != corner_points_.end();
       ++i)
  {
    Point* p = *i;
    p->Coordinates(x);
    x += 3;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::GetAllPoints(Mesh& mesh, PointSet& cut_points, bool dotriangulate)
{
  if (is_planar(mesh, dotriangulate))
  {
    std::copy(points_.begin(), points_.end(), std::inserter(cut_points, cut_points.begin()));
    for (plain_facet_set::iterator i = holes_.begin(); i != holes_.end(); ++i)
    {
      Facet* h = *i;
      h->GetAllPoints(mesh, cut_points);
    }
  }
  else
  {
    for (std::vector<std::vector<Point*>>::iterator i = triangulation_.begin();
         i != triangulation_.end(); ++i)
    {
      std::vector<Point*>& tri = *i;
      std::copy(tri.begin(), tri.end(), std::inserter(cut_points, cut_points.begin()));
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::AddHole(Facet* hole)
{
  double dot = 1.0;
  while (dot > 0)
  {
    std::vector<double> eqn_plane_hole = Kernel::EqnPlaneOfPolygon(hole->points_);
    std::vector<double> eqn_plane_facet = Kernel::EqnPlaneOfPolygon(corner_points_);
    dot = 0.0;
    for (unsigned int dim = 0; dim < 3; ++dim) dot += eqn_plane_hole[dim] * eqn_plane_facet[dim];
    if (dot >
        0)  // if this happens we get problems with the ear clipping algorithm for triangulation
    {
      std::reverse(hole->corner_points_.begin(), hole->corner_points_.end());
      std::reverse(hole->points_.begin(), hole->points_.end());
      std::cout << "WARNING:: Facet::AddHole: points_ and hole->points_ are not alighned correct "
                   "(other direction)! "
                << std::endl;
    }
  }

  for (std::vector<Point*>::iterator i = hole->points_.begin(); i != hole->points_.end(); ++i)
  {
    Point* p = *i;
    p->Register(this);
  }
  holes_.insert(hole);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::is_planar(Mesh& mesh, bool dotriangulate)
{
  if (dotriangulate)
  {
    // if ( not IsTriangulated() and points_.size() > 3 )
    if (not IsTriangulated() and corner_points_.size() > 3)
    {
      create_triangulation(mesh, points_);
      planar_known_ = true;
      planar_ = false;
    }
  }

  if (not planar_known_)
  {
    planar_ = is_planar(mesh, points_);
    planar_known_ = true;
  }

  return planar_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::is_planar(Mesh& mesh, const std::vector<Point*>& points)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::Cut::Facet::is_planar" );

  if (is_planar_computed_) return is_planar_;


  Core::LinAlg::Matrix<3, 1> x1;
  Core::LinAlg::Matrix<3, 1> x2;
  Core::LinAlg::Matrix<3, 1> x3;

  Core::LinAlg::Matrix<3, 1> b1;
  Core::LinAlg::Matrix<3, 1> b2;
  Core::LinAlg::Matrix<3, 1> b3;

  unsigned i = normal(points, x1, x2, x3, b1, b2, b3);
  if (i == 0)  // all on one line is ok
  {
    is_planar_computed_ = true;
    is_planar_ = true;
    return true;
  }

  Core::LinAlg::Matrix<3, 3> A;
  std::copy(b1.data(), b1.data() + 3, A.data());
  std::copy(b2.data(), b2.data() + 3, A.data() + 3);
  std::copy(b3.data(), b3.data() + 3, A.data() + 6);

  // std::copy( points.begin(), points.end(), std::ostream_iterator<Point*>( std::cout, "; " ) );
  // std::cout << "\n";

  for (++i; i < points.size(); ++i)
  {
    Point* p = points[i];
    p->Coordinates(x3.data());

    x3.update(-1, x1, 1);

    Core::LinAlg::Matrix<3, 3> B;
    B = A;
    x2 = 0;
    double det = Core::LinAlg::gaussElimination<true, 3, double>(B, x3, x2);
    if (fabs(det) < LINSOLVETOL)
    {
      FOUR_C_THROW("failed to find point position");
    }

    if (fabs(x2(2)) > PLANARTOL)
    {
      // there is one point that is not within the plain

      create_triangulation(mesh, points);

      is_planar_computed_ = true;
      is_planar_ = false;
      return false;
    }
  }
  is_planar_computed_ = true;
  is_planar_ = true;
  return true;
}

/*----------------------------------------------------------------------------*
 * Only used in tetmesh debugging!!
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::ShareSameCutSide(Facet* f) /* Does the facet share the same cut-side? */
{
  if (this->parentside_->Id() == f->ParentSide()->Id())
  {
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
  Perform triangulation of the facet.                            sudhakar 06/15
  When convex, centre-point-triangulation is used that can better approximate
  the geometry. For concave facets, Ear-clipping procedure is utilized --
  slightly less accurate due to deleting the inline points but no other choice
*-----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::create_triangulation(Mesh& mesh, const std::vector<Point*>& points)
{
  // Perform centre-point triangulation
  // Find the middle point (M) and join them with the corners of the facet to create triangles
  // Works only for convex facets, because of the middle point calculation
  //                                                                  |
  //                                       4 _______ 3                |
  //                                        /.     .\                 |
  //                                       /  .   .  \                |
  //                                      /    . .    \               |
  //                                  5  /      . M    \ 2            |
  //                                     \     . .     /              |
  //                                      \   .   .   /               |
  //                                       \._______./                |
  //                                      0          1                |
  //

  // TODO: We force the facets created from levelset side to use centre-pt-tri
  // This is not completely correct but since we use a single levelset side for each element
  // (even though the levelset can break into two), we get some weird facet shapes.
  // This procedure seems to give better result even if the facet is concave in such cases
  // Also we assume that the levelset facet never contains a hole --> May be in future changes?
  Core::Geo::Cut::FacetShape geoType;
  std::vector<int> concave_ids = Kernel::CheckConvexity(corner_points_, geoType, false, false);
  if (((geoType == Core::Geo::Cut::Convex) or belongs_to_level_set_side()) and (not HasHoles()))
  {
    std::vector<Point*> pts(points);
    // Find the middle point
    Core::LinAlg::Matrix<3, 1> cur;
    Core::LinAlg::Matrix<3, 1> avg(true);  // fill with zeros
    for (std::vector<Point*>::iterator i = pts.begin(); i != pts.end(); i++)
    {
      Point* p1 = *i;
      p1->Coordinates(cur.data());
      avg.update(1.0, cur, 1.0);
    }
    avg.scale(1.0 / pts.size());
    // One could create a better approx for the "midpoint" here parsing information from the
    // levelset
    //  and thus finding the zero level set easier. Might be worthwhile testing.
    Point* p_mid = mesh.NewPoint(avg.data(), nullptr, ParentSide(),
        0.0);  // change tolerance here intelligently !!! - basically
               // there is no reason why I'd like to merge here!
    p_mid->Position(Position());
    p_mid->Register(this);

    // form triangles and store them in triangulation_
    triangulation_.clear();
    std::vector<Point*> newtri(3);
    newtri[0] = p_mid;
    for (int i = 0; i < (int)pts.size(); i++)
    {
      Point* pt1 = pts[i];
      Point* pt2 = pts[(i + 1) % pts.size()];

      // For LevelSet, make sure to orient a concave point on the 2nd entry of the triangulation.
      //  This assumption is necessary for Direct Divergence to work correctly.
      if (!belongs_to_level_set_side())
      {
        newtri[1] = pt1;
        newtri[2] = pt2;
      }
      else
      {
        bool reverse_order = false;
        for (std::vector<int>::iterator j = concave_ids.begin(); j != concave_ids.end(); ++j)
        {
          if ((*j) == i) reverse_order = true;
        }
        if (!reverse_order)
        {
          newtri[1] = pt1;
          newtri[2] = pt2;
        }
        else
        {
          newtri[2] = pt1;
          newtri[1] = pt2;
        }
      }

      triangulation_.push_back(newtri);
    }
  }
  // When the facet is concave, perform triangulation by applying Ear-clipping procedure
  // This is less accurate than the centre-point-triangulation, as it involves deleting in-line
  // points upto a tolerance But we have no choice in the case of concave polyhedra
  else
  {
    triangulation_.clear();
    std::vector<Point*> facetpts(points);

    if (not HasHoles())
    {
      std::vector<int> ptc;
      Core::Geo::Cut::TriangulateFacet tf(facetpts);
      tf.EarClipping(ptc, true, true);
      triangulation_ = tf.GetSplitCells();
    }
    else
    {
      std::vector<std::vector<Point*>> holepts;
      for (plain_facet_set::iterator ifacet = holes_.begin(); ifacet != holes_.end(); ++ifacet)
      {
        holepts.push_back((*ifacet)->corner_points_);
      }
      Core::Geo::Cut::TriangulateFacet tf(facetpts, holepts);
      tf.ear_clipping_with_holes(parentside_);
      triangulation_ = tf.GetSplitCells();
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::get_nodal_ids(
    Mesh& mesh, const std::vector<Point*>& points, std::vector<int>& nids)
{
  for (std::vector<Point*>::const_iterator i = points.begin(); i != points.end(); ++i)
  {
    Point* p = *i;
    Node* n = p->CutNode();
    if (n == nullptr)
    {
      plain_int_set point_id;
      point_id.insert(p->Id());
      n = mesh.GetNode(point_id, p->X());
    }
    nids.push_back(n->Id());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::equals(
    const std::vector<Point*>& my_points, const std::vector<Point*>& facet_points)
{
  if (my_points.size() != facet_points.size() or holes_.size() > 0) return false;

  unsigned size = my_points.size();
  unsigned shift =
      std::find(my_points.begin(), my_points.end(), facet_points.front()) - my_points.begin();

  bool forward_match = true;
  for (unsigned i = 0; i < size; ++i)
  {
    unsigned j = (i + shift) % size;
    if (my_points[j] != facet_points[i])
    {
      forward_match = false;
      break;
    }
  }
  if (not forward_match)
  {
    for (unsigned i = 0; i < size; ++i)
    {
      unsigned j = (shift + size - i) % size;
      if (my_points[j] != facet_points[i])
      {
        return false;
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::IsCutSide(Side* side)
{
  if (parentside_ == side) return false;

  int count = 0;
  for (std::vector<Point*>::iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    if (p->IsCut(side))
    {
      count += 1;
      if (count > 1)
      {
        return true;
      }
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::Position(Point::PointPosition pos)
{
  FOUR_C_ASSERT(IsCutPositionUnchanged(position_, pos),
      "Are you sure that you want to change the facet-position from inside to outside or vice "
      "versa?");

  if (position_ == Point::undecided)
  {
    if (position_ != pos)
    {
      position_ = pos;
      if (pos == Point::outside or pos == Point::inside)
      {
        for (std::vector<Point*>::iterator i = points_.begin(); i != points_.end(); ++i)
        {
          Point* p = *i;
          if (p->Position() == Point::undecided)
          {
            p->Position(pos);
          }
        }
        for (plain_volumecell_set::iterator i = cells_.begin(); i != cells_.end(); ++i)
        {
          VolumeCell* c = *i;
          c->Position(pos);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::GetLines(std::map<std::pair<Point*, Point*>, plain_facet_set>& lines)
{
  GetLines(points_, lines);

  // add hole lines but do not connect with parent facet
  for (plain_facet_set::iterator i = holes_.begin(); i != holes_.end(); ++i)
  {
    Facet* hole = *i;
    hole->GetLines(lines);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::GetLines(
    const std::vector<Point*>& points, std::map<std::pair<Point*, Point*>, plain_facet_set>& lines)
{
  unsigned length = points.size();
  for (unsigned i = 0; i < length; ++i)
  {
    unsigned j = (i + 1) % length;

    Point* p1 = points[i];
    Point* p2 = points[j];

    if (p1->Id() < p2->Id())
    {
      lines[std::make_pair(p1, p2)].insert(this);
    }
    else if (p1->Id() > p2->Id())
    {
      lines[std::make_pair(p2, p1)].insert(this);
    }
    else if (p1->Id() == p2->Id() and length == 1)
    {
      lines[std::make_pair(p1, p2)].insert(this);
    }
    else
      FOUR_C_THROW("line creation with identical begin and end points\n");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::is_line(Point* p1, Point* p2)
{
  if (IsTriangulated())
  {
    for (std::vector<std::vector<Point*>>::iterator i = triangulation_.begin();
         i != triangulation_.end(); ++i)
    {
      std::vector<Point*>& points = *i;
      if (is_line(points, p1, p2)) return true;
    }
  }
  else
  {
    if (is_line(points_, p1, p2)) return true;
    for (plain_facet_set::iterator i = holes_.begin(); i != holes_.end(); ++i)
    {
      Facet* hole = *i;
      if (hole->is_line(p1, p2)) return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::is_line(const std::vector<Point*>& points, Point* p1, Point* p2)
{
  std::vector<Point*>::const_iterator i1 = std::find(points.begin(), points.end(), p1);
  if (i1 != points.end())
  {
    std::vector<Point*>::const_iterator i2 = i1 + 1;
    if (i2 != points.end())
    {
      if (*i2 == p2)
      {
        return true;
      }
    }
    else
    {
      i2 = points.begin();
      if (*i2 == p2)
      {
        return true;
      }
    }
    if (i1 != points.begin())
    {
      i2 = i1 - 1;
      if (*i2 == p2)
      {
        return true;
      }
    }
    else
    {
      i2 = points.end() - 1;
      if (*i2 == p2)
      {
        return true;
      }
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::Contains(Point* p) const
{
  if (IsTriangulated())
  {
    for (std::vector<std::vector<Point*>>::const_iterator i = triangulation_.begin();
         i != triangulation_.end(); ++i)
    {
      const std::vector<Point*>& points = *i;
      if (std::find(points.begin(), points.end(), p) != points.end()) return true;
    }
  }
  else
  {
    if (std::find(points_.begin(), points_.end(), p) != points_.end()) return true;
    for (plain_facet_set::const_iterator i = holes_.begin(); i != holes_.end(); ++i)
    {
      Facet* hole = *i;
      if (hole->Contains(p)) return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::Contains(const std::vector<Point*>& side) const
{
  for (std::vector<Point*>::const_iterator i = side.begin(); i != side.end(); ++i)
  {
    Point* p = *i;
    if (not Contains(p))
    {
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::Contains(const plain_facet_set& vcell) const
{
  for (plain_volumecell_set::const_iterator i_vc = cells_.begin(); i_vc != cells_.end(); ++i_vc)
  {
    const VolumeCell& vc = **i_vc;
    if (vc.IsEqual(vcell)) return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::ContainsSome(const std::vector<Point*>& side) const
{
  for (std::vector<Point*>::const_iterator i = side.begin(); i != side.end(); ++i)
  {
    Point* p = *i;
    if (Contains(p))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::Touches(Facet* f)
{
  for (std::vector<Point*>::iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = *i;
    if (f->Contains(p))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------------------------------*
 * If this Facet has a CommonEdge with another facet, based on this edge the point ordering is
 *checked ager 08/15
 *----------------------------------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::have_consistant_normal(Facet* f, bool& result)
{
  // 1// create map which stores references to all points cycles, which should be matched!
  //(this function uses the fact, that ordering of the points in holes is opposite to the ordering
  // of the outer points cycle!)
  std::map<std::vector<Point*>*, std::vector<Point*>*> matchingcycles;
  matchingcycles[&(corner_points_)] = &(f->corner_points_);  // outer points cycles
  if (HasHoles())  // outer point cycle with inner hole point cycles
  {
    for (plain_facet_set::iterator hole = holes_.begin(); hole != holes_.end(); ++hole)
    {
      matchingcycles[&((*hole)->corner_points_)] = &(f->corner_points_);
    }
  }
  if (f->HasHoles())
  {
    for (plain_facet_set::iterator hole = f->holes_.begin(); hole != f->holes_.end(); ++hole)
    {
      matchingcycles[&(corner_points_)] = &((*hole)->corner_points_);
    }
  }
  if (HasHoles() && f->HasHoles())  // inner hole point cycles
  {
    for (plain_facet_set::iterator hole = holes_.begin(); hole != holes_.end(); ++hole)
    {
      for (plain_facet_set::iterator fhole = f->holes_.begin(); fhole != f->holes_.end(); ++fhole)
      {
        matchingcycles[&((*hole)->corner_points_)] = &((*fhole)->corner_points_);
      }
    }
  }

  for (std::map<std::vector<Point*>*, std::vector<Point*>*>::iterator cyclepair =
           matchingcycles.begin();
       cyclepair != matchingcycles.end(); ++cyclepair)
  {
    std::vector<Point*>& points1 = *(cyclepair->first);
    std::vector<Point*>& points2 = *(cyclepair->second);

    // 2// find the corner points which are matching and store the local idx in a map<int,int>
    std::map<unsigned int, unsigned int> matchingppoints;

    for (unsigned int i = 0; i < points1.size(); i++)
    {
      for (unsigned int j = 0; j < points2.size(); j++)
      {
        if (points1[i] == points2[j])
        {
          matchingppoints[i] = j;
        }
      }
    }

    // 3// if there are more than 2 matching points, remove points at the end of the common edge(s)!
    if (matchingppoints.size() > 2)
    {
      std::cout << "WARNING: Got " << matchingppoints.size()
                << " common points of two facets! --> but just two points are used to check the "
                   "orientation!"
                << "-- if you see any problems, look into this!" << std::endl;

      std::map<unsigned int, unsigned int>::iterator bmatchingppoints1 = matchingppoints.begin();
      std::map<unsigned int, unsigned int>::iterator bmatchingppoints2 = matchingppoints.begin();
      unsigned int inlinecounter = 0;
      bmatchingppoints2++;
      while (matchingppoints.size() > 2)
      {
        if ((((bmatchingppoints1->first + 1) % (points1.size()) != bmatchingppoints2->first) &&
                ((bmatchingppoints2->first + 1) % (points1.size()) != bmatchingppoints1->first)) ||
            inlinecounter > matchingppoints.size())
        {
          matchingppoints.erase(bmatchingppoints1);
          bmatchingppoints1 = matchingppoints.begin();
          bmatchingppoints2 = matchingppoints.begin();
          bmatchingppoints2++;
          inlinecounter = 0;
          continue;
        }
        else  // in case all points are anyway on a line, it does not matter which point we delete!
        {
          inlinecounter++;
        }
        bmatchingppoints1++;
        bmatchingppoints2++;
        if (bmatchingppoints1 == matchingppoints.end()) bmatchingppoints1 = matchingppoints.begin();
        if (bmatchingppoints2 == matchingppoints.end()) bmatchingppoints2 = matchingppoints.begin();
      }
    }

    // 4// we have two common points left and check based on this information the orientation of the
    // facets!
    if (matchingppoints.size() == 2)
    {
      std::map<unsigned int, unsigned int>::iterator bmatchingppoints = matchingppoints.begin();
      bmatchingppoints++;


      bool firstordering = false;  // ordering of this facet
      if ((matchingppoints.begin()->first + 1) % (points1.size()) == bmatchingppoints->first)
        firstordering = true;
      else if ((matchingppoints.begin()->first) == (bmatchingppoints->first + 1) % (points1.size()))
        firstordering = false;
      else
        FOUR_C_THROW("Is there a point between the matched edges?");

      bool secondordering = false;  // ordering of the facet given as parameter
      if ((matchingppoints.begin()->second + 1) % (points2.size()) == bmatchingppoints->second)
        secondordering = true;
      else if ((matchingppoints.begin()->second) ==
               (bmatchingppoints->second + 1) % (points2.size()))
        secondordering = false;
      else
      {
        // it can be almost-hole facet, but touching at one point
        // in that case pointgraph produces two no-hole facets
        Side* parentside = f->ParentSide();
        const auto& facets = parentside->Facets();
        std::set<Point*> my_points(f->points_.begin(), f->points_.end());

        for (auto& facet : facets)
        {
          if (facet != f)
          {
            std::set<Point*> common;
            std::set<Point*> other_points(facet->points_.begin(), facet->points_.end());
            std::set_intersection(my_points.begin(), my_points.end(), other_points.begin(),
                other_points.end(), std::inserter(common, common.begin()));
            // if one is fully inside another
            if (common == other_points)
            {
              // we cannot find anything meaningful on this cycle -> it
              // includes all points of other facet, that is not a hole..
              // better try that facet instead
              return false;
            }
            else
            {
              FOUR_C_THROW("Check this case!");
            }
          }
        }
        FOUR_C_THROW("Is there a point between the matched edges?");
      }

      // ordering of the edge hast to be different, if the facets are ordered equal and vice versa!
      result = !(firstordering == secondordering);
      return true;
    }
  }
  // nowhere 2 common points in any cycle found --> facets have no common edge!
  return false;
}

Core::Geo::Cut::VolumeCell* Core::Geo::Cut::Facet::Neighbor(VolumeCell* cell)
{
  if (cells_.size() > 2) FOUR_C_THROW("can only have two neighbors");

  if (cells_.count(cell) == 0) FOUR_C_THROW("not my neighbor");

  for (plain_volumecell_set::iterator i = cells_.begin(); i != cells_.end(); ++i)
  {
    VolumeCell* vc = *i;
    if (vc != cell)
    {
      return vc;
    }
  }
  return nullptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::Neighbors(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& done, plain_volumecell_set& connected, plain_element_set& elements)
{
  for (plain_volumecell_set::iterator i = cells_.begin(); i != cells_.end(); ++i)
  {
    VolumeCell* c = *i;
    if (cells.count(c) > 0)
    {
      if (done.count(c) == 0 and connected.count(c) == 0 and
          elements.count(c->parent_element()) == 0)
      {
        connected.insert(c);
        elements.insert(c->parent_element());
        c->Neighbors(p, cells, done, connected, elements);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::Neighbors(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& done, plain_volumecell_set& connected)
{
  for (plain_volumecell_set::iterator i = cells_.begin(); i != cells_.end(); ++i)
  {
    VolumeCell* c = *i;
    if (cells.count(c) > 0)
    {
      if (done.count(c) == 0 and connected.count(c) == 0)
      {
        connected.insert(c);
        c->Neighbors(p, cells, done, connected);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::Equals(Core::FE::CellType distype)
{
  if (holes_.size() == 0)
  {
    find_corner_points();
    switch (distype)
    {
      case Core::FE::CellType::point1:
        return Kernel::IsValidPoint1(corner_points_);
        break;
      case Core::FE::CellType::line2:
        return Kernel::IsValidLine2(corner_points_);
        break;
      case Core::FE::CellType::quad4:
        return Kernel::IsValidQuad4(corner_points_);
        break;
      case Core::FE::CellType::tri3:
        return Kernel::IsValidTri3(corner_points_);
        break;
      default:
        FOUR_C_THROW("unsupported distype requested");
        break;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned Core::Geo::Cut::Facet::normal(const std::vector<Point*>& points,
    Core::LinAlg::Matrix<3, 1>& x1, Core::LinAlg::Matrix<3, 1>& x2, Core::LinAlg::Matrix<3, 1>& x3,
    Core::LinAlg::Matrix<3, 1>& b1, Core::LinAlg::Matrix<3, 1>& b2, Core::LinAlg::Matrix<3, 1>& b3)
{
  unsigned pointsize = points.size();
  if (pointsize < 3) return 0;

  points[0]->Coordinates(x1.data());
  points[1]->Coordinates(x2.data());

  b1.update(1, x2, -1, x1, 0);
  b1.scale(1. / b1.norm2());

  if (b1.norm2() < std::numeric_limits<double>::min())
    FOUR_C_THROW("same point in facet not supported");

  bool found = false;
  unsigned i = 2;
  for (; i < pointsize; ++i)
  {
    Point* p = points[i];
    p->Coordinates(x3.data());

    b2.update(1, x3, -1, x1, 0);
    b2.scale(1. / b2.norm2());

    // cross product to get the normal at the point
    b3(0) = b1(1) * b2(2) - b1(2) * b2(1);
    b3(1) = b1(2) * b2(0) - b1(0) * b2(2);
    b3(2) = b1(0) * b2(1) - b1(1) * b2(0);

    if (b3.norm2() > PLANARTOL)
    {
      found = true;
      break;
    }
  }
  if (not found)
  {
    // all on one line, no normal
    return 0;
  }

  b3.scale(1. / b3.norm2());
  return i;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::TriangulationPoints(PointSet& points)
{
  for (std::vector<std::vector<Point*>>::const_iterator i = triangulation_.begin();
       i != triangulation_.end(); ++i)
  {
    const std::vector<Point*>& tri = *i;
    std::copy(tri.begin(), tri.end(), std::inserter(points, points.begin()));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::NewPoint1Cell(Mesh& mesh, VolumeCell* volume,
    const std::vector<Point*>& points, plain_boundarycell_set& bcells)
{
  BoundaryCell* bc = mesh.NewPoint1Cell(volume, this, points);
  bcells.insert(bc);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::NewLine2Cell(Mesh& mesh, VolumeCell* volume,
    const std::vector<Point*>& points, plain_boundarycell_set& bcells)
{
  BoundaryCell* bc = mesh.NewLine2Cell(volume, this, points);
  bcells.insert(bc);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::NewTri3Cell(Mesh& mesh, VolumeCell* volume,
    const std::vector<Point*>& points, plain_boundarycell_set& bcells)
{
  BoundaryCell* bc = mesh.NewTri3Cell(volume, this, points);
  bcells.insert(bc);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::NewQuad4Cell(Mesh& mesh, VolumeCell* volume,
    const std::vector<Point*>& points, plain_boundarycell_set& bcells)
{
  if (mesh.CreateOptions().GenQuad4())
  {
    BoundaryCell* bc = mesh.NewQuad4Cell(volume, this, points);
    bcells.insert(bc);
  }
  else
  {  // split into two tri3 cells
    std::vector<Point*> tri3_points = points;
    tri3_points.pop_back();  // erase the last (fourth) point to obtain points (1,2,3)
    BoundaryCell* bc = mesh.NewTri3Cell(volume, this, tri3_points);
    bcells.insert(bc);
    tri3_points.erase(tri3_points.begin() + 1);  // erase the second point (1,3)
    tri3_points.push_back(points.back());        // add the last point again (1,3,4)
    bc = mesh.NewTri3Cell(volume, this, tri3_points);
    bcells.insert(bc);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::NewArbitraryCell(Mesh& mesh, VolumeCell* volume,
    const std::vector<Point*>& points, plain_boundarycell_set& bcells,
    const Core::FE::GaussIntegration& gp, const Core::LinAlg::Matrix<3, 1>& normal)
{
  BoundaryCell* bc = mesh.NewArbitraryCell(volume, this, points, gp, normal);
  bcells.insert(bc);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::GetBoundaryCells(plain_boundarycell_set& bcells)
{
  for (plain_volumecell_set::iterator vit = cells_.begin(); vit != cells_.end(); ++vit)
  {
    VolumeCell* vc = *vit;
    const plain_boundarycell_set& vbcells = vc->BoundaryCells();
    for (plain_boundarycell_set::const_iterator i = vbcells.begin(); i != vbcells.end(); ++i)
    {
      BoundaryCell* bc = *i;
      if (bc->GetFacet() == this) bcells.insert(bc);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::TestFacetArea(double tolerance, bool istetmeshintersection)
{
  if (OnCutSide() and cells_.size() > 1)
  {
    std::vector<double> area;
    area.reserve(2);
    // Does not necessarily have to be same on both sides of cutside.
    // Especially true for tetmeshintersection.
    for (plain_volumecell_set::iterator i = cells_.begin(); i != cells_.end(); ++i)
    {
      VolumeCell* vc = *i;

      area.push_back(0.);
      double& a = area.back();

      const plain_boundarycell_set& vbcells = vc->BoundaryCells();
      for (plain_boundarycell_set::const_iterator i = vbcells.begin(); i != vbcells.end(); ++i)
      {
        BoundaryCell* bc = *i;
        if (bc->GetFacet() == this)
        {
          a += bc->Area();
        }
      }
    }
    if (area.size() != 2)
    {
      FOUR_C_THROW("expect two volume cells at facet");
    }
    double diff = area[0] - area[1];
    if (fabs(diff) >= tolerance)
    {
      int eleId = (*cells_.begin())->parent_element()->Id();
      std::stringstream str;
      if (istetmeshintersection)
      {
        str << "In eleId=" << eleId << ", in TetMeshIntersection,"
            << " area mismatch: a1=" << area[0] << " a2=" << area[1] << " diff=" << diff;
      }
      else
      {
        str << "In eleId=" << eleId << " area mismatch: a1=" << area[0] << " a2=" << area[1]
            << " diff=" << diff;
      }

      std::cout << "WARNING: " << str.str() << "\n";
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::find_corner_points()
{
  if (corner_points_.size() == 0)
  {
    corner_points_ = points_;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::print(std::ostream& stream) const
{
  stream << "--- Facet ( address: " << this << " )\n";
  stream << "# Facet: "
         << "numpoints " << points_.size() << "\n# sideID: " << this->SideId()
         << "\n# OnCutSide: " << (OnCutSide() ? "TRUE" : "FALSE")
         << "\n# registered volume cells: ";
  std::copy(cells_.begin(), cells_.end(), std::ostream_iterator<VolumeCell*>(stream, " "));
  stream << "\n# points: ";
  std::copy(points_.begin(), points_.end(), std::ostream_iterator<Point*>(stream, " "));
  stream << "\n";
  for (std::vector<Point*>::const_iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point& p = **i;
    p.print(stream);
    stream << "\n";
  }
  stream << "\n";

  for (plain_facet_set::const_iterator i = holes_.begin(); i != holes_.end(); ++i)
  {
    Facet* hole = *i;
    hole->print(stream);
  }

  stream << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::IsTriangle(const std::vector<Point*>& tri) const
{
  if (tri.size() != 3) FOUR_C_THROW("three points expected");

  return points_.size() == 3 and not IsTriangulated() and not HasHoles() and Contains(tri);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::IsTriangulatedSide(const std::vector<Point*>& tri) const
{
  if (tri.size() != 3) FOUR_C_THROW("three points expected");

  for (std::vector<std::vector<Point*>>::const_iterator i = triangulation_.begin();
       i != triangulation_.end(); ++i)
  {
    const std::vector<Point*>& t = *i;
    bool found = true;
    for (std::vector<Point*>::const_iterator i = tri.begin(); i != tri.end(); ++i)
    {
      Point* p = *i;
      if (std::find(t.begin(), t.end(), p) == t.end())
      {
        found = false;
        break;
      }
    }
    if (found)
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned Core::Geo::Cut::Facet::NumPoints()
{
  unsigned numpoints = points_.size();
  if (IsTriangulated()) return numpoints + 1;
  for (plain_facet_set::iterator i = holes_.begin(); i != holes_.end(); ++i)
  {
    Facet* hole = *i;
    numpoints += hole->NumPoints();
  }
  return numpoints;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::Point* Core::Geo::Cut::Facet::OtherPoint(Point* p1, Point* p2)
{
  Point* result = nullptr;
  if (not HasHoles() and not IsTriangulated() and Points().size() == 3)
  {
    for (std::vector<Point*>::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = *i;
      if (p != p1 and p != p2)
      {
        if (result == nullptr)
        {
          result = p;
        }
        else
        {
          FOUR_C_THROW("point not unique");
        }
      }
    }
  }
  else
  {
    FOUR_C_THROW("plain triangular facet required");
  }
  return result;
}

/*----------------------------------------------------------------------------*
  return the local coordinates of corner points with respect to the given element
  if shadow=true, then the mapping is w.r. to the parent quad element from which
  this element is derived
*-----------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::CornerPointsLocal(
    Element* elem1, std::vector<std::vector<double>>& cornersLocal, bool shadow)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::Cut::Facet::CornerPointsLocal" );

  // TODO: this routine is expensive
  cornersLocal.clear();

  const std::vector<Point*>& corners = CornerPoints();

  for (std::vector<Point*>::const_iterator k = corners.begin(); k != corners.end(); k++)
  {
    std::vector<double> pt_local;
    const Point* po = *k;
    const double* coords = po->X();
    Core::LinAlg::Matrix<3, 1> glo, loc;

    glo(0, 0) = coords[0];
    glo(1, 0) = coords[1];
    glo(2, 0) = coords[2];

    // TODO: do we really need call a Newton here?
    if (shadow and elem1->isShadow())
      elem1->local_coordinates_quad(glo, loc);
    else
      elem1->local_coordinates(glo, loc);

    pt_local.push_back(loc(0, 0));
    pt_local.push_back(loc(1, 0));
    pt_local.push_back(loc(2, 0));

    cornersLocal.push_back(pt_local);
  }
}

/*----------------------------------------------------------------------------*
 * Return the global coordinates all of its corner points in order
 *                                                              sudhakar 05/15
 *----------------------------------------------------------------------------*/
std::vector<std::vector<double>> Core::Geo::Cut::Facet::CornerPointsGlobal(
    Element* elem1, bool shadow)
{
  const std::vector<Point*>& corners = CornerPoints();
  std::vector<std::vector<double>> cornersLocal;
  for (std::vector<Point*>::const_iterator k = corners.begin(); k != corners.end(); k++)
  {
    std::vector<double> pt_local(3);
    const Point* po = *k;
    const double* coords = po->X();

    for (unsigned dim = 0; dim < 3; dim++) pt_local[dim] = coords[dim];

    cornersLocal.push_back(pt_local);
  }
  return cornersLocal;
}

/*------------------------------------------------------------------------*
 *        Split the facet into a number of tri and quad cells
 *------------------------------------------------------------------------*/
void Core::Geo::Cut::Facet::SplitFacet(const std::vector<Point*>& facetpts)
{
  if (not this->HasHoles())
  {
    TriangulateFacet tf(facetpts);
    tf.SplitFacet();
    split_cells_ = tf.GetSplitCells();
  }
  else
  {
    std::vector<std::vector<Point*>> holepts;
    for (plain_facet_set::iterator ifacet = holes_.begin(); ifacet != holes_.end(); ++ifacet)
    {
      holepts.push_back((*ifacet)->corner_points_);
    }
    Core::Geo::Cut::TriangulateFacet tf(facetpts, holepts);
    tf.ear_clipping_with_holes(parentside_);
    split_cells_ = tf.GetSplitCells();
  }
}

/*----------------------------------------------------------------------------*
          Returns true if the facet is convex shaped
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::isConvex()
{
  if (this->HasHoles()) return false;

  Core::Geo::Cut::FacetShape geoType;
  std::vector<int> ptConcavity = Kernel::CheckConvexity(corner_points_, geoType, false, false);

  if (geoType == Core::Geo::Cut::Convex) return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Facet::belongs_to_level_set_side() { return parentside_->IsLevelSetSide(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, Core::Geo::Cut::Facet& f)
{
  stream << "facet: {";
  if (f.IsTriangulated())
  {
    const std::vector<std::vector<Core::Geo::Cut::Point*>>& triangulation = f.Triangulation();
    for (std::vector<std::vector<Core::Geo::Cut::Point*>>::const_iterator i = triangulation.begin();
         i != triangulation.end(); ++i)
    {
      const std::vector<Core::Geo::Cut::Point*>& tri = *i;
      stream << "{";
      for (std::vector<Core::Geo::Cut::Point*>::const_iterator i = tri.begin(); i != tri.end(); ++i)
      {
        Core::Geo::Cut::Point* p = *i;
        p->print(stream);
        stream << ",";
      }
      stream << "},";
    }
  }
  else
  {
    const std::vector<Core::Geo::Cut::Point*>& points = f.Points();
    for (std::vector<Core::Geo::Cut::Point*>::const_iterator i = points.begin(); i != points.end();
         ++i)
    {
      Core::Geo::Cut::Point* p = *i;
      p->print(stream);
      stream << ",";
    }
    if (f.HasHoles())
    {
      const Core::Geo::Cut::plain_facet_set& holes = f.Holes();
      for (Core::Geo::Cut::plain_facet_set::const_iterator i = holes.begin(); i != holes.end(); ++i)
      {
        Core::Geo::Cut::Facet& h = **i;
        stream << h;
      }
    }
  }
  stream << "}";
  return stream;
}

FOUR_C_NAMESPACE_CLOSE
