/*---------------------------------------------------------------------*/
/*! \file

\brief cut boundary cell

\level 3


*----------------------------------------------------------------------*/

#include "4C_cut_boundarycell.hpp"

#include "4C_cut_cycle.hpp"
#include "4C_cut_element.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_discretization_geometry_element_volume.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundaryCell::BoundaryCell(
    const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points)
    : facet_(facet), points_(Teuchos::rcp(new Cycle(points)))
{
  xyz_.shape(3, xyz.numCols());

  /* This is necessary, because it is possible, that the given
   * LinAlg::SerialDenseMatrix has a total row number smaller than 3 (equal the
   * actual problem dimension).                              hiermeier 11/16 */
  for (unsigned c = 0; c < static_cast<unsigned>(xyz.numCols()); ++c)
  {
    std::copy(&xyz(0, c), &xyz(0, c) + xyz.numRows(), &xyz_(0, c));
    std::fill(&xyz_(0, c) + xyz.numRows(), &xyz_(0, c) + 3, 0.0);
  }

  // Assign reference position of boundary cell
  xyz_ref_ = xyz_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundaryCell::Clear() { points_->clear(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::BoundaryCell::IsValid() const { return points_->size() > 0; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <Core::FE::CellType celldistype>
void Core::Geo::Cut::BoundaryCell::transform_local_coords(Element* elem1,
    const Core::LinAlg::Matrix<2, 1>& eta, Core::LinAlg::Matrix<3, 1>& x_gp_lin,
    Core::LinAlg::Matrix<3, 1>& normal, double& drs, bool shadow)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::Cut::BoundaryCell::transform_local_coords" );


  const int numnodes = Core::FE::num_nodes<celldistype>;
  Core::LinAlg::Matrix<3, numnodes> xyzeGlo(this->xyz_, true), xyze;

  for (int i = 0; i < numnodes; i++)
  {
    Core::LinAlg::Matrix<3, 1> glo(&xyzeGlo(0, i)), loc;

    // map w.r to linear shadow element
    if (not shadow)
    {
      elem1->local_coordinates(glo, loc);
    }
    // map w.r to parent quad element
    else
    {
      if (not elem1->isShadow())
        elem1->local_coordinates(glo, loc);
      else
        elem1->local_coordinates_quad(glo, loc);
    }
    xyze(0, i) = loc(0);
    xyze(1, i) = loc(1);
    xyze(2, i) = loc(2);
  }
  Core::LinAlg::Matrix<numnodes, 1> funct;
  Core::LinAlg::Matrix<2, numnodes> deriv;
  Core::LinAlg::Matrix<2, 2> metrictensor;

  Core::FE::shape_function_2D(funct, eta(0), eta(1), celldistype);
  Core::FE::shape_function_2D_deriv1(deriv, eta(0), eta(1), celldistype);
  Core::FE::ComputeMetricTensorForBoundaryEle<celldistype>(xyze, deriv, metrictensor, drs, &normal);

  x_gp_lin.Multiply(xyze, funct);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Core::Geo::Cut::Point*>& Core::Geo::Cut::BoundaryCell::Points() const
{
  return (*points_)();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Core::Geo::Cut::Line2BoundaryCell::Area() { return Core::Geo::ElementArea(Shape(), xyz_); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Core::Geo::Cut::Tri3BoundaryCell::Area()
{
  const int numnodes = Core::FE::num_nodes<Core::FE::CellType::tri3>;

  const std::vector<Point*> points = this->Points();

  // create planes consisting of 3 nodes each
  Core::LinAlg::Matrix<numnodes, 1> p0(points[0]->X());
  Core::LinAlg::Matrix<numnodes, 1> p1(points[1]->X());
  Core::LinAlg::Matrix<numnodes, 1> p2(points[2]->X());

  Core::LinAlg::Matrix<numnodes, 1> v01;
  Core::LinAlg::Matrix<numnodes, 1> v02;
  Core::LinAlg::Matrix<numnodes, 1> v12;

  v01.Update(1, p1, -1, p0, 0);
  v02.Update(1, p2, -1, p0, 0);
  v12.Update(1, p1, -1, p2, 0);


  double a = v01.Norm2();
  double b = v02.Norm2();
  double c = v12.Norm2();

  // Herons formula:
  //    sqrt((a+(b+c))(c-(a-b))(c+(a-b))(a+(b-c)))/4
  //    double s=0.5*(v01.Norm2() + v02.Norm2() + v12.Norm2());
  //    double area = sqrt( s*(s-v01.Norm2())*(s-v02.Norm2())*(s-v12.Norm2()) );

  double areasqr = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
  if (areasqr < 0.0)
  {
    if (std::abs(areasqr) > AREA_BC_TOL)
    {
      std::ofstream file("invalid_boundary_cell.pos");
      // Write BCs for outside VolumeCell:
      Core::Geo::Cut::Output::GmshNewSection(file, "BoundaryCells");
      DumpGmsh(file);
      Core::Geo::Cut::Output::GmshEndSection(file, false);
      Core::Geo::Cut::Output::GmshNewSection(file, "BoundaryCellsNormal");
      DumpGmshNormal(file);
      Core::Geo::Cut::Output::GmshEndSection(file);
      FOUR_C_THROW("Boundary Cell not valid! Written GMSH output in Invalid_boundary_cell.pos!");
    }
    else
    {
#if EXTENDED_CUT_DEBUG_OUTPUT
      std::cout << "NOTE: Boundary cell area calculated with Heron's formula is negative"
                << std::endl;
#endif
      // TODO: What should we really do here?
      // For now we create really small boundary cells
      areasqr *= -1;
    }
  }
  double area = 0.25 * sqrt(areasqr);

  // Alternative options:
  //  - cross-product (see cut_tetmesh.cpp in is_valid_tet())

  return area;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundaryCell::DumpGmsh(std::ofstream& file, int* value)
{
  int default_value = facet_->SideId();
  if (not value) value = &default_value;

  Output::GmshCellDump(file, Shape(), xyz_, nullptr, value);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Point1BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  // there is no normal for one point
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Tri3BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  file.precision(16);

  file << "VP(";
  Core::LinAlg::Matrix<3, 1> midpoint_triag(true);
  for (int i = 0; i < 3; ++i)
  {
    Core::LinAlg::Matrix<3, 1> cur(true);
    cur(0, 0) = xyz_(0, i);
    cur(1, 0) = xyz_(1, i);
    cur(2, 0) = xyz_(2, i);
    midpoint_triag.Update(1.0, cur, 1.0);
  }

  midpoint_triag.Scale(1.0 / 3.0);
  file << midpoint_triag(0, 0) << "," << midpoint_triag(1, 0) << "," << midpoint_triag(2, 0);
  file << "){";

  // Choose midpoint of triangle as normal for now. Not best choice possibly.
  Core::LinAlg::Matrix<2, 1> eta;
  eta(0, 0) = 0.5;
  eta(1, 0) = 0.5;
  Core::LinAlg::Matrix<3, 1> normal;
  this->Normal(eta, normal);
  file << normal(0) << "," << normal(1) << "," << normal(2);

  file << "};\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Quad4BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  file.precision(16);

  file << "VP(";
  Core::LinAlg::Matrix<3, 1> midpoint_quad(true);
  for (int i = 0; i < 4; ++i)
  {
    Core::LinAlg::Matrix<3, 1> cur(true);
    cur(0, 0) = xyz_(0, i);
    cur(1, 0) = xyz_(1, i);
    cur(2, 0) = xyz_(2, i);
    midpoint_quad.Update(1.0, cur, 1.0);
  }

  midpoint_quad.Scale(1.0 / 4.0);
  file << midpoint_quad(0, 0) << "," << midpoint_quad(1, 0) << "," << midpoint_quad(2, 0);
  file << "){";

  // Choose midpoint of triangle as normal for now. Not best choice possibly.
  Core::LinAlg::Matrix<2, 1> eta;
  eta(0, 0) = 0.5;
  eta(1, 0) = 0.5;
  Core::LinAlg::Matrix<3, 1> normal;
  this->Normal(eta, normal);
  file << normal(0) << "," << normal(1) << "," << normal(2);

  file << "};\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Line2BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  file.precision(16);

  const unsigned num_nodes = Core::FE::num_nodes<Core::FE::CellType::line2>;

  file << "VP(";
  Core::LinAlg::Matrix<3, 1> midpoint(true);
  for (unsigned i = 0; i < num_nodes; ++i)
  {
    Core::LinAlg::Matrix<3, 1> xyz(&xyz_(0, i), true);
    midpoint.Update(1.0, xyz, 1.0);
  }

  midpoint.Scale(1.0 / num_nodes);
  for (unsigned j = 0; j < 3; ++j)
  {
    if (j > 0) file << ",";
    file << midpoint(j, 0);
  }
  file << "){";

  // Choose midpoint of line2
  Core::LinAlg::Matrix<2, 1> eta;
  eta = 0.0;
  Core::LinAlg::Matrix<3, 1> normal;
  this->Normal(eta, normal);
  for (unsigned j = 0; j < 3; ++j)
  {
    if (j > 0) file << ",";
    file << normal(j, 0);
  }

  file << "};\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::ArbitraryBoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  // TO DO: implement gmsh output for arbitrarily shaped bcell
  //  FOUR_C_THROW("not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Point1BoundaryCell::Normal(
    const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const
{
  normal.putScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Line2BoundaryCell::Normal(
    const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const
{
  const unsigned probdim = facet_->ParentSide()->ProbDim();
  switch (probdim)
  {
    case 2:
      EvalNormalVectors<2, Core::FE::CellType::line2>(xyz_, xsi, normal);
      break;
    case 3:
      EvalNormalVectors<3, Core::FE::CellType::line2>(xyz_, xsi, normal);
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Tri3BoundaryCell::Normal(
    const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const
{
  // get derivatives at pos
  Core::LinAlg::Matrix<3, 3> side_xyze(xyz_.values(), true);

  Core::LinAlg::Matrix<2, 3> deriv;
  Core::LinAlg::Matrix<2, 3> A;

  Core::FE::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), Core::FE::CellType::tri3);
  A.MultiplyNT(deriv, side_xyze);

  // cross product to get the normal at the point
  normal(0) = A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1);
  normal(1) = A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2);
  normal(2) = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

  double norm = normal.Norm2();
  normal.Scale(1. / norm);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Quad4BoundaryCell::Normal(
    const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const
{
  // get derivatives at pos
  Core::LinAlg::Matrix<3, 4> side_xyze(xyz_.values(), true);
  // Position2d<Core::FE::CellType::quad4> position( side_xyze, xsi );
  // position.Normal( xsi, normal );

  Core::LinAlg::Matrix<2, 4> deriv;
  Core::LinAlg::Matrix<2, 3> A;

  Core::FE::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), Core::FE::CellType::quad4);
  A.MultiplyNT(deriv, side_xyze);

  // cross product to get the normal at the point
  normal(0) = A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1);
  normal(1) = A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2);
  normal(2) = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

  double norm = normal.Norm2();
  normal.Scale(1. / norm);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::ArbitraryBoundaryCell::Normal(
    const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const
{
  FOUR_C_THROW("Call GetNormalVector() to get normal for arbitrary boundarycells");
  exit(1);
  /*// cross product to get the normal at the point
  normal( 0 ) = A( 0, 1 )*A( 1, 2 ) - A( 0, 2 )*A( 1, 1 );
  normal( 1 ) = A( 0, 2 )*A( 1, 0 ) - A( 0, 0 )*A( 1, 2 );
  normal( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 0, 1 )*A( 1, 0 );*/

  //   double norm = normal.Norm2();
  //   normal.Scale( 1./norm );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::GaussIntegration Core::Geo::Cut::Point1BoundaryCell::gaussRule(int cubaturedegree)
{
  Core::FE::GaussIntegration gi(Core::FE::CellType::point1, cubaturedegree);
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::GaussIntegration Core::Geo::Cut::Line2BoundaryCell::gaussRule(int cubaturedegree)
{
  Core::FE::GaussIntegration gi(Core::FE::CellType::line2, cubaturedegree);
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::GaussIntegration Core::Geo::Cut::Tri3BoundaryCell::gaussRule(int cubaturedegree)
{
  Core::FE::GaussIntegration gi(Core::FE::CellType::tri3, cubaturedegree);
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::GaussIntegration Core::Geo::Cut::Quad4BoundaryCell::gaussRule(int cubaturedegree)
{
  Core::FE::GaussIntegration gi(Core::FE::CellType::quad4, cubaturedegree);
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::GaussIntegration Core::Geo::Cut::ArbitraryBoundaryCell::gaussRule(int cubaturedegree)
{
  return gauss_rule_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::ArbitraryBoundaryCell::element_center(Core::LinAlg::Matrix<3, 1>& midpoint)
{
  FOUR_C_THROW("Element Center for ArbitraryBoundaryCells not implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Point1BoundaryCell::element_center(Core::LinAlg::Matrix<3, 1>& midpoint)
{
  std::copy(xyz_.values(), xyz_.values() + 3, midpoint.A());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Line2BoundaryCell::element_center(Core::LinAlg::Matrix<3, 1>& midpoint)
{
  Core::LinAlg::Matrix<3, 1> center_rst(true);
  my_element_center<Core::FE::CellType::line2>(center_rst, midpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Tri3BoundaryCell::element_center(Core::LinAlg::Matrix<3, 1>& midpoint)
{
  Core::LinAlg::Matrix<3, 1> center;
  center(0, 0) = 0.25;
  center(1, 0) = 0.25;
  center(2, 0) = 0.0;
  my_element_center<Core::FE::CellType::tri3>(center, midpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Quad4BoundaryCell::element_center(Core::LinAlg::Matrix<3, 1>& midpoint)
{
  Core::LinAlg::Matrix<3, 1> center;
  center(0, 0) = 0.0;
  center(1, 0) = 0.0;
  center(2, 0) = 0.0;
  my_element_center<Core::FE::CellType::quad4>(center, midpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::Matrix<3, 1> Core::Geo::Cut::Point1BoundaryCell::GetNormalVector()
{
  FOUR_C_THROW("There is no normal for Point1 boundarycell");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::Matrix<3, 1> Core::Geo::Cut::Line2BoundaryCell::GetNormalVector()
{
  Core::LinAlg::Matrix<3, 1> normal(true);
  Core::LinAlg::Matrix<2, 1> xsi(true);

  Normal(xsi, normal);

  return normal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::Matrix<3, 1> Core::Geo::Cut::Tri3BoundaryCell::GetNormalVector()
{
  FOUR_C_THROW("Call Transform function to get normal for Tri3 boundarycell");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::Matrix<3, 1> Core::Geo::Cut::Quad4BoundaryCell::GetNormalVector()
{
  FOUR_C_THROW("Call Transform function to get normal for Quad4 boundarycell");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::LinAlg::Matrix<3, 1> Core::Geo::Cut::ArbitraryBoundaryCell::GetNormalVector()
{
  return normal_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundaryCell::Print(std::ostream& stream)
{
  stream << "--- boundary cell ( address: " << std::setw(10) << this << " )\n";
  for (unsigned i = 0; i < points_->size(); i++)
  {
    (*points_)()[i]->Print(stream);
    stream << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<std::vector<double>> Core::Geo::Cut::BoundaryCell::CoordinatesV()
{
  std::vector<std::vector<double>> corners;
  for (std::vector<Point*>::const_iterator j = (*points_)().begin(); j != (*points_)().end(); ++j)
  {
    std::vector<double> cornerLocal(3);
    cornerLocal.assign((*j)->X(), (*j)->X() + 3);
    corners.push_back(cornerLocal);
  }
  return corners;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Tri3BoundaryCell::IsValidBoundaryCell()
{
  const int numnodes = Core::FE::num_nodes<Core::FE::CellType::tri3>;

  const std::vector<Point*> points = this->Points();

  // create planes consisting of 3 nodes each
  Core::LinAlg::Matrix<numnodes, 1> p0(points[0]->X());
  Core::LinAlg::Matrix<numnodes, 1> p1(points[1]->X());
  Core::LinAlg::Matrix<numnodes, 1> p2(points[2]->X());

  Core::LinAlg::Matrix<numnodes, 1> v01;
  Core::LinAlg::Matrix<numnodes, 1> v02;
  Core::LinAlg::Matrix<numnodes, 1> v12;

  v01.Update(1, p1, -1, p0, 0);
  v02.Update(1, p2, -1, p0, 0);
  v12.Update(1, p1, -1, p2, 0);

  // Get distance to origin
  Core::LinAlg::Matrix<numnodes, 1> temp(true);
  temp(0, 0) = p0.Norm2();  // Distance of points to origin
  temp(1, 0) = p1.Norm2();
  temp(2, 0) = p2.Norm2();

  // This is to scale the tolerance, it determines our maximum precision (i.e. machine precision)
  //  double max_dist_to_orgin = temp.NormInf();

  // Distance between points in triangle
  temp(0, 0) = v01.Norm2();
  temp(1, 0) = v02.Norm2();
  temp(2, 0) = v12.Norm2();

  double min_dist_in_tri = temp.MinValue();
  // We want to test with this one I think... But might lead to problems.
  //  double tolerance = LINSOLVETOL*max_dist_to_orgin;

  temp(0, 0) = points[0]->Tolerance();
  temp(1, 0) = points[1]->Tolerance();
  temp(2, 0) = points[2]->Tolerance();
  double max_tol_points = temp.MaxValue();


  //  std::cout << "min_dist_in_tri: " << min_dist_in_tri << std::endl;
  //  std::cout << "tolerance: " << tolerance << std::endl;
  //  std::cout << "max_tol_points: " << max_tol_points << std::endl;

  if (min_dist_in_tri < max_tol_points)
  {
    return false;
  }

  return true;
}

// function specializations
template void Core::Geo::Cut::BoundaryCell::transform_local_coords<Core::FE::CellType::tri3>(
    Element* elem1, const Core::LinAlg::Matrix<2, 1>& eta, Core::LinAlg::Matrix<3, 1>& x_gp_lin,
    Core::LinAlg::Matrix<3, 1>& normal, double& drs, bool shadow);
template void Core::Geo::Cut::BoundaryCell::transform_local_coords<Core::FE::CellType::quad4>(
    Element* elem1, const Core::LinAlg::Matrix<2, 1>& eta, Core::LinAlg::Matrix<3, 1>& x_gp_lin,
    Core::LinAlg::Matrix<3, 1>& normal, double& drs, bool shadow);

FOUR_C_NAMESPACE_CLOSE
