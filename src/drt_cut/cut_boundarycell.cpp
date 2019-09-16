/*---------------------------------------------------------------------*/
/*! \file

\brief cut boundary cell

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include "cut_boundarycell.H"
#include "cut_volumecell.H"
#include "cut_cycle.H"
#include "cut_element.H"
#include "cut_output.H"
#include "cut_side.H"

#include "../drt_geometry/element_volume.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::BoundaryCell::BoundaryCell(
    const Epetra_SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points)
    : facet_(facet), points_(Teuchos::rcp(new Cycle(points)))
{
  xyz_.Shape(3, xyz.N());

  /* This is necessary, because it is possible, that the given
   * Epetra_SerialDenseMatrix has a total row number smaller than 3 (equal the
   * actual problem dimension).                              hiermeier 11/16 */
  for (unsigned c = 0; c < static_cast<unsigned>(xyz.N()); ++c)
  {
    std::copy(&xyz(0, c), &xyz(0, c) + xyz.M(), &xyz_(0, c));
    std::fill(&xyz_(0, c) + xyz.M(), &xyz_(0, c) + 3, 0.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::BoundaryCell::Clear() { points_->clear(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::CUT::BoundaryCell::IsValid() const { return points_->size() > 0; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType celldistype>
void GEO::CUT::BoundaryCell::TransformLocalCoords(Element* elem1, const LINALG::Matrix<2, 1>& eta,
    LINALG::Matrix<3, 1>& x_gp_lin, LINALG::Matrix<3, 1>& normal, double& drs, bool shadow)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::BoundaryCell::TransformLocalCoords" );


  const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement;
  LINALG::Matrix<3, numnodes> xyzeGlo(this->xyz_, true), xyze;

  for (int i = 0; i < numnodes; i++)
  {
    LINALG::Matrix<3, 1> glo(&xyzeGlo(0, i)), loc;

    // map w.r to linear shadow element
    if (not shadow)
    {
      elem1->LocalCoordinates(glo, loc);
    }
    // map w.r to parent quad element
    else
    {
      if (not elem1->isShadow())
        elem1->LocalCoordinates(glo, loc);
      else
        elem1->LocalCoordinatesQuad(glo, loc);
    }
    xyze(0, i) = loc(0);
    xyze(1, i) = loc(1);
    xyze(2, i) = loc(2);
  }
  LINALG::Matrix<numnodes, 1> funct;
  LINALG::Matrix<2, numnodes> deriv;
  LINALG::Matrix<2, 2> metrictensor;

  DRT::UTILS::shape_function_2D(funct, eta(0), eta(1), celldistype);
  DRT::UTILS::shape_function_2D_deriv1(deriv, eta(0), eta(1), celldistype);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<celldistype>(
      xyze, deriv, metrictensor, drs, &normal);

  x_gp_lin.Multiply(xyze, funct);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<GEO::CUT::Point*>& GEO::CUT::BoundaryCell::Points() const { return (*points_)(); }

#if 0
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Tri3BoundaryCell::CollectCoordinates( const std::vector<Point*> & side,
                                                     std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  if ( side.size()!=3 )
  {
    throw std::runtime_error( "wrong number of points" );
  }

  plain_facet_set facets;
  FindCommonFacets( side[0], side[1], side[2], facets );

  Facet * f;
  if ( facets.size()==1 )
  {
    f = *facets.begin();
  }
  else
  {
    Facet * found = NULL;
    for ( plain_facet_set::iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->IsTriangle( side ) )
      {
        if ( found==NULL )
        {
          found = f;
        }
        else
        {
          throw std::runtime_error( "not unique" );
        }
      }
    }

    if ( found==NULL )
    {
      throw std::runtime_error( "not a valid cut side" );
    }

    f = found;
  }

  if ( not f->OnCutSide() )
  {
    return;
  }

  std::vector<Epetra_SerialDenseMatrix> & side_xyz = sides_xyz[f];
  side_xyz.push_back( Epetra_SerialDenseMatrix( 3, 3 ) );
  Epetra_SerialDenseMatrix & xyz = side_xyz.back();

  for ( int j=0; j<3; ++j )
  {
    Point * p = side[j];
    p->Coordinates( &xyz( 0, j ) );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Quad4BoundaryCell::CollectCoordinates(
    const std::vector<Point*> & side,
    std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz )
{
  if ( side.size()!=4 )
  {
    throw std::runtime_error( "wrong number of points" );
  }

  plain_facet_set facets;
  FindCommonFacets( side[0], side[1], side[2], side[3], facets );

  if ( facets.size()!=1 )
  {
    throw std::runtime_error( "not a valid cut side" );
  }

  Facet * f = *facets.begin();
  if ( not f->OnCutSide() )
  {
    return;
  }

  std::vector<Epetra_SerialDenseMatrix> & side_xyz = sides_xyz[f];
  side_xyz.push_back( Epetra_SerialDenseMatrix( 3, 4 ) );
  Epetra_SerialDenseMatrix & xyz = side_xyz.back();

  for ( int j=0; j<4; ++j )
  {
    Point * p = side[j];
    p->Coordinates( &xyz( 0, j ) );
  }
}
#endif

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double GEO::CUT::Line2BoundaryCell::Area() { return GEO::ElementArea(Shape(), xyz_); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double GEO::CUT::Tri3BoundaryCell::Area()
{
  const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

  const std::vector<Point*> points = this->Points();

  // create planes consisting of 3 nodes each
  LINALG::Matrix<numnodes, 1> p0(points[0]->X());
  LINALG::Matrix<numnodes, 1> p1(points[1]->X());
  LINALG::Matrix<numnodes, 1> p2(points[2]->X());

  LINALG::Matrix<numnodes, 1> v01;
  LINALG::Matrix<numnodes, 1> v02;
  LINALG::Matrix<numnodes, 1> v12;

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
      GEO::CUT::OUTPUT::GmshNewSection(file, "BoundaryCells");
      DumpGmsh(file);
      GEO::CUT::OUTPUT::GmshEndSection(file, false);
      GEO::CUT::OUTPUT::GmshNewSection(file, "BoundaryCellsNormal");
      DumpGmshNormal(file);
      GEO::CUT::OUTPUT::GmshEndSection(file);
      dserror("Boundary Cell not valid! Written GMSH output in Invalid_boundary_cell.pos!");
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
  //  - cross-product (see cut_tetmesh.cpp in IsValidTet())

  return area;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::BoundaryCell::DumpGmsh(std::ofstream& file, int* value)
{
  int default_value = facet_->SideId();
  if (not value) value = &default_value;

  OUTPUT::GmshCellDump(file, Shape(), xyz_, NULL, value);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Point1BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  // there is no normal for one point
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Tri3BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  file.precision(16);

  file << "VP(";
  LINALG::Matrix<3, 1> midpoint_triag(true);
  for (int i = 0; i < 3; ++i)
  {
    LINALG::Matrix<3, 1> cur(true);
    cur(0, 0) = xyz_(0, i);
    cur(1, 0) = xyz_(1, i);
    cur(2, 0) = xyz_(2, i);
    midpoint_triag.Update(1.0, cur, 1.0);
  }

  midpoint_triag.Scale(1.0 / 3.0);
  file << midpoint_triag(0, 0) << "," << midpoint_triag(1, 0) << "," << midpoint_triag(2, 0);
  file << "){";

  // Choose midpoint of triangle as normal for now. Not best choice possibly.
  LINALG::Matrix<2, 1> eta;
  eta(0, 0) = 0.5;
  eta(1, 0) = 0.5;
  LINALG::Matrix<3, 1> normal;
  this->Normal(eta, normal);
  file << normal(0) << "," << normal(1) << "," << normal(2);

  file << "};\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Quad4BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  file.precision(16);

  file << "VP(";
  LINALG::Matrix<3, 1> midpoint_quad(true);
  for (int i = 0; i < 4; ++i)
  {
    LINALG::Matrix<3, 1> cur(true);
    cur(0, 0) = xyz_(0, i);
    cur(1, 0) = xyz_(1, i);
    cur(2, 0) = xyz_(2, i);
    midpoint_quad.Update(1.0, cur, 1.0);
  }

  midpoint_quad.Scale(1.0 / 4.0);
  file << midpoint_quad(0, 0) << "," << midpoint_quad(1, 0) << "," << midpoint_quad(2, 0);
  file << "){";

  // Choose midpoint of triangle as normal for now. Not best choice possibly.
  LINALG::Matrix<2, 1> eta;
  eta(0, 0) = 0.5;
  eta(1, 0) = 0.5;
  LINALG::Matrix<3, 1> normal;
  this->Normal(eta, normal);
  file << normal(0) << "," << normal(1) << "," << normal(2);

  file << "};\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Line2BoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  file.precision(16);

  const unsigned num_nodes =
      DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::line2>::numNodePerElement;

  file << "VP(";
  LINALG::Matrix<3, 1> midpoint(true);
  for (unsigned i = 0; i < num_nodes; ++i)
  {
    LINALG::Matrix<3, 1> xyz(&xyz_(0, i), true);
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
  LINALG::Matrix<2, 1> eta;
  eta = 0.0;
  LINALG::Matrix<3, 1> normal;
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
void GEO::CUT::ArbitraryBoundaryCell::DumpGmshNormal(std::ofstream& file)
{
  // TO DO: implement gmsh output for arbitrarily shaped bcell
  //  dserror("not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Point1BoundaryCell::Normal(
    const LINALG::Matrix<2, 1>& xsi, LINALG::Matrix<3, 1>& normal) const
{
  normal.Scale(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Line2BoundaryCell::Normal(
    const LINALG::Matrix<2, 1>& xsi, LINALG::Matrix<3, 1>& normal) const
{
  const unsigned probdim = facet_->ParentSide()->ProbDim();
  switch (probdim)
  {
    case 2:
      EvalNormalVectors<2, DRT::Element::line2>(xyz_, xsi, normal);
      break;
    case 3:
      EvalNormalVectors<3, DRT::Element::line2>(xyz_, xsi, normal);
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Tri3BoundaryCell::Normal(
    const LINALG::Matrix<2, 1>& xsi, LINALG::Matrix<3, 1>& normal) const
{
#if 1
  // get derivatives at pos
  LINALG::Matrix<3, 3> side_xyze(xyz_.A(), true);

  LINALG::Matrix<2, 3> deriv;
  LINALG::Matrix<2, 3> A;

  DRT::UTILS::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), DRT::Element::tri3);
  A.MultiplyNT(deriv, side_xyze);

  // cross product to get the normal at the point
  normal(0) = A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1);
  normal(1) = A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2);
  normal(2) = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

  double norm = normal.Norm2();
  normal.Scale(1. / norm);
#else
  LINALG::Matrix<3, 1> x1(&xyz_(0, 0));
  LINALG::Matrix<3, 1> x2(&xyz_(0, 1));
  LINALG::Matrix<3, 1> x3(&xyz_(0, 2));

  LINALG::Matrix<3, 1> r1(true);
  LINALG::Matrix<3, 1> r2(true);

  r1.Update(1., x1, -1., x3, 0.);
  r2.Update(1., x2, -1., x1, 0.);

  //  r1.Update( +1., x3, -1., x2, 0.);
  //  r2.Update( +1., x1, -1., x3, 0.);

  //  r1.Update( +1., x2, -1., x1, 0.);
  //  r2.Update( +1., x3, -1., x2, 0.);


  r2.Scale(1.0 / r2.Norm2());
  r1.Scale(1.0 / r1.Norm2());

  // cross product to get the normal at the point
  normal(0) = r1(1) * r2(2) - r1(2) * r2(1);
  normal(1) = r1(2) * r2(0) - r1(0) * r2(2);
  normal(2) = r1(0) * r2(1) - r1(1) * r2(0);

  double norm = normal.Norm2();
  normal.Scale(1. / norm);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Quad4BoundaryCell::Normal(
    const LINALG::Matrix<2, 1>& xsi, LINALG::Matrix<3, 1>& normal) const
{
  // get derivatives at pos
  LINALG::Matrix<3, 4> side_xyze(xyz_.A(), true);
  // Position2d<DRT::Element::quad4> position( side_xyze, xsi );
  // position.Normal( xsi, normal );

  LINALG::Matrix<2, 4> deriv;
  LINALG::Matrix<2, 3> A;

  DRT::UTILS::shape_function_2D_deriv1(deriv, xsi(0), xsi(1), DRT::Element::quad4);
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
void GEO::CUT::ArbitraryBoundaryCell::Normal(
    const LINALG::Matrix<2, 1>& xsi, LINALG::Matrix<3, 1>& normal) const
{
  dserror("Call GetNormalVector() to get normal for arbitrary boundarycells");
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
DRT::UTILS::GaussIntegration GEO::CUT::Point1BoundaryCell::gaussRule()
{
  DRT::UTILS::GaussIntegration gi(DRT::Element::point1, CubatureDegree(DRT::Element::point1));
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::UTILS::GaussIntegration GEO::CUT::Line2BoundaryCell::gaussRule()
{
  DRT::UTILS::GaussIntegration gi(DRT::Element::line2, CubatureDegree(DRT::Element::line2));
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::UTILS::GaussIntegration GEO::CUT::Tri3BoundaryCell::gaussRule()
{
  DRT::UTILS::GaussIntegration gi(DRT::Element::tri3, CubatureDegree(DRT::Element::tri3));
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::UTILS::GaussIntegration GEO::CUT::Quad4BoundaryCell::gaussRule()
{
  DRT::UTILS::GaussIntegration gi(DRT::Element::quad4, CubatureDegree(DRT::Element::quad4));
  return gi;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::UTILS::GaussIntegration GEO::CUT::ArbitraryBoundaryCell::gaussRule() { return gaussRule_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::ArbitraryBoundaryCell::ElementCenter(LINALG::Matrix<3, 1>& midpoint)
{
  dserror("Element Center for ArbitraryBoundaryCells not implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Point1BoundaryCell::ElementCenter(LINALG::Matrix<3, 1>& midpoint)
{
  std::copy(xyz_.A(), xyz_.A() + 3, midpoint.A());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Line2BoundaryCell::ElementCenter(LINALG::Matrix<3, 1>& midpoint)
{
  LINALG::Matrix<3, 1> center_rst(true);
  MyElementCenter<DRT::Element::line2>(center_rst, midpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Tri3BoundaryCell::ElementCenter(LINALG::Matrix<3, 1>& midpoint)
{
  LINALG::Matrix<3, 1> center;
  center(0, 0) = 0.25;
  center(1, 0) = 0.25;
  center(2, 0) = 0.0;
  MyElementCenter<DRT::Element::tri3>(center, midpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::Quad4BoundaryCell::ElementCenter(LINALG::Matrix<3, 1>& midpoint)
{
  LINALG::Matrix<3, 1> center;
  center(0, 0) = 0.0;
  center(1, 0) = 0.0;
  center(2, 0) = 0.0;
  MyElementCenter<DRT::Element::quad4>(center, midpoint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::CUT::Point1BoundaryCell::GetNormalVector()
{
  dserror("There is no normal for Point1 boundarycell");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::CUT::Line2BoundaryCell::GetNormalVector()
{
  LINALG::Matrix<3, 1> normal(true);
  LINALG::Matrix<2, 1> xsi(true);

  Normal(xsi, normal);

  return normal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::CUT::Tri3BoundaryCell::GetNormalVector()
{
  dserror("Call Transform function to get normal for Tri3 boundarycell");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::CUT::Quad4BoundaryCell::GetNormalVector()
{
  dserror("Call Transform function to get normal for Quad4 boundarycell");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::CUT::ArbitraryBoundaryCell::GetNormalVector() { return normal_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::BoundaryCell::Print(std::ostream& stream)
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
std::vector<std::vector<double>> GEO::CUT::BoundaryCell::CoordinatesV()
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
bool GEO::CUT::Tri3BoundaryCell::IsValidBoundaryCell()
{
  const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

  const std::vector<Point*> points = this->Points();

  // create planes consisting of 3 nodes each
  LINALG::Matrix<numnodes, 1> p0(points[0]->X());
  LINALG::Matrix<numnodes, 1> p1(points[1]->X());
  LINALG::Matrix<numnodes, 1> p2(points[2]->X());

  LINALG::Matrix<numnodes, 1> v01;
  LINALG::Matrix<numnodes, 1> v02;
  LINALG::Matrix<numnodes, 1> v12;

  v01.Update(1, p1, -1, p0, 0);
  v02.Update(1, p2, -1, p0, 0);
  v12.Update(1, p1, -1, p2, 0);

  // Get distance to origin
  LINALG::Matrix<numnodes, 1> temp(true);
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
template void GEO::CUT::BoundaryCell::TransformLocalCoords<DRT::Element::tri3>(Element* elem1,
    const LINALG::Matrix<2, 1>& eta, LINALG::Matrix<3, 1>& x_gp_lin, LINALG::Matrix<3, 1>& normal,
    double& drs, bool shadow);
template void GEO::CUT::BoundaryCell::TransformLocalCoords<DRT::Element::quad4>(Element* elem1,
    const LINALG::Matrix<2, 1>& eta, LINALG::Matrix<3, 1>& x_gp_lin, LINALG::Matrix<3, 1>& normal,
    double& drs, bool shadow);
