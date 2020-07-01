/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to surface geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_H_
#define UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_H_


#include <cxxtest/TestSuite.h>

#include "src/drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "src/drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.H"
#include "src/drt_geometry_pair/geometry_pair_element_functions.H"
#include "src/drt_geometry_pair/geometry_pair_utility_classes.H"
#include "src/drt_so3/so_surface.H"
#include "src/drt_beam3/beam3r.H"
#include "src/drt_inpar/inpar_beam_to_solid.H"

#include "src/drt_geometry_pair/unit_geometry_pair_line_to_surface_geometry.H"


namespace GEOMETRYPAIR
{
  class GeometryPairLineToSurface_TestSuite;
}

/**
 * Class to test the line to volume geometry pair segmentation algorithm.
 */
class GEOMETRYPAIR::GeometryPairLineToSurface_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * Set up the testing environment.
   */
  void setUp()
  {
    // Set up the evaluation data container for the geometry pairs.
    Teuchos::ParameterList line_to_surface_params_list;
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_surface_params_list);
    evaluation_data_ =
        Teuchos::rcp(new GEOMETRYPAIR::LineToSurfaceEvaluationData(line_to_surface_params_list));
  }

  /**
   * Delete pointers and other class variables.
   */
  void tearDown()
  {
    // Dereference the pointers.
    evaluation_data_ = Teuchos::null;
  }

  /**
   * Test the projection of a point to a tri3 surface, with default normals on the surface.
   */
  void test_point_to_surface_projection_tri3()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    xtest_setup_tri3(q_solid);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.3;
    point(1) = 0.1;
    point(2) = 0.2;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3436484045755569;
    xi_result(1) = 0.2877784467188441;
    xi_result(2) = 0.03189763881277458;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri3 surface, with given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_tri3()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    LINALG::Matrix<9, 1, double> nodal_normals(true);
    xtest_setup_tri3(q_solid, &nodal_normals);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.3;
    point(1) = 0.1;
    point(2) = 0.2;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, &nodal_normals);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3457692493957274;
    xi_result(1) = 0.2853120425437799;
    xi_result(2) = 0.03218342274405913;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri6 surface, with default normals on the surface.
   */
  void test_point_to_surface_projection_tri6()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    xtest_setup_tri6(q_solid);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.3;
    point(1) = 0.1;
    point(2) = 0.2;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.1935801417994475;
    xi_result(1) = 0.1678155116663445;
    xi_result(2) = 0.236826220497202;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri3 surface, with given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_tri6()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    LINALG::Matrix<18, 1, double> nodal_normals(true);
    xtest_setup_tri6(q_solid, &nodal_normals);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.3;
    point(1) = 0.1;
    point(2) = 0.2;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, &nodal_normals);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3274411842809972;
    xi_result(1) = 0.1649919700896869;
    xi_result(2) = 0.2749865824042791;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad4 surface, with default normals on the surface.
   */
  void test_point_to_surface_projection_quad4()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the normals on the solid element.
    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid(true);
    xtest_setup_quad4(q_solid);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.5856297224156624;
    xi_result(1) = -0.2330351551569786;
    xi_result(2) = 0.1132886291998745;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad4 surface, with given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_quad4()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<12, 1, double> q_solid(true);
    LINALG::Matrix<12, 1, double> nodal_normals(true);
    xtest_setup_quad4(q_solid, &nodal_normals);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, &nodal_normals);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.6306816217205055;
    xi_result(1) = -0.2391123963538002;
    xi_result(2) = 0.1168739495183324;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad8 surface, with default normals on the surface.
   */
  void test_point_to_surface_projection_quad8()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    xtest_setup_quad8(q_solid);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.4869140501387866;
    xi_result(1) = -0.6545313748232923;
    xi_result(2) = 0.4772682324027889;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad8 surface, with given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_quad8()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    LINALG::Matrix<24, 1, double> nodal_normals;
    xtest_setup_quad8(q_solid, &nodal_normals);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, &nodal_normals);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = -0.167932271257968;
    xi_result(1) = 0.1593451990533972;
    xi_result(2) = 0.6729448863050194;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad9 surface, with default normals on the surface.
   */
  void test_point_to_surface_projection_quad9()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    xtest_setup_quad9(q_solid);

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.4374951399531939;
    xi_result(1) = -0.4006486973745378;
    xi_result(2) = 0.2412946023554158;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad9 surface, with given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_quad9()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    LINALG::Matrix<27, 1, double> nodal_normals;
    xtest_setup_quad9(q_solid, &nodal_normals);


    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, q_solid, xi, projection_result, &nodal_normals);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3784195771508677;
    xi_result(1) = -0.436333510864013;
    xi_result(2) = 0.2483249147920992;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the intersection of a line with a tri3 surface, with default normals on the surface.
   */
  void test_line_to_surface_intersection_tri3()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    xtest_setup_tri3(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.0;
    xi_result(0, 1) = 0.5441734719700435;
    xi_result(1, 0) = 0.1074360140351795;
    xi_result(1, 1) = 0.4558265280299565;
    xi_result(2, 0) = 0.1140207710811362;
    xi_result(2, 1) = 0.00821450263257107;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.959535845440973;
    eta_result(1) = -0.2754895911921936;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri3 surface, with given normals on the nodes.
   */
  void test_line_to_surface_intersection_normal_interpolation_tri3()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    LINALG::Matrix<9, 1, double> nodal_normals;
    xtest_setup_tri3(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.;
    xi_result(0, 1) = 0.5449151431901401;
    xi_result(1, 0) = 0.0892976752542103;
    xi_result(1, 1) = 0.4550848568098599;
    xi_result(2, 0) = 0.1071908576829917;
    xi_result(2, 1) = 0.00852036464820085;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.933108361186308;
    eta_result(1) = -0.2769233373990823;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri6 surface, with default normals on the surface.
   */
  void test_line_to_surface_intersection_tri6()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    xtest_setup_tri6(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.0;
    xi_result(0, 1) = 0.661306368091275;
    xi_result(1, 0) = 0.1351724121757158;
    xi_result(1, 1) = 0.338693631908725;
    xi_result(2, 0) = 0.1130371451881858;
    xi_result(2, 1) = 0.133409588649314;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.845498535448603;
    eta_result(1) = -0.1960742371555871;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri6 surface, with given normals on the nodes.
   */
  void test_line_to_surface_intersection_normal_interpolation_tri6()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    LINALG::Matrix<18, 1, double> nodal_normals;
    xtest_setup_tri6(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.;
    xi_result(0, 1) = 0.6584629848688872;
    xi_result(1, 0) = 0.1326786387805501;
    xi_result(1, 1) = 0.3415370151311128;
    xi_result(2, 0) = 0.1167772617143948;
    xi_result(2, 1) = 0.117654537323362;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.834928214700044;
    eta_result(1) = -0.1707134503670001;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad4 surface, with default normals on the surface.
   */
  void test_line_to_surface_intersection_quad4()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid;
    xtest_setup_quad4(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.785985513536155;
    xi_result(1, 1) = 0.0135117312962169;
    xi_result(2, 0) = 0.113108951013877;
    xi_result(2, 1) = 0.1177337444785567;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.957076826689831;
    eta_result(1) = 0.4600569936643898;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad4 surface, with given normals on the nodes.
   */
  void test_line_to_surface_intersection_normal_interpolation_quad4()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid;
    LINALG::Matrix<12, 1, double> nodal_normals;
    xtest_setup_quad4(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.825474249880623;
    xi_result(1, 1) = -0.01145366341249682;
    xi_result(2, 0) = 0.107340226468075;
    xi_result(2, 1) = 0.119547807682323;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.932605711413924;
    eta_result(1) = 0.4202318513645913;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad8 surface, with default normals on the surface.
   */
  void test_line_to_surface_intersection_quad8()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    xtest_setup_quad8(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.7289003389787947;
    xi_result(1, 1) = -0.2401689430824591;
    xi_result(2, 0) = 0.1151116342572037;
    xi_result(2, 1) = 0.3985715991803625;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.839363025185972;
    eta_result(1) = 0.5611477338536844;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad8 surface, with given normals on the nodes.
   */
  void test_line_to_surface_intersection_normal_interpolation_quad8()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    LINALG::Matrix<24, 1, double> nodal_normals;
    xtest_setup_quad8(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.6839738851708264;
    xi_result(1, 1) = -0.3051161431281305;
    xi_result(2, 0) = 0.1455754614884382;
    xi_result(2, 1) = 0.5364371832797651;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.7799706383258106;
    eta_result(1) = 0.2729951612552455;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad9 surface, with default normals on the surface.
   */
  void test_line_to_surface_intersection_quad9()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    xtest_setup_quad9(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.7317907464850744;
    xi_result(1, 1) = -0.02799989440327506;
    xi_result(2, 0) = 0.1080035769948319;
    xi_result(2, 1) = 0.3188357119982439;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.865583933012948;
    eta_result(1) = 0.926806412303738;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad9 surface, with given normals on the nodes.
   */
  void test_line_to_surface_intersection_normal_interpolation_quad9()
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    xtest_setup_beam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    pair.is_unit_test_ = true;
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    LINALG::Matrix<27, 1, double> nodal_normals;
    xtest_setup_quad9(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    TS_ASSERT_EQUALS(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.6516378999140468;
    xi_result(1, 1) = -0.03862428489685134;
    xi_result(2, 0) = 0.111426072236278;
    xi_result(2, 1) = 0.33200167129208;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.869816485526844;
    eta_result(1) = 0.808011110533093;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      TS_ASSERT_DELTA(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        TS_ASSERT_DELTA(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

 private:
  //! Evaluation data container for geometry pairs.
  Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData> evaluation_data_;
};

#endif
