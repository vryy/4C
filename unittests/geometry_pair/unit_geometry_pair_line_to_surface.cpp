/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to surface geometry pairs.

\level 1
*/
// End doxygen header.


#include <gtest/gtest.h>

#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_line_to_surface_evaluation_data.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_utility_classes.H"
#include "so3_surface.H"
#include "beam3_reissner.H"
#include "inpar_beam_to_solid.H"

#include "unit_geometry_pair_line_to_surface_geometry.H"

using namespace GEOMETRYPAIR;

namespace
{
  /**
   * Class to test the line to volume geometry pair segmentation algorithm.
   */
  class GeometryPairLineToSurfaceTest : public ::testing::Test
  {
   protected:
    /**
     * Set up the testing environment.
     */
    GeometryPairLineToSurfaceTest()
    {
      // Set up the evaluation data container for the geometry pairs.
      Teuchos::ParameterList line_to_surface_params_list;
      INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_surface_params_list);
      INPAR::GEOMETRYPAIR::SetValidParametersLineToSurface(line_to_surface_params_list);
      evaluation_data_ =
          Teuchos::rcp(new GEOMETRYPAIR::LineToSurfaceEvaluationData(line_to_surface_params_list));
    }

    /**
     * Set that the pair is a unit test pair. This has to be done here sine otherwise a gtest
     * specific macro has to be used to define the friend class.
     */
    template <typename A, typename B>
    void SetIsUnitTest(
        GEOMETRYPAIR::GeometryPairLineToSurface<double, A, B>& pair, const int is_unit_test)
    {
      pair.is_unit_test_ = is_unit_test;
    }

    //! Evaluation data container for geometry pairs.
    Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData> evaluation_data_;
  };  // namespace

  /**
   * Test the projection of a point to a tri3 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionTri3)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    XtestSetupTri3(q_solid);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri3 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationTri3)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    LINALG::Matrix<9, 1, double> nodal_normals(true);
    XtestSetupTri3(q_solid, &nodal_normals);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri6 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionTri6)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    XtestSetupTri6(q_solid);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri3 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationTri6)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    LINALG::Matrix<18, 1, double> nodal_normals(true);
    XtestSetupTri6(q_solid, &nodal_normals);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad4 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionQuad4)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the normals on the solid element.
    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid(true);
    XtestSetupQuad4(q_solid);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad4 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationQuad4)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<12, 1, double> q_solid(true);
    LINALG::Matrix<12, 1, double> nodal_normals(true);
    XtestSetupQuad4(q_solid, &nodal_normals);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad8 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionQuad8)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    XtestSetupQuad8(q_solid);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad8 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationQuad8)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    LINALG::Matrix<24, 1, double> nodal_normals;
    XtestSetupQuad8(q_solid, &nodal_normals);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad9 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionQuad9)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    XtestSetupQuad9(q_solid);

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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad9 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationQuad9)
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(beam.get(), NULL);
    pair.Setup();

    // Define the coordinates and normals for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    LINALG::Matrix<27, 1, double> nodal_normals;
    XtestSetupQuad9(q_solid, &nodal_normals);


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
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the intersection of a line with a tri3 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionTri3)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    XtestSetupTri3(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.0;
    xi_result(0, 1) = 0.5442036547500066;
    xi_result(1, 0) = 0.1074351908646558;
    xi_result(1, 1) = 0.4557963452499935;
    xi_result(2, 0) = 0.1140198081425712;
    xi_result(2, 1) = 0.00817486544751219;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.959558791949879;
    eta_result(1) = -0.2755154609844958;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri3 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationTri3)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<9, 1, double> q_solid;
    LINALG::Matrix<9, 1, double> nodal_normals;
    XtestSetupTri3(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.;
    xi_result(0, 1) = 0.5449417178907668;
    xi_result(1, 0) = 0.0892959494780295;
    xi_result(1, 1) = 0.4550582821092332;
    xi_result(2, 0) = 0.1071885022364286;
    xi_result(2, 1) = 0.00847913923973852;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.933144945731849;
    eta_result(1) = -0.2769426498640102;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri6 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionTri6)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    XtestSetupTri6(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.0;
    xi_result(0, 1) = 0.6613364114808006;
    xi_result(1, 0) = 0.1351609040799041;
    xi_result(1, 1) = 0.3386635885191994;
    xi_result(2, 0) = 0.1130261935455847;
    xi_result(2, 1) = 0.133359925421189;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.845578728729394;
    eta_result(1) = -0.1960488586704193;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri6 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationTri6)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<18, 1, double> q_solid;
    LINALG::Matrix<18, 1, double> nodal_normals;
    XtestSetupTri6(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.0;
    xi_result(0, 1) = 0.6584930664718799;
    xi_result(1, 0) = 0.1326645539292272;
    xi_result(1, 1) = 0.3415069335281202;
    xi_result(2, 0) = 0.1167653982105834;
    xi_result(2, 1) = 0.1176127305359125;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.835011952021638;
    eta_result(1) = -0.1706930518170407;
    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad4 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionQuad4)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid;
    XtestSetupQuad4(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.7859873421054778;
    xi_result(1, 1) = 0.01350316176774645;
    xi_result(2, 0) = 0.1131078893969968;
    xi_result(2, 1) = 0.1177634472685727;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.957101099360353;
    eta_result(1) = 0.4601648421155885;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad4 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationQuad4)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid;
    LINALG::Matrix<12, 1, double> nodal_normals;
    XtestSetupQuad4(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.825477582092077;
    xi_result(1, 1) = -0.01147092212982951;
    xi_result(2, 0) = 0.1073378316553475;
    xi_result(2, 1) = 0.11957685229373;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.932642561926223;
    eta_result(1) = 0.4203246988865186;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad8 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionQuad8)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    XtestSetupQuad8(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.7289250572942389;
    xi_result(1, 1) = -0.24019526509169;
    xi_result(2, 0) = 0.1150995008049619;
    xi_result(2, 1) = 0.3986046875693287;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.839447041376278;
    eta_result(1) = 0.5612360819815785;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad8 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationQuad8)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    LINALG::Matrix<24, 1, double> nodal_normals;
    XtestSetupQuad8(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.684038067025361;
    xi_result(1, 1) = -0.3051565576006982;
    xi_result(2, 0) = 0.1455461215481048;
    xi_result(2, 1) = 0.5364549622511008;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.7800927748816529;
    eta_result(1) = 0.273043123755964;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad9 surface, with default normals on the surface.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionQuad9)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    XtestSetupQuad9(q_solid);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, NULL);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.7318086192269153;
    xi_result(1, 1) = -0.02799916735239928;
    xi_result(2, 0) = 0.107995627440907;
    xi_result(2, 1) = 0.3188379218922715;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.865652166867077;
    eta_result(1) = 0.92684679887125;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad9 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationQuad9)
  {
    // Set up the beam.
    Teuchos::RCP<DRT::Element> element_1;
    LINALG::Matrix<12, 1, double> q_beam;
    XtestSetupBeam(element_1, q_beam);

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(evaluation_data_);
    SetIsUnitTest(pair, true);
    pair.Init(element_1.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    LINALG::Matrix<27, 1, double> nodal_normals;
    XtestSetupQuad9(q_solid, &nodal_normals);

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    LINALG::Matrix<3, 1, double> xi_start(true);
    pair.IntersectLineWithOther(q_beam, q_solid, intersection_points, 0., xi_start, &nodal_normals);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    LINALG::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.6516666236508463;
    xi_result(1, 1) = -0.0386216926927474;
    xi_result(2, 0) = 0.1114178909512789;
    xi_result(2, 1) = 0.3320135883392168;

    LINALG::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.869880338147344;
    eta_result(1) = 0.808094949905998;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

}  // namespace
