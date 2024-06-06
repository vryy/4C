/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to surface geometry pairs.

\level 1
*/
// End doxygen header.


#include <gtest/gtest.h>

#include "4C_geometry_pair_line_to_surface.hpp"

#include "4C_beam3_reissner.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_geometry_test.hpp"
#include "4C_geometry_pair_utility_classes.hpp"
#include "4C_inpar_beam_to_solid.hpp"
#include "4C_so3_surface.hpp"

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
      Inpar::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_surface_params_list);
      Inpar::GEOMETRYPAIR::SetValidParametersLineToSurface(line_to_surface_params_list);
      evaluation_data_ =
          Teuchos::rcp(new GEOMETRYPAIR::LineToSurfaceEvaluationData(line_to_surface_params_list));
    }

    /**
     * Set that the pair is a unit test pair. This has to be done here sine otherwise a gtest
     * specific macro has to be used to define the friend class.
     */
    template <typename A, typename B>
    void set_is_unit_test(
        GEOMETRYPAIR::GeometryPairLineToSurface<double, A, B>& pair, const int is_unit_test)
    {
      pair.is_unit_test_ = is_unit_test;
    }

    //! Evaluation data container for geometry pairs.
    Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData> evaluation_data_;
  };  // namespace

  /**
   * Test the projection of a point to a tri3 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationTri3)
  {
    // Set up the pair.
    Teuchos::RCP<Core::Elements::Element> beam = Teuchos::rcp(new Discret::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(beam.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates and normals for the solid element.
    const auto element_data_solid = XtestSetupTri3();

    // Point to project to.
    Core::LinAlg::Matrix<3, 1, double> point(true);
    point(0) = 0.3;
    point(1) = 0.1;
    point(2) = 0.2;

    // Project the point to the surface.
    Core::LinAlg::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, element_data_solid, xi, projection_result);

    // Check the results.
    Core::LinAlg::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3457692493957274;
    xi_result(1) = 0.2853120425437799;
    xi_result(2) = 0.03218342274405913;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a tri3 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationTri6)
  {
    // Set up the pair.
    Teuchos::RCP<Core::Elements::Element> beam = Teuchos::rcp(new Discret::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(beam.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates and normals for the solid element.
    const auto element_data_solid = XtestSetupTri6();

    // Point to project to.
    Core::LinAlg::Matrix<3, 1, double> point(true);
    point(0) = 0.3;
    point(1) = 0.1;
    point(2) = 0.2;

    // Project the point to the surface.
    Core::LinAlg::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, element_data_solid, xi, projection_result);

    // Check the results.
    Core::LinAlg::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3274411842809972;
    xi_result(1) = 0.1649919700896869;
    xi_result(2) = 0.2749865824042791;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad4 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationQuad4)
  {
    // Set up the pair.
    Teuchos::RCP<Core::Elements::Element> beam = Teuchos::rcp(new Discret::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(beam.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates and normals for the solid element.
    const auto element_data_solid = XtestSetupQuad4();

    // Point to project to.
    Core::LinAlg::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    Core::LinAlg::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, element_data_solid, xi, projection_result);

    // Check the results.
    Core::LinAlg::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.6306816217205055;
    xi_result(1) = -0.2391123963538002;
    xi_result(2) = 0.1168739495183324;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad8 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationQuad8)
  {
    // Set up the pair.
    Teuchos::RCP<Core::Elements::Element> beam = Teuchos::rcp(new Discret::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(beam.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates for the solid element.
    const auto element_data_solid = XtestSetupQuad8();

    // Point to project to.
    Core::LinAlg::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    Core::LinAlg::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, element_data_solid, xi, projection_result);

    // Check the results.
    Core::LinAlg::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = -0.167932271257968;
    xi_result(1) = 0.1593451990533972;
    xi_result(2) = 0.6729448863050194;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad9 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestPointToSurfaceProjectionNormalInterpolationQuad9)
  {
    // Set up the pair.
    Teuchos::RCP<Core::Elements::Element> beam = Teuchos::rcp(new Discret::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(beam.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates and normals for the solid element.
    const auto element_data_solid = XtestSetupQuad9();

    // Point to project to.
    Core::LinAlg::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    Core::LinAlg::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToOther(point, element_data_solid, xi, projection_result);

    // Check the results.
    Core::LinAlg::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3784195771508677;
    xi_result(1) = -0.436333510864013;
    xi_result(2) = 0.2483249147920992;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      EXPECT_NEAR(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
  }

  /**
   * Test the intersection of a line with a tri3 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationTri3)
  {
    // Set up the beam.
    const auto [element_line, element_data_line] = XtestSetupBeam();

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>
        pair(element_line.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates for the solid element.
    const auto element_data_solid = XtestSetupTri3();

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    Core::LinAlg::Matrix<3, 1, double> xi_start(true);
    pair.intersect_line_with_other(
        element_data_line, element_data_solid, intersection_points, 0., xi_start);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    Core::LinAlg::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.;
    xi_result(0, 1) = 0.5449417178907668;
    xi_result(1, 0) = 0.0892959494780295;
    xi_result(1, 1) = 0.4550582821092332;
    xi_result(2, 0) = 0.1071885022364286;
    xi_result(2, 1) = 0.00847913923973852;

    Core::LinAlg::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.933144945731849;
    eta_result(1) = -0.2769426498640102;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::Constants::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a tri6 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationTri6)
  {
    // Set up the beam.
    const auto [element_line, element_data_line] = XtestSetupBeam();

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>
        pair(element_line.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates for the solid element.
    const auto element_data_solid = XtestSetupTri6();

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    Core::LinAlg::Matrix<3, 1, double> xi_start(true);
    pair.intersect_line_with_other(
        element_data_line, element_data_solid, intersection_points, 0., xi_start);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    Core::LinAlg::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = 0.0;
    xi_result(0, 1) = 0.6584930664718799;
    xi_result(1, 0) = 0.1326645539292272;
    xi_result(1, 1) = 0.3415069335281202;
    xi_result(2, 0) = 0.1167653982105834;
    xi_result(2, 1) = 0.1176127305359125;

    Core::LinAlg::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.835011952021638;
    eta_result(1) = -0.1706930518170407;
    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::Constants::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad4 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationQuad4)
  {
    // Set up the beam.
    const auto [element_line, element_data_line] = XtestSetupBeam();

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair(element_line.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates for the solid element.
    const auto element_data_solid = XtestSetupQuad4();

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    Core::LinAlg::Matrix<3, 1, double> xi_start(true);
    pair.intersect_line_with_other(
        element_data_line, element_data_solid, intersection_points, 0., xi_start);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    Core::LinAlg::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.825477582092077;
    xi_result(1, 1) = -0.01147092212982951;
    xi_result(2, 0) = 0.1073378316553475;
    xi_result(2, 1) = 0.11957685229373;

    Core::LinAlg::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.932642561926223;
    eta_result(1) = 0.4203246988865186;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::Constants::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad8 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationQuad8)
  {
    // Set up the beam.
    const auto [element_line, element_data_line] = XtestSetupBeam();

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair(element_line.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);

    // Define the coordinates for the solid element.
    const auto element_data_solid = XtestSetupQuad8();

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    Core::LinAlg::Matrix<3, 1, double> xi_start(true);
    pair.intersect_line_with_other(
        element_data_line, element_data_solid, intersection_points, 0., xi_start);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    Core::LinAlg::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.684038067025361;
    xi_result(1, 1) = -0.3051565576006982;
    xi_result(2, 0) = 0.1455461215481048;
    xi_result(2, 1) = 0.5364549622511008;

    Core::LinAlg::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.7800927748816529;
    eta_result(1) = 0.273043123755964;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::Constants::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
    }
  }

  /**
   * Test the intersection of a line with a quad9 surface, with given normals on the nodes.
   */
  TEST_F(GeometryPairLineToSurfaceTest, TestLineToSurfaceIntersectionNormalInterpolationQuad9)
  {
    // Set up the beam.
    const auto [element_line, element_data_line] = XtestSetupBeam();

    // Set up the pair.
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair(element_line.get(), nullptr, evaluation_data_);
    set_is_unit_test(pair, true);
    // Define the coordinates for the solid element.
    const auto element_data_solid = XtestSetupQuad9();

    // Intersect the beam with the surface.
    std::vector<ProjectionPoint1DTo3D<double>> intersection_points;
    Core::LinAlg::Matrix<3, 1, double> xi_start(true);
    pair.intersect_line_with_other(
        element_data_line, element_data_solid, intersection_points, 0., xi_start);

    // Check the results.
    EXPECT_EQ(intersection_points.size(), 2);

    Core::LinAlg::Matrix<3, 2, double> xi_result;
    xi_result(0, 0) = -1.;
    xi_result(0, 1) = 1.;
    xi_result(1, 0) = -0.6516666236508463;
    xi_result(1, 1) = -0.0386216926927474;
    xi_result(2, 0) = 0.1114178909512789;
    xi_result(2, 1) = 0.3320135883392168;

    Core::LinAlg::Matrix<2, 1, double> eta_result;
    eta_result(0) = -0.869880338147344;
    eta_result(1) = 0.808094949905998;

    for (unsigned int i_intersection = 0; i_intersection < 2; i_intersection++)
    {
      EXPECT_NEAR(intersection_points[i_intersection].GetEta(), eta_result(i_intersection),
          GEOMETRYPAIR::Constants::projection_xi_eta_tol);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        EXPECT_NEAR(intersection_points[i_intersection].GetXi()(i_dir),
            xi_result(i_dir, i_intersection), GEOMETRYPAIR::Constants::projection_xi_eta_tol);
    }
  }

}  // namespace
