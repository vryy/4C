/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to surface geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_H_
#define UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_H_


#include <cxxtest/TestSuite.h>

#include "src/drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "src/drt_geometry_pair/geometry_pair_line_to_volume_evaluation_data.H"
#include "src/drt_geometry_pair/geometry_pair_element_functions.H"
#include "src/drt_geometry_pair/geometry_pair_utility_classes.H"
#include "src/drt_so3/so_surface.H"
#include "src/drt_beam3/beam3r.H"
#include "src/drt_lib/drt_globalproblem.H"


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
    // Add an empty parameter list, then the default values are set.
    DRT::Problem& problem = (*DRT::Problem::Instance());
    Teuchos::RCP<Teuchos::ParameterList> global_parameter =
        Teuchos::rcp(new Teuchos::ParameterList());
    problem.setParameterList(global_parameter);

    // Set up the evaluation data container for the geometry pairs.
    evaluation_data_ = Teuchos::rcp(new GEOMETRYPAIR::GeometryEvaluationDataGlobal());
    evaluation_data_->BuildLineToVolumeEvaluationData();
    evaluation_data_->LineToVolumeEvaluationData()->SetGaussRule(
        DRT::UTILS::GaussRule1D::intrule_line_6point);
  }

  /**
   * Delete pointers and other class variables.
   */
  void tearDown()
  {
    // Dereference the pointers.
    evaluation_data_ = Teuchos::null;

    // We need to make sure the DRT::Problem instance created in setUp is deleted again. If this is
    // not done, some troubles arise where unit tests influence each other on some configurations.
    // We suspect that missing singleton destruction might be the reason for that.
    DRT::Problem::Done();
  }

  /**
   * Test the projection of a point to a quad4 surface, with default normals on the surface.
   */
  void test_point_to_surface_projection_quad4()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair;
    pair.Init(evaluation_data_, beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid(true);
    q_solid(0) = 0;
    q_solid(1) = 0;
    q_solid(2) = 0;
    q_solid(3) = 1.0000000000000000000;
    q_solid(4) = -0.50000000000000000000;
    q_solid(5) = 0.50000000000000000000;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0;

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToSurface(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.5856297224156624;
    xi_result(1) = -0.2330351551569786;
    xi_result(2) = 0.1132886291998745;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad4 surface, with default given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_quad4()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>
        pair;
    pair.Init(evaluation_data_, beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<12, 1, double> q_solid(true);
    q_solid(0) = 0;
    q_solid(1) = 0;
    q_solid(2) = 0;
    q_solid(3) = 1.0000000000000000000;
    q_solid(4) = -0.50000000000000000000;
    q_solid(5) = 0.50000000000000000000;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0;

    // Define the normals on the solid element.
    LINALG::Matrix<12, 1, double> nodal_normals(true);
    nodal_normals(0) = -0.2627627396383057;
    nodal_normals(1) = 0.0814482045510598;
    nodal_normals(2) = 0.961416628019913;
    nodal_normals(3) = -0.696398235712595;
    nodal_normals(4) = -0.00250557448550925;
    nodal_normals(5) = 0.7176511822556153;
    nodal_normals(6) = -0.5381757744868267;
    nodal_normals(7) = 0.2523235485556319;
    nodal_normals(8) = 0.804176387740773;
    nodal_normals(9) = -0.5695583742587585;
    nodal_normals(10) = 0.3920519545501438;
    nodal_normals(11) = 0.7224254447658469;

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToSurface(point, q_solid, xi, projection_result, &nodal_normals);

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
        pair;
    pair.Init(evaluation_data_, beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0.;
    q_solid(12) = 0.7;
    q_solid(13) = -0.1;
    q_solid(14) = 0.;
    q_solid(15) = 1.5;
    q_solid(16) = 0.5;
    q_solid(17) = 0.4285714285714286;
    q_solid(18) = 0.6;
    q_solid(19) = 1.;
    q_solid(20) = 0.;
    q_solid(21) = 0.;
    q_solid(22) = 0.4;
    q_solid(23) = 0.;

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToSurface(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.4869140501387866;
    xi_result(1) = -0.6545313748232923;
    xi_result(2) = 0.4772682324027889;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad8 surface, with default given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_quad8()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>
        pair;
    pair.Init(evaluation_data_, beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<24, 1, double> q_solid;
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0.;
    q_solid(12) = 0.7;
    q_solid(13) = -0.1;
    q_solid(14) = 0.;
    q_solid(15) = 1.5;
    q_solid(16) = 0.5;
    q_solid(17) = 0.4285714285714286;
    q_solid(18) = 0.6;
    q_solid(19) = 1.;
    q_solid(20) = 0.;
    q_solid(21) = 0.;
    q_solid(22) = 0.4;
    q_solid(23) = 0.;

    // Define the normals on the solid element.
    LINALG::Matrix<24, 1, double> nodal_normals;
    nodal_normals(0) = 0.4647142046388283;
    nodal_normals(1) = 0.07560665247369116;
    nodal_normals(2) = 0.882226922117334;
    nodal_normals(3) = -0.836343726511398;
    nodal_normals(4) = 0.4199949186410029;
    nodal_normals(5) = 0.3523257575607623;
    nodal_normals(6) = -0.685419065478165;
    nodal_normals(7) = -0.2637809647101298;
    nodal_normals(8) = 0.6786901408858333;
    nodal_normals(9) = 0.084751268287718;
    nodal_normals(10) = 0.5628096661844073;
    nodal_normals(11) = 0.822230200231674;
    nodal_normals(12) = -0.4965263011922144;
    nodal_normals(13) = -0.01413350221776395;
    nodal_normals(14) = 0.867906605770136;
    nodal_normals(15) = -0.979140049557007;
    nodal_normals(16) = 0.05236426002737823;
    nodal_normals(17) = 0.1963230695188068;
    nodal_normals(18) = 0.2410702660375203;
    nodal_normals(19) = -0.349981608608912;
    nodal_normals(20) = 0.905206054149064;
    nodal_normals(21) = 0.5745894251085984;
    nodal_normals(22) = -0.2036034941085485;
    nodal_normals(23) = 0.792712185941506;

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToSurface(point, q_solid, xi, projection_result, &nodal_normals);

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
        pair;
    pair.Init(evaluation_data_, beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0.;
    q_solid(12) = 0.7;
    q_solid(13) = -0.1;
    q_solid(14) = 0.;
    q_solid(15) = 1.5;
    q_solid(16) = 0.5;
    q_solid(17) = 0.4285714285714286;
    q_solid(18) = 0.6;
    q_solid(19) = 1.;
    q_solid(20) = 0.;
    q_solid(21) = 0.;
    q_solid(22) = 0.4;
    q_solid(23) = 0.;
    q_solid(24) = 0.5;
    q_solid(25) = 0.5;
    q_solid(26) = 0.2;

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToSurface(point, q_solid, xi, projection_result, NULL);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.4374951399531939;
    xi_result(1) = -0.4006486973745378;
    xi_result(2) = 0.2412946023554158;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

  /**
   * Test the projection of a point to a quad9 surface, with default given normals on the nodes.
   */
  void test_point_to_surface_projection_normal_interpolation_quad9()
  {
    // Set up the pair.
    Teuchos::RCP<DRT::Element> beam = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>
        pair;
    pair.Init(evaluation_data_, beam.get(), NULL);
    pair.Setup();

    // Define the coordinates for the solid element.
    LINALG::Matrix<27, 1, double> q_solid;
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0.;
    q_solid(12) = 0.7;
    q_solid(13) = -0.1;
    q_solid(14) = 0.;
    q_solid(15) = 1.5;
    q_solid(16) = 0.5;
    q_solid(17) = 0.4285714285714286;
    q_solid(18) = 0.6;
    q_solid(19) = 1.;
    q_solid(20) = 0.;
    q_solid(21) = 0.;
    q_solid(22) = 0.4;
    q_solid(23) = 0.;
    q_solid(24) = 0.5;
    q_solid(25) = 0.5;
    q_solid(26) = 0.2;

    // Define the normals on the solid element.
    LINALG::Matrix<27, 1, double> nodal_normals;
    nodal_normals(0) = 0.4647142046388283;
    nodal_normals(1) = 0.07560665247369116;
    nodal_normals(2) = 0.882226922117334;
    nodal_normals(3) = -0.836343726511398;
    nodal_normals(4) = 0.4199949186410029;
    nodal_normals(5) = 0.3523257575607623;
    nodal_normals(6) = -0.685419065478165;
    nodal_normals(7) = -0.2637809647101298;
    nodal_normals(8) = 0.6786901408858333;
    nodal_normals(9) = 0.084751268287718;
    nodal_normals(10) = 0.5628096661844073;
    nodal_normals(11) = 0.822230200231674;
    nodal_normals(12) = -0.6102068773340782;
    nodal_normals(13) = -0.5672458156179746;
    nodal_normals(14) = 0.5530639669315768;
    nodal_normals(15) = -0.477142251285473;
    nodal_normals(16) = -0.03320684296047792;
    nodal_normals(17) = 0.878198484181582;
    nodal_normals(18) = 0.0853062914327283;
    nodal_normals(19) = 0.898922710199504;
    nodal_normals(20) = 0.4297217678097922;
    nodal_normals(21) = -0.3719814588140768;
    nodal_normals(22) = -0.3844673823348876;
    nodal_normals(23) = 0.844875509302472;
    nodal_normals(24) = 0.2831165159078166;
    nodal_normals(25) = 0.1617600702632027;
    nodal_normals(26) = 0.945345819310935;

    // Point to project to.
    LINALG::Matrix<3, 1, double> point(true);
    point(0) = 0.8;
    point(1) = 0.2;
    point(2) = 0.5;

    // Project the point to the surface.
    LINALG::Matrix<3, 1, double> xi(true);
    ProjectionResult projection_result;
    pair.ProjectPointToSurface(point, q_solid, xi, projection_result, &nodal_normals);

    // Check the results.
    LINALG::Matrix<3, 1, double> xi_result(true);
    xi_result(0) = 0.3784195771508677;
    xi_result(1) = -0.436333510864013;
    xi_result(2) = 0.2483249147920992;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      TS_ASSERT_DELTA(xi(i_dim), xi_result(i_dim), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
  }

 private:
  //! Evaluation data container for geometry pairs.
  Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_;
};

#endif
