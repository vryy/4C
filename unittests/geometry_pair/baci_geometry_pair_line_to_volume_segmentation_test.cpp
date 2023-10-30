/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to volume geometry pairs with the segmentation algorithm.

\level 1
*/
// End doxygen header.


#include <gtest/gtest.h>

#include "baci_geometry_pair_line_to_volume_segmentation.H"

#include "baci_beam3_reissner.H"
#include "baci_geometry_pair_element_functions.H"
#include "baci_geometry_pair_line_to_3D_evaluation_data.H"
#include "baci_geometry_pair_line_to_volume_segmentation_geometry_functions_test.H"
#include "baci_geometry_pair_utility_classes.H"
#include "baci_so3_hex27.H"


namespace
{
  /**
   * Class to test the line to volume geometry pair segmentation algorithm.
   */
  class GeometryPairLineToVolumeSegmentationTest : public ::testing::Test
  {
   protected:
    /**
     * Set up the testing environment.
     */
    GeometryPairLineToVolumeSegmentationTest()
    {
      // Set up the evaluation data container for the geometry pairs.
      Teuchos::ParameterList line_to_volume_params_list;
      INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_volume_params_list);
      evaluation_data_ =
          Teuchos::rcp(new GEOMETRYPAIR::LineTo3DEvaluationData(line_to_volume_params_list));
    }

    /**
     * This function set up the geometry pairs, line elements and calls evaluate on the pairs. The
     * geometry of the test case is given in the various input parameters. The calculated segments
     * are returned which can be checked for results.
     * @param geometry_pairs (out) Vector with the geometry pairs.
     * @param q_line_elements (in) Vector of DOF vectors for the positions and tangents of the line
     * element.
     * @param q_rot_line_elements (in) Vector with rotation vectors of the line nodes, in following
     * order: 0, 2, 1
     * @param q_volume_elements (in) Vector of DOF vectors for the positions of the volume element.
     * @param segments_vector (out) Vector with found segments for each pair.
     */
    template <typename el1, typename el2>
    void CreateEvaluatePairs(
        std::vector<Teuchos::RCP<
            GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, el1, el2>>>& geometry_pairs,
        const std::vector<CORE::LINALG::Matrix<el1::n_dof_, 1, double>>& q_line_elements,
        const std::vector<CORE::LINALG::Matrix<9, 1, double>>& q_rot_line_elements,
        const std::vector<CORE::LINALG::Matrix<el2::n_dof_, 1, double>>& q_volume_elements,
        std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>>& segments_vector)
    {
      // Check that the vectors have the right size.
      if (line_elements_.size() != q_line_elements.size())
        dserror("Size for line elements and line q does not match!");
      if (volume_elements_.size() != q_volume_elements.size())
        dserror("Size for volume elements and volume q does not match!");

      // Set the reference length for the beam.
      for (unsigned int i_beam = 0; i_beam < line_elements_.size(); i_beam++)
      {
        // Split up the rotational and positional DOF, needed to calculate the reference length.
        std::vector<double> xrefe(6);
        for (unsigned int j = 0; j < 2; j++)
          for (unsigned int i = 0; i < 3; i++)
            xrefe[i + 3 * j] = q_line_elements[i_beam](i + j * 6);

        // Get the rotational vector.
        std::vector<double> rotrefe(9);
        for (unsigned int i = 0; i < 9; i++) rotrefe[i] = q_rot_line_elements[i_beam](i);

        // Cast beam element.
        Teuchos::RCP<DRT::ELEMENTS::Beam3r> beam_element =
            Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Beam3r>(line_elements_[i_beam], true);

        // Set the hermitian interpolation.
        beam_element->SetCenterlineHermite(true);

        // Calculate the reference length.
        beam_element->SetUpReferenceGeometry<3, 2, 2>(xrefe, rotrefe);
      }

      // Create the geometry pairs.
      for (auto& line : line_elements_)
      {
        // Loop over each solid with this beam and create a pair.
        for (auto& volume : volume_elements_)
        {
          geometry_pairs.push_back(
              Teuchos::rcp(new GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, el1, el2>(
                  line.get(), volume.get(), evaluation_data_)));
        }
      }
      segments_vector.resize(geometry_pairs.size());

      // Evaluate the segmentation.
      unsigned int counter = 0;
      for (unsigned int i_line = 0; i_line < q_line_elements.size(); i_line++)
      {
        for (unsigned int i_volume = 0; i_volume < q_volume_elements.size(); i_volume++)
        {
          geometry_pairs[counter]->Evaluate(
              q_line_elements[i_line], q_volume_elements[i_volume], segments_vector[counter]);
          counter++;
        }
      }
    }

    //! Evaluation data container for geometry pairs.
    Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> evaluation_data_;

    //! Vector of line elements.
    std::vector<Teuchos::RCP<DRT::Element>> line_elements_;

    //! Vector of volume elements.
    std::vector<Teuchos::RCP<DRT::Element>> volume_elements_;
  };

  /**
   * Test a non straight beam that lies exactly between two solid. The segmentation should only be
   * done on one of the pairs.
   */
  TEST_F(GeometryPairLineToVolumeSegmentationTest, TestLineAlongElementSurface)
  {
    // Definition of variables for this test case.
    std::vector<CORE::LINALG::Matrix<12, 1, double>> q_line_elements;
    std::vector<CORE::LINALG::Matrix<9, 1, double>> q_rot_line_elements;
    std::vector<CORE::LINALG::Matrix<24, 1, double>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>>>
        geometry_pairs;

    // Get the geometry.
    XtestLineAlongElementSurfaceGeometry(
        line_elements_, volume_elements_, q_line_elements, q_rot_line_elements, q_volume_elements);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // The segment is found on both pairs, but only evaluated on the first one that found it.
      EXPECT_EQ(segments_vector[0].size(), 1);
      EXPECT_EQ(segments_vector[1].size(), 0);

      // The first pair contains the full beam.
      EXPECT_NEAR(
          -1., segments_vector[0][0].GetEtaA(), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(
          1., segments_vector[0][0].GetEtaB(), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test a non straight beam that lies between two volume elements. The elements are relatively
   * small. This is to check that the parameter coordinates are converged in the local Newton
   * iterations.
   */
  TEST_F(GeometryPairLineToVolumeSegmentationTest, TestLineInSmallElements)
  {
    // Definition of variables for this test case.
    std::vector<CORE::LINALG::Matrix<12, 1, double>> q_line_elements;
    std::vector<CORE::LINALG::Matrix<9, 1, double>> q_rot_line_elements;
    std::vector<CORE::LINALG::Matrix<24, 1, double>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>>>
        geometry_pairs;

    // Get the geometry.
    XtestLineInSmallElementsGeometry(
        line_elements_, volume_elements_, q_line_elements, q_rot_line_elements, q_volume_elements);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      EXPECT_EQ(segments_vector[0].size(), 1);
      EXPECT_EQ(segments_vector[1].size(), 1);

      // Check the parameter coordinates.
      EXPECT_NEAR(
          -1., segments_vector[0][0].GetEtaA(), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(0.36285977578126655, segments_vector[0][0].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(0.36285977578126744, segments_vector[1][0].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(
          1.0, segments_vector[1][0].GetEtaB(), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test a beam that has multiple intersections with a HEX27 element.
   */
  TEST_F(GeometryPairLineToVolumeSegmentationTest, TestMultipleIntersectionsHex27)
  {
    // Definition of variables for this test case.
    std::vector<CORE::LINALG::Matrix<12, 1, double>> q_line_elements;
    std::vector<CORE::LINALG::Matrix<9, 1, double>> q_rot_line_elements;
    std::vector<CORE::LINALG::Matrix<81, 1, double>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>>>
        geometry_pairs;

    // Get the geometry.
    XtestMultipleIntersectionsHex27Geometry(
        line_elements_, volume_elements_, q_line_elements, q_rot_line_elements, q_volume_elements);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // Two segments should be found.
      EXPECT_EQ(segments_vector[0].size(), 2);

      // Check the segment coordinates on the line.
      EXPECT_NEAR(-0.7495456134309243, segments_vector[0][0].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(-0.44451080329628256, segments_vector[0][0].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(0.076870238957896297, segments_vector[0][1].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(0.95200105410689462, segments_vector[0][1].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test a beam that has multiple intersections with a TET10 element.
   */
  TEST_F(GeometryPairLineToVolumeSegmentationTest, TestMultipleIntersectionsTet10)
  {
    // Definition of variables for this test case.
    std::vector<CORE::LINALG::Matrix<12, 1, double>> q_line_elements;
    std::vector<CORE::LINALG::Matrix<9, 1, double>> q_rot_line_elements;
    std::vector<CORE::LINALG::Matrix<30, 1, double>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>>>
        geometry_pairs;

    // Get the geometry.
    XtestMultipleIntersectionsTet10Geometry(
        line_elements_, volume_elements_, q_line_elements, q_rot_line_elements, q_volume_elements);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // Two segments should be found.
      EXPECT_EQ(segments_vector[0].size(), 2);

      // Check the segment coordinates on the line.
      EXPECT_NEAR(-0.40853230756138476, segments_vector[0][0].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(-0.079518054716933712, segments_vector[0][0].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(0.88282786408473413, segments_vector[0][1].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(0.97917616983415101, segments_vector[0][1].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }


  /**
   * Test a beam that has multiple intersections with a NURBS27 element.
   */
  TEST_F(GeometryPairLineToVolumeSegmentationTest, TestMultipleIntersectionsNurbs27)
  {
    // Definition of variables for this test case.
    std::vector<CORE::LINALG::Matrix<12, 1, double>> q_line_elements;
    std::vector<CORE::LINALG::Matrix<9, 1, double>> q_rot_line_elements;
    std::vector<CORE::LINALG::Matrix<81, 1, double>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27>>>
        geometry_pairs;

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Add the relevant nurbs information to the discretization.
    Teuchos::RCP<DRT::NURBS::NurbsDiscretization> structdis =
        Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", Teuchos::null));
    DRT::Problem::Instance()->AddDis("structure", structdis);

    // Get the geometry.
    XtestMultipleIntersectionsNurbs27Geometry(line_elements_, volume_elements_, q_line_elements,
        q_rot_line_elements, q_volume_elements, structdis);

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // Two segments should be found.
      EXPECT_EQ(segments_vector[0].size(), 1);

      // Check the segment coordinates on the line.
      EXPECT_NEAR(-0.40769440465702655, segments_vector[0][0].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      EXPECT_NEAR(-0.0090552523537554153, segments_vector[0][0].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

}  // namespace
