/*!

\brief Unit tests for line to volume geometry pairs with the segmentation algorithm.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_H_
#define UNIT_GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_H_


#include <cxxtest/TestSuite.h>

#include "src/drt_beam3/beam3r.H"
#include "src/drt_so3/so_hex8.H"
#include "src/drt_so3/so_hex27.H"
#include "src/drt_so3/so_tet10.H"
#include "src/drt_geometry_pair/geometry_pair_element_types.H"
#include "src/drt_geometry_pair/geometry_pair_evaluation_data_global.H"
#include "src/drt_geometry_pair/geometry_pair_line_to_volume_segmentation.H"
#include "src/drt_geometry_pair/geometry_pair_line_to_volume_evaluation_data.H"
#include "src/drt_geometry_pair/geometry_pair_utility_classes.H"
#include "src/drt_lib/drt_globalproblem.H"

#include "../drt_geometry_pair/unit_geometry_pair_line_to_volume_segmentation_geometry_functions.H"


namespace GEOMETRYPAIR
{
  class GeometryPairLineToVolumeSegmentation_TestSuite;
}

/**
 * Class to test the line to volume geometry pair segmentation algorithm.
 */
class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation_TestSuite : public CxxTest::TestSuite
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

    // Clear the vectors.
    line_elements_.clear();
    volume_elements_.clear();

    // We need to make sure the DRT::Problem instance created in setUp is deleted again. If this is
    // not done, some troubles arise where unit tests influence each other on some configurations.
    // We suspect that missing singleton destruction might be the reason for that.
    DRT::Problem::Done();
  }

  /**
   * Test a non straight beam that lies exactly between two solid. The segmentation should only be
   * done on one of the pairs.
   */
  void test_line_along_element_surface()
  {
    // Definition of variables for this test case.
    std::vector<LINALG::TMatrix<double, 12, 1>> q_line_elements;
    std::vector<LINALG::TMatrix<double, 9, 1>> q_rot_line_elements;
    std::vector<LINALG::TMatrix<double, 24, 1>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>>>
        geometry_pairs;

    // Get the geometry.
    xtest_line_along_element_surface_geometry(line_elements_, volume_elements_, q_line_elements,
        q_rot_line_elements, q_volume_elements, geometry_pairs);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // The segment is found on both pairs, but only evaluated on the first one that found it.
      TS_ASSERT_EQUALS(segments_vector[0].size(), 1);
      TS_ASSERT_EQUALS(segments_vector[1].size(), 0);

      // The first pair contains the full beam.
      TS_ASSERT_DELTA(
          -1., segments_vector[0][0].GetEtaA(), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(
          1., segments_vector[0][0].GetEtaB(), GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test a beam that has multiple intersections with a HEX27 element.
   */
  void test_multiple_intersections_hex27()
  {
    // Definition of variables for this test case.
    std::vector<LINALG::TMatrix<double, 12, 1>> q_line_elements;
    std::vector<LINALG::TMatrix<double, 9, 1>> q_rot_line_elements;
    std::vector<LINALG::TMatrix<double, 81, 1>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>>>
        geometry_pairs;

    // Get the geometry.
    xtest_multiple_intersections_hex27_geometry(line_elements_, volume_elements_, q_line_elements,
        q_rot_line_elements, q_volume_elements, geometry_pairs);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // Two segments should be found.
      TS_ASSERT_EQUALS(segments_vector[0].size(), 2);

      // Check the segment coordinates on the line.
      TS_ASSERT_DELTA(-0.748962623795456, segments_vector[0][0].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(-0.441636318022276, segments_vector[0][0].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(0.0775865904109607, segments_vector[0][1].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(0.951758620200254, segments_vector[0][1].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

  /**
   * Test a beam that has multiple intersections with a TET10 element.
   */
  void test_multiple_intersections_tet10()
  {
    // Definition of variables for this test case.
    std::vector<LINALG::TMatrix<double, 12, 1>> q_line_elements;
    std::vector<LINALG::TMatrix<double, 9, 1>> q_rot_line_elements;
    std::vector<LINALG::TMatrix<double, 30, 1>> q_volume_elements;
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
        GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>>>
        geometry_pairs;

    // Get the geometry.
    xtest_multiple_intersections_tet10_geometry(line_elements_, volume_elements_, q_line_elements,
        q_rot_line_elements, q_volume_elements, geometry_pairs);

    // Vector with vector of segments for Evaluate.
    std::vector<std::vector<GEOMETRYPAIR::LineSegment<double>>> segments_vector;

    // Create and evaluate the geometry pairs.
    CreateEvaluatePairs(
        geometry_pairs, q_line_elements, q_rot_line_elements, q_volume_elements, segments_vector);

    // Check results.
    {
      // Two segments should be found.
      TS_ASSERT_EQUALS(segments_vector[0].size(), 2);

      // Check the segment coordinates on the line.
      TS_ASSERT_DELTA(-0.412044388317299, segments_vector[0][0].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(-0.0831807800943276, segments_vector[0][0].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(0.879733788402666, segments_vector[0][1].GetEtaA(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
      TS_ASSERT_DELTA(0.978734285704557, segments_vector[0][1].GetEtaB(),
          GEOMETRYPAIR::CONSTANTS::projection_xi_eta_tol);
    }
  }

 private:
  /**
   * This function set up the geometry pairs, line elements and calls evaluate on the pairs. The
   * geometry of the test case is given in the various input parameters. The calculated segments are
   * returned which can be checked for results.
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
      const std::vector<LINALG::TMatrix<double, el1::n_dof_, 1>>& q_line_elements,
      const std::vector<LINALG::TMatrix<double, 9, 1>>& q_rot_line_elements,
      const std::vector<LINALG::TMatrix<double, el2::n_dof_, 1>>& q_volume_elements,
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
        for (unsigned int i = 0; i < 3; i++) xrefe[i + 3 * j] = q_line_elements[i_beam](i + j * 6);

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
        geometry_pairs.push_back(Teuchos::rcp(
            new GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, el1, el2>()));
        geometry_pairs.back()->Init(evaluation_data_, line.get(), volume.get());
        geometry_pairs.back()->Setup();
      }
    }
    segments_vector.resize(geometry_pairs.size());

    // Evaluate the segmentation.
    unsigned int counter = 0;
    for (unsigned int i_line = 0; i_line < q_line_elements.size(); i_line++)
    {
      for (unsigned int i_volume = 0; i_volume < q_volume_elements.size(); i_volume++)
      {
        std::vector<GEOMETRYPAIR::LineSegment<double>> segments;
        geometry_pairs[counter]->Evaluate(
            q_line_elements[i_line], q_volume_elements[i_volume], segments_vector[counter]);
        counter++;
      }
    }
  }

 private:
  //! Evaluation data container for geometry pairs.
  Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_;

  //! Vector of line elements.
  std::vector<Teuchos::RCP<DRT::Element>> line_elements_;

  //! Vector of volume elements.
  std::vector<Teuchos::RCP<DRT::Element>> volume_elements_;
};

#endif
