/*!
\file unit_beaminteraction_beam_to_solid_volume_meshtying_mortar.H

\brief Unit tests for the beam to volume meshtying mortar pair.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_BEAMINTERACTION_BEAM_TO_SOLIV_VOLUME_MESHTYING_MORTAR_H_
#define UNIT_BEAMINTERACTION_BEAM_TO_SOLIV_VOLUME_MESHTYING_MORTAR_H_


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
#include "src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_mortar.H"
#include "src/drt_lib/drt_globalproblem.H"
#include "src/linalg/linalg_serialdensevector.H"


#include "unit_beaminteraction_beam_to_solid_volume_meshtying_mortar_geometry.H"


namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingPairMortar_TestSuite;
}

/**
 * Class to test the local mortar matrices calculated by the beam to volume meshtying mortar pair.
 */
class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * \brief Set up the testing environment.
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
   * \brief Delete pointers and other class variables.
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
   * \brief Test a non straight beam in a hex8 element, with line2 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex8_line2()
  {
    // Definition of variables for this test case.
    std::vector<LINALG::TMatrix<double, 12, 1>> q_line_elements;
    std::vector<LINALG::TMatrix<double, 9, 1>> q_rot_line_elements;
    std::vector<LINALG::TMatrix<double, 24, 1>> q_volume_elements;

    // Get the geometry.
    xtest_beam_in_hex8_geometry(
        line_elements_, volume_elements_, q_line_elements, q_rot_line_elements, q_volume_elements);

    // Set up the beam reference geometry.
    SetUpReferenceGeometry(q_line_elements, q_rot_line_elements);

    // Vector with the two element for this pair.
    std::vector<const DRT::Element*> pair_elements;
    pair_elements.push_back(&(*line_elements_[0]));
    pair_elements.push_back(&(*volume_elements_[0]));

    // Create the meshtying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>
        beam_contact_pair =
            BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>();
    beam_contact_pair.Init(Teuchos::null, evaluation_data_, pair_elements);

    // Matrices that will be filled in local evaluation.
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Evaluate the local matrices.
    beam_contact_pair.CastGeometryPair()->Setup();
    beam_contact_pair.ele1posref_ = q_line_elements[0];
    beam_contact_pair.ele2posref_ = q_volume_elements[0];
    beam_contact_pair.EvaluateDM(local_D, local_M, local_kappa);

    // Check the results for kappa.
    TS_ASSERT_EQUALS(local_kappa.M(), GEOMETRYPAIR::t_line2::n_dof_);
    const double result_kappa[] = {0.3197459060512249, 0.3197459060512249, 0.3197459060512249,
        0.3333977460013357, 0.3333977460013357, 0.3333977460013357};
    for (int i_row = 0; i_row < local_kappa.M(); i_row++)
      TS_ASSERT_DELTA(local_kappa(i_row), result_kappa[i_row], 1e-14);

    // Check the results for D.
    TS_ASSERT_EQUALS(local_D.M(), GEOMETRYPAIR::t_line2::n_dof_);
    TS_ASSERT_EQUALS(local_D.N(), GEOMETRYPAIR::t_hermite::n_dof_);
    const double result_D[][12] = {{0.2190511901271037, 0, 0, 0.0206808975254317, 0, 0,
                                       0.1006947159241212, 0, 0, -0.01441371675994627, 0, 0},
        {0, 0.2190511901271037, 0, 0, 0.0206808975254317, 0, 0, 0.1006947159241212, 0, 0,
            -0.01441371675994627, 0},
        {0, 0, 0.2190511901271037, 0, 0, 0.0206808975254317, 0, 0, 0.1006947159241212, 0, 0,
            -0.01441371675994627},
        {0.0987608135029245, 0, 0, 0.01441371675994627, 0, 0, 0.2346369324984112, 0, 0,
            -0.02193259846463873, 0, 0},
        {0, 0.0987608135029245, 0, 0, 0.01441371675994627, 0, 0, 0.2346369324984112, 0, 0,
            -0.02193259846463873, 0},
        {0, 0, 0.0987608135029245, 0, 0, 0.01441371675994627, 0, 0, 0.2346369324984112, 0, 0,
            -0.02193259846463873}};
    for (int i_row = 0; i_row < local_D.M(); i_row++)
      for (int i_col = 0; i_col < local_D.N(); i_col++)
        TS_ASSERT_DELTA(local_D(i_row, i_col), result_D[i_row][i_col], 1e-14);

    // Check the results for M.
    TS_ASSERT_EQUALS(local_M.M(), GEOMETRYPAIR::t_line2::n_dof_);
    TS_ASSERT_EQUALS(local_M.N(), GEOMETRYPAIR::t_hex8::n_dof_);
    const double result_M[][24] = {
        {0.06017774373419256, 0, 0, 0.05408827084348954, 0, 0, 0.06839405790991376, 0, 0,
            0.07628556578127677, 0, 0, 0.01344127602858902, 0, 0, 0.01252108498227593, 0, 0,
            0.01678164537837932, 0, 0, 0.01805626139310799, 0, 0},
        {0, 0.06017774373419256, 0, 0, 0.05408827084348954, 0, 0, 0.06839405790991376, 0, 0,
            0.07628556578127677, 0, 0, 0.01344127602858902, 0, 0, 0.01252108498227593, 0, 0,
            0.01678164537837932, 0, 0, 0.01805626139310799, 0},
        {0, 0, 0.06017774373419256, 0, 0, 0.05408827084348954, 0, 0, 0.06839405790991376, 0, 0,
            0.07628556578127677, 0, 0, 0.01344127602858902, 0, 0, 0.01252108498227593, 0, 0,
            0.01678164537837932, 0, 0, 0.01805626139310799},
        {0.04493081627696245, 0, 0, 0.04298226440679511, 0, 0, 0.0555958294619839, 0, 0,
            0.05837011361724712, 0, 0, 0.02938122210054368, 0, 0, 0.0293593399569368, 0, 0,
            0.03630140942759987, 0, 0, 0.03647675075326669, 0, 0},
        {0, 0.04493081627696245, 0, 0, 0.04298226440679511, 0, 0, 0.0555958294619839, 0, 0,
            0.05837011361724712, 0, 0, 0.02938122210054368, 0, 0, 0.0293593399569368, 0, 0,
            0.03630140942759987, 0, 0, 0.03647675075326669, 0},
        {0, 0, 0.04493081627696245, 0, 0, 0.04298226440679511, 0, 0, 0.0555958294619839, 0, 0,
            0.05837011361724712, 0, 0, 0.02938122210054368, 0, 0, 0.0293593399569368, 0, 0,
            0.03630140942759987, 0, 0, 0.03647675075326669}};
    for (int i_row = 0; i_row < local_M.M(); i_row++)
      for (int i_col = 0; i_col < local_M.N(); i_col++)
        TS_ASSERT_DELTA(local_M(i_row, i_col), result_M[i_row][i_col], 1e-14);
  }

  /**
   * \brief Test a non straight beam in a hex8 element, with line3 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex8_line3()
  {
    // Definition of variables for this test case.
    std::vector<LINALG::TMatrix<double, 12, 1>> q_line_elements;
    std::vector<LINALG::TMatrix<double, 9, 1>> q_rot_line_elements;
    std::vector<LINALG::TMatrix<double, 24, 1>> q_volume_elements;

    // Get the geometry.
    xtest_beam_in_hex8_geometry(
        line_elements_, volume_elements_, q_line_elements, q_rot_line_elements, q_volume_elements);

    // Set up the beam reference geometry.
    SetUpReferenceGeometry(q_line_elements, q_rot_line_elements);

    // Vector with the two element for this pair.
    std::vector<const DRT::Element*> pair_elements;
    pair_elements.push_back(&(*line_elements_[0]));
    pair_elements.push_back(&(*volume_elements_[0]));

    // Create the meshtying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>
        beam_contact_pair =
            BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>();
    beam_contact_pair.Init(Teuchos::null, evaluation_data_, pair_elements);

    // Matrices that will be filled in local evaluation.
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Evaluate the local matrices.
    beam_contact_pair.CastGeometryPair()->Setup();
    beam_contact_pair.ele1posref_ = q_line_elements[0];
    beam_contact_pair.ele2posref_ = q_volume_elements[0];
    beam_contact_pair.EvaluateDM(local_D, local_M, local_kappa);

    // Check the results for kappa.
    TS_ASSERT_EQUALS(local_kappa.M(), GEOMETRYPAIR::t_line3::n_dof_);
    const double result_kappa[] = {0.0989906405971891, 0.0989906405971891, 0.0989906405971891,
        0.11264248054729985, 0.11264248054729985, 0.11264248054729985, 0.4415105309080716,
        0.4415105309080716, 0.4415105309080716};
    for (int i_row = 0; i_row < local_kappa.M(); i_row++)
      TS_ASSERT_DELTA(local_kappa(i_row), result_kappa[i_row], 1e-14);

    // Check the results for D.
    TS_ASSERT_EQUALS(local_D.M(), GEOMETRYPAIR::t_line3::n_dof_);
    TS_ASSERT_EQUALS(local_D.N(), GEOMETRYPAIR::t_hermite::n_dof_);
    const double result_D[][12] = {{0.1112459262865433, 0, 0, 0.006473801602571051, 0, 0,
                                       -0.01225528568935419, 0, 0, 0.0002066208370856273, 0, 0},
        {0, 0.1112459262865433, 0, 0, 0.006473801602571051, 0, 0, -0.01225528568935419, 0, 0,
            0.0002066208370856273, 0},
        {0, 0, 0.1112459262865433, 0, 0, 0.006473801602571051, 0, 0, -0.01225528568935419, 0, 0,
            0.0002066208370856273},
        {-0.009044450337635913, 0, 0, 0.0002066208370856275, 0, 0, 0.1216869308849358, 0, 0,
            -0.007312260867606827, 0, 0},
        {0, -0.009044450337635913, 0, 0, 0.0002066208370856275, 0, 0, 0.1216869308849358, 0, 0,
            -0.007312260867606827, 0},
        {0, 0, -0.009044450337635913, 0, 0, 0.0002066208370856275, 0, 0, 0.1216869308849358, 0, 0,
            -0.007312260867606827},
        {0.2156105276811209, 0, 0, 0.02841419184572129, 0, 0, 0.2259000032269507, 0, 0,
            -0.0292406751940638, 0, 0},
        {0, 0.2156105276811209, 0, 0, 0.02841419184572129, 0, 0, 0.2259000032269507, 0, 0,
            -0.0292406751940638, 0},
        {0, 0, 0.2156105276811209, 0, 0, 0.02841419184572129, 0, 0, 0.2259000032269507, 0, 0,
            -0.0292406751940638}};
    for (int i_row = 0; i_row < local_D.M(); i_row++)
      for (int i_col = 0; i_col < local_D.N(); i_col++)
        TS_ASSERT_DELTA(local_D(i_row, i_col), result_D[i_row][i_col], 1e-14);

    // Check the results for M.
    TS_ASSERT_EQUALS(local_M.M(), GEOMETRYPAIR::t_line3::n_dof_);
    TS_ASSERT_EQUALS(local_M.N(), GEOMETRYPAIR::t_hex8::n_dof_);
    const double result_M[][24] = {
        {0.02523442610928899, 0, 0, 0.02219440099004521, 0, 0, 0.02563540648028574, 0, 0,
            0.02937626210868358, 0, 0, -0.0006039567482481712, 0, 0, -0.0009567576673671217, 0, 0,
            -0.001170602343818399, 0, 0, -0.0007185383316807301, 0, 0},
        {0, 0.02523442610928899, 0, 0, 0.02219440099004521, 0, 0, 0.02563540648028574, 0, 0,
            0.02937626210868358, 0, 0, -0.0006039567482481712, 0, 0, -0.0009567576673671217, 0, 0,
            -0.001170602343818399, 0, 0, -0.0007185383316807301, 0},
        {0, 0, 0.02523442610928899, 0, 0, 0.02219440099004521, 0, 0, 0.02563540648028574, 0, 0,
            0.02937626210868358, 0, 0, -0.0006039567482481712, 0, 0, -0.0009567576673671217, 0, 0,
            -0.001170602343818399, 0, 0, -0.0007185383316807301},
        {0.009987498652058875, 0, 0, 0.01108839455335078, 0, 0, 0.01283717803235589, 0, 0,
            0.01146080994465393, 0, 0, 0.01533598932370649, 0, 0, 0.01588149730729375, 0, 0,
            0.01834916170540214, 0, 0, 0.01770195102847798, 0, 0},
        {0, 0.009987498652058875, 0, 0, 0.01108839455335078, 0, 0, 0.01283717803235589, 0, 0,
            0.01146080994465393, 0, 0, 0.01533598932370649, 0, 0, 0.01588149730729375, 0, 0,
            0.01834916170540214, 0, 0, 0.01770195102847798, 0},
        {0, 0, 0.009987498652058875, 0, 0, 0.01108839455335078, 0, 0, 0.01283717803235589, 0, 0,
            0.01146080994465393, 0, 0, 0.01533598932370649, 0, 0, 0.01588149730729375, 0, 0,
            0.01834916170540214, 0, 0, 0.01770195102847798},
        {0.06988663524980716, 0, 0, 0.06378773970688865, 0, 0, 0.08551730285925604, 0, 0,
            0.09381860734518636, 0, 0, 0.02809046555367439, 0, 0, 0.0269556852992861, 0, 0,
            0.03590449544439545, 0, 0, 0.03754959944957743, 0, 0},
        {0, 0.06988663524980716, 0, 0, 0.06378773970688865, 0, 0, 0.08551730285925604, 0, 0,
            0.09381860734518636, 0, 0, 0.02809046555367439, 0, 0, 0.0269556852992861, 0, 0,
            0.03590449544439545, 0, 0, 0.03754959944957743, 0},
        {0, 0, 0.06988663524980716, 0, 0, 0.06378773970688865, 0, 0, 0.08551730285925604, 0, 0,
            0.09381860734518636, 0, 0, 0.02809046555367439, 0, 0, 0.0269556852992861, 0, 0,
            0.03590449544439545, 0, 0, 0.03754959944957743}};
    for (int i_row = 0; i_row < local_M.M(); i_row++)
      for (int i_col = 0; i_col < local_M.N(); i_col++)
        TS_ASSERT_DELTA(local_M(i_row, i_col), result_M[i_row][i_col], 1e-14);
  }

 private:
  /**
   * \brief This function sets up the beam elements in a way where the reference length can be
   * calculated (is needed for shale function evaluation).
   * @param q_line_elements (in) Vector of DOF vectors for the positions and tangents of the line
   * element.
   * @param q_rot_line_elements (in) Vector with rotation vectors of the line nodes, in following
   * order: 0, 2, 1
   */
  template <unsigned int n_dof_el1>
  void SetUpReferenceGeometry(
      const std::vector<LINALG::TMatrix<double, n_dof_el1, 1>>& q_line_elements,
      const std::vector<LINALG::TMatrix<double, 9, 1>>& q_rot_line_elements)
  {
    // Check that the vector has the right size.
    if (line_elements_.size() != q_line_elements.size())
      dserror("Size for line elements and line q does not match!");

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
