/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the beam to fluid meshtying gpts pair.

\level 3
\maintainer Nora Hagmeyer
*/
// End doxygen header.


#ifndef UNIT_BEAM_TO_FLUID_MESHTYING_PAIR_GAUSS_POINT_H
#define UNIT_BEAM_TO_FLUID_MESHTYING_PAIR_GAUSS_POINT_H


#include <cxxtest/TestSuite.h>
#include <vector>

#include "../../../src/drt_beam3/beam3eb.H"
#include "../../../src/drt_fluid_ele/fluid_ele.H"
#include "../../../src/drt_geometry_pair/geometry_pair_element_functions.H"
#include "../../../src/drt_geometry_pair/geometry_pair_line_to_volume_segmentation.H"
#include "../../../src/drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../../../src/drt_geometry_pair/geometry_pair_utility_classes.H"
#include "../../../src/drt_fbi/beam_to_fluid_meshtying_pair_gauss_point.H"
#include "../../../src/drt_lib/drt_globalproblem.H"
#include "../../../src/drt_lib/drt_element.H"
#include "../../../src/linalg/linalg_serialdensevector.H"
#include "../../../src/drt_beaminteraction/beam_contact_pair.H"
#include "../../../src/drt_fbi/beam_to_fluid_meshtying_params.H"


namespace BEAMINTERACTION
{
  class BeamToFluidMeshtyingPairGPTS_TestSuite;
}

/**
 * Class to test the local mortar matrices calculated by the beam to volume meshtying mortar pair.
 */
class BEAMINTERACTION::BeamToFluidMeshtyingPairGPTS_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * \brief Set up the testing environment.
   */
  void setUp()
  {
    // Set up the evaluation data container for the geometry pairs.
    Teuchos::ParameterList line_to_volume_params_list;
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_volume_params_list);
    evaluation_data_ =
        Teuchos::rcp(new GEOMETRYPAIR::LineTo3DEvaluationData(line_to_volume_params_list));
  }

  /**
   * \brief Delete pointers and other class variables.
   */
  void tearDown()
  {
    // Dereference the pointers.
    evaluation_data_ = Teuchos::null;
  }

  /**
   * \brief Test a moving straight beam in a hex8 element with hermite line2  shape functions.
   */
  void test_beam_to_fluid_meshtying_hex8_moving_beam()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex8 fluid_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam_type, fluid_type> pair =
        BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam_type, fluid_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<beam_type::n_dof_, 1, double> v_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<fluid_type::n_dof_, 1, double> q_fluid;
    LINALG::Matrix<fluid_type::n_dof_, 1, double> v_fluid;
    std::vector<double> beam_centerline_dofvec;
    std::vector<double> fluid_dofvec;

    // Matrices for the results.
    LINALG::Matrix<fluid_type::n_dof_, fluid_type::n_dof_, double> results_kff(true);
    LINALG::Matrix<fluid_type::n_dof_, beam_type::n_dof_, double> results_kfs(true);
    LINALG::Matrix<beam_type::n_dof_, fluid_type::n_dof_, double> results_ksf(true);
    LINALG::SerialDenseVector results_fs(beam_type::n_dof_, true);
    LINALG::SerialDenseVector results_ff(fluid_type::n_dof_, true);
    results_fs.Zero();
    results_ff.Zero();


    // Define the geometry of the two elements.
    q_beam(0) = 0.5;
    q_beam(1) = -1.0;
    q_beam(2) = 0.5;
    q_beam(3) = 0.0;
    q_beam(4) = 1.0;
    q_beam(5) = 0.0;
    q_beam(6) = 0.5;
    q_beam(7) = 1.0;
    q_beam(8) = 0.5;
    q_beam(9) = 0.0;
    q_beam(10) = 1.0;
    q_beam(11) = 0.0;

    v_beam(0) = 1.0;
    v_beam(1) = 0.0;
    v_beam(2) = 0.0;
    v_beam(3) = 0.0;
    v_beam(4) = 0.0;
    v_beam(5) = 0.0;
    v_beam(6) = 1.0;
    v_beam(7) = 0.0;
    v_beam(8) = 0.0;
    v_beam(9) = 0.0;
    v_beam(10) = 0.0;
    v_beam(11) = 0.0;

    for (unsigned int i = 0; i < beam_type::n_dof_; i++)
    {
      beam_centerline_dofvec.push_back(q_beam(i));
    }
    for (unsigned int i = 0; i < beam_type::n_dof_; i++)
    {
      beam_centerline_dofvec.push_back(v_beam(i));
    }

    // Positional DOFs of the solid.
    q_fluid(2) = -1.0;
    q_fluid(0) = -1.0;
    q_fluid(1) = -1.0;

    q_fluid(5) = 1.0;
    q_fluid(3) = -1.0;
    q_fluid(4) = -1.0;

    q_fluid(8) = 1.0;
    q_fluid(6) = -1.0;
    q_fluid(7) = 1.0;

    q_fluid(11) = -1.0;
    q_fluid(9) = -1.0;
    q_fluid(10) = 1.0;

    q_fluid(14) = -1.0;
    q_fluid(12) = 1.0;
    q_fluid(13) = -1.0;

    q_fluid(17) = 1.0;
    q_fluid(15) = 1.0;
    q_fluid(16) = -1.0;

    q_fluid(20) = 1.0;
    q_fluid(18) = 1.0;
    q_fluid(19) = 1.0;

    q_fluid(23) = -1.0;
    q_fluid(21) = 1.0;
    q_fluid(22) = 1.0;

    v_fluid(0) = 1.0;
    v_fluid(1) = 0.0;
    v_fluid(2) = 0.0;
    v_fluid(3) = 1.0;
    v_fluid(4) = 0.0;
    v_fluid(5) = 0.0;
    v_fluid(6) = 1.0;
    v_fluid(7) = 0.0;
    v_fluid(8) = 0.0;
    v_fluid(9) = 1.0;
    v_fluid(10) = 0.0;
    v_fluid(11) = 0.0;
    v_fluid(12) = 1.0;
    v_fluid(13) = 0.0;
    v_fluid(14) = 0.0;
    v_fluid(15) = 1.0;
    v_fluid(16) = 0.0;
    v_fluid(17) = 0.0;
    v_fluid(18) = 1.0;
    v_fluid(19) = 0.0;
    v_fluid(20) = 0.0;
    v_fluid(21) = 1.0;
    v_fluid(22) = 0.0;
    v_fluid(23) = 0.0;

    for (unsigned int i = 0; i < fluid_type::n_dof_; i++)
    {
      fluid_dofvec.push_back(q_fluid(i));
    }
    for (unsigned int i = 0; i < fluid_type::n_dof_; i++)
    {
      fluid_dofvec.push_back(v_fluid(i));
    }

    results_kff(0, 0) = 0.0026041666666667;
    results_kff(1, 1) = 0.0026041666666667;
    results_kff(2, 2) = 0.0026041666666667;
    results_kff(3, 3) = 0.0234375000000000;
    results_kff(4, 4) = 0.0234375000000000;
    results_kff(5, 5) = 0.0234375000000000;
    results_kff(6, 6) = 0.0234375000000000;
    results_kff(7, 7) = 0.0234375000000000;
    results_kff(8, 8) = 0.0234375000000000;
    results_kff(9, 9) = 0.0026041666666667;
    results_kff(10, 10) = 0.0026041666666667;
    results_kff(11, 11) = 0.0026041666666667;
    results_kff(12, 12) = 0.0234375000000000;
    results_kff(13, 13) = 0.0234375000000000;
    results_kff(14, 14) = 0.0234375000000000;
    results_kff(15, 15) = 0.2109375000000001;
    results_kff(16, 16) = 0.2109375000000001;
    results_kff(17, 17) = 0.2109375000000001;
    results_kff(18, 18) = 0.2109375000000001;
    results_kff(19, 19) = 0.2109375000000001;
    results_kff(20, 20) = 0.2109375000000001;
    results_kff(21, 21) = 0.0234375000000000;
    results_kff(22, 22) = 0.0234375000000000;
    results_kff(23, 23) = 0.0234375000000000;

    results_kfs(0, 0) = 0.0039062500000000;
    results_kfs(0, 1) = 0.0039062500000000;
    results_kfs(0, 2) = 0.0039062500000000;
    results_kfs(0, 3) = 0.0039062500000000;
    results_kfs(0, 4) = 0.0039062500000000;
    results_kfs(0, 5) = 0.0039062500000000;
    results_kfs(0, 6) = 0.0039062500000000;
    results_kfs(0, 7) = 0.0039062500000000;
    results_kfs(0, 8) = 0.0039062500000000;
    results_kfs(0, 9) = 0.0039062500000000;
    results_kfs(0, 10) = 0.0039062500000000;
    results_kfs(0, 11) = 0.0039062500000000;

    // Perform the unit tests.
    PerformGPTSPairUnitTest<beam_type, fluid_type>(pair, q_beam, beam_centerline_dofvec, q_beam_rot,
        q_fluid, fluid_dofvec, results_fs, results_ff, results_ksf, results_kfs, results_kff);
  }

 private:
  /**
   * \brief Set up the contact pair so it can be evaluated and compare the results.
   */
  template <typename beam_type, typename fluid_type>
  void PerformGPTSPairUnitTest(
      BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam_type, fluid_type>& pair,
      const LINALG::Matrix<beam_type::n_dof_, 1, double>& q_beam,
      const std::vector<double>& beam_dofvec, const LINALG::Matrix<9, 1, double>& q_beam_rot,
      const LINALG::Matrix<fluid_type::n_dof_, 1, double>& q_fluid,
      const std::vector<double>& fluid_dofvec, LINALG::SerialDenseVector results_fs,
      LINALG::SerialDenseVector results_ff,
      const LINALG::Matrix<beam_type::n_dof_, fluid_type::n_dof_, double> results_ksf,
      const LINALG::Matrix<fluid_type::n_dof_, beam_type::n_dof_, double> results_kfs,
      const LINALG::Matrix<fluid_type::n_dof_, fluid_type::n_dof_, double> results_kff)
  {
    // Create the elements.
    const int dummy_node_ids[2] = {0, 1};
    Teuchos::RCP<DRT::Element> beam_element = Teuchos::rcp(new DRT::ELEMENTS::Beam3eb(0, 0));
    beam_element->SetNodeIds(2, dummy_node_ids);
    Teuchos::RCP<DRT::ELEMENTS::Fluid> fluid_element = Teuchos::rcp(new DRT::ELEMENTS::Fluid(1, 0));
    fluid_element->SetDisType(DRT::Element::hex8);

    // Set up the beam element.
    std::vector<double> xrefe(6);
    for (unsigned int n = 0; n < 2; n++)
    {
      for (unsigned int i = 0; i < 3; i++)
      {
        xrefe[3 * n + i] = q_beam(6 * n + i);
      }
    }
    // Cast beam element and set the geometry.
    Teuchos::RCP<DRT::ELEMENTS::Beam3eb> beam_element_cast =
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Beam3eb>(beam_element, true);
    beam_element_cast->SetUpReferenceGeometry(xrefe);

    Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> intersection_params =
        Teuchos::rcp(new FBI::BeamToFluidMeshtyingParams());

    // Call Init on the beam contact pair.
    std::vector<const DRT::Element*> pair_elements;
    pair_elements.push_back(&(*beam_element));
    pair_elements.push_back(&(*fluid_element));
    pair.CreateGeometryPair(evaluation_data_);
    pair.Init(intersection_params, pair_elements);
    pair.CastGeometryPair()->Setup();
    pair.ele1posref_ = q_beam;
    pair.ele2posref_ = q_fluid;

    pair.ResetState(beam_dofvec, fluid_dofvec);

    // Evaluate the local matrices.
    LINALG::SerialDenseMatrix local_kff;
    LINALG::SerialDenseMatrix local_kfs;
    LINALG::SerialDenseMatrix local_ksf;
    LINALG::SerialDenseMatrix local_kss;
    LINALG::SerialDenseVector local_fs;
    LINALG::SerialDenseVector local_ff;
    pair.PreEvaluate();
    bool projects =
        pair.Evaluate(&local_fs, &local_ff, &local_kss, &local_ksf, &local_kfs, &local_kff);

    TS_ASSERT(projects);
    TS_ASSERT_EQUALS(local_kff.M(), fluid_type::n_dof_);
    TS_ASSERT_EQUALS(local_kff.N(), fluid_type::n_dof_);
    TS_ASSERT_EQUALS(local_kfs.M(), fluid_type::n_dof_);
    TS_ASSERT_EQUALS(local_kfs.N(), beam_type::n_dof_);
    TS_ASSERT_EQUALS(local_fs.Length(), beam_type::n_dof_);
    for (int i_row = 0; i_row < local_kff.M(); i_row++)
    {
      TS_ASSERT_DELTA((local_kff)(i_row, i_row), results_kff(i_row, i_row), 1e-11);
      TS_ASSERT_DELTA(local_ff(i_row), results_ff(i_row), 1e-11);
      for (int i_col = 0; i_col < local_kfs.N(); i_col++)
        TS_ASSERT_DELTA((local_kfs)(i_row, i_col), local_ksf(i_col, i_row), 1e-11);
    }
    for (unsigned int i_col = 0; i_col < beam_type::n_dof_; i_col++)
      TS_ASSERT_DELTA(local_fs(i_col), results_fs(i_col), 1e-11);
  }


 private:
  //! Evaluation data container for geometry pairs.
  Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> evaluation_data_;
};

#endif
