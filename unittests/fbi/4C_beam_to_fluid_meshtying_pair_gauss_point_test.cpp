/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the beam to fluid meshtying gpts pair.

\level 2
*/


#include <gtest/gtest.h>

#include "4C_beam3_euler_bernoulli.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_pair_gauss_point.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_volume_segmentation.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <vector>

namespace
{
  using namespace FourC;

  /**
   * Class to test the local coupling matrices calculated by the beam to fluid meshtying gpts pair.
   */
  class BeamToFluidMeshtyingPairGPTSTest : public ::testing::Test

  {
   public:
    /**
     * \brief Set up the testing environment.
     */
    BeamToFluidMeshtyingPairGPTSTest()
    {
      // Set up the evaluation data container for the geometry pairs.
      Teuchos::ParameterList line_to_volume_params_list;
      Inpar::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_volume_params_list);
      evaluation_data_ =
          Teuchos::rcp(new GEOMETRYPAIR::LineTo3DEvaluationData(line_to_volume_params_list));
    }

    /**
     * \brief Set up the pair so it can be evaluated and compare the results.
     */
    template <typename beam_type, typename fluid_type>
    void perform_gpts_pair_unit_test(
        const Core::LinAlg::Matrix<beam_type::n_dof_, 1, double>& q_beam,
        const std::vector<double>& beam_dofvec,
        const Core::LinAlg::Matrix<fluid_type::n_dof_, 1, double>& q_fluid,
        const std::vector<double>& fluid_dofvec, Core::LinAlg::SerialDenseVector results_fs,
        Core::LinAlg::SerialDenseVector results_ff,
        const Core::LinAlg::Matrix<beam_type::n_dof_, fluid_type::n_dof_, double> results_ksf,
        const Core::LinAlg::Matrix<fluid_type::n_dof_, beam_type::n_dof_, double> results_kfs,
        const Core::LinAlg::Matrix<fluid_type::n_dof_, fluid_type::n_dof_, double> results_kff)
    {
      using scalar_type = GEOMETRYPAIR::line_to_volume_scalar_type<beam_type, fluid_type>;

      // Create the mesh tying mortar pair.
      BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam_type, fluid_type> pair =
          BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam_type, fluid_type>();

      // Create the elements.
      const int dummy_node_ids[2] = {0, 1};
      Teuchos::RCP<Core::Elements::Element> beam_element =
          Teuchos::rcp(new Discret::ELEMENTS::Beam3eb(0, 0));
      beam_element->SetNodeIds(2, dummy_node_ids);
      Teuchos::RCP<Discret::ELEMENTS::Fluid> fluid_element =
          Teuchos::rcp(new Discret::ELEMENTS::Fluid(1, 0));
      fluid_element->SetDisType(Core::FE::CellType::hex8);

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
      Teuchos::RCP<Discret::ELEMENTS::Beam3eb> beam_element_cast =
          Teuchos::rcp_dynamic_cast<Discret::ELEMENTS::Beam3eb>(beam_element, true);
      beam_element_cast->set_up_reference_geometry(xrefe);

      Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> intersection_params =
          Teuchos::rcp(new FBI::BeamToFluidMeshtyingParams());

      // Call Init on the beam contact pair.
      std::vector<const Core::Elements::Element*> pair_elements;
      pair_elements.push_back(&(*beam_element));
      pair_elements.push_back(&(*fluid_element));
      pair.CreateGeometryPair(pair_elements[0], pair_elements[1], evaluation_data_);
      pair.Init(intersection_params, pair_elements);

      pair.ele1pos_ = GEOMETRYPAIR::InitializeElementData<beam_type, scalar_type>::Initialize(
          beam_element.get());
      pair.ele1posref_ =
          GEOMETRYPAIR::InitializeElementData<beam_type, double>::Initialize(beam_element.get());
      pair.ele1poscur_ =
          GEOMETRYPAIR::InitializeElementData<beam_type, double>::Initialize(beam_element.get());
      pair.ele1vel_ = GEOMETRYPAIR::InitializeElementData<beam_type, scalar_type>::Initialize(
          beam_element.get());
      pair.ele1posref_.element_position_ = q_beam;
      pair.ele2posref_.element_position_ = q_fluid;

      pair.ResetState(beam_dofvec, fluid_dofvec);

      const int fluid_dofs = fluid_type::n_dof_;
      const int beam_dofs = beam_type::n_dof_;

      // Evaluate the local matrices.
      Core::LinAlg::SerialDenseMatrix local_kff;
      Core::LinAlg::SerialDenseMatrix local_kfs;
      Core::LinAlg::SerialDenseMatrix local_ksf;
      Core::LinAlg::SerialDenseMatrix local_kss;
      Core::LinAlg::SerialDenseVector local_fs;
      Core::LinAlg::SerialDenseVector local_ff;
      pair.pre_evaluate();
      bool projects =
          pair.evaluate(&local_fs, &local_ff, &local_kss, &local_ksf, &local_kfs, &local_kff);

      EXPECT_TRUE(projects);
      EXPECT_EQ(local_kff.numRows(), fluid_dofs);
      EXPECT_EQ(local_kff.numCols(), fluid_dofs);
      EXPECT_EQ(local_kfs.numRows(), fluid_dofs);
      EXPECT_EQ(local_kfs.numCols(), beam_dofs);
      EXPECT_EQ(local_fs.length(), beam_dofs);


      for (int i_row = 0; i_row < local_kff.numRows(); i_row++)
      {
        EXPECT_NEAR((local_kff)(i_row, i_row), results_kff(i_row, i_row), 1e-11)
            << " for i_row = " << i_row;
        EXPECT_NEAR(local_ff(i_row), results_ff(i_row), 1e-11) << " for i_row = " << i_row;
        for (int i_col = 0; i_col < local_kfs.numCols(); i_col++)
          EXPECT_NEAR((local_kfs)(i_row, i_col), local_ksf(i_col, i_row), 1e-11)
              << " for i_row = " << i_row << ", i_col = " << i_col;
      }
      for (unsigned int i_col = 0; i_col < beam_type::n_dof_; i_col++)
        EXPECT_NEAR(local_fs(i_col), results_fs(i_col), 1e-11) << " for i_col = " << i_col;
    }


   private:
    //! Evaluation data container for geometry pairs.
    Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> evaluation_data_;
  };

  /**
   * \brief Test a moving straight beam in a hex8 element with hermite line2 shape functions.
   */
  TEST_F(BeamToFluidMeshtyingPairGPTSTest, TestBeamToFluidMeshtyingHex8MovingBeam)
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex8 fluid_type;

    // Definition of variables for this test case.
    Core::LinAlg::Matrix<beam_type::n_dof_, 1, double> q_beam;
    Core::LinAlg::Matrix<beam_type::n_dof_, 1, double> v_beam;
    Core::LinAlg::Matrix<9, 1, double> q_beam_rot;
    Core::LinAlg::Matrix<fluid_type::n_dof_, 1, double> q_fluid;
    Core::LinAlg::Matrix<fluid_type::n_dof_, 1, double> v_fluid;
    std::vector<double> beam_centerline_dofvec;
    std::vector<double> fluid_dofvec;

    // Matrices for the results.
    Core::LinAlg::Matrix<fluid_type::n_dof_, fluid_type::n_dof_, double> results_kff(true);
    Core::LinAlg::Matrix<fluid_type::n_dof_, beam_type::n_dof_, double> results_kfs(true);
    Core::LinAlg::Matrix<beam_type::n_dof_, fluid_type::n_dof_, double> results_ksf(true);
    Core::LinAlg::SerialDenseVector results_fs(beam_type::n_dof_, true);
    Core::LinAlg::SerialDenseVector results_ff(fluid_type::n_dof_, true);
    results_fs.putScalar(0.0);
    results_ff.putScalar(0.0);


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

    // Coordinates and velocity DOFs of the fluid element
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
    perform_gpts_pair_unit_test<beam_type, fluid_type>(q_beam, beam_centerline_dofvec, q_fluid,
        fluid_dofvec, results_fs, results_ff, results_ksf, results_kfs, results_kff);
  }
}  // namespace
