// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_solid_scatra_ele_calc.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_solid_ele_calc_displacement_based.hpp"
#include "4C_solid_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_ele_calc_fbar.hpp"
#include "4C_solid_scatra_ele.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <array>
#include <cstdlib>
#include <memory>

namespace
{
  using namespace FourC;

  constexpr Core::FE::CellType celltype = Core::FE::CellType::hex8;
  constexpr int num_nodes = Core::FE::num_nodes(celltype);
  constexpr int num_dim = Core::FE::dim<celltype>;
  constexpr int num_dof = num_nodes * num_dim;

  //! thermo-elastic material with the temperature-dependent Young's modulus E(T) = 2.1e5 - 50 T
  std::shared_ptr<Mat::ThermoStVenantKirchhoff> create_material()
  {
    Core::IO::InputParameterContainer material_data;
    material_data.add("YOUNG", std::vector<double>{2.1e5, -50.0});
    material_data.add("NUE", 0.3);
    material_data.add("DENS", 1.0);
    material_data.add("THEXPANS", 1.0e-3);
    material_data.add("INITTEMP", 290.0);
    auto parameters = std::dynamic_pointer_cast<Mat::PAR::ThermoStVenantKirchhoff>(std::shared_ptr(
        Mat::make_parameter(1, Core::Materials::MaterialType::m_thermostvenant, material_data)));

    return {new Mat::ThermoStVenantKirchhoff(parameters.get()),
        [parameters](Mat::ThermoStVenantKirchhoff* material) { delete material; }};
  }

  /*!
   * @brief A single, slightly distorted solid-scatra element with a thermo-elastic material
   *
   * The element is part of a discretization with the displacements in the first and the
   * temperatures in the second dofset, i.e., the setup of the structural discretization in TSI.
   */
  class SolidScatraElementTest : public ::testing::Test
  {
   protected:
    // finite deformation, such that the volume change varies within the element
    const std::array<std::array<double, num_dim>, num_nodes> displacements_ = {{
        {0.03, -0.02, 0.04},
        {-0.01, 0.05, -0.03},
        {0.04, 0.01, 0.02},
        {-0.05, -0.04, 0.01},
        {0.02, 0.03, -0.05},
        {-0.03, 0.02, 0.04},
        {0.01, -0.05, -0.02},
        {0.05, 0.04, 0.03},
    }};

    // nodal temperatures
    const std::array<double, num_nodes> temperatures_ = {
        305.0, 296.0, 309.0, 292.0, 301.0, 307.0, 294.0, 299.0};

    std::shared_ptr<Mat::So3Material> material_ = create_material();

    // empty solid discretization
    Core::FE::Discretization discretization_{"structure", MPI_COMM_WORLD, num_dim};

    // the solid-scatra element
    std::shared_ptr<Discret::Elements::SolidScatra<num_dim>> element_ =
        std::make_shared<Discret::Elements::SolidScatra<num_dim>>(0, 0);
    Core::Elements::LocationArray location_array_{2};

    // set the nodal displacements and temperatures in the discretization
    void set_states(const std::array<std::array<double, num_dim>, num_nodes>& displacements,
        const std::array<double, num_nodes>& temperatures)
    {
      Core::LinAlg::Vector<double> displacement_state(*discretization_.dof_row_map(0), true);
      for (int dof = 0; dof < num_dof; ++dof)
      {
        // sort the displacements into the state vector
        const auto result = std::div(dof, num_dim);
        displacement_state.replace_global_value(
            location_array_[0].lm_[dof], displacements[result.quot][result.rem]);
      }
      discretization_.set_state(0, "displacement", displacement_state);

      Core::LinAlg::Vector<double> temperature_state(*discretization_.dof_row_map(1), true);
      for (int node = 0; node < num_nodes; ++node)
        temperature_state.replace_global_value(location_array_[1].lm_[node], temperatures[node]);
      discretization_.set_state(1, "temperature", temperature_state);
    }

    void SetUp() override
    {
      // distorted reference configuration of the element
      const std::array<std::array<double, num_dim>, num_nodes> reference_coordinates = {{
          {-1.02, -0.97, -1.04},
          {1.03, -1.01, -0.96},
          {0.98, 1.04, -1.02},
          {-0.96, 0.99, -0.98},
          {-1.01, -1.03, 1.02},
          {0.97, -0.98, 1.03},
          {1.04, 1.02, 0.97},
          {-0.99, 0.96, 1.01},
      }};

      // add the element to the discretization
      std::array<int, num_nodes> node_ids{};
      for (int node = 0; node < num_nodes; ++node)
      {
        discretization_.add_node(reference_coordinates[node], node, nullptr);
        node_ids[node] = node;
      }
      element_->set_node_ids(num_nodes, node_ids.data());
      discretization_.add_element(element_);

      // add the dofsets for the displacements and temperatures to the discretization and fill it
      discretization_.add_dof_set(
          std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(1, 0, 0, true));
      discretization_.fill_complete();
      element_->location_vector(discretization_, location_array_);
    }
  };

  //! evaluation of the element with the given solid formulation
  template <typename SolidFormulation>
  class SolidScatraEleCalcTest : public SolidScatraElementTest
  {
   protected:
    // the calculation class for the solid-scatra element with the given solid formulation
    Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation> element_calc_;

    void SetUp() override
    {
      SolidScatraElementTest::SetUp();
      element_calc_.setup(*material_, Core::IO::InputParameterContainer{});
    }

    //! internal force vector of the element for the given nodal displacements and temperatures
    Core::LinAlg::SerialDenseVector internal_force(
        const std::array<std::array<double, num_dim>, num_nodes>& displacements,
        const std::array<double, num_nodes>& temperatures)
    {
      set_states(displacements, temperatures);
      Teuchos::ParameterList params;
      Core::LinAlg::SerialDenseVector force(num_dof);
      element_calc_.evaluate_nonlinear_force_stiffness_mass(*element_, *material_, discretization_,
          location_array_, params, &force, nullptr, nullptr);
      return force;
    }

    //! analytic linearization of the internal force vector wrt the nodal temperatures
    Core::LinAlg::SerialDenseMatrix d_internal_force_d_temperatures(
        const std::array<std::array<double, num_dim>, num_nodes>& displacements,
        const std::array<double, num_nodes>& temperatures)
    {
      set_states(displacements, temperatures);
      Teuchos::ParameterList params;
      Core::LinAlg::SerialDenseMatrix d_force_d_temperatures(num_dof, num_nodes);
      element_calc_.evaluate_d_stress_d_scalar(
          *element_, *material_, discretization_, location_array_, params, d_force_d_temperatures);
      return d_force_d_temperatures;
    }
  };

  using Formulations = ::testing::Types<Discret::Elements::DisplacementBasedFormulation<celltype>,
      Discret::Elements::DisplacementBasedLinearKinematicsFormulation<celltype>,
      Discret::Elements::FBarFormulation<celltype>>;
  TYPED_TEST_SUITE(SolidScatraEleCalcTest, Formulations);

  TYPED_TEST(SolidScatraEleCalcTest, InternalForceWrtTemperatureMatchesFiniteDifferences)
  {
    const auto d_force_d_temperatures =
        this->d_internal_force_d_temperatures(this->displacements_, this->temperatures_);

    // Here, The internal force is quadratic in the temperature since the Young's modulus is linear
    // in it. Thus, the central difference is exact and the tolerance only accounts for round-off.
    const double perturbation_size = 1.0e-4;
    Core::LinAlg::SerialDenseMatrix d_force_d_temperatures_FD(num_dof, num_nodes);
    for (int node = 0; node < num_nodes; ++node)
    {
      auto temperatures_plus = this->temperatures_;
      auto temperatures_minus = this->temperatures_;
      temperatures_plus[node] += perturbation_size;
      temperatures_minus[node] -= perturbation_size;
      const auto force_plus = this->internal_force(this->displacements_, temperatures_plus);
      const auto force_minus = this->internal_force(this->displacements_, temperatures_minus);
      for (int dof = 0; dof < num_dof; ++dof)
        d_force_d_temperatures_FD(dof, node) =
            (force_plus[dof] - force_minus[dof]) / (2.0 * perturbation_size);
    }

    FOUR_C_EXPECT_NEAR(d_force_d_temperatures, d_force_d_temperatures_FD, 1.0e-6);
  }
}  // namespace
