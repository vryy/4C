// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_beaminteraction_contact_beam_to_solid_mortar_shape_functions_dual_hermite.hpp"

namespace
{
  using namespace FourC;

  TEST(BeamInteractionMortarShapeFunctionsTest, TestDualHermiteLine2Length03)
  {
    const double xi = 0.123456789;
    const double ref_length = 0.3;

    Core::LinAlg::Matrix<4, 1> phi(Core::LinAlg::Initialization::zero);

    GeometryPair::ShapeFunctionData<BeamInteraction::t_hermite_dual> shape_function_data;
    shape_function_data.ref_length_ = ref_length;

    GeometryPair::EvaluateShapeFunction<BeamInteraction::t_hermite_dual>::evaluate(
        phi, xi, shape_function_data);

    const std::vector<double> phi_ref = {
        -0.24634578918994065, -0.029314912689680812, -1.1393423701836304, 7.981608777268251};

    for (std::size_t i = 0; i < phi_ref.size(); ++i) EXPECT_NEAR(phi(i), phi_ref[i], 1.0e-10);
  }

}  // namespace
