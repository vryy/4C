// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_element.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_meshintersection.hpp"

#include "cut_test_utils.hpp"

void test_hex8_quad4_double_cut()
{
  Cut::MeshIntersection intersection(2);
  intersection.get_options().init_for_cuttests();  // use full cln

  Cut::Mesh& mesh = intersection.normal_mesh();

  // Cut::Element * hex8 =
  create_hex8(mesh);

  Core::LinAlg::SerialDenseMatrix xyze(3, 4);

  Cut::Mesh& cut_mesh1 = intersection.cut_mesh(0);

  xyze(0, 0) = 0.25;
  xyze(1, 0) = -0.2;
  xyze(2, 0) = -0.2;

  xyze(0, 1) = 0.75;
  xyze(1, 1) = -0.2;
  xyze(2, 1) = 1.2;

  xyze(0, 2) = 0.75;
  xyze(1, 2) = 1.2;
  xyze(2, 2) = 1.2;

  xyze(0, 3) = 0.25;
  xyze(1, 3) = 1.2;
  xyze(2, 3) = -0.2;

  // Cut::Side* quad4_1 =
  create_quad4(cut_mesh1, xyze);

  Cut::Mesh& cut_mesh2 = intersection.cut_mesh(1);

  xyze(0, 0) = 0.75;
  xyze(1, 0) = -0.2;
  xyze(2, 0) = -0.2;

  xyze(0, 1) = 0.25;
  xyze(1, 1) = -0.2;
  xyze(2, 1) = 1.2;

  xyze(0, 2) = 0.25;
  xyze(1, 2) = 1.2;
  xyze(2, 2) = 1.2;

  xyze(0, 3) = 0.75;
  xyze(1, 3) = 1.2;
  xyze(2, 3) = -0.2;

  // Cut::Side* quad4_2 =
  create_quad4(cut_mesh2, xyze);

  // intersection.SelfCutTest_Cut();

  // OutputGenerator generator;
  intersection.cut_test_cut(true, Cut::VCellGaussPts_DirectDivergence);
}
