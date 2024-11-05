// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_intersection.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"

#include "cut_test_utils.hpp"

void test_unit_intersection_touch()
{
  double scale = 1.1e-6;
  for (int i = 0; i < 7; ++i)
  {
    double x = std::pow(0.1, i);
    Cut::Options options;
    options.init_for_cuttests();  // use cln
    Cut::Mesh mesh(options, x);

    Core::LinAlg::SerialDenseMatrix xyze(3, 4);

    xyze(0, 0) = 0;
    xyze(1, 0) = 0;
    xyze(2, 0) = 0;

    xyze(0, 1) = x;
    xyze(1, 1) = 0;
    xyze(2, 1) = 0;

    xyze(0, 2) = x;
    xyze(1, 2) = 0;
    xyze(2, 2) = x;

    xyze(0, 3) = 0;
    xyze(1, 3) = 0;
    xyze(2, 3) = x;

    Cut::Side* s1 = create_quad4(mesh, xyze);

    xyze(0, 0) = 0;
    xyze(1, 0) = -scale * x;
    xyze(2, 0) = 0;

    xyze(0, 1) = 0;
    xyze(1, 1) = x;
    xyze(2, 1) = 0;

    xyze(0, 2) = 0;
    xyze(1, 2) = x;
    xyze(2, 2) = x;

    xyze(0, 3) = 0;
    xyze(1, 3) = scale * x;
    xyze(2, 3) = x;

    Cut::Side* s2 = create_quad4(mesh, xyze);

    const std::vector<Cut::Edge*>& edges = s2->edges();

    Cut::Edge* e = edges[3];

    if (e->nodes()[0]->point()->id() != 7 or e->nodes()[1]->point()->id() != 4)
    {
      FOUR_C_THROW("unexpected nodal id");
    }


    std::shared_ptr<Cut::IntersectionBase> intersection =
        Cut::IntersectionBase::create(Core::FE::CellType::line2, Core::FE::CellType::quad4);
    intersection->init(&mesh, e, s1, false, false, false);

    Cut::PointSet cuts;
    intersection->intersect(cuts);

    for (Cut::PointSet::iterator i = cuts.begin(); i != cuts.end(); ++i)
    {
      Cut::Point* p = *i;
      if (p->id() != 8)
      {
        FOUR_C_THROW("unexpected nodal id");
      }
    }
  }
}
