/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

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
    Core::Geo::Cut::Options options;
    options.Init_for_Cuttests();  // use cln
    Core::Geo::Cut::Mesh mesh(options, x);

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

    Core::Geo::Cut::Side* s1 = create_quad4(mesh, xyze);

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

    Core::Geo::Cut::Side* s2 = create_quad4(mesh, xyze);

    const std::vector<Core::Geo::Cut::Edge*>& edges = s2->Edges();

    Core::Geo::Cut::Edge* e = edges[3];

    if (e->Nodes()[0]->point()->Id() != 7 or e->Nodes()[1]->point()->Id() != 4)
    {
      FOUR_C_THROW("unexpected nodal id");
    }


    Teuchos::RCP<Core::Geo::Cut::IntersectionBase> intersection =
        Core::Geo::Cut::IntersectionBase::Create(
            Core::FE::CellType::line2, Core::FE::CellType::quad4);
    intersection->Init(&mesh, e, s1, false, false, false);

    Core::Geo::Cut::PointSet cuts;
    intersection->Intersect(cuts);

    for (Core::Geo::Cut::PointSet::iterator i = cuts.begin(); i != cuts.end(); ++i)
    {
      Core::Geo::Cut::Point* p = *i;
      if (p->Id() != 8)
      {
        FOUR_C_THROW("unexpected nodal id");
      }
    }
  }
}
