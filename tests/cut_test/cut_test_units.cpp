/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#include "baci_cut_intersection.H"
#include "baci_cut_mesh.H"
#include "baci_cut_options.H"
#include "baci_cut_side.H"

#include "cut_test_utils.H"

void test_unit_intersection_touch()
{
  double scale = 1.1e-6;
  for (int i = 0; i < 7; ++i)
  {
    double x = std::pow(0.1, i);
    CORE::GEO::CUT::Options options;
    options.Init_for_Cuttests();  // use cln
    CORE::GEO::CUT::Mesh mesh(options, x);

    CORE::LINALG::SerialDenseMatrix xyze(3, 4);

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

    CORE::GEO::CUT::Side* s1 = create_quad4(mesh, xyze);

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

    CORE::GEO::CUT::Side* s2 = create_quad4(mesh, xyze);

    const std::vector<CORE::GEO::CUT::Edge*>& edges = s2->Edges();

    CORE::GEO::CUT::Edge* e = edges[3];

    if (e->Nodes()[0]->point()->Id() != 7 or e->Nodes()[1]->point()->Id() != 4)
    {
      dserror("unexpected nodal id");
    }


    Teuchos::RCP<CORE::GEO::CUT::IntersectionBase> intersection =
        CORE::GEO::CUT::IntersectionBase::Create(
            CORE::FE::CellType::line2, CORE::FE::CellType::quad4);
    intersection->Init(&mesh, e, s1, false, false, false);

    CORE::GEO::CUT::PointSet cuts;
    intersection->Intersect(cuts);

    for (CORE::GEO::CUT::PointSet::iterator i = cuts.begin(); i != cuts.end(); ++i)
    {
      CORE::GEO::CUT::Point* p = *i;
      if (p->Id() != 8)
      {
        dserror("unexpected nodal id");
      }
    }
  }
}
