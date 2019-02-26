/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_units.cpp

\level 1

\maintainer Christoph Ager
*----------------------------------------------------------------------*/

#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_intersection.H"
#include "../../src/drt_cut/cut_side.H"
#include "cut_test_utils.H"

void test_unit_intersection_touch()
{
  double scale = 1.1e-6;
  for (int i = 0; i < 7; ++i)
  {
    double x = std::pow(0.1, i);
    GEO::CUT::Options options;
    options.Init_for_Cuttests();  // use cln
    GEO::CUT::Mesh mesh(options, x);

    Epetra_SerialDenseMatrix xyze(3, 4);

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

    GEO::CUT::Side* s1 = create_quad4(mesh, xyze);

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

    GEO::CUT::Side* s2 = create_quad4(mesh, xyze);

    const std::vector<GEO::CUT::Edge*>& edges = s2->Edges();

    GEO::CUT::Edge* e = edges[3];

    if (e->Nodes()[0]->point()->Id() != 7 or e->Nodes()[1]->point()->Id() != 4)
    {
      throw std::runtime_error("unexpected nodal id");
    }


    Teuchos::RCP<GEO::CUT::IntersectionBase> intersection =
        GEO::CUT::IntersectionBase::Create(DRT::Element::line2, DRT::Element::quad4);
    intersection->Init(&mesh, e, s1, false, false, false);

    GEO::CUT::PointSet cuts;
    intersection->Intersect(cuts);

    for (GEO::CUT::PointSet::iterator i = cuts.begin(); i != cuts.end(); ++i)
    {
      GEO::CUT::Point* p = *i;
      if (p->Id() != 8)
      {
        run_time_error("unexpected nodal id");
      }
    }
  }
}
