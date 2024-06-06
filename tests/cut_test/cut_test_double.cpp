/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#include "4C_cut_element.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_meshintersection.hpp"

#include "cut_test_utils.hpp"

void test_hex8_quad4_double_cut()
{
  Core::Geo::Cut::MeshIntersection intersection(2);
  intersection.GetOptions().Init_for_Cuttests();  // use full cln

  Core::Geo::Cut::Mesh& mesh = intersection.NormalMesh();

  // Geo::Cut::Element * hex8 =
  create_hex8(mesh);

  Core::LinAlg::SerialDenseMatrix xyze(3, 4);

  Core::Geo::Cut::Mesh& cut_mesh1 = intersection.CutMesh(0);

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

  // Geo::Cut::Side* quad4_1 =
  create_quad4(cut_mesh1, xyze);

  Core::Geo::Cut::Mesh& cut_mesh2 = intersection.CutMesh(1);

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

  // Geo::Cut::Side* quad4_2 =
  create_quad4(cut_mesh2, xyze);

  // intersection.SelfCutTest_Cut();

  // OutputGenerator generator;
  intersection.CutTest_Cut(true, Inpar::Cut::VCellGaussPts_DirectDivergence);
}
