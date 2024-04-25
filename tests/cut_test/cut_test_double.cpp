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
  CORE::GEO::CUT::MeshIntersection intersection(2);
  intersection.GetOptions().Init_for_Cuttests();  // use full cln

  CORE::GEO::CUT::Mesh& mesh = intersection.NormalMesh();

  // GEO::CUT::Element * hex8 =
  create_hex8(mesh);

  CORE::LINALG::SerialDenseMatrix xyze(3, 4);

  CORE::GEO::CUT::Mesh& cut_mesh1 = intersection.CutMesh(0);

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

  // GEO::CUT::Side* quad4_1 =
  create_quad4(cut_mesh1, xyze);

  CORE::GEO::CUT::Mesh& cut_mesh2 = intersection.CutMesh(1);

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

  // GEO::CUT::Side* quad4_2 =
  create_quad4(cut_mesh2, xyze);

  // intersection.SelfCutTest_Cut();

  // OutputGenerator generator;
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
