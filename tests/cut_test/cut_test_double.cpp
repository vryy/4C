/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "cut_test_utils.H"

void test_quad4_surface_mesh_cut()
{
#if 0
  GEO::CUT::MeshIntersection intersection;

  create_quad4_cylinder_mesh( intersection, 0., 0., 60, 10 );
  create_quad4_cylinder_mesh( intersection, 0.5, 0., 60, 10 );

  intersection.SelfCutTest_Cut();
  intersection.Status();
#endif
}

void test_hex8_quad4_double_cut()
{
  GEO::CUT::MeshIntersection intersection(2);
  intersection.GetOptions().Init_for_Cuttests();  // use full cln

  GEO::CUT::Mesh& mesh = intersection.NormalMesh();

  // GEO::CUT::Element * hex8 =
  create_hex8(mesh);

  Epetra_SerialDenseMatrix xyze(3, 4);

  GEO::CUT::Mesh& cut_mesh1 = intersection.CutMesh(0);

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

  GEO::CUT::Mesh& cut_mesh2 = intersection.CutMesh(1);

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
  intersection.Status();

  // OutputGenerator generator;
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
