/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_simple.cpp

\level 1

\maintainer Ager Christoph
*----------------------------------------------------------------------*/

#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_position.H"
#include "cut_test_utils.H"
#include "../../src/drt_cut/cut_triangulateFacet.H"

GEO::CUT::Element* create_tet4(GEO::CUT::Mesh& mesh)
{
  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 2;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 2;
  xyze(1, 1) = 0;
  xyze(2, 1) = 1;

  xyze(0, 2) = 2;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 0.5;

  return create_tet4(mesh, xyze);
}

GEO::CUT::Side* create_quad4(
    GEO::CUT::Mesh& mesh, double x, double dx, double dz, bool reverse = false)
{
  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = x - dx;
  xyze(1, 0) = -0.5;
  xyze(2, 0) = -0.5 - dz;

  xyze(0, 1) = x + dx;
  xyze(1, 1) = -0.5;
  xyze(2, 1) = 1.5 + dz;

  xyze(0, 2) = x + dx;
  xyze(1, 2) = 1.5;
  xyze(2, 2) = 1.5 + dz;

  xyze(0, 3) = x - dx;
  xyze(1, 3) = 1.5;
  xyze(2, 3) = -0.5 - dz;

  if (reverse)
  {
    std::swap(xyze(0, 1), xyze(0, 3));
    std::swap(xyze(1, 1), xyze(1, 3));
    std::swap(xyze(2, 1), xyze(2, 3));
  }

  return create_quad4(mesh, xyze);
}

void test_hex8_simple()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);
  GEO::CUT::Element* e = create_hex8(mesh);
  GEO::CUT::Side* s = create_quad4(mesh, 0.5, 0.1, 0);

  mesh.Status();

  e->Cut(mesh, *(s));

  cutmesh(mesh);
}

void test_tet4_simple()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 0;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 1;

  GEO::CUT::Element* e = create_tet4(mesh, xyze);
  GEO::CUT::Side* s = create_quad4(mesh, 0.5, 0.1, 0);

  mesh.Status();

  e->Cut(mesh, *(s));

  cutmesh(mesh);
}

void test_pyramid5_simple()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  Epetra_SerialDenseMatrix xyze(3, 5);

  xyze(0, 0) = 0;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0;
  xyze(1, 3) = 1;
  xyze(2, 3) = 0;

  xyze(0, 4) = 0.5;
  xyze(1, 4) = 1;
  xyze(2, 4) = 1;

  GEO::CUT::Element* e = create_pyramid5(mesh, xyze);
  GEO::CUT::Side* s = create_quad4(mesh, 0.5, 0.1, 0);

  mesh.Status();

  e->Cut(mesh, *(s));

  cutmesh(mesh);
}

void test_wedge6_simple()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  Epetra_SerialDenseMatrix xyze(3, 6);

  xyze(0, 0) = 0;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0;
  xyze(1, 3) = 0;
  xyze(2, 3) = 1;

  xyze(0, 4) = 1;
  xyze(1, 4) = 0;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = 1;
  xyze(2, 5) = 1;

  GEO::CUT::Element* e = create_wedge6(mesh, xyze);
  GEO::CUT::Side* s = create_quad4(mesh, 0.5, 0.1, 0);

  mesh.Status();

  e->Cut(mesh, *(s));

  cutmesh(mesh);
}


void test_hex8_fullside()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);
  GEO::CUT::Element* e = create_hex8(mesh);
  GEO::CUT::Side* s = create_quad4(mesh, 1, 0, 0);

  mesh.Status();

  e->Cut(mesh, *(s));

  cutmesh(mesh);
}

void test_hex8_diagonal()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);
  GEO::CUT::Element* e = create_hex8(mesh);
  GEO::CUT::Side* s = create_quad4(mesh, 0.5, 1, 0);

  mesh.Status();

  e->Cut(mesh, *(s));

  cutmesh(mesh);
}


void test_hex8_tet4()
{
  SimpleWrapper w;

  w.CreateHex8();
  w.CreateTet4Sides();

  Epetra_SerialDenseMatrix xyze(3, 4);

  // add second cut to be able to find nodal positions

  xyze(0, 0) = -0.1;
  xyze(1, 0) = 0.1;
  xyze(2, 0) = -0.1;

  xyze(0, 1) = 1.1;
  xyze(1, 1) = 0.1;
  xyze(2, 1) = -0.1;

  xyze(0, 2) = 1.1;
  xyze(1, 2) = -0.1;
  xyze(2, 2) = 0.1;

  xyze(0, 3) = -0.1;
  xyze(1, 3) = -0.1;
  xyze(2, 3) = 0.1;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_hex8()
{
  SimpleWrapper w;

  w.CreateHex8();
  w.CreateHex8Sides(0.5, 0.5, 0.5);
  w.Status();
  w.CutTest_Cut();
}

void test_hex8_touch()
{
  SimpleWrapper w;

  w.CreateHex8();
  w.CreateHex8Sides(1, 0, 0);

  Epetra_SerialDenseMatrix xyze(3, 4);

  // add second cut to be able to find nodal positions

  xyze(0, 0) = 0.1;
  xyze(1, 0) = -0.1;
  xyze(2, 0) = -0.1;

  xyze(0, 1) = 0.1;
  xyze(1, 1) = 1.1;
  xyze(2, 1) = -0.1;

  xyze(0, 2) = -0.1;
  xyze(1, 2) = 1.1;
  xyze(2, 2) = 0.1;

  xyze(0, 3) = -0.1;
  xyze(1, 3) = -0.1;
  xyze(2, 3) = 0.1;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_touch2()
{
  SimpleWrapper w;

  w.CreateHex8();
  w.CreateHex8Sides(1, 0.5, 0.5);
  w.Status();
  w.CutTest_Cut();
}

void test_hex8_schraeg()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 0.5;
  xyze(1, 1) = 1;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = 1;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0;
  xyze(2, 3) = 1;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_woelbung()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = -0.5;
  xyze(1, 0) = -0.5;
  xyze(2, 0) = -1.5;

  xyze(0, 1) = 2.5;
  xyze(1, 1) = -0.5;
  xyze(2, 1) = 1.5;

  xyze(0, 2) = 2.5;
  xyze(1, 2) = 1.5;
  xyze(2, 2) = -1.5;

  xyze(0, 3) = -0.5;
  xyze(1, 3) = 1.5;
  xyze(2, 3) = 1.5;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
  w.AssumeVolumeCells(2);
}

void test_hex8_tet4_touch()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 2;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 2;
  xyze(1, 1) = 0;
  xyze(2, 1) = 1;

  xyze(0, 2) = 2;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 0.5;

  w.CreateTet4Sides(xyze);

  // add second cut to be able to find nodal positions

  xyze(0, 0) = -0.1;
  xyze(1, 0) = 0.1;
  xyze(2, 0) = -0.1;

  xyze(0, 1) = 1.1;
  xyze(1, 1) = 0.1;
  xyze(2, 1) = -0.1;

  xyze(0, 2) = 1.1;
  xyze(1, 2) = -0.1;
  xyze(2, 2) = 0.1;

  xyze(0, 3) = -0.1;
  xyze(1, 3) = -0.1;
  xyze(2, 3) = 0.1;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut(true, true);  // as cut_sides are just touching!!
}

void test_hex8_tet4_touch2()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = -0.5;

  xyze(0, 1) = 1.5;
  xyze(1, 1) = 0;
  xyze(2, 1) = -0.5;

  xyze(0, 2) = 1.5;
  xyze(1, 2) = 1;
  xyze(2, 2) = -0.5;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 1.5;

  w.CreateTet4Sides(xyze);

  w.Status();
  w.CutTest_Cut(true, true);
}

void test_hex8_mesh()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  create_hex8_mesh(mesh, 10, 10, 10);

  GEO::CUT::Side* s = create_quad4(mesh, 0.5, 0.5, 0);

  mesh.Status();

  GEO::CUT::plain_element_set done;
  GEO::CUT::plain_element_set elements_done;
  mesh.Cut(*(s), done, elements_done);

  cutmesh(mesh);
}

void test_hex8_double()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);
  GEO::CUT::Element* e = create_hex8(mesh);
  GEO::CUT::Side* s1 = create_quad4(mesh, 0.4, 0.1, 0);
  GEO::CUT::Side* s2 = create_quad4(mesh, 0.6, 0.1, 0);

  mesh.Status();

  e->Cut(mesh, *(s1));
  e->Cut(mesh, *(s2));

  cutmesh(mesh);
}

void test_hex8_bad1()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = 0.7291666666666666;
  xyze(1, 0) = 0.5208333333332368;
  xyze(2, 0) = 0.02500000000896939;

  xyze(0, 1) = 0.7291666666666667;
  xyze(1, 1) = 0.5208333333333334;
  xyze(2, 1) = 0;

  xyze(0, 2) = 0.75;
  xyze(1, 2) = 0.5208333333333334;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0.7499999999999999;
  xyze(1, 3) = 0.5208333333332476;
  xyze(2, 3) = 0.02500000000797485;

  xyze(0, 4) = 0.7291666666666667;
  xyze(1, 4) = 0.5;
  xyze(2, 4) = 0.025;

  xyze(0, 5) = 0.7291666666666667;
  xyze(1, 5) = 0.5;
  xyze(2, 5) = 0;

  xyze(0, 6) = 0.75;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = 0;

  xyze(0, 7) = 0.75;
  xyze(1, 7) = 0.5;
  xyze(2, 7) = 0.025;

  GEO::CUT::Element* e = create_hex8(mesh, xyze);

  xyze(0, 0) = 0.75;
  xyze(1, 0) = 0.5010108360343256;
  xyze(2, 0) = 0;

  xyze(0, 1) = 0.7435592801990288;
  xyze(1, 1) = 0.5208333333333334;
  xyze(2, 1) = 0;

  xyze(0, 2) = 0.7435592801990578;
  xyze(1, 2) = 0.5208333333332442;
  xyze(2, 2) = 0.02500000000828232;

  xyze(0, 3) = 0.75;
  xyze(1, 3) = 0.5010108360343257;
  xyze(2, 3) = 0.02500000000038694;

  GEO::CUT::Side* quad4 = create_quad4(mesh, xyze);

  mesh.Status();

  e->Cut(mesh, *(quad4));

  cutmesh(mesh);
}

void test_hex8_bad2()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = 1.05556;
  xyze(1, 0) = 0.444444;
  xyze(2, 0) = -4.82103e-20;

  xyze(0, 1) = 1.05556;
  xyze(1, 1) = 0.444444;
  xyze(2, 1) = -0.05;

  xyze(0, 2) = 1.05556;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = -0.05;

  xyze(0, 3) = 1.05556;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 0;

  xyze(0, 4) = 1.11111;
  xyze(1, 4) = 0.444444;
  xyze(2, 4) = 1.41172e-22;

  xyze(0, 5) = 1.11111;
  xyze(1, 5) = 0.444444;
  xyze(2, 5) = -0.05;

  xyze(0, 6) = 1.11111;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = -0.05;

  xyze(0, 7) = 1.11111;
  xyze(1, 7) = 0.5;
  xyze(2, 7) = 0;

  GEO::CUT::Element* e = create_hex8(mesh, xyze);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = -0.0505;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0.5;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1.05714;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = -9.3343e-19;

  xyze(0, 3) = 1.05714;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = -0.0505;

  GEO::CUT::Side* quad4 = create_quad4(mesh, xyze);

  mesh.Status();

  e->Cut(mesh, *(quad4));

  cutmesh(mesh);
}

void test_hex8_bad3()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = 1.05556;
  xyze(1, 0) = 0.444444;
  xyze(2, 0) = 0.05;

  xyze(0, 1) = 1.05556;
  xyze(1, 1) = 0.444444;
  xyze(2, 1) = -4.82103e-20;

  xyze(0, 2) = 1.05556;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = 0;

  xyze(0, 3) = 1.05556;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 0.05;

  xyze(0, 4) = 1.11111;
  xyze(1, 4) = 0.444444;
  xyze(2, 4) = 0.05;

  xyze(0, 5) = 1.11111;
  xyze(1, 5) = 0.444444;
  xyze(2, 5) = 1.41172e-22;

  xyze(0, 6) = 1.11111;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = 0;

  xyze(0, 7) = 1.11111;
  xyze(1, 7) = 0.5;
  xyze(2, 7) = 0.05;

  GEO::CUT::Element* e = create_hex8(mesh, xyze);

  xyze(0, 0) = 1.05714;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = -9.3343e-19;

  xyze(0, 1) = 1.05714;
  xyze(1, 1) = 0.5;
  xyze(2, 1) = 0.0505;

  xyze(0, 2) = 1.11429;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = 0.0505;

  xyze(0, 3) = 1.11429;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 1.60089e-18;

  GEO::CUT::Side* quad4 = create_quad4(mesh, xyze);

  mesh.Status();

  e->Cut(mesh, *(quad4));

  cutmesh(mesh);
}

/*
 * (0.944444,0,0.05), (0.944444,0,6.55604e-19),
 * (0.944444,0.0555556,1.29045e-18), (0.944444,0.0555556,0.05), (1,0,0.05),
 * (1,0,0), (1,0.0555556,3.40507e-19), (1,0.0555556,0.05)
 *
 * {(1,0,0.0505), (1,0.0555556,0.0505), (1,0.0555556,3.40507e-19), (1,0,0)}
 */

void test_hex8_bad4()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  double hex8_xyz[24] = {0.944444, 0, 0.05, 0.944444, 0, 6.55604e-19, 0.944444, 0.0555556,
      1.29045e-18, 0.944444, 0.0555556, 0.05, 1, 0, 0.05, 1, 0, 0, 1, 0.0555556, 3.40507e-19, 1,
      0.0555556, 0.05};

  double quad4_xyz[12] = {1, 0, 0.0505, 1, 0.0555556, 0.0505, 1, 0.0555556, 3.40507e-19, 1, 0, 0};

  Epetra_SerialDenseMatrix xyze(3, 8);

  std::copy(hex8_xyz, hex8_xyz + 24, &xyze(0, 0));

  GEO::CUT::Element* e = create_hex8(mesh, xyze);

  std::copy(quad4_xyz, quad4_xyz + 12, &xyze(0, 0));

  GEO::CUT::Side* quad4 = create_quad4(mesh, xyze);

  mesh.Status();

  e->Cut(mesh, *(quad4));

  cutmesh(mesh);
}

void test_hex8_wedge6()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 6);

  xyze(0, 0) = 0.5;
  xyze(1, 0) = 2;
  xyze(2, 0) = -0.5;

  xyze(0, 1) = 0.5;
  xyze(1, 1) = 0.5;
  xyze(2, 1) = -0.5;

  xyze(0, 2) = 3;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = -0.5;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 2;
  xyze(2, 3) = 1.5;

  xyze(0, 4) = 0.5;
  xyze(1, 4) = 0.5;
  xyze(2, 4) = 1.5;

  xyze(0, 5) = 3;
  xyze(1, 5) = 0.5;
  xyze(2, 5) = 1.5;

  w.CreateWedge6Sides(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_touch()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 1.5;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1.5;
  xyze(2, 2) = 1.5;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0;
  xyze(2, 3) = 1.5;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut(true, true);
}

void test_hex8_quad4_touch2()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 1.5;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1.5;
  xyze(2, 2) = 1.5;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 1.5;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut(true, true);
}

void test_hex8_quad4_touch3()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = -0.5;

  xyze(0, 1) = 1;
  xyze(1, 1) = 1.5;
  xyze(2, 1) = -0.5;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1.5;
  xyze(2, 2) = 1.5;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 1.5;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut(true, true);
}

void test_hex8_quad4_cut()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 0.5;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 0.5;
  xyze(1, 1) = 1;
  xyze(2, 1) = 0;

  xyze(0, 2) = 0.5;
  xyze(1, 2) = 1;
  xyze(2, 2) = 1;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 0;
  xyze(2, 3) = 1;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_gedreht()
{
  SimpleWrapper w;

  w.CreateHex8();

#if 0
  Epetra_SerialDenseMatrix xyze( 3, 4 );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = -0.2;

  xyze( 0, 1 ) = 0.5;
  xyze( 1, 1 ) = 1.2;
  xyze( 2, 1 ) = 0.5;

  xyze( 0, 2 ) = 0.5;
  xyze( 1, 2 ) = 0.5;
  xyze( 2, 2 ) = 1.2;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = -0.2;
  xyze( 2, 3 ) = 0.5;

  GEO::CUT::Side* quad4 = create_quad4( mesh, xyze );
#endif

  w.CreateQuad4Mesh(2, 2);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_hex8_durchstoss()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = -0.5;
  xyze(1, 0) = 0.2;
  xyze(2, 0) = 0.2;

  xyze(0, 1) = -0.5;
  xyze(1, 1) = 0.8;
  xyze(2, 1) = 0.2;

  xyze(0, 2) = -0.5;
  xyze(1, 2) = 0.8;
  xyze(2, 2) = 0.8;

  xyze(0, 3) = -0.5;
  xyze(1, 3) = 0.2;
  xyze(2, 3) = 0.8;

  xyze(0, 4) = 1.5;
  xyze(1, 4) = 0.2;
  xyze(2, 4) = 0.2;

  xyze(0, 5) = 1.5;
  xyze(1, 5) = 0.8;
  xyze(2, 5) = 0.2;

  xyze(0, 6) = 1.5;
  xyze(1, 6) = 0.8;
  xyze(2, 6) = 0.8;

  xyze(0, 7) = 1.5;
  xyze(1, 7) = 0.2;
  xyze(2, 7) = 0.8;

  w.CreateHex8Sides(xyze);

#if 1
  // add second cut to be able to find nodal positions

  xyze(0, 0) = -0.1;
  xyze(1, 0) = 0.1;
  xyze(2, 0) = -0.1;

  xyze(0, 1) = 1.1;
  xyze(1, 1) = 0.1;
  xyze(2, 1) = -0.1;

  xyze(0, 2) = 1.1;
  xyze(1, 2) = -0.1;
  xyze(2, 2) = 0.1;

  xyze(0, 3) = -0.1;
  xyze(1, 3) = -0.1;
  xyze(2, 3) = 0.1;

  w.CreateQuad4(xyze);
#endif

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_hex8_onside()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = 0.5;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = 0.2;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0.1;
  xyze(2, 1) = 0.2;

  xyze(0, 2) = 1.5;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = 0.2;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0.9;
  xyze(2, 3) = 0.2;

  xyze(0, 4) = 0.5;
  xyze(1, 4) = 0.5;
  xyze(2, 4) = 0.8;

  xyze(0, 5) = 1;
  xyze(1, 5) = 0.1;
  xyze(2, 5) = 0.8;

  xyze(0, 6) = 1.5;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = 0.8;

  xyze(0, 7) = 1;
  xyze(1, 7) = 0.9;
  xyze(2, 7) = 0.8;

  w.CreateHex8Sides(xyze);

  // add second cut to be able to find nodal positions

  xyze(0, 0) = -0.1;
  xyze(1, 0) = 0.1;
  xyze(2, 0) = -0.1;

  xyze(0, 1) = 1.1;
  xyze(1, 1) = 0.1;
  xyze(2, 1) = -0.1;

  xyze(0, 2) = 1.1;
  xyze(1, 2) = -0.1;
  xyze(2, 2) = 0.1;

  xyze(0, 3) = -0.1;
  xyze(1, 3) = -0.1;
  xyze(2, 3) = 0.1;


  w.CreateQuad4(xyze);
  w.Status();
  w.CutTest_Cut();
}

void test_hex8_hex8_internal()
{
  SimpleWrapper w;

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = -1;
  xyze(1, 0) = -1;
  xyze(2, 0) = -1;

  xyze(0, 1) = 1;
  xyze(1, 1) = -1;
  xyze(2, 1) = -1;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = -1;

  xyze(0, 3) = -1;
  xyze(1, 3) = 1;
  xyze(2, 3) = -1;

  xyze(0, 4) = -1;
  xyze(1, 4) = -1;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = -1;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = -1;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  w.CreateHex8(xyze);

  xyze(0, 0) = -1.5;
  xyze(1, 0) = -0.5;
  xyze(2, 0) = 0.707107;

  xyze(0, 1) = -0.5;
  xyze(1, 1) = -1.5;
  xyze(2, 1) = -0.707107;

  xyze(0, 2) = -0.207107;
  xyze(1, 2) = 0.207107;
  xyze(2, 2) = -1.70711;

  xyze(0, 3) = -1.20711;
  xyze(1, 3) = 1.20711;
  xyze(2, 3) = -0.292893;

  xyze(0, 4) = 0.207107;
  xyze(1, 4) = -0.207107;
  xyze(2, 4) = 1.70711;

  xyze(0, 5) = 1.20711;
  xyze(1, 5) = -1.20711;
  xyze(2, 5) = 0.292893;

  xyze(0, 6) = 1.5;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = -0.707107;

  xyze(0, 7) = 0.5;
  xyze(1, 7) = 1.5;
  xyze(2, 7) = 0.707107;

  w.CreateHex8Sides(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_hex8_sideintersection()
{
  SimpleWrapper w;

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = -1;
  xyze(1, 0) = -1;
  xyze(2, 0) = -1;

  xyze(0, 1) = 1;
  xyze(1, 1) = -1;
  xyze(2, 1) = -1;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = -1;

  xyze(0, 3) = -1;
  xyze(1, 3) = 1;
  xyze(2, 3) = -1;

  xyze(0, 4) = -1;
  xyze(1, 4) = -1;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = -1;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = -1;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  w.CreateHex8(xyze);

  xyze(0, 0) = 0.5;
  xyze(1, 0) = -0.5;
  xyze(2, 0) = -0.5;

  xyze(0, 1) = 1.5;
  xyze(1, 1) = -0.5;
  xyze(2, 1) = -0.5;

  xyze(0, 2) = 1.5;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = -0.5;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = -0.5;

  xyze(0, 4) = 0.5;
  xyze(1, 4) = -0.5;
  xyze(2, 4) = 0.5;

  xyze(0, 5) = 1.5;
  xyze(1, 5) = -0.5;
  xyze(2, 5) = 0.5;

  xyze(0, 6) = 1.5;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = 0.5;

  xyze(0, 7) = 0.5;
  xyze(1, 7) = 0.5;
  xyze(2, 7) = 0.5;

  w.CreateHex8Sides(xyze);

  // add second cut to be able to find nodal positions

  xyze(0, 0) = -1.1;
  xyze(1, 0) = -0.9;
  xyze(2, 0) = -1.1;

  xyze(0, 1) = 1.1;
  xyze(1, 1) = -0.9;
  xyze(2, 1) = -1.1;

  xyze(0, 2) = 1.1;
  xyze(1, 2) = -1.1;
  xyze(2, 2) = -0.9;

  xyze(0, 3) = -1.1;
  xyze(1, 3) = -1.1;
  xyze(2, 3) = -0.9;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_hex8_inside()
{
  SimpleWrapper w;

  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = -1;
  xyze(1, 0) = -1;
  xyze(2, 0) = -1;

  xyze(0, 1) = 1;
  xyze(1, 1) = -1;
  xyze(2, 1) = -1;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = -1;

  xyze(0, 3) = -1;
  xyze(1, 3) = 1;
  xyze(2, 3) = -1;

  xyze(0, 4) = -1;
  xyze(1, 4) = -1;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = -1;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = -1;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  w.CreateHex8(xyze);

  xyze(0, 0) = -0.5;
  xyze(1, 0) = -0.5;
  xyze(2, 0) = -0.5;

  xyze(0, 1) = 0.5;
  xyze(1, 1) = -0.5;
  xyze(2, 1) = -0.5;

  xyze(0, 2) = 0.5;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = -0.5;

  xyze(0, 3) = -0.5;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = -0.5;

  xyze(0, 4) = -0.5;
  xyze(1, 4) = -0.5;
  xyze(2, 4) = 0.5;

  xyze(0, 5) = 0.5;
  xyze(1, 5) = -0.5;
  xyze(2, 5) = 0.5;

  xyze(0, 6) = 0.5;
  xyze(1, 6) = 0.5;
  xyze(2, 6) = 0.5;

  xyze(0, 7) = -0.5;
  xyze(1, 7) = 0.5;
  xyze(2, 7) = 0.5;

  w.CreateHex8Sides(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_schnitt()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 0.5;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = -0.2;

  xyze(0, 1) = 1.5;
  xyze(1, 1) = 0.5;
  xyze(2, 1) = -0.2;

  xyze(0, 2) = 1.5;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = 1.2;

  xyze(0, 3) = 0.5;
  xyze(1, 3) = 0.5;
  xyze(2, 3) = 1.2;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_touch4()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 0.2;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1.5;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1.5;
  xyze(1, 2) = 0;
  xyze(2, 2) = 1.2;

  xyze(0, 3) = 0.2;
  xyze(1, 3) = 0;
  xyze(2, 3) = 1.2;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_touch5()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 0.2;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1.5;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1.5;
  xyze(1, 2) = 0;
  xyze(2, 2) = 1.2;

  xyze(0, 3) = 1.2;
  xyze(1, 3) = 0;
  xyze(2, 3) = 1.2;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_touch6()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 1;
  xyze(2, 1) = 0.5;

  xyze(0, 2) = 1;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = 1;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0;
  xyze(2, 3) = 0.5;

  w.CreateQuad4(xyze);

  // add second cut to be able to find nodal positions

  xyze(0, 0) = 0.1;
  xyze(1, 0) = -0.1;
  xyze(2, 0) = -0.1;

  xyze(0, 1) = 0.1;
  xyze(1, 1) = 1.1;
  xyze(2, 1) = -0.1;

  xyze(0, 2) = -0.1;
  xyze(1, 2) = 1.1;
  xyze(2, 2) = 0.1;

  xyze(0, 3) = -0.1;
  xyze(1, 3) = -0.1;
  xyze(2, 3) = 0.1;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

void test_hex8_quad4_touch7()
{
  SimpleWrapper w;

  w.CreateHex8();

  Epetra_SerialDenseMatrix xyze(3, 4);

  xyze(0, 0) = 1;
  xyze(1, 0) = 0.5;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0.8;
  xyze(2, 1) = 0.5;

  xyze(0, 2) = 1;
  xyze(1, 2) = 0.5;
  xyze(2, 2) = 1;

  xyze(0, 3) = 1;
  xyze(1, 3) = 0.2;
  xyze(2, 3) = 0.5;

  w.CreateQuad4(xyze);

  w.Status();
  w.CutTest_Cut();
}

#if 0
void test_quad4_quad4_simple()
{
  GEO::CUT::Mesh mesh;

  Epetra_SerialDenseMatrix xyze( 3, 8 );

  xyze( 0, 0 ) = 0;
  xyze( 1, 0 ) = 0;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1;
  xyze( 1, 1 ) = 0;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1;
  xyze( 1, 2 ) = 1;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0;
  xyze( 1, 3 ) = 1;
  xyze( 2, 3 ) = 0;

  GEO::CUT::Side* quad4_1 = create_quad4( mesh, xyze );

  xyze( 0, 0 ) = 0.5;
  xyze( 1, 0 ) = 0.5;
  xyze( 2, 0 ) = 0;

  xyze( 0, 1 ) = 1.5;
  xyze( 1, 1 ) = 0.5;
  xyze( 2, 1 ) = 0;

  xyze( 0, 2 ) = 1.5;
  xyze( 1, 2 ) = 1.5;
  xyze( 2, 2 ) = 0;

  xyze( 0, 3 ) = 0.5;
  xyze( 1, 3 ) = 1.5;
  xyze( 2, 3 ) = 0;

  GEO::CUT::Side* quad4_2 = create_quad4( mesh, xyze );

  mesh.Status();

  quad4_1->Cut( mesh, *quad4_2, NULL );

  mesh.Status();
}
#endif

void test_hex8_quad4_mesh()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  GEO::CUT::Mesh mesh(options);

  create_hex8_mesh(mesh, 2, 2, 2);

  std::vector<GEO::CUT::Side*> sides;
  create_quad4_mesh(mesh, 3, 3, sides);

  mesh.Status();

  for (std::vector<GEO::CUT::Side*>::iterator i = sides.begin(); i != sides.end(); ++i)
  {
    GEO::CUT::Side* quad4 = *i;
    GEO::CUT::plain_element_set done;
    GEO::CUT::plain_element_set elements_done;
    mesh.Cut(*(quad4), done, elements_done);
  }

  cutmesh(mesh);
}

void test_position2d()
{
  LINALG::Matrix<3, 3> side_xyze;
  LINALG::Matrix<3, 1> xyz;
  LINALG::Matrix<3, 1> shift;

  side_xyze(0, 0) = -0.20710678118654757;
  side_xyze(1, 0) = 0;
  side_xyze(2, 0) = 0.62132034355964261;
  side_xyze(0, 1) = -0.20710678118654757;
  side_xyze(1, 1) = 0;
  side_xyze(2, 1) = -0.62132034355964261;
  side_xyze(0, 2) = 0.41421356237309503;
  side_xyze(1, 2) = 0;
  side_xyze(2, 2) = 0;

  xyz(0) = -0.20710678118654757;
  xyz(1) = -0.62132046378538341;
  xyz(2) = -0.62132034355964261;

  shift(0) = -0.41421356237309503;
  shift(1) = 1.2022574075492253e-07;
  shift(2) = 0;

  for (int i = 0; i < 3; ++i)
  {
    LINALG::Matrix<3, 1> x1(&side_xyze(0, i), true);
    x1.Update(1, shift, 1);
  }
  xyz.Update(1, shift, 1);

  double scale = 1.6094757082487299;

  side_xyze.Scale(scale);
  xyz.Scale(scale);

  GEO::CUT::PositionFactory::SpecifyGeneralDistFloattype(INPAR::CUT::floattype_cln);    // use cln
  GEO::CUT::PositionFactory::SpecifyGeneralPosFloattype(INPAR::CUT::floattype_double);  // use
                                                                                        // double
  Teuchos::RCP<GEO::CUT::Position> pos =
      GEO::CUT::Position::Create(side_xyze, xyz, DRT::Element::tri3);
  pos->Compute();
}
