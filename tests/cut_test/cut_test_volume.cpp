/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_volume.cpp

\level 1

\maintainer Ager Christoph
*----------------------------------------------------------------------*/

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_mesh.H"

void test_cut_volumes()
{
  GEO::CUT::Options options;
  options.Init_for_Cuttests();
  // this is meant to be used with matching boundaries. Thus, no
  // inside/outside positions.
  options.SetFindPositions(false);

  GEO::CUT::Mesh mesh1(options);
  GEO::CUT::Mesh mesh2(options, 1, mesh1.Points());

  create_hex8_mesh(mesh1, 4, 4, 4);
  create_hex8_mesh(mesh2, 3, 5, 2);

  mesh2.CreateSideIds_CutTest();

  mesh1.Status();
  mesh2.Status();

  GEO::CUT::plain_element_set elements_done;

  mesh2.Cut(mesh1, elements_done);

  cutmesh(mesh1);

  mesh2.AssignOtherVolumeCells_CutTest(mesh1);
}

void test_cut_volumes2()
{
  for (int i = 2; i < 5; ++i)
  {
    for (int j = 2; j < 5; ++j)
    {
      for (int k = 2; k < 5; ++k)
      {
        GEO::CUT::Options options;
        options.Init_for_Cuttests();
        // this is meant to be used with matching boundaries. Thus, no
        // inside/outside positions.
        options.SetFindPositions(false);

        GEO::CUT::Mesh mesh1(options);
        GEO::CUT::Mesh mesh2(options, 1, mesh1.Points());

#if 1
        create_hex8_mesh(mesh1, 1, 1, 1);
        create_hex8_mesh(mesh2, i, j, k);
#else
        create_hex8_mesh(mesh1, 4, 4, 4);
        create_hex8_mesh(mesh2, 8, 9, 10);
#endif

        mesh2.CreateSideIds_CutTest();

        mesh1.Status();
        mesh2.Status();

        GEO::CUT::plain_element_set elements_done;

        mesh2.Cut(mesh1, elements_done);

        cutmesh(mesh1);

        mesh2.AssignOtherVolumeCells_CutTest(mesh1);
      }
    }
  }
}

void test_cut_volumes3()
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

  xyze(0, 0) = 0;
  xyze(1, 0) = -1;
  xyze(2, 0) = -1;

  xyze(0, 1) = 0;
  xyze(1, 1) = -1;
  xyze(2, 1) = 1;

  xyze(0, 2) = 0;
  xyze(1, 2) = 1;
  xyze(2, 2) = 1;

  xyze(0, 3) = 0;
  xyze(1, 3) = 1;
  xyze(2, 3) = -1;

  xyze(0, 4) = -0.5;
  xyze(1, 4) = 0;
  xyze(2, 4) = 0;

  w.CreatePyramid5Sides(xyze);

#if 0
  xyze( 0, 0 ) =  0;
  xyze( 1, 0 ) = -1;
  xyze( 2, 0 ) = -1;

  xyze( 0, 1 ) =  0;
  xyze( 1, 1 ) =  1;
  xyze( 2, 1 ) = -1;

  xyze( 0, 2 ) =  0;
  xyze( 1, 2 ) =  1;
  xyze( 2, 2 ) =  1;

  xyze( 0, 3 ) =  0;
  xyze( 1, 3 ) = -1;
  xyze( 2, 3 ) =  1;

  xyze( 0, 4 ) = 0.5;
  xyze( 1, 4 ) = 0;
  xyze( 2, 4 ) = 0;

  w.CreatePyramid5Sides( xyze );
#endif

  w.Status();
  w.CutTest_Cut();
}
