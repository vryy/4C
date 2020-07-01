/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <list>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_side.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_cut/cut_tetmeshintersection.H"
#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_volumecell.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_hex8quad4selfcut20()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(14);
    nids.push_back(13);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut21()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut22()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.0;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.0;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.0;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut23()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.0;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.0;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut24()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut30()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(14);
    nids.push_back(13);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(13);
    nids.push_back(12);
    nids.push_back(33);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut31()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut32()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut33()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut34()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut35()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut36()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut37()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut38()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.4;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.4;

    nids.clear();
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(27);
    nids.push_back(28);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut39()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.0;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.0;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.0;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.0;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(27);
    nids.push_back(28);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut41()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut42()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);

  std::vector<double> dirdivVol;

  GEO::CUT::Mesh mesh = intersection.NormalMesh();
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = mesh.VolumeCells();
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    dirdivVol.push_back(vc->Volume());
  }

  for (uint i = 0; i < dirdivVol.size(); ++i)
    std::cout << "volume is: " << dirdivVol[i] << std::endl;
}

void test_hex8quad4selfcut43()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -2.3;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -2.3;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut51()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(27);
    nids.push_back(28);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(27);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(28);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut52()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(27);
    nids.push_back(28);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(27);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(28);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut53()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -2.3;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -2.3;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(27);
    nids.push_back(28);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(27);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(28);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut61()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut62()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut63()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut64()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut65()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut66()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut71()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.7;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.7;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.7;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.7;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);

  std::vector<double> tessVol, momFitVol, dirDivVol;

  GEO::CUT::Mesh mesh = intersection.NormalMesh();
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = mesh.VolumeCells();
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    tessVol.push_back(vc->Volume());
  }

  intersection.Status();

  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    vc->MomentFitGaussWeights(
        vc->ParentElement(), mesh, true, INPAR::CUT::BCellGaussPts_Tessellation);
    momFitVol.push_back(vc->Volume());
  }

  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    vc->DirectDivergenceGaussRule(
        vc->ParentElement(), mesh, true, INPAR::CUT::BCellGaussPts_Tessellation);
    dirDivVol.push_back(vc->Volume());
  }

  std::cout << "the volumes predicted by\ntessellation \tMomentFitting \tDirectDivergence\n";
  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    std::cout << tessVol[i] << "\t\t" << momFitVol[i] << "\t\t" << dirDivVol[i] << "\n";
    if (fabs(tessVol[i] - momFitVol[i]) > 1e-9 || fabs(dirDivVol[i] - momFitVol[i]) > 1e-9)
      std::cout << "WARNING: volume predicted by either one of the method is wrong\n";
  }
}

void test_hex8quad4selfcut72()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.7;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.7;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.7;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.7;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut81()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.45;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.45;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.45;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.45;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut82()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.45;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.45;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.45;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.45;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut83()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }


  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut84()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }


  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut85()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }


  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut86()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }


  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut91()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(40);
    nids.push_back(39);
    nids.push_back(30);
    nids.push_back(41);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut92()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(40);
    nids.push_back(39);
    nids.push_back(30);
    nids.push_back(41);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}

/*------------------------------------------------------------------------------------*
 * hex8 background element is cut by two quad cut surfaces. An edge of one cut surface
 * coincides with an edge of another cut surface
 *------------------------------------------------------------------------------------*/
void test_hex8quad4alignedEdges()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.5;  // 0.8;
    quad4_xyze(1, 0) = -0.1;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = 0.5;  // 0.8;
    quad4_xyze(1, 1) = -0.1;
    quad4_xyze(2, 1) = -0.1;

    quad4_xyze(0, 2) = 0.5;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = -0.1;

    quad4_xyze(0, 3) = 0.5;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  {
    Epetra_SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.5;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = 0.5;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = -0.1;

    quad4_xyze(0, 2) = 0.5;  // 0.6;
    quad4_xyze(1, 2) = 1.1;
    quad4_xyze(2, 2) = -0.1;

    quad4_xyze(0, 3) = 0.5;  // 0.6;
    quad4_xyze(1, 3) = 1.1;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(15);
    nids.push_back(16);
    nids.push_back(17);
    nids.push_back(18);
    intersection.AddCutSide(++sidecount, nids, quad4_xyze, DRT::Element::quad4);
  }

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.0;
  hex8_xyze(1, 0) = 1.0;
  hex8_xyze(2, 0) = 1.0;
  hex8_xyze(0, 1) = 1.0;
  hex8_xyze(1, 1) = 0.0;
  hex8_xyze(2, 1) = 1.0;
  hex8_xyze(0, 2) = 0.0;
  hex8_xyze(1, 2) = 0.0;
  hex8_xyze(2, 2) = 1.0;
  hex8_xyze(0, 3) = 0.0;
  hex8_xyze(1, 3) = 1.0;
  hex8_xyze(2, 3) = 1.0;
  hex8_xyze(0, 4) = 1.0;
  hex8_xyze(1, 4) = 1.0;
  hex8_xyze(2, 4) = 0.0;
  hex8_xyze(0, 5) = 1.0;
  hex8_xyze(1, 5) = 0.0;
  hex8_xyze(2, 5) = 0.0;
  hex8_xyze(0, 6) = 0.0;
  hex8_xyze(1, 6) = 0.0;
  hex8_xyze(2, 6) = 0.0;
  hex8_xyze(0, 7) = 0.0;
  hex8_xyze(1, 7) = 1.0;
  hex8_xyze(2, 7) = 0.0;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();

  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
