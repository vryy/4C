/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_alex_50.cpp

\level 1

\maintainer Ager Christoph
*----------------------------------------------------------------------*/

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.H"

#include "../../src/drt_cut/cut_side.H"
#include "../../src/drt_cut/cut_meshintersection.H"
#include "../../src/drt_cut/cut_tetmeshintersection.H"
#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_volumecell.H"

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_alex50()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.154469e-01;
    tri3_xyze(1, 0) = 5.803694e-02;
    tri3_xyze(2, 0) = 1.925507e-01;
    tri3_xyze(0, 1) = 9.154385e-01;
    tri3_xyze(1, 1) = 5.805569e-02;
    tri3_xyze(2, 1) = 2.247699e-01;
    tri3_xyze(0, 2) = 9.138705e-01;
    tri3_xyze(1, 2) = 4.345965e-02;
    tri3_xyze(2, 2) = 2.086410e-01;
    nids.clear();
    nids.push_back(541);
    nids.push_back(537);
    nids.push_back(542);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.154469e-01;
    tri3_xyze(1, 0) = 5.803694e-02;
    tri3_xyze(2, 0) = 1.925507e-01;
    tri3_xyze(0, 1) = 9.122981e-01;
    tri3_xyze(1, 1) = 2.886343e-02;
    tri3_xyze(2, 1) = 1.925326e-01;
    tri3_xyze(0, 2) = 9.138608e-01;
    tri3_xyze(1, 2) = 4.345167e-02;
    tri3_xyze(2, 2) = 1.764694e-01;
    nids.clear();
    nids.push_back(541);
    nids.push_back(539);
    nids.push_back(546);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.154385e-01;
    tri3_xyze(1, 0) = 5.805569e-02;
    tri3_xyze(2, 0) = 2.247699e-01;
    tri3_xyze(0, 1) = 9.154469e-01;
    tri3_xyze(1, 1) = 5.803694e-02;
    tri3_xyze(2, 1) = 1.925507e-01;
    tri3_xyze(0, 2) = 9.178737e-01;
    tri3_xyze(1, 2) = 7.270872e-02;
    tri3_xyze(2, 2) = 2.086677e-01;
    nids.clear();
    nids.push_back(537);
    nids.push_back(541);
    nids.push_back(572);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.122981e-01;
    tri3_xyze(1, 0) = 2.886343e-02;
    tri3_xyze(2, 0) = 1.925326e-01;
    tri3_xyze(0, 1) = 9.154469e-01;
    tri3_xyze(1, 1) = 5.803694e-02;
    tri3_xyze(2, 1) = 1.925507e-01;
    tri3_xyze(0, 2) = 9.138705e-01;
    tri3_xyze(1, 2) = 4.345965e-02;
    tri3_xyze(2, 2) = 2.086410e-01;
    nids.clear();
    nids.push_back(539);
    nids.push_back(541);
    nids.push_back(542);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.154179e-01;
    tri3_xyze(1, 0) = 5.804047e-02;
    tri3_xyze(2, 0) = 1.603961e-01;
    tri3_xyze(0, 1) = 9.154469e-01;
    tri3_xyze(1, 1) = 5.803694e-02;
    tri3_xyze(2, 1) = 1.925507e-01;
    tri3_xyze(0, 2) = 9.138608e-01;
    tri3_xyze(1, 2) = 4.345167e-02;
    tri3_xyze(2, 2) = 1.764694e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(541);
    nids.push_back(546);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.122986e-01;
    tri3_xyze(1, 0) = 2.888255e-02;
    tri3_xyze(2, 0) = 2.247107e-01;
    tri3_xyze(0, 1) = 9.122981e-01;
    tri3_xyze(1, 1) = 2.886343e-02;
    tri3_xyze(2, 1) = 1.925326e-01;
    tri3_xyze(0, 2) = 9.138705e-01;
    tri3_xyze(1, 2) = 4.345965e-02;
    tri3_xyze(2, 2) = 2.086410e-01;
    nids.clear();
    nids.push_back(535);
    nids.push_back(539);
    nids.push_back(542);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.203177e-01;
    tri3_xyze(1, 0) = 8.736618e-02;
    tri3_xyze(2, 0) = 1.925568e-01;
    tri3_xyze(0, 1) = 9.154469e-01;
    tri3_xyze(1, 1) = 5.803694e-02;
    tri3_xyze(2, 1) = 1.925507e-01;
    tri3_xyze(0, 2) = 9.178674e-01;
    tri3_xyze(1, 2) = 7.270475e-02;
    tri3_xyze(2, 2) = 1.764735e-01;
    nids.clear();
    nids.push_back(571);
    nids.push_back(541);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.122981e-01;
    tri3_xyze(1, 0) = 2.886343e-02;
    tri3_xyze(2, 0) = 1.925326e-01;
    tri3_xyze(0, 1) = 9.122801e-01;
    tri3_xyze(1, 1) = 2.886584e-02;
    tri3_xyze(2, 1) = 1.603981e-01;
    tri3_xyze(0, 2) = 9.138608e-01;
    tri3_xyze(1, 2) = 4.345167e-02;
    tri3_xyze(2, 2) = 1.764694e-01;
    nids.clear();
    nids.push_back(539);
    nids.push_back(543);
    nids.push_back(546);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.154469e-01;
    tri3_xyze(1, 0) = 5.803694e-02;
    tri3_xyze(2, 0) = 1.925507e-01;
    tri3_xyze(0, 1) = 9.154179e-01;
    tri3_xyze(1, 1) = 5.804047e-02;
    tri3_xyze(2, 1) = 1.603961e-01;
    tri3_xyze(0, 2) = 9.178674e-01;
    tri3_xyze(1, 2) = 7.270475e-02;
    tri3_xyze(2, 2) = 1.764735e-01;
    nids.clear();
    nids.push_back(541);
    nids.push_back(545);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.154469e-01;
    tri3_xyze(1, 0) = 5.803694e-02;
    tri3_xyze(2, 0) = 1.925507e-01;
    tri3_xyze(0, 1) = 9.203177e-01;
    tri3_xyze(1, 1) = 8.736618e-02;
    tri3_xyze(2, 1) = 1.925568e-01;
    tri3_xyze(0, 2) = 9.178737e-01;
    tri3_xyze(1, 2) = 7.270872e-02;
    tri3_xyze(2, 2) = 2.086677e-01;
    nids.clear();
    nids.push_back(541);
    nids.push_back(571);
    nids.push_back(572);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 9.253731e-01;
  hex8_xyze(1, 0) = 2.962963e-02;
  hex8_xyze(2, 0) = 1.764706e-01;
  hex8_xyze(0, 1) = 9.253731e-01;
  hex8_xyze(1, 1) = 5.925926e-02;
  hex8_xyze(2, 1) = 1.764706e-01;
  hex8_xyze(0, 2) = 8.955224e-01;
  hex8_xyze(1, 2) = 5.925926e-02;
  hex8_xyze(2, 2) = 1.764706e-01;
  hex8_xyze(0, 3) = 8.955224e-01;
  hex8_xyze(1, 3) = 2.962963e-02;
  hex8_xyze(2, 3) = 1.764706e-01;
  hex8_xyze(0, 4) = 9.253731e-01;
  hex8_xyze(1, 4) = 2.962963e-02;
  hex8_xyze(2, 4) = 2.058824e-01;
  hex8_xyze(0, 5) = 9.253731e-01;
  hex8_xyze(1, 5) = 5.925926e-02;
  hex8_xyze(2, 5) = 2.058824e-01;
  hex8_xyze(0, 6) = 8.955224e-01;
  hex8_xyze(1, 6) = 5.925926e-02;
  hex8_xyze(2, 6) = 2.058824e-01;
  hex8_xyze(0, 7) = 8.955224e-01;
  hex8_xyze(1, 7) = 2.962963e-02;
  hex8_xyze(2, 7) = 2.058824e-01;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();

  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
