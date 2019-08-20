/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

\maintainer Christoph Ager
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

void test_alex43()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.157270e-01;
    tri3_xyze(1, 0) = 5.772857e-02;
    tri3_xyze(2, 0) = 1.925566e-01;
    tri3_xyze(0, 1) = 9.156979e-01;
    tri3_xyze(1, 1) = 5.773074e-02;
    tri3_xyze(2, 1) = 1.603983e-01;
    tri3_xyze(0, 2) = 9.182831e-01;
    tri3_xyze(1, 2) = 7.231715e-02;
    tri3_xyze(2, 2) = 1.764784e-01;
    nids.clear();
    nids.push_back(541);
    nids.push_back(545);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.156979e-01;
    tri3_xyze(1, 0) = 5.773074e-02;
    tri3_xyze(2, 0) = 1.603983e-01;
    tri3_xyze(0, 1) = 9.123871e-01;
    tri3_xyze(1, 1) = 2.870955e-02;
    tri3_xyze(2, 1) = 1.603993e-01;
    tri3_xyze(0, 2) = 9.140278e-01;
    tri3_xyze(1, 2) = 4.322599e-02;
    tri3_xyze(2, 2) = 1.443401e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(543);
    nids.push_back(550);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.156979e-01;
    tri3_xyze(1, 0) = 5.773074e-02;
    tri3_xyze(2, 0) = 1.603983e-01;
    tri3_xyze(0, 1) = 9.156601e-01;
    tri3_xyze(1, 1) = 5.774456e-02;
    tri3_xyze(2, 1) = 1.282767e-01;
    tri3_xyze(0, 2) = 9.182471e-01;
    tri3_xyze(1, 2) = 7.232722e-02;
    tri3_xyze(2, 2) = 1.443339e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(549);
    nids.push_back(576);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.208392e-01;
    tri3_xyze(1, 0) = 8.690832e-02;
    tri3_xyze(2, 0) = 1.603939e-01;
    tri3_xyze(0, 1) = 9.156979e-01;
    tri3_xyze(1, 1) = 5.773074e-02;
    tri3_xyze(2, 1) = 1.603983e-01;
    tri3_xyze(0, 2) = 9.182471e-01;
    tri3_xyze(1, 2) = 7.232722e-02;
    tri3_xyze(2, 2) = 1.443339e-01;
    nids.clear();
    nids.push_back(573);
    nids.push_back(545);
    nids.push_back(576);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.123871e-01;
    tri3_xyze(1, 0) = 2.870955e-02;
    tri3_xyze(2, 0) = 1.603993e-01;
    tri3_xyze(0, 1) = 9.123661e-01;
    tri3_xyze(1, 1) = 2.871913e-02;
    tri3_xyze(2, 1) = 1.282860e-01;
    tri3_xyze(0, 2) = 9.140278e-01;
    tri3_xyze(1, 2) = 4.322599e-02;
    tri3_xyze(2, 2) = 1.443401e-01;
    nids.clear();
    nids.push_back(543);
    nids.push_back(547);
    nids.push_back(550);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.156601e-01;
    tri3_xyze(1, 0) = 5.774456e-02;
    tri3_xyze(2, 0) = 1.282767e-01;
    tri3_xyze(0, 1) = 9.156979e-01;
    tri3_xyze(1, 1) = 5.773074e-02;
    tri3_xyze(2, 1) = 1.603983e-01;
    tri3_xyze(0, 2) = 9.140278e-01;
    tri3_xyze(1, 2) = 4.322599e-02;
    tri3_xyze(2, 2) = 1.443401e-01;
    nids.clear();
    nids.push_back(549);
    nids.push_back(545);
    nids.push_back(550);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.123871e-01;
    tri3_xyze(1, 0) = 2.870955e-02;
    tri3_xyze(2, 0) = 1.603993e-01;
    tri3_xyze(0, 1) = 9.156979e-01;
    tri3_xyze(1, 1) = 5.773074e-02;
    tri3_xyze(2, 1) = 1.603983e-01;
    tri3_xyze(0, 2) = 9.140544e-01;
    tri3_xyze(1, 2) = 4.321917e-02;
    tri3_xyze(2, 2) = 1.764726e-01;
    nids.clear();
    nids.push_back(543);
    nids.push_back(545);
    nids.push_back(546);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.124055e-01;
    tri3_xyze(1, 0) = 2.870781e-02;
    tri3_xyze(2, 0) = 1.925361e-01;
    tri3_xyze(0, 1) = 9.123871e-01;
    tri3_xyze(1, 1) = 2.870955e-02;
    tri3_xyze(2, 1) = 1.603993e-01;
    tri3_xyze(0, 2) = 9.140544e-01;
    tri3_xyze(1, 2) = 4.321917e-02;
    tri3_xyze(2, 2) = 1.764726e-01;
    nids.clear();
    nids.push_back(539);
    nids.push_back(543);
    nids.push_back(546);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.156979e-01;
    tri3_xyze(1, 0) = 5.773074e-02;
    tri3_xyze(2, 0) = 1.603983e-01;
    tri3_xyze(0, 1) = 9.208392e-01;
    tri3_xyze(1, 1) = 8.690832e-02;
    tri3_xyze(2, 1) = 1.603939e-01;
    tri3_xyze(0, 2) = 9.182831e-01;
    tri3_xyze(1, 2) = 7.231715e-02;
    tri3_xyze(2, 2) = 1.764784e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(573);
    nids.push_back(574);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.156979e-01;
    tri3_xyze(1, 0) = 5.773074e-02;
    tri3_xyze(2, 0) = 1.603983e-01;
    tri3_xyze(0, 1) = 9.157270e-01;
    tri3_xyze(1, 1) = 5.772857e-02;
    tri3_xyze(2, 1) = 1.925566e-01;
    tri3_xyze(0, 2) = 9.140544e-01;
    tri3_xyze(1, 2) = 4.321917e-02;
    tri3_xyze(2, 2) = 1.764726e-01;
    nids.clear();
    nids.push_back(545);
    nids.push_back(541);
    nids.push_back(546);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 9.253731e-01;
  hex8_xyze(1, 0) = 2.962963e-02;
  hex8_xyze(2, 0) = 1.470588e-01;
  hex8_xyze(0, 1) = 9.253731e-01;
  hex8_xyze(1, 1) = 5.925926e-02;
  hex8_xyze(2, 1) = 1.470588e-01;
  hex8_xyze(0, 2) = 8.955224e-01;
  hex8_xyze(1, 2) = 5.925926e-02;
  hex8_xyze(2, 2) = 1.470588e-01;
  hex8_xyze(0, 3) = 8.955224e-01;
  hex8_xyze(1, 3) = 2.962963e-02;
  hex8_xyze(2, 3) = 1.470588e-01;
  hex8_xyze(0, 4) = 9.253731e-01;
  hex8_xyze(1, 4) = 2.962963e-02;
  hex8_xyze(2, 4) = 1.764706e-01;
  hex8_xyze(0, 5) = 9.253731e-01;
  hex8_xyze(1, 5) = 5.925926e-02;
  hex8_xyze(2, 5) = 1.764706e-01;
  hex8_xyze(0, 6) = 8.955224e-01;
  hex8_xyze(1, 6) = 5.925926e-02;
  hex8_xyze(2, 6) = 1.764706e-01;
  hex8_xyze(0, 7) = 8.955224e-01;
  hex8_xyze(1, 7) = 2.962963e-02;
  hex8_xyze(2, 7) = 1.764706e-01;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);


  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
