/*!----------------------------------------------------------------------
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

void test_alex48()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.183851e-01;
    tri3_xyze(1, 0) = 7.454914e-01;
    tri3_xyze(2, 0) = 1.604224e-01;
    tri3_xyze(0, 1) = 9.134352e-01;
    tri3_xyze(1, 1) = 7.728957e-01;
    tri3_xyze(2, 1) = 1.604116e-01;
    tri3_xyze(0, 2) = 9.159278e-01;
    tri3_xyze(1, 2) = 7.591951e-01;
    tri3_xyze(2, 2) = 1.765066e-01;
    nids.clear();
    nids.push_back(1452);
    nids.push_back(1471);
    nids.push_back(1472);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.134352e-01;
    tri3_xyze(1, 0) = 7.728957e-01;
    tri3_xyze(2, 0) = 1.604116e-01;
    tri3_xyze(0, 1) = 9.134630e-01;
    tri3_xyze(1, 1) = 7.728987e-01;
    tri3_xyze(2, 1) = 1.925733e-01;
    tri3_xyze(0, 2) = 9.159278e-01;
    tri3_xyze(1, 2) = 7.591951e-01;
    tri3_xyze(2, 2) = 1.765066e-01;
    nids.clear();
    nids.push_back(1471);
    nids.push_back(1469);
    nids.push_back(1472);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.134630e-01;
    tri3_xyze(1, 0) = 7.728987e-01;
    tri3_xyze(2, 0) = 1.925733e-01;
    tri3_xyze(0, 1) = 9.134352e-01;
    tri3_xyze(1, 1) = 7.728957e-01;
    tri3_xyze(2, 1) = 1.604116e-01;
    tri3_xyze(0, 2) = 9.118745e-01;
    tri3_xyze(1, 2) = 7.864986e-01;
    tri3_xyze(2, 2) = 1.764712e-01;
    nids.clear();
    nids.push_back(1469);
    nids.push_back(1471);
    nids.push_back(1486);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.134352e-01;
    tri3_xyze(1, 0) = 7.728957e-01;
    tri3_xyze(2, 0) = 1.604116e-01;
    tri3_xyze(0, 1) = 9.103000e-01;
    tri3_xyze(1, 1) = 8.001000e-01;
    tri3_xyze(2, 1) = 1.604000e-01;
    tri3_xyze(0, 2) = 9.118745e-01;
    tri3_xyze(1, 2) = 7.864986e-01;
    tri3_xyze(2, 2) = 1.764712e-01;
    nids.clear();
    nids.push_back(1471);
    nids.push_back(1262);
    nids.push_back(1486);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.134030e-01;
    tri3_xyze(1, 0) = 7.728806e-01;
    tri3_xyze(2, 0) = 1.282873e-01;
    tri3_xyze(0, 1) = 9.134352e-01;
    tri3_xyze(1, 1) = 7.728957e-01;
    tri3_xyze(2, 1) = 1.604116e-01;
    tri3_xyze(0, 2) = 9.158883e-01;
    tri3_xyze(1, 2) = 7.591843e-01;
    tri3_xyze(2, 2) = 1.443511e-01;
    nids.clear();
    nids.push_back(1473);
    nids.push_back(1471);
    nids.push_back(1474);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.134352e-01;
    tri3_xyze(1, 0) = 7.728957e-01;
    tri3_xyze(2, 0) = 1.604116e-01;
    tri3_xyze(0, 1) = 9.183851e-01;
    tri3_xyze(1, 1) = 7.454914e-01;
    tri3_xyze(2, 1) = 1.604224e-01;
    tri3_xyze(0, 2) = 9.158883e-01;
    tri3_xyze(1, 2) = 7.591843e-01;
    tri3_xyze(2, 2) = 1.443511e-01;
    nids.clear();
    nids.push_back(1471);
    nids.push_back(1452);
    nids.push_back(1474);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.134352e-01;
    tri3_xyze(1, 0) = 7.728957e-01;
    tri3_xyze(2, 0) = 1.604116e-01;
    tri3_xyze(0, 1) = 9.134030e-01;
    tri3_xyze(1, 1) = 7.728806e-01;
    tri3_xyze(2, 1) = 1.282873e-01;
    tri3_xyze(0, 2) = 9.118595e-01;
    tri3_xyze(1, 2) = 7.864941e-01;
    tri3_xyze(2, 2) = 1.443497e-01;
    nids.clear();
    nids.push_back(1471);
    nids.push_back(1473);
    nids.push_back(1487);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.103000e-01;
    tri3_xyze(1, 0) = 8.001000e-01;
    tri3_xyze(2, 0) = 1.283000e-01;
    tri3_xyze(0, 1) = 9.103000e-01;
    tri3_xyze(1, 1) = 8.001000e-01;
    tri3_xyze(2, 1) = 1.604000e-01;
    tri3_xyze(0, 2) = 9.118595e-01;
    tri3_xyze(1, 2) = 7.864941e-01;
    tri3_xyze(2, 2) = 1.443497e-01;
    nids.clear();
    nids.push_back(1264);
    nids.push_back(1262);
    nids.push_back(1487);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.103000e-01;
    tri3_xyze(1, 0) = 8.001000e-01;
    tri3_xyze(2, 0) = 1.604000e-01;
    tri3_xyze(0, 1) = 9.134352e-01;
    tri3_xyze(1, 1) = 7.728957e-01;
    tri3_xyze(2, 1) = 1.604116e-01;
    tri3_xyze(0, 2) = 9.118595e-01;
    tri3_xyze(1, 2) = 7.864941e-01;
    tri3_xyze(2, 2) = 1.443497e-01;
    nids.clear();
    nids.push_back(1262);
    nids.push_back(1471);
    nids.push_back(1487);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.103000e-01;
    tri3_xyze(1, 0) = 8.001000e-01;
    tri3_xyze(2, 0) = 1.604000e-01;
    tri3_xyze(0, 1) = 9.103000e-01;
    tri3_xyze(1, 1) = 8.001000e-01;
    tri3_xyze(2, 1) = 1.925000e-01;
    tri3_xyze(0, 2) = 9.118745e-01;
    tri3_xyze(1, 2) = 7.864986e-01;
    tri3_xyze(2, 2) = 1.764712e-01;
    nids.clear();
    nids.push_back(1262);
    nids.push_back(1260);
    nids.push_back(1486);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 9.253731e-01;
  hex8_xyze(1, 0) = 7.703703e-01;
  hex8_xyze(2, 0) = 1.470588e-01;
  hex8_xyze(0, 1) = 9.253731e-01;
  hex8_xyze(1, 1) = 8.000000e-01;
  hex8_xyze(2, 1) = 1.470588e-01;
  hex8_xyze(0, 2) = 8.955224e-01;
  hex8_xyze(1, 2) = 8.000000e-01;
  hex8_xyze(2, 2) = 1.470588e-01;
  hex8_xyze(0, 3) = 8.955224e-01;
  hex8_xyze(1, 3) = 7.703703e-01;
  hex8_xyze(2, 3) = 1.470588e-01;
  hex8_xyze(0, 4) = 9.253731e-01;
  hex8_xyze(1, 4) = 7.703703e-01;
  hex8_xyze(2, 4) = 1.764706e-01;
  hex8_xyze(0, 5) = 9.253731e-01;
  hex8_xyze(1, 5) = 8.000000e-01;
  hex8_xyze(2, 5) = 1.764706e-01;
  hex8_xyze(0, 6) = 8.955224e-01;
  hex8_xyze(1, 6) = 8.000000e-01;
  hex8_xyze(2, 6) = 1.764706e-01;
  hex8_xyze(0, 7) = 8.955224e-01;
  hex8_xyze(1, 7) = 7.703703e-01;
  hex8_xyze(2, 7) = 1.764706e-01;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
