/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

\maintainer Martin Kronbichler
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

void test_alex42()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.792e-01;
    tri3_xyze(1, 0) = 2.604e-01;
    tri3_xyze(2, 0) = -1.000e-04;
    tri3_xyze(0, 1) = 9.906e-01;
    tri3_xyze(1, 1) = 2.897e-01;
    tri3_xyze(2, 1) = -1.000e-04;
    tri3_xyze(0, 2) = 9.849e-01;
    tri3_xyze(1, 2) = 2.750e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(494);
    nids.push_back(498);
    nids.push_back(716);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.906e-01;
    tri3_xyze(1, 0) = 2.897e-01;
    tri3_xyze(2, 0) = -1.000e-04;
    tri3_xyze(0, 1) = 9.906e-01;
    tri3_xyze(1, 1) = 2.896e-01;
    tri3_xyze(2, 1) = 3.196e-02;
    tri3_xyze(0, 2) = 9.849e-01;
    tri3_xyze(1, 2) = 2.750e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(498);
    nids.push_back(714);
    nids.push_back(716);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.906e-01;
    tri3_xyze(1, 0) = 2.896e-01;
    tri3_xyze(2, 0) = 3.196e-02;
    tri3_xyze(0, 1) = 9.792e-01;
    tri3_xyze(1, 1) = 2.603e-01;
    tri3_xyze(2, 1) = 3.196e-02;
    tri3_xyze(0, 2) = 9.849e-01;
    tri3_xyze(1, 2) = 2.750e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(714);
    nids.push_back(695);
    nids.push_back(716);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.792e-01;
    tri3_xyze(1, 0) = 2.603e-01;
    tri3_xyze(2, 0) = 3.196e-02;
    tri3_xyze(0, 1) = 9.792e-01;
    tri3_xyze(1, 1) = 2.604e-01;
    tri3_xyze(2, 1) = -1.000e-04;
    tri3_xyze(0, 2) = 9.849e-01;
    tri3_xyze(1, 2) = 2.750e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(695);
    nids.push_back(494);
    nids.push_back(716);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.906e-01;
    tri3_xyze(1, 0) = 2.896e-01;
    tri3_xyze(2, 0) = 3.196e-02;
    tri3_xyze(0, 1) = 9.906e-01;
    tri3_xyze(1, 1) = 2.897e-01;
    tri3_xyze(2, 1) = -1.000e-04;
    tri3_xyze(0, 2) = 9.964e-01;
    tri3_xyze(1, 2) = 3.044e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(714);
    nids.push_back(498);
    nids.push_back(735);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.002e+00;
    tri3_xyze(1, 0) = 3.190e-01;
    tri3_xyze(2, 0) = 3.195e-02;
    tri3_xyze(0, 1) = 9.906e-01;
    tri3_xyze(1, 1) = 2.896e-01;
    tri3_xyze(2, 1) = 3.196e-02;
    tri3_xyze(0, 2) = 9.964e-01;
    tri3_xyze(1, 2) = 3.044e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(733);
    nids.push_back(714);
    nids.push_back(735);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 9.906e-01;
    tri3_xyze(1, 0) = 2.897e-01;
    tri3_xyze(2, 0) = -1.000e-04;
    tri3_xyze(0, 1) = 1.002e+00;
    tri3_xyze(1, 1) = 3.191e-01;
    tri3_xyze(2, 1) = -1.000e-04;
    tri3_xyze(0, 2) = 9.964e-01;
    tri3_xyze(1, 2) = 3.044e-01;
    tri3_xyze(2, 2) = 1.593e-02;
    nids.clear();
    nids.push_back(498);
    nids.push_back(502);
    nids.push_back(735);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.000e+00;
  hex8_xyze(1, 0) = 3.000e-01;
  hex8_xyze(2, 0) = 2.500e-02;
  hex8_xyze(0, 1) = 1.000e+00;
  hex8_xyze(1, 1) = 2.750e-01;
  hex8_xyze(2, 1) = 2.500e-02;
  hex8_xyze(0, 2) = 9.750e-01;
  hex8_xyze(1, 2) = 2.750e-01;
  hex8_xyze(2, 2) = 2.500e-02;
  hex8_xyze(0, 3) = 9.750e-01;
  hex8_xyze(1, 3) = 3.000e-01;
  hex8_xyze(2, 3) = 2.500e-02;
  hex8_xyze(0, 4) = 1.000e+00;
  hex8_xyze(1, 4) = 3.000e-01;
  hex8_xyze(2, 4) = 0.000e+00;
  hex8_xyze(0, 5) = 1.000e+00;
  hex8_xyze(1, 5) = 2.750e-01;
  hex8_xyze(2, 5) = 0.000e+00;
  hex8_xyze(0, 6) = 9.750e-01;
  hex8_xyze(1, 6) = 2.750e-01;
  hex8_xyze(2, 6) = 0.000e+00;
  hex8_xyze(0, 7) = 9.750e-01;
  hex8_xyze(1, 7) = 3.000e-01;
  hex8_xyze(2, 7) = 0.000e+00;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);


  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
