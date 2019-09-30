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

void test_axel9()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950328;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.427257;
    tri3_xyze(0, 1) = 0.905384;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.438462;
    tri3_xyze(0, 2) = 0.935517;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.463586;
    nids.clear();
    nids.push_back(431);
    nids.push_back(432);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.905384;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.438462;
    tri3_xyze(0, 1) = 0.920706;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.499914;
    tri3_xyze(0, 2) = 0.935517;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.463586;
    nids.clear();
    nids.push_back(432);
    nids.push_back(433);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.499914;
    tri3_xyze(0, 1) = 0.965649;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.488709;
    tri3_xyze(0, 2) = 0.935517;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.463586;
    nids.clear();
    nids.push_back(433);
    nids.push_back(434);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.488709;
    tri3_xyze(0, 1) = 0.950328;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.427257;
    tri3_xyze(0, 2) = 0.935517;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.463586;
    nids.clear();
    nids.push_back(434);
    nids.push_back(431);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950328;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.427257;
    tri3_xyze(0, 1) = 0.952148;
    tri3_xyze(1, 1) = 0.69277;
    tri3_xyze(2, 1) = 0.368503;
    tri3_xyze(0, 2) = 0.928766;
    tri3_xyze(1, 2) = 0.679742;
    tri3_xyze(2, 2) = 0.403483;
    nids.clear();
    nids.push_back(431);
    nids.push_back(436);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.907205;
    tri3_xyze(1, 0) = 0.649577;
    tri3_xyze(2, 0) = 0.379708;
    tri3_xyze(0, 1) = 0.905384;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.438462;
    tri3_xyze(0, 2) = 0.928766;
    tri3_xyze(1, 2) = 0.679742;
    tri3_xyze(2, 2) = 0.403483;
    nids.clear();
    nids.push_back(437);
    nids.push_back(432);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.905384;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.438462;
    tri3_xyze(0, 1) = 0.950328;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.427257;
    tri3_xyze(0, 2) = 0.928766;
    tri3_xyze(1, 2) = 0.679742;
    tri3_xyze(2, 2) = 0.403483;
    nids.clear();
    nids.push_back(432);
    nids.push_back(431);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.488709;
    tri3_xyze(0, 1) = 1.0213;
    tri3_xyze(1, 1) = 0.73135;
    tri3_xyze(2, 1) = 0.474834;
    tri3_xyze(0, 2) = 0.985813;
    tri3_xyze(1, 2) = 0.720628;
    tri3_xyze(2, 2) = 0.451046;
    nids.clear();
    nids.push_back(434);
    nids.push_back(406);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.0213;
    tri3_xyze(1, 0) = 0.73135;
    tri3_xyze(2, 0) = 0.474834;
    tri3_xyze(0, 1) = 1.00598;
    tri3_xyze(1, 1) = 0.73135;
    tri3_xyze(2, 1) = 0.413382;
    tri3_xyze(0, 2) = 0.985813;
    tri3_xyze(1, 2) = 0.720628;
    tri3_xyze(2, 2) = 0.451046;
    nids.clear();
    nids.push_back(406);
    nids.push_back(407);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.00598;
    tri3_xyze(1, 0) = 0.73135;
    tri3_xyze(2, 0) = 0.413382;
    tri3_xyze(0, 1) = 0.950328;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.427257;
    tri3_xyze(0, 2) = 0.985813;
    tri3_xyze(1, 2) = 0.720628;
    tri3_xyze(2, 2) = 0.451046;
    nids.clear();
    nids.push_back(407);
    nids.push_back(431);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950328;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.427257;
    tri3_xyze(0, 1) = 0.965649;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.488709;
    tri3_xyze(0, 2) = 0.985813;
    tri3_xyze(1, 2) = 0.720628;
    tri3_xyze(2, 2) = 0.451046;
    nids.clear();
    nids.push_back(431);
    nids.push_back(434);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950328;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.427257;
    tri3_xyze(0, 1) = 1.00598;
    tri3_xyze(1, 1) = 0.73135;
    tri3_xyze(2, 1) = 0.413382;
    tri3_xyze(0, 2) = 0.975689;
    tri3_xyze(1, 2) = 0.712676;
    tri3_xyze(2, 2) = 0.394458;
    nids.clear();
    nids.push_back(431);
    nids.push_back(407);
    nids.push_back(440);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.952148;
    tri3_xyze(1, 0) = 0.69277;
    tri3_xyze(2, 0) = 0.368503;
    tri3_xyze(0, 1) = 0.950328;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.427257;
    tri3_xyze(0, 2) = 0.975689;
    tri3_xyze(1, 2) = 0.712676;
    tri3_xyze(2, 2) = 0.394458;
    nids.clear();
    nids.push_back(436);
    nids.push_back(431);
    nids.push_back(440);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.875763;
    tri3_xyze(1, 0) = 0.62352;
    tri3_xyze(2, 0) = 0.51112;
    tri3_xyze(0, 1) = 0.920706;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.499914;
    tri3_xyze(0, 2) = 0.890574;
    tri3_xyze(1, 2) = 0.645117;
    tri3_xyze(2, 2) = 0.474791;
    nids.clear();
    nids.push_back(442);
    nids.push_back(433);
    nids.push_back(443);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.499914;
    tri3_xyze(0, 1) = 0.905384;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.438462;
    tri3_xyze(0, 2) = 0.890574;
    tri3_xyze(1, 2) = 0.645117;
    tri3_xyze(2, 2) = 0.474791;
    nids.clear();
    nids.push_back(433);
    nids.push_back(432);
    nids.push_back(443);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.488709;
    tri3_xyze(0, 1) = 0.920706;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.499914;
    tri3_xyze(0, 2) = 0.950839;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.525038;
    nids.clear();
    nids.push_back(434);
    nids.push_back(433);
    nids.push_back(448);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.499914;
    tri3_xyze(0, 1) = 0.936028;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.561367;
    tri3_xyze(0, 2) = 0.950839;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.525038;
    nids.clear();
    nids.push_back(433);
    nids.push_back(446);
    nids.push_back(448);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.980971;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.550161;
    tri3_xyze(0, 1) = 0.965649;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.488709;
    tri3_xyze(0, 2) = 0.950839;
    tri3_xyze(1, 2) = 0.68831;
    tri3_xyze(2, 2) = 0.525038;
    nids.clear();
    nids.push_back(447);
    nids.push_back(434);
    nids.push_back(448);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.0213;
    tri3_xyze(1, 0) = 0.73135;
    tri3_xyze(2, 0) = 0.474834;
    tri3_xyze(0, 1) = 0.965649;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.488709;
    tri3_xyze(0, 2) = 1.00113;
    tri3_xyze(1, 2) = 0.720628;
    tri3_xyze(2, 2) = 0.512498;
    nids.clear();
    nids.push_back(406);
    nids.push_back(434);
    nids.push_back(449);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649;
    tri3_xyze(1, 0) = 0.709907;
    tri3_xyze(2, 0) = 0.488709;
    tri3_xyze(0, 1) = 0.980971;
    tri3_xyze(1, 1) = 0.709907;
    tri3_xyze(2, 1) = 0.550161;
    tri3_xyze(0, 2) = 1.00113;
    tri3_xyze(1, 2) = 0.720628;
    tri3_xyze(2, 2) = 0.512498;
    nids.clear();
    nids.push_back(434);
    nids.push_back(447);
    nids.push_back(449);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.499914;
    tri3_xyze(0, 1) = 0.875763;
    tri3_xyze(1, 1) = 0.62352;
    tri3_xyze(2, 1) = 0.51112;
    tri3_xyze(0, 2) = 0.905895;
    tri3_xyze(1, 2) = 0.645117;
    tri3_xyze(2, 2) = 0.536243;
    nids.clear();
    nids.push_back(433);
    nids.push_back(442);
    nids.push_back(451);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.936028;
    tri3_xyze(1, 0) = 0.666714;
    tri3_xyze(2, 0) = 0.561367;
    tri3_xyze(0, 1) = 0.920706;
    tri3_xyze(1, 1) = 0.666714;
    tri3_xyze(2, 1) = 0.499914;
    tri3_xyze(0, 2) = 0.905895;
    tri3_xyze(1, 2) = 0.645117;
    tri3_xyze(2, 2) = 0.536243;
    nids.clear();
    nids.push_back(446);
    nids.push_back(433);
    nids.push_back(451);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.916667;
  hex8_xyze(1, 0) = 0.666667;
  hex8_xyze(2, 0) = 0.5;
  hex8_xyze(0, 1) = 0.916667;
  hex8_xyze(1, 1) = 0.666667;
  hex8_xyze(2, 1) = 0.416667;
  hex8_xyze(0, 2) = 0.916667;
  hex8_xyze(1, 2) = 0.75;
  hex8_xyze(2, 2) = 0.416667;
  hex8_xyze(0, 3) = 0.916667;
  hex8_xyze(1, 3) = 0.75;
  hex8_xyze(2, 3) = 0.5;
  hex8_xyze(0, 4) = 1;
  hex8_xyze(1, 4) = 0.666667;
  hex8_xyze(2, 4) = 0.5;
  hex8_xyze(0, 5) = 1;
  hex8_xyze(1, 5) = 0.666667;
  hex8_xyze(2, 5) = 0.416667;
  hex8_xyze(0, 6) = 1;
  hex8_xyze(1, 6) = 0.75;
  hex8_xyze(2, 6) = 0.416667;
  hex8_xyze(0, 7) = 1;
  hex8_xyze(1, 7) = 0.75;
  hex8_xyze(2, 7) = 0.5;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
