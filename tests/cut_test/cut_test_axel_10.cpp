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

void test_axel10()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950327694;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.427256823;
    tri3_xyze(0, 1) = 0.905384481;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.438462406;
    tri3_xyze(0, 2) = 0.9355169535;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.463585645;
    nids.clear();
    nids.push_back(431);
    nids.push_back(432);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.905384481;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.438462406;
    tri3_xyze(0, 1) = 0.920706213;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.499914467;
    tri3_xyze(0, 2) = 0.9355169535;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.463585645;
    nids.clear();
    nids.push_back(432);
    nids.push_back(433);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706213;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.499914467;
    tri3_xyze(0, 1) = 0.965649426;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.488708884;
    tri3_xyze(0, 2) = 0.9355169535;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.463585645;
    nids.clear();
    nids.push_back(433);
    nids.push_back(434);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649426;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.488708884;
    tri3_xyze(0, 1) = 0.950327694;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.427256823;
    tri3_xyze(0, 2) = 0.9355169535;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.463585645;
    nids.clear();
    nids.push_back(434);
    nids.push_back(431);
    nids.push_back(435);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950327694;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.427256823;
    tri3_xyze(0, 1) = 0.95214808;
    tri3_xyze(1, 1) = 0.692770064;
    tri3_xyze(2, 1) = 0.368502647;
    tri3_xyze(0, 2) = 0.9287662802;
    tri3_xyze(1, 2) = 0.6797418;
    tri3_xyze(2, 2) = 0.4034825265;
    nids.clear();
    nids.push_back(431);
    nids.push_back(436);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.907204866;
    tri3_xyze(1, 0) = 0.649576843;
    tri3_xyze(2, 0) = 0.37970823;
    tri3_xyze(0, 1) = 0.905384481;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.438462406;
    tri3_xyze(0, 2) = 0.9287662802;
    tri3_xyze(1, 2) = 0.6797418;
    tri3_xyze(2, 2) = 0.4034825265;
    nids.clear();
    nids.push_back(437);
    nids.push_back(432);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.905384481;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.438462406;
    tri3_xyze(0, 1) = 0.950327694;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.427256823;
    tri3_xyze(0, 2) = 0.9287662802;
    tri3_xyze(1, 2) = 0.6797418;
    tri3_xyze(2, 2) = 0.4034825265;
    nids.clear();
    nids.push_back(432);
    nids.push_back(431);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649426;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.488708884;
    tri3_xyze(0, 1) = 1.02129769;
    tri3_xyze(1, 1) = 0.731349707;
    tri3_xyze(2, 1) = 0.474834234;
    tri3_xyze(0, 2) = 0.9858126925;
    tri3_xyze(1, 2) = 0.720628232;
    tri3_xyze(2, 2) = 0.451045521;
    nids.clear();
    nids.push_back(434);
    nids.push_back(406);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.02129769;
    tri3_xyze(1, 0) = 0.731349707;
    tri3_xyze(2, 0) = 0.474834234;
    tri3_xyze(0, 1) = 1.00597596;
    tri3_xyze(1, 1) = 0.731349707;
    tri3_xyze(2, 1) = 0.413382143;
    tri3_xyze(0, 2) = 0.9858126925;
    tri3_xyze(1, 2) = 0.720628232;
    tri3_xyze(2, 2) = 0.451045521;
    nids.clear();
    nids.push_back(406);
    nids.push_back(407);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.00597596;
    tri3_xyze(1, 0) = 0.731349707;
    tri3_xyze(2, 0) = 0.413382143;
    tri3_xyze(0, 1) = 0.950327694;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.427256823;
    tri3_xyze(0, 2) = 0.9858126925;
    tri3_xyze(1, 2) = 0.720628232;
    tri3_xyze(2, 2) = 0.451045521;
    nids.clear();
    nids.push_back(407);
    nids.push_back(431);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950327694;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.427256823;
    tri3_xyze(0, 1) = 0.965649426;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.488708884;
    tri3_xyze(0, 2) = 0.9858126925;
    tri3_xyze(1, 2) = 0.720628232;
    tri3_xyze(2, 2) = 0.451045521;
    nids.clear();
    nids.push_back(431);
    nids.push_back(434);
    nids.push_back(439);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.950327694;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.427256823;
    tri3_xyze(0, 1) = 1.00597596;
    tri3_xyze(1, 1) = 0.731349707;
    tri3_xyze(2, 1) = 0.413382143;
    tri3_xyze(0, 2) = 0.9756891572;
    tri3_xyze(1, 2) = 0.7126762273;
    tri3_xyze(2, 2) = 0.3944579215;
    nids.clear();
    nids.push_back(431);
    nids.push_back(407);
    nids.push_back(440);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.95214808;
    tri3_xyze(1, 0) = 0.692770064;
    tri3_xyze(2, 0) = 0.368502647;
    tri3_xyze(0, 1) = 0.950327694;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.427256823;
    tri3_xyze(0, 2) = 0.9756891572;
    tri3_xyze(1, 2) = 0.7126762273;
    tri3_xyze(2, 2) = 0.3944579215;
    nids.clear();
    nids.push_back(436);
    nids.push_back(431);
    nids.push_back(440);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.875762999;
    tri3_xyze(1, 0) = 0.623520255;
    tri3_xyze(2, 0) = 0.511120081;
    tri3_xyze(0, 1) = 0.920706213;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.499914467;
    tri3_xyze(0, 2) = 0.890573755;
    tri3_xyze(1, 2) = 0.6451168955;
    tri3_xyze(2, 2) = 0.4747912435;
    nids.clear();
    nids.push_back(442);
    nids.push_back(433);
    nids.push_back(443);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706213;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.499914467;
    tri3_xyze(0, 1) = 0.905384481;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.438462406;
    tri3_xyze(0, 2) = 0.890573755;
    tri3_xyze(1, 2) = 0.6451168955;
    tri3_xyze(2, 2) = 0.4747912435;
    nids.clear();
    nids.push_back(433);
    nids.push_back(432);
    nids.push_back(443);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649426;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.488708884;
    tri3_xyze(0, 1) = 0.920706213;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.499914467;
    tri3_xyze(0, 2) = 0.9508386853;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.5250377133;
    nids.clear();
    nids.push_back(434);
    nids.push_back(433);
    nids.push_back(448);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706213;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.499914467;
    tri3_xyze(0, 1) = 0.936027944;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.561366558;
    tri3_xyze(0, 2) = 0.9508386853;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.5250377133;
    nids.clear();
    nids.push_back(433);
    nids.push_back(446);
    nids.push_back(448);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.980971158;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.550160944;
    tri3_xyze(0, 1) = 0.965649426;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.488708884;
    tri3_xyze(0, 2) = 0.9508386853;
    tri3_xyze(1, 2) = 0.6883101465;
    tri3_xyze(2, 2) = 0.5250377133;
    nids.clear();
    nids.push_back(447);
    nids.push_back(434);
    nids.push_back(448);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.02129769;
    tri3_xyze(1, 0) = 0.731349707;
    tri3_xyze(2, 0) = 0.474834234;
    tri3_xyze(0, 1) = 0.965649426;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.488708884;
    tri3_xyze(0, 2) = 1.001134423;
    tri3_xyze(1, 2) = 0.720628232;
    tri3_xyze(2, 2) = 0.512497589;
    nids.clear();
    nids.push_back(406);
    nids.push_back(434);
    nids.push_back(449);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.965649426;
    tri3_xyze(1, 0) = 0.709906757;
    tri3_xyze(2, 0) = 0.488708884;
    tri3_xyze(0, 1) = 0.980971158;
    tri3_xyze(1, 1) = 0.709906757;
    tri3_xyze(2, 1) = 0.550160944;
    tri3_xyze(0, 2) = 1.001134423;
    tri3_xyze(1, 2) = 0.720628232;
    tri3_xyze(2, 2) = 0.512497589;
    nids.clear();
    nids.push_back(434);
    nids.push_back(447);
    nids.push_back(449);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.920706213;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.499914467;
    tri3_xyze(0, 1) = 0.875762999;
    tri3_xyze(1, 1) = 0.623520255;
    tri3_xyze(2, 1) = 0.511120081;
    tri3_xyze(0, 2) = 0.9058954717;
    tri3_xyze(1, 2) = 0.6451168955;
    tri3_xyze(2, 2) = 0.5362433045;
    nids.clear();
    nids.push_back(433);
    nids.push_back(442);
    nids.push_back(451);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.936027944;
    tri3_xyze(1, 0) = 0.666713536;
    tri3_xyze(2, 0) = 0.561366558;
    tri3_xyze(0, 1) = 0.920706213;
    tri3_xyze(1, 1) = 0.666713536;
    tri3_xyze(2, 1) = 0.499914467;
    tri3_xyze(0, 2) = 0.9058954717;
    tri3_xyze(1, 2) = 0.6451168955;
    tri3_xyze(2, 2) = 0.5362433045;
    nids.clear();
    nids.push_back(446);
    nids.push_back(433);
    nids.push_back(451);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.916666687;
  hex8_xyze(1, 0) = 0.666666687;
  hex8_xyze(2, 0) = 0.5;
  hex8_xyze(0, 1) = 0.916666687;
  hex8_xyze(1, 1) = 0.666666687;
  hex8_xyze(2, 1) = 0.416666657;
  hex8_xyze(0, 2) = 0.916666687;
  hex8_xyze(1, 2) = 0.75;
  hex8_xyze(2, 2) = 0.416666657;
  hex8_xyze(0, 3) = 0.916666687;
  hex8_xyze(1, 3) = 0.75;
  hex8_xyze(2, 3) = 0.5;
  hex8_xyze(0, 4) = 1;
  hex8_xyze(1, 4) = 0.666666687;
  hex8_xyze(2, 4) = 0.5;
  hex8_xyze(0, 5) = 1;
  hex8_xyze(1, 5) = 0.666666687;
  hex8_xyze(2, 5) = 0.416666657;
  hex8_xyze(0, 6) = 1;
  hex8_xyze(1, 6) = 0.75;
  hex8_xyze(2, 6) = 0.416666657;
  hex8_xyze(0, 7) = 1;
  hex8_xyze(1, 7) = 0.75;
  hex8_xyze(2, 7) = 0.5;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
