/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_alex_54.cpp

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

void test_alex54()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.03082333333;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.06165666667;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.8966666667;
    tri3_xyze(1, 2) = 0.04624;
    tri3_xyze(2, 2) = 0.31999;

    {
      int data[] = {1374389534, 1072504504, -495238368, 1067421734, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1259857073, 1072476542, -276252298, 1067953348, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(442);
    nids.push_back(446);
    nids.push_back(447);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.8833333333;
    tri3_xyze(1, 1) = 0.06165666667;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.8966666667;
    tri3_xyze(1, 2) = 0.04624;
    tri3_xyze(2, 2) = 0.31999;

    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1145324612, 1072448580, 1842598235, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1259857073, 1072476542, -276252298, 1067953348, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(446);
    nids.push_back(29);
    nids.push_back(447);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.8833333333;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.06165666667;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.8966666667;
    tri3_xyze(1, 2) = 0.07707333333;
    tri3_xyze(2, 2) = 0.31999;

    {
      int data[] = {1145324612, 1072448580, 1842598235, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1259857073, 1072476542, -166759264, 1068743443, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(29);
    nids.push_back(446);
    nids.push_back(469);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.09249;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.8966666667;
    tri3_xyze(1, 2) = 0.07707333333;
    tri3_xyze(2, 2) = 0.31999;

    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -1254817646, 1069002092, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1259857073, 1072476542, -166759264, 1068743443, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(446);
    nids.push_back(468);
    nids.push_back(469);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.28799;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.06165666667;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.91;
    tri3_xyze(1, 2) = 0.04624;
    tri3_xyze(2, 2) = 0.30399;

    {
      int data[] = {1374389534, 1072504504, 1842598235, 1068470646, -1679504012, 1070755437};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -276252298, 1067953348, 2031347732, 1070822546};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(531);
    nids.push_back(446);
    nids.push_back(532);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.03082333333;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.91;
    tri3_xyze(1, 2) = 0.04624;
    tri3_xyze(2, 2) = 0.30399;

    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -495238368, 1067421734, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -276252298, 1067953348, 2031347732, 1070822546};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(446);
    nids.push_back(442);
    nids.push_back(532);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.06165666667;
    tri3_xyze(2, 1) = 0.28799;
    tri3_xyze(0, 2) = 0.91;
    tri3_xyze(1, 2) = 0.07707333333;
    tri3_xyze(2, 2) = 0.30399;

    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, 1842598235, 1068470646, -1679504012, 1070755437};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -166759264, 1068743443, 2031347732, 1070822546};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(446);
    nids.push_back(531);
    nids.push_back(566);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.06165666667;
    tri3_xyze(2, 0) = 0.28799;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.09249;
    tri3_xyze(2, 1) = 0.28799;
    tri3_xyze(0, 2) = 0.91;
    tri3_xyze(1, 2) = 0.07707333333;
    tri3_xyze(2, 2) = 0.30399;

    {
      int data[] = {1374389534, 1072504504, 1842598235, 1068470646, -1679504012, 1070755437};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -1254817646, 1069002092, -1679504012, 1070755437};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -166759264, 1068743443, 2031347732, 1070822546};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(531);
    nids.push_back(565);
    nids.push_back(566);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.91;
    tri3_xyze(1, 0) = 0.09249;
    tri3_xyze(2, 0) = 0.31999;
    tri3_xyze(0, 1) = 0.91;
    tri3_xyze(1, 1) = 0.06165666667;
    tri3_xyze(2, 1) = 0.31999;
    tri3_xyze(0, 2) = 0.91;
    tri3_xyze(1, 2) = 0.07707333333;
    tri3_xyze(2, 2) = 0.30399;

    {
      int data[] = {1374389534, 1072504504, -1254817646, 1069002092, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, 1842598236, 1068470646, 1447232180, 1070889655};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1374389534, 1072504504, -166759264, 1068743443, 2031347732, 1070822546};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(468);
    nids.push_back(446);
    nids.push_back(566);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.9104477612;
  hex8_xyze(1, 0) = 0.06153846154;
  hex8_xyze(2, 0) = 0.2941176471;
  hex8_xyze(0, 1) = 0.9104477612;
  hex8_xyze(1, 1) = 0.06923076923;
  hex8_xyze(2, 1) = 0.2941176471;
  hex8_xyze(0, 2) = 0.9029850746;
  hex8_xyze(1, 2) = 0.06923076923;
  hex8_xyze(2, 2) = 0.2941176471;
  hex8_xyze(0, 3) = 0.9029850746;
  hex8_xyze(1, 3) = 0.06153846154;
  hex8_xyze(2, 3) = 0.2941176471;
  hex8_xyze(0, 4) = 0.9104477612;
  hex8_xyze(1, 4) = 0.06153846154;
  hex8_xyze(2, 4) = 0.3235294118;
  hex8_xyze(0, 5) = 0.9104477612;
  hex8_xyze(1, 5) = 0.06923076923;
  hex8_xyze(2, 5) = 0.3235294118;
  hex8_xyze(0, 6) = 0.9029850746;
  hex8_xyze(1, 6) = 0.06923076923;
  hex8_xyze(2, 6) = 0.3235294118;
  hex8_xyze(0, 7) = 0.9029850746;
  hex8_xyze(1, 7) = 0.06153846154;
  hex8_xyze(2, 7) = 0.3235294118;

  {
    int data[] = {1474391756, 1072505443, 528611326, 1068466680, -757935383, 1070781138};
    std::memcpy(&hex8_xyze(0, 0), data, 3 * sizeof(double));
  }
  {
    int data[] = {1474391755, 1072505443, -1850139775, 1068611867, -757935383, 1070781138};
    std::memcpy(&hex8_xyze(0, 1), data, 3 * sizeof(double));
  }
  {
    int data[] = {-192311971, 1072489792, -1850139775, 1068611867, -757935383, 1070781138};
    std::memcpy(&hex8_xyze(0, 2), data, 3 * sizeof(double));
  }
  {
    int data[] = {-192311970, 1072489792, 528611326, 1068466680, -757935383, 1070781138};
    std::memcpy(&hex8_xyze(0, 3), data, 3 * sizeof(double));
  }
  {
    int data[] = {1474391756, 1072505443, 528611328, 1068466680, -1263225680, 1070904500};
    std::memcpy(&hex8_xyze(0, 4), data, 3 * sizeof(double));
  }
  {
    int data[] = {1474391755, 1072505443, -1850139774, 1068611867, -1263225680, 1070904500};
    std::memcpy(&hex8_xyze(0, 5), data, 3 * sizeof(double));
  }
  {
    int data[] = {-192311971, 1072489792, -1850139774, 1068611867, -1263225680, 1070904500};
    std::memcpy(&hex8_xyze(0, 6), data, 3 * sizeof(double));
  }
  {
    int data[] = {-192311970, 1072489792, 528611328, 1068466680, -1263225680, 1070904500};
    std::memcpy(&hex8_xyze(0, 7), data, 3 * sizeof(double));
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
