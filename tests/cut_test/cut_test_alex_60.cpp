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

void test_alex60()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3587888889;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3693444444;
    tri3_xyze(2, 2) = 0.3209;

    {
      int data[] = {597859447, 1072463033, -1409512822, 1071052389, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, -512341876, 1071096662, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(980);
    nids.push_back(984);
    nids.push_back(985);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.87015;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3693444444;
    tri3_xyze(2, 2) = 0.3209;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, -512341876, 1071096662, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(984);
    nids.push_back(787);
    nids.push_back(985);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.87015;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.31086875;

    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, 384829070, 1071140936, 178670640, 1070851398};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(787);
    nids.push_back(984);
    nids.push_back(989);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3008375;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.31086875;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, 384829070, 1071140936, 178670640, 1070851398};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(984);
    nids.push_back(988);
    nids.push_back(989);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3008375;
    tri3_xyze(0, 1) = 0.87015;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3008375;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.31086875;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, 384829070, 1071140936, 178670640, 1070851398};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(988);
    nids.push_back(825);
    nids.push_back(989);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.87015;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3008375;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3008375;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.29080625;

    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, 384829070, 1071140936, -783402035, 1070767249};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(825);
    nids.push_back(988);
    nids.push_back(991);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3008375;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.280775;
    tri3_xyze(0, 2) = 0.8801875;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.29080625;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -1264438372, 1070725175};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-103079216, 1072441982, 384829070, 1071140936, -783402035, 1070767249};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(988);
    nids.push_back(990);
    nids.push_back(991);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3693444444;
    tri3_xyze(2, 2) = 0.3209;

    {
      int data[] = {1999736773, 1072505133, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, -512341876, 1071096662, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(1120);
    nids.push_back(984);
    nids.push_back(1121);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3587888889;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3693444444;
    tri3_xyze(2, 2) = 0.3209;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, -1409512822, 1071052389, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, -512341876, 1071096662, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(984);
    nids.push_back(980);
    nids.push_back(1121);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3008375;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.31086875;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, 384829070, 1071140936, 178670640, 1070851398};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(988);
    nids.push_back(984);
    nids.push_back(1125);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3209;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3209;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.31086875;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1999736773, 1072505133, 384829070, 1071140936, 659706977, 1070893472};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, 384829070, 1071140936, 178670640, 1070851398};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(984);
    nids.push_back(1120);
    nids.push_back(1125);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3008375;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3008375;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.31086875;

    {
      int data[] = {1999736773, 1072505133, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, 384829070, 1071140936, 178670640, 1070851398};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(1124);
    nids.push_back(988);
    nids.push_back(1125);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.280775;
    tri3_xyze(0, 1) = 0.890225;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3008375;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.29080625;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -1264438372, 1070725175};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, 384829070, 1071140936, -783402035, 1070767249};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(990);
    nids.push_back(988);
    nids.push_back(1127);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.890225;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.3008375;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.3008375;
    tri3_xyze(0, 2) = 0.9002625;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.29080625;

    {
      int data[] = {597859447, 1072463033, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1999736773, 1072505133, 384829070, 1071140936, -302365698, 1070809323};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {1298798110, 1072484083, 384829070, 1071140936, -783402035, 1070767249};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(988);
    nids.push_back(1124);
    nids.push_back(1127);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.8955223881;
  hex8_xyze(1, 0) = 0.3851851852;
  hex8_xyze(2, 0) = 0.2941176471;
  hex8_xyze(0, 1) = 0.8880597015;
  hex8_xyze(1, 1) = 0.3851851852;
  hex8_xyze(2, 1) = 0.2941176471;
  hex8_xyze(0, 2) = 0.8880597015;
  hex8_xyze(1, 2) = 0.3777777778;
  hex8_xyze(2, 2) = 0.2941176471;
  hex8_xyze(0, 3) = 0.8955223881;
  hex8_xyze(1, 3) = 0.3777777778;
  hex8_xyze(2, 3) = 0.2941176471;
  hex8_xyze(0, 4) = 0.8955223881;
  hex8_xyze(1, 4) = 0.3851851852;
  hex8_xyze(2, 4) = 0.3235294118;
  hex8_xyze(0, 5) = 0.8880597015;
  hex8_xyze(1, 5) = 0.3851851852;
  hex8_xyze(2, 5) = 0.3235294118;
  hex8_xyze(0, 6) = 0.8880597015;
  hex8_xyze(1, 6) = 0.3777777778;
  hex8_xyze(2, 6) = 0.3235294118;
  hex8_xyze(0, 7) = 0.8955223881;
  hex8_xyze(1, 7) = 0.3777777778;
  hex8_xyze(2, 7) = 0.3235294118;

  {
    int data[] = {-1859015696, 1072474142, -1018066322, 1071163103, -757935390, 1070781138};
    std::memcpy(&hex8_xyze(0, 0), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247873, 1072458492, -1018066322, 1071163103, -757935390, 1070781138};
    std::memcpy(&hex8_xyze(0, 1), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247873, 1072458492, -668106024, 1071132034, -757935390, 1070781138};
    std::memcpy(&hex8_xyze(0, 2), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1859015696, 1072474142, -668106024, 1071132034, -757935390, 1070781138};
    std::memcpy(&hex8_xyze(0, 3), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1859015696, 1072474142, -1018066322, 1071163103, -1263225696, 1070904500};
    std::memcpy(&hex8_xyze(0, 4), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247873, 1072458492, -1018066322, 1071163103, -1263225696, 1070904500};
    std::memcpy(&hex8_xyze(0, 5), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247873, 1072458492, -668106024, 1071132034, -1263225696, 1070904500};
    std::memcpy(&hex8_xyze(0, 6), data, 3 * sizeof(double));
  }
  {
    int data[] = {-1859015696, 1072474142, -668106024, 1071132034, -1263225696, 1070904500};
    std::memcpy(&hex8_xyze(0, 7), data, 3 * sizeof(double));
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
