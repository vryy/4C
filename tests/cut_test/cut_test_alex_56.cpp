/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

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

void test_alex56()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.8567666667;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.8835333333;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.87015;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {-306946996, 1072392865, 384829070, 1071140936, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(433);
    nids.push_back(435);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.8835333333;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.8835333333;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.87015;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(435);
    nids.push_back(399);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.8835333333;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.8567666667;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.87015;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-306946996, 1072392865, 384829070, 1071140936, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-804017878, 1072420932, 384829070, 1071140936, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(399);
    nids.push_back(397);
    nids.push_back(438);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.8835333333;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.8835333333;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.8969166667;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1798159642, 1072477066, 384829070, 1071140936, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(399);
    nids.push_back(435);
    nids.push_back(526);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.8835333333;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.8969166667;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1999736773, 1072505133, 384829070, 1071140936, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1798159642, 1072477066, 384829070, 1071140936, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(435);
    nids.push_back(524);
    nids.push_back(526);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.8835333333;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.8969166667;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {1999736773, 1072505133, 384829070, 1071140936, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1301088760, 1072448999, 384829070, 1071140936, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1798159642, 1072477066, 384829070, 1071140936, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(506);
    nids.push_back(399);
    nids.push_back(526);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.8805970149;
  hex8_xyze(1, 0) = 0.3851851852;
  hex8_xyze(2, 0) = 0;
  hex8_xyze(0, 1) = 0.8805970149;
  hex8_xyze(1, 1) = 0.3777777778;
  hex8_xyze(2, 1) = 0;
  hex8_xyze(0, 2) = 0.8880597015;
  hex8_xyze(1, 2) = 0.3777777778;
  hex8_xyze(2, 2) = 0;
  hex8_xyze(0, 3) = 0.8880597015;
  hex8_xyze(1, 3) = 0.3851851852;
  hex8_xyze(2, 3) = 0;
  hex8_xyze(0, 4) = 0.8805970149;
  hex8_xyze(1, 4) = 0.3851851852;
  hex8_xyze(2, 4) = 0.02941176471;
  hex8_xyze(0, 5) = 0.8805970149;
  hex8_xyze(1, 5) = 0.3777777778;
  hex8_xyze(2, 5) = 0.02941176471;
  hex8_xyze(0, 6) = 0.8880597015;
  hex8_xyze(1, 6) = 0.3777777778;
  hex8_xyze(2, 6) = 0.02941176471;
  hex8_xyze(0, 7) = 0.8880597015;
  hex8_xyze(1, 7) = 0.3851851852;
  hex8_xyze(2, 7) = 0.02941176471;

  {
    int data[] = {-897455855, 1072442841, -1018066322, 1071163103, 0, 0};
    std::memcpy(&hex8_xyze(0, 0), data, 3 * sizeof(double));
  }
  {
    int data[] = {-897455855, 1072442841, -668106024, 1071132034, 0, 0};
    std::memcpy(&hex8_xyze(0, 1), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247872, 1072458492, -668106024, 1071132034, 0, 0};
    std::memcpy(&hex8_xyze(0, 2), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247872, 1072458492, -1018066322, 1071163103, 0, 0};
    std::memcpy(&hex8_xyze(0, 3), data, 3 * sizeof(double));
  }
  {
    int data[] = {-897455855, 1072442841, -1018066322, 1071163103, 505290247, 1067327006};
    std::memcpy(&hex8_xyze(0, 4), data, 3 * sizeof(double));
  }
  {
    int data[] = {-897455855, 1072442841, -668106024, 1071132034, 505290247, 1067327006};
    std::memcpy(&hex8_xyze(0, 5), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247872, 1072458492, -668106024, 1071132034, 505290247, 1067327006};
    std::memcpy(&hex8_xyze(0, 6), data, 3 * sizeof(double));
  }
  {
    int data[] = {769247872, 1072458492, -1018066322, 1071163103, 505290247, 1067327006};
    std::memcpy(&hex8_xyze(0, 7), data, 3 * sizeof(double));
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
