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

#include "../../src/drt_fem_general/drt_utils_local_connectivity_matrices.H"

void test_shadan6()
{
  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.2354503677e-01;
    tri3_xyze(1, 0) = 4.7072755304e-01;
    tri3_xyze(2, 0) = -3.0000000000e-02;
    tri3_xyze(0, 1) = 1.1182310292e-01;
    tri3_xyze(1, 1) = 4.4311240744e-01;
    tri3_xyze(2, 1) = -3.0000000000e-02;
    tri3_xyze(0, 2) = 1.1768406985e-01;
    tri3_xyze(1, 2) = 4.5691998024e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-997849845, 1069523109, 1969820200, 1071521894, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {255813208, 1069424779, 1291933380, 1071463981, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(18);
    nids.push_back(15);
    nids.push_back(19);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.1182310292e-01;
    tri3_xyze(1, 0) = 4.4311240744e-01;
    tri3_xyze(2, 0) = -3.0000000000e-02;
    tri3_xyze(0, 1) = 1.1182310292e-01;
    tri3_xyze(1, 1) = 4.4311240744e-01;
    tri3_xyze(2, 1) = 0.0000000000e+00;
    tri3_xyze(0, 2) = 1.1768406985e-01;
    tri3_xyze(1, 2) = 4.5691998024e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, 0, 0};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {255813208, 1069424779, 1291933380, 1071463981, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(15);
    nids.push_back(14);
    nids.push_back(19);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.0010116906e-01;
    tri3_xyze(1, 0) = 4.1549726183e-01;
    tri3_xyze(2, 0) = -3.0000000000e-02;
    tri3_xyze(0, 1) = 8.8379235209e-02;
    tri3_xyze(1, 1) = 3.8788211623e-01;
    tri3_xyze(2, 1) = -3.0000000000e-02;
    tri3_xyze(0, 2) = 9.4240202136e-02;
    tri3_xyze(1, 2) = 4.0168968903e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-2065806119, 1068933125, -2097500724, 1071174415, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {975498124, 1069031456, -1419613903, 1071232328, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(9);
    nids.push_back(2);
    nids.push_back(10);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.0010116906e-01;
    tri3_xyze(1, 0) = 4.1549726183e-01;
    tri3_xyze(2, 0) = 0.0000000000e+00;
    tri3_xyze(0, 1) = 1.0010116906e-01;
    tri3_xyze(1, 1) = 4.1549726183e-01;
    tri3_xyze(2, 1) = -3.0000000000e-02;
    tri3_xyze(0, 2) = 9.4240202136e-02;
    tri3_xyze(1, 2) = 4.0168968903e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, 0, 0};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {975498124, 1069031456, -1419613903, 1071232328, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(8);
    nids.push_back(9);
    nids.push_back(10);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.1182310292e-01;
    tri3_xyze(1, 0) = 4.4311240744e-01;
    tri3_xyze(2, 0) = -3.0000000000e-02;
    tri3_xyze(0, 1) = 1.0010116906e-01;
    tri3_xyze(1, 1) = 4.1549726183e-01;
    tri3_xyze(2, 1) = -3.0000000000e-02;
    tri3_xyze(0, 2) = 1.0596213599e-01;
    tri3_xyze(1, 2) = 4.2930483464e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1531827982, 1069228117, -63840262, 1071348154, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(15);
    nids.push_back(9);
    nids.push_back(16);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.0010116906e-01;
    tri3_xyze(1, 0) = 4.1549726183e-01;
    tri3_xyze(2, 0) = -3.0000000000e-02;
    tri3_xyze(0, 1) = 1.0010116906e-01;
    tri3_xyze(1, 1) = 4.1549726183e-01;
    tri3_xyze(2, 1) = 0.0000000000e+00;
    tri3_xyze(0, 2) = 1.0596213599e-01;
    tri3_xyze(1, 2) = 4.2930483464e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, 0, 0};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1531827982, 1069228117, -63840262, 1071348154, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(9);
    nids.push_back(8);
    nids.push_back(16);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.1182310292e-01;
    tri3_xyze(1, 0) = 4.4311240744e-01;
    tri3_xyze(2, 0) = 0.0000000000e+00;
    tri3_xyze(0, 1) = 1.1182310292e-01;
    tri3_xyze(1, 1) = 4.4311240744e-01;
    tri3_xyze(2, 1) = -3.0000000000e-02;
    tri3_xyze(0, 2) = 1.0596213599e-01;
    tri3_xyze(1, 2) = 4.2930483464e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, 0, 0};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1531827982, 1069228117, -63840262, 1071348154, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(14);
    nids.push_back(15);
    nids.push_back(16);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8379235209e-02;
    tri3_xyze(1, 0) = 3.8788211623e-01;
    tri3_xyze(2, 0) = -3.0000000000e-02;
    tri3_xyze(0, 1) = 8.8379235209e-02;
    tri3_xyze(1, 1) = 3.8788211623e-01;
    tri3_xyze(2, 1) = 0.0000000000e+00;
    tri3_xyze(0, 2) = 9.4240202136e-02;
    tri3_xyze(1, 2) = 4.0168968903e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-2065806119, 1068933125, -2097500724, 1071174415, -343597384, -1080117167};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-2065806119, 1068933125, -2097500724, 1071174415, 0, 0};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {975498124, 1069031456, -1419613903, 1071232328, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(2);
    nids.push_back(1);
    nids.push_back(10);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.1182310292e-01;
    tri3_xyze(1, 0) = 4.4311240744e-01;
    tri3_xyze(2, 0) = 0.0000000000e+00;
    tri3_xyze(0, 1) = 1.2354503677e-01;
    tri3_xyze(1, 1) = 4.7072755304e-01;
    tri3_xyze(2, 1) = 0.0000000000e+00;
    tri3_xyze(0, 2) = 1.1768406985e-01;
    tri3_xyze(1, 2) = 4.5691998024e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, 0, 0};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-997849845, 1069523109, 1969820200, 1071521894, 0, 0};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {255813208, 1069424779, 1291933380, 1071463981, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(14);
    nids.push_back(17);
    nids.push_back(19);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 8.8379235209e-02;
    tri3_xyze(1, 0) = 3.8788211623e-01;
    tri3_xyze(2, 0) = 0.0000000000e+00;
    tri3_xyze(0, 1) = 1.0010116906e-01;
    tri3_xyze(1, 1) = 4.1549726183e-01;
    tri3_xyze(2, 1) = 0.0000000000e+00;
    tri3_xyze(0, 2) = 9.4240202136e-02;
    tri3_xyze(1, 2) = 4.0168968903e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-2065806119, 1068933125, -2097500724, 1071174415, 0, 0};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, 0, 0};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {975498124, 1069031456, -1419613903, 1071232328, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(1);
    nids.push_back(8);
    nids.push_back(10);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 1.0010116906e-01;
    tri3_xyze(1, 0) = 4.1549726183e-01;
    tri3_xyze(2, 0) = 0.0000000000e+00;
    tri3_xyze(0, 1) = 1.1182310292e-01;
    tri3_xyze(1, 1) = 4.4311240744e-01;
    tri3_xyze(2, 1) = 0.0000000000e+00;
    tri3_xyze(0, 2) = 1.0596213599e-01;
    tri3_xyze(1, 2) = 4.2930483464e-01;
    tri3_xyze(2, 2) = -1.5000000000e-02;

    {
      int data[] = {-278164929, 1069129786, -741727082, 1071290241, 0, 0};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {1509476261, 1069326448, 614046558, 1071406068, 0, 0};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {-1531827982, 1069228117, -63840262, 1071348154, -343597384, -1081165743};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(8);
    nids.push_back(14);
    nids.push_back(16);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, DRT::Element::tri3);
  }
  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 5.7142857143e-02;
  hex8_xyze(1, 0) = 4.0000000000e-01;
  hex8_xyze(2, 0) = 0.0000000000e+00;
  hex8_xyze(0, 1) = 5.7142857143e-02;
  hex8_xyze(1, 1) = 4.0000000000e-01;
  hex8_xyze(2, 1) = -3.0000000000e-02;
  hex8_xyze(0, 2) = 5.7142857143e-02;
  hex8_xyze(1, 2) = 4.5714285714e-01;
  hex8_xyze(2, 2) = -3.0000000000e-02;
  hex8_xyze(0, 3) = 5.7142857143e-02;
  hex8_xyze(1, 3) = 4.5714285714e-01;
  hex8_xyze(2, 3) = 0.0000000000e+00;
  hex8_xyze(0, 4) = 1.1428571429e-01;
  hex8_xyze(1, 4) = 4.0000000000e-01;
  hex8_xyze(2, 4) = 0.0000000000e+00;
  hex8_xyze(0, 5) = 1.1428571429e-01;
  hex8_xyze(1, 5) = 4.0000000000e-01;
  hex8_xyze(2, 5) = -3.0000000000e-02;
  hex8_xyze(0, 6) = 1.1428571429e-01;
  hex8_xyze(1, 6) = 4.5714285714e-01;
  hex8_xyze(2, 6) = -3.0000000000e-02;
  hex8_xyze(0, 7) = 1.1428571429e-01;
  hex8_xyze(1, 7) = 4.5714285714e-01;
  hex8_xyze(2, 7) = 0.0000000000e+00;

  {
    int data[] = {490853408, 1068319188, -1717986918, 1071225241, 0, 0};
    std::memcpy(&hex8_xyze(0, 0), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1068319188, -1717986918, 1071225241, -343597384, -1080117167};
    std::memcpy(&hex8_xyze(0, 1), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1068319188, 490853405, 1071464916, -343597385, -1080117167};
    std::memcpy(&hex8_xyze(0, 2), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1068319188, 490853405, 1071464916, 0, 0};
    std::memcpy(&hex8_xyze(0, 3), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1069367764, -1717986918, 1071225241, 0, 0};
    std::memcpy(&hex8_xyze(0, 4), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1069367764, -1717986918, 1071225241, -343597384, -1080117167};
    std::memcpy(&hex8_xyze(0, 5), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1069367764, 490853406, 1071464916, -343597384, -1080117167};
    std::memcpy(&hex8_xyze(0, 6), data, 3 * sizeof(double));
  }
  {
    int data[] = {490853408, 1069367764, 490853406, 1071464916, 0, 0};
    std::memcpy(&hex8_xyze(0, 7), data, 3 * sizeof(double));
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.AddElement(1, nids, hex8_xyze, DRT::Element::hex8);


  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_DirectDivergence);
}
