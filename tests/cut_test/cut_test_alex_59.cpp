/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_alex59()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.03156666667;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.06323333333;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -1846263273, 1067460993, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(85);
    nids.push_back(90);
    nids.push_back(100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.06323333333;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.06323333333;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(90);
    nids.push_back(99);
    nids.push_back(100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.06323333333;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.03156666667;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.0474;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1846263272, 1067460993, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 329853492, 1067992272, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(99);
    nids.push_back(94);
    nids.push_back(100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.06323333333;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.0949;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.07906666667;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -865865406, 1069042525, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1953923787, 1068776886, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(90);
    nids.push_back(132);
    nids.push_back(135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.0949;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.06323333333;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.07906666667;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, -865865406, 1069042525, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1953923787, 1068776886, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(134);
    nids.push_back(99);
    nids.push_back(135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.83;
    tri3_xyze(1, 0) = 0.06323333333;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.83;
    tri3_xyze(1, 1) = 0.06323333333;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.83;
    tri3_xyze(1, 2) = 0.07906666667;
    tri3_xyze(2, 2) = 0.01595;

    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -350469331, -1088801054};
      std::memcpy(&tri3_xyze(0, 0), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, 1252985128, 1068511247, -755914248, 1067475533};
      std::memcpy(&tri3_xyze(0, 1), data, 3 * sizeof(double));
    }
    {
      int data[] = {687194767, 1072336732, -1953923787, 1068776886, 1635523542, 1066423602};
      std::memcpy(&tri3_xyze(0, 2), data, 3 * sizeof(double));
    }

    nids.clear();
    nids.push_back(99);
    nids.push_back(90);
    nids.push_back(135);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.8358208955;
  hex8_xyze(1, 0) = 0.05925925926;
  hex8_xyze(2, 0) = 0;
  hex8_xyze(0, 1) = 0.8358208955;
  hex8_xyze(1, 1) = 0.06666666667;
  hex8_xyze(2, 1) = 0;
  hex8_xyze(0, 2) = 0.828358209;
  hex8_xyze(1, 2) = 0.06666666667;
  hex8_xyze(2, 2) = 0;
  hex8_xyze(0, 3) = 0.828358209;
  hex8_xyze(1, 3) = 0.05925925926;
  hex8_xyze(2, 3) = 0;
  hex8_xyze(0, 4) = 0.8358208955;
  hex8_xyze(1, 4) = 0.05925925926;
  hex8_xyze(2, 4) = 0.02941176471;
  hex8_xyze(0, 5) = 0.8358208955;
  hex8_xyze(1, 5) = 0.06666666667;
  hex8_xyze(2, 5) = 0.02941176471;
  hex8_xyze(0, 6) = 0.828358209;
  hex8_xyze(1, 6) = 0.06666666667;
  hex8_xyze(2, 6) = 0.02941176471;
  hex8_xyze(0, 7) = 0.828358209;
  hex8_xyze(1, 7) = 0.05925925926;
  hex8_xyze(2, 7) = 0.02941176471;

  {
    int data[] = {1987223673, 1072348939, -922622608, 1068390202, 0, 0};
    std::memcpy(&hex8_xyze(0, 0), data, 3 * sizeof(double));
  }
  {
    int data[] = {1987223673, 1072348939, 286331152, 1068568849, 0, 0};
    std::memcpy(&hex8_xyze(0, 1), data, 3 * sizeof(double));
  }
  {
    int data[] = {320519946, 1072333289, 286331152, 1068568849, 0, 0};
    std::memcpy(&hex8_xyze(0, 2), data, 3 * sizeof(double));
  }
  {
    int data[] = {320519946, 1072333289, -922622608, 1068390202, 0, 0};
    std::memcpy(&hex8_xyze(0, 3), data, 3 * sizeof(double));
  }
  {
    int data[] = {1987223673, 1072348939, -922622612, 1068390202, 505290245, 1067327006};
    std::memcpy(&hex8_xyze(0, 4), data, 3 * sizeof(double));
  }
  {
    int data[] = {1987223673, 1072348939, 286331150, 1068568849, 505290245, 1067327006};
    std::memcpy(&hex8_xyze(0, 5), data, 3 * sizeof(double));
  }
  {
    int data[] = {320519947, 1072333289, 286331150, 1068568849, 505290245, 1067327006};
    std::memcpy(&hex8_xyze(0, 6), data, 3 * sizeof(double));
  }
  {
    int data[] = {320519947, 1072333289, -922622612, 1068390202, 505290245, 1067327006};
    std::memcpy(&hex8_xyze(0, 7), data, 3 * sizeof(double));
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}
