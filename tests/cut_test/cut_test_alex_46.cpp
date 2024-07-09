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

void test_alex46()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.856767;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.87015;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;
    nids.clear();
    nids.push_back(433);
    nids.push_back(435);
    nids.push_back(438);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.87015;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;
    nids.clear();
    nids.push_back(435);
    nids.push_back(399);
    nids.push_back(438);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.856767;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.87015;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;
    nids.clear();
    nids.push_back(399);
    nids.push_back(397);
    nids.push_back(438);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;
    nids.clear();
    nids.push_back(399);
    nids.push_back(435);
    nids.push_back(526);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = 0.032;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = 0.032;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;
    nids.clear();
    nids.push_back(435);
    nids.push_back(524);
    nids.push_back(526);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.3799;
    tri3_xyze(2, 0) = -0.0001;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.3799;
    tri3_xyze(2, 1) = -0.0001;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.3799;
    tri3_xyze(2, 2) = 0.01595;
    nids.clear();
    nids.push_back(506);
    nids.push_back(399);
    nids.push_back(526);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.880597;
  hex8_xyze(1, 0) = 0.385185;
  hex8_xyze(2, 0) = 0;
  hex8_xyze(0, 1) = 0.880597;
  hex8_xyze(1, 1) = 0.377778;
  hex8_xyze(2, 1) = 0;
  hex8_xyze(0, 2) = 0.88806;
  hex8_xyze(1, 2) = 0.377778;
  hex8_xyze(2, 2) = 0;
  hex8_xyze(0, 3) = 0.88806;
  hex8_xyze(1, 3) = 0.385185;
  hex8_xyze(2, 3) = 0;
  hex8_xyze(0, 4) = 0.880597;
  hex8_xyze(1, 4) = 0.385185;
  hex8_xyze(2, 4) = 0.0294118;
  hex8_xyze(0, 5) = 0.880597;
  hex8_xyze(1, 5) = 0.377778;
  hex8_xyze(2, 5) = 0.0294118;
  hex8_xyze(0, 6) = 0.88806;
  hex8_xyze(1, 6) = 0.377778;
  hex8_xyze(2, 6) = 0.0294118;
  hex8_xyze(0, 7) = 0.88806;
  hex8_xyze(1, 7) = 0.385185;
  hex8_xyze(2, 7) = 0.0294118;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}
