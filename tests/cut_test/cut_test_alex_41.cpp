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

void test_alex41()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.1283;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.4201;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(811);
    nids.push_back(820);
    nids.push_back(1203);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.4201;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(820);
    nids.push_back(1202);
    nids.push_back(1203);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.1283;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.4201;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(1202);
    nids.push_back(1200);
    nids.push_back(1203);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.883533;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0641;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.4201;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(820);
    nids.push_back(829);
    nids.push_back(1205);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0641;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.4201;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1204);
    nids.push_back(1202);
    nids.push_back(1205);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.883533;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.896917;
    tri3_xyze(1, 2) = 0.4201;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1202);
    nids.push_back(820);
    nids.push_back(1205);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.1283;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.435933;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(1200);
    nids.push_back(1202);
    nids.push_back(1298);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.451767;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.435933;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(1202);
    nids.push_back(1297);
    nids.push_back(1298);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.451767;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.451767;
    tri3_xyze(2, 1) = 0.1283;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.435933;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(1297);
    nids.push_back(1293);
    nids.push_back(1298);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.451767;
    tri3_xyze(2, 0) = 0.1283;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.451767;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.4676;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(1293);
    nids.push_back(1297);
    nids.push_back(1300);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.451767;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.483433;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.4676;
    tri3_xyze(2, 2) = 0.11225;
    nids.clear();
    nids.push_back(1297);
    nids.push_back(1299);
    nids.push_back(1300);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.4201;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0641;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.435933;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1202);
    nids.push_back(1204);
    nids.push_back(1302);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.451767;
    tri3_xyze(2, 0) = 0.0641;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.451767;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.435933;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1301);
    nids.push_back(1297);
    nids.push_back(1302);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.451767;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.4201;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.435933;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1297);
    nids.push_back(1202);
    nids.push_back(1302);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.451767;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.451767;
    tri3_xyze(2, 1) = 0.0641;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.4676;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1297);
    nids.push_back(1301);
    nids.push_back(1304);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 0) = 0.9103;
    tri3_xyze(1, 0) = 0.483433;
    tri3_xyze(2, 0) = 0.0962;
    tri3_xyze(0, 1) = 0.9103;
    tri3_xyze(1, 1) = 0.451767;
    tri3_xyze(2, 1) = 0.0962;
    tri3_xyze(0, 2) = 0.9103;
    tri3_xyze(1, 2) = 0.4676;
    tri3_xyze(2, 2) = 0.08015;
    nids.clear();
    nids.push_back(1299);
    nids.push_back(1297);
    nids.push_back(1304);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 0.923077;
  hex8_xyze(1, 0) = 0.457143;
  hex8_xyze(2, 0) = 0.111111;
  hex8_xyze(0, 1) = 0.923077;
  hex8_xyze(1, 1) = 0.419048;
  hex8_xyze(2, 1) = 0.111111;
  hex8_xyze(0, 2) = 0.884615;
  hex8_xyze(1, 2) = 0.419048;
  hex8_xyze(2, 2) = 0.111111;
  hex8_xyze(0, 3) = 0.884615;
  hex8_xyze(1, 3) = 0.457143;
  hex8_xyze(2, 3) = 0.111111;
  hex8_xyze(0, 4) = 0.923077;
  hex8_xyze(1, 4) = 0.457143;
  hex8_xyze(2, 4) = 0.0833333;
  hex8_xyze(0, 5) = 0.923077;
  hex8_xyze(1, 5) = 0.419048;
  hex8_xyze(2, 5) = 0.0833333;
  hex8_xyze(0, 6) = 0.884615;
  hex8_xyze(1, 6) = 0.419048;
  hex8_xyze(2, 6) = 0.0833333;
  hex8_xyze(0, 7) = 0.884615;
  hex8_xyze(1, 7) = 0.457143;
  hex8_xyze(2, 7) = 0.0833333;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);

  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}
