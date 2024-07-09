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
#include <list>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_hex8quad4selfcut20()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(14);
    nids.push_back(13);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut21()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut22()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.0;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.0;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.0;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut23()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.0;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.0;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut24()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut30()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(14);
    nids.push_back(13);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(13);
    nids.push_back(12);
    nids.push_back(33);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut31()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut32()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut33()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut34()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut35()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = -0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = -0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut36()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut37()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut38()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.4;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.4;

    nids.clear();
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(27);
    nids.push_back(28);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut39()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.0;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.0;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.0;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.0;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(27);
    nids.push_back(28);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut41()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut42()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);

  std::vector<double> dirdivVol;

  Core::Geo::Cut::Mesh mesh = intersection.normal_mesh();
  const std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell>>& other_cells = mesh.volume_cells();
  for (std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    Core::Geo::Cut::VolumeCell* vc = &**i;
    dirdivVol.push_back(vc->volume());
  }

  for (uint i = 0; i < dirdivVol.size(); ++i)
    std::cout << "volume is: " << dirdivVol[i] << std::endl;
}

void test_hex8quad4selfcut43()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -2.3;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -2.3;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut51()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(27);
    nids.push_back(28);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(27);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(28);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut52()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(27);
    nids.push_back(28);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(27);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(28);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut53()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -2.3;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -2.3;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(27);
    nids.push_back(28);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(27);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(28);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut61()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut62()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut63()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut64()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut65()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut66()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut71()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.7;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.7;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.7;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.7;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);

  std::vector<double> tessVol, momFitVol, dirDivVol;

  Core::Geo::Cut::Mesh mesh = intersection.normal_mesh();
  const std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell>>& other_cells = mesh.volume_cells();
  for (std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    Core::Geo::Cut::VolumeCell* vc = &**i;
    tessVol.push_back(vc->volume());
  }

  for (std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    Core::Geo::Cut::VolumeCell* vc = &**i;
    vc->moment_fit_gauss_weights(
        vc->parent_element(), mesh, true, Core::Geo::Cut::BCellGaussPts_Tessellation);
    momFitVol.push_back(vc->volume());
  }

  for (std::list<Teuchos::RCP<Core::Geo::Cut::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    Core::Geo::Cut::VolumeCell* vc = &**i;
    vc->direct_divergence_gauss_rule(
        vc->parent_element(), mesh, true, Core::Geo::Cut::BCellGaussPts_Tessellation);
    dirDivVol.push_back(vc->volume());
  }

  std::cout << "the volumes predicted by\ntessellation \tMomentFitting \tDirectDivergence\n";
  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    std::cout << tessVol[i] << "\t\t" << momFitVol[i] << "\t\t" << dirDivVol[i] << "\n";
    if (fabs(tessVol[i] - momFitVol[i]) > 1e-9 || fabs(dirDivVol[i] - momFitVol[i]) > 1e-9)
      std::cout << "WARNING: volume predicted by either one of the method is wrong\n";
  }
}

void test_hex8quad4selfcut72()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.7;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.7;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.7;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.7;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut81()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.45;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.45;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.45;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.45;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut82()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.45;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.45;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.45;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.45;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.95;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.95;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.95;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.95;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut83()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }


  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut84()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }


  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = -0.2;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(41);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(44);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = -0.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(45);
    nids.push_back(46);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut85()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }


  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut86()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.5;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.5;

    nids.clear();
    nids.push_back(22);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(25);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(21);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(32);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(33);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }


  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.7;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.5;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.5;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.7;

    nids.clear();
    nids.push_back(42);
    nids.push_back(41);
    nids.push_back(44);
    nids.push_back(43);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.7;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.7;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(45);
    nids.push_back(42);
    nids.push_back(43);
    nids.push_back(46);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut91()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -0.2;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -0.2;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(40);
    nids.push_back(39);
    nids.push_back(30);
    nids.push_back(41);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

void test_hex8quad4selfcut92()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.2;
    quad4_xyze(1, 0) = -0.2;
    quad4_xyze(2, 0) = -3.0;

    quad4_xyze(0, 1) = 0.2;
    quad4_xyze(1, 1) = -0.2;
    quad4_xyze(2, 1) = 1.2;

    quad4_xyze(0, 2) = 0.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 1.2;

    quad4_xyze(0, 3) = 0.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = -3.0;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.2;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = 0.2;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = 0.2;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.2;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(21);
    nids.push_back(22);
    nids.push_back(23);
    nids.push_back(24);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 0.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.8;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.8;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 0.2;

    nids.clear();
    nids.push_back(22);
    nids.push_back(25);
    nids.push_back(26);
    nids.push_back(23);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.8;
    quad4_xyze(2, 0) = 0.8;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.2;
    quad4_xyze(2, 1) = 0.8;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.2;
    quad4_xyze(2, 2) = 0.8;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.8;
    quad4_xyze(2, 3) = 0.8;

    nids.clear();
    nids.push_back(25);
    nids.push_back(21);
    nids.push_back(24);
    nids.push_back(26);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.3;
    quad4_xyze(2, 0) = 0.9;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 1.2;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 1.2;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.3;
    quad4_xyze(2, 3) = 0.9;

    nids.clear();
    nids.push_back(31);
    nids.push_back(32);
    nids.push_back(33);
    nids.push_back(34);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.4;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.3;
    quad4_xyze(2, 1) = 0.9;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.3;
    quad4_xyze(2, 2) = 0.9;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.4;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(35);
    nids.push_back(31);
    nids.push_back(34);
    nids.push_back(36);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.0;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.4;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.4;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.0;

    nids.clear();
    nids.push_back(37);
    nids.push_back(35);
    nids.push_back(36);
    nids.push_back(38);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 0.6;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.0;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.0;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 0.6;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(39);
    nids.push_back(37);
    nids.push_back(38);
    nids.push_back(30);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = -0.2;
    quad4_xyze(1, 0) = 1.2;
    quad4_xyze(2, 0) = 1.2;

    quad4_xyze(0, 1) = -0.2;
    quad4_xyze(1, 1) = 0.6;
    quad4_xyze(2, 1) = 1.1;

    quad4_xyze(0, 2) = 1.2;
    quad4_xyze(1, 2) = 0.6;
    quad4_xyze(2, 2) = 1.1;

    quad4_xyze(0, 3) = 1.2;
    quad4_xyze(1, 3) = 1.2;
    quad4_xyze(2, 3) = 1.2;

    nids.clear();
    nids.push_back(40);
    nids.push_back(39);
    nids.push_back(30);
    nids.push_back(41);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);
  {
    hex8_xyze(0, 0) = 1.0;
    hex8_xyze(1, 0) = 1.0;
    hex8_xyze(2, 0) = 1.0;
    hex8_xyze(0, 1) = 1.0;
    hex8_xyze(1, 1) = 0.0;
    hex8_xyze(2, 1) = 1.0;
    hex8_xyze(0, 2) = 0.0;
    hex8_xyze(1, 2) = 0.0;
    hex8_xyze(2, 2) = 1.0;
    hex8_xyze(0, 3) = 0.0;
    hex8_xyze(1, 3) = 1.0;
    hex8_xyze(2, 3) = 1.0;
    hex8_xyze(0, 4) = 1.0;
    hex8_xyze(1, 4) = 1.0;
    hex8_xyze(2, 4) = 0.0;
    hex8_xyze(0, 5) = 1.0;
    hex8_xyze(1, 5) = 0.0;
    hex8_xyze(2, 5) = 0.0;
    hex8_xyze(0, 6) = 0.0;
    hex8_xyze(1, 6) = 0.0;
    hex8_xyze(2, 6) = 0.0;
    hex8_xyze(0, 7) = 0.0;
    hex8_xyze(1, 7) = 1.0;
    hex8_xyze(2, 7) = 0.0;
  }

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);
  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}

/*------------------------------------------------------------------------------------*
 * hex8 background element is cut by two quad cut surfaces. An edge of one cut surface
 * coincides with an edge of another cut surface
 *------------------------------------------------------------------------------------*/
void test_hex8quad4alignedEdges()
{
  Core::Geo::Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.5;  // 0.8;
    quad4_xyze(1, 0) = -0.1;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = 0.5;  // 0.8;
    quad4_xyze(1, 1) = -0.1;
    quad4_xyze(2, 1) = -0.1;

    quad4_xyze(0, 2) = 0.5;
    quad4_xyze(1, 2) = 0.5;
    quad4_xyze(2, 2) = -0.1;

    quad4_xyze(0, 3) = 0.5;
    quad4_xyze(1, 3) = 0.5;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(11);
    nids.push_back(12);
    nids.push_back(13);
    nids.push_back(14);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  {
    Core::LinAlg::SerialDenseMatrix quad4_xyze(3, 4);

    quad4_xyze(0, 0) = 0.5;
    quad4_xyze(1, 0) = 0.5;
    quad4_xyze(2, 0) = 1.1;

    quad4_xyze(0, 1) = 0.5;
    quad4_xyze(1, 1) = 0.5;
    quad4_xyze(2, 1) = -0.1;

    quad4_xyze(0, 2) = 0.5;  // 0.6;
    quad4_xyze(1, 2) = 1.1;
    quad4_xyze(2, 2) = -0.1;

    quad4_xyze(0, 3) = 0.5;  // 0.6;
    quad4_xyze(1, 3) = 1.1;
    quad4_xyze(2, 3) = 1.1;

    nids.clear();
    nids.push_back(15);
    nids.push_back(16);
    nids.push_back(17);
    nids.push_back(18);
    intersection.add_cut_side(++sidecount, nids, quad4_xyze, Core::FE::CellType::quad4);
  }

  Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

  hex8_xyze(0, 0) = 1.0;
  hex8_xyze(1, 0) = 1.0;
  hex8_xyze(2, 0) = 1.0;
  hex8_xyze(0, 1) = 1.0;
  hex8_xyze(1, 1) = 0.0;
  hex8_xyze(2, 1) = 1.0;
  hex8_xyze(0, 2) = 0.0;
  hex8_xyze(1, 2) = 0.0;
  hex8_xyze(2, 2) = 1.0;
  hex8_xyze(0, 3) = 0.0;
  hex8_xyze(1, 3) = 1.0;
  hex8_xyze(2, 3) = 1.0;
  hex8_xyze(0, 4) = 1.0;
  hex8_xyze(1, 4) = 1.0;
  hex8_xyze(2, 4) = 0.0;
  hex8_xyze(0, 5) = 1.0;
  hex8_xyze(1, 5) = 0.0;
  hex8_xyze(2, 5) = 0.0;
  hex8_xyze(0, 6) = 0.0;
  hex8_xyze(1, 6) = 0.0;
  hex8_xyze(2, 6) = 0.0;
  hex8_xyze(0, 7) = 0.0;
  hex8_xyze(1, 7) = 1.0;
  hex8_xyze(2, 7) = 0.0;

  nids.clear();
  for (int i = 0; i < 8; ++i) nids.push_back(i);

  intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);

  intersection.cut_test_cut(true, Core::Geo::Cut::VCellGaussPts_DirectDivergence);
}
