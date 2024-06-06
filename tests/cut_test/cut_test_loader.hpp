/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#ifndef CUT_TEST_LOADER_HPP
#define CUT_TEST_LOADER_HPP

#include "4C_cut_meshintersection.hpp"

#include <map>
#include <vector>

using namespace FourC;

class MeshLoader
{
 public:
  void GetCutNode(int nid, double x, double y, double z, double lsv);

  void GetNode(int nid, double x, double y, double z, double lsv);

  void create_side(int sid, int nid1, int nid2, int nid3, int nid4, Core::FE::CellType shape);

  void create_element(int eid, int nid1, int nid2, int nid3, int nid4, int nid5, int nid6, int nid7,
      int nid8, Core::FE::CellType);

  void CutTest_Cut(bool include_inner, bool do_Cut_Positions_Dofsets = false)
  {
    mesh_.GetOptions().Init_for_Cuttests();  // use full cln
    mesh_.CutTest_Cut(include_inner, Inpar::Cut::VCellGaussPts_DirectDivergence,
        Inpar::Cut::BCellGaussPts_Tessellation, true, true, do_Cut_Positions_Dofsets);
    mesh_.print_cell_stats();
  }

 private:
  void fill(std::map<int, std::vector<double>>& nodes, int nid, double* values);

  Core::Geo::Cut::MeshIntersection mesh_;
  std::map<int, std::vector<double>> nodes_;
  std::map<int, std::vector<double>> cut_nodes_;
};
#endif
