/*----------------------------------------------------------------------*/
/*! \file
\brief Test for the CUT Library

\level 1

*----------------------------------------------------------------------*/

#ifndef CUT_TEST_LOADER_H
#define CUT_TEST_LOADER_H

#include <map>
#include <vector>

#include "../../src/drt_cut/cut_meshintersection.H"

class MeshLoader
{
 public:
  void GetCutNode(int nid, double x, double y, double z, double lsv);

  void GetNode(int nid, double x, double y, double z, double lsv);

  void CreateSide(
      int sid, int nid1, int nid2, int nid3, int nid4, DRT::Element::DiscretizationType shape);

  void CreateElement(int eid, int nid1, int nid2, int nid3, int nid4, int nid5, int nid6, int nid7,
      int nid8, DRT::Element::DiscretizationType);

  void CutTest_Cut(bool include_inner, bool do_Cut_Positions_Dofsets = false)
  {
    mesh_.GetOptions().Init_for_Cuttests();  // use full cln
    mesh_.CutTest_Cut(include_inner, INPAR::CUT::VCellGaussPts_DirectDivergence,
        INPAR::CUT::BCellGaussPts_Tessellation, true, true, do_Cut_Positions_Dofsets);
    mesh_.PrintCellStats();
  }

 private:
  void Fill(std::map<int, std::vector<double>>& nodes, int nid, double* values);

  GEO::CUT::MeshIntersection mesh_;
  std::map<int, std::vector<double>> nodes_;
  std::map<int, std::vector<double>> cut_nodes_;
};
#endif
