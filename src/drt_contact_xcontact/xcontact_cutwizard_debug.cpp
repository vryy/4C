/*----------------------------------------------------------------------------*/
/**
\brief debug output methods for the xcontact cutwizard

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------------*/


#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_boundarycell.H"
#include "xcontact_cutwizard.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::CutWizard::Cut_Debug(const DRT::Element* cutted_ele) const
{
#ifdef DEBUG_XCONTACT_INTERSECTION
  std::cout << "\n\nDRT::Element #" << cutted_ele->Id() << std::endl;
  cutted_ele->Print(std::cout);

  GEO::CUT::ElementHandle* ehandle = GetElement(cutted_ele);

  std::cout << "\n\nRelated GEO::CUT::VolumeCell " << std::endl;
  std::vector<GEO::CUT::plain_volumecell_set> cell_sets;
  std::vector<std::vector<int>> nds_sets;
  std::vector<std::vector<DRT::UTILS::GaussIntegration>> intpoints_sets;
  bool has_xfem_integration_rule =
      ehandle->GetCellSets_DofSets_GaussPoints(cell_sets, nds_sets, intpoints_sets, true);
  std::cout << "has_xfem_integration_rule = " << (has_xfem_integration_rule ? "TRUE" : "FALSE")
            << "\n";

  for (std::vector<GEO::CUT::plain_volumecell_set>::const_iterator cit_sets = cell_sets.begin();
       cit_sets != cell_sets.end(); ++cit_sets)
  {
    const GEO::CUT::plain_volumecell_set& cells = *cit_sets;
    for (GEO::CUT::plain_volumecell_set::const_iterator cit = cells.begin(); cit != cells.end();
         ++cit)
    {
      GEO::CUT::VolumeCell& cell = **cit;
      cell.Print(std::cout);
    }
  }

  std::cout << "--- Nodal XFEM DoF-Sets\n";
  for (unsigned i = 0; i < nds_sets.size(); ++i)
  {
    std::cout << "Nodal DoF-set [" << i << "]:\n";
    for (unsigned j = 0; j < nds_sets[i].size(); ++j)
    {
      std::cout << "node LID [" << j << "] = " << nds_sets[i][j] << std::endl;
    }
  }

  std::cout << "\n\n--- Boundary Integrations Cells with inside position ---\n\n";
  std::vector<GEO::CUT::plain_boundarycell_set> bcellsets;
  ehandle->GetBoundaryCellSets(GEO::CUT::Point::inside, bcellsets);

  std::cout << "bcellsets.size() = " << bcellsets.size() << std::endl;
  for (GEO::CUT::plain_boundarycell_set::const_iterator ibc = bcellsets[0].begin();
       ibc != bcellsets[0].end(); ++ibc)
  {
    (*ibc)->Print();
  }
#endif
}
