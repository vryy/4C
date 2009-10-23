/*!
\file xdofmapcreation.cpp

\brief defines unknowns based on the intersection pattern from the xfem intersection

this is related to the physics of the fluid problem and therefore should not be part of the standard xfem routines

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include <algorithm>
#include <set>
#include <iterator>


#include "xdofmapcreation.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_utils.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::EnrichmentInDofSet(
    const XFEM::Enrichment::EnrType     testenr,
    const std::set<XFEM::FieldEnr>&     fieldenrset)
{
  bool enr_in_set = false;
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
  {
    if (fieldenr->getEnrichment().Type() == testenr)
    {
      enr_in_set = true;
      break;
    }
  }
  return enr_in_set;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool XFEM::EnrichmentInNodalDofSet(
    const int                                           gid,
    const XFEM::Enrichment::EnrType                     testenr,
    const std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet)
{
  bool enr_at_node = false;
  //check for testenrichment in the given nodalDofSet
  std::map<int, std::set<XFEM::FieldEnr> >::const_iterator setiter = nodalDofSet.find(gid);
  if (setiter != nodalDofSet.end())
  {
    const std::set<XFEM::FieldEnr>& fieldenrset = setiter->second;
    enr_at_node = XFEM::EnrichmentInDofSet(testenr, fieldenrset);
  }
  return enr_at_node;
}




//! check whether to sets of NodeIds are overlapping or not
bool XFEM::ConnectedElements(
    const set<int>& patchset,
    const set<int>& testset
    )
{
  std::set<int> result;
  std::set_intersection(patchset.begin(), patchset.end(), testset.begin(), testset.end(),
    std::insert_iterator<set<int> >(result, result.begin()));
  return result.size() > 0;
}

void XFEM::separateByLabel(
    const XFEM::InterfaceHandleXFSI&             ih,
    const std::set<int>&                         beles,
    map<int, set<int> >&                         elementsByLabel,
    map<int, set<int> >&                         nodesByLabel
    )
{
    for (std::set<int>::const_iterator iele=beles.begin();iele != beles.end(); ++iele)
    {
      const int label = ih.GetLabelPerBoundaryElementId(*iele);
      elementsByLabel[label].insert(*iele);
    }

    // translate to node set
    for (map<int, set<int> >::const_iterator entry = elementsByLabel.begin(); entry != elementsByLabel.end();++entry)
    {
      const int label = entry->first;
      const set<int> eleGiDs = entry->second;
      for (set<int>::const_iterator ibele = eleGiDs.begin(); ibele != eleGiDs.end(); ++ibele)
      {
        const DRT::Element* bele = ih.cutterdis()->gElement(*ibele);
        const int* bnodeids = bele->NodeIds();
        for (int inode = 0; inode < bele->NumNode(); ++inode)
        {
          nodesByLabel[label].insert(bnodeids[inode]);
        }
      }
    }
}




#endif  // #ifdef CCADISCRET
