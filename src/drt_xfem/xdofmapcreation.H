/*!
\file xdofmapcreation.H

\brief defines unknowns based on the intersection pattern from the xfem intersection


<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XDOFMAPCREATION_H
#define XDOFMAPCREATION_H

#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/interface.H"
#include "../drt_xfem/element_ansatz.H"

namespace XFEM
{

  bool EnrichmentInDofSet(
    const XFEM::Enrichment::EnrType     testenr,
    const std::set<XFEM::FieldEnr>&     fieldenrset);


  bool EnrichmentInNodalDofSet(
    const int                                           gid,
    const XFEM::Enrichment::EnrType                     testenr,
    const std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet);

  //! check whether to sets of NodeIds are overlapping or not
  bool ConnectedElements(
      const set<int>& patchset,
      const set<int>& testset
      );

  void separateByLabel(
      const XFEM::InterfaceHandleXFSI&             ih,
      const std::set<int>&                         beles,
      map<int, set<int> >&                         elementsByLabel,
      map<int, set<int> >&                         nodesByLabel
      );

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
