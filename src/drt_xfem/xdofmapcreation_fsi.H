/*!
\file xdofmapcreation_fsi.H

\brief defines unknowns based on the intersection pattern from the xfem intersection


<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XDOFMAPCREATION_FSI_H
#define XDOFMAPCREATION_FSI_H

#include "../drt_inpar/inpar_combust.H"

namespace XFEM
{
  void processVoidEnrichmentForElement(
      const DRT::Element*                           xfemele,
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
      const XFEM::InterfaceHandleXFSI&              ih,
      const set<int>&                               labels,
      const double                                  boundaryRatioLimit,
      const double                                  volumeRatioLimit,
      std::set<XFEM::FieldEnr>&                     elementDofs,
      bool&                                         skipped_elem_enr
      );

  /*!
   * \brief create a map that for each node (gid) and element (gid) on this processor
   * \note this is 3d fluid specific and should somehow be generalized
   */
  void createDofMapFSI(
      const XFEM::InterfaceHandleXFSI&                  ih,                  ///< interface information
      std::map<int, const std::set<XFEM::FieldEnr> >&   nodalDofSetFinal,    ///< enriched fields per node
      std::map<int, const std::set<XFEM::FieldEnr> >&   elementalDofsFinal,  ///< enriched fields per element
      const std::set<XFEM::PHYSICS::Field>&             fieldset,
      const XFEM::ElementAnsatz&                        elementAnsatz,
      const Teuchos::ParameterList&                     params,
      const vector<int>                                 MovingFluidnodeGIDs
      );

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
