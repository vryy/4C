/*!
\file xdofmapcreation_combust.H

\brief defines unknowns based on the intersection pattern from the xfem intersection


<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XDOFMAPCREATION_COMBUST_H
#define XDOFMAPCREATION_COMBUST_H

#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/interface.H"
#include "../drt_xfem/element_ansatz.H"
#include "../drt_inpar/inpar_combust.H"
// Der ist hier nur drin, damit ich (Florian Henke) auf die FlameFront im InterfaceHandle zugreifen kann
#include "../drt_combust/combust_interface.H"

namespace XFEM
{

  //! apply nodal enrichments to intersected element: used for combustion problem only
  bool ApplyJumpEnrichment(
      const DRT::Element*                           xfemele,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      const double                                  volumeRatioLimit,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
      );

  //! apply nodal enrichments to intersected element: used for two phase problem only
  bool ApplyKinkEnrichment(
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      const double                              volumeRatioLimit,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  /*!
   * \brief build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
   * \note preliminary split-off function from createDofMap() for combustion problems
   */
  void createDofMapCombust(
      const COMBUST::InterfaceHandleCombust&    ih,            ///< interface information
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap,    ///< enriched fields per node
      std::map<int, std::set<XFEM::FieldEnr> >& elementDofMap, ///< enriched fields per element
      const std::set<XFEM::PHYSICS::Field>&     fieldset,      ///< physical fields to enrich
      const Teuchos::ParameterList&             params         ///< XFEM parameter list
      );

  /*!
   * \brief apply standard dofs to nodes
   * \note preliminary split-off function from applyStandardEnrichment() for combustion problems
   */
  void ApplyStandardEnrichmentCombust(
      const DRT::Element*                           xfemele,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
      );

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
