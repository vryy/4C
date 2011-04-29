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
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );
  
  // schott Aug 3, 2010
  //! apply nodal enrichments to touched element: used for combustion and Nitsche problem
  bool ApplyJumpEnrichmentToTouched(
      const COMBUST::InterfaceHandleCombust&    ih,
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  //! apply nodal enrichments to intersected element: used for two phase problem only
  bool ApplyKinkEnrichment(
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );
  
  // schott Jun 15, 2010
  //! apply nodal enrichments to intersected element: used for two phase problem with surface tension only
  bool ApplyKinkJumpEnrichment(
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );
  
  // schott Aug 3, 2010
  //! apply nodal enrichments to touched element: used for two phase problem with surface tension only
  bool ApplyKinkJumpEnrichmentToTouched(
      const COMBUST::InterfaceHandleCombust&    ih,
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
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
      const XFEM::ElementAnsatz&                elementAnsatz, ///< approach for extra stress field
      const Teuchos::ParameterList&             params         ///< XFEM parameter list
      );

  /*!
   * \brief apply standard dofs to nodes
   * \note preliminary split-off function from applyStandardEnrichment() for combustion problems
   */
  void ApplyStandardEnrichmentCombust(
      const DRT::Element*                           xfemele,   ///<
      const std::set<XFEM::PHYSICS::Field>&         fieldset,  ///<
      std::map<int, std::set<XFEM::FieldEnr> >&     nodeDofMap ///<
      );

    /*!
   * \brief apply standard dofs to nodes
   * \note preliminary split-off function from applyStandardEnrichment() for combustion problems
   */
  void ApplyElementEnrichmentCombust(
      const DRT::Element*                                                xfemele,           ///<
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,    ///< element dofs
      std::map<int, std::set<XFEM::FieldEnr> >&                          elementDofMap,     ///<
      const COMBUST::InterfaceHandleCombust&                             ih                 ///<
      );

  void ApplyElementEnrichmentCombust(
      const DRT::Element*                                                xfemele,
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
      std::set<XFEM::FieldEnr>&                                          elementFieldEnrSet,
      const COMBUST::InterfaceHandleCombust&                             ih
      );

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
