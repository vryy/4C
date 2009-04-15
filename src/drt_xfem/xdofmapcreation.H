/*!
\file xdofmapcreation.H

\brief defines unknowns based on the intersection pattern from the xfem intersection

this is related to the physics of the fluid problem and therefore not part of the standard xfem routines
somehow this should work for each physical problem differently

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
// Der ist hier nur drin, damit ich auf die FlameFront im InterfaceHandle zugreifen kann
#include "../drt_combust/combust_interface.H"

namespace XFEM
{

  bool EnrichmentInDofSet(
    const XFEM::Enrichment::EnrType     testenr,
    const std::set<XFEM::FieldEnr>&     fieldenrset);


  bool EnrichmentInNodalDofSet(
    const int                                           gid,
    const XFEM::Enrichment::EnrType                     testenr,
    const std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet);


  //! apply nodal enrichments to intersected element
  bool ApplyNodalEnrichments(
      const DRT::Element*                           xfemele,
      const XFEM::InterfaceHandle&                  ih,
      const int&                                    label,
      const XFEM::Enrichment::EnrType               enrtype,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      const double                                  volumeRatioLimit,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
      );

  //! apply nodal enrichments to intersected element: used for combustion problem only
  bool ApplyJumpEnrichmentCombust(
      const DRT::Element*                           xfemele,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      const double                                  volumeRatioLimit,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
      );

  //! apply nodal enrichments to intersected element
  bool ApplyNodalEnrichmentsNodeWise(
      const DRT::Element*                           xfemele,
      const XFEM::InterfaceHandle&                  ih,
      const int&                                    label,
      const XFEM::Enrichment::EnrType               enrtype,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      const double                                  volumeRatioLimit,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
      );


  //! apply element enrichments to intersected element
  bool ApplyElementEnrichments(
      const DRT::Element*                           xfemele,
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
      const XFEM::InterfaceHandle&                  ih,
      const int&                                    label,
      const XFEM::Enrichment::EnrType               enrtype,
      const double                                  boundaryRatioLimit,
      std::set<XFEM::FieldEnr>&                     enrfieldset
      );


  //! apply void enrichments and add this info to maps
  void ApplyVoidEnrichmentForElement(
      const DRT::Element*                           xfemele,
      const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>&  element_ansatz,
      const XFEM::InterfaceHandle&                  ih,
      const int&                                    label,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      const double                                  volumeRatioLimit,
      const double                                  boundaryRatioLimit,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
      std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs,
      bool&                                         skipped_node_enr,
      bool&                                         skipped_elem_enr
      );

  /*!
   * \brief create a map that for each node (gid) and element (gid) on this processor
   * \note this is 3d fluid specific and should somehow be generalized
   */
  void createDofMap(
      const XFEM::InterfaceHandle&                      ih,                  ///< interface information
      std::map<int, const std::set<XFEM::FieldEnr> >&   nodalDofSetFinal,    ///< enriched fields per node
      std::map<int, const std::set<XFEM::FieldEnr> >&   elementalDofsFinal,  ///< enriched fields per element
      const std::set<XFEM::PHYSICS::Field>&             fieldset,
      const XFEM::ElementAnsatz&                        elementAnsatz,
      const Teuchos::ParameterList&                     params
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

  /*!
   * \brief apply standard dofs to nodes
   * \note this is 3d fluid specific and should somehow be generalized
   */
  void applyStandardEnrichment(
      const XFEM::InterfaceHandle&                  ih,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
      std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
      );


  /*!
   * \brief apply standard dofs to nodes
   * \note this is 3d fluid specific and should somehow be generalized
   */
  void applyStandardEnrichmentNodalBasedApproach(
      const XFEM::InterfaceHandle&                  ih,
      const std::set<XFEM::PHYSICS::Field>&         fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
      );


} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
