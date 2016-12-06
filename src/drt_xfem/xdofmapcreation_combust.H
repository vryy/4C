/*!
\file xdofmapcreation_combust.H

\brief defines unknowns based on the intersection pattern from the xfem intersection

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#ifndef XDOFMAPCREATION_COMBUST_H
#define XDOFMAPCREATION_COMBUST_H


#include "../drt_inpar/inpar_combust.H"
#include "../drt_lib/drt_element.H"
#include "physics.H"


namespace COMBUST
{
  class InterfaceHandleCombust;
}
namespace DRT
{
  class Element;
}


namespace XFEM
{
  class ElementAnsatz;
  class FieldEnr;

  //! apply nodal enrichments to intersected element: used for combustion problem only
  bool ApplyJumpEnrichment(
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  //! apply nodal enrichments to touched element: used for combustion and Nitsche problem
  bool ApplyJumpEnrichmentToTouched(
      const COMBUST::InterfaceHandleCombust&    ih,
      const Epetra_Vector*                      phinp,
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  //! apply nodal enrichments to intersected element: used for two phase problem only
  bool ApplyKinkEnrichment(
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      const INPAR::COMBUST::SelectedEnrichment& selectedEnrichment,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  //! apply nodal enrichments to intersected element: used for two phase problem with surface tension only
  bool ApplyKinkJumpEnrichment(
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      const INPAR::COMBUST::SelectedEnrichment& selectedEnrichment,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  // schott Aug 3, 2010
  //! apply nodal enrichments to touched element: used for two phase problem with surface tension only
  bool ApplyKinkJumpEnrichmentToTouched(
      const COMBUST::InterfaceHandleCombust&    ih,
      const Epetra_Vector*                      phinp,
      const DRT::Element*                       xfemele,
      const std::set<XFEM::PHYSICS::Field>&     fieldset,
      const INPAR::COMBUST::SelectedEnrichment& selectedEnrichment,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap
      );

  /*!
   * \brief build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
   * \note preliminary split-off function from createDofMap() for combustion problems
   */
  void createDofMapCombust(
      const COMBUST::InterfaceHandleCombust&    ih,            ///< interface information
      const Epetra_Vector*                      phinp,
      std::map<int, std::set<XFEM::FieldEnr> >& nodeDofMap,    ///< enriched fields per node
      const std::set<XFEM::PHYSICS::Field>&     fieldset,      ///< physical fields to enrich
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

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
