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

#include "dof_management.H"
#include "interfacexfsi.H"


namespace XFEM
{
/// apply nodal enrichments to intersected element
void ApplyNodalEnrichments(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const XFEM::Enrichment&                       voidenr,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
    );

/// apply element enrichments to intersected element
void ApplyElementEnrichments(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const XFEM::Enrichment&                       voidenr,
    const bool                                    DLM_condensation,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
    );

/// apply enrichments and add this info to maps
void ApplyVoidEnrichmentForElement(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const std::map<int,int>&                      labelPerElementId,
    const XFEM::Enrichment&                       voidenr,
    const bool                                    DLM_condensation,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
    );

/*! 
 * \brief create a map that for each node (gid) and element (gid) on this processor
 * \note this is 3d fluid specific and should somehow be generalized
 */
void createDofMap(
    const XFEM::InterfaceHandle&                      ih,                  ///< interface information
    std::map<int, const std::set<XFEM::FieldEnr> >&   nodalDofSetFinal,    ///< enriched fields per node
    std::map<int, const std::set<XFEM::FieldEnr> >&   elementalDofsFinal,  ///< enriched fields per element
    const bool DLM_condensation
);

/*! 
 * \brief apply standard dofs to nodes
 * \note this is 3d fluid specific and should somehow be generalized
 */
void applyStandardEnrichment(
    const XFEM::InterfaceHandle&                  ih,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
);

/*! 
 * \brief apply standard dofs to nodes
 * \note this is 3d fluid specific and should somehow be generalized
 */
void applyStandardEnrichmentNodalBasedApproach(
    const XFEM::InterfaceHandle&                  ih,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet,
    std::map<int, std::set<XFEM::FieldEnr> >&     elementalDofs
);

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
