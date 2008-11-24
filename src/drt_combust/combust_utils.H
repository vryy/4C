/*!---------------------------------------------------------------------
\file combust_utils.H

\brief collection of functions in namespace COMBUST

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_UTILS_H
#define COMBUST_UTILS_H

#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/interfacexfsi.H"
#include "../drt_xfem/field_enriched.H"

namespace COMBUST
{

/*------------------------------------------------------------------------------------------------*
 | print COMBUST module logo on screen                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void printlogo();

/*! 
 * \brief create a map which, for each node (gid) and element (gid) on this processor, holds the
 * corresponding enriched fields 
 * 
 * remark: Sind das nicht eigentlich einfach die dofs? Es sind ja auch Standard
 * Anreicherungen (also eben auch nicht angereicherte dofs!) dabei              .   henke 10/08
 */
void createDofMap(
    const XFEM::InterfaceHandle&                      ih,                  ///< interface information
    std::map<int, const std::set<XFEM::FieldEnr> >&   nodalDofSetFinal,    ///< enriched fields per node
    std::map<int, const std::set<XFEM::FieldEnr> >&   elementalDofsFinal   ///< enriched fields per element
);

/// apply nodal enrichments to intersected element
void ApplyNodalEnrichments(
    const DRT::Element*                           xfemele,
    const XFEM::InterfaceHandle&                  ih,
    const int&                                    label,
    std::map<int, std::set<XFEM::FieldEnr> >&     nodalDofSet
    );

/// apply standard dofs to nodes

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

} //namespace COMBUST

#endif  // #ifndef COMBUST_UTILS_H
#endif  // #ifdef CCADISCRET
