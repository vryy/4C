/*!
\file xdofmapcreation.H

\brief defines unknowns based on the intersection pappern from the xfem intersection

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

#include "xfem.H"
//#include "physics.H"
//#include "enrichment.H"
#include "integrationcell.H"
//#include "interface.H"
#include "../drt_lib/drt_discret.H"


// namespace for all xfem related stuff
namespace XFEM
{

// create a map that for each node (gid) and element (gid) on this processor
// this is fluid specific and should somehow be generalized
void createDofMap(
        const RefCountPtr<DRT::Discretization>    xfemdis,
        const RefCountPtr<DRT::Discretization>    cutterdis,
        const map<int, DomainIntCells >&          elementDomainIntCellMap,
        map<int, const set<XFEM::FieldEnr> >&     nodalDofSetFinal,        ///< the map for nodal dofs
        map<int, const set<XFEM::FieldEnr> >&     elementalDofsFinal       ///< the map for elemental dofs        
);

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
