/*!
 \file xfem_condition.H

 \brief support routines

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */

#ifdef CCADISCRET
#ifndef XFEM_CONDITION_H
#define XFEM_CONDITION_H

#include <vector>
#include <set>
#include <map>

#include "../drt_lib/drt_discret.H"

namespace XFEM
{

/*!
 \brief It knows everything about the XFEM enrichments for one particular interface setup

 */
void CollectElementsByXFEMCouplingLabel(
    const DRT::Discretization&           cutterdis,
    std::map<int,set<int> >&             elementsByLabel
    );

}

#endif  // #ifndef INTERFACE_H
#endif  // #ifdef CCADISCRET
