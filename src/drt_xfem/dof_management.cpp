/*!
\file dof_management.cpp

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "dof_management.H"

//
// ctor 
// ag 08/07
//
XFEM::EnrPhysVar::EnrPhysVar(
        const PhysVar var,
        const Enrichment enr) :
    var_(var), enr_(enr)
{
    return;
}

//
// copy-ctor
// ag 08/07
//
XFEM::EnrPhysVar::EnrPhysVar(
        const EnrPhysVar& other) :
var_(other.var_), 
enr_(other.enr_)
{
    assert(&other != this);
    return;
}

//
// dtor
// ag 08/07
//
XFEM::EnrPhysVar::~EnrPhysVar()
{
    return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM
