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
#ifdef CCADISCRET

#include "dof_management.H"

//
// ctor 
// ag 08/07
//
XFEM::EnrPhysVar::EnrPhysVar(
        const PhysVar physvar,
        const Enrichment enr) :
            physvar_(physvar), enr_(enr)
{
    return;
}

//
// copy-ctor
// ag 08/07
//
XFEM::EnrPhysVar::EnrPhysVar(
        const EnrPhysVar& other) :
            physvar_(other.physvar_), 
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

string XFEM::EnrPhysVar::toString() const
{
    stringstream s;
    s << "Enriched PhysVar: " << Physics::physVarToString(this->physvar_) << ", Enrichment: " << enr_.toString();
    return s.str();
}

#endif  // #ifdef CCADISCRET
