/*!
\file field_enriched.cpp

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include "field_enriched.H"



using namespace std;

/*----------------------------------------------------------------------*
 |  default ctor                                                ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr::FieldEnr() :
            field_(XFEM::PHYSICS::undefinedField),
            enr_(Enrichment())
{
    dserror("FieldEnr() -> please don't call me!");
    return;
}

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr::FieldEnr(
        const XFEM::PHYSICS::Field field,
        const Enrichment enr) :
            field_(field), enr_(enr)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                   ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr::FieldEnr(
        const FieldEnr& other) :
            field_(other.field_), 
            enr_(other.enr_)
{
    assert(&other != this);
    return;
}

/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr::~FieldEnr()
{
    return;
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
string XFEM::FieldEnr::toString() const
{
    stringstream s;
    s << "Enriched Field: " << PHYSICS::physVarToString(field_) << ", Enrichment: " << enr_.toString();
    return s.str();
}



#endif  // #ifdef CCADISCRET
