/*!
\file enrichment.cpp

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "enrichment.H"
#include "../drt_lib/drt_dserror.H"
#include <string>
#include <sstream>

using namespace XFEM;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const int id,
        const EnrType type) :
    id_(id), type_(type)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const Enrichment& other) :
    id_(other.id_), type_(other.type_)
{
    assert(&other != this);
    return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::~Enrichment()
{
    return;
}

/*----------------------------------------------------------------------*
 |  create string                                                       |
 *----------------------------------------------------------------------*/
string Enrichment::toString() const
{
    stringstream s;
    s << "Enrichment id: " << this->id_ << ", type: " << enrTypeToString(this->type_);
    return s.str();
}

string Enrichment::enrTypeToString(const EnrType type) const
{
    string typetext;
    switch (type){
        case typeStandard:  typetext = "Standard"; break;
        case typeJump:      typetext = "Jump    "; break;
        case typeVoid:      typetext = "Void    "; break;
        default: dserror("no string defined for EnrType");
    };
    return typetext;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
