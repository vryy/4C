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

using namespace XFEM;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const int id,
        const Type type) :
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
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void Enrichment::Print(
        ostream& os) const
{
    os << "Enrichment ";
    cout << "id: "<< this->id_ << "type: "<< this->type_ << endl;
    return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
