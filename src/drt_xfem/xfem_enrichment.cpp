/*!
\file xfem_enrichment.cpp

\brief contains information about enrichments at nodes and elements

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "xfem_enrichment.H"

using namespace Enrichments;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const EnrId enr_id,
        const EnrType enr_type) :
    enr_id_(enr_id), enr_type_(enr_type)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
Enrichment::Enrichment(
        const Enrichment& other) :
    enr_id_(other.enr_id_), enr_type_(other.enr_type_)
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
    cout << "enr_id: "<< enr_id_ << "enr_type: "<< enr_type_ << endl;
    return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef XFEM
