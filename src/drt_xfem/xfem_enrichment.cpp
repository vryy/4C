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
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "xfem_enrichment.H"
#include "../drt_lib/drt_dserror.H"

using namespace Enrichments;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
ElementEnrichment::ElementEnrichment(const int enr_id, const EnrichmentType enr_type) :
enr_id_(enr_id),
enr_type_(enr_type)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
ElementEnrichment::ElementEnrichment(const ElementEnrichment& other) :
enr_id_(other.enr_id_),
enr_type_(other.enr_type_)
{
    assert(&other != this);
    return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
ElementEnrichment::~ElementEnrichment()
{
    return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void ElementEnrichment::Print(ostream& os) const
{
    os << "ElementEnrichment ";
    cout << "enr_id: " << enr_id_ << "enr_type: " << enr_type_ << endl;
    return;
}






/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
NodeEnrichment::NodeEnrichment(const int enr_id, const EnrichmentType enr_type) :
enr_id_(enr_id),
enr_type_(enr_type)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
NodeEnrichment::NodeEnrichment(const NodeEnrichment& other) :
enr_id_(other.enr_id_),
enr_type_(other.enr_type_)
{
    assert(&other != this);
    return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
NodeEnrichment::~NodeEnrichment()
{
    return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void NodeEnrichment::Print(ostream& os) const
{
    os << "NodeEnrichment ";
    cout << "enr_id: " << enr_id_ << "enr_type: " << enr_type_ << endl;
    return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
