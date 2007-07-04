/*!----------------------------------------------------------------------
\file xfem_enrichment.H

\brief contains information about enrichments at nodes and elements

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_FLUID3_XFEM_ENRICHMENT_H
#define DRT_FLUID3_XFEM_ENRICHMENT_H

#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;

// namespace for all enrichment related stuff
namespace Enrichments
{

enum EnrichmentType
{
    enrichment_type_standard = 0,
    enrichment_type_jump     = 1
};


/*!
\class ElementEnrichment

\brief Corresponds to a specific set of nodal enrichments with unique id and type

*/
class ElementEnrichment
{
public:

    ElementEnrichment(const int enr_id, const EnrichmentType enr_type);
    ElementEnrichment(const ElementEnrichment& old);
    ~ElementEnrichment();
    void Print(ostream& os) const;
   
    inline int getId() const {return enr_id_;};
    inline int getType() const {return enr_type_;};

private:
    // a global id used to distuinguish several surfaces within one element
    int             enr_id_; 
    // the specific type of enrichment
    EnrichmentType  enr_type_;
}; // class ElementEnrichment

/*!
\class NodeEnrichment

\brief Corresponds to a specific nodal enrichment with unique id and type

*/
class NodeEnrichment
{
public:

    NodeEnrichment(const int enr_id, const EnrichmentType enr_type);
    NodeEnrichment(const NodeEnrichment& other);
    ~NodeEnrichment();
    void Print(ostream& os) const;
   
    inline int getId() const {return enr_id_;};
    inline int getType() const {return enr_type_;};

private:
    // a global id used to distuinguish several surfaces within one element
    int             enr_id_; 
    // the specific type of enrichment
    EnrichmentType  enr_type_;
}; // class NodeEnrichment

}; // namespace Enrichments

#endif  // #ifndef FLUID3_XFEM_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
