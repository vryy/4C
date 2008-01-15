/*!
\file xfem.H

\brief general xfem routines without a class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFEM_H
#define XFEM_H

#include <blitz/array.h>
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "physics.H"
#include "dof_management.H"
#include "interface.H"

using namespace Teuchos;


namespace XFEM
{

typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix

enum IntegrationType
{
    sub_integration,
    standard_integration
};

inline BlitzVec toBlitzArray(const double* x)
{
    BlitzVec blitz_x(3);
    blitz_x(0) = x[0];
    blitz_x(1) = x[1];
    blitz_x(2) = x[2];
    
    // doesn't work... why?
    // BlitzVec blitz_x(x, blitz::shape(3), blitz::duplicateData);
    return blitz_x;
}

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const BlitzVec& funct,
        BlitzVec& enr_funct
        );

//
// For a given situation compute the enrichement values
// 
inline BlitzVec ComputeEnrichmentArrayForNodalDofs(
        DRT::Element&  ele,
        const int numparamvelx,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos
        )
{
    
    BlitzVec enrvals(numparamvelx);
    
    DRT::Node** const nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); inode++)
    {
        const DRT::Node* node = nodes[inode]; 
        const int gid = node->Id();
        const BlitzVec nodalpos(toBlitzArray(node->X()));

        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const XFEM::Enrichment enr = enrfield->getEnrichment();
                enrvals(dofcounter) = enr.ModifiedEnrValue(actpos, nodalpos, ih->cutterdis());
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
    return enrvals;
}

void ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        const BlitzMat& derxy2,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy,
        BlitzMat& enr_derxy2
        );

void ComputeEnrichedStressShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy
        );

inline std::map<Enrichment,double> getEnrichmentValues(
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const BlitzVec& actpos
        )
{
    std::map<Enrichment,double> result;
    
    std::set<Enrichment> enrset= dofman.getUniqueEnrichmentSet();
    for (std::set<Enrichment>::const_iterator enr = enrset.begin(); enr != enrset.end(); ++enr)
    {
        std::pair<Enrichment,double> p((*enr), enr->EnrValue(actpos, ih->cutterdis()));
        result.insert(p);
    }
        
    return result;
}

}


#endif  // #ifndef XFEM_H
#endif  // #ifdef CCADISCRET
