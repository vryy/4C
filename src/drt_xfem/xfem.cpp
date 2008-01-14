/*!
\file xfem.cpp

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "xfem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_node.H"
#include <string>
#include <sstream>


inline blitz::Array<double,1> XFEM::toBlitzArray(const double* x)
{
    blitz::Array<double,1> nodalpos(3);
    nodalpos(0) = x[0];
    nodalpos(1) = x[1];
    nodalpos(2) = x[2];
    return nodalpos;
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& funct,
        blitz::Array<double,1>& enr_funct
        )
{
    
    DRT::Node** const nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); inode++)
    {
        const int gid = nodes[inode]->Id();
        const blitz::Array<double,1> nodalpos(toBlitzArray(nodes[inode]->X()));

        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const XFEM::Enrichment enr = enrfield->getEnrichment();
                const double enrval = enr.ModifiedEnrValue(actpos, nodalpos, ih->cutterdis());
                enr_funct(dofcounter) = funct(inode) * enrval;
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& funct,
        const blitz::Array<double,2>& derxy,
        const blitz::Array<double,2>& derxy2,
        blitz::Array<double,1>& enr_funct,
        blitz::Array<double,2>& enr_derxy,
        blitz::Array<double,2>& enr_derxy2
        )
{
    blitz::Range _  = blitz::Range::all();
    
    DRT::Node** const nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); inode++)
    {
        const int gid = nodes[inode]->Id();
        const blitz::Array<double,1> nodalpos(toBlitzArray(nodes[inode]->X()));

        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const XFEM::Enrichment enr = enrfield->getEnrichment();
                const double enrval = enr.ModifiedEnrValue(actpos, nodalpos, ih->cutterdis());
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                enr_derxy2(_,dofcounter) = derxy2(_,inode) * enrval;
                dofcounter += 1;
            }
        }
    }
}


#endif  // #ifdef CCADISCRET
