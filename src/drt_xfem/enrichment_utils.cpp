/*!
\file enrichment_utils.cpp

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include <string>
#include <sstream>
#include <blitz/array.h>
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "intersection_service.H"
#include "xfem.H"
#include "physics.H"
#include "enrichment_utils.H"
#include "dof_management.H"
#include "interface.H"


using namespace XFEM;


std::map<XFEM::Enrichment, double> computeEnrvalMap(
        const RCP<XFEM::InterfaceHandle>  ih,
        const std::set<XFEM::Enrichment> enrset,
        const BlitzVec& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection
        )
{
    std::map<XFEM::Enrichment, double> enrvals;
    
    for (std::set<XFEM::Enrichment>::const_iterator enriter =
        enrset.begin(); enriter != enrset.end(); ++enriter)
    {
        const double enrval = enriter->EnrValue(actpos, ih->cutterdis(), approachdirection);
        const XFEM::Enrichment enr = *enriter;
        enrvals.insert(make_pair(enr, enrval));
    }
    return enrvals;
}



//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const blitz::Array<double,1>& funct,
        BlitzVec& enr_funct
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals = computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection);
    
    const DRT::Node*const* nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
        //const blitz::Array<double,1> nodalpos(toBlitzArray(nodes[inode]->X()));

        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                //const XFEM::Enrichment enr = enrfield->getEnrichment();
                //const double enrval = enr.ModifiedEnrValue(actpos, nodalpos, ih->cutterdis());
                //const double enrval = enr.EnrValue(actpos, ih->cutterdis(), approachdirection);
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
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
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const blitz::Array<double,1>& funct,
        const blitz::Array<double,2>& derxy,
        blitz::Array<double,1>& enr_funct,
        blitz::Array<double,2>& enr_derxy
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals = computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection);
    
    const blitz::Range _  = blitz::Range::all();
    
    const DRT::Node*const* nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
//        const blitz::Array<double,1> nodalpos(toBlitzArray(nodes[inode]->X()));

        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
//                const XFEM::Enrichment enr = enrfield->getEnrichment();
//                //const double enrval = enr.ModifiedEnrValue(actpos, nodalpos, ih->cutterdis());
//                const double enrval = enr.EnrValue(actpos, ih->cutterdis(), approachdirection);
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                dofcounter += 1;
            }
        }
    }
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const blitz::Array<double,1>& funct,
        const blitz::Array<double,2>& derxy,
        const blitz::Array<double,2>& derxy2,
        blitz::Array<double,1>& enr_funct,
        blitz::Array<double,2>& enr_derxy,
        blitz::Array<double,2>& enr_derxy2
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals = computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection);
    
    const blitz::Range _  = blitz::Range::all();
    
    const DRT::Node*const* nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
//        const blitz::Array<double,1> nodalpos(toBlitzArray(nodes[inode]->X()));

        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
//                const XFEM::Enrichment enr = enrfield->getEnrichment();
//                //const double enrval = enr.ModifiedEnrValue(actpos, nodalpos, ih->cutterdis());
//                const double enrval = enr.EnrValue(actpos, ih->cutterdis(), approachdirection);
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                enr_derxy2(_,dofcounter) = derxy2(_,inode) * enrval;
                dofcounter += 1;
            }
        }
    }
}


//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedStressShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const blitz::Array<double,1>& funct,
        blitz::Array<double,1>& enr_funct
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals = computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection);
    
    int dofcounter = 0;
    for (int inode = 0; inode < dofman.NumVirtualNodes(); ++inode)
    {
        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerElement();
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
//                const XFEM::Enrichment enr = enrfield->getEnrichment();
//                const double enrval = enr.EnrValue(actpos, ih->cutterdis(), approachdirection);
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                dofcounter += 1;                
            }
        }
    }
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedStressShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const blitz::Array<double,1>& funct,
        const blitz::Array<double,2>& derxy,
        blitz::Array<double,1>& enr_funct,
        blitz::Array<double,2>& enr_derxy
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals = computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection);
    
    const blitz::Range _  = blitz::Range::all();
    
    int dofcounter = 0;
    for (int inode = 0; inode < dofman.NumVirtualNodes(); ++inode)
    {
        const std::set<XFEM::FieldEnr> enrfieldset = dofman.FieldEnrSetPerElement();
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                //const XFEM::Enrichment enr = enrfield->getEnrichment();
                //const double enrval = enr.EnrValue(actpos, ih->cutterdis(), approachdirection);
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                dofcounter += 1;                
            }
        }
    }
}

#endif  // #ifdef CCADISCRET
