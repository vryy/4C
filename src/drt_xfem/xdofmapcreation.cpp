/*!
\file xdofmapcreation.cpp

\brief defines unknowns based on the intersection pappern from the xfem intersection

this is related to the physics of the fluid problem and therefore not part of the standard xfem routines

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "dof_management.H"
#include "xfem.H"
#include "xdofmapcreation.H"

using namespace std;



/*----------------------------------------------------------------------*
 |  construct dofmap                                            ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::createDofMap(
        const RCP<DRT::Discretization>            xfemdis,
        const RCP<DRT::Discretization>            cutterdis,
        const map<int, XFEM::DomainIntCells >&    elementDomainIntCellMap,
        map<int, const set<XFEM::FieldEnr> >&     nodalDofSetFinal,
        map<int, const set<XFEM::FieldEnr> >&     elementalDofsFinal
        )
{
    // temporary assembly
    map<int, set<XFEM::FieldEnr> >  nodalDofSet;
    map<int, set<XFEM::FieldEnr> >  elementalDofs;
    
//    // loop my row nodes and add standard degrees of freedom to nodes
//    for (int i=0; i<xfemdis->NumMyColNodes(); ++i)
//    {
//        // standard enrichment used for all nodes (for now -> we can remove them from holes in the fluid)
//        const int standard_label = 0;    
//        const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
//        const DRT::Node* actnode = xfemdis->lColNode(i);
//        const int gid = actnode->Id();
//        nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
//        nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
//        nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
//        nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
//    };
    
//    const int standard_label = 0;
//    const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
//    for (int i=0; i<xfemdis->NumMyColElements(); ++i)
//    {
//        DRT::Element* actele = xfemdis->lColElement(i);
//        if ( not (elementDomainIntCellMap.count(actele->Id()) >= 1))
//        {
//            const int* nodeidptrs = actele->NodeIds();
//            const Epetra_SerialDenseVector nodalpos = toEpetraArray(actele->Nodes()[0]->X());
//            
//            const bool in_solid = XFEM::PositionWithinDiscretization(cutterdis, nodalpos);
//            
//            if (not in_solid)
//            {
//                for (int inen = 0; inen<actele->NumNode(); ++inen)
//                {
//                    const int node_gid = nodeidptrs[inen];
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
//                };
//            }
////             // add continuous stress unknowns
////                const int numstressparam = XFLUID::getNumberOfStressDofs(actele->Shape());
////                const int element_gid = actele->Id();
////                for (int inen = 0; inen<numstressparam; ++inen)
////                {
////                    elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxx, voidenr));
////                    elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyy, voidenr));
////                    elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauzz, voidenr));
////                    elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxy, voidenr));
////                    elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxz, voidenr));
////                    elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyz, voidenr));
////                }
//        }
//    };
    
    // loop xfem conditions and add appropriate enrichments
    vector< DRT::Condition * >              xfemConditions;
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
    cout << "numcondition = " << xfemConditions.size() << endl;
    
    for (vector< DRT::Condition * >::const_iterator condition_iter = xfemConditions.begin(); condition_iter != xfemConditions.end(); ++condition_iter)
    {
        const DRT::Condition * condition = *condition_iter;
        const int label = condition->Getint("label");
        cout << "condition with label = " << label << endl;
    
        // for surface 1, loop my col elements and add void enrichments to each elements member nodes
        const XFEM::Enrichment jumpenr(label, XFEM::Enrichment::typeJump);
        const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);
        for (int i=0; i<xfemdis->NumMyColElements(); ++i)
        {
            const DRT::Element* actele = xfemdis->lColElement(i);
            const int element_gid = actele->Id();
            if (elementDomainIntCellMap.count(element_gid) >= 1)
            {
                //TODO: check if element is intersected by the CURRENT condition label
                
                const int nen = actele->NumNode();
                const int* nodeidptrs = actele->NodeIds();
                for (int inen = 0; inen<nen; ++inen)
                {
                    const int node_gid = nodeidptrs[inen];
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, jumpenr));
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, jumpenr));
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, jumpenr));
//                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, jumpenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, voidenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, voidenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, voidenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, voidenr));
                };
                // add discontinuous stress unknowns
                // the number of each of these parameters will be determined later
                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxx, voidenr));
                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyy, voidenr));
                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauzz, voidenr));
                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxy, voidenr));
                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxz, voidenr));
                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyz, voidenr));
            }
        };
    };
    
    
    const int standard_label = 0;
    const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
    for (int i=0; i<xfemdis->NumMyColElements(); ++i)
    {
        DRT::Element* actele = xfemdis->lColElement(i);
        if ( not (elementDomainIntCellMap.count(actele->Id()) >= 1))
        {
            const int* nodeidptrs = actele->NodeIds();
            const Epetra_SerialDenseVector nodalpos = toEpetraArray(actele->Nodes()[0]->X());
            
            const bool in_solid = XFEM::PositionWithinDiscretization(cutterdis, nodalpos);
            
            if (not in_solid)
            {
                for (int inen = 0; inen<actele->NumNode(); ++inen)
                {
                    const int node_gid = nodeidptrs[inen];
                    bool voidenrichment_in_set = false;
                    //check for void enrichement in a given set, if such set already exists for this node_gid
                    std::map<int, std::set<FieldEnr> >::const_iterator setiter = nodalDofSet.find(node_gid);
                    if (setiter != nodalDofSet.end())
                    {
                        
                        std::set<FieldEnr> fieldenrset = setiter->second;
                        for (std::set<FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
                        {
                            if (fieldenr->getEnrichment().Type() == Enrichment::typeVoid)
                            {
                                voidenrichment_in_set = true;
                                break;
                            }
                        }
                    }
                    if (not voidenrichment_in_set)
                    {
                        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
                        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
                        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
                        nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
                    }
                };

                // add continuous stress unknowns
//                const int element_gid = actele->Id();
//                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxx, enr_std));
//                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyy, enr_std));
//                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauzz, enr_std));
//                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxy, enr_std));
//                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauxz, enr_std));
//                elementalDofs[element_gid].insert(XFEM::FieldEnr(PHYSICS::Tauyz, enr_std));
            }
        }
    };
//    if (stress_unknowns_used)
//    {
//        cout << (elementalDofs[0].size() / 6) << " continuous stress unknown(s) applied to " << elementalDofs.size() << " elements." << endl;
////        const set<XFEM::FieldEnr> fieldset = elementalDofs[0];
////        for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldset.begin(); fieldenr != fieldset.end(); ++fieldenr)
////        {
////            cout << fieldenr->toString() << endl;
////        }
//    }

//    for (vector< DRT::Condition * >::const_iterator condition_iter = xfemConditions.begin(); condition_iter != xfemConditions.end(); ++condition_iter)
//    {
//        const DRT::Condition * condition = *condition_iter;
//        const int label = condition->Getint("label");
//        cout << "condition with label = " << label << endl;
//    
//        // for surface 1, loop my col elements and add void enrichments to each elements member nodes
//        const XFEM::Enrichment jumpenr(label, XFEM::Enrichment::typeJump);
//        const XFEM::Enrichment voidenr(label, XFEM::Enrichment::typeVoid);
//        for (int i=0; i<xfemdis->NumMyColElements(); ++i)
//        {
//            const DRT::Element* actele = xfemdis->lColElement(i);
//            if (elementDomainIntCellMap.count(actele->Id()) >= 1)
//            {
//                const int nen = actele->NumNode();
//                //const int* nodeidptrs = actele->NodeIds();
//                for (int inen = 0; inen<nen; ++inen)
//                {
////                    const int node_gid = nodeidptrs[inen];
////                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, jumpenr));
////                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, jumpenr));
////                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, jumpenr));
////                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, jumpenr));
//                };
//            }
//        };
//    };
    
    // create const sets from standard sets, so the sets cannot be accidentily changed
    // could be removed later, if this is a performance bottleneck
    for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = nodalDofSet.begin(); oneset != nodalDofSet.end(); ++oneset )
    {
        nodalDofSetFinal.insert( make_pair(oneset->first, oneset->second));
    };
    
    for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator onevec = elementalDofs.begin(); onevec != elementalDofs.end(); ++onevec )
    {
        elementalDofsFinal.insert( make_pair(onevec->first, onevec->second));
    };
}


#endif  // #ifdef CCADISCRET
