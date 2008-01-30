/*!
\file dof_management.cpp

\brief provides a class that represents an enriched physical scalar field

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



using namespace std;

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
    s << "Enriched Field: " << PHYSICS::physVarToString(this->field_) << ", Enrichment: " << enr_.toString();
    return s.str();
}





/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager()
{
	return;
}

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
        DRT::Element& ele,
		const map<int, const set<XFEM::FieldEnr> >& nodalDofSet,
        const set<XFEM::FieldEnr>& elementDofs,
        const int numeleparam
		) :
			nodalDofSet_(nodalDofSet),
			elementDofs_(elementDofs),
			numeleparam_(numeleparam)
{
	// count number of dofs for each node
	map<int, const set<XFEM::FieldEnr> >::const_iterator tmp;
	for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp) {
		const int gid = tmp->first;
		const set<XFEM::FieldEnr> enrfieldset = tmp->second;
		const int numdof = enrfieldset.size();
		//dsassert(numdof > 0, "sollte jetzt noch nicht sein");
		nodalNumDof_[gid] = numdof;
	}
	
	// set number of parameters per field to zero
	for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp) {
		const set<XFEM::FieldEnr> enrfieldset = tmp->second;
		for (set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
			const XFEM::PHYSICS::Field field = enrfield->getField();
			numParamsPerField_[field] = 0;
			paramsLocalEntries_[field] = vector<int>();
		}
	}
    for (set<XFEM::FieldEnr>::const_iterator enrfield = elementDofs.begin(); enrfield != elementDofs.end(); ++enrfield) {
        const XFEM::PHYSICS::Field field = enrfield->getField();
        numParamsPerField_[field] = 0;
        paramsLocalEntries_[field] = vector<int>();
    }
	
	
	
	// count number of parameters per field
	// define local position of unknown by looping first over nodes and then over its unknowns!
	int counter = 0;
	DRT::Node** const nodes = ele.Nodes();
	for (int inode=0; inode<ele.NumNode(); inode++)
	{
	    const int gid = nodes[inode]->Id();
	    map<int, const set <XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
	    if (entry == nodalDofSet_.end())
	        dserror("impossible ;-)");
	    const set<XFEM::FieldEnr> enrfieldset = entry->second;
	    
	    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
	    {
	        const XFEM::PHYSICS::Field field = enrfield->getField();
	        numParamsPerField_[field] += 1;
	        paramsLocalEntries_[field].push_back(counter);
	        counter++;
	    }
	}
	// loop now over element dofs
	// for that we have to loop over the "virtual element nodes" asanstz function
	for (int i = 0; i < numeleparam; ++i)
    {
	    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = elementDofs.begin(); enrfield != elementDofs.end(); ++enrfield)
        {
            const XFEM::PHYSICS::Field field = enrfield->getField();
            numParamsPerField_[field] += 1;
            paramsLocalEntries_[field].push_back(counter);
            counter++;
        }
    }
	
	return;
}

/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(const ElementDofManager& old)
{
    dserror("no copying");
    return;
}
        
/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::~ElementDofManager()
{
    return;
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
std::string XFEM::ElementDofManager::toString() const
{
	stringstream s;
	map<int, const set<XFEM::FieldEnr> >::const_iterator tmp;
	for (tmp = nodalDofSet_.begin(); tmp != nodalDofSet_.end(); ++tmp)
	{
		const int gid = tmp->first;
		const set <XFEM::FieldEnr> actset = tmp->second;
		for ( set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
		{
			s << "Node: " << gid << ", " << var->toString() << endl;
		};
	};
	return s.str();
}






/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::DofManager(const RCP<XFEM::InterfaceHandle> ih) :
        	xfemdis_(ih->xfemdis())
{
    XFEM::createDofMap(
            ih->xfemdis(), ih->cutterdis(), ih->elementalDomainIntCells(),
            nodalDofSet_, elementalDofs_);
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                   ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::DofManager(const XFEM::DofManager& dofman)
{
    dserror("A DofManager shall not be copied!");
}
		
/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::~DofManager()
{
    return;
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
string XFEM::DofManager::toString() const
{
	stringstream s;
	for (int i=0; i<xfemdis_->NumMyRowNodes(); ++i)
	{
		const int gid = xfemdis_->lRowNode(i)->Id();
		const set <XFEM::FieldEnr> actset = nodalDofSet_.find(gid)->second;
		for ( set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
		{
			s << "Node: " << gid << ", " << var->toString() << endl;
		};
	};
	return s.str();
}

/*----------------------------------------------------------------------*
 |  construct element dof manager                               ag 11/07|
 *----------------------------------------------------------------------*/
const XFEM::ElementDofManager XFEM::DofManager::constructElementDofManager(DRT::Element& ele, const int numeleparam) const
{
    // create a list with number of dofs per local node 
    const int numnode = ele.NumNode();
    const int* nodegids = ele.NodeIds();
    
    // nodal dofs for ele
    map<int, const set <XFEM::FieldEnr> > nodaldofset; 
    for (int inode = 0; inode < numnode; ++inode) {
        const int gid = nodegids[inode];
        nodaldofset.insert(std::pair<int, const set<XFEM::FieldEnr> >(gid,this->getNodeDofSet(gid)));
    }

    // element dofs for ele
    const std::set<XFEM::FieldEnr> elementdofs = this->getElementDofs(ele.Id());
    
    // create a local dofmanager
    XFEM::ElementDofManager eleDofManager = XFEM::ElementDofManager(ele, nodaldofset, elementdofs, numeleparam);

    return eleDofManager;
}

/*----------------------------------------------------------------------*
 |  sanity check                                                ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::DofManager::checkForConsistency(
        DRT::Element& ele,
        const XFEM::ElementDofManager& stored_eledofman,
        const int numeleparam
        ) const
{
    // create local copy of current information about dofs
    const XFEM::ElementDofManager current_eledofman = this->constructElementDofManager(ele, numeleparam);
    
    // compare with given and report error  
    if (current_eledofman != stored_eledofman)
    {
        dserror("given elementdofmanager is not consistent with global dofmanger");
    }
    return;        
}

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

XFEM::AssemblyType XFEM::CheckForStandardEnrichmentsOnly(
        const ElementDofManager&   eleDofManager_,
        const int                  numnode,
        const int*                 nodeids
        )
{
    // find out whether we can use standard assembly or need xfem assembly
    XFEM::AssemblyType assembly_type = XFEM::standard_assembly;
    for (int inode = 0; inode < numnode; ++inode)
    {
        if (assembly_type == XFEM::xfem_assembly)
        {
            break;
        }
        const int gid = nodeids[inode];
        std::set<XFEM::FieldEnr> fields = eleDofManager_.FieldEnrSetPerNode(gid);
        if (fields.size() != 4)
        {
            assembly_type = XFEM::xfem_assembly;
            break;
        };
        for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fields.begin(); fieldenr != fields.end(); ++fieldenr)
        {
            if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
            {
                assembly_type = XFEM::xfem_assembly;
                break;
            };
        };
    };
    const int eledof = eleDofManager_.NumDofPerElement();
    if (eledof != 0)
    {
        assembly_type = XFEM::xfem_assembly;
    }
    
    return assembly_type;
}


#endif  // #ifdef CCADISCRET
