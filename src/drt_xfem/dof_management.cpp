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
#include "integrationcell.H"
#include "dof_management.H"

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr::FieldEnr(
        const PHYSICS::Field field,
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
        const vector<XFEM::FieldEnr>& elementDofs
		) :
			nodalDofSet_(nodalDofSet),
			elementDofs_(elementDofs)
{
	// count number of dofs for each node
	std::map<int, const set<XFEM::FieldEnr> >::const_iterator tmp;
	for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp) {
		const int gid = tmp->first;
		const set<XFEM::FieldEnr> enrfieldset = tmp->second;
		const int numdof = enrfieldset.size();
		//dsassert(numdof > 0, "sollte jetzt noch nicht sein");
		nodalNumDof_[gid] = numdof;
	}
	
	// set number of parameters per field to zero
	for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp) {
		const std::set<XFEM::FieldEnr> enrfieldset = tmp->second;
		for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
			const XFEM::PHYSICS::Field field = enrfield->getField();
			numParamsPerField_[field] = 0;
			paramsLocalEntries_[field] = vector<int>();
		}
	}
    for (std::vector<XFEM::FieldEnr>::const_iterator enrfield = elementDofs.begin(); enrfield != elementDofs.end(); ++enrfield) {
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
	    
	    set<XFEM::FieldEnr>::const_iterator enrfield;
	    for (enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
	    {
	        const XFEM::PHYSICS::Field field = enrfield->getField();
	        numParamsPerField_[field] += 1;
	        paramsLocalEntries_[field].push_back(counter);
	        counter++;
	    }
	}
	// loop now over element dofs
    for (std::vector<XFEM::FieldEnr>::const_iterator enrfield = elementDofs.begin(); enrfield != elementDofs.end(); ++enrfield) {
        const XFEM::PHYSICS::Field field = enrfield->getField();
        numParamsPerField_[field] += 1;
        paramsLocalEntries_[field].push_back(counter);
    }
	
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
    XFEM::DofManager::createDofMap(
            ih->xfemdis(), ih->cutterdis(), ih->elementalDomainIntCells(),
            nodalDofSet_, elementalDofs_);
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
const XFEM::ElementDofManager XFEM::DofManager::constructElementDofManager(DRT::Element& ele) const
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
    const std::vector<XFEM::FieldEnr> elementdofs = this->getElementDofs(ele.Id());
    
    // create a local dofmanager
    XFEM::ElementDofManager eleDofManager = XFEM::ElementDofManager(ele, nodaldofset, elementdofs);

    return eleDofManager;
}

/*----------------------------------------------------------------------*
 |  sanity check                                                ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::DofManager::checkForConsistency(
        DRT::Element& ele,
        const XFEM::ElementDofManager& stored_eledofman) const
{
    // create local copy of current information about dofs
    const XFEM::ElementDofManager current_eledofman = this->constructElementDofManager(ele);
    
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
void XFEM::DofManager::createDofMap(
        const RCP<DRT::Discretization>            xfemdis,
        const RCP<DRT::Discretization>            cutterdis,
        const map<int, XFEM::DomainIntCells >&    elementDomainIntCellMap,
        map<int, const set<XFEM::FieldEnr> >&     nodalDofSetFinal,
        map<int, const vector<XFEM::FieldEnr> >&  elementalDofsFinal
        ) const
{
	// temporary assembly
    map<int, set<XFEM::FieldEnr> >  nodalDofSet;
    map<int, vector<XFEM::FieldEnr> >  elementalDofs;
    
    // standard enrichment used for all nodes (for now -> we can remove them from holes in the fluid)
    const int standard_label = 0;
    const XFEM::Enrichment enr_std(standard_label, XFEM::Enrichment::typeStandard);
    // loop my row nodes and add standard degrees of freedom
    for (int i=0; i<xfemdis->NumMyColNodes(); ++i)
    {
        const DRT::Node* actnode = xfemdis->lColNode(i);
        const int gid = actnode->Id();
        
        //if (actnode->X()[0] < 0.9 or actnode->X()[0] > 2.1)
        //{
            nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
            nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
            nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
            nodalDofSet[gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
        //}
    }

    vector< DRT::Condition * >              xfemConditions;
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
    cout << "numcondition = " << xfemConditions.size() << endl;
    
    vector< DRT::Condition * >::const_iterator condition_iter;
    for (condition_iter = xfemConditions.begin(); condition_iter != xfemConditions.end(); ++condition_iter)
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
            if (elementDomainIntCellMap.count(actele->Id()) >= 1)
            {
                const int nen = actele->NumNode();
                const int* nodeidptrs = actele->NodeIds();
                for (int inen = 0; inen<nen; ++inen)
                {
                    const int node_gid = nodeidptrs[inen];
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, jumpenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, jumpenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, jumpenr));
                    nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, jumpenr));
                    //              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::LMPLambdax, enr_void1));
                    //              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::LMPLambday, enr_void1));
                    //              nodalDofSet[node_gid].insert(XFEM::FieldEnr(PHYSICS::LMPLambdaz, enr_void1));              
                };
                // add discontinuous stress unknowns
//                int numstressparam = 0;
//                switch (actele->Shape())
//                {
//                    case DRT::Element::hex20: case DRT::Element::hex27: case DRT::Element::tet10:
//                        numstressparam = 4;
//                        break;
//                    case DRT::Element::hex8: case DRT::Element::tet4:
//                        numstressparam = 1;
//                        break;
//                    default:
//                        dserror("nope, not yet defined in createDofMap");
//                };
//                
//                const int element_gid = actele->NumNode();
//                for (int inen = 0; inen<numstressparam; ++inen)
//                {
//                    elementalDofs[element_gid].push_back(XFEM::FieldEnr(PHYSICS::Tauxx, enr));
//                    elementalDofs[element_gid].push_back(XFEM::FieldEnr(PHYSICS::Tauyy, enr));
//                    elementalDofs[element_gid].push_back(XFEM::FieldEnr(PHYSICS::Tauzz, enr));
//                    elementalDofs[element_gid].push_back(XFEM::FieldEnr(PHYSICS::Tauxy, enr));
//                    elementalDofs[element_gid].push_back(XFEM::FieldEnr(PHYSICS::Tauxz, enr));
//                    elementalDofs[element_gid].push_back(XFEM::FieldEnr(PHYSICS::Tauyz, enr));
//                }
            }
        };
    };
    
    // create const sets from standard sets, so the sets cannot be accidentily changed
    // could be removed later, if this is a performance bottleneck
    for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = nodalDofSet.begin(); oneset != nodalDofSet.end(); ++oneset )
    {
        nodalDofSetFinal.insert( make_pair(oneset->first, oneset->second));
    };
    
    for ( std::map<int, std::vector<XFEM::FieldEnr> >::const_iterator oneset = elementalDofs.begin(); oneset != elementalDofs.end(); ++oneset )
    {
        elementalDofsFinal.insert( make_pair(oneset->first, oneset->second));
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
    return assembly_type;
}


#endif  // #ifdef CCADISCRET
