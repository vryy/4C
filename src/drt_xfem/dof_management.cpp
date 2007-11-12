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

#include "integrationcell.H"
#include "dof_management.H"

//
// ctor 
// ag 08/07
//
XFEM::EnrPhysVar::EnrPhysVar(
        const PhysVar physvar,
        const Enrichment enr) :
            physvar_(physvar), enr_(enr)
{
    return;
}

//
// copy-ctor
// ag 08/07
//
XFEM::EnrPhysVar::EnrPhysVar(
        const EnrPhysVar& other) :
            physvar_(other.physvar_), 
            enr_(other.enr_)
{
    assert(&other != this);
    return;
}

//
// dtor
// ag 08/07
//
XFEM::EnrPhysVar::~EnrPhysVar()
{
    return;
}

string XFEM::EnrPhysVar::toString() const
{
    stringstream s;
    s << "Enriched PhysVar: " << Physics::physVarToString(this->physvar_) << ", Enrichment: " << enr_.toString();
    return s.str();
}


//
// ctor 
// ag 08/07
//
XFEM::DofManager::DofManager(
		RefCountPtr<DRT::Discretization> xfemdis,
        const map<int, DomainIntCells >&  elementalDomainIntCells) :
        	xfemdis_(xfemdis)
{
	nodalDofMap_ = XFEM::DofManager::createNodalDofMap(xfemdis, elementalDomainIntCells);
}
		
//
// dtor
// ag 08/07
//
XFEM::DofManager::~DofManager()
{
    return;
}




// debug: print enrichments to screen
string XFEM::DofManager::toString() const
{
	stringstream s;
	for (int i=0; i<xfemdis_->NumMyRowNodes(); ++i)
	{
		const int gid = xfemdis_->lRowNode(i)->Id();
		const set <XFEM::EnrPhysVar> actset = nodalDofMap_.find(gid)->second;
		for ( set<XFEM::EnrPhysVar>::const_iterator var = actset.begin(); var != actset.end(); ++var )
		{
			s << "Node: " << gid << ", " << var->toString() << endl;
		};
	};
	return s.str();
}


map<int, const set <XFEM::EnrPhysVar> > XFEM::DofManager::createNodalDofMap(
        const RefCountPtr<DRT::Discretization>        xfemdis,
        const map<int, XFEM::DomainIntCells >&  elementDomainIntCellMap) const
{
	// the final map
	map<int, const set <XFEM::EnrPhysVar> >  nodalDofMapFinal;

	// temporary assembly
    map<int, set <XFEM::EnrPhysVar> >  nodalDofMap;
    
    // loop my row nodes and add standard degrees of freedom
    const XFEM::Enrichment enr_std(0, XFEM::Enrichment::typeStandard);
    for (int i=0; i<xfemdis->NumMyRowNodes(); ++i)
    {
        const DRT::Node* actnode = xfemdis->lRowNode(i);
        nodalDofMap[actnode->Id()].insert(XFEM::EnrPhysVar(Physics::Velx, enr_std));
        nodalDofMap[actnode->Id()].insert(XFEM::EnrPhysVar(Physics::Vely, enr_std));
        nodalDofMap[actnode->Id()].insert(XFEM::EnrPhysVar(Physics::Velz, enr_std));
        nodalDofMap[actnode->Id()].insert(XFEM::EnrPhysVar(Physics::Pres, enr_std));
    }

    // for surface 1, loop my col elements and add void enrichments to each elements member nodes
    const XFEM::Enrichment enr_void1(1, XFEM::Enrichment::typeVoid);
    for (int i=0; i<xfemdis->NumMyColElements(); ++i)
    {
        const DRT::Element* actele = xfemdis->lColElement(i);
        if (elementDomainIntCellMap.count(actele->Id()))
        {
            const int nen = actele->NumNode();
            const int* nodeidptrs = actele->NodeIds();
            for (int inen = 0; inen<nen; ++inen)
            {
                const int node_gid = nodeidptrs[inen];
                nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::Velx, enr_void1));
                nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::Vely, enr_void1));
                nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::Velz, enr_void1));
                nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::Pres, enr_void1));
                //              nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::LMPLambdax, enr_void1));
                //              nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::LMPLambday, enr_void1));
                //              nodalDofMap[node_gid].insert(XFEM::EnrPhysVar(Physics::LMPLambdaz, enr_void1));              
            };
        }
    };
    
    // create const sets from standard sets, so the sets cannot be changed by accident
    // could be removed later, if this is aperformance bottleneck
    for ( map<int, set <XFEM::EnrPhysVar> >::const_iterator oneset = nodalDofMap.begin(); oneset != nodalDofMap.end(); ++oneset )
    {
        nodalDofMapFinal.insert( make_pair(oneset->first, oneset->second));
    };
    return nodalDofMapFinal;
}

#endif  // #ifdef CCADISCRET
