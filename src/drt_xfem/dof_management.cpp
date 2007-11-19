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

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::EnrField::EnrField(
        const Physics::Field field,
        const Enrichment enr) :
        	field_(field), enr_(enr)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                   ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::EnrField::EnrField(
        const EnrField& other) :
        	field_(other.field_), 
            enr_(other.enr_)
{
    assert(&other != this);
    return;
}

/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::EnrField::~EnrField()
{
    return;
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
string XFEM::EnrField::toString() const
{
    stringstream s;
    s << "Enriched Field: " << Physics::physVarToString(this->field_) << ", Enrichment: " << enr_.toString();
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
		map<int, const set <XFEM::EnrField> >& nodalDofMap) :
			nodalDofMap_(nodalDofMap)
{
	set<XFEM::EnrField>::const_iterator enrfield;
	map<int, const set<XFEM::EnrField> >::const_iterator tmp;
	
	for (tmp = nodalDofMap.begin(); tmp != nodalDofMap.end(); ++tmp) {
		const int gid = tmp->first;
		const set<XFEM::EnrField> enrfieldset = tmp->second;
		const int numdof = enrfieldset.size();
		nodalNumDofMap_[gid] = numdof;
	}
	
	// set number of parameters per field to zero
	for (tmp = nodalDofMap.begin(); tmp != nodalDofMap.end(); ++tmp) {
		const set<XFEM::EnrField> enrfieldset = tmp->second;
		for (enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
			const XFEM::Physics::Field field = enrfield->getField();
			numParamsPerFieldMap_[field] = 0;
			paramsLocalEntries_[field] = vector<int>();
		}
	}
	
	// count number of parameters per field
	int counter = 0;
	for (tmp = nodalDofMap.begin(); tmp != nodalDofMap.end(); ++tmp) {
		const set<XFEM::EnrField> enrfieldset = tmp->second;
		
		for (enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
			const XFEM::Physics::Field field = enrfield->getField();
			numParamsPerFieldMap_[field] += 1;
			paramsLocalEntries_[field].push_back(counter);
			counter++;
		}
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                   ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
        const ElementDofManager& other) :
        	nodalDofMap_(other.nodalDofMap_), 
        	nodalNumDofMap_(other.nodalNumDofMap_)
{
    assert(&other != this);
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
	map<int, const set<XFEM::EnrField> >::const_iterator tmp;
	for (tmp = nodalDofMap_.begin(); tmp != nodalDofMap_.end(); ++tmp)
	{
		const int gid = tmp->first;
		const set <XFEM::EnrField> actset = tmp->second;
		for ( set<XFEM::EnrField>::const_iterator var = actset.begin(); var != actset.end(); ++var )
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
	nodalDofMap_ = XFEM::DofManager::createNodalDofMap(ih->xfemdis(), ih->elementalDomainIntCells());
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
		const set <XFEM::EnrField> actset = nodalDofMap_.find(gid)->second;
		for ( set<XFEM::EnrField>::const_iterator var = actset.begin(); var != actset.end(); ++var )
		{
			s << "Node: " << gid << ", " << var->toString() << endl;
		};
	};
	return s.str();
}

/*----------------------------------------------------------------------*
 |  construct dofmap                                            ag 11/07|
 *----------------------------------------------------------------------*/
const map<int, const set <XFEM::EnrField> > XFEM::DofManager::createNodalDofMap(
        const RCP<DRT::Discretization>        xfemdis,
        const map<int, XFEM::DomainIntCells >&  elementDomainIntCellMap) const
{
	// the final map
	map<int, const set <XFEM::EnrField> >  nodalDofMapFinal;

	// temporary assembly
    map<int, set <XFEM::EnrField> >  nodalDofMap;
    
    // standard enrichment used for all nodes (for now -> we can remove them from holes in the fluid)
    const XFEM::Enrichment enr_std(0, XFEM::Enrichment::typeStandard);
    // loop my row nodes and add standard degrees of freedom
    for (int i=0; i<xfemdis->NumMyColNodes(); ++i)
    {
        const DRT::Node* actnode = xfemdis->lColNode(i);
        const int gid = actnode->Id();
        nodalDofMap[gid].insert(XFEM::EnrField(Physics::Velx, enr_std));
        nodalDofMap[gid].insert(XFEM::EnrField(Physics::Vely, enr_std));
        nodalDofMap[gid].insert(XFEM::EnrField(Physics::Velz, enr_std));
        nodalDofMap[gid].insert(XFEM::EnrField(Physics::Pres, enr_std));
    }

    // for surface 1, loop my col elements and add void enrichments to each elements member nodes
//    const XFEM::Enrichment enr_void1(1, XFEM::Enrichment::typeVoid);
//    for (int i=0; i<xfemdis->NumMyColElements(); ++i)
//    {
//        const DRT::Element* actele = xfemdis->lColElement(i);
//        if (elementDomainIntCellMap.count(actele->Id()))
//        {
//            const int nen = actele->NumNode();
//            const int* nodeidptrs = actele->NodeIds();
//            for (int inen = 0; inen<nen; ++inen)
//            {
//                const int node_gid = nodeidptrs[inen];
//                nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::Velx, enr_void1));
//                nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::Vely, enr_void1));
//                nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::Velz, enr_void1));
//                nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::Pres, enr_void1));
//                //              nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::LMPLambdax, enr_void1));
//                //              nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::LMPLambday, enr_void1));
//                //              nodalDofMap[node_gid].insert(XFEM::EnrField(Physics::LMPLambdaz, enr_void1));              
//            };
//        }
//    };
    
    // create const sets from standard sets, so the sets cannot be changed by accident
    // could be removed later, if this is aperformance bottleneck
    for ( map<int, set <XFEM::EnrField> >::const_iterator oneset = nodalDofMap.begin(); oneset != nodalDofMap.end(); ++oneset )
    {
        nodalDofMapFinal.insert( make_pair(oneset->first, oneset->second));
    };
    return nodalDofMapFinal;
}

#endif  // #ifdef CCADISCRET
