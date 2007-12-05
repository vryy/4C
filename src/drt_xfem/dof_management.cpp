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
		map<int, const set <XFEM::FieldEnr> >& nodalDofMap) :
			nodalDofMap_(nodalDofMap)
{
	// count number of dos for each node
	map<int, const set<XFEM::FieldEnr> >::const_iterator tmp;
	for (tmp = nodalDofMap.begin(); tmp != nodalDofMap.end(); ++tmp) {
		const int gid = tmp->first;
		const set<XFEM::FieldEnr> enrfieldset = tmp->second;
		const int numdof = enrfieldset.size();
		dsassert(numdof > 0, "sollte jetzt noch nicht sein");
		nodalNumDofMap_[gid] = numdof;
	}
	
	// set number of parameters per field to zero
	set<XFEM::FieldEnr>::const_iterator enrfield;
	for (tmp = nodalDofMap.begin(); tmp != nodalDofMap.end(); ++tmp) {
		const set<XFEM::FieldEnr> enrfieldset = tmp->second;
		for (enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
			const XFEM::PHYSICS::Field field = enrfield->getField();
			numParamsPerFieldMap_[field] = 0;
			paramsLocalEntries_[field] = vector<int>();
		}
	}
	
	// count number of parameters per field
	int counter = 0;
	for (tmp = nodalDofMap.begin(); tmp != nodalDofMap.end(); ++tmp) {
		const set<XFEM::FieldEnr> enrfieldset = tmp->second;
		
		for (enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
			const XFEM::PHYSICS::Field field = enrfield->getField();
			numParamsPerFieldMap_[field] += 1;
			paramsLocalEntries_[field].push_back(counter);
			counter++;
		}
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
	for (tmp = nodalDofMap_.begin(); tmp != nodalDofMap_.end(); ++tmp)
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
		const set <XFEM::FieldEnr> actset = nodalDofMap_.find(gid)->second;
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
    
    map<int, const set <XFEM::FieldEnr> > nodaldofset; 
    for (int inode = 0; inode < numnode; ++inode) {
        const int gid = nodegids[inode];
        nodaldofset.insert(this->getDofsAsPair(gid));
    }

    // create a local dofmanager
    XFEM::ElementDofManager eleDofManager = XFEM::ElementDofManager(nodaldofset);
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
const map<int, const set <XFEM::FieldEnr> > XFEM::DofManager::createNodalDofMap(
        const RCP<DRT::Discretization>        xfemdis,
        const map<int, XFEM::DomainIntCells >&  elementDomainIntCellMap) const
{
	// the final map
	map<int, const set <XFEM::FieldEnr> >  nodalDofMapFinal;

	// temporary assembly
    map<int, set <XFEM::FieldEnr> >  nodalDofMap;
    
    // standard enrichment used for all nodes (for now -> we can remove them from holes in the fluid)
    const XFEM::Enrichment enr_std(0, XFEM::Enrichment::typeStandard);
    // loop my row nodes and add standard degrees of freedom
    for (int i=0; i<xfemdis->NumMyColNodes(); ++i)
    {
        const DRT::Node* actnode = xfemdis->lColNode(i);
        const int gid = actnode->Id();
        nodalDofMap[gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_std));
        nodalDofMap[gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_std));
        nodalDofMap[gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_std));
        nodalDofMap[gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_std));
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
                //nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velx, enr_void1));
                //nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::Vely, enr_void1));
                //nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::Velz, enr_void1));
                //nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::Pres, enr_void1));
                //              nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::LMPLambdax, enr_void1));
                //              nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::LMPLambday, enr_void1));
                //              nodalDofMap[node_gid].insert(XFEM::FieldEnr(PHYSICS::LMPLambdaz, enr_void1));              
            };
        }
    };
    
    // create const sets from standard sets, so the sets cannot be changed by accident
    // could be removed later, if this is aperformance bottleneck
    for ( map<int, set <XFEM::FieldEnr> >::const_iterator oneset = nodalDofMap.begin(); oneset != nodalDofMap.end(); ++oneset )
    {
        nodalDofMapFinal.insert( make_pair(oneset->first, oneset->second));
    };
    return nodalDofMapFinal;
}



bool inCircleCylinder(
        const blitz::Array<double,1>& pos,
        const blitz::Array<double,1>& center,
        const double cylinder_radius
        )
{
    blitz::Range _  = blitz::Range::all();
    const blitz::Array<double,1> origincircle(pos(_) - center(_));
    
    const double circle_radius = sqrt(origincircle(0)*origincircle(0) + origincircle(1)*origincircle(1));
    
    bool in_circle = false;
    if (circle_radius <= cylinder_radius){
        in_circle = true;
    } else {
        in_circle = false;
    }
    return in_circle;
}

/*----------------------------------------------------------------------*
 |  get enrichment value                                        ag 11/07|
 *----------------------------------------------------------------------*/
double XFEM::enrValue(
        const XFEM::Enrichment enr,
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& nodalpos
        )
{
    // return value
    double enrval = 1.0;
    
    switch (enr.Type()){
    case XFEM::Enrichment::typeStandard:
    {
        enrval = 1.0;
        break;
    }
//    case XFEM::Enrichment::typeVoid:
//    {
//        // TODO: generalize
//        blitz::Array<double,1> center(3);
//        center(0) = 0.6; center(1) = 0.5; center(2) = 0.0;
//        const double cylinder_radius = 0.2;
//       
//        double actpos_enr_val = 0.0;
//        if (inCircleCylinder(actpos, center, cylinder_radius)) {
//            actpos_enr_val = 0.0;
//        } else {
//            actpos_enr_val = 1.0;
//        }
//        
//        double nodepos_enr_val = 0.0;
//        if (inCircleCylinder(actpos, center, cylinder_radius)) {
//            nodepos_enr_val = 0.0;
//        } else {
//            nodepos_enr_val = 1.0;
//        }
//        
//        enrval = actpos_enr_val - nodepos_enr_val;
//        
//        break;
//    }
    case XFEM::Enrichment::typeVoid:
    {
        // TODO: generalize
        double actpos_enr_val = 0.0;
        if (actpos(0) > 1.525) {
            actpos_enr_val = 0.0;
        } else {
            actpos_enr_val = 1.0;
        }
        
        double nodepos_enr_val = 0.0;
        if (nodalpos(0) > 1.525) {
            nodepos_enr_val = 0.0;
        } else {
            nodepos_enr_val = 1.0;
        }
        
        enrval = actpos_enr_val - nodepos_enr_val;
        enrval = 1.0;
        dserror("not yet");
        
        break;
    }
    default:
        dserror("unsupported enrichment!");
    }
    return enrval;
}



#endif  // #ifdef CCADISCRET
