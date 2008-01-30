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
#include "xdofmapcreation.H"



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
