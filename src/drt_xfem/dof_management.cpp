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
#include "../io/gmsh.H"



using namespace std;

/*----------------------------------------------------------------------*
 |  default ctor                                                ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::FieldEnr::FieldEnr() :
            field_(XFEM::PHYSICS::undefinedField),
            enr_(Enrichment())
{
    dserror("FieldEnr() -> please don't call me!");
    return;
}

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
    s << "Enriched Field: " << PHYSICS::physVarToString(field_) << ", Enrichment: " << enr_.toString();
    return s.str();
}





/*----------------------------------------------------------------------*
 |  default ctor                                                ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager() :
    nodalDofSet_(),
    elementDofs_(),
    numVirtualNodes_(0)
{
    return;
}

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
        const DRT::Element& ele,
        const map<int, const set<XFEM::FieldEnr> >& nodalDofSet,
        const set<XFEM::FieldEnr>& elementDofs,
        const int numeleparam
        ) :
            nodalDofSet_(nodalDofSet),
            elementDofs_(elementDofs),
            numVirtualNodes_(numeleparam)
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
    
    
    unique_enrichments_.clear();
    // count number of parameters per field
    // define local position of unknown by looping first over nodes and then over its unknowns!
    int counter = 0;
    const DRT::Node*const* nodes = ele.Nodes();
    for (int inode=0; inode<ele.NumNode(); ++inode)
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
            unique_enrichments_.insert(enrfield->getEnrichment());
            counter++;
        }
    }
    // loop now over element dofs
    // for that we have to loop over the "virtual element nodes" ansatz function
    for (int i = 0; i < numeleparam; ++i)
    {
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield = elementDofs.begin(); enrfield != elementDofs.end(); ++enrfield)
        {
            const XFEM::PHYSICS::Field field = enrfield->getField();
            numParamsPerField_[field] += 1;
            paramsLocalEntries_[field].push_back(counter);
            unique_enrichments_.insert(enrfield->getEnrichment());
            counter++;
        }
    }
    
    return;
}

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
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
    
    unique_enrichments_.clear();
    for (map<int, const set<XFEM::FieldEnr> >::const_iterator fieldenriter=nodalDofSet_.begin();
            fieldenriter!=nodalDofSet_.end(); ++fieldenriter)
    {
        const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            unique_enrichments_.insert(enrfield->getEnrichment());
        }
    }
    for (map<int, const set<XFEM::FieldEnr> >::const_iterator fieldenriter=elementalDofs_.begin();
            fieldenriter!=elementalDofs_.end(); ++fieldenriter)
    {
        const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            unique_enrichments_.insert(enrfield->getEnrichment());
        }
    }
    cout << " Enrichments available:" << endl;
    for (std::set<XFEM::Enrichment>::const_iterator enr =
        unique_enrichments_.begin(); enr != unique_enrichments_.end(); ++enr)
    {
        cout << "  - " << enr->toString() << endl;
    }
    
    
    
    std::ofstream f_system("numdof_coupled_system.pos");
    //f_system << IO::GMSH::disToString("Fluid", 0.0, ih->xfemdis(), ih->elementalDomainIntCells());
    f_system << IO::GMSH::disToString("Solid", 1.0, ih->cutterdis());
    {
        // draw elements with associated gid
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "Element->Id() \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
        {
            DRT::Element* actele = ih->xfemdis()->lColElement(i);
            gmshfilecontent << IO::GMSH::elementToString(double(actele->Id()), actele);
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << " Stress unknowns in element \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
        {
            DRT::Element* actele = ih->xfemdis()->lColElement(i);
            const int ele_gid = actele->Id();
            double val = 0.0;
            std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = elementalDofs_.find(ele_gid);
            
            if (blub != elementalDofs_.end())
            {
                const set<XFEM::FieldEnr> schnapp = blub->second;
                val = schnapp.size();
                gmshfilecontent << IO::GMSH::elementToString(val, actele);
            }
            
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "NumDof per node \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
        {
            //DRT::Element* actele = ih->xfemdis()->lColElement(i);
            const DRT::Node* actnode = ih->xfemdis()->lColNode(i);
            const blitz::Array<double,1> pos(toBlitzArray(actnode->X()));
            const int node_gid = actnode->Id();
            
            double val = 0.0;
            std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
            
            if (blub != nodalDofSet_.end())
            {
                const set<XFEM::FieldEnr> schnapp = blub->second;
                val = schnapp.size();
            

            gmshfilecontent << "SP(";
            gmshfilecontent << scientific << pos(0) << ",";
            gmshfilecontent << scientific << pos(1) << ",";
            gmshfilecontent << scientific << pos(2);
            gmshfilecontent << "){";
            gmshfilecontent << val << "};" << endl;
            }
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    
    {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "NumDof Jump enriched nodes \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
        {
            const DRT::Node* actnode = ih->xfemdis()->lColNode(i);
            const blitz::Array<double,1> pos(toBlitzArray(actnode->X()));
            const int node_gid = actnode->Id();
            
            double val = 0.0;
            std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
            if (blub != nodalDofSet_.end())
            {
                const std::set<XFEM::FieldEnr> fields = blub->second;
                for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
                {
                    if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeJump)
                    {
                        val = val+1.0;
                    }
                }
                if (val > 0.5)
                {
                    gmshfilecontent << "SP(";
                    gmshfilecontent << scientific << pos(0) << ",";
                    gmshfilecontent << scientific << pos(1) << ",";
                    gmshfilecontent << scientific << pos(2);
                    gmshfilecontent << "){";
                    gmshfilecontent << val << "};" << endl;
                }
            }
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    
    {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "NumDof" << " standard enriched nodes \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
        {
            const DRT::Node* actnode = ih->xfemdis()->lColNode(i);
            const blitz::Array<double,1> pos(toBlitzArray(actnode->X()));
            const int node_gid = actnode->Id();
            
            double val = 0.0;
            std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
            if (blub != nodalDofSet_.end())
            {
                const std::set<XFEM::FieldEnr> fields = blub->second;
                for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
                {
                    if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeStandard)
                    {
                        val = val+1.0;
                    }
                }
                if (val > 0.5)
                {
                    gmshfilecontent << "SP(";
                    gmshfilecontent << scientific << pos(0) << ",";
                    gmshfilecontent << scientific << pos(1) << ",";
                    gmshfilecontent << scientific << pos(2);
                    gmshfilecontent << "){";
                    gmshfilecontent << val << "};" << endl;
                }
            }
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    
    {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << "NumDof" << " Void enriched nodes \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
        {
            const DRT::Node* actnode = ih->xfemdis()->lColNode(i);
            const blitz::Array<double,1> pos(toBlitzArray(actnode->X()));
            const int node_gid = actnode->Id();
            
            double val = 0.0;
            std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
            if (blub != nodalDofSet_.end())
            {
                const std::set<XFEM::FieldEnr> fields = blub->second;
                for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
                {
                    if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeVoid)
                    {
                        val = val+1.0;
                    }
                }
                if (val > 0.5)
                {
                    gmshfilecontent << "SP(";
                    gmshfilecontent << scientific << pos(0) << ",";
                    gmshfilecontent << scientific << pos(1) << ",";
                    gmshfilecontent << scientific << pos(2);
                    gmshfilecontent << "){";
                    gmshfilecontent << val << "};" << endl;
                }
            }
        };
        gmshfilecontent << "};" << endl;
        f_system << gmshfilecontent.str();
    }
    
    
    //f_system << IO::GMSH::getConfigString(2);
    f_system.close();
    
    
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
const XFEM::ElementDofManager XFEM::DofManager::constructElementDofManager(
        const DRT::Element&  ele,
        const int      numeleparam
        ) const
{
    // create a list with number of dofs per local node 
    const int numnode = ele.NumNode();
    const int* nodegids = ele.NodeIds();
    
    // nodal dofs for ele
    std::map<int, const set <XFEM::FieldEnr> > nodaldofset;
    for (int inode = 0; inode < numnode; ++inode) {
        const int gid = nodegids[inode];
        nodaldofset.insert(make_pair(gid,this->getNodeDofSet(gid)));
    }

    // element dofs for ele
    const std::set<XFEM::FieldEnr> elementdofs = this->getElementDofs(ele.Id());
    
    // create a local dofmanager
    XFEM::ElementDofManager eleDofManager(ele, nodaldofset, elementdofs, numeleparam);

    return eleDofManager;
}

/*----------------------------------------------------------------------*
 |  sanity check                                                ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::DofManager::checkForConsistency(
        const DRT::Element& ele,
        const XFEM::ElementDofManager& stored_eledofman
        ) const
{
    // create local copy of current information about dofs
    XFEM::ElementDofManager current_eledofman = this->constructElementDofManager(ele, stored_eledofman.NumVirtualNodes());
    
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
