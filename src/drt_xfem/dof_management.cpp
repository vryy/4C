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



map<int, set <XFEM::EnrPhysVar> > XFEM::createNodalDofMap(
        RefCountPtr<DRT::Discretization>            xfemdis,
        const map<int, vector <XFEM::Integrationcell> >&  elementIntCellMap)
{
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
        if (elementIntCellMap.count(actele->Id()))
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
    return nodalDofMap;
}

#endif  // #ifdef CCADISCRET
