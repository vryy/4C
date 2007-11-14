/*!
\file dof_management.H

\brief provides a class that represents an enriched physical scalar field

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef DOF_MANAGEMENT_H
#define DOF_MANAGEMENT_H

#include "physics.H"
#include "enrichment.H"
#include "integrationcell.H"
#include "../drt_lib/drt_discret.H"


// namespace for all xfem related stuff
namespace XFEM
{

using namespace Physics;

/*!
 \class EnrPhysVar

 \brief Corresponds to a specific enriched physical field variable

 type EnrPhysVar = (PhysVar, Enrichment)

 */
class EnrPhysVar
{
public:

    EnrPhysVar(
            const PhysVar physvar,
            const Enrichment enr);
    EnrPhysVar(
            const EnrPhysVar& other);
    ~EnrPhysVar();

    string toString() const;

    inline PhysVar getPhysVar() const
    {
        return physvar_;
    };

    inline Enrichment getEnrichment() const
    {
        return enr_;
    };
    
    inline bool operator <(const EnrPhysVar& rhs) const
    {
        if (getPhysVar() < rhs.getPhysVar())
            return true;
        else if (getPhysVar() > rhs.getPhysVar())
            return false;
        else
        {
            if (getEnrichment() < rhs.getEnrichment())
                return true;
            else
                return false;
        }
    } 
    

private:
    // physical field variable to be enriched
    PhysVar physvar_;
    // specific enrichment for this physical field
    Enrichment enr_;

}; // class EnrPhysVar


/*!
 \class ElementDofManager

 \brief He knows everything about the XFEM enrichments for a particular interface setup for one given element

 */
class ElementDofManager
{
public:

	ElementDofManager();
	
	ElementDofManager(
			map<int, const set <XFEM::EnrPhysVar> >&  nodalDofMap);
	ElementDofManager(
            const ElementDofManager& other);
    ~ElementDofManager();

    std::string toString() const;

    inline const std::set<EnrPhysVar> getDofs(const int gid) const
    {
    	map<int, const std::set<EnrPhysVar> >::const_iterator tmp = nodalDofMap_.find(gid);
    	if (tmp == nodalDofMap_.end()){
    	    //dserror("global node id not found!");
    	    cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
    	    std::set<EnrPhysVar> x;
    	    return x;
    	}
        return tmp->second;
    };
    
    inline int getNumDof(const int gid) const
    {
    	map<int,int>::const_iterator tmp = nodalNumDofMap_.find(gid);
    	if (tmp == nodalNumDofMap_.end()){
    		//dserror("global node id not found!");
    		cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
    		return 0;
    	}
        return tmp->second;
    };

private:
	// unknowns for each global node id (gid)
	map<int, const set <XFEM::EnrPhysVar> > nodalDofMap_;

	// number of unknowns for each global node id (gid)
	map<int, int> nodalNumDofMap_;
};



/*!
 \class DofManager

 \brief He knows everything about the XFEM enrichments for one particular interface setup

 */
class DofManager
{
public:

	DofManager(
			RCP<DRT::Discretization> xfemdis,
	        const map<int, DomainIntCells >&  elementalDomainIntCells);
	DofManager(
            const DofManager& other);
    ~DofManager();

    string toString() const;

    inline const std::set<EnrPhysVar> getDofs(const int gid) const
    {
    	map<int, const std::set<EnrPhysVar> >::const_iterator tmp = nodalDofMap_.find(gid);
    	if (tmp == nodalDofMap_.end())
    	{	
    		//dserror("global node id not found!");
    		cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
    		std::set<EnrPhysVar> x;
    		return x;
    	}
        return tmp->second;
    };
    
    inline std::pair<int, const set<XFEM::EnrPhysVar> > getDofsAsPair(const int gid) const
    {
    	map<int, const std::set<EnrPhysVar> >::const_iterator tmp = nodalDofMap_.find(gid);
    	if (tmp == nodalDofMap_.end()){
    		//dserror("global node id not found!");
    		cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
//    		std::set<EnrPhysVar> s;
//    		std::pair<int, const set<XFEM::EnrPhysVar> > x(-1,s);
    		return std::pair<int, const set<XFEM::EnrPhysVar> >(tmp->first,tmp->second);
    	}
        return std::pair<int, const set<XFEM::EnrPhysVar> >(tmp->first,tmp->second);
    };

private:
	// pointer to the xfem discretization of which we hold the unknowns
	RefCountPtr<DRT::Discretization> xfemdis_;
	
	// unknowns for each global node id (gid) on this processor
	map<int, const set <XFEM::EnrPhysVar> > nodalDofMap_;
	
	// create a map that for each node (gid) on this processor has its unknowns
	map<int, const set<EnrPhysVar> > createNodalDofMap(
			const RefCountPtr<DRT::Discretization> xfemdis,
	        const map<int, DomainIntCells >&  elementDomainIntCellMap) const;
};



}; // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
