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
#include "interface.H"
#include "../drt_lib/drt_discret.H"


// namespace for all xfem related stuff
namespace XFEM
{

using namespace Physics;

/*!
 \class EnrField

 \brief Corresponds to a specific enriched physical field variable

 type EnrField = (Field, Enrichment)

 */
class EnrField
{
public:

    explicit EnrField(
            const Field physvar,
            const Enrichment enr);
    EnrField(
            const EnrField& other);
    ~EnrField();

    string toString() const;

    inline Field getField() const
    {
        return field_;
    };

    inline Enrichment getEnrichment() const
    {
        return enr_;
    };
    
    inline bool operator <(const EnrField& rhs) const
    {
        if (getField() < rhs.getField())
            return true;
        else if (getField() > rhs.getField())
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
	Field field_;
    // specific enrichment for this physical field
    Enrichment enr_;

}; // class EnrField


/*!
 \class ElementDofManager

 \brief He knows everything about the XFEM enrichments for a particular interface setup for one given element

 */
class ElementDofManager
{
public:

	explicit ElementDofManager();
	
	explicit ElementDofManager(
			map<int, const set <XFEM::EnrField> >&  nodalDofMap);
	ElementDofManager(
            const ElementDofManager& other);
    ~ElementDofManager();

    std::string toString() const;

    inline const std::set<EnrField> getDofs(const int gid) const
    {
        //TODO: return proper default set given as second argument
    	map<int, const std::set<EnrField> >::const_iterator tmp = nodalDofMap_.find(gid);
    	if (tmp == nodalDofMap_.end()){
    	    //dserror("global node id not found!");
    	    cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
    	    std::set<EnrField> x;
    	    return x;
    	}
        return tmp->second;
    };
    
    /// function that returns node specific number of dofs or the given default number
    inline int getNumDof(
            const int gid,            /// global node id
            const int defaultnumdof   /// default number of degrees of freedom to be returned
            ) const
    {
    	map<int,int>::const_iterator tmp = nodalNumDofMap_.find(gid);
    	if (tmp == nodalNumDofMap_.end()){
    		//dserror("global node id not found!");
    		// "warning: global node id not found! Could be due to not properly initialized ElementDofManager"
    	    cout << "EleDM warning! ";
    		return defaultnumdof;
    	}
        return tmp->second;
    };
    
    inline int NumDof(const XFEM::Physics::Field field) const
    {
    	map<XFEM::Physics::Field, int>::const_iterator tmp = numParamsPerFieldMap_.find(field);
    	if (tmp == numParamsPerFieldMap_.end()){
    		cout << XFEM::Physics::physVarToString(field) << endl;
    		dserror("field not found!");
    	}
        return tmp->second;
    }
    
    inline vector<int> Dof(const XFEM::Physics::Field field) const
    {
    	map<XFEM::Physics::Field, vector<int> >::const_iterator tmp = paramsLocalEntries_.find(field);
    	if (tmp == paramsLocalEntries_.end()){
    		cout << XFEM::Physics::physVarToString(field) << endl;
    		dserror("field not found!");
    	}
        return tmp->second;
    }

private:
	//! unknowns for each global node id (gid)
	map<int, const set <XFEM::EnrField> > nodalDofMap_;

	//! number of unknowns for each global node id (gid)
	map<int, int> nodalNumDofMap_;
	
	//! number of parameters for each field
	map<XFEM::Physics::Field, int> numParamsPerFieldMap_;
	
	//! local (elemental) dofpos
	map<XFEM::Physics::Field, vector<int> > paramsLocalEntries_;
};



/*!
 \class DofManager

 \brief He knows everything about the XFEM enrichments for one particular interface setup

 */
class DofManager
{
public:

	explicit DofManager(const RCP<XFEM::InterfaceHandle> ih);
    ~DofManager();

    string toString() const;

    inline const std::set<EnrField> getDofs(const int gid) const
    {
    	map<int, const std::set<EnrField> >::const_iterator tmp = nodalDofMap_.find(gid);
    	if (tmp == nodalDofMap_.end())
    	{	
    		//dserror("global node id not found!");
    		cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
    		std::set<EnrField> x;
    		return x;
    	}
        return tmp->second;
    };
    
    inline std::pair<int, const set<XFEM::EnrField> > getDofsAsPair(const int gid) const
    {
    	map<int, const std::set<EnrField> >::const_iterator tmp = nodalDofMap_.find(gid);
    	if (tmp == nodalDofMap_.end()){
    		//dserror("global node id not found!");
    		cout << "warning: global node id not found! Could be due to not properly initialized ElementDofManager";
    		return std::pair<int, const set<XFEM::EnrField> >(tmp->first,tmp->second);
    	}
        return std::pair<int, const set<XFEM::EnrField> >(tmp->first,tmp->second);
    };

private:
	// pointer to the xfem discretization of which we hold the unknowns
	RefCountPtr<DRT::Discretization> xfemdis_;
	
	// unknowns for each global node id (gid) on this processor
	map<int, const set <XFEM::EnrField> > nodalDofMap_;
	
	// create a map that for each node (gid) on this processor has its unknowns
	const map<int, const set<EnrField> > createNodalDofMap(
			const RefCountPtr<DRT::Discretization> xfemdis,
	        const map<int, DomainIntCells >&  elementDomainIntCellMap) const;
};

typedef pair<int,EnrField> DofKey;

}; // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
