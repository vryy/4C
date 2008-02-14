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

#include "xfem.H"
#include "physics.H"
#include "enrichment.H"
#include "integrationcell.H"
#include "interface.H"
#include "../drt_lib/drt_discret.H"


// namespace for all xfem related stuff
namespace XFEM
{

/*!
 \brief Allows to distuingish between normal assembly and xfem assembly
 */
enum AssemblyType
{
    xfem_assembly,
    standard_assembly
};

/*!
 \class FieldEnr

 \brief Corresponds to a specific enriched physical field variable

 type FieldEnr = (Field, Enrichment)

 */
class FieldEnr
{
public:

    explicit FieldEnr();
    
    explicit FieldEnr(
            const PHYSICS::Field physvar,
            const Enrichment enr);
    FieldEnr(
            const FieldEnr& other);
    ~FieldEnr();

    string toString() const;

    PHYSICS::Field getField() const
    {
        return field_;
    };

    Enrichment getEnrichment() const
    {
        return enr_;
    };
    
    bool operator <(const FieldEnr& rhs) const
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
    
    bool operator ==(const FieldEnr& rhs) const
    {
        if (getField() == rhs.getField() and getEnrichment() == rhs.getEnrichment())
            return true;
        else
            return false;
    }
    
    bool operator !=(const FieldEnr& rhs) const
    {
        return (!(*this == rhs));
    }
    

private:
    // physical field variable to be enriched
    PHYSICS::Field field_;
    // specific enrichment for this physical field
    Enrichment enr_;

}; // class FieldEnr


/*!
 \class ElementDofManager

 \brief He knows everything about the XFEM enrichments for a particular interface setup for one given element

 */
class ElementDofManager
{
public:

	explicit ElementDofManager();
	
	explicit ElementDofManager(
	        const DRT::Element& ele,
			const map<int, const set<XFEM::FieldEnr> >&  nodalDofSet,
		    const set<XFEM::FieldEnr>& elementDofs,
		    const int numeleparam
		    );
	ElementDofManager(const ElementDofManager& old);
	
    ~ElementDofManager();

    std::string toString() const;

    /*!
     * \brief get a set with the enriched fields on this node
     * 
     * this function shall always fail, if the dofmanager not set properly
     * the number of unknowns has a default, but this one not!!!
     * 
     * \return set with enriched fields
     * 
     */
    std::set<FieldEnr> FieldEnrSetPerNode(const int gid) const
    {
    	map<int, const set<FieldEnr> >::const_iterator tmp = nodalDofSet_.find(gid);
    	if (tmp == nodalDofSet_.end()){
    	    cout << gid << endl;
    	    // let this always be here so we can recognize errors early
    	    dserror("FieldEnrSetPerNode: global node id not found!");
    	    return std::set<FieldEnr>();
    	}
        return tmp->second;
    };
    
    std::set<FieldEnr> FieldEnrSetPerElement() const
    {
        return elementDofs_;
    };
    
    //! function that returns node specific number of dofs or the given default number
    int NumDofPerNode(
            const int gid,            /// global node id
            const int defaultnumdof   /// default number of degrees of freedom to be returned
            ) const
    {
    	map<int,int>::const_iterator tmp = nodalNumDof_.find(gid);
    	if (tmp == nodalNumDof_.end()){
    		return defaultnumdof;
    	}
        return tmp->second;
    };
    
    /// function that returns element specific number of dofs or the given default number
    int NumDofPerElement() const
    {
        return elementDofs_.size()*numVirtualNodes_;
    };
    
    /// function that returns element specific number of dofs or the given default number
    int NumVirtualNodes() const
    {
        return numVirtualNodes_;
    };
    
    /// return number of dofs on this element for a given field. E.g., if we have 4 nodes with 
    /// standard x displacements at each node and 2 enriched x displacements, this function will return 6
    /// we do not distuingish between nodal and element dofs
    int NumDofPerField(const XFEM::PHYSICS::Field field) const
    {
    	map<XFEM::PHYSICS::Field, int>::const_iterator tmp = numParamsPerField_.find(field);
    	if (tmp == numParamsPerField_.end()){
            //cout << XFEM::PHYSICS::physVarToString(field) << endl;
    		return 0;
    	}
        return tmp->second;
    };
    
    /// return a list of local positions in a array of dofs
    std::vector<int> LocalDofPosPerField(const XFEM::PHYSICS::Field field) const
    {
    	std::map<XFEM::PHYSICS::Field, vector<int> >::const_iterator tmp = paramsLocalEntries_.find(field);
    	if (tmp == paramsLocalEntries_.end()){
    		//cout << XFEM::PHYSICS::physVarToString(field) << endl;
    		return vector<int>(0);
    	}
        return tmp->second;
    };
    
    bool operator !=(const ElementDofManager& rhs) const
    {
        // in principle, they can be only not equal, if the nodalDofSet is not equal
        // everything else is derived from that in the constructor
        bool not_equal = false;
        if (nodalDofSet_ != rhs.nodalDofSet_)
        {
            not_equal = true;
        }
        else
        {
            not_equal = false;
        }
        return not_equal;
    };
            
private:
	//! unknowns for each global node id (gid)
	std::map<int, const set<XFEM::FieldEnr> > nodalDofSet_;

	//! unknown fields for the element
    std::set<XFEM::FieldEnr> elementDofs_;
    
    /*! multiply the size of elementDofs_ with this number and you get the number of overall
     * element unknowns.
     */
    int numVirtualNodes_;
	
	//! number of unknowns for each global node id (gid)
	std::map<int, int> nodalNumDof_;
	
	//! number of parameters for each field
	std::map<XFEM::PHYSICS::Field, int> numParamsPerField_;
	
	//! local (elemental) dofpos
	std::map<XFEM::PHYSICS::Field, std::vector<int> > paramsLocalEntries_;
};



/*!
 \class DofManager

 \brief He knows everything about the XFEM enrichments for one particular interface setup

 */
class DofManager
{
public:
    explicit DofManager();
	explicit DofManager(const RCP<XFEM::InterfaceHandle> ih);
	explicit DofManager(const XFEM::DofManager& dofman);
    ~DofManager();

    string toString() const;

    const std::set<FieldEnr> getNodeDofSet(const int gid) const
    {
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
        if (entry == nodalDofSet_.end())
            if (entry == nodalDofSet_.end())
            {
                return std::set<XFEM::FieldEnr>();
            }
            return entry->second;
    };
    
    const std::set<XFEM::FieldEnr> getElementDofs(const int gid) const
    {
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = elementalDofs_.find(gid);
        if (entry == elementalDofs_.end())
        {
            return std::set<XFEM::FieldEnr>();
        }
        return entry->second;
    };
    
    const XFEM::ElementDofManager constructElementDofManager(
            const DRT::Element& ele,
            const int numeleparam) const;
    
    void checkForConsistency(
            const DRT::Element& ele,
            const XFEM::ElementDofManager& eledofman
            ) const;
    
private:
	// pointer to the xfem discretization of which we hold the unknowns
	RCP<DRT::Discretization> xfemdis_;
	
	// unknowns for each global node id (gid) on this processor
	// if it does not contain the key, the number of unknowns is supposed to be zero
	// meaning, the the dofmanager will not crash but give back an empty set of unknowns
	// see also getNodeDofSet()
	map<int, const set<XFEM::FieldEnr> > nodalDofSet_;
	
	// unknowns for each global element id (gid) on this processor
	// same as above, see getElementDofs()
	map<int, const set<XFEM::FieldEnr> > elementalDofs_;
	
};


XFEM::AssemblyType CheckForStandardEnrichmentsOnly(
        const ElementDofManager&   eleDofManager_,
        const int                  numnode,
        const int*                 nodeids
        );

} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
