/*!
\file dof_management.H

\brief provides the dofmanager classes

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
#include "field_enriched.H"
#include "dofkey.H"
#include "integrationcell.H"
#include "interface.H"
#include "../drt_lib/drt_discret.H"


namespace XFEM
{

/*!
 \brief Allows to distuingish between normal assembly and xfem assembly
 */
enum AssemblyType
{
    standard_assembly,
    xfem_assembly
};


/*!
 \brief He knows everything about the XFEM enrichments for a particular interface setup for one given element

 */
class ElementDofManager
{
public:
    //! default constructor
    explicit ElementDofManager();
    
    //! standard constructor
    explicit ElementDofManager(
            const DRT::Element& ele,
            const map<int, const set<XFEM::FieldEnr> >&  nodalDofSet,
            const map<int, const set<XFEM::FieldEnr> >&  elementnodalDofSet,
            const int numvirtualnodes
            );
    //! copy constructor
    ElementDofManager(const ElementDofManager& old);
    //! destructor
    ~ElementDofManager();
    //! return a string representation of this class
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
    const std::set<FieldEnr>& FieldEnrSetPerNode(
            const int  gid    ///< unique global node id
            ) const
    {
        map<int, const set<FieldEnr> >::const_iterator tmp = nodalDofSet_.find(gid);
        if (tmp == nodalDofSet_.end()){
            cout << gid << endl;
            // let this always be here so we can recognize errors early
            dserror("FieldEnrSetPerNode: global node id not found!");
            //return std::set<FieldEnr>();
        }
        return tmp->second;
    };
    
    /*!
     * \brief get a set with the enriched fields on this element
     * 
     * \return set with enriched fields
     * 
     */
    const std::set<FieldEnr>& FieldEnrSetPerVirtualElementNode(
        const int  lid    ///< local node id
        ) const
    {
      map<int, const set<FieldEnr> >::const_iterator tmp = elementnodalDofSet_.find(lid);
      if (tmp == elementnodalDofSet_.end()){
          cout << lid << endl;
          // let this always be here so we can recognize errors early
          dserror("FieldEnrSetPerVirtualElementNode: local node id not found!");
          //return std::set<FieldEnr>();
      }
      return tmp->second;
    };
    
    //! function that returns node specific number of dofs or the given default number
    int NumDofPerNode(
            const int gid,            ///< unique global node id
            const int defaultnumdof   ///< default number of degrees of freedom to be returned
            ) const
    {
        map<int,int>::const_iterator tmp = nodalNumDof_.find(gid);
        if (tmp == nodalNumDof_.end())
        {
            return defaultnumdof;
        }
        return tmp->second;
    };
    
    //! function that returns element specific number of dofs or the given default number
    int NumDofPerElement() const
    {
        int numdofele = 0;
        const int numvirtualnodes = elementnodalDofSet_.size();
        if (numvirtualnodes > 0)
        {
          numdofele = numvirtualnodes * elementnodalDofSet_.find(0)->second.size();
        }
        return numdofele;
    };
    
    //! function that returns element specific number of dofs or the given default number
    int NumVirtualNodes() const
    {
        return elementnodalDofSet_.size();
    };
    
    /*! return number of dofs on this element for a given field. E.g., if we have 4 nodes with
     *  standard x displacements at each node and 2 enriched x displacements, this function will return 6
     *  we do not distuingish between nodal and element dofs
     */
    int NumDofPerField(
            const XFEM::PHYSICS::Field  field  ///< field for which we seek the number of DOFs
            ) const
    {
        map<XFEM::PHYSICS::Field, int>::const_iterator tmp = numParamsPerField_.find(field);
        if (tmp == numParamsPerField_.end()){
            //cout << XFEM::PHYSICS::physVarToString(field) << endl;
            return 0;
        }
        return tmp->second;
    };
    
    //! return reference to list of local positions in a array of dofs
    const std::vector<int>& LocalDofPosPerField(const XFEM::PHYSICS::Field field) const
    {
        return paramsLocalEntries_.find(field)->second;
    };
    
    //! return unique enrichments for this proc
    const std::set<XFEM::Enrichment>& getUniqueEnrichments() const
    {
        return unique_enrichments_;
    }
    
    //! not equal operator (for use in STL iterators)
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

    //! unknowns for each virtual node id (lid)
    std::map<int, const set<XFEM::FieldEnr> > elementnodalDofSet_;
    
    //! number of unknowns for each global node id (gid)
    std::map<int, int> nodalNumDof_;
    
    //! number of parameters for each field
    std::map<XFEM::PHYSICS::Field, int> numParamsPerField_;
    
    //! local (elemental) dofpos
    std::map<XFEM::PHYSICS::Field, std::vector<int> > paramsLocalEntries_;
    
    //! set of all enrichments that are available on this proc.
    std::set<XFEM::Enrichment> unique_enrichments_;
};



/*!
 \brief He knows everything about the XFEM enrichments for one particular interface setup

 */
class DofManager
{
public:
    //! default constructor
    explicit DofManager();
    //! constructor
    explicit DofManager(const RCP<XFEM::InterfaceHandle> ih);
    //! copy constructor
    explicit DofManager(const XFEM::DofManager& dofman);
    //! destructor
    ~DofManager();
    //! return string representation of the dofmanager
    std::string toString() const;
    //! print a Gmsh representation of this class
    void toGmsh(const Teuchos::RCP<XFEM::InterfaceHandle> ih, const int step) const;

    const std::set<FieldEnr> getNodeDofSet(const int gid) const
    {
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
        if (entry == nodalDofSet_.end())
        {
            return std::set<XFEM::FieldEnr>();
        }
        return entry->second;
    };
    
    const std::set<XFEM::FieldEnr> getElementDofSet(const int gid) const
    {
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = elementalDofs_.find(gid);
        if (entry == elementalDofs_.end())
        {
            return std::set<XFEM::FieldEnr>();
        }
        return entry->second;
    };
    
    XFEM::ElementDofManager constructElementDofManager(
            const DRT::Element& ele,
            const int numeleparam) const;
    
    void checkForConsistency(
            const DRT::Element& ele,
            const XFEM::ElementDofManager& eledofman
            ) const;
    
    void fillDofDistributionMap(
            DofPosMap&   NodalDofDistributionMap
            ) const;
    
private:
    //! pointer to the xfem discretization of which we hold the unknowns
    Teuchos::RCP<DRT::Discretization> xfemdis_;
    
    /*! unknowns for each global node id (gid) on this processor
     * if it does not contain the key, the number of unknowns is supposed to be zero
     * meaning, the the dofmanager will not crash but give back an empty set of unknowns
     * see also getNodeDofSet()
     */
    std::map<int, const set<XFEM::FieldEnr> > nodalDofSet_;
    
    /*! unknowns for each global element id (gid) on this processor
     * same as above, see getElementDofSet()
     */
    std::map<int, const set<XFEM::FieldEnr> > elementalDofs_;
    
    //! set of all enrichments that are available on this proc.
    std::set<XFEM::Enrichment> unique_enrichments_;
    
};

//! return assembly type depending on the DOF distribution
XFEM::AssemblyType CheckForStandardEnrichmentsOnly(
        const ElementDofManager&   eleDofManager_,
        const int                  numnode,
        const int*                 nodeids
        );

} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
