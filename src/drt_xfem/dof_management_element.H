/*!
\file dof_management_element.H

\brief provides the element dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET
#ifndef DOF_MANAGEMENT_ELEMENT_H
#define DOF_MANAGEMENT_ELEMENT_H


#include "dof_management.H"

namespace XFEM
{
  
  //! Allows to distinguish between normal assembly and xfem assembly
  enum AssemblyType
  {
    standard_assembly,
    xfem_assembly
  };
  
  
  /*!
   * \brief He knows everything about the XFEM enrichments for a particular interface setup for one given element
   * 
   */
  class ElementDofManager
  {
    public:
      
      //! default constructor
      explicit ElementDofManager();
      
      //! standard constructor
      explicit ElementDofManager(
          const DRT::Element& ele,
          const map<int, const std::set<XFEM::FieldEnr> >&  nodalDofSet,
          const std::set<XFEM::FieldEnr>& enrfieldset,
          const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
      );
      
      //! for a given element, collect all information about Dofs for this element and return it
      explicit ElementDofManager(
          const DRT::Element&                                                     ele,            ///< construct EleDofMan for this element
          const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz, ///< information about the requested element unknowns
          const XFEM::DofManager& dofman
          );
      
      //! return a string representation of this class
      std::string toString() const;
      
      /*!
       * \brief get a set with the enriched fields on this node
       * 
       * this function shall always fail, if the dofmanager not set properly
       * the number of unknowns has a default, but this one not!!!
       * 
       * \return std::set with enriched fields
       */
      const std::set<FieldEnr>& FieldEnrSetPerNode(
          const int  gid    ///< unique global node id
      ) const
      {
        map<int, const std::set<FieldEnr> >::const_iterator tmp = nodalDofSet_.find(gid);
        if (tmp == nodalDofSet_.end())
        {
          std::cout << gid << endl;
          // let this always be here so we can recognize errors early
          dserror("FieldEnrSetPerNode: global node id not found!");
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
        return numElemDof_;
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
      
      //! return reference to list of local positions in a array of dofs
      template <XFEM::PHYSICS::Field field>
      const std::vector<int>& LocalDofPosPerField() const
      {
        return paramsLocalEntries_.find(field)->second;
      }
      
      //! return unique enrichments for this proc
      const std::set<XFEM::Enrichment>& getUniqueEnrichments() const
      {
        return unique_enrichments_;
      }
      
      //! return unique enrichments for this proc
      std::set<int> getUniqueEnrichmentLabels() const
      {
        std::set<int> xlabelset;
        for(std::set<XFEM::Enrichment>::const_iterator enr = unique_enrichments_.begin();
        enr != unique_enrichments_.end();
        ++enr)
        {
          xlabelset.insert(enr->XFEMConditionLabel());
        }
        return xlabelset;
      }
      
      //! access to map
      std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> getDisTypePerFieldMap() const
      {
        return DisTypePerElementField_;
      }
      
      //! access to map
      DRT::Element::DiscretizationType getDisTypePerField(XFEM::PHYSICS::Field field) const
      {
        return DisTypePerElementField_.find(field)->second;
      }
      
      //! access to map
      const std::set<XFEM::FieldEnr> getEnrichedFieldsPerEleField(XFEM::PHYSICS::Field field) const
      {
        return enrichedFieldperPhysField_.find(field)->second;
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

      //! number of element degrees of freedom
      int numElemDof_;
      
      //! unknowns for each global node id (gid)
      std::map<int, const std::set<XFEM::FieldEnr> > nodalDofSet_;
      
      //! number of unknowns for each global node id (gid)
      std::map<int, int> nodalNumDof_;
      
      //! number of parameters for each field
      std::map<XFEM::PHYSICS::Field, int> numParamsPerField_;
      
      //! number of parameters for each field
      std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> DisTypePerElementField_;
      
      //! local (elemental) dofpos
      std::map<XFEM::PHYSICS::Field, std::vector<int> > paramsLocalEntries_;
      
      //! local (elemental) dofpos
      std::map<XFEM::PHYSICS::Field, std::set<XFEM::FieldEnr> > enrichedFieldperPhysField_;
      
      //! set of all enrichments that are available on this proc.
      std::set<XFEM::Enrichment> unique_enrichments_;
      
      //! used in constructor
      void ComputeDependendInfo(
          const DRT::Element& ele,
          const map<int, const std::set<XFEM::FieldEnr> >&  nodalDofSet,
          const std::set<XFEM::FieldEnr>& enrfieldset,
          const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
      );
      
      /// no copy of this object is allowed
      
      //! private =
      ElementDofManager operator = (const ElementDofManager& old);
      
      //! private copy constructor
      ElementDofManager(const ElementDofManager& old);
  };
  
  
  //! return assembly type depending on the DOF distribution
  XFEM::AssemblyType CheckForStandardEnrichmentsOnly(
      const ElementDofManager&   eleDofManager,
      const int                  numnode,
      const int*                 nodeids
  );
  
} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
