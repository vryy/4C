/*!
\file dof_management_element.H

\brief provides the element dofmanager classes

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 */


#ifndef DOF_MANAGEMENT_ELEMENT_H
#define DOF_MANAGEMENT_ELEMENT_H


#include "physics.H"
#include "../drt_lib/drt_element.H"


namespace XFEM
{
  class DofManager;
  class Enrichment;
  class FieldEnr;

  /*!
   * \brief It knows everything about the XFEM enrichments for a particular interface setup for one given element
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
          const std::map<int, const std::set<XFEM::FieldEnr> >&  nodalDofSet,
          const std::set<XFEM::FieldEnr>& enrfieldset,
          const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
      );

      //! for a given element, collect all information about Dofs for this element and return it
      explicit ElementDofManager(
          const DRT::Element&                                                     ele,            ///< construct EleDofMan for this element
          const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz, ///< information about the requested element unknowns
          const XFEM::DofManager& dofman
          );

      //! return a std::string representation of this class
      std::string toString() const;

      /*!
       * \brief get a set with the enriched fields on this node
       *
       * this function shall always fail, if the dofmanager not set properly
       * the number of unknowns has a default, but this one not!!!
       *
       * \return std::set with enriched fields
       */
      const std::set<XFEM::FieldEnr>& FieldEnrSetPerNode(
          const int  gid    ///< unique global node id
          ) const;

      //! return enriched field for local dof id
      const XFEM::FieldEnr& FieldEnrSetPerDof(
          const int  lid    ///< local dof id
          ) const;

      //! function that returns node specific number of dofs or the given default number
      std::size_t NumDofPerNode(
          const int gid             ///< unique global node id
      ) const;

      //! return number of dofs at nodes
      std::size_t NumNodeDof() const { return numNodeDof_; };

      //! return number of dofs at elements
      std::size_t NumElemDof() const { return numElemDof_; };

      //! return number of dofs at nodes + elements
      std::size_t NumDofElemAndNode() const { return numNodeDof_ + numElemDof_; };

      //! function that returns dofs for each node of this particular element
      const std::map<int, const std::set<XFEM::FieldEnr> >& getNodalDofSet() const;

      /*! return number of dofs on this element for a given field. E.g., if we have 4 nodes with
       *  standard x displacements at each node and 2 enriched x displacements, this function will return 6
       *  we do not distinguish between nodal and element dofs
       */
      std::size_t NumDofPerField(
          const XFEM::PHYSICS::Field  field  ///< field for which we seek the number of DOFs
      ) const;

      // get number of fields
      int NumFields() const
      {
        return paramsLocalEntries_.size();
      };

      // get vector of available fields
      const std::vector<XFEM::PHYSICS::Field> GetFields() const
      {
        std::vector<XFEM::PHYSICS::Field> fields;
        std::map<XFEM::PHYSICS::Field, std::vector<int> >::const_iterator it;
        for (it=paramsLocalEntries_.begin(); it!=paramsLocalEntries_.end();it++)
          fields.push_back(it->first);
        return fields;
      };

      //! return reference to list of local positions in a array of dofs
      const std::vector<int>& LocalDofPosPerField(const XFEM::PHYSICS::Field field) const;

      //! return reference to list of local positions in a array of dofs
      template <XFEM::PHYSICS::Field field>
      const std::vector<int>& LocalDofPosPerField() const
      {
        std::map<XFEM::PHYSICS::Field, std::vector<int> >::const_iterator tmp = paramsLocalEntries_.find(field);
        if (tmp == paramsLocalEntries_.end()){
          dserror( "no local dofs for field %d", field );
        }
        return tmp->second;
      }

      //! return unique enrichments for this proc
      const std::set<XFEM::Enrichment>& getUniqueEnrichments() const;

      //! return unique enrichments for this proc
      std::set<int> getUniqueEnrichmentLabels() const;

      //! access to map
      DRT::Element::DiscretizationType getDisTypePerField(XFEM::PHYSICS::Field field) const;

      //! access to map
      const std::set<XFEM::FieldEnr>& getEnrichedFieldsPerEleField(XFEM::PHYSICS::Field field) const;

      //! not equal operator (for use in STL iterators)
      bool operator !=(const ElementDofManager& rhs) const;

    private:

      //! unknowns for each global node id (gid)
      std::map<int, const std::set<XFEM::FieldEnr> > nodalDofSet_;

      //! number of unknowns for each global node id (gid)
      //std::map<int, std::size_t> nodalNumDof_;

      //! number of parameters for each field
      //std::map<XFEM::PHYSICS::Field, std::size_t> numParamsPerField_;

      //! number of parameters for each field
      std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> DisTypePerElementField_;

      //! local (elemental) dofpos
      std::map<XFEM::PHYSICS::Field, std::vector<int> > paramsLocalEntries_;

      //! local (elemental) dofpos
      std::map<XFEM::PHYSICS::Field, std::set<XFEM::FieldEnr> > enrichedFieldperPhysField_;

      //! field and enrichment for local dof (nodal dofs)
      std::vector<XFEM::FieldEnr> lidtofieldenr_;

      //! set of all enrichments that are available on this proc.
      std::set<XFEM::Enrichment> unique_enrichments_;

      //! total number of nodal degrees of freedom
      std::size_t numNodeDof_;

      //! total number of element degrees of freedom
      std::size_t numElemDof_;

      //! used in constructor
      void ComputeDependentInfo(
          const DRT::Element& ele,
          const std::map<int, const std::set<XFEM::FieldEnr> >&  nodalDofSet, ///< std::map<NodeGid,std::set<FieldEnr>>
          const std::set<XFEM::FieldEnr>& enrfieldset,
          const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
      );

      /// no copy of this object is allowed

      //! private =
      ElementDofManager operator = (const ElementDofManager& old);

      //! private copy constructor
      ElementDofManager(const ElementDofManager& old);
  };

} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
