/*!
\file dof_management.H

\brief provides the dofmanager class

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

#include "dofkey.H"
#include "element_ansatz.H"
// Der ist hier nur drin, damit ich auf die FlameFront im InterfaceHandle zugreifen kann
#include "../drt_combust/combust_interface.H"

namespace XFEM
{

  /*!
 \brief He knows everything about the XFEM enrichments for one particular interface setup

   */
  class DofManager
  {
    public:

      //! default constructor
      explicit DofManager();

      //! constructor
      explicit DofManager(
          const RCP<XFEM::InterfaceHandle>&      ih,
          const std::set<XFEM::PHYSICS::Field>&  fieldset,
          const XFEM::ElementAnsatz&             element_ansatz,
          const Teuchos::ParameterList&          params
          );

      //! constructor (for combustion problems)
      explicit DofManager(
          const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle, ///< geometrical interface information
          const std::set<XFEM::PHYSICS::Field>&      fieldset,        ///< physical fields to be enriched
          const Teuchos::ParameterList&              params           ///< XFEM parameters
          );

      //! copy constructor
      explicit DofManager(const XFEM::DofManager& dofman);

      //! destructor
      ~DofManager();

      //! return string representation of the dofmanager
      std::string toString() const;


      //! print a Gmsh representation of this class
      void toGmsh(
          const int step                                 ///< step for file numbering
          ) const;

      //! return global (i.e. from all procs) number of nodal dofs only
      int NumNodalDof() const;


      //! return set with enriched fields for requested node
      const std::set<FieldEnr>& getNodeDofSet(
          const int gid      ///< global node id
          ) const
      {
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
        if (entry == nodalDofSet_.end())
        {
          return emptyset_;
        }
        return entry->second;
      };


      //! return set with enriched fields for requested element
      const std::set<XFEM::FieldEnr>& getElementDofSet(
          const int gid      ///< global element id
          ) const
      {
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = elementalDofs_.find(gid);
        if (entry == elementalDofs_.end())
        {
          return emptyset_;
        }
        return entry->second;
      };

//      //! check whether given ElementDofManager fits the given element
//      void checkForConsistency(
//          const DRT::Element&            ele,        ///< element for which one wants to check the eledofmanager
//          const XFEM::ElementDofManager& eledofman   ///< eledofmanager in question
//      ) const;

      //! fill maps with nodal and element dof positions in EpetraVectors for all DofKeys on this proc
      void fillDofDistributionMaps(
          map<DofKey<onNode>, DofGID>&  NodalDofDistributionMap,        ///< nodal dofs
          map<DofKey<onElem>, DofGID>&  ElementalDofDistributionMap     ///< element dofs
      ) const;


      //! for output XFEM is turned of and a uniform field is drawn. This method returns a vector with uniform number of unknowns per node
      Teuchos::RCP<Epetra_Vector> fillPhysicalOutputVector(
          const Epetra_Vector&                    original_vector,
          const DRT::DofSet&                      dofset_out,
          const map<DofKey<onNode>, DofGID>&      nodalDofDistributionMap,
          const std::set<XFEM::PHYSICS::Field>&   fields_out
      ) const;

      //! return pointer to interface handle
      const Teuchos::RCP<XFEM::InterfaceHandle> getInterfaceHandle() const { return ih_; }

    private:
      
      DofManager operator = (const DofManager& old);

      /// find unique enrichments among all dofs in this DofManager
      void GatherUniqueEnrichments() const;

      //! pointer to the xfem discretization of which we hold the unknowns
      const Teuchos::RCP<XFEM::InterfaceHandle>         ih_;

      /*! unknowns for each global node id (gid) on this processor
       * if it does not contain the key, the number of unknowns is supposed to be zero
       * meaning, the the dofmanager will not crash but give back an empty set of unknowns
       * see also getNodeDofSet()
       */
      std::map<int, const std::set<XFEM::FieldEnr> >    nodalDofSet_;

      /*! unknowns for each global element id (gid) on this processor
       * same as above, see getElementDofSet()
       */
      std::map<int, const std::set<XFEM::FieldEnr> >    elementalDofs_;

      //! only used for combustion! henke 03/09
      std::map<int, std::set<XFEM::FieldEnr> >    nodeDofMap_;

      //! only used for combustion! henke 03/09
      std::map<int, std::set<XFEM::FieldEnr> >    elementDofMap_;

      /// when returning a reference to a set, this set can be used when an empty set is needed
      /// avoids creating a lot of empty sets -> speed
      const std::set<XFEM::FieldEnr>                    emptyset_;

  };

} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
