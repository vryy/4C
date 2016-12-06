/*!
\file dof_management.H

\brief provides the dofmanager class

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 */


#ifndef DOF_MANAGEMENT_H
#define DOF_MANAGEMENT_H


#include "enrichment.H"
#include "physics.H"


namespace COMBUST
{
  class InterfaceHandleCombust;
}

namespace DRT
{
  class DofSet;
}

namespace XFEM
{
  class DofKey;
  class ElementAnsatz;
  class FieldEnr;

  /*!
 \brief He knows everything about the XFEM enrichments for one particular interface setup

   */
  class DofManager
  {
    public:

      //! default constructor
      explicit DofManager();

      //! constructor (for combustion problems)
      explicit DofManager(
          const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle, ///< geometrical interface information
          const Teuchos::RCP<Epetra_Vector>&                   phinp,
          const std::set<XFEM::PHYSICS::Field>&                fieldset,        ///< physical fields to be enriched
          const Teuchos::ParameterList&                        params,          ///< XFEM parameters
          const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap           ///< nodes on periodic boundaries
          );

      //! copy constructor
      explicit DofManager(const XFEM::DofManager& dofman);

      //! destructor
      ~DofManager();

      //! return std::string representation of the dofmanager
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
        return emptyset_;
      };

      //! fill maps with nodal and element dof positions in EpetraVectors for all DofKeys on this proc
      void fillDofRowDistributionMaps(
          std::map<XFEM::DofKey, XFEM::DofGID>&  NodalDofDistributionMap        ///< nodal dofs
      ) const;

      // fill column map with nodal dof positions in EpetraVectors for all DofKeys on this proc
      void fillNodalDofColDistributionMap(
          std::map<XFEM::DofKey, XFEM::DofGID>&  NodalDofColDistributionMap
      ) const;

      /*!
       \brief transform XFEM vector to (standard FEM) output vector

       This function assumes that a so-called modified enrichment strategy was used, which guarantees
       that the values of standard dofs represent the enriched XFEM field. That is, all enrichment
       values at nodes are identically 0!

       An XFEM vector holds a different number of unknowns for each node. To write the output, a
       uniform number of degrees of freedom is expected for each node. The XFEM vector is trans-
       formed to the original (standard FEM) layout (dofsetout) with one dof per physical field per
       node. The desired physical fields for output can be specified (fieldsout).

       This function writes the value of the XFEM standard enrichment for every physical field to
       the corresponding standard FEM dof in the output vector.

       \author henke 07/09
      */
      Teuchos::RCP<Epetra_Vector> transformXFEMtoStandardVector(
          const Epetra_Vector&                        originalvector,
          const DRT::DofSet&                          dofsetout,
          const std::map<XFEM::DofKey, XFEM::DofGID>& nodalDofDistributionMap,
          const std::set<XFEM::PHYSICS::Field>&       fieldsout
      ) const;

      /*!
        \author tansfor from Standard to XFEM
       */
       Teuchos::RCP<Epetra_Vector> transformStandardToXFEMVector(
           const Epetra_Vector&                        originalvector,
           const DRT::DofSet&                          dofsetout,
           const std::map<XFEM::DofKey, XFEM::DofGID>& nodalDofDistributionMap,
           const std::set<XFEM::PHYSICS::Field>&       fieldsout,
           const Epetra_BlockMap*                      xfemdofrowmap
       ) const;



      /*!
       \brief overwrite the values of a certain physical field and enrichment
      */
      void overwritePhysicalField(
          Teuchos::RCP<Epetra_Vector>&                vector,
          const std::map<XFEM::DofKey, XFEM::DofGID>& nodalDofDistributionMap,
          const XFEM::PHYSICS::Field&                 physfield,
          const XFEM::Enrichment::EnrType&            enrichment,
          const double                                value,
          const bool                                  strict = false
      ) const;

    private:

      DofManager operator = (const DofManager& old);

      /// find unique enrichments among all dofs in this DofManager
      void GatherUniqueEnrichments() const;

      //! pointer to the xfem discretization of which we hold the unknowns
      const Teuchos::RCP<COMBUST::InterfaceHandleCombust> ih_;

      //! connection between master and slave row nodes
      const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap_;

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

      /// when returning a reference to a set, this set can be used when an empty set is needed
      /// avoids creating a lot of empty sets -> speed
      const std::set<XFEM::FieldEnr>                    emptyset_;

  };

} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
