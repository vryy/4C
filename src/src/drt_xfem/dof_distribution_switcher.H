/*!
\file dof_distribution_switcher.H

\brief provides the dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef DOF_DISTRIBUTION_SWITCHER_H
#define DOF_DISTRIBUTION_SWITCHER_H

#include "interfacexfsi.H"
#include "dofkey.H"


namespace XFEM
{

  /*!
  \brief this class is used to switch from one distribution to a second

  */
  class DofDistributionSwitcher
  {
  public:

    //! constructor
    explicit DofDistributionSwitcher(
        const RCP<XFEM::InterfaceHandle>     ih,
        const RCP<XFEM::DofManager>          dofman,
        const Epetra_Map&                    olddofrowmap,
        const Epetra_Map&                    newdofrowmap,
        const map<DofKey<onNode>, DofGID>&   oldNodalDofDistrib,
        const map<DofKey<onNode>, DofGID>&   newNodalDofDistrib,
        const map<DofKey<onElem>, DofGID>&   oldElementalDofDistrib,
        const map<DofKey<onElem>, DofGID>&   newElementalDofDistrib
    ) :
      ih_(ih),
      dofman_(dofman),
      olddofrowmap_(olddofrowmap),
      newdofrowmap_(newdofrowmap),
      oldNodalDofDistrib_(oldNodalDofDistrib),
      newNodalDofDistrib_(newNodalDofDistrib),
      oldElementalDofDistrib_(oldElementalDofDistrib),
      newElementalDofDistrib_(newElementalDofDistrib)
    {
      extractDofKeysForInitialization(unknownFieldEnr_);
      cout << "Number of nodes with missing nodal unknowns: " << unknownFieldEnr_.size() << endl;
      return;
    }


    //! destructor
    ~DofDistributionSwitcher()
    {
      return;
    }


    //! switch from old to new dof distribution
    void mapVectorToNewDofDistribution(
        Teuchos::RCP<Epetra_Vector>&     vector
    ) const;

    //! switch from old to new dof distribution for jump/kink enrichment problems
    void mapVectorToNewDofDistributionCombust(
        Teuchos::RCP<Epetra_Vector>&     vector
    ) const;

    void extractDofKeysForInitialization(
        std::map<int, set<XFEM::FieldEnr> >& unknownFieldEnr
    ) const;

    void extrapolateOldTimeStepValues(
        const RCP<DRT::Discretization>    bdis,
        const std::map<int,LINALG::Matrix<3,1> >&  cutterposnp,
        const RCP<const Epetra_Vector>    ivector,
        RCP<Epetra_Vector>                state_vector
    ) const;

    //! switch from old to new dof distribution
    void generateTransferInformation(
        const RCP<Epetra_Vector>&     vector
    ) const;

  private:

    //! disabled copy constructor
    explicit DofDistributionSwitcher(const XFEM::DofDistributionSwitcher& dofdist);

    //! disabled assignment operator
    DofDistributionSwitcher operator = (const DofDistributionSwitcher& old);

    //! disabled default constructor
    explicit DofDistributionSwitcher();


    const RCP<XFEM::InterfaceHandle>     ih_;                      ///< interface information
    const RCP<XFEM::DofManager>          dofman_;                  ///< dofmanager for this proc
    const Epetra_Map                     olddofrowmap_;            ///< dofrowmap of the old state vector
    const Epetra_Map                     newdofrowmap_;            ///< dofrowmap of the new state vector
    const map<DofKey<onNode>, DofGID>    oldNodalDofDistrib_;      ///< old positions of dof in Epetra Vector per nodal dofkey
    const map<DofKey<onNode>, DofGID>    newNodalDofDistrib_;      ///< new positions of dof in Epetra Vector per nodal dofkey
    const map<DofKey<onElem>, DofGID>    oldElementalDofDistrib_;  ///< old positions of dof in Epetra Vector per element dofkey
    const map<DofKey<onElem>, DofGID>    newElementalDofDistrib_;  ///< new positions of dof in Epetra Vector per element dofkey

    std::map<int, set<XFEM::FieldEnr> >      unknownFieldEnr_;
    //    map<int,std::set<XFEM::FieldEnr> >   oldEnrFieldSetNode_;      ///< old enriched fields per node
  };


  //! try to find another enrichment for this physical field
  XFEM::Enrichment genAlternativeEnrichment(
      const int                    gnodeid,
      const XFEM::PHYSICS::Field   oldphysvar,
      const RCP<XFEM::DofManager>& dofman
  );


} // namespace XFEM

#endif
#endif  // #ifdef CCADISCRET
