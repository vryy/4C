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
        const Teuchos::RCP<InterfaceHandle>       ih,
        const Teuchos::RCP<DofManager>            dofman,
        const Epetra_Map&                         olddofrowmap,
        const Epetra_Map&                         newdofrowmap,
        const std::map<DofKey<onNode>, DofGID>&   oldNodalDofDistrib,
        const std::map<DofKey<onNode>, DofGID>&   newNodalDofDistrib,
        const std::map<DofKey<onElem>, DofGID>&   oldElementalDofDistrib,
        const std::map<DofKey<onElem>, DofGID>&   newElementalDofDistrib
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
        Teuchos::RCP<Epetra_Vector>&     vector,
        const bool                       quasi_static_enr
    ) const;

    void extractDofKeysForInitialization(
        std::map<int, std::set<FieldEnr> >& unknownFieldEnr
    ) const;

    void extrapolateOldTimeStepValues(
        const Teuchos::RCP<DRT::Discretization>    bdis,
        const std::map<int,LINALG::Matrix<3,1> >&  cutterposnp,
        const Teuchos::RCP<const Epetra_Vector>    ivector,
        Teuchos::RCP<Epetra_Vector>                state_vector
    ) const;

    //! switch from old to new dof distribution
    void generateTransferInformation(
        const Teuchos::RCP<Epetra_Vector>&     vector
    ) const;

  private:

    //! disabled copy constructor
    explicit DofDistributionSwitcher(const DofDistributionSwitcher& dofdist);

    //! disabled assignment operator
    DofDistributionSwitcher operator = (const DofDistributionSwitcher& old);

    //! disabled default constructor
    explicit DofDistributionSwitcher();


    const Teuchos::RCP<InterfaceHandle>       ih_;                      ///< interface information
    const Teuchos::RCP<DofManager>            dofman_;                  ///< dofmanager for this proc
    const Epetra_Map                          olddofrowmap_;            ///< dofrowmap of the old state vector
    const Epetra_Map                          newdofrowmap_;            ///< dofrowmap of the new state vector
    const std::map<DofKey<onNode>, DofGID>    oldNodalDofDistrib_;      ///< old positions of dof in Epetra Vector per nodal dofkey
    const std::map<DofKey<onNode>, DofGID>    newNodalDofDistrib_;      ///< new positions of dof in Epetra Vector per nodal dofkey
    const std::map<DofKey<onElem>, DofGID>    oldElementalDofDistrib_;  ///< old positions of dof in Epetra Vector per element dofkey
    const std::map<DofKey<onElem>, DofGID>    newElementalDofDistrib_;  ///< new positions of dof in Epetra Vector per element dofkey

    std::map<int, std::set<FieldEnr> >  unknownFieldEnr_;
  };


  //! try to find another enrichment for this physical field
  Enrichment genAlternativeEnrichment(
      const int                             gnodeid,
      const PHYSICS::Field            oldphysvar,
      const Teuchos::RCP<DofManager>& dofman
  );


} // namespace XFEM

#endif
#endif  // #ifdef CCADISCRET
