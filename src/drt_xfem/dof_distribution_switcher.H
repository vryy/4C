/*!
\file dof_distribution_switcher.H

\brief provides the dofmanager classes

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/


#ifndef DOF_DISTRIBUTION_SWITCHER_H
#define DOF_DISTRIBUTION_SWITCHER_H


#include <Teuchos_RCP.hpp>
#include "xfem_utils_combust.H"


class Epetra_Map;
class Epetra_Vector;

namespace COMBUST
{
  class InterfaceHandleCombust;
}

namespace DRT
{
  class Discretization;
}

namespace XFEM
{
  class DofKey;
  class DofManager;
  class FieldEnr;

  /*!
  \brief this class is used to switch from one distribution to a second

  */
  class DofDistributionSwitcher
  {
  public:

    //! constructor
    explicit DofDistributionSwitcher(
        const Teuchos::RCP<COMBUST::InterfaceHandleCombust> ih,
        const Teuchos::RCP<DofManager>             dofman,
        const Epetra_Map&                          olddofrowmap,
        const Epetra_Map&                          newdofrowmap,
        const std::map<DofKey, DofGID>&    oldNodalDofDistrib,
        const std::map<DofKey, DofGID>&    newNodalDofDistrib
    ) :
      ih_(ih),
      dofman_(dofman),
      olddofrowmap_(olddofrowmap),
      newdofrowmap_(newdofrowmap),
      oldNodalDofDistrib_(oldNodalDofDistrib),
      newNodalDofDistrib_(newNodalDofDistrib)
    {
      extractDofKeysForInitialization(unknownFieldEnr_);
#ifdef DEBUG
      if (unknownFieldEnr_.size()>0)
        std::cout << "Number of nodes with missing nodal unknowns: " << unknownFieldEnr_.size() << std::endl;
#endif
      return;
    }

    //! destructor
    ~DofDistributionSwitcher()
    {
      return;
    }

    //! switch from old to new dof distribution for jump/kink enrichment problems
    void mapVectorToNewDofDistributionCombust(
        Teuchos::RCP<Epetra_Vector>&     vector,
        const bool                       quasi_static_enr
    ) const;

    void extractDofKeysForInitialization(
        std::map<int, std::set<FieldEnr> >& unknownFieldEnr
    ) const;

    void GmshOutput(
      int                               step,
      const Teuchos::RCP<DRT::Discretization>    patchdis,
      std::map<int,int>                mapforoutput
    ) const;

  private:

    //! disabled copy constructor
    explicit DofDistributionSwitcher(const DofDistributionSwitcher& dofdist);

    //! disabled assignment operator
    DofDistributionSwitcher operator = (const DofDistributionSwitcher& old);

    //! disabled default constructor
    explicit DofDistributionSwitcher();

    const Teuchos::RCP<COMBUST::InterfaceHandleCombust> ih_;                      ///< interface information
    const Teuchos::RCP<DofManager>                      dofman_;                  ///< dofmanager for this proc
    const Epetra_Map&                                   olddofrowmap_;            ///< dofrowmap of the old state vector
    const Epetra_Map&                                   newdofrowmap_;            ///< dofrowmap of the new state vector
    const std::map<DofKey, DofGID>              oldNodalDofDistrib_;      ///< old positions of dof in Epetra Vector per nodal dofkey
    const std::map<DofKey, DofGID>              newNodalDofDistrib_;      ///< new positions of dof in Epetra Vector per nodal dofkey

    std::map<int, std::set<FieldEnr> >  unknownFieldEnr_;
  };
} // namespace XFEM

#endif
