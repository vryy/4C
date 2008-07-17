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

#include "interface.H"


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
            const map<DofKey<onNode>, DofPos>&   oldNodalDofDistrib,
            const map<DofKey<onNode>, DofPos>&   newNodalDofDistrib,
            const map<DofKey<onElem>, DofPos>&   oldElementalDofDistrib,
            const map<DofKey<onElem>, DofPos>&   newElementalDofDistrib
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
        // derive information for dofkeys per node
        map<DofKey<onNode>, DofPos>::const_iterator dofkeyiter;
        for (dofkeyiter = oldNodalDofDistrib.begin(); dofkeyiter != oldNodalDofDistrib.end(); ++dofkeyiter)
        {
          const DofKey<onNode> dofkey = dofkeyiter->first;
          const int gid = dofkey.getGid();
          const XFEM::FieldEnr fieldenr = dofkey.getFieldEnr();
          oldEnrFieldSetNode_[gid].insert(fieldenr);
        }
        return;
    }
            
    //! destructor
    ~DofDistributionSwitcher()
    {
        return;
    }
    
    //! switch from old to new dof distribution
    void mapVectorToNewDofDistribution(
            RCP<Epetra_Vector>&     vector
            ) const
    {
      BlitzVec3 dummy;
      dummy(0) = 0.0;
      dummy(1) = 0.0;
      dummy(2) = 0.0;
      mapVectorToNewDofDistribution(vector,dummy);
    }
    
    //! switch from old to new dof distribution
    void mapVectorToNewDofDistribution(
            RCP<Epetra_Vector>&     vector,
            BlitzVec3               ivelnrigid_body
            ) const;
    
private:
    //! disabled copy constructor
    explicit DofDistributionSwitcher(const XFEM::DofDistributionSwitcher& dofdist);
    //! disabled default constructor
    explicit DofDistributionSwitcher();

    
    const RCP<XFEM::InterfaceHandle>     ih_;                      ///< interface information
    const RCP<XFEM::DofManager>          dofman_;                  ///< dofmanager for this proc
    const Epetra_Map                     olddofrowmap_;            ///< dofrowmap of the old state vector
    const Epetra_Map                     newdofrowmap_;            ///< dofrowmap of the new state vector
    const map<DofKey<onNode>, DofPos>&   oldNodalDofDistrib_;      ///< old positions of dof in Epetra Vector per nodal dofkey
    const map<DofKey<onNode>, DofPos>&   newNodalDofDistrib_;      ///< new positions of dof in Epetra Vector per nodal dofkey
    const map<DofKey<onElem>, DofPos>&   oldElementalDofDistrib_;  ///< old positions of dof in Epetra Vector per element dofkey
    const map<DofKey<onElem>, DofPos>&   newElementalDofDistrib_;  ///< new positions of dof in Epetra Vector per element dofkey
    map<int,std::set<XFEM::FieldEnr> >   oldEnrFieldSetNode_;      ///< old enriched fields per node 
};

} // namespace XFEM

#endif
#endif  // #ifdef CCADISCRET
