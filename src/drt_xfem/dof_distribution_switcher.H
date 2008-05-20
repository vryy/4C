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

//#include "xfem.H"
//#include "physics.H"
//#include "enrichment.H"
//#include "field_enriched.H"
//#include "dofkey.H"
//#include "integrationcell.H"
#include "interface.H"
//#include "../drt_lib/drt_discret.H"


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
        return;
    }
            
    //! destructor
    ~DofDistributionSwitcher()
    {
        return;
    }
    
    void mapVectorToNewDofDistribution(
            RCP<Epetra_Vector>&     vector
            ) const;
    
    
private:
    //! disabled copy constructor
    explicit DofDistributionSwitcher(const XFEM::DofDistributionSwitcher& dofdist);
    //! disabled default constructor
    explicit DofDistributionSwitcher();

    
    const RCP<XFEM::InterfaceHandle>     ih_;
    const RCP<XFEM::DofManager>          dofman_;
    const Epetra_Map                     olddofrowmap_;
    const Epetra_Map                     newdofrowmap_;
    const map<DofKey<onNode>, DofPos>&   oldNodalDofDistrib_;
    const map<DofKey<onNode>, DofPos>&   newNodalDofDistrib_;
    const map<DofKey<onElem>, DofPos>&   oldElementalDofDistrib_;
    const map<DofKey<onElem>, DofPos>&   newElementalDofDistrib_;
};

} // namespace XFEM

#endif
#endif  // #ifdef CCADISCRET
