/*!
\file dof_distribution_switcher.cpp

\brief provides the dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include "xfem.H"
#include "dof_management.H"
#include "dof_distribution_switcher.H"
#include "dofkey.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"


static XFEM::Enrichment genAlternativeEnrichment(
        const int                    gnodeid,
        const XFEM::PHYSICS::Field   oldphysvar,
        const RCP<XFEM::DofManager>  dofman
        )
{
    std::set<XFEM::FieldEnr> fieldset = dofman->getNodeDofSet(gnodeid);
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenriter = fieldset.begin(); fieldenriter != fieldset.end(); ++fieldenriter)
    {
        const XFEM::PHYSICS::Field physvar = fieldenriter->getField();
        if (oldphysvar == physvar)
        {
            return fieldenriter->getEnrichment();
            break;
        }
    }
    return XFEM::Enrichment();
}

void XFEM::DofDistributionSwitcher::mapVectorToNewDofDistribution(
        RCP<Epetra_Vector>&             vector
        ) const
{
    // create new vector with new number of dofs 
    RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);
    
    if (vector == null)
    {
      std::cout << "created new vector with all zeros" << endl;
    }
    else
    {
        const RCP<Epetra_Vector> oldVector = vector;
        const Epetra_BlockMap& oldmap = oldVector->Map();
//        std::cout << "olddofrowmap_" << endl;
//        std::cout << (olddofrowmap_) << endl;
//        std::cout << "newdofrowmap_" << endl;
//        std::cout << (newdofrowmap_) << endl;
        
        if (not oldmap.SameAs(olddofrowmap_)) dserror("bug!");
        
        // step 1: find predecessor of new nodal dofkey
        for (NodalDofPosMap::const_iterator newdof = newNodalDofDistrib_.begin();
                                       newdof != newNodalDofDistrib_.end();
                                       ++newdof)
        {
            const XFEM::DofKey<XFEM::onNode> newdofkey = newdof->first;
            const int newdofpos = newdof->second;
            
            NodalDofPosMap::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
            if (olddof != oldNodalDofDistrib_.end())  // if dofkey has existed before, use old value
            {
                const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
                const int olddofpos = olddof->second;
                //cout << newdofkey.toString() << " -> init to old value" << endl;
                (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
            }
            else // if dofkey has not been existed before, initialize to zero
            {
                //cout << newdofkey.toString() << " -> init to zero" << endl;
                (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
            }
        }

        // step 2: find sucessor of old nodal dofkey to summ up values
        for (NodalDofPosMap::const_iterator olddof = oldNodalDofDistrib_.begin();
                                       olddof != oldNodalDofDistrib_.end();
                                       ++olddof)
        {
            const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
            const int olddofpos = olddof->second;
            const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();
            
            // try to find successor
            NodalDofPosMap::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
            if (newdof == newNodalDofDistrib_.end())  // if no successor found
            {
                const int gnodeid = olddofkey.getGid();
                const BlitzVec3 actpos(toBlitzArray(ih_->xfemdis()->gNode(gnodeid)->X()));
                const XFEM::FieldEnr oldfieldenr = olddofkey.getFieldEnr();
                const XFEM::Enrichment oldenr = oldfieldenr.getEnrichment();
                const double enrval = oldenr.EnrValue(actpos, *ih_, Enrichment::approachUnknown);
                
                // create alternative dofkey
                XFEM::Enrichment altenr(genAlternativeEnrichment(gnodeid, oldphysvar, dofman_));
                
                if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
                {
                    // find dof position of alternative key
                    const XFEM::FieldEnr altfieldenr(oldfieldenr.getField(), altenr);
                    const XFEM::DofKey<XFEM::onNode> altdofkey(gnodeid, altfieldenr);
                    const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;
                    
                    //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
                    if (newdofpos < 0)
                    {
                      std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
                      std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
                      dserror("bug!");
                    }
                    
                    // add old value to already existing values
                    (*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofpos];
                }
                else // if not alternative is found
                {
                    // this can only happen in the void enrichment case and in that case,
                    // the dof value is zero anyway, which coincides with the fact that we have no place,
                    // where we could store it ;-) 
                }
            }
            else
            {
              // do nothing, this case was handled in step 1
            }
        }
        
        // step 3: find predecessor of new elemetal dofkey
        for (ElementalDofPosMap::const_iterator newdof = newElementalDofDistrib_.begin();
                                       newdof != newElementalDofDistrib_.end();
                                       ++newdof)
        {
            const XFEM::DofKey<XFEM::onElem> newdofkey = newdof->first;
            const int newdofpos = newdof->second;
            
            ElementalDofPosMap::const_iterator olddof = oldElementalDofDistrib_.find(newdofkey);
            if (olddof != oldElementalDofDistrib_.end())  // if dofkey has existed before, use old value
            {
                const XFEM::DofKey<XFEM::onElem> olddofkey = olddof->first;
                const int olddofpos = olddof->second;
                //cout << "init to old value" << endl;
                (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
            }
            else // if dofkey has not been existed before, initialize to zero
            {
                //cout << "init to zero" << endl;
                (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
            }
        }
        //exit(1);
    }
    
    
    // set vector to zero or initialized vector
    vector = newVector;
}






#endif  // #ifdef CCADISCRET
