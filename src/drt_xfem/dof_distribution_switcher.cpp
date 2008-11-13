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

#include "../drt_geometry/vector_definitions.H"
#include "dof_management.H"
#include "dof_distribution_switcher.H"
#include "dofkey.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/drt_colors.H"

//! try to find another enrichment for this physical field
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
        RCP<Epetra_Vector>&             vector,
        LINALG::Vec3                    ivalrigid_body
        ) const
{
    // create new vector with new number of dofs 
    RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);
    
    if (vector == null)
    {
#ifdef DEBUG
      std::cout << "created new vector with all zeros" << endl;
#endif
    }
    else
    {
        bool completely_unchanged = true;
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
            else // if dofkey has not been existed before, check for other dofs on the dofkeys node
            {
//              const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();
              
                // initialize to zero
                (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
            }
        }

        // step 2: find sucessor of old nodal dofkey to sum up values
        for (NodalDofPosMap::const_iterator olddof = oldNodalDofDistrib_.begin();
                                       olddof != oldNodalDofDistrib_.end();
                                       ++olddof)
        {
            const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
            const int olddofpos = olddof->second;
            const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();
            
            // try to find successor
            NodalDofPosMap::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
            if (newdof == newNodalDofDistrib_.end())  // if no successor found (was handled already in step 1)
            {
                // try to find another usefull value
                // current assumption: there is only one type of enrichment per node
                // no overlapping enrichments allowed for now
                const int nodegid = olddofkey.getGid();
                const BlitzVec3 actpos(toBlitzArray(ih_->xfemdis()->gNode(nodegid)->X()));
                const XFEM::Enrichment oldenr(olddofkey.getFieldEnr().getEnrichment());
                //const double enrval = oldenr.EnrValue(actpos, *ih_, Enrichment::approachUnknown);
                
                // create alternative dofkey
                XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, oldphysvar, dofman_));
                
                if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
                {
                    // find dof position of alternative key
                    const XFEM::FieldEnr altfieldenr(olddofkey.getFieldEnr().getField(), altenr);
                    const XFEM::DofKey<XFEM::onNode> altdofkey(nodegid, altfieldenr);
                    const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;
                    
                    //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
                    if (newdofpos < 0)
                    {
                      std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
                      std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
                      dserror("bug!");
                    }
                    
                    // add old value to already existing values
                    //(*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofrowmap_.LID(olddofpos)];
                    (*newVector)[newdofrowmap_.LID(newdofpos)] += (*oldVector)[olddofrowmap_.LID(olddofpos)];
                    completely_unchanged = false;
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
#if 0
        // step 3: PaveMaker - level all values by using interface velocities
        for (NodalDofPosMap::const_iterator newdof = newNodalDofDistrib_.begin();
                                       newdof != newNodalDofDistrib_.end();
                                       ++newdof)
        {
          const XFEM::DofKey<XFEM::onNode> newdofkey = newdof->first;
          const int newdofpos = newdof->second;
          const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();
          const int node_gid = newdofkey.getGid();
          
          const BlitzVec3 nodalpos(toBlitzArray(ih_->xfemdis()->gNode(node_gid)->X()));
          const int label = ih_->PositionWithinConditionNP(nodalpos);
          bool in_fluid = false;
          if (label == 0)
          {
            in_fluid = true;
          }
          else
          {
            in_fluid = false;
          }
          if (not in_fluid)
          {
            // reset with interface velocity - how to get velocity in non-rigid case?
            //cout << "reset" << endl;
            if (field == XFEM::PHYSICS::Velx)
              (*newVector)[newdofrowmap_.LID(newdofpos)] = ivalrigid_body(0);
            else if (field == XFEM::PHYSICS::Vely)
              (*newVector)[newdofrowmap_.LID(newdofpos)] = ivalrigid_body(1);
            else if (field == XFEM::PHYSICS::Velz)
              (*newVector)[newdofrowmap_.LID(newdofpos)] = ivalrigid_body(2);
            else
            {
              // keep previously set values
            }
          }
        }
#endif
#if 1
        
        // step 4: find predecessor of new elemental dofkey
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
                completely_unchanged = false;
            }
        }
#endif
        
#ifdef DEBUG
        if (completely_unchanged)
          cout << "completely unchanged vector" << endl;
        else
          cout << "modified vector" << endl;
#endif
    }
    
    
    // set vector to zero or initialized vector
    vector = newVector;
}


////////////////////////////////////////////
void XFEM::DofDistributionSwitcher::generateTransferInformation(
        const RCP<Epetra_Vector>&             vector
        ) const
{
    // create new vector with new number of dofs 
    RCP<Epetra_Vector> newVector = LINALG::CreateVector(newdofrowmap_,true);
    
    if (vector == null)
    {
#ifdef DEBUG
      std::cout << "created new vector with all zeros" << endl;
#endif
    }
    else
    {
        bool completely_unchanged = true;
        const RCP<Epetra_Vector> oldVector = vector;
        const Epetra_BlockMap& oldmap = oldVector->Map();
//        std::cout << "olddofrowmap_" << endl;
//        std::cout << (olddofrowmap_) << endl;
//        std::cout << "newdofrowmap_" << endl;
//        std::cout << (newdofrowmap_) << endl;
        
        if (not oldmap.SameAs(olddofrowmap_)) dserror("bug!");
        
        map<XFEM::DofKey<XFEM::onNode>,set<XFEM::DofKey<XFEM::onNode> > > transferop;
        
        // step 1: find predecessor of new nodal dofkey
        for (NodalDofPosMap::const_iterator newdof = newNodalDofDistrib_.begin();
                                       newdof != newNodalDofDistrib_.end();
                                       ++newdof)
        {
            const XFEM::DofKey<XFEM::onNode> newdofkey = newdof->first;
//            const int newdofpos = newdof->second;
            if (newdofkey.getFieldEnr().getField() == XFEM::PHYSICS::Velx)
            {
              transferop[newdofkey] = set<XFEM::DofKey<XFEM::onNode> >();
            }
            NodalDofPosMap::const_iterator olddof = oldNodalDofDistrib_.find(newdofkey);
            if (olddof != oldNodalDofDistrib_.end())  // if dofkey has existed before, use old value
            {
                const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
//                const int olddofpos = olddof->second;
                //cout << newdofkey.toString() << " -> init to old value" << endl;
//                (*newVector)[newdofrowmap_.LID(newdofpos)] = (*oldVector)[olddofrowmap_.LID(olddofpos)];
                
                if (newdofkey.getFieldEnr().getField() == XFEM::PHYSICS::Velx)
                {
                  transferop[newdofkey].insert(olddofkey);
                }
                
                
            }
            else // if dofkey has not been existed before, check for other dofs on the dofkeys node
            {
//              const XFEM::PHYSICS::Field field = newdofkey.getFieldEnr().getField();
              
                // initialize to zero
//                (*newVector)[newdofrowmap_.LID(newdofpos)] = 0.0;
            }
        }

        // step 2: find sucessor of old nodal dofkey to sum up values
        for (NodalDofPosMap::const_iterator olddof = oldNodalDofDistrib_.begin();
                                       olddof != oldNodalDofDistrib_.end();
                                       ++olddof)
        {
            const XFEM::DofKey<XFEM::onNode> olddofkey = olddof->first;
//            const int olddofpos = olddof->second;
            const XFEM::PHYSICS::Field oldphysvar = olddofkey.getFieldEnr().getField();
            
            // try to find successor
            NodalDofPosMap::const_iterator newdof = newNodalDofDistrib_.find(olddofkey);
            if (newdof == newNodalDofDistrib_.end())  // if no successor found
            {
                // try to find another usefull value
                // current assumption: there is only one type of enrichment per node
                // no overlapping enrichments allowed for now
                const int nodegid = olddofkey.getGid();
                const BlitzVec3 actpos(toBlitzArray(ih_->xfemdis()->gNode(nodegid)->X()));
                const XFEM::Enrichment oldenr = olddofkey.getFieldEnr().getEnrichment();
                //const double enrval = oldenr.EnrValue(actpos, *ih_, Enrichment::approachUnknown);
                
                // create alternative dofkey
                XFEM::Enrichment altenr(genAlternativeEnrichment(nodegid, oldphysvar, dofman_));
                
                if (altenr.Type() != XFEM::Enrichment::typeUndefined) // if alternative key found, add old solution to it
                {
                    // find dof position of alternative key
                    const XFEM::FieldEnr altfieldenr(olddofkey.getFieldEnr().getField(), altenr);
                    const XFEM::DofKey<XFEM::onNode> altdofkey(nodegid, altfieldenr);
                    const int newdofpos = newNodalDofDistrib_.find(altdofkey)->second;
                    
                    //std::cout << olddofkey.toString() << " -> " << altdofkey.toString() << endl;
                    if (newdofpos < 0)
                    {
                      std::cout << "old Dofkey" << endl << olddofkey.toString() << endl;
                      std::cout << "alt Dofkey" << endl << altdofkey.toString() << endl;
                      dserror("bug!");
                    }
                    
                    // add old value to already existing values
                    //(*newVector)[newdofrowmap_.LID(newdofpos)] += enrval*(*oldVector)[olddofrowmap_.LID(olddofpos)];
//                    (*newVector)[newdofrowmap_.LID(newdofpos)] += (*oldVector)[olddofrowmap_.LID(olddofpos)];
                    completely_unchanged = false;
                    if (altdofkey.getFieldEnr().getField() == XFEM::PHYSICS::Velx)
                    {
                      transferop[altdofkey].insert(olddofkey);
                    }
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
        
#ifdef DEBUG
        if (completely_unchanged)
          cout << "completely unchanged vector" << endl;
        else
          cout << "modified vector" << endl;
#endif
        map<XFEM::DofKey<XFEM::onNode>,set<XFEM::DofKey<XFEM::onNode> > >::const_iterator entry;
        for (entry = transferop.begin(); entry!=transferop.end(); ++entry)
        {
          if (entry->second.size() > 1)
            dserror("only one successor allowed for now");
          
          if (entry->second.empty())
          {
            cout << RED << entry->first << END_COLOR << endl;
          }
          else
          {
            set<XFEM::DofKey<XFEM::onNode> >::const_iterator first_set_entry = entry->second.begin();
            XFEM::DofKey<XFEM::onNode> old = *entry->second.begin();
            if (entry->first != old)
            {
              cout << GREEN << entry->first << END_COLOR << " <-- " << GREEN << old << END_COLOR << endl;
            }
          }
        }
        //exit(1);
    }
}



#endif  // #ifdef CCADISCRET
