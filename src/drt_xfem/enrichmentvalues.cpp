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


#include "../drt_combust/combust_defines.H"
#include "dof_management_element.H"
#include "dof_management.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "enrichmentvalues.H"
#include <iostream>



XFEM::Enrichmentvalues::Enrichmentvalues(
    const RCP<DRT::Discretization> discret,
    const RCP<DofManager> olddofman,
    const RCP<DofManager> dofman,
    vector<RCP<Epetra_Vector> > oldVectors,
    vector<RCP<Epetra_Vector> > newVectors,
    const RCP<COMBUST::FlameFront> flamefront,
    const RCP<InterfaceHandle> ihold,
    const RCP<InterfaceHandle> ih,
    const Epetra_Map& olddofcolmap,
    const map<DofKey<onNode>, DofGID>& oldNodalDofColDistrib,
    const Epetra_Map& newdofrowmap,
    const map<DofKey<onNode>, DofGID>& newNodalDofRowDistrib
) : 
#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
  critTol_(1.0e-02),
#endif
#endif
  oldVectors_(oldVectors),
  newVectors_(newVectors),
  ih_old_(ihold),
  ih_(ih),
  discret_(discret),
  olddofman_(olddofman),
  dofman_(dofman),
  olddofcolmap_(olddofcolmap),
  oldNodalDofColDistrib_(oldNodalDofColDistrib),
  newdofrowmap_(newdofrowmap),
  newNodalDofRowDistrib_(newNodalDofRowDistrib),
  exporter_(discret_->Comm()),
  myrank_(discret_->Comm().MyPID()),
  numproc_(discret_->Comm().NumProc())
{
  // remark: in flamefront the phi vectors allways have to fit to the current
  // discretization so no mapping from node col id to phi dof col id is needed
  phiOldCol_ = flamefront->Phin(); // last phi vector in column map
  phiNewCol_ = flamefront->Phinp(); // current phi vector in column map
  
#ifdef ENR_FULL
#ifdef MIN_ENR_VALUES
  dserror("fully and minimal recomputation at same time impossible! Fix defines!");
#endif
#endif
  
  setEnrichmentValues();
  return;
}




void XFEM::Enrichmentvalues::setEnrichmentValues(
)
{
  // Initialization
#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
  computeDomainCellVolumes();
#endif
#endif
  
  
  
/*----------------------------------------------------------------*
 * first part: compute jump and kink height for every intersected *
 * element (touched elements are not counted as intersected)      *
 *----------------------------------------------------------------*/
  oldJumpAndKinkValues();
  
  
  
/*-----------------------------------------------------------*
 * second part: compute enrichment values for new interface  *
 * according to the old values. this works when around an    *
 * enriched node at new interface an intersected element at  *
 * old interface has been. else node is added to failedNodes *
 *-----------------------------------------------------------*/
  computeNewEnrichments();
  
  
  
/*---------------------------------------------------------------*
 * third part: handle the nodes which didn't get an sensible     *
 * enrichment value by the above step. look for nearest enriched *
 * node at old interface and take its surround elements          *
 *---------------------------------------------------------------*/
  handleFailedNodes();
  
  
  
  return;
} // end function setEnrichmentValues



void XFEM::Enrichmentvalues::oldJumpAndKinkValues(
)
{
  const int nsd = 3;
  const int numnode = 8;
  // velocity and pressure enrichment values of the elements nodes
  vector<LINALG::Matrix<nsd+1,numnode> > enrValues(oldVectors_.size(),LINALG::Matrix<nsd+1,numnode>(true));
  
  // remark: Finally values shall shall be computed for row nodes.
  // Therefore we loop here over col elements. If nodes of this col
  // element are not in the current procs col node map, the according
  // col element is not needed finally (test within this function)
  for (int iele=0;iele<discret_->NumMyColElements();iele++)
  {
    const DRT::Element* ele = discret_->lColElement(iele);
    
    if (ih_old_->ElementBisected(ele) || ih_old_->ElementTouchedMinus(ele)) // element is intersected
    {
      bool allEleNodesOnProc = true;
      LINALG::Matrix<1,numnode> dmin(true);
      LINALG::Matrix<nsd,1> normal(true);
      
      getDataOfIntersectedEle(ele,dmin,allEleNodesOnProc,enrValues,normal);
      
      int numnodeused = numnode;
      analyseEnrichments(ele,dmin,numnodeused);
//      cout << numnodeused << " nodes are used in " << *ele << " and dmins are\n" << dmin;

      if(allEleNodesOnProc) // if this is false element is not needed
      {
#ifndef ENR_VEL_JUMP_SCALAR
        if (ih_old_->ElementBisected(ele)) // compute jump and kink in bisected ele
        {
          // compute the kink and jump height at the interface in this element
          // a least squares approach is used for the computation
          LINALG::Matrix<2,2> sysmat(true);
          LINALG::Matrix<2,1> rhs(true);
          
          LINALG::Matrix<2,1> currJumpAndKink(true);
          LINALG::Matrix<2,4> currentJumpsAndKinks(true);
          vector<LINALG::Matrix<2,4> > eleJumpAndKinks;
          for (size_t field=0;field<oldVectors_.size();field++)
          {
            currentJumpsAndKinks.Clear();
            for (int ientry=0;ientry<nsd+1;ientry++) // number of enrichments for every node and field
            {
              sysmat.Clear();
              rhs.Clear();
              currJumpAndKink.Clear();
              
              sysmat(0,0) = numnodeused;
              for (int index=0;index<numnode;index++)
              {
                if (dmin(index)!=INFINITY) // if infinity then dont used node
                {
                  rhs(0) -= enrValues[field](ientry,index);
                  rhs(1) -= enrValues[field](ientry,index)*dmin(index);
                  sysmat(1,0) += dmin(index);
                  sysmat(1,1) += dmin(index)*dmin(index);
                }
              }
              sysmat(0,1) = sysmat(1,0);
              sysmat.Scale(0.5);
//              cout << "sysmat is " << sysmat << " and rhs is " << rhs << endl;
              
              sysmat.Invert();
              currJumpAndKink.Multiply(sysmat,rhs);
              
              currentJumpsAndKinks(0,ientry) = currJumpAndKink(0);
              currentJumpsAndKinks(1,ientry) = currJumpAndKink(1);
            } // end loop over enrichment number
            eleJumpAndKinks.push_back(currentJumpsAndKinks);
//            cout << "jump and kinks are " << currentJumpsAndKinks;
          } // end loop over vectorfield size
          cout << "final old jumps and kinks at bisected element are " << eleJumpAndKinks[2] << endl;
          eleJumpsAndKinks_.insert(make_pair(ele->Id(),eleJumpAndKinks));
        } // end if element bisected
        else // in touchedminus ele just compute jump since kink value = 0 (singulary sysmat)
        {
          // compute the jump height at the interface in this element as a mean value
          double sysmat = 0;
          double rhs = 0;
//          cout << "enrvalues of touched ele are " << enrValues[2];
          double currJump = 0;
          LINALG::Matrix<2,4> currentJumpsAndKinks(true);
          vector<LINALG::Matrix<2,4> > eleJumpAndKinks;
          for (size_t field=0;field<oldVectors_.size();field++)
          {
            currentJumpsAndKinks.Clear();
            for (int ientry=0;ientry<nsd+1;ientry++) // number of enrichments for every node and field
            {
              sysmat = 0.5*numnodeused;
              rhs = 0;
              for (int index=0;index<numnode;index++)
              {
                if (dmin(index)!=INFINITY) // if infinity then dont used node
                  rhs -= enrValues[field](ientry,index);
              }
//              cout << "rhs in touched ele is " << rhs << endl;
              
              currJump = rhs/sysmat;
//              cout << "curr jump in touched ele is " << currJump << endl;
              currentJumpsAndKinks(0,ientry) = currJump;
              
            } // end loop over enrichment number
            eleJumpAndKinks.push_back(currentJumpsAndKinks);
          } // end loop over vectorfield size
          cout << "final old jumps and kinks at touched element are " << eleJumpAndKinks[2] << endl;
          eleJumpsAndKinks_.insert(make_pair(ele->Id(),eleJumpAndKinks));
        } // end if element touchedminus
#else
        if (ih_old_->ElementBisected(ele)) // compute jump and kink in bisected ele
        {
          // compute the kink and jump height at the interface in this element
          // a least squares approach is used for the computation
          LINALG::Matrix<2,2> sysmat(true);
          LINALG::Matrix<2,1> rhsvel(true);
          LINALG::Matrix<2,1> rhspres(true);
          
          LINALG::Matrix<2,1> velJumpAndKink(true);
          LINALG::Matrix<2,1> presJumpAndKink(true);
          
          LINALG::Matrix<2,4> currentJumpsAndKinks(true);
          vector<LINALG::Matrix<2,4> > eleJumpAndKinks;
          for (size_t field=0;field<oldVectors_.size();field++)
          {
            currentJumpsAndKinks.Clear();
            sysmat.Clear();
            rhsvel.Clear();
            rhspres.Clear();
            velJumpAndKink.Clear();
            presJumpAndKink.Clear();
            
            sysmat(0,0) = numnodeused;
            
            for (int index=0;index<numnode;index++)
            {
              if (dmin(index)!=INFINITY) // if infinity then dont used node
              {
                sysmat(1,0) += dmin(index);
                sysmat(1,1) += dmin(index)*dmin(index);
                
                LINALG::Matrix<nsd,1> currVelEnrValues(true);
                for (int i=0;i<nsd;i++)
                  currVelEnrValues(i) = enrValues[field](i,index);
                
                double currVelEnrScalar = currVelEnrValues.Norm2();
                rhsvel(0) -= currVelEnrScalar;
                rhsvel(1) -= currVelEnrScalar*dmin(index);
                
                rhspres(0) -= enrValues[field](nsd,index);
                rhspres(1) -= enrValues[field](nsd,index)*dmin(index);
              }
            }
            sysmat(0,1) = sysmat(1,0);
            sysmat.Scale(0.5);
//              cout << "sysmat is " << sysmat << " and rhs is " << rhs << endl;
            
            sysmat.Invert();
            velJumpAndKink.Multiply(sysmat,rhsvel);
            presJumpAndKink.Multiply(sysmat,rhspres);
            
            currentJumpsAndKinks(0,0) = velJumpAndKink(0);
            currentJumpsAndKinks(1,0) = velJumpAndKink(1);
            currentJumpsAndKinks(0,1) = presJumpAndKink(0);
            currentJumpsAndKinks(1,1) = presJumpAndKink(1);
              
            eleJumpAndKinks.push_back(currentJumpsAndKinks);
          } // end loop over vectorfield size
          cout << "final old jumps and kinks at bisected element are " << eleJumpAndKinks[2] << endl;
          eleJumpsAndKinks_.insert(make_pair(ele->Id(),eleJumpAndKinks));
        } // end if element bisected
        else // in touchedminus ele just compute jump since kink value = 0 (singulary sysmat)
        {
          // compute the jump height at the interface in this element as a mean value
          double sysmat = 0;
          double rhsvel = 0;
          double rhspres = 0;
//          cout << "enrvalues of touched ele are " << enrValues[2];
          double velJump = 0;
          double presJump = 0;
          
          LINALG::Matrix<2,4> currentJumpsAndKinks(true);
          vector<LINALG::Matrix<2,4> > eleJumpAndKinks;
          for (size_t field=0;field<oldVectors_.size();field++)
          {
            currentJumpsAndKinks.Clear();
            
            sysmat = 0.5*numnodeused;
            rhsvel = 0;
            rhspres = 0;
            
            for (int index=0;index<numnode;index++)
            {
              if (dmin(index)!=INFINITY) // if infinity then dont used node
              {
                LINALG::Matrix<nsd,1> currVelEnrValues;
                for (int i=0;i<nsd;i++)
                  currVelEnrValues(i) = enrValues[field](i,index);
                
                double currVelEnrScalar = currVelEnrValues.Dot(normal);
                rhsvel -= currVelEnrScalar;
                
                rhspres -= enrValues[field](4,index);
              }
            }
//            cout << "rhs in touched ele is " << rhs << endl;
            
            velJump = rhsvel/sysmat;
            presJump = rhspres/sysmat;
            
//            cout << "curr jump in touched ele is " << currJump << endl;
            currentJumpsAndKinks(0,0) = velJump;
            currentJumpsAndKinks(0,1) = presJump;
            
            eleJumpAndKinks.push_back(currentJumpsAndKinks);
          } // end loop over vectorfield size
          cout << "final old jumps and kinks at touched element are " << eleJumpAndKinks[2] << endl;
          eleJumpsAndKinks_.insert(make_pair(ele->Id(),eleJumpAndKinks));
        } // end if element touchedminus
#endif
      } // end if all nodes of element on proc
    } // end if element intersected
  } // end loop over col ele nodes
} // end function jumpKinkHeight



void XFEM::Enrichmentvalues::computeNewEnrichments(
)
{
  const int nsd = 3;
  for (int inode=0;inode<discret_->NumMyRowNodes();inode++)
  {
    const DRT::Node* currnode = discret_->lRowNode(inode); // current node
    
    if (newEnrValueNeeded(currnode))
    {
      const DRT::Element* const* eles = currnode->Elements();
      const int numeles = currnode->NumElement();
      vector<LINALG::Matrix<2,nsd+1> > currJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true));
      vector<LINALG::Matrix<2,nsd+1> > averageJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true));
      
      // compute average jump and kink value around node if elements are enriched
      int numOldIntersectedEle=0; // index for number of bisected elements
      for (int iele=0;iele<numeles;iele++)
      {
        const int elegid = eles[iele]->Id(); // global id of current element
        
        if (ih_old_->ElementBisected(eles[iele]) || ih_old_->ElementTouchedMinus(eles[iele])) // element intersected
        {
          numOldIntersectedEle++;
          currJumpsAndKinks = eleJumpsAndKinks_.find(elegid)->second;
//          cout << "for ele " << *eles[iele] << "currjumpsandkinks are " << currJumpsAndKinks[2];
          for (size_t ivector=0;ivector<newVectors_.size();ivector++)
            averageJumpsAndKinks[ivector] += currJumpsAndKinks[ivector];
        } // end if element bisected
      } // end loop over elements containing the node
//      cout << "summed jumps and kinks are " << averageJumpsAndKinks[2] << endl;
      
      vector<LINALG::Matrix<1,nsd+1> > finalEnrichmentValues(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true));
      if (numOldIntersectedEle > 0) // >=1 elements around the node where enriched -> computed value can be used
      {
        for (size_t ivector=0;ivector<newVectors_.size();ivector++)
          averageJumpsAndKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
//        cout << "average jumps and kinks are " << averageJumpsAndKinks[2] << endl;
        
        int numNewIntersectedEle = 0;
        
        for (int iele=0;iele<numeles;iele++)
        {
          if (ih_->ElementBisected(eles[iele]) || ih_->ElementTouchedMinus(eles[iele])) // element intersected
          {
            numNewIntersectedEle++;
            double dist = 0.0;
            LINALG::Matrix<nsd,1> normal(true);
            SignedDistance(currnode,eles[iele]->Id(),ih_,dist,normal);
            
            for (size_t ivector=0;ivector<newVectors_.size();ivector++)
            {
#ifndef ENR_VEL_JUMP_SCALAR
              for (int entry=0;entry<nsd+1;entry++)
                finalEnrichmentValues[ivector](entry) += 
                   -0.5*averageJumpsAndKinks[ivector](0,entry); // -0.5*jump
//                   +0.5*dist*averageJumpsAndKinks[ivector](1,entry); // +0.5*kink*signeddistance
#else
              // velocity entries
              for (int entry=0;entry<nsd;entry++)
                finalEnrichmentValues[ivector](entry) += 
                   -normal(entry)*(-0.5)*averageJumpsAndKinks[ivector](0,0); // -0.5*jump
//                   +0.5*dist*averageJumpsAndKinks[ivector](1,0)); // +0.5*kink*signeddistance
              
              // pressure entry
              finalEnrichmentValues[ivector](nsd) += 
                 -0.5*averageJumpsAndKinks[ivector](0,1); // -0.5*jump
//                 +0.5*dist*averageJumpsAndKinks[ivector](1,1); // +0.5*kink*signeddistance
#endif
            } // end loop over vector size
          } // end if element bisected
        } // end loop over elements containing the node
        for (size_t ivector=0;ivector<newVectors_.size();ivector++)
          finalEnrichmentValues[ivector].Scale(1.0/static_cast<double>(numNewIntersectedEle));
        cout << *currnode << " has final enr values in std case are " << finalEnrichmentValues[2] << endl;
        
        
        int i=0; // index which entry has to be used
        // set nodal velocities and pressures with help of the field set of node
        const std::set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(currnode->Id()));
        for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
            fieldenr != fieldenrset.end();++fieldenr)
        {
          if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
          {
            const DofKey<onNode> newdofkey(currnode->Id(), *fieldenr);
            const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
            const int lid = newdofrowmap_.LID(newdofpos);
            
            for (size_t index=0;index<newVectors_.size();index++)
              (*newVectors_[index])[lid] = finalEnrichmentValues[index](0,i);
            i++;
          } // end if enrichment type not standard
        } // end loop over fieldenr
      }
      else
      {
        failedData failed(LINALG::Matrix<nsd,1>(currnode->X()),INFINITY,vector<LINALG::Matrix<2,4> >(newVectors_.size(),LINALG::Matrix<2,4>(true)));
        failed_.insert(make_pair(currnode->Id(),failed));
      }
    } // end if node is enriched
  } // end loop over row nodes
}



void XFEM::Enrichmentvalues::handleFailedNodes(
)
{
  const int nsd = 3;
  
  // evaluate all enriched nodes at old interface position
  set<int> oldEnrNodes;
  for (int inode=0;inode<discret_->NumMyRowNodes();inode++)
  {
    const int nodegid = discret_->lRowNode(inode)->Id(); // current node
    
    // set nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldenrset(olddofman_->getNodeDofSet(nodegid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
        oldEnrNodes.insert(nodegid);
        break;
      }
    }
  }
  
  
#ifdef PARALLEL
  // get the nearest enriched node at old timestep
  for (int iproc=0;iproc<numproc_;iproc++)
  {
#endif // PARALLEL
    for (map<int,failedData>::iterator newnode=failed_.begin();
        newnode!=failed_.end(); newnode++)
    {
      failedData& failed = newnode->second;
      
      for (set<int>::const_iterator enrnode=oldEnrNodes.begin();
          enrnode!=oldEnrNodes.end(); enrnode++)
      {
        LINALG::Matrix<nsd,1> oldCoords(discret_->gNode(*enrnode)->X());
        LINALG::Matrix<nsd,1> diff(true);
        diff.Update(1.0,failed.coords_,-1.0,oldCoords);
        
        if (diff.Norm2() < failed.dist_) // new nearest node
        {
          failed.dist_ = diff.Norm2();
          
          DRT::Node* currOldNode = discret_->gNode(*enrnode);
          const DRT::Element* const* eles = currOldNode->Elements();
          const int numeles = currOldNode->NumElement();
          
          vector<LINALG::Matrix<2,nsd+1> > finalJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true));
          vector<LINALG::Matrix<2,nsd+1> > currJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true));
          
          int numOldIntersectedEle=0; // index how much elements are intersected
          for (int iele=0;iele<numeles;iele++)
          {
            const int elegid = eles[iele]->Id(); // global id of current element
            
            if (ih_old_->ElementBisected(eles[iele]) || ih_old_->ElementTouchedMinus(eles[iele])) // element intersected
            {
              numOldIntersectedEle++;
              currJumpsAndKinks=eleJumpsAndKinks_.find(elegid)->second;
              for(size_t ivector=0;ivector<oldVectors_.size();ivector++)
              {
                finalJumpsAndKinks[ivector] += currJumpsAndKinks[ivector];
              }
            } // end if element bisected
          } // end loop over elements around the node
          
          if (numOldIntersectedEle==0)
          {
            cout << *currOldNode << " is enriched at old timestep" << endl;
            dserror("Enriched node shall have an intersected element around");
          }
          
          for(size_t ivector=0;ivector<oldVectors_.size();ivector++)
            finalJumpsAndKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
          failed.enrValues_ = finalJumpsAndKinks;
        } // end if new nearest node
      } // loop over enriched nodes on current proc at old interface
    } // loop over failed nodes
    
#ifdef PARALLEL
    exportEnrichmentData();
#endif // PARALLEL
    
  } // loop over processors
  
  
//  for (map<int,failedData>::iterator i=failed_.begin();
//      i!=failed_.end();i++)
//  {
//    failedData failed = i->second;
//    cout << "nodegid is " << i->first << ", coords are " << failed.coords_ << ", dist is "
//        << failed.dist_ << " and enrvalues are " << failed.enrValues_[2] << endl;
//  }
  
  
  // compute new enrichment values
  for (map<int,failedData>::iterator newnode=failed_.begin();
      newnode!=failed_.end(); newnode++)
  {
    const DRT::Node* currnode = discret_->gNode(newnode->first);
    
    vector<LINALG::Matrix<2,nsd+1> > currJumpsAndKinks = newnode->second.enrValues_;
    vector<LINALG::Matrix<1,nsd+1> > finalEnrichmentValues(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true));
    int numNewIntersectedEle=0; // index for number of intersected elements
    
    
    const DRT::Element* const* eles = currnode->Elements();
    const int numele = currnode->NumElement();
//    cout << "current jumps and kinks " << currJumpsAndKinks[2] << endl;
    
    for (int iele=0;iele<numele;iele++)
    {
      if (ih_->ElementBisected(eles[iele]) || ih_->ElementTouchedMinus(eles[iele]))
      {
        numNewIntersectedEle++;
        double dist = 0;
        LINALG::Matrix<nsd,1> normal(true);
        SignedDistance(currnode,eles[iele]->Id(),ih_,dist,normal);
        
        for (size_t ivector=0;ivector<newVectors_.size();ivector++)
        {
#ifndef ENR_VEL_JUMP_SCALAR
          for (int entry=0;entry<nsd+1;entry++)
            finalEnrichmentValues[ivector](entry) += 
               -0.5*currJumpsAndKinks[ivector](0,entry); // -0.5*jump
//               +0.5*dist*currJumpsAndKinks[ivector](1,entry); // +0.5*kink*signeddistance
#else
              // velocity entries
              for (int entry=0;entry<nsd;entry++)
                finalEnrichmentValues[ivector](entry) += 
                   -normal(entry)*(-0.5)*currJumpsAndKinks[ivector](0,0); // -0.5*jump
//                   +0.5*dist*currJumpsAndKinks[ivector](1,0)); // +0.5*kink*signeddistance
              
              // pressure entry
              finalEnrichmentValues[ivector](nsd) += 
                 -0.5*currJumpsAndKinks[ivector](0,1); // -0.5*jump
//                 +0.5*dist*currJumpsAndKinks[ivector](1,1); // +0.5*kink*signeddistance
#endif
        }
      } // end loop over vector size
    } // end loop over enriched elements around a node
    
    for (size_t ivector=0;ivector<newVectors_.size();ivector++)
    {
      finalEnrichmentValues[ivector].Scale(1.0/static_cast<double>(numNewIntersectedEle));
    }
    cout << *currnode << " has final enrichment values in alternative case are " << finalEnrichmentValues[2] << endl;
    
    int i=0; // index which entry has to be used
    const std::set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(currnode->Id()));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
      {
        const DofKey<onNode> newdofkey(newnode->first, *fieldenr);
        const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
        const int lid = newdofrowmap_.LID(newdofpos);
        
        for (size_t index=0;index<newVectors_.size();index++)
          (*newVectors_[index])[lid] = finalEnrichmentValues[index](0,i);
        i++;
      } // end if enrichment type not standard
    } // end loop over fieldenr
  } // end loop over failed nodes
}



void XFEM::Enrichmentvalues::getDataOfIntersectedEle(
    const DRT::Element* element,
    LINALG::Matrix<1,8>& dmin,
    bool& allEleNodesOnProc,
    vector<LINALG::Matrix<4,8> >& enrValues,
    LINALG::Matrix<3,1>& normal
)
{
  const int numnode = 8;
  const DRT::Node* const* elenodes = element->Nodes();
  // clear content of enrValues, but hold its size
  for (size_t ivector=0;ivector<enrValues.size();ivector++)
    enrValues[ivector].Clear();
  
  for (int inode=0;inode<numnode;inode++) // loop over element nodes
  {
    const DRT::Node* currnode = elenodes[inode]; // current node
    
    // evaluate nodal velocity and pressure enrichment values
    const std::set<XFEM::FieldEnr>& fieldenrset(olddofman_->getNodeDofSet(currnode->Id()));
    
    if (fieldenrset.empty()) // node data not on processor
    {
      allEleNodesOnProc = false;
      break; // col element not needed, see remark above
    }
    else
    {
      // evaluate the minimal distance of the nodes to the interface segments in this element
      // this part is adopted from src/drt_combust/combust_reinitializer.cpp
      // more information can be found there in function 'signedDistanceFunction'
      SignedDistance(currnode,element->Id(),ih_old_,dmin(inode),normal);
      
      int i=0; // index which entry has to be set
      for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
          fieldenr != fieldenrset.end();++fieldenr)
      {
        const DofKey<onNode> olddofkey(currnode->Id(), *fieldenr);
        const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
        
        // remark: in touched elements nodes may be not enriched
        // they will not give entries here and therefore the enrvalues
        // stay zero in some entries what is ok!
        if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
        {
          const int lid = olddofcolmap_.LID(olddofpos);
          
          for (size_t index=0;index<oldVectors_.size();index++)
            enrValues[index](i,inode) = (*oldVectors_[index])[lid];
          i++;
        }
      } // end loop over fieldenr
    } // end if fieldenrset is empty
  } // end loop over nodes
//  cout << "for element " << *element << " enrvalues are " << enrValues[2] << endl;
}



void XFEM::Enrichmentvalues::analyseEnrichments(
    const DRT::Element* element,
    LINALG::Matrix<1,8>& dmin,
    int& numnodeused
)
{
  const int numnode = 8;
  const int nsd = 3;
  
  // remark: this function shall handle critical enrichments
  // therefore it shall sort out enrichments nearly as far away
  // from the elements interface as the element diameter is.
  // But still enrichment values for both element sides shall be
  // available for computation as far as we have no touched element
  
  // compute element diameter
  double diameter = -1;
  const DRT::Node* const* nodes = element->Nodes();
  for (int inode=0;inode<numnode;inode++)
  {
    LINALG::Matrix<nsd,1> icoords(nodes[inode]->X());
    
    for (int jnode=0;jnode<numnode;jnode++)
    {
      LINALG::Matrix<nsd,1> jcoords(nodes[jnode]->X());
      jcoords.Update(1.0,icoords,-1.0);
      
      if (jcoords.Norm2()>diameter) diameter = jcoords.Norm2();
    } // end loop over element nodes
  } // end loop over element nodes
  if (diameter < 0)
    dserror("element diameter shall be greater than zero");
  
  // handle touched elements (the enrichment value influences the touchedminuselement)
  if (ih_old_->ElementTouchedMinus(element))
  {
    for (int inode=0;inode<numnode;inode++)
    {
      if (fabs(dmin(inode)/diameter) > 1e-6) // no touching node is not enriched by this element and shall therefore not be used
      {
        numnodeused--;
        dmin(inode) = INFINITY;
      }
    }
  }
  
  // handle bisected elements
  else if (ih_old_->ElementBisected(element))
  {
    // sort out nodes far away since their enrichment values are critical
    // but guarantee that both interface sides still give entries
    LINALG::Matrix<1,numnode> tmpdmin = dmin;
    double Tol = 0.5;
    numnodeused = numnode;
    bool plusside = false;
    bool minusside = false;
    
    while(true)
    {
      for (int inode=0;inode<numnode;inode++)
      {
        if ((diameter-fabs(tmpdmin(inode)))/diameter < Tol) // node shall not be used
        {
          numnodeused--;
          tmpdmin(inode) = INFINITY;
        }
        else
        {
          if (tmpdmin(inode)>=0) plusside = true;
          else minusside = true;
        }
      } // end loop over nodes
      
      // check if enrichment values for both sides exist
      if (plusside && minusside) // both sides exist
      {
//        cout << "final tolerance is " << Tol << endl;
        break;
      }
      else
      {
        // reset values
        tmpdmin = dmin;
        Tol = 0.8*Tol; // get tol smaller
        numnodeused = numnode;
        plusside = false;
        minusside = false;
      }
    } // end while loop
    
    dmin = tmpdmin;
  }
  else
    dserror("BUG! This function shall just handle bisected or touched minus elements!");
}



bool XFEM::Enrichmentvalues::newEnrValueNeeded(
    const DRT::Node* node
)
{
  // case 1: the node was not enriched in old timestep and therefore needs a new enrichment
  const int gid = node->Id();
//  cout << "here with node " << *node << endl;
  
  const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gid));
  for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
      fieldenr != fieldenrset.end();++fieldenr)
  {
    const DofKey<onNode> newdofkey(gid, *fieldenr);
    
    if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
    {
#ifdef ENR_FULL
      return true;
#endif
      map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofColDistrib_.find(newdofkey);
      
      // olddof not found -> enrichment did not exist before and has to be set
      if (olddof == oldNodalDofColDistrib_.end())
        return true;
    }
  }
  
  
#ifndef MIN_ENR_VALUES
#ifndef ENR_FULL
  // case 2: the node was enriched but the according support is very small.
  // Then the enrichment value is potentially much too high and therefore needs a new value.
  bool critCut = false; // true if intersected ele around and all intersected eles have small support for enr shape fcn
  const DRT::Element* const* eles = node->Elements();
  
  // determine interface side of node
  bool domainPlus;
  if ((*phiOldCol_)[node->LID()] >= 0)
    domainPlus = true;
  else
    domainPlus = false;
  
  // loop over elements around node
  for (int iele=0;iele<node->NumElement();iele++)
  {
    if (ih_old_->ElementBisected(eles[iele])) // element intersected
    {
      if (domainPlus) // when node is in plus domain, support of enrichment is the minus part of the element
      {
        set<int>::const_iterator tmp = critElesMinus_.find(eles[iele]->Id());
      
        if (tmp == critElesMinus_.end()) // volumes just saved in critical cases
        {
          critCut = false; // one element has a big enough support for the node -> no problem with values
          break;
        }
        else
          critCut = true;
      } // end if node in plus domain
      else
      {
        set<int>::const_iterator tmp = critElesPlus_.find(eles[iele]->Id());
      
        if (tmp == critElesPlus_.end()) // volumes just saved in critical cases
        {
          critCut = false; // one element has a big enough support for the node -> no problem with values
          break;
        }
        else
          critCut = true;
      } // end if node in minus domain
    } // end if element bisected or touched
  } // end loop over elements around node
//  cout << "bool is " << critCut << endl;
  return critCut;
#endif
#endif
  return false; // if true shall be returned it was done before
}



#ifndef ENR_FULL
#ifndef MIN_ENR_VALUES
void XFEM::Enrichmentvalues::computeDomainCellVolumes(
)
{
  double plusVol;
  double minusVol;
  double currVol;
  double eleVol;
  
  DRT::Element* currEle = NULL;
  
  for (int iele=0;iele<discret_->NumMyColElements();iele++)
  {
    currEle = discret_->lColElement(iele);
    
    if (ih_old_->ElementBisected(currEle))
    {
      plusVol = 0.0;
      minusVol = 0.0;
      
      // get domain integration cells for this element
      const GEO::DomainIntCells&  domainIntCells(ih_old_->GetDomainIntCells(currEle));
      
      // loop over domain integration cells
      for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
      {
        currVol = cell->VolumeInPhysicalDomain();
        
        if (cell->getDomainPlus()) plusVol += currVol;
        else minusVol += currVol;
      } // end loop over domain integration cells
      
      if (plusVol == 0.0 || minusVol == 0.0)
      {
        cout << "element " << *currEle << " shall be bisected" << endl;
        cout << "element volume in plus domain is " << plusVol << endl;
        cout << "element volume in minus domain is " << minusVol << endl;
        cout << "WARNING!!! Bisected domain shall have integration cells on both interface sides!";
      }
      
      eleVol = plusVol + minusVol;
      
      if (plusVol/eleVol < critTol_)
        critElesPlus_.insert(currEle->Id());
      
      if (minusVol/eleVol < critTol_)
        critElesMinus_.insert(currEle->Id());
    } // end if element bisected
  } // end loop over col elements
}
#endif
#endif



void XFEM::Enrichmentvalues::SignedDistance(
    const DRT::Node* node,
    const int elegid,
    const RCP<InterfaceHandle> curr_ih,
    double& dist,
    LINALG::Matrix<3,1>& normal
)
{
  const int nsd = 3;
  LINALG::Matrix<nsd,1> nodecoord(node->X()); // coordinates of this node
  
  //-----------------------------------------------------------
  // compute smallest distance to the flame front for this node
  //-----------------------------------------------------------
  // smallest distance to the vertex of a flame front patch
  double vertexdist = 7777.7; // default value
  // smallest distance to flame front
  double mindist = 5555.5; // default value
  
  // number of flamefront patches for this element
  const std::vector<GEO::BoundaryIntCell> patches = curr_ih->GetBoundaryIntCells(elegid);
  const int numpatch = patches.size();
  
  // loop flame front patches of this element
  for(int ipatch=0; ipatch<numpatch; ++ipatch)
  {
    // get a single patch from group of flamefront patches
    const GEO::BoundaryIntCell patch = patches[ipatch];
    // only triangles and quadrangles are allowed as flame front patches (boundary cells)
    if (!(patch.Shape() == DRT::Element::tri3 or
          patch.Shape() == DRT::Element::quad4))
      dserror("invalid type of boundary integration cell for reinitialization");
    
    // get coordinates of vertices defining flame front patch
    const LINALG::SerialDenseMatrix& patchcoord = patch.CellNodalPosXYZ();
    
    // compute normal vector to flame front patch
    normal.Clear();
    ComputeNormalVectorToFlameFront(patch,patchcoord,normal);
    
    //-----------------------------------------
    // find flame front patches facing the node
    //-----------------------------------------
    // boolean indicating if facing patch was found
    bool facenode = false;
    // distance to the facing patch
    double patchdist = 7777.7; // default value
    // check if this patch faces the node
    FindFacingPatchProjCellSpace(nodecoord,patch,patchcoord,normal,facenode,patchdist);
    
    // a facing patch was found
    if (facenode == true)
    {
      // overwrite smallest distance if computed patch distance is smaller
      if (fabs(patchdist) < fabs(mindist))
        mindist = patchdist;
    }
    
    // compute smallest distance to vertices of this flame front patch
    ComputeDistanceToPatch(nodecoord,patch,patchcoord,normal,vertexdist);
  }
  
  if (fabs(vertexdist) < fabs(mindist))
  {
    // if the sign has been changed by mistake in ComputeDistanceToPatch(), this has to be corrected here
    if ((*phiOldCol_)[node->LID()] * vertexdist < 0.0 )
      mindist = -vertexdist;
    else
      mindist = vertexdist;
  }
  if (mindist == 5555.5) // case 3
  {
    dserror ("computation of minimal distance failed");
    return;
  }
  else
  {
    dist = mindist;
    return;
  }
}



/*------------------------------------------------------------------------------------------------*
 | private: find a facing flame front patch by projecton of node into boundary cell space         |
 |                                                                                    henke 12/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::Enrichmentvalues::FindFacingPatchProjCellSpace(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    bool&                            facenode,
    double&                          patchdist)
{
  // indicator
  facenode = false;

  static LINALG::Matrix<2,1> eta(true);
  double alpha = 0.0;

  //-------------------------------------------------------
  // perform Newton-Raphson method to project node on patch
  //-------------------------------------------------------
  bool converged = false;
  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    converged = ProjectNodeOnPatch<DRT::Element::tri3>(node, patch, patchcoord, normal, eta, alpha);
    break;
  }
  case DRT::Element::quad4:
  {
    converged = ProjectNodeOnPatch<DRT::Element::quad4>(node, patch, patchcoord, normal, eta, alpha);
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
  }

  // Newton iteration converged
  //  cout << "Newton iteration converged in " << iter << " steps!" << endl;

  //----------------------------------------------------
  // check if projection lies within boundary cell space
  //----------------------------------------------------
  // remark: - tolerance has to be of same order as the tolerance that coordinates of projected nodes
  //           differ from an exact position on edges of patches (e.g. 1.0E-7 ~ 1.0E-8 -> 1.0E-6)
  //         - if this is not the case, the level set function can become tilted, since valid
  //           patches are ignored
  double TOL= 1e-6;
  
  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    // criteria for tri3 patch
    if ((eta(0) > -TOL) and (eta(0) < 1.0+TOL) and
        (eta(1) > -TOL) and (eta(1) < 1.0+TOL) and
        (1.0-eta(0)-eta(1) > -TOL) and (1.0-eta(0)-eta(1) < 1.0+TOL) and
        converged)
    {
      facenode = true;
      patchdist = alpha;
//      cout << "facing patch found (tri3 patch)! coordinates eta(0): " << eta(0) << " eta(1) " << eta(1) << endl;
    }
    break;
  }
  case DRT::Element::quad4:
  {
    // criteria for quad4 patch
    if ((eta(0) > -1.0-TOL) and (eta(0) < 1.0+TOL) and
        (eta(1) > -1.0-TOL) and (eta(1) < 1.0+TOL) and
        converged)
    {
      facenode = true;
      patchdist = alpha;
//      cout << "facing patch found (quad4 patch)!" << endl;
    }
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
  }
//  if (!converged)
//  {
//    cout << "node x component " << node(0,0) << endl;
//    cout << "node y component " << node(1,0) << endl;
//    cout << "node z component " << node(2,0) << endl;
//    cout << "eta1 " << eta(0) << endl;
//    cout << "eta2 " << eta(1) << endl;
//    cout << "alpha " << alpha << endl;
//    cout << "patch vertices x component " << patchcoord(0,0) << " " << patchcoord(0,1) << " " << patchcoord(0,2) << endl;
//    cout << "patch vertices y component " << patchcoord(1,0) << " " << patchcoord(1,1) << " " << patchcoord(1,2) << endl;
//    cout << "patch vertices z component " << patchcoord(2,0) << " " << patchcoord(2,1) << " " << patchcoord(2,2) << endl;
//  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | private: compute distance to vertex of patch                                       henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::Enrichmentvalues::ComputeDistanceToPatch(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    double&                          vertexdist
)
{
  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch
  static LINALG::Matrix<3,1> vertex(true);
  // distance vector from patch to node
  static LINALG::Matrix<3,1> dist(true);

  // compute distance to all vertices of patch
  for(size_t ivert = 0; ivert<numvertices; ++ivert)
  {
    // vertex of flame front patch
    vertex(0) = patchcoord(0,ivert);
    vertex(1) = patchcoord(1,ivert);
    vertex(2) = patchcoord(2,ivert);

    // compute distance vector from flame front to node
    dist.Update(1.0, node, -1.0, vertex);

    // compute L2-norm of distance vector
    double normdist = sqrt(dist(0)*dist(0) + dist(1)*dist(1) + dist(2)*dist(2));
    if(normdist < fabs(vertexdist))
    {
      // determine sign of distance vector from node to flame front
      double tmp = normal(0)*dist(0) + normal(1)*dist(1) + normal(2)*dist(2);
      // tmp < 0 if node in burnt domain (G>0) and vice versa
      if (tmp <= 0.0) // 'normal' and 'dist' point in different directions
        vertexdist = normdist;
      else // 'normal' and 'dist' point in the same direction
        vertexdist = -normdist;
      //cout << "distance to vertex is overwritten by: " << vertexdist << endl;
    }
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 | private: compute normal vector to flame front patch                                henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void XFEM::Enrichmentvalues::ComputeNormalVectorToFlameFront(
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      LINALG::Matrix<3,1>&             normal)
{
  // first point of flame front patch
  LINALG::Matrix<3,1> point1;
  point1(0) = patchcoord(0,0);
  point1(1) = patchcoord(1,0);
  point1(2) = patchcoord(2,0);

  // second point of flame front patch
  LINALG::Matrix<3,1> point2;
  point2(0) = patchcoord(0,1);
  point2(1) = patchcoord(1,1);
  point2(2) = patchcoord(2,1);

  // first edge of flame front patch
  LINALG::Matrix<3,1> edge1;
  edge1.Update(1.0, point2, -1.0, point1);

  // third point of flame front patch
  point2(0) = patchcoord(0,2);
  point2(1) = patchcoord(1,2);
  point2(2) = patchcoord(2,2);

  // second edge of flame front patch (if patch is triangle; if not: edge 2 is secant of polygon)
  LINALG::Matrix<3,1> edge2;
  edge2.Update(1.0, point2, -1.0, point1);

  // compute normal vector of patch (cross product: edge1 x edge2)
  // remark: normal vector points into unburnt domain (G<0)
  normal(0) = (edge1(1)*edge2(2) - edge1(2)*edge2(1));
  normal(1) = (edge1(2)*edge2(0) - edge1(0)*edge2(2));
  normal(2) = (edge1(0)*edge2(1) - edge1(1)*edge2(0));
  // TODO remove restriction to 2D
  // schott
  normal(2) = 0.0;

//  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
//  cout << "proc " << comm.MyPID() << " normal " <<  normal << endl;
//  cout << "proc " << comm.MyPID() << " patch " <<  patchcoord << endl;

  // compute unit (normed) normal vector
  double norm = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
  if (norm == 0.0) dserror("norm of normal vector is zero!");
  normal.Scale(1.0/norm);

#ifdef DEBUG
  if (!((normal(2) > 0.0-1.0E-8) and (normal(2) < 0.0+1.0E-8)))
  {
    cout << "z-component of normal: " << normal(2) << endl;
    dserror ("pseudo-3D problem not symmetric anymore!");
  }
#endif

  return;
}



/// project node into the boundary cell space (2D)
template<DRT::Element::DiscretizationType DISTYPE>
bool XFEM::Enrichmentvalues::ProjectNodeOnPatch(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    LINALG::Matrix<2,1>&             eta,
    double&                          alpha)
{
  // indicator for convergence of Newton-Raphson scheme
  bool converged = false;
  // number space dimensions for 3d combustion problems
  const size_t nsd = 3;
  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // get coordinates of vertices of flame front patch
  // remark: here we only get a view (bool true) on the SerialDenseMatrix returned by CellNodalPosXYZ()
  LINALG::Matrix<nsd,numvertices> patchcoordfix(patchcoord.A(),true);

  static LINALG::Matrix<numvertices,1> funct(true);
  static LINALG::Matrix<2,numvertices> deriv(true);
  static LINALG::Matrix<nsd,1> projX(true);
  static LINALG::Matrix<nsd,2> gradprojX(true);

  //----------------------------------
  // start values for iterative scheme
  //----------------------------------
  // start position (barycenter of triangular boundary cell)
  eta(0) = 1.0/3.0;
  eta(1) = 1.0/3.0;
  // auxiliary variable
  // remark: third unknown to close system of equations; arbitrary value
  alpha = 0.0;

  // function F (system of equations)
  static LINALG::Matrix<nsd,1> f(true);
  // gradient of function F (dF/deta(0), dF/deta(1), dF/dalpha)
  static LINALG::Matrix<nsd,nsd> gradf(true);
  // increment in Newton iteration (unknown to be solved for)
  static LINALG::Matrix<nsd,1> incr(true);

  // maximum number Newton iterations
  size_t maxiter = 3;
  // convergence tolerance
  double conv = 0.0;

  //------------------------------------------------------
  // Newton-Raphson loop for non-linear projection problem
  //------------------------------------------------------
  for (size_t iter=0;iter<maxiter;++iter)
  {
    // evaluate shape functions in boundary cell space at current position \eta_1,\eta_2 on the patch
    funct.Clear();
    DRT::UTILS::shape_function_2D(funct,eta(0),eta(1),patch.Shape());
    // evaluate derivatives of shape functions in boundary cell space at current position \eta_1,\eta_2 on the patch
    deriv.Clear();
    DRT::UTILS::shape_function_2D_deriv1(deriv,eta(0),eta(1),patch.Shape());

    // evaluate projection X of node P at current position \eta_1,\eta_2 on the patch
    // projX(i,j) = patchcoord(i,k)*funct(k,1)
    projX.Clear();
    projX.MultiplyNN(patchcoordfix,funct);

    // evaluate gradient of projection X of node P at current position \eta_1,\eta_2 on the patch
    // gradprojX(i,j) = patchcoord(i,k)*deriv(j,k)
    gradprojX.Clear();
    gradprojX.MultiplyNT(patchcoordfix,deriv);

    //---------------------------------------------------
    // build system of equations F and its gradient gradF
    //---------------------------------------------------
    f.Clear();
    gradf.Clear();
    incr.Clear();
    for (size_t icoord=0;icoord<nsd;++icoord)
    {
      // evaluate function f
      f(icoord) = projX(icoord) + alpha*normal(icoord) - node(icoord);
      // evaluate gradient of function at current position on patch in boundary cell space
      gradf(icoord,0) = gradprojX(icoord,0);
      gradf(icoord,1) = gradprojX(icoord,1);
      gradf(icoord,2) = normal(icoord);
    }

    // check convergence
    conv = sqrt(f(0)*f(0)+f(1)*f(1)+f(2)*f(2));
    //cout << "iteration " << iter << ": -> |f|=" << conv << endl;
    if (conv <= 1.0E-12) break;

    //----------------------------------------------------
    // solve linear system of equations: gradF * incr = -F
    //----------------------------------------------------
    // F = F*-1.0
    f.Scale(-1.0);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<nsd,nsd,1> solver;
    solver.SetMatrix(gradf);              // set A=gradF
    solver.SetVectors(incr, f);           // set X=incr, B=F
    solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();           // ?
    int err = solver.Solve();             // incr = gradF^-1.F
    if ((err != 0) || (err2!=0))
      dserror("solving linear system in Newton-Raphson method for projection failed");

    // update eta and alpha
    eta(0) += incr(0);
    eta(1) += incr(1);
    alpha  += incr(2);
    //cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << " alpha: " << alpha << endl;
  }
  // change sign to preserve sign of G-function
  alpha = -alpha;

  // Newton iteration unconverged
  if (conv > 1.0E-12)
  {
    alpha = 7777.7;
//        cout << "projection did not converge" << endl;
    //dserror("projection did not converge!");
  }
  else
  {
    converged = true;
    //cout << "convergence criterion " << conv << endl;
    //cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << " alpha: " << alpha << endl;
  }

  return converged;
}




#ifdef PARALLEL

/*------------------------------------------------------------------------------------------------*
 * export enrichment data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Enrichmentvalues::exportEnrichmentData()
{
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  
  // destination proc (the "next" one)
  int dest = myrank_+1;
  if(myrank_ == (numproc_-1))
    dest = 0;
  
  // source proc (the "last" one)
  int source = myrank_-1;
  if(myrank_ == 0)
    source = numproc_-1;
  
  // vector including all data that has to be send to the next proc
  vector<char> dataSend;
  
  // packing the data
  DRT::ParObject::AddtoPack(dataSend,failed_.size());
  for (map<int,failedData>::iterator newnode=failed_.begin();
      newnode!=failed_.end(); newnode++)
  {
    failedData& failed(newnode->second);
    
    DRT::ParObject::AddtoPack(dataSend,newnode->first);
    DRT::ParObject::AddtoPack(dataSend,failed.coords_);
    DRT::ParObject::AddtoPack(dataSend,failed.dist_);
    DRT::ParObject::AddtoPack(dataSend,failed.enrValues_);
  }
  
  vector<int> lengthSend(1,0);
  lengthSend[0] = dataSend.size();
  int size_one = 1;
  
#ifdef DEBUG
  cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
#endif
  
  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter_.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
  cout << "here1" << endl;// ... and receive length
  vector<int> lengthRecv(1,0);
  exporter_.Receive(source, length_tag, lengthRecv, size_one);
  exporter_.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);

  discret_->Comm().Barrier();

  // ... and receive data
  vector<char> dataRecv(lengthRecv[0]);
  exporter_.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
  exporter_.Wait(req_data);
  
#ifdef DEBUG
  cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << endl;
#endif
  
  
  // pointer to current position of group of cells in global string (counts bytes)
  size_t posinData = 0;
  
  // initialize temporary vectors that should be filled
  size_t numberOfNodes = 0;
  
  // clear vector that should be filled
  failed_.clear();
  
  // unpack received data
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
  
  for (size_t inode=0;inode<numberOfNodes;inode++)
  {
    int gid;
    LINALG::Matrix<nsd,1> coords;
    double dist;
    vector<LINALG::Matrix<2,4> > enrValues;

    DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,dist);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,enrValues);
    
    failedData failed(coords,dist,enrValues);
    failed_.insert(make_pair(gid,failed));
  }
} // end exportEnrichmentData

#endif // parallel

#endif // CCADISCRET

