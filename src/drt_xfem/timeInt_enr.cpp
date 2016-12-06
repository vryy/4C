/*!-----------------------------------------------------------------------------------------------*
\file timeInt_enr.cpp

\brief provides the enrichment computation class(es)

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#include "timeInt_enr.H"


/*------------------------------------------------------------------------------------------------*
 * constructor of the enrichment recomputation with projection approach          winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
XFEM::EnrichmentProjection::EnrichmentProjection(
    XFEM::TIMEINT& timeInt,
    const double& veljump,
    INPAR::COMBUST::XFEMTimeIntegrationEnr timeIntEnr,
    INPAR::COMBUST::XFEMTimeIntegrationEnrComp timeIntEnrType
) : ENR(timeInt,timeIntEnr,timeIntEnrType),
veljump_(veljump)
{
#ifndef COMBUST_SETJUMP
  callOldValues(); // compute the old enrichment data
#endif
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * call computation of the projection approach for enrichment values             winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::compute(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
)
{
  if (FGIType_==FRSNot1_)
    return;
  else if (FGIType_==FRS1FGINot1_)
    timeIntData_->clear();

  handleVectors(newRowVectorsn,newRowVectorsnp);

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



/*------------------------------------------------------------------------------------------------*
 * call computation of old values at old interface position                      winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::callOldValues(
)
{
  // remark: Finally values shall be computed for row nodes.
  // Therefore we loop here over col elements. If nodes of this col
  // element are not in the current procs col node map, the according
  // col element is not required finally (test within this function)
  for (int iele=0;iele<discret_->NumMyColElements();iele++)
  {
    const DRT::Element* ele = discret_->lColElement(iele); // current element

      if (intersectionStatus(ele)!=XFEM::TIMEINT::uncut_) // element is intersected
      {
        switch (ele->Shape())
        {
        case DRT::Element::hex8:
          oldValues<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement>(ele);
          break;
        case DRT::Element::hex20:
          oldValues<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement>(ele);
          break;
        default:
          dserror("xfem assembly type not yet implemented in time integration");
        };
      } // end if element intersected
  } // end loop over col ele nodes
}



/*------------------------------------------------------------------------------------------------*
 * compute values for different enrichments for the old interface position       winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
template <const int numnode>
void XFEM::EnrichmentProjection::oldValues(
    const DRT::Element*& ele
)
{
  const int nsd = 3;

  // velocity and pressure enrichment values of the elements nodes
  std::vector<LINALG::Matrix<nsd+1,numnode> > jumpEnrValues(oldVectors_.size(),LINALG::Matrix<nsd+1,numnode>(true));
  std::vector<LINALG::Matrix<nsd+1,numnode> > kinkEnrValues(oldVectors_.size(),LINALG::Matrix<nsd+1,numnode>(true));

  bool allEleNodesOnProc = true; // true if all nodes of an element are on current proc
  LINALG::Matrix<1,numnode> dmin(true); // distances from node to interface
  LINALG::Matrix<nsd,1> normal(true); // normal vector into interface normal direction
  int numnodeused = numnode; // number of nodes used for computation (=uncritical nodes)

  // compute required data
  getDataOfIntersectedEle<numnode>(ele,dmin,allEleNodesOnProc,jumpEnrValues,kinkEnrValues,normal);

  // sort out critical nodes
  analyseEnrichments<numnode>(ele,dmin,numnodeused);

  if(allEleNodesOnProc) // if this is false element is not needed
  {
    oldJumpAndKinkValues<numnode>(ele,dmin,normal,jumpEnrValues,numnodeused);
    oldKinkValues<numnode>(ele,dmin,kinkEnrValues,numnodeused);
  } // end if all nodes of element on proc
} // end function oldJumpAndKinkValues



/*------------------------------------------------------------------------------------------------*
 * compute values for jump enriched dofs for the old interface position          winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
template <const int numnode>
void XFEM::EnrichmentProjection::oldJumpAndKinkValues(
    const DRT::Element* ele,
    LINALG::Matrix<1,numnode>& dmin,
    LINALG::Matrix<3,1> normal,
    std::vector<LINALG::Matrix<4,numnode> >& jumpEnrValues,
    int& numnodeused
)
{
  const int nsd = 3;

  if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project)
  {
    if (intersectionStatus(ele)==XFEM::TIMEINT::cut_) // compute jump and kink in bisected ele
    {
      // compute the kink and jump height at the interface in this element
      // a least squares approach is used for the computation
      LINALG::Matrix<2,2> sysmat(true); // system matrix
      LINALG::Matrix<2,1> rhs(true); // right hand side

      LINALG::Matrix<2,1> currJumpAndKink(true); // jump and kink value
      LINALG::Matrix<2,nsd+1> currentJumpsAndKinks(true); // jump and kink values for all enrichments of current node
      std::vector<LINALG::Matrix<2,nsd+1> > eleJumpAndKinks; // all jump and kink values for a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over fields that shall be set
      {
        currentJumpsAndKinks.Clear();
        for (int ientry=0;ientry<nsd+1;ientry++) // number of enrichments for every node and field
        {
          // set up system matrix and right hand side
          sysmat.Clear();
          rhs.Clear();
          currJumpAndKink.Clear();

          sysmat(0,0) = numnodeused;
          for (int index=0;index<numnode;index++) // loop over element nodes
          {
            if (dmin(index)!=INFINITY) // if infinity then dont used node
            {
              rhs(0) -= jumpEnrValues[field](ientry,index);
              rhs(1) -= jumpEnrValues[field](ientry,index)*dmin(index);
              sysmat(1,0) += dmin(index);
              sysmat(1,1) -= dmin(index)*dmin(index);
            }
          }
          sysmat(0,1) = -sysmat(1,0);
          sysmat.Scale(0.5);

          // compute and set jump and kink value
          sysmat.Invert();
          currJumpAndKink.Multiply(sysmat,rhs);
          currentJumpsAndKinks(0,ientry) = currJumpAndKink(0);
          currentJumpsAndKinks(1,ientry) = currJumpAndKink(1);
        } // end loop over enrichment number

        eleJumpAndKinks.push_back(currentJumpsAndKinks);
      } // end loop over vectorfield size

      //      std::cout << "in bisected " << *ele << " final old jumps and kinks are " << eleJumpAndKinks[0];
      eleJumpsAndKinks_.insert(std::make_pair(ele->Id(),eleJumpAndKinks));
    } // end if element bisected
    else // in touchedminus ele just compute jump since kink value = 0 (singulary sysmat)
    {
      // compute the jump height at the interface in this element as a mean value
      double sysmat = 0; // system matrix
      double rhs = 0; // right hand side
      double currJump = 0; // jump value
      LINALG::Matrix<2,nsd+1> currentJumpsAndKinks(true); // jump and kink values for all enrichments of current node
      std::vector<LINALG::Matrix<2,nsd+1> > eleJumpAndKinks; // all jump and kink values for a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over fields that shall be set
      {
        currentJumpsAndKinks.Clear();
        for (int ientry=0;ientry<nsd+1;ientry++) // number of enrichments for every node and field
        {
          // set up system matrix and right hand side
          sysmat = 0.5*numnodeused;
          rhs = 0;
          for (int index=0;index<numnode;index++) // loop over element nodes
          {
            if (dmin(index)!=INFINITY) // if infinity then dont used node
              rhs -= jumpEnrValues[field](ientry,index);
          }

          // compute and set jump value
          currJump = rhs/sysmat;
          currentJumpsAndKinks(0,ientry) = currJump;

        } // end loop over enrichment number

        eleJumpAndKinks.push_back(currentJumpsAndKinks);
      } // end loop over vectorfield size

      //      std::cout << "final old jumps and kinks at touched " << *ele << " are " << eleJumpAndKinks[0];
      eleJumpsAndKinks_.insert(std::make_pair(ele->Id(),eleJumpAndKinks));
    } // end if element touchedminus
  }
  else if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project_scalar)
  {
    dserror("normal vector erroneous");
    if (intersectionStatus(ele)==XFEM::TIMEINT::cut_) // compute jump and kink in bisected ele
    {
      // compute the kink and jump height at the interface in this element
      // a least squares approach is used for the computation
      LINALG::Matrix<2,2> sysmat(true); // system matrix
      LINALG::Matrix<2,1> rhsvel(true); // right hand side for the velocity
      LINALG::Matrix<2,1> rhspres(true); // right hand side for the pressure
      LINALG::Matrix<2,1> velJumpAndKink(true); // jump and kink value for the velocity
      LINALG::Matrix<2,1> presJumpAndKink(true); // jump and kink value for the pressure
      LINALG::Matrix<2,nsd+1> currentJumpsAndKinks(true); // jump and kink values for all physical fields of a node
      std::vector<LINALG::Matrix<2,nsd+1> > eleJumpAndKinks; // all jump and kink values of a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over vectors that shall be set
      {
        // clear vectors and matrices
        currentJumpsAndKinks.Clear();
        sysmat.Clear();
        rhsvel.Clear();
        rhspres.Clear();
        velJumpAndKink.Clear();
        presJumpAndKink.Clear();

        // compute system matrix and right hand side
        sysmat(0,0) = numnodeused;
        for (int index=0;index<numnode;index++) // loop over element nodes
        {
          if (dmin(index)!=INFINITY) // if infinity then dont used node
          {
            sysmat(1,0) += dmin(index);
            sysmat(1,1) -= dmin(index)*dmin(index);

            LINALG::Matrix<nsd,1> currVelEnrValues(true);
            for (int i=0;i<nsd;i++)
              currVelEnrValues(i) = jumpEnrValues[field](i,index);

            double currVelEnrScalar = currVelEnrValues.Norm2();
            rhsvel(0) -= currVelEnrScalar;
            rhsvel(1) -= currVelEnrScalar*dmin(index);

            rhspres(0) -= jumpEnrValues[field](nsd,index);
            rhspres(1) -= jumpEnrValues[field](nsd,index)*dmin(index);
          }
        }
        sysmat(0,1) = -sysmat(1,0);
        sysmat.Scale(0.5);

        // compute and set jump and kink values
        sysmat.Invert();
        velJumpAndKink.Multiply(sysmat,rhsvel);
        presJumpAndKink.Multiply(sysmat,rhspres);
        currentJumpsAndKinks(0,0) = velJumpAndKink(0);
        currentJumpsAndKinks(1,0) = velJumpAndKink(1);
        currentJumpsAndKinks(0,1) = presJumpAndKink(0);
        currentJumpsAndKinks(1,1) = presJumpAndKink(1);
        eleJumpAndKinks.push_back(currentJumpsAndKinks);
      } // end loop over vectorfield size

//      std::cout << "final old jumps and kinks at bisected element are " << eleJumpAndKinks[0] << std::endl;
      eleJumpsAndKinks_.insert(std::make_pair(ele->Id(),eleJumpAndKinks));
    } // end if element bisected
    else // in touchedminus ele just compute jump since kink value = 0 (singulary sysmat)
    {
      // compute the jump height at the interface in this element as a mean value
      double sysmat = 0; // system matrix
      double rhsvel = 0; // right hand side for the velocity
      double rhspres = 0; // right hand side for the pressure
      double velJump = 0; // velocity jump
      double presJump = 0; // pressure jump
      LINALG::Matrix<2,nsd+1> currentJumpsAndKinks(true); // jump and kink values of the current node
      std::vector<LINALG::Matrix<2,nsd+1> > eleJumpAndKinks; // all jump and kink values of the current node

      for (size_t field=0;field<oldVectors_.size();field++) // lop over vectors that shall be set
      {
        currentJumpsAndKinks.Clear();

        // compute system matrix and right hand side
        sysmat = 0.5*numnodeused;
        rhsvel = 0;
        rhspres = 0;

        for (int index=0;index<numnode;index++) // loop over element nodes
        {
          if (dmin(index)!=INFINITY) // if infinity then dont used node
          {
            LINALG::Matrix<nsd,1> currVelEnrValues;
            for (int i=0;i<nsd;i++)
              currVelEnrValues(i) = jumpEnrValues[field](i,index);

            double currVelEnrScalar = currVelEnrValues.Dot(normal);
            rhsvel -= currVelEnrScalar;

            rhspres -= jumpEnrValues[field](3,index);
          }
        }

        velJump = rhsvel/sysmat;
        presJump = rhspres/sysmat;

        currentJumpsAndKinks(0,0) = velJump;
        currentJumpsAndKinks(0,1) = presJump;

        eleJumpAndKinks.push_back(currentJumpsAndKinks);
      } // end loop over vectorfield size
      //      std::cout << "final old jumps and kinks at touched " << *ele << " are " << eleJumpAndKinks[0] << std::endl;
      eleJumpsAndKinks_.insert(std::make_pair(ele->Id(),eleJumpAndKinks));
    } // end if element touchedminus
  }
}



/*------------------------------------------------------------------------------------------------*
 * compute values for kink enriched dofs for the old interface position          winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
template <const int numnode>
void XFEM::EnrichmentProjection::oldKinkValues(
    const DRT::Element* ele,
    LINALG::Matrix<1,numnode>& dmin,
    std::vector<LINALG::Matrix<4,numnode> >& kinkEnrValues,
    int& numnodeused
)
{
  const int nsd = 3;

  if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project)
  {
    if (intersectionStatus(ele)==XFEM::TIMEINT::cut_) // compute kink in bisected ele
    {
      double enrSum; // sum of enrichment values of element nodes
      LINALG::Matrix<1,nsd+1> currKinks(true); // kink values for all enrichments of current node
      std::vector<LINALG::Matrix<1,nsd+1> > eleKinks; // all kink values for a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over fields that shall be set
      {
        currKinks.Clear();
        for (int ientry=0;ientry<nsd+1;ientry++) // number of enrichments for every node and field
        {
          enrSum = 0.0;
          for (int index=0;index<numnode;index++) // loop over element nodes
          {
            if (dmin(index)!=INFINITY) // if infinity then don't use node
              enrSum += kinkEnrValues[field](ientry,index);
          }
          currKinks(ientry) = enrSum/numnodeused;
        } // end loop over enrichment number
        eleKinks.push_back(currKinks);
      } // end loop over vectorfield size
      //      std::cout << "final old kink values at bisected " << *ele << " are " << eleKinks[0];
      eleKinks_.insert(std::make_pair(ele->Id(),eleKinks));
    } // end if element bisected
    else
    {
      LINALG::Matrix<1,nsd+1> currKinks(true); // kink values for all enrichments of current node
      std::vector<LINALG::Matrix<1,nsd+1> > eleKinks; // all kink values for a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over fields that shall be set
        eleKinks.push_back(currKinks);

      eleKinks_.insert(std::make_pair(ele->Id(),eleKinks));
    }
  }
  else if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project_scalar)
  {
    if (intersectionStatus(ele)==XFEM::TIMEINT::cut_) // compute kink in bisected ele
    {
      double enrSumVel; // sum of scalar velocity enrichment values of element nodes
      double enrSumPres; // sum of pressure enrichment values

      LINALG::Matrix<nsd,1> enrVel;

      LINALG::Matrix<1,nsd+1> currKinks(true); // kink values for all enrichments of current node
      std::vector<LINALG::Matrix<1,nsd+1> > eleKinks; // all kink values for a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over fields that shall be set
      {
        currKinks.Clear();
        enrSumPres = 0;

        for (int index=0;index<numnode;index++) // loop over element nodes
        {
          if (dmin(index)!=INFINITY) // if infinity then don't use node
          {
            enrSumPres += kinkEnrValues[field](3,index);
            for (int ientry=0;ientry<nsd;ientry++) // number of enrichments for every node and field
              enrVel(0,ientry) += kinkEnrValues[field](ientry,index);
            enrSumVel = enrVel.Norm2();
          } // end if distance infinity
        } // end loop over element nodes
        enrSumVel = enrVel.Norm2();

        currKinks(0) = enrSumVel/numnodeused;
        currKinks(1) = enrSumPres/numnodeused;

        eleKinks.push_back(currKinks);
      } // end loop over vectorfield size
      //      std::cout << "final old kink values at bisected element are " << eleKinks[0];
      eleKinks_.insert(std::make_pair(ele->Id(),eleKinks));
    } // end if element bisected
    else
    {
      LINALG::Matrix<1,nsd+1> currKinks(true); // kink values for all enrichments of current node
      std::vector<LINALG::Matrix<1,nsd+1> > eleKinks; // all kink values for a node

      for (size_t field=0;field<oldVectors_.size();field++) // loop over fields that shall be set
        eleKinks.push_back(currKinks);

      //      std::cout << "final old kink values at touched " << *ele << " are " << eleKinks[0];
      eleKinks_.insert(std::make_pair(ele->Id(),eleKinks));
    } // end if element touched
  }
}



/*------------------------------------------------------------------------------------------------*
 * standard computation of the new enrichment values                             winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::computeNewEnrichments(
)
{
  const int nsd = 3; // dimension
  for (int inode=0;inode<discret_->NumMyRowNodes();inode++) // loop over nodes
  {
    const DRT::Node* currnode = discret_->lRowNode(inode); // current node

    if (newEnrValueNeeded(currnode)) // current node needs new enrichment value
    {
      // vector of elements located around this node
      std::vector<const DRT::Element*> eles;

      addPBCelements(currnode, eles);
      const int numeles = eles.size(); // number of elements located around this node

      std::vector<LINALG::Matrix<2,nsd+1> > currJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true)); // currently computed jump and kink value
      std::vector<LINALG::Matrix<2,nsd+1> > averageJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true)); // jump and kink values for current node

      std::vector<LINALG::Matrix<1,nsd+1> > currKinks(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // currently computed jump and kink value
      std::vector<LINALG::Matrix<1,nsd+1> > averageKinks(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // jump and kink values for current node

#ifndef COMBUST_SETJUMP
      // compute average jump and kink value around node if elements are enriched
      int numOldIntersectedEle=0; // index for number of bisected elements at old interface position
      for (int iele=0;iele<numeles;iele++)
      {
        const int elegid = eles[iele]->Id(); // global id of current element

        if (intersectionStatus(eles[iele])!=XFEM::TIMEINT::uncut_) // element intersected
        {
          numOldIntersectedEle++;
          currJumpsAndKinks = eleJumpsAndKinks_.find(elegid)->second;
          currKinks = eleKinks_.find(elegid)->second;
          for (size_t ivector=0;ivector<newVectors_.size();ivector++)
          {
            averageJumpsAndKinks[ivector] += currJumpsAndKinks[ivector];
            averageKinks[ivector] += currKinks[ivector];
          }
        } // end if element bisected
      } // end loop over elements containing the node
      //      std::cout << "summed jumps and kinks are " << averageJumpsAndKinks[0] << std::endl;

      std::vector<LINALG::Matrix<1,nsd+1> > finalEnrichmentValues(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true));
      if (numOldIntersectedEle > 0) // >=1 elements around the node where enriched -> computed value can be used
      {
        for (size_t ivector=0;ivector<newVectors_.size();ivector++) // scaling
        {
          averageJumpsAndKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
          averageKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
        }
        computeJumpEnrichmentValues(currnode,averageJumpsAndKinks);
        computeKinkEnrichmentValues(currnode,averageKinks);
      } // end if old bisected element around new enriched node
      else // no old enriched element around new enriched node -> other way of computation required
      {
        TimeIntData data(
            *currnode,
            INFINITY,
            1,
            std::vector<LINALG::Matrix<2,4> >(newVectors_.size(),LINALG::Matrix<2,4>(true)),
            std::vector<LINALG::Matrix<1,4> >(newVectors_.size(),LINALG::Matrix<1,4>(true)));
        timeIntData_->push_back(data);
      }
#else
      setJumpEnrichmentValues(currnode);
      setKinkEnrichmentValues(currnode);
#endif
    } // end if node is enriched
  } // end loop over row nodes
}// end function computeNewEnrichments



/*------------------------------------------------------------------------------------------------*
 * alternative computation of the new enrichment values                          winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::handleFailedNodes(
)
{
  const int nsd = 3;
  const double TOL = 1.0e-3;

  // evaluate all enriched nodes at old interface position
  std::set<int> oldEnrNodes; // set of all enriched nodes at old interface position
  for (int inode=0;inode<discret_->NumMyRowNodes();inode++) // loop over processor nodes
  {
    const int nodegid = discret_->lRowNode(inode)->Id(); // current node

    // set nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldenrset(olddofman_->getNodeDofSet(nodegid));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump
          || fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeKink)
      {
        oldEnrNodes.insert(nodegid);
        break;
      }
    }
  } // end loop over processor nodes

#ifdef PARALLEL
  // get the nearest enriched node at old timestep
  for (int iproc=0;iproc<numproc_;iproc++) // loop over processors
  {
#endif // PARALLEL
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      for (std::set<int>::const_iterator enrnode=oldEnrNodes.begin();
          enrnode!=oldEnrNodes.end(); enrnode++)
      {
        DRT::Node* oldnode = discret_->gNode(*enrnode);
        double oldphivalue = (*phin_)[oldnode->LID()]; // distance from old node to interface as second indicator
        LINALG::Matrix<nsd,1> oldCoords(oldnode->X()); // coordinates of old enriched node
        LINALG::Matrix<nsd,1> diff(true); // difference between new and old enriched node
        diff.Update(1.0,LINALG::Matrix<nsd,1>(data->node_.X()),-1.0,oldCoords);

        double finalDiff = diff.Norm2()+fabs(oldphivalue); // diff is the main indicator, the phi-value the second if two nodes have (nearly) equal distance to the new node

        if (finalDiff-data->phiValue_+TOL*finalDiff < 0) // new nearest node and more nearer than the old + some tolerance (both nodes not equal good -> prevents special case)
        {
          data->phiValue_ = finalDiff; // new smallest difference
          data->counter_ = 1; // one single nearest node

          const DRT::Node* currOldNode = discret_->gNode(*enrnode);

          // vector of elements located around this node
          std::vector<const DRT::Element*> eles;
          addPBCelements(currOldNode,eles);

          const int numeles = eles.size();

          std::vector<LINALG::Matrix<2,nsd+1> > finalJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true)); // final jump and kink values for node
          std::vector<LINALG::Matrix<2,nsd+1> > currJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true)); // jump and kink values of one element around the node

          std::vector<LINALG::Matrix<1,nsd+1> > finalKinks(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // final jump and kink values for node
          std::vector<LINALG::Matrix<1,nsd+1> > currKinks(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // jump and kink values of one element around the node

          int numOldIntersectedEle=0; // index how much elements are intersected
          for (int iele=0;iele<numeles;iele++)
          {
            const int elegid = eles[iele]->Id(); // global id of current element
            if (intersectionStatus(eles[iele])!=XFEM::TIMEINT::uncut_) // element intersected
            {
              numOldIntersectedEle++;

              currJumpsAndKinks=eleJumpsAndKinks_.find(elegid)->second;
              currKinks=eleKinks_.find(elegid)->second;

              for(size_t ivector=0;ivector<oldVectors_.size();ivector++)
              {
                finalJumpsAndKinks[ivector] += currJumpsAndKinks[ivector];
                finalKinks[ivector]+=currKinks[ivector];
              }
            } // end if element bisected
          } // end loop over elements around the node

          if (numOldIntersectedEle==0)
          {
            std::cout << *currOldNode << " is enriched at old timestep" << std::endl;
            dserror("Enriched node shall have an intersected element around");
          }

          for(size_t ivector=0;ivector<oldVectors_.size();ivector++) // scaling loop
          {
            finalJumpsAndKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
            finalKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
          }

          data->jumpAndKinkValues_ = finalJumpsAndKinks;
          data->kinkValues_ = finalKinks;
        } // end if new nearest node
        else if ((finalDiff-data->phiValue_ > -TOL*finalDiff) and
            (finalDiff-data->phiValue_ < TOL*finalDiff)) // handles special case that two nodes are very similar near
        {
          data->phiValue_ = (finalDiff+data->phiValue_)/2.0; // average between the two equal near nodes
          data->counter_ += 1; // increase number of nearest nodes by 1

          DRT::Node* currOldNode = discret_->gNode(*enrnode);

          // vector of elements located around this node
          std::vector<const DRT::Element*> eles;
          addPBCelements(currOldNode,eles);

          const int numeles = eles.size();

          std::vector<LINALG::Matrix<2,nsd+1> > finalJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true)); // final jump and kink values for node
          std::vector<LINALG::Matrix<2,nsd+1> > currJumpsAndKinks(newVectors_.size(),LINALG::Matrix<2,nsd+1>(true)); // jump and kink values of one element around the node

          std::vector<LINALG::Matrix<1,nsd+1> > finalKinks(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // final jump and kink values for node
          std::vector<LINALG::Matrix<1,nsd+1> > currKinks(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // jump and kink values of one element around the node


          int numOldIntersectedEle=0; // index how much elements are intersected
          for (int iele=0;iele<numeles;iele++)
          {
            const int elegid = eles[iele]->Id(); // global id of current element

            if (intersectionStatus(eles[iele])!=XFEM::TIMEINT::uncut_) // element intersected
            {
              numOldIntersectedEle++;

              currJumpsAndKinks=eleJumpsAndKinks_.find(elegid)->second;
              currKinks=eleKinks_.find(elegid)->second;

              for(size_t ivector=0;ivector<oldVectors_.size();ivector++)
              {
                finalJumpsAndKinks[ivector] += currJumpsAndKinks[ivector];
                finalKinks[ivector]+=currKinks[ivector];
              }
            } // end if element bisected
          } // end loop over elements around the node

          if (numOldIntersectedEle==0)
          {
            std::cout << *currOldNode << " is enriched at old timestep" << std::endl;
            dserror("Enriched node shall have an intersected element around");
          }

          for(size_t ivector=0;ivector<oldVectors_.size();ivector++) // scaling loop
          {
            finalJumpsAndKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
            data->jumpAndKinkValues_[ivector] += finalJumpsAndKinks[ivector];

            finalKinks[ivector].Scale(1.0/static_cast<double>(numOldIntersectedEle));
            data->kinkValues_[ivector] += finalKinks[ivector];
          }
        } // end if two equal near nodes
      } // loop over enriched nodes on current proc at old interface
    } // loop over failed nodes

#ifdef PARALLEL
    exportEnrichmentData();
  } // loop over processors
#endif // PARALLEL

  // remark: after n loops over the n processors the failed data is again at the
  // processor where it should be. It contains the nearest node that was enriched
  // at the old interface position and all required according data at this node

  // compute new enrichment values
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    for(size_t ivector=0;ivector<oldVectors_.size();ivector++) // scaling loop
    {
      data->jumpAndKinkValues_[ivector].Scale(1.0/static_cast<double>(data->counter_));
      data->kinkValues_[ivector].Scale(1.0/static_cast<double>(data->counter_));
    }
    computeJumpEnrichmentValues(&data->node_,data->jumpAndKinkValues_);
    computeKinkEnrichmentValues(&data->node_,data->kinkValues_);
  } // end loop over failed nodes
} // end function handleFailedNodes



/*------------------------------------------------------------------------------------------------*
 * compute and set values for jump enriched dofs for the new interface position  winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::computeJumpEnrichmentValues(
    const DRT::Node* node,
    std::vector<LINALG::Matrix<2,4> > jumpsAndKinks
)
{
  const int nsd = 3; // dimension
  int numNewIntersectedEle = 0; // number of intersected elements at new interface position

  // vector of elements located around this node
  std::vector<const DRT::Element*> eles;
  addPBCelements(node,eles);
  const int numeles = eles.size();

  std::vector<LINALG::Matrix<1,nsd+1> > finalEnrichmentValues(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // enrichment values for the node

  for (int iele=0;iele<numeles;iele++) // loop over elements around enriched node
  {
    if (intersectionStatus(eles[iele],false)!=XFEM::TIMEINT::uncut_) // element intersected
    {
      numNewIntersectedEle++;
      double dist = 0.0; // minimal distance from node to interface segment of current element
      LINALG::Matrix<nsd,1> normal(true); // normal vector from node to interface segment of current element
      LINALG::Matrix<nsd,1> dummy(true);
      SignedDistance(node,eles[iele]->Id(),dist,normal,dummy,false);
      double phiDiff = dist - (*phin_)[node->LID()]; // difference of the signed distances of the node to interface between this and last time step

      for (size_t ivector=0;ivector<newVectors_.size();ivector++) // loop over global vectors
      {
        if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project)
        {
          for (int entry=0;entry<nsd+1;entry++)
            finalEnrichmentValues[ivector](entry) +=
                -0.5*(jumpsAndKinks[ivector](0,entry)+phiDiff*jumpsAndKinks[ivector](1,entry)) // -0.5*new_jump
                +0.5*dist*jumpsAndKinks[ivector](1,entry); // +0.5*kink*signeddistance
        }
        else if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project_scalar)
        {
          dserror("normal vector erroneous");
          // velocity entry
          for (int entry=0;entry<nsd;entry++)
            finalEnrichmentValues[ivector](entry) +=
                -normal(entry)*((-0.5)*(jumpsAndKinks[ivector](0,0)+phiDiff*jumpsAndKinks[ivector](1,0)) // -0.5*new_jump
                    +0.5*dist*jumpsAndKinks[ivector](1,0)); // +0.5*kink*signeddistance

          // pressure entry
          finalEnrichmentValues[ivector](nsd) +=
              -0.5*(jumpsAndKinks[ivector](0,1)+phiDiff*jumpsAndKinks[ivector](1,1)) // -0.5*jump
              +0.5*dist*jumpsAndKinks[ivector](1,1); // +0.5*kink*signeddistance
        }
      } // end loop over vector size
    } // end if element bisected
  } // end loop over elements containing the node
  if (numNewIntersectedEle==0)
  {
    std::cout << std::endl << std::endl << *node << std::endl;
    dserror("node is enriched at new time step but has no intersected or touched element around");
  }
  for (size_t ivector=0;ivector<newVectors_.size();ivector++)
    finalEnrichmentValues[ivector].Scale(1.0/static_cast<double>(numNewIntersectedEle));


  int indexJumpEnr=0; // index which entry has to be used
  // set nodal velocities and pressures with help of the field set of node
  const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(node->Id()));
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
      fieldenr != fieldenrset.end();++fieldenr)
  {
    if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
    {
      const DofKey newdofkey(node->Id(), *fieldenr);
      const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
      const int lid = newdofrowmap_.LID(newdofpos);

      //      std::cout << *node << ": enrichment value "<< (*newVectors_[0])[lid] << " becomes " << finalEnrichmentValues[0](0,indexJumpEnr) << std::endl;
      for (size_t index=0;index<vectorSize(TimeIntData::predictor_);index++)
        (*newVectors_[index])[lid] = finalEnrichmentValues[index](0,indexJumpEnr);
      indexJumpEnr++;
    } // end if enrichment type jump enrichment
  } // end loop over fieldenr
}



/*------------------------------------------------------------------------------------------------*
 * compute and set values for kink enriched dofs for the new interface position  winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::computeKinkEnrichmentValues(
    const DRT::Node* node,
    std::vector<LINALG::Matrix<1,4> > kinks
)
{
  int numNewIntersectedEle = 0; // number of intersected elements at new interface position
  const int nsd = 3; // dimension

  // vector of elements located around this node
  std::vector<const DRT::Element*> eles;
  addPBCelements(node,eles);

  const int numeles=eles.size();
  std::vector<LINALG::Matrix<1,nsd+1> > finalEnrichmentValues(newVectors_.size(),LINALG::Matrix<1,nsd+1>(true)); // enrichment values for the node

  for (int iele=0;iele<numeles;iele++) // loop over elements around enriched node
  {
    if (intersectionStatus(eles[iele],false)==XFEM::TIMEINT::cut_) // element intersected
    {
      numNewIntersectedEle++;
      double dist = 0.0; // currently dummy
      LINALG::Matrix<nsd,1> normal(true); // normal vector from node to interface segment of current element
      LINALG::Matrix<nsd,1> dummy(true);
      SignedDistance(node,eles[iele]->Id(),dist,normal,dummy,false);

      for (size_t ivector=0;ivector<newVectors_.size();ivector++) // loop over global vectors
      {
        for (int entry=0;entry<nsd+1;entry++)
          finalEnrichmentValues[ivector](entry) += kinks[ivector](entry);

        if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project)
        {
          for (int entry=0;entry<nsd+1;entry++)
            finalEnrichmentValues[ivector](entry) += 0.5*kinks[ivector](entry);
        }
        else if (timeIntEnr_==INPAR::COMBUST::xfemtimeintenr_project_scalar)
        {
          dserror("normal vector erroneous");
          // velocity entry
          for (int entry=0;entry<nsd;entry++)
            finalEnrichmentValues[ivector](0) +=
                -normal(entry)*(0.5*kinks[ivector](0)); // +0.5*kink*signeddistance

          // pressure entry
          finalEnrichmentValues[ivector](nsd) += 0.5*kinks[ivector](1); // +0.5*kink*signeddistance
        }
      } // end loop over vector size
    } // end if element intersected
  } // end loop over elements containing the node
  for (size_t ivector=0;ivector<newVectors_.size();ivector++)
    finalEnrichmentValues[ivector].Scale(1.0/static_cast<double>(numNewIntersectedEle));
  //  std::cout << *node << " has final enr values " << finalEnrichmentValues[0] << std::endl;

  int indexKinkEnr=0; // index which entry has to be used
  // set nodal velocities and pressures with help of the field set of node
  const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(node->Id()));
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
      fieldenr != fieldenrset.end();++fieldenr)
  {
    if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeKink)
    {
      const DofKey newdofkey(node->Id(), *fieldenr);
      const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
      const int lid = newdofrowmap_.LID(newdofpos);

      for (size_t index=0;index<vectorSize(TimeIntData::predictor_);index++)
        (*newVectors_[index])[lid] = finalEnrichmentValues[index](0,indexKinkEnr);
      indexKinkEnr++;
    } // end if enrichment type jump enrichment
  } // end loop over fieldenr
}



/*------------------------------------------------------------------------------------------------*
 * compute and set values for jump enriched dofs for the new interface position  winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::setJumpEnrichmentValues(
    const DRT::Node* node
)
{
  // get the global id for current node
  const int gid = node->Id();

  // get local processor id according to global node id
  const int nodelid = (*gradphi_).Map().LID(gid);
  if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphi_).Comm().MyPID(),gid);

  const int numcol = (*gradphi_).NumVectors();
  if( numcol != 3) dserror("number of columns in Epetra_MultiVector is not 3");

  //-------------------------------------------------------------
  // get (smoothed) gradient of the G-function field at this node
  //-------------------------------------------------------------
  LINALG::Matrix<3,1> mygradphi(true);

  // loop over dimensions (= number of columns in multivector)
  for(int col=0; col< numcol; col++)
  {
    // get columns vector of multivector
    double* globalcolumn = (*gradphi_)[col];
    // set smoothed gradient entry of phi into column of global multivector
    mygradphi(col) = globalcolumn[nodelid];
  }

  //------------------------------------
  // compute smoothed unit normal vector
  // n = - grad phi / |grad phi|
  //------------------------------------
  // smoothed normal vector at this node
  LINALG::Matrix<3,1> nvec = mygradphi;
  // compute norm of smoothed normal vector
  const double ngradphi = mygradphi.Norm2();

  // divide vector by its norm to get unit normal vector
  if (fabs(ngradphi < 1.0E-12))// 'ngradnorm' == 0.0
  {
    // length of smoothed normal is zero at this node -> node must be on a singularity of the
    // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
    // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
    //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
    std::cout << "\n/!\\ phi gradient too small at node " << gid << " -> flame velocity is only the convective velocity" << std::endl;
    nvec.PutScalar(0.0);
  }
  else
  {
    nvec.Scale(-1.0/ngradphi);
  }

  LINALG::Matrix<3,1> veljump(true);
  veljump.Update(veljump_,nvec);

  double wallfac = 1.0;
#ifdef ORACLES
  //--------------------------------------------------------
  // check physical coordinates of node
  // remark: we want to blend the flame speed close to walls
  //--------------------------------------------------------

  //    wall
  // 1.0 |     ______
  //     |    /
  // 0.0 |___/ ,
  //     |     H/6
  const double wallzone = 0.0299/6.0;
  if (node->X()[0] > 0.0) // inside combustion chamber
  {
    if ( (0.0653-abs(node->X()[1])) < wallzone or // close to top or bottom wall
                       node->X()[0] < wallzone )  // close to step
    {
      // wall factor is 0 at the wall and 1 at H/6 or further away from the wall
      wallfac = 6.0/0.0299 * std::min(0.0653-abs(node->X()[1]),node->X()[0]);
      if (wallfac < 0.1) // cut off the last 10% to guarantee a zero jump at the wall
        wallfac = 0.0;
    }
  }
#endif
  veljump.Scale(wallfac);

  // set nodal velocities and pressures with help of the field set of node
  const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(node->Id()));
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
      fieldenr != fieldenrset.end();++fieldenr)
  {
    if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeJump)
    {
      const DofKey newdofkey(node->Id(), *fieldenr);
      const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
      const int lid = newdofrowmap_.LID(newdofpos);

      // velocity set for old and new time step
      for (size_t i=0;i<vectorSize(TimeIntData::predictor_);i=i+2)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          (*newVectors_[i])[lid] = -0.5*veljump(0);
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          (*newVectors_[i])[lid] = -0.5*veljump(1);
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          (*newVectors_[i])[lid] = -0.5*veljump(2);
      }

      // acceleration set to zero for old and new time step
      for (size_t i=1;i<vectorSize(TimeIntData::predictor_);i=i+2)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          (*newVectors_[i])[lid] = 0.0;
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          (*newVectors_[i])[lid] = 0.0;
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          (*newVectors_[i])[lid] = 0.0;
      }
    } // end if enrichment type jump enrichment
  } // end loop over fieldenr
}



/*------------------------------------------------------------------------------------------------*
 * compute and set values for kink enriched dofs for the new interface position  winklmaier 08/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::setKinkEnrichmentValues(
    const DRT::Node* node
)
{
  // set nodal velocities and pressures with help of the field set of node
  const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(node->Id()));
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
      fieldenr != fieldenrset.end();++fieldenr)
  {
    if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeKink)
    {
      const DofKey newdofkey(node->Id(), *fieldenr);
      const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;
      const int lid = newdofrowmap_.LID(newdofpos);

      for (size_t index=0;index<vectorSize(TimeIntData::predictor_);index++)
        (*newVectors_[index])[lid] = 0.0;
    } // end if enrichment type jump enrichment
  } // end loop over fieldenr
}



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export enrichment data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::EnrichmentProjection::exportEnrichmentData()
{
  const int nsd = 3;
  // destination proc (the "next" one)
  int dest = myrank_+1;
  if(myrank_ == (numproc_-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank_-1;
  if(myrank_ == 0)
    source = numproc_-1;

  // vector including all data that has to be send to the next proc
  DRT::PackBuffer dataSend;

  // packing the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::failedEnr_)
    {
      packNode(dataSend,data->node_);
      DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
      DRT::ParObject::AddtoPack(dataSend,data->counter_);
      DRT::ParObject::AddtoPack(dataSend,data->jumpAndKinkValues_);
      DRT::ParObject::AddtoPack(dataSend,data->kinkValues_);
    }
    else
      dserror("In enrichment computation all data should have enrichment status");
  }

  dataSend.StartPacking();

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::failedEnr_)
    {
      packNode(dataSend,data->node_);
      DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
      DRT::ParObject::AddtoPack(dataSend,data->counter_);
      DRT::ParObject::AddtoPack(dataSend,data->jumpAndKinkValues_);
      DRT::ParObject::AddtoPack(dataSend,data->kinkValues_);
    }
  }


  std::vector<char> dataRecv;
  sendData(dataSend,dest,source,dataRecv);

  // pointer to current position of group of cells in global std::string (counts bytes)
  std::vector<char>::size_type posinData = 0;

  // clear vector that should be filled
  clearState(TimeIntData::failedEnr_);

  // unpack received data
  while (posinData < dataRecv.size())
  {
    double coords[nsd] = {0.0};
    DRT::Node node(0,(double*)coords,0);
    double phiValue;
    int counter;
    std::vector<LINALG::Matrix<2,4> > jumpAndKinkValues;
    std::vector<LINALG::Matrix<1,4> > kinkValues;

    unpackNode(posinData,dataRecv,node);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,counter);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,jumpAndKinkValues);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,kinkValues);

    TimeIntData data(node,phiValue,counter,jumpAndKinkValues,kinkValues);
    timeIntData_->push_back(data);
  }

  discret_->Comm().Barrier();
} // end exportEnrichmentData

#endif // parallel

