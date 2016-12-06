/*!------------------------------------------------------------------------------------------------*
\file timeInt_std_SemiLagrange.cpp

\brief provides the SemiLagrange class

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#include "timeInt_std_extrapolation.H"
#include "timeInt_std_SemiLagrange.H"
#include "../linalg/linalg_utils.H"


/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrange Back-Tracking algorithm constructor                             winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
XFEM::SemiLagrange::SemiLagrange(
    XFEM::TIMEINT& timeInt,
    INPAR::COMBUST::XFEMTimeIntegration& timeIntType,
    const Teuchos::RCP<Epetra_Vector> veln,
    const double& dt,
    const double& theta,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const double& veljump,
    bool initialize
) :
STD(timeInt,timeIntType,veln,dt,flamefront,initialize),
theta_default_(theta),
veljump_(veljump)
{
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrangian Back-Tracking main algorithm                                  winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::compute(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  handleVectors(newRowVectorsn,newRowVectorsnp);
  newIteration_prepare(newVectors_);

  switch (FGIType_)
  {
  case FRS1FGI1_:
  {
    resetState(TimeIntData::basicStd_,TimeIntData::currSL_);
    break;
  }
  case FRSNot1_:
  {
    resetState(TimeIntData::doneStd_,TimeIntData::currSL_);
    break;
  }
  case FRS1FGINot1_:
  {
    reinitializeData();
    resetState(TimeIntData::basicStd_,TimeIntData::currSL_);
    resetState(TimeIntData::doneStd_,TimeIntData::currSL_);
    break;
  }
  default:
  {
    dserror("not implemented");
    break;
  }
  } // end switch

  /*----------------------------------------------------*
   * first part: get the correct origin for the node    *
   * in a lagrangian point of view using a newton loop  *
   *----------------------------------------------------*/
#ifdef PARALLEL
  int counter = 0; // loop counter to avoid infinite loops

  // loop over nodes which still don't have and may get a good startvalue
  while (true)
  {
    counter += 1;

    // counter limit because maximal max_iter newton iterations with maximal
    // numproc processor changes per iteration (avoids infinite loop)
    if (!globalNewtonFinished(counter))
    {
#endif
      // loop over all nodes which changed interface side
      // remark: negative loop so that deleted elements don't influence other elements position
      for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
          data!=timeIntData_->end(); data++)
      {
        if (data->state_ == TimeIntData::currSL_)
        {
          //cout << endl << "on proc " << myrank_ << " iter " << data->counter_ <<
          //    " starts for " << data->node_ << " and init startpoint " << data->startpoint_;
          // Initialization
          DRT::Element* ele = NULL; // pointer to the element where start point lies in
          LINALG::Matrix<nsd,1> xi(true); // local transformed coordinates of x
          LINALG::Matrix<nsd,1> vel(true); // velocity of the start point approximation
          double phin = 0.0; // phi-value of the start point approximation
          bool elefound = false;             // true if an element for a point was found on the processor
          bool stdBackTracking = true; // true if standard back-tracking shall be done

          // search for an element where the current startpoint lies in
          // if found, give out all data at the startpoint
          elementSearch(ele,data->startpoint_,xi,vel,phin,elefound);

          // if element is not found, look at another processor and so add all
          // according data to the vectors which will be sent to the next processor
          if (!elefound)
          {
            if (data->searchedProcs_ < numproc_+10) // safety + 10
            {
              data->state_ = TimeIntData::nextSL_;
              data->searchedProcs_ += 1;
            }
            else // all procs searched -> point not in domain
            {
              data->state_ = TimeIntData::failedSL_;
              std::cout << "WARNING! Lagrangian start point not in domain!" << std::endl;
            }
          } // end if elefound
          else  // if element is found, the newton iteration to find a better startpoint can start
          {
            if (interfaceSideCompare(ele,data->startpoint_,0,data->phiValue_) == false) // Lagrangian origin and original node on different interface sides
            {
              if (data->counter_==0)
                data->state_ = TimeIntData::initfailedSL_;
              else
                data->state_ = TimeIntData::failedSL_;
            }
            else  // Newton loop just for sensible points
              NewtonLoop(ele,&*data,xi,vel,phin,elefound,stdBackTracking);

            // if iteration can go on (that is when startpoint is on
            // correct interface side and iter < max_iter)
            if ((data->counter_<newton_max_iter_) and (data->state_==TimeIntData::currSL_))
            {
              // if element is not found in a newton step, look at another processor and so add
              // all according data to the vectors which will be sent to the next processor
              if (!elefound)
              {
                data->searchedProcs_ = 2;
                data->state_ = TimeIntData::nextSL_;
              }
              else // newton iteration converged to a good startpoint and so the data can be used to go on
                callBackTracking(ele,&*data,xi,"standard");
            } // end if

            if (data->counter_ == newton_max_iter_) // maximum number of iterations reached
            { // do not use the lagrangian origin since this case is strange and potential dangerous
              std::cout << "WARNING: newton iteration to find start value did not converge!" << std::endl;
              data->state_ = TimeIntData::failedSL_;
            }
          } // end if over nodes which changed interface

          findNearProcs(&*data,ele);
          //std::cout << "on proc " << myrank_ << " after " << data->counter_ << " iterations "
          //    << data->node_ << " has startpoint " << data->startpoint_ <<
          //    " and state " << data->stateToString() << std::endl;
        }
      } // end loop over all nodes with changed interface side
    } // end if movenodes is empty or max. iter reached
    else //
      resetState(TimeIntData::currSL_,TimeIntData::failedSL_);

#ifdef PARALLEL
    // export nodes and according data for which the startpoint isn't still found (next_ vector) to next proc
    bool procDone = globalNewtonFinished();
    exportIterDataNew(procDone);

    // convergencecheck: procfinished == 1 just if all procs have finished
    if (procDone)
      break;
  } // end while loop over searched nodes
#endif



  /*-----------------------------------------------------------------------------*
   * second part: get sensible startvalues for nodes where the algorithm failed, *
   * using another algorithm, and combine the "Done" and the "Failed" - vectors  *
   *-----------------------------------------------------------------------------*/
  if (FGIType_==FRSNot1_) // failed nodes stay equal after the first computation
    clearState(TimeIntData::failedSL_);
  else
  {
    getDataForNotConvergedNodes(); // compute final data for failed nodes
  }



  /*-----------------------------------------------------------*
   * third part: set the computed values into the state vector *
   *-----------------------------------------------------------*/
#ifdef PARALLEL
  // send the computed startvalues for every node which needs
  // new start data to the processor where the node is
  exportFinalData();
#endif

  // now every proc has the whole data for the nodes and so the data can be set to the right place now
  setFinalData();

#ifdef DEBUG
#ifdef PARALLEL
  if (counter > 8*numproc_) // too much loops shouldnt be if all this works
    std::cout << "WARNING: semiLagrangeExtrapolation seems to run an infinite loop!" << std::endl;
#endif
#endif

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
  }
} // end semiLagrangeExtrapolation



/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::NewtonLoop(
    DRT::Element*& ele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& elefound,
    bool& stdBackTracking
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  stdBackTracking = true; // standard Newton loop -> standard back tracking

  // Initialization
  LINALG::Matrix<nsd,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<nsd,1> incr(true);                 // increment of the newton system

  const double relTolIncr = 1.0e-10;   // tolerance for the increment
  const double relTolRes = 1.0e-10;    // tolerance for the residual

  LINALG::Matrix<nsd,1> origNodeCoords(true); // coordinates of endpoint of Lagrangian characteristics
  for (int i=0;i<nsd;i++)
    origNodeCoords(i) = data->node_.X()[i];

  // initialize residual (Theta = 0 at predictor step)
  residuum.Clear();
  residuum.Update((1.0-Theta(data)),vel,Theta(data),data->vel_); // dt*v(data->startpoint_)
  residuum.Update(1.0,data->startpoint_,-1.0,origNodeCoords,dt_);  // R = data->startpoint_ - data->movNode_ + dt*v(data->startpoint_)

  while(data->counter_ < newton_max_iter_)  // newton loop
  {
    data->counter_ += 1;

    switch (ele->Shape())
    {
    case DRT::Element::hex8:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      NewtonIter<numnode,DRT::Element::hex8>(ele,data,xi,vel,residuum,incr,phi,elefound);
    }
    break;
    case DRT::Element::hex20:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
      NewtonIter<numnode,DRT::Element::hex20>(ele,data,xi,vel,residuum,incr,phi,elefound);
    }
    break;
    default:
    {
      dserror("xfem assembly type not yet implemented in time integration");
      break;
    }
    }; // end switch element type

    if (elefound) // element of data->startpoint_ at this processor
    {
      if (interfaceSideCompare(ele,data->startpoint_,0,data->phiValue_) == false)
      {
        data->state_ = TimeIntData::failedSL_;
        break; // leave newton loop if point is on wrong domain side
      }
      else
      {
        // reset residual
        residuum.Clear();
        residuum.Update((1.0-Theta(data)),vel,Theta(data),data->vel_); // dt*v(data->startpoint_)
        residuum.Update(1.0,data->startpoint_,-1.0,origNodeCoords,dt_);  // R = data->startpoint_ - data->movNode_ + dt*v(data->startpoint_)

        // convergence criterion
        if (data->startpoint_.Norm2()>1e-3)
        {
          if (incr.Norm2()/data->startpoint_.Norm2() < relTolIncr && residuum.Norm2()/data->startpoint_.Norm2() < relTolRes)
            break;
        }
        else
        {
          if (incr.Norm2() < relTolIncr && residuum.Norm2() < relTolRes)
            break;
        }
      } // end if interface side is the same
    } // end if elefound is true
    else // element of data->startpoint_ not at this processor
      break; // stop newton loop on this proc
  } // end while Newton loop
#ifdef DEBUG
  // did newton iteration converge?
  if(data->counter_ == newton_max_iter_){std::cout << "WARNING: newton iteration for finding start value not converged for point\n" << std::endl;}
  //    std::cout << "after " << data->iter_ << " iterations the endpoint is\n" << xAppr << std::endl;
#endif
} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * One Newton iteration of the Semi-Lagrangian Back-Tracking algorithm           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::SemiLagrange::NewtonIter(
    DRT::Element*& ele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    LINALG::Matrix<3,1>& residuum,
    LINALG::Matrix<3,1>& incr,
    double& phi,
    bool& elefound
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  LINALG::Matrix<numnode,1> shapeFcn(true);          // shape functions
  LINALG::Matrix<nsd,nsd> xji(true);        // invers of jacobian
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true); // dummy
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true); // dummy

  LINALG::Matrix<nsd,nsd> sysmat(true); // matrix for the newton system

  LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
  LINALG::Matrix<1,2*numnode> nodepres(true); // node pressures, just for function call

  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true); // dummy
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);

  // build systemMatrix
  pointdataXFEM<numnode,DISTYPE>(
      ele,
      xi,
      xji,
      shapeFcn,
      enrShapeFcnVel,
      enrShapeFcnPres,
      enrShapeXYVelDeriv1,
      enrShapeXYPresDeriv1,
      false
  );

  elementsNodalData<numnode>(
      ele,
      veln_,
      olddofman_,
      olddofcolmap_,
      oldNodalDofColDistrib_,
      nodevel,
      nodepres); // nodal data of the element nodes

  sysmat.MultiplyNT((1.0-Theta(data))*dt_,nodevel,enrShapeXYVelDeriv1); // (1-theta) * dt * v_nodes * dN/dx

  for (int i=0;i<nsd;i++)
    sysmat(i,i) += 1.0; // I + dt*velDerivXY
  sysmat.Invert();
  // invers system Matrix built

  //solve Newton iteration
  incr.Clear();
  incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

  // update iteration
  for (int i=0;i<nsd;i++)
    data->startpoint_(i) += incr(i);
  //std::cout << "in newton loop: approximate startvalue is " << data->startpoint_(0) << " " << data->startpoint_(1) << " " << data->startpoint_(2) << std::endl;

  //=============== update residuum================
  elementSearch(ele, data->startpoint_, xi, vel,phi,elefound);
} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * Computing final data where Semi-Lagrangian approach failed                    winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::getDataForNotConvergedNodes(
)
{
#ifdef PARALLEL
  if (timeIntType_!=INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew)
    exportDataToStartpointProc();
  // otherwise export is done in called function
#endif

  switch (timeIntType_)
  {
  case INPAR::COMBUST::xfemtimeint_mixedSLExtrapol:
  case INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew:
  {
    Teuchos::RCP<XFEM::TIMEINT> timeIntData = Teuchos::rcp(new XFEM::TIMEINT(
        discret_,
        olddofman_,
        newdofman_,
        oldVectors_,
        flamefront_,
        olddofcolmap_,
        newdofrowmap_,
        oldNodalDofColDistrib_,
        newNodalDofRowDistrib_,
        pbcmap_));

    Teuchos::RCP<XFEM::STD> extrapol = Teuchos::null;

    if (timeIntType_==INPAR::COMBUST::xfemtimeint_mixedSLExtrapol)
    {
      extrapol = Teuchos::rcp(new XFEM::ExtrapolationOld(
          *timeIntData,
          timeIntType_,
          veln_,
          dt_,
          flamefront_,
          veljump_,
          false));
    }
    else
    {
      extrapol = Teuchos::rcp(new XFEM::ExtrapolationNew(
          *timeIntData,
          timeIntType_,
          veln_,
          dt_,
          flamefront_,
          false));
    }

    // set state for extrapolation
    resetState(TimeIntData::failedSL_,TimeIntData::extrapolateStd_);

    // vector with data for extrapolation
    extrapol->timeIntData_ = Teuchos::rcp(new std::vector<TimeIntData>);

    // add data for extrapolation
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::extrapolateStd_)
        extrapol->timeIntData_->push_back(*data);
    }

    // clear data which is handled by extrapolation
    clearState(TimeIntData::extrapolateStd_);

    // call computation
    extrapol->compute(newVectors_);

    // add extrapolation data again in Semi-Lagrange data
    timeIntData_->insert(timeIntData_->end(),
        extrapol->timeIntData_->begin(),
        extrapol->timeIntData_->end());

    if (timeIntType_==INPAR::COMBUST::xfemtimeint_mixedSLExtrapol)
      break;

    if (timeIntType_==INPAR::COMBUST::xfemtimeint_mixedSLExtrapolNew)
    {
      bool done = true;
      for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
          data!=timeIntData_->end(); data++)
      {
        if (data->state_!=TimeIntData::doneStd_)
          done = false;
      }

      if (done==true)
        break;
      else
        exportDataToStartpointProc();
    }

    break;
  }
  case INPAR::COMBUST::xfemtimeint_semilagrange:
  {
    const int nsd = 3; // 3 dimensions for a 3d fluid element

    // remark: all data has to be sent to the processor where
    //         the startpoint lies before calling this function
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::failedSL_)
      {
        if (data->startGid_.size() != 1)
          dserror("data for alternative nodes shall be computed for one node here");

        DRT::Node* node = discret_->gNode(data->startGid_[0]); // failed node
        DRT::Element* ele = node->Elements()[0]; // element of failed node
        data->startpoint_=LINALG::Matrix<nsd,1>(node->X());

        LINALG::Matrix<nsd,1> x(node->X()); // coordinates of failed node
        LINALG::Matrix<nsd,1> xi(true); // local coordinates of failed node
        LINALG::Matrix<nsd,1> vel(true); // velocity at pseudo-Lagrangian origin
        bool elefound = false;

        /*----------------------------------------------------------*
         * element data at pseudo-Lagrangian origin                 *
         * remark: an element must be found since the intersected   *
         *         node must be on the current processor            *
         *----------------------------------------------------------*/
        callXToXiCoords(ele,x,xi,elefound);

        if (!elefound) // possibly slave node looked for element of master node or vice versa
        {
          // get pbcnode
          bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
          DRT::Node* pbcnode = NULL;
          findPBCNode(node,pbcnode,pbcnodefound);

          // get local coordinates
          LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
          callXToXiCoords(ele,pbccoords,xi,elefound);

          if (!elefound) // now something is really wrong...
            dserror("element of a row node not on same processor as node?! BUG!");
        }

        callBackTracking(ele,&*data,xi,static_cast<const char*>("failing"));
      }

    } // end loop over nodes

    break;
  }
  default:
  {
    dserror("Unknown XFEM time-integration scheme");
    break;
  }
  }
} // end getDataForNotConvergedNodes



/*------------------------------------------------------------------------------------------------*
 * call back-tracking of data at final Lagrangian origin of a point              winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::callBackTracking(
    DRT::Element*& ele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    const char* backTrackingType
)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    backTracking<numnode,DRT::Element::hex8>(ele,data,xi,backTrackingType);
  }
  break;
  case DRT::Element::hex20:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    backTracking<numnode,DRT::Element::hex20>(ele,data,xi,backTrackingType);
  }
  break;
  default:
  {
    dserror("xfem assembly type not yet implemented in time integration");
    break;
  }
  };
} // end backTracking



/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::SemiLagrange::backTracking(
    DRT::Element*& fittingele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    const char* backTrackingType
)
{
  const int nsd = 3; // dimension

  if ((strcmp(backTrackingType,static_cast<const char*>("standard"))!=0) and
      (strcmp(backTrackingType,static_cast<const char*>("failing"))!=0))
    dserror("backTrackingType not implemented");

  //  std::cout << data->node_ << "has lagrange origin " << data->startpoint_ << "with xi-coordinates "
  //      << xi << "in element " << *fittingele << std::endl;
  LINALG::Matrix<nsd,nsd> xji(true); // invers of jacobian
  LINALG::Matrix<numnode,1> shapeFcn(true); // shape function
  double deltaT = 0; // pseudo time-step size

  // enriched shape functions and there derivatives in local coordinates (N * \Psi)
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);

  // data for the final back-tracking
  LINALG::Matrix<nsd,1> vel(true); // velocity data
  std::vector<LINALG::Matrix<nsd,1> > veln(oldVectors_.size(),LINALG::Matrix<nsd,1>(true)); // velocity at t^n
  std::vector<LINALG::Matrix<nsd,nsd> > velnDeriv1(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)); // first derivation of velocity data
  LINALG::Matrix<1,1> pres(true); // pressure data
  std::vector<LINALG::Matrix<1,nsd> > presnDeriv1(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)); // first derivation of pressure data
  LINALG::Matrix<nsd,1> transportVeln(true); // transport velocity at Lagrangian origin (x_Lagr(t^n))

  int numele; // number of elements
  std::vector<const DRT::Element*> nodeeles;

  if ((data->startGid_.size() != 1) and
      (strcmp(backTrackingType,static_cast<const char*>("failing")) == 0))
    dserror("back-tracking shall be done only for one node here!");

  DRT::Node* node = discret_->gNode(data->startGid_[0]); // current node
  if (strcmp(backTrackingType,"failing")==0)
  {
    addPBCelements(node,nodeeles);
    numele=nodeeles.size();
  }
  else
    numele=1;

  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd,2*numnode> nodevel(true);
  // node velocities of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,2*numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size(),LINALG::Matrix<1,2*numnode>(true));
  std::vector<LINALG::Matrix<nsd,1> > velValues(oldVectors_.size(),LINALG::Matrix<nsd,1>(true)); // velocity of the data that should be changed
  std::vector<double> presValues(oldVectors_.size(),0); // pressures of the data that should be changed

  for (int iele=0;iele<numele;iele++) // loop over elements containing the startpoint (usually one)
  {
    for (size_t index=0;index<oldVectors_.size();index++)
    {
      nodeveldata[index].Clear();
      nodepresdata[index].Clear();
    }

    DRT::Element* ele = NULL; // current element
    if (numele>1)
    {
      ele = (DRT::Element*)nodeeles[iele];
      LINALG::Matrix<nsd,1> coords(node->X());
      LINALG::Matrix<nsd,1> vel(true); // dummy
      bool elefound = false;
      callXToXiCoords(ele,coords,xi,elefound);

      if (!elefound) // possibly slave node looked for element of master node or vice versa
      {
        // get pbcnode
        bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
        DRT::Node* pbcnode = NULL;
        findPBCNode(node,pbcnode,pbcnodefound);

        // get local coordinates
        LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
        callXToXiCoords(ele,pbccoords,xi,elefound);

        if (!elefound) // now something is really wrong...
          dserror("element of a row node not on same processor as node?! BUG!");
      }
    }
    else
      ele = fittingele;

    // evaluate data for the given point
    LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
    LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true); // dummy

    pointdataXFEM<numnode,DISTYPE>(
        ele,
        xi,
        xji,
        shapeFcn,
        enrShapeFcnVel,
        enrShapeFcnPres,
        enrShapeXYVelDeriv1,
        enrShapeXYPresDeriv1,
        false
    );

    const int* elenodeids = ele->NodeIds();

    int dofcounterVelx = 0;
    int dofcounterVely = 0;
    int dofcounterVelz = 0;
    int dofcounterPres = 0;

    for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
    {
      // get nodal velocities and pressures with help of the field set of node
      const std::set<XFEM::FieldEnr>& fieldEnrSet(olddofman_->getNodeDofSet(elenodeids[nodeid]));
      for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
          fieldenr != fieldEnrSet.end();++fieldenr)
      {
        const DofKey olddofkey(elenodeids[nodeid], *fieldenr);
        const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
        switch (fieldenr->getEnrichment().Type())
        {
        case XFEM::Enrichment::typeStandard :
        case XFEM::Enrichment::typeJump :
        case XFEM::Enrichment::typeVoid :
        case XFEM::Enrichment::typeKink :
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            if (iele==0)
              nodevel(0,dofcounterVelx) = (*veln_)[olddofcolmap_.LID(olddofpos)];
            for (size_t index=0;index<oldVectors_.size();index++)
              nodeveldata[index](0,dofcounterVelx) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            dofcounterVelx++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          {
            if (iele==0)
              nodevel(1,dofcounterVely) = (*veln_)[olddofcolmap_.LID(olddofpos)];
            for (size_t index=0;index<oldVectors_.size();index++)
              nodeveldata[index](1,dofcounterVely) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            dofcounterVely++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            if (iele==0)
              nodevel(2,dofcounterVelz) = (*veln_)[olddofcolmap_.LID(olddofpos)];
            for (size_t index=0;index<oldVectors_.size();index++)
              nodeveldata[index](2,dofcounterVelz) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            dofcounterVelz++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            for (size_t index=0;index<oldVectors_.size();index++)
              nodepresdata[index](0,dofcounterPres) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            dofcounterPres++;
          }
          else
          {
            std::cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << std::endl;
            dserror("not implemented physical field!");
          }
          break;
        }
        case XFEM::Enrichment::typeUndefined :
        default :
        {
          std::cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << std::endl;
          dserror("unknown enrichment type");
          break;
        }
        } // end switch enrichment
      } // end loop over fieldenr
    } // end loop over element nodes

    if (iele==0) // compute transportvel just once!
    {
      // interpolate velocity and pressure values at starting point
      transportVeln.Multiply(nodevel, enrShapeFcnVel);

      // computing pseudo time-step deltaT
      // remark: if x is the Lagrange-origin of node, deltaT = dt with respect to small errors.
      // if its not, deltaT estimates the time x needs to move to node)
      if (data->type_==TimeIntData::predictor_)
      {
        LINALG::Matrix<nsd,1> diff(data->node_.X());
        diff -= data->startpoint_; // diff = x_Node - x_Appr

        double numerator = transportVeln.Dot(diff); // numerator = v^T*(x_Node-x_Appr)
        double denominator = transportVeln.Dot(transportVeln); // denominator = v^T*v

        if (denominator>1e-15) deltaT = numerator/denominator; // else deltaT = 0 as initialized
      }
      else
        deltaT = dt_;
    } // end if first ele

    // interpolate velocity and pressure gradients for all fields at starting point and get final values
    for (size_t index=0;index<oldVectors_.size();index++)
    {
      if (iele==0)
        veln[index].Multiply(nodeveldata[index],enrShapeFcnVel);
      velnDeriv1[index].MultiplyNT(1.0,nodeveldata[index],enrShapeXYVelDeriv1,1.0);
      //      if (index==0)
      //        std::cout << *ele << "with nodevels " << nodeveldata[index] << ", shapefcnderiv " <<
      //            enrShapeXYVelDeriv1 << " and summed up velnderiv currently is " << velnDeriv1[index] << std::endl;

      presnDeriv1[index].MultiplyNT(1.0,nodepresdata[index],enrShapeXYPresDeriv1,1.0);
    } // end loop over vectors to be read from
  } // end loop over elements containing the point (usually one)

  for (size_t index=0;index<oldVectors_.size();index++)
  {
    velnDeriv1[index].Scale(1.0/numele);
    presnDeriv1[index].Scale(1.0/numele);

    vel.Multiply(1.0-Theta(data),velnDeriv1[index],transportVeln); // v = (1-theta)*Dv^n/Dx*v^n
    vel.Multiply(Theta(data),data->velDeriv_[index],data->vel_,1.0); // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
    vel.Update(1.0,veln[index],deltaT); // v = v_n + dt*(theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n)
    velValues[index]=vel;

    pres.Multiply(1.0-Theta(data),presnDeriv1[index],transportVeln); // p = (1-theta)*Dp^n/Dx*v^n
    pres.Multiply(Theta(data),data->presDeriv_[index],data->vel_,1.0); // p = theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n
    pres.Multiply(1.0,nodepresdata[index],enrShapeFcnPres,deltaT); // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
    presValues[index] = pres(0);
  } // loop over vectors to be set

  data->startOwner_ = std::vector<int>(1,myrank_);
  data->velValues_ = velValues;
  data->presValues_ = presValues;
  data->state_ = TimeIntData::doneStd_;
} // end backTracking



/*------------------------------------------------------------------------------------------------*
 * rewrite data for new computation                                              winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::newIteration_prepare(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors
)
{
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    data->searchedProcs_ = 1;
    data->counter_ = 0;
    data->velValues_.clear();
    data->presValues_.clear();
  }

  newIteration_nodalData(newRowVectors); // data at t^n+1 not used in predictor
  newRowVectors.clear(); // no more needed
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::newIteration_nodalData(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors
)
{
  const int nsd = 3;

  // data about column vectors required
  const Epetra_Map& newdofcolmap = *discret_->DofColMap();
  std::map<XFEM::DofKey,XFEM::DofGID> newNodalDofColDistrib;
  newdofman_->fillNodalDofColDistributionMap(newNodalDofColDistrib);

  std::vector<Teuchos::RCP<Epetra_Vector> > newColVectors;

  for (size_t index=0;index<newRowVectors.size();index++)
  {
    Teuchos::RCP<Epetra_Vector> tmpColVector = Teuchos::rcp(new Epetra_Vector(newdofcolmap,true));
    newColVectors.push_back(tmpColVector);
    LINALG::Export(*newRowVectors[index],*newColVectors[index]);
  }

  // computed data
  std::vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1(static_cast<int>(oldVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  std::vector<LINALG::Matrix<1,nsd> > presnpDeriv1(static_cast<int>(oldVectors_.size()),LINALG::Matrix<1,nsd>(true));

  std::vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1Tmp(static_cast<int>(oldVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  std::vector<LINALG::Matrix<1,nsd> > presnpDeriv1Tmp(static_cast<int>(oldVectors_.size()),LINALG::Matrix<1,nsd>(true));

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    DRT::Node& node = data->node_;

    std::vector<const DRT::Element*> eles;
    addPBCelements(&node,eles);
    const int numeles=eles.size();

    for (size_t i=0;i<newColVectors.size();i++)
    {
      velnpDeriv1[i].Clear();
      presnpDeriv1[i].Clear();
    }

    for (int iele=0;iele<numeles;iele++)
    {
      const DRT::Element* currele = eles[iele]; // current element

      switch (currele->Shape())
      {
      case DRT::Element::hex8:
      {
        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
        computeNodalGradient<numnode,DRT::Element::hex8>(newColVectors,newdofcolmap,newNodalDofColDistrib,currele,&node,velnpDeriv1Tmp,presnpDeriv1Tmp);
      }
      break;
      case DRT::Element::hex20:
      {
        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
        computeNodalGradient<numnode,DRT::Element::hex20>(newColVectors,newdofcolmap,newNodalDofColDistrib,currele,&node,velnpDeriv1Tmp,presnpDeriv1Tmp);
      }
      break;
      default:
      {
        dserror("xfem assembly type not yet implemented in time integration");
        break;
      }
      };

      for (size_t i=0;i<newColVectors.size();i++)
      {
        velnpDeriv1[i]+=velnpDeriv1Tmp[i];
        presnpDeriv1[i]+=presnpDeriv1Tmp[i];
      }
    } // end loop over elements around node

    // set transport velocity at this node
    const int gid = node.Id();
    const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(gid));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey newdofkey(gid, *fieldenr);

      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
          data->vel_(0) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
          data->vel_(1) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
          data->vel_(2) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
        }
      }
    } // end loop over fieldenr

    for (size_t i=0;i<newColVectors.size();i++)
    {
      velnpDeriv1[i].Scale(1.0/numeles);
      presnpDeriv1[i].Scale(1.0/numeles);
    }

    data->velDeriv_ = velnpDeriv1;
    data->presDeriv_ = presnpDeriv1;
    //    std::cout << "after setting transportvel is " << data->vel_ << ", velderiv is " << velnpDeriv1[0]
    //         << " and presderiv is " << presnpDeriv1[0] << std::endl;
  }
}



/*------------------------------------------------------------------------------------------------*
 * reinitialize data for new computation                                         winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::reinitializeData()
{
  const int nsd = 3; // dimension
  LINALG::Matrix<nsd,1> dummyStartpoint; // dummy startpoint for comparison
  for (int i=0;i<nsd;i++) dummyStartpoint(i) = 777.777;

  // fill curr_ structure with the data for the nodes which changed interface side
  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
  {
    DRT::Node* currnode = discret_->lColNode(lnodeid);
    // node on current processor which changed interface side
    if ((currnode->Owner() == myrank_) &&
        (interfaceSideCompare((*phinp_)[lnodeid],(*phinpi_)[lnodeid])==false))
    {
      if (interfaceSideCompare((*phinp_)[lnodeid],(*phin_)[lnodeid]) == false) // real new side
      {
        timeIntData_->push_back(TimeIntData(
            *currnode,
            LINALG::Matrix<nsd,1>(true),
            std::vector<LINALG::Matrix<nsd,nsd> >(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),
            std::vector<LINALG::Matrix<1,nsd> >(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)),
            dummyStartpoint,
            (*phinp_)[lnodeid],
            1,
            0,
            std::vector<int>(1,-1),
            std::vector<int>(1,-1),
            INFINITY,
            TimeIntData::predictor_));
      }
      else // other side than last FSI, but same side as old solution at last time step
      {
        const int nodeid = currnode->Id();

          // 1) delete data, loop backward so that deleting an vector-entry does not disturbe the iterator procedure
        for (std::vector<TimeIntData>::iterator data=timeIntData_->end()-1;
            data>=timeIntData_->begin(); data--)
        {
          if (data->node_.Id()==nodeid)
            timeIntData_->erase(data);
        }

        // 2) reset value of old solution
        // get nodal velocities and pressures with help of the field set of node
        const std::set<XFEM::FieldEnr>& fieldEnrSet(newdofman_->getNodeDofSet(nodeid));
        for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
            fieldenr != fieldEnrSet.end();++fieldenr)
        {
          const DofKey dofkey(nodeid, *fieldenr);
          const int newdofpos = newNodalDofRowDistrib_.find(dofkey)->second;
          const int olddofpos = oldNodalDofColDistrib_.find(dofkey)->second;

          switch (fieldenr->getEnrichment().Type())
          {
          case XFEM::Enrichment::typeJump :
          case XFEM::Enrichment::typeKink : break; // just standard dofs
          case XFEM::Enrichment::typeStandard :
          case XFEM::Enrichment::typeVoid :
          {
            for (size_t index=0;index<newVectors_.size();index++) // reset standard dofs due to old solution
              (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] =
                  (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            break;
          }
          case XFEM::Enrichment::typeUndefined : break;
          default :
          {
            std::cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << std::endl;
            dserror("unknown enrichment type");
            break;
          }
          } // end switch enrichment
        } // end loop over fieldenr
      }
    }
  } // end loop over processor nodes

  startpoints();

  // test loop if all initial startpoints have been computed
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->startpoint_==dummyStartpoint)
      dserror("WARNING! No enriched node on one interface side found!\nThis "
          "indicates that the whole area is at one side of the interface!");
  } // end loop over nodes
} // end function reinitializeData



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::SemiLagrange::computeNodalGradient(
    std::vector<Teuchos::RCP<Epetra_Vector> >& newColVectors,
    const Epetra_Map& newdofcolmap,
    std::map<XFEM::DofKey,XFEM::DofGID>& newNodalDofColDistrib,
    const DRT::Element* ele,
    DRT::Node* node,
    std::vector<LINALG::Matrix<3,3> >& velnpDeriv1,
    std::vector<LINALG::Matrix<1,3> >& presnpDeriv1
) const
{
  const int nsd = 3;

  for (size_t i=0;i<newColVectors.size();i++)
  {
    velnpDeriv1[i].Clear();
    presnpDeriv1[i].Clear();
  }

  // shape fcn data
  LINALG::Matrix<nsd,1> xi(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);

  // nodal data
  LINALG::Matrix<nsd,2*numnode> nodevel(true);
  LINALG::Matrix<1,2*numnode> nodepres(true);
  LINALG::Matrix<nsd,1> nodecoords(true);

  // dummies for function call
  LINALG::Matrix<nsd,nsd> jacobiDummy(true);
  LINALG::Matrix<numnode,1> shapeFcnDummy(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnDummy(true);

  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);

  { // get local coordinates
    bool elefound = false;
    LINALG::Matrix<nsd,1> coords(node->X());
    callXToXiCoords(ele,coords,xi,elefound);

    if (!elefound) // possibly slave node looked for element of master node or vice versa
    {
      // get pbcnode
      bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
      DRT::Node* pbcnode = NULL;
      findPBCNode(node,pbcnode,pbcnodefound);

      // get local coordinates
      LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
      callXToXiCoords(ele,pbccoords,xi,elefound);

      if (!elefound) // now something is really wrong...
        dserror("element of a row node not on same processor as node?! BUG!");
    }
  }

  pointdataXFEM<numnode,DISTYPE>(
      (DRT::Element*)ele,
      xi,
      jacobiDummy,
      shapeFcnDummy,
      enrShapeFcnVel,
      enrShapeFcnDummy,
      enrShapeXYVelDeriv1,
      enrShapeXYPresDeriv1,
      true);

  //std::cout << "shapefcnvel is " << enrShapeFcnVel << ", velderiv is " << enrShapeXYVelDeriv1 << " and presderiv is " << enrShapeXYPresDeriv1 << std::endl;
  for (size_t i=0;i<newColVectors.size();i++)
  {
    elementsNodalData<numnode>(
        (DRT::Element*)ele,
        newColVectors[i],
        newdofman_,
        newdofcolmap,
        newNodalDofColDistrib,
        nodevel,
        nodepres);

    velnpDeriv1[i].MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1,1.0);
    presnpDeriv1[i].MultiplyNT(1.0,nodepres,enrShapeXYPresDeriv1,1.0);
  }
} // end function compute nodal gradient



/*------------------------------------------------------------------------------------------------*
 * get the time integration factor theta fitting to the computation type         winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
double XFEM::SemiLagrange::Theta(TimeIntData* data) const
{
  double theta = -1.0;
  switch (data->type_)
  {
  case TimeIntData::predictor_: theta = 0.0; break;
  case TimeIntData::standard_ : theta = theta_default_; break;
  default: dserror("type not implemented"); break;
  }
  if (theta < 0.0) dserror("something wrong");
  return theta;
} // end function theta



/*------------------------------------------------------------------------------------------------*
 * ind processors around the current startpoint                                  winklmaier 09/13 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::findNearProcs(
    TimeIntData* data,
    const DRT::Element* const & ele
)
{
  switch (data->state_)
  {
  case TimeIntData::nextSL_:
  {
    if (data->searchedProcs_==2)
    {
      std::vector<const DRT::Element*> eles;

      // in this initial case currently procs are searched as follows:
      // initial node -> surrounding elements -> procs of their nodes
      if (data->counter_==0)
      {
        addPBCelements(&data->node_,eles);
      }
      else // counter>0 -> new iter was on curr_proc -> point in ele of curr_proc
      {
        // in this case procs are searched as follows:
        // (Lagrangian) point -> it's element -> proc's nodes -> their elements -> procs of their nodes
        eles.push_back(ele);

        for (int i=0;i<ele->NumNode();i++)
        {
          if (ele->Nodes()[i]->Owner()==myrank_)
          {
            const DRT::Node* node = ele->Nodes()[i];

            for (int j=0;j<node->NumElement();j++)
            {
              bool newproc = true;
              for (size_t k=0;k<eles.size();k++)
              {
                if (eles[k]->Owner()==node->Elements()[j]->Owner())
                {
                  newproc=false;
                  break;
                }
              }
              if (newproc) eles.push_back(node->Elements()[j]);
            }
          }
        }
      }

      if (data->node_.Owner()!=myrank_ and data->counter_==0)
        dserror("something wrong");

      for (size_t i=0;i<eles.size();i++)
      {
        const DRT::Node* const* nodes = eles[i]->Nodes();

        for (int j=0;j<eles[i]->NumNode();j++)
        {
          if ((*nodes)->Owner()!=myrank_)
          {
            bool newproc = true;
            for (std::vector<int>::iterator it=data->near_procs_.begin();it!=data->near_procs_.end();it++)
            {
              if ((*nodes)->Owner()==(*it))
              {
                newproc=false;
                break;
              }
            }

            if (newproc) data->near_procs_.push_back((*nodes)->Owner());
          }

          nodes++;
        }
      }
    } // else do-nothing
    break;
  }
  case TimeIntData::initfailedSL_:
  {
    data->near_procs_.push_back(data->startOwner_[0]);
    break;
  }
  default: break;
  }
}



/*------------------------------------------------------------------------------------------------*
 * check if newton iteration searching for the Lagrangian origin has finished    winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::SemiLagrange::globalNewtonFinished(
    int counter
) const
{
  if (counter == newton_max_iter_*numproc_)
    return true; // maximal number of iterations reached
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if ((data->state_==TimeIntData::currSL_) or
        (data->state_==TimeIntData::nextSL_) or
        (data->state_==TimeIntData::initfailedSL_))
      return false; // one node requires more data
  }
  return true; // if no more node requires data, we are done
}



# ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export alternative algo data to neighbour proc                                winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::exportDataToStartpointProc()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  std::vector<std::vector<TimeIntData> > dataVec(numproc_);

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::failedSL_)
    {
      if (data->startOwner_.size()!=1)
      {
        std::vector<int> gids = data->startGid_;
        std::vector<int> owners = data->startOwner_;

        for (size_t i=0;i<owners.size();i++)
        {
          data->startGid_.clear();
          data->startGid_.push_back(gids[i]);
          data->startOwner_.clear();
          data->startOwner_.push_back(owners[i]);
          dataVec[data->startOwner_[i]].push_back(*data);
        }
      }
      else // this case is handled explicit since it will happen (nearly) always
        dataVec[data->startOwner_[0]].push_back(*data);
    }
  }

  clearState(TimeIntData::failedSL_);
  timeIntData_->insert(timeIntData_->end(),
      dataVec[myrank_].begin(),
      dataVec[myrank_].end());

  dataVec[myrank_].clear(); // clear the set data from the vector

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    DRT::PackBuffer dataSend; // vector including all data that has to be send to dest proc

    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;

    // pack data to be sent
    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      if (data->state_==TimeIntData::failedSL_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend,data->startGid_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      if (data->state_==TimeIntData::failedSL_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend,data->startGid_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    // clear the no more needed data
    dataVec[dest].clear();

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      LINALG::Matrix<nsd,1> vel;
      std::vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      std::vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      std::vector<int> startGid;
      int newtype;

      unpackNode(posinData,dataRecv,node);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

      timeIntData_->push_back(TimeIntData(
          node,
          vel,
          velDeriv,
          presDeriv,
          startpoint,
          phiValue,
          startGid,
          std::vector<int>(1,myrank_),
          (TimeIntData::type)newtype)); // startOwner is current proc
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors
} // end exportDataToStartpointProc



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::exportIterData(
    bool& procDone
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  int dest = myrank_+1; // destination proc (the "next" one)
  if(myrank_ == (numproc_-1))
    dest = 0;

  int source = myrank_-1; // source proc (the "last" one)
  if(myrank_ == 0)
    source = numproc_-1;


  /*-------------------------------------------*
   * first part: send procfinished in order to *
   * check whether all procs have finished     *
   *-------------------------------------------*/
  for (int iproc=0;iproc<numproc_-1;iproc++)
  {
    DRT::PackBuffer dataSend;

    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));
    dataSend.StartPacking();
    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    size_t posinData = 0;
    int allProcsDone;

    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,allProcsDone);

    if (allProcsDone==0)
      procDone = 0;

    // processors wait for each other
    discret_->Comm().Barrier();
  }


  /*--------------------------------------*
   * second part: if not all procs have   *
   * finished send data to neighbour proc *
   *--------------------------------------*/
  if (!procDone)
  {
    DRT::PackBuffer dataSend;

    // fill vectors with the data
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      switch (data->state_)
      {
      case TimeIntData::initfailedSL_:
      {
        if (data->startOwner_[0]==myrank_)
        {
          data->state_ = TimeIntData::currSL_;
          DRT::Node* startnode = discret_->gNode(data->startGid_[0]);
          data->startpoint_ = LINALG::Matrix<nsd,1>(startnode->X());

          break;
        } // else -> do same as in nextsl
      }
      case TimeIntData::nextSL_:
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend,data->counter_);
        DRT::ParObject::AddtoPack(dataSend,data->startGid_);
        DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->state_);

        break;
      }
      default:
        break; // do nothing
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if ((data->state_==TimeIntData::nextSL_) or
          (data->state_==TimeIntData::initfailedSL_ and data->startOwner_[0]!=myrank_))
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend,data->counter_);
        DRT::ParObject::AddtoPack(dataSend,data->startGid_);
        DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->state_);
      }
    }

    clearState(TimeIntData::nextSL_);
    clearState(TimeIntData::initfailedSL_);

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      LINALG::Matrix<nsd,1> vel;
      std::vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      std::vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      int searchedProcs;
      int iter;
      std::vector<int> startGid;
      std::vector<int> startOwner;
      int newtype;
      int newstate;

      unpackNode(posinData,dataRecv,node);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,iter);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newstate);

      // reset data of initially failed nodes which get a second chance if on correct proc
      switch ((TimeIntData::state)newstate)
      {
      case TimeIntData::initfailedSL_:
      {
        if (startOwner[0]==myrank_)
        {
          newstate = (int)TimeIntData::currSL_;
          DRT::Node* startnode = discret_->gNode(startGid[0]);
          startpoint = LINALG::Matrix<nsd,1>(startnode->X());
        }
        break;
      }
      case TimeIntData::nextSL_:
      {
        newstate = (int)TimeIntData::currSL_;
        break;
      }
      default:
      {
        dserror("data sent to new proc which should not be sent!");
        break;
      }
      }

      timeIntData_->push_back(TimeIntData(
          node,
          vel,
          velDeriv,
          presDeriv,
          startpoint,
          phiValue,
          searchedProcs,
          iter,
          startGid,
          startOwner,
          (TimeIntData::type)newtype,
          (TimeIntData::state)newstate));
    } // end loop over number of points to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end if procfinished == false
} // end exportIterData



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::SemiLagrange::exportIterDataNew(
    bool& procDone
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  int dest = myrank_+1; // destination proc (the "next" one)
  if(myrank_ == (numproc_-1))
    dest = 0;

  int source = myrank_-1; // source proc (the "last" one)
  if(myrank_ == 0)
    source = numproc_-1;


  /*-------------------------------------------*
   * first part: send procfinished in order to *
   * check whether all procs have finished     *
   *-------------------------------------------*/
  for (int iproc=0;iproc<numproc_-1;iproc++)
  {
    DRT::PackBuffer dataSend;

    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));
    dataSend.StartPacking();
    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    size_t posinData = 0;
    int allProcsDone;

    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,allProcsDone);

    if (allProcsDone==0)
      procDone = 0;

    // processors wait for each other
    discret_->Comm().Barrier();
  }


  /*--------------------------------------*
   * second part: if not all procs have   *
   * finished send data to neighbour proc *
   *--------------------------------------*/
  if (!procDone)
  {
    std::vector<DRT::PackBuffer> dataSend(numproc_);

    // fill vectors with the data
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      int targetproc = -1;
      if (data->near_procs_.size()!=0)
        targetproc = data->near_procs_[data->near_procs_.size()-1];
      else
      {
        targetproc = myrank_+1;
        if(myrank_ == (numproc_-1))
            targetproc = 0;
      }

      switch (data->state_)
      {
      case TimeIntData::initfailedSL_:
      {
        if (data->startOwner_[0]==myrank_)
        {
          data->state_ = TimeIntData::currSL_;
          DRT::Node* startnode = discret_->gNode(data->startGid_[0]);
          data->startpoint_ = LINALG::Matrix<nsd,1>(startnode->X());

          break;
        } // else -> do same as in nextsl
      }
      case TimeIntData::nextSL_:
      {
        packNode(dataSend[targetproc],data->node_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->vel_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->counter_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->near_procs_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->startGid_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->startOwner_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],(int)data->type_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],(int)data->state_);

        break;
      }
      default:
        break; // do nothing
      }
    }

    for (int i=0;i<numproc_;i++)
      dataSend[i].StartPacking();

    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if ((data->state_==TimeIntData::nextSL_) or
          (data->state_==TimeIntData::initfailedSL_ and data->startOwner_[0]!=myrank_))
      {
        int targetproc = -1;
        if (data->near_procs_.size()!=0)
        {
          targetproc = data->near_procs_[data->near_procs_.size()-1];
          data->near_procs_.pop_back();
        }
        else
        {
          targetproc = myrank_+1;
          if(myrank_ == (numproc_-1))
              targetproc = 0;
        }

        packNode(dataSend[targetproc],data->node_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->vel_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->counter_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->near_procs_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->startGid_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],data->startOwner_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],(int)data->type_);
        DRT::ParObject::AddtoPack(dataSend[targetproc],(int)data->state_);
      }
    }

    clearState(TimeIntData::nextSL_);
    clearState(TimeIntData::initfailedSL_);

    std::vector<std::vector<char> > dataRecv(numproc_);

    for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
    {
      // Initialization
      int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
      if (source<0)       source+=numproc_;
      else if (source>=numproc_)  source-=numproc_;

      sendData(dataSend[dest],dest,source,dataRecv[source]);
      discret_->Comm().Barrier(); // processors wait for each other
    }

    for (int i=0;i<numproc_;i++)
    {
      // pointer to current position of group of cells in global std::string (counts bytes)
      std::vector<char>::size_type posinData = 0;

      // unpack received data
      while (posinData < dataRecv[i].size())
      {
        double coords[nsd] = {0.0};
        DRT::Node node(0,(double*)coords,0);
        LINALG::Matrix<nsd,1> vel;
        std::vector<LINALG::Matrix<nsd,nsd> > velDeriv;
        std::vector<LINALG::Matrix<1,nsd> > presDeriv;
        LINALG::Matrix<nsd,1> startpoint;
        double phiValue;
        int searchedProcs;
        int iter;
        std::vector<int> near_procs;
        std::vector<int> startGid;
        std::vector<int> startOwner;
        int newtype;
        int newstate;

        unpackNode(posinData,dataRecv[i],node);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],vel);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],velDeriv);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],presDeriv);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],startpoint);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],phiValue);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],searchedProcs);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],iter);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],near_procs);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],startGid);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],startOwner);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],newtype);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv[i],newstate);

        // reset data of initially failed nodes which get a second chance if on correct proc
        switch ((TimeIntData::state)newstate)
        {
        case TimeIntData::initfailedSL_:
        {
          if (startOwner[0]==myrank_)
          {
            newstate = (int)TimeIntData::currSL_;
            DRT::Node* startnode = discret_->gNode(startGid[0]);
            startpoint = LINALG::Matrix<nsd,1>(startnode->X());
          }
          else
            dserror("in this new version processor has to fit");

          break;
        }
        case TimeIntData::nextSL_:
        {
          newstate = (int)TimeIntData::currSL_;
          break;
        }
        default:
          dserror("data sent to new proc which should not be sent!");
          break;
        }

        timeIntData_->push_back(TimeIntData(
            node,
            vel,
            velDeriv,
            presDeriv,
            startpoint,
            phiValue,
            searchedProcs,
            iter,
            near_procs,
            startGid,
            startOwner,
            (TimeIntData::type)newtype,
            (TimeIntData::state)newstate));
      } // end loop over number of points to get
    }
  } // end if procfinished == false
} // end exportIterData
#endif // parallel


