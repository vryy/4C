/*!------------------------------------------------------------------------------------------------*
\file startvalues.cpp

\brief provides the SemiLagrange class

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "xfluid_timeInt_std_SemiLagrange.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../linalg/linalg_utils.H"

#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_volumecell.H"


/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrange Back-Tracking algorithm constructor                             winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_SemiLagrange::XFLUID_SemiLagrange(
    XFEM::XFLUID_TIMEINT_BASE& timeInt,
    const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,
    INPAR::XFEM::XFluidTimeInt& timeIntType,
    const RCP<Epetra_Vector> veln,
    const double& dt,
    const double& theta,
    const RCP<COMBUST::FlameFront> flamefront,
    bool initialize
) :
XFLUID_STD(timeInt, reconstr_method,
timeIntType,veln,dt,flamefront,initialize),
theta_default_(theta)
{ cout << "in constructor of XFLUID_SemiLagrange" << endl;
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrangian Back-Tracking main algorithm                                  winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::compute(
    vector<RCP<Epetra_Vector> >& newRowVectorsn
)
{cout << "in XFLUID_SemiLagrange::compute" << endl;
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  handleVectors(newRowVectorsn);

  // REMARK: in case of a new FGI iteration we have values! at new position
  //newIteration_prepare(newVectors_);

  switch (FGIType_)
  {
  case FRS1FGI1_:
  { cout << "case FRS1FGI1_" << endl;
    resetState(TimeIntData::basicStd_,TimeIntData::currSL_);
    break;
  }
  case FRSNot1_:
  { cout << "case FRSNot1_" << endl;
    resetState(TimeIntData::doneStd_,TimeIntData::currSL_);
    break;
  }
  case FRS1FGINot1_:
  { cout << "case FRS1FGINot1_" << endl;
    reinitializeData();
    resetState(TimeIntData::basicStd_,TimeIntData::currSL_);
    resetState(TimeIntData::doneStd_,TimeIntData::currSL_);
    break;
  }
  default: dserror("not implemented");
  } // end switch


  cout << "vector<TimeIntData>.size()" << timeIntData_->size() << endl;


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

      cout << "vector<TimeIntData>.size()" << timeIntData_->size() << endl;

      // loop over all nodes which changed interface side
      // remark: negative loop so that deleted elements don't influence other elements position
      for (vector<TimeIntData>::iterator data=timeIntData_->begin();
          data!=timeIntData_->end(); data++)
      {
        cout << "iter for node " << data->node_.Id() << endl;
        cout << "\t start with initial point" << data->initialpoint_ << endl;

        if (data->state_ == TimeIntData::currSL_)
        {
//          cout << endl << "on proc " << myrank_ << " iteration starts for " <<
//              data->node_ << " with initial startpoint " << data->startpoint_;

          // Initialization
          DRT::Element* ele = NULL; // pointer to the element where start point lies in
          LINALG::Matrix<nsd,1> xi(true); // local transformed coordinates of x
          LINALG::Matrix<nsd,1> vel(true); // velocity of the start point approximation
          double phin = 0.0; // phi-value of the start point approximation
          bool elefound = false;             // true if an element for a point was found on the processor
          bool stdBackTracking = true; // true if standard back-tracking shall be done

          // search for an element where the current startpoint lies in
          // if found, give out all data at the startpoint
          elementSearch(ele,data->startpoint_,xi,elefound);

          if(elefound)
          {

            bool step_np = false; // new timestep or old timestep
            getNodalDofSet(ele, data->startpoint_,data->nds_, step_np);

            LINALG::Matrix<nsd,nsd> vel_deriv_tmp(true); // dummy matrix

            // compute the velocity at startpoint
            getGPValues(ele,xi,data->nds_,step_np,vel,vel_deriv_tmp,true);

cout << "\t velocity" << vel << endl;
          }

          // if element is not found, look at another processor and so add all
          // according data to the vectors which will be sent to the next processor
          if (!elefound)
          {
            if (data->searchedProcs_ < numproc_)
            {
              data->state_ = TimeIntData::nextSL_;
              data->searchedProcs_ += 1;
            }
            else // all procs searched -> point not in domain
            {
              data->state_ = TimeIntData::failedSL_;
              cout << "WARNING! Lagrangian start point not in domain!" << endl;
            }
          } // end if elefound
          else  // if element is found, the newton iteration to find a better startpoint can start
          {
//            if (interfaceSideCompare(ele,data->startpoint_,0,data->phiValue_) == false) // Lagrangian origin and original node on different interface sides
            if (ChangedSide(ele, data->startpoint_,false, ele, data->initialpoint_, false) == true) // Lagrangian origin and original node on different interface sides
            {
              data->state_ = TimeIntData::failedSL_;
            }
            else  // Newton loop just for sensible points
            {
              NewtonLoop(ele,&*data,xi,vel,phin,elefound,stdBackTracking);
            }

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
              cout << "WARNING: newton iteration to find start value did not converge!" << endl;
              data->state_ = TimeIntData::failedSL_;
            }
          } // end if over nodes which changed interface

//          cout << "on proc " << myrank_ << " after " << data->counter_ << " iterations "
//              << data->node_ << " has startpoint " << data->startpoint_ << endl;
        }
      } // end loop over all nodes with changed interface side
    } // end if movenodes is empty or max. iter reached
    else //
      resetState(TimeIntData::currSL_,TimeIntData::failedSL_);

#ifdef PARALLEL
    // export nodes and according data for which the startpoint isn't still found (next_ vector) to next proc
    bool procDone = globalNewtonFinished();

    if(procDone) cout << " procDone on processor " << myrank_ << endl;

    exportIterData(procDone);

    // convergencecheck: procfinished == 1 just if all procs have finished
    if (procDone)
    {
      cout << "!!!!!!!!!! procDone!!!!!!!!"<< endl;
      break;
    }
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
#ifdef PARALLEL
    exportAlternativAlgoData(); // export data of failed nodes
#endif
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
    cout << "WARNING: semiLagrangeExtrapolation seems to run an infinite loop!" << endl;
#endif
#endif
} // end semiLagrangeExtrapolation


void XFEM::XFLUID_SemiLagrange::getNodalDofSet(DRT::Element* ele,
                                               LINALG::Matrix<3,1>& x ,
                                               std::vector<int>& nds,
                                               bool step_np)
{

  RCP<XFEM::FluidWizard> wizard = step_np ? wizard_new_ : wizard_old_;

  GEO::CUT::ElementHandle* e = wizard->GetElement(ele);

  if ( e!=NULL )
  {
    GEO::CUT::plain_volumecell_set cells;
    e->VolumeCells(cells);

    for(GEO::CUT::plain_volumecell_set::iterator cell_it=cells.begin(); cell_it!=cells.end(); cell_it++)
    {
      GEO::CUT::VolumeCell* cell = *cell_it;
      if(cell->Contains(x))
      {
        nds = cell->NodalDofSet();
        return;
      }
    }
  }
  else
  {
    int numnode = ele->NumNode();

    for(int inode=0; inode <numnode; inode++)
    {
      nds.push_back(0);
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::NewtonLoop(
    DRT::Element*& ele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& elefound,
    bool& stdBackTracking
)
{ cout << "NewtonLoop" << endl;
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

  // data->vel_ = vel^(n+1) for FGI>1, vel = vel^n
  residuum.Update((1.0-Theta(data)),vel,Theta(data),data->vel_); // dt*v(data->startpoint_)
  residuum.Update(1.0,data->startpoint_,-1.0,origNodeCoords,dt_);  // R = data->startpoint_ - data->movNode_ + dt*v(data->startpoint_)

  cout << "residuum " << residuum << endl;

  while(data->counter_ < newton_max_iter_)  // newton loop
  {
    data->counter_ += 1;
cout << "Newton iteration" << data->counter_ << endl;
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
      dserror("xfem assembly type not yet implemented in time integration");
    }; // end switch element type

    if (elefound) // element of data->startpoint_ at this processor
    {
      cout << "elefound" << endl;

//      if (interfaceSideCompare(ele,data->startpoint_,0,data->phiValue_) == false)
      if(ChangedSide(ele, data->startpoint_,false, ele, data->initialpoint_, false) == true)
      { cout << "changed side!!!!" << endl;
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
          {
            cout << "converged!" << endl;
            break;
          }
        }
        else
        {
          if (incr.Norm2() < relTolIncr && residuum.Norm2() < relTolRes)
          {
            cout << "converged!" << endl;
            break;
          }
        }
      } // end if interface side is the same
    } // end if elefound is true
    else // element of data->startpoint_ not at this processor
    {
      cout << "!elefound" << endl;
      break; // stop newton loop on this proc
    }

  } // end while Newton loop
#ifdef DEBUG
  // did newton iteration converge?
  if(data->counter_ == newton_max_iter_){cout << "WARNING: newton iteration for finding start value not converged for point\n" << endl;}
  //    cout << "after " << data->iter_ << " iterations the endpoint is\n" << xAppr << endl;
#endif
} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * One Newton iteration of the Semi-Lagrangian Back-Tracking algorithm           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::NewtonIter(
    DRT::Element*& ele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    LINALG::Matrix<3,1>& residuum,
    LINALG::Matrix<3,1>& incr,
    double& phi,
    bool& elefound
)
{ cout << "NewtonIter" << endl;
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  LINALG::Matrix<nsd,nsd> vel_deriv(true); // matrix for the newton system
  LINALG::Matrix<nsd,nsd> sysmat(true); // matrix for the newton system

  int step_np = false;

  // compute the velocity at startpoint
  getGPValues(ele,xi,data->nds_,step_np,vel,vel_deriv,true);


  cout << "velocity in NewtonIter" << vel << endl;

  cout << "vel_deriv in NewtonIter" << vel_deriv << endl;

  // build sysmat
  cout << "Theta" << Theta(data) << endl;
  cout << "dt" << dt_ << endl;

  sysmat.Update((1.0-Theta(data))*dt_,vel_deriv); // (1-theta) * dt * v_nodes * dN/dx

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
  cout << "in newton loop: approximate startvalue is " << data->startpoint_(0) << " " << data->startpoint_(1) << " " << data->startpoint_(2) << endl;

  //=============== update residuum================
  elementSearch(ele, data->startpoint_, xi,elefound);

  if(elefound)
  {
    bool step_np = false; // new timestep or old timestep

    getNodalDofSet(ele, data->startpoint_, data->nds_, step_np);

    LINALG::Matrix<nsd,nsd> vel_deriv_tmp(true); // dummy matrix

    // compute the velocity at startpoint
    getGPValues(ele,xi,data->nds_,step_np,vel,vel_deriv_tmp,true);

  }

} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * Computing final data where Semi-Lagrangian approach failed                    schott 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::getDataForNotConvergedNodes(
)
{
  if(true) // switch alternative algos!
  {
    const int nsd = 3; // 3 dimensions for a 3d fluid element

    // remark: all data has to be sent to the processor where
    //         the startpoint lies before calling this function
    for (vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::failedSL_)
      {

        dserror(" what to do in case that failedSL -> implement alternative algo!");

//        if (data->startGid_.size() != 1)
//          dserror("data for alternative nodes shall be computed for one node here");
//
//        DRT::Node* node = discret_->gNode(data->startGid_[0]); // failed node
//        DRT::Element* ele = node->Elements()[0]; // element of failed node
//        data->startpoint_=LINALG::Matrix<nsd,1>(node->X());
//
//        LINALG::Matrix<nsd,1> x(node->X()); // coordinates of failed node
//        LINALG::Matrix<nsd,1> xi(true); // local coordinates of failed node
//        LINALG::Matrix<nsd,1> vel(true); // velocity at pseudo-Lagrangian origin
//        bool elefound = false;
//
//        /*----------------------------------------------------------*
//         * element data at pseudo-Lagrangian origin                 *
//         * remark: an element must be found since the intersected   *
//         *         node must be on the current processor            *
//         *----------------------------------------------------------*/
//        callXToXiCoords(ele,x,xi,elefound);
//
//        if (!elefound) // possibly slave node looked for element of master node or vice versa
//        {
//          // get pbcnode
//          bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
//          DRT::Node* pbcnode = NULL;
//          findPBCNode(node,pbcnode,pbcnodefound);
//
//          // get local coordinates
//          LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
//          callXToXiCoords(ele,pbccoords,xi,elefound);
//
//          if (!elefound) // now something is really wrong...
//            dserror("element of a row node not on same processor as node?! BUG!");
//        }
//
//        callBackTracking(ele,&*data,xi,static_cast<const char*>("failing"));
      }

    } // end loop over nodes
  }












//  if (timeIntType_==INPAR::COMBUST::xfemtimeint_mixedSLExtrapol) // use Extrapolation as alternative
//  {
//    RCP<XFEM::TIMEINT> timeIntData = rcp(new XFEM::TIMEINT(
//        discret_,
//        olddofman_,
//        newdofman_,
//        oldVectors_,
//        flamefront_,
//        olddofcolmap_,
//        newdofrowmap_,
//        oldNodalDofColDistrib_,
//        newNodalDofRowDistrib_,
//        pbcmap_));
//
//    RCP<XFEM::Extrapolation> extrapol = rcp(new XFEM::Extrapolation(
//        *timeIntData,
//        timeIntType_,
//        veln_,
//        dt_,
//        flamefront_,
//        false));
//
//    // set state for extrapolation
//    resetState(TimeIntData::failedSL_,TimeIntData::extrapolateStd_);
//
//    // vector with data for extrapolation
//    extrapol->timeIntData_ = rcp(new vector<TimeIntData>);
//
//    // add data for extrapolation
//    for (vector<TimeIntData>::iterator data=timeIntData_->begin();
//        data!=timeIntData_->end(); data++)
//    {
//      if (data->state_==TimeIntData::extrapolateStd_)
//        extrapol->timeIntData_->push_back(*data);
//    }
//
//    // clear data which is handled by extrapolation
//    clearState(TimeIntData::extrapolateStd_);
//
//    // call computation
//    extrapol->compute(newVectors_);
//
//    // add extrapolation data again in Semi-Lagrange data
//    timeIntData_->insert(timeIntData_->end(),
//        extrapol->timeIntData_->begin(),
//        extrapol->timeIntData_->end());
//
//  }
//  else
//  {
//    const int nsd = 3; // 3 dimensions for a 3d fluid element
//
//    // remark: all data has to be sent to the processor where
//    //         the startpoint lies before calling this function
//    for (vector<TimeIntData>::iterator data=timeIntData_->begin();
//        data!=timeIntData_->end(); data++)
//    {
//      if (data->state_==TimeIntData::failedSL_)
//      {
//        if (data->startGid_.size() != 1)
//          dserror("data for alternative nodes shall be computed for one node here");
//
//        DRT::Node* node = discret_->gNode(data->startGid_[0]); // failed node
//        DRT::Element* ele = node->Elements()[0]; // element of failed node
//        data->startpoint_=LINALG::Matrix<nsd,1>(node->X());
//
//        LINALG::Matrix<nsd,1> x(node->X()); // coordinates of failed node
//        LINALG::Matrix<nsd,1> xi(true); // local coordinates of failed node
//        LINALG::Matrix<nsd,1> vel(true); // velocity at pseudo-Lagrangian origin
//        bool elefound = false;
//
//        /*----------------------------------------------------------*
//         * element data at pseudo-Lagrangian origin                 *
//         * remark: an element must be found since the intersected   *
//         *         node must be on the current processor            *
//         *----------------------------------------------------------*/
//        callXToXiCoords(ele,x,xi,elefound);
//
//        if (!elefound) // possibly slave node looked for element of master node or vice versa
//        {
//          // get pbcnode
//          bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
//          DRT::Node* pbcnode = NULL;
//          findPBCNode(node,pbcnode,pbcnodefound);
//
//          // get local coordinates
//          LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
//          callXToXiCoords(ele,pbccoords,xi,elefound);
//
//          if (!elefound) // now something is really wrong...
//            dserror("element of a row node not on same processor as node?! BUG!");
//        }
//
//        callBackTracking(ele,&*data,xi,static_cast<const char*>("failing"));
//      }
//
//    } // end loop over nodes
//  }
} // end getDataForNotConvergedNodes



/*------------------------------------------------------------------------------------------------*
 * call back-tracking of data at final Lagrangian origin of a point              winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::callBackTracking(
    DRT::Element*& ele,
    TimeIntData* data,
    LINALG::Matrix<3,1>& xi,
    const char* backTrackingType
)
{
  cout << "in callBackTracking" << endl;

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
    dserror("xfem assembly type not yet implemented in time integration");
  };
} // end backTracking





/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::extractNodalValuesFromVector(
    LINALG::Matrix<3,numnode>& evel,
    LINALG::Matrix<numnode,1>& epre,
    RCP<Epetra_Vector> vel_vec,
    std::vector<int>& lm
    )
{
  const int nsd = 3;
  const int numdofpernode = nsd +1;

  evel.Clear();
  epre.Clear();

  if(vel_vec == null)
    dserror("vector is null");

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*vel_vec,mymatrix,lm);

  for (int inode=0; inode<numnode; ++inode)  // number of nodes
  {
    for(int idim=0; idim<nsd; ++idim) // number of dimensions
    {
      (evel)(idim,inode) = mymatrix[idim+(inode*numdofpernode)];
    }
    (epre)(inode,0) = mymatrix[nsd+(inode*numdofpernode)];
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::backTracking(
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

    cout << data->node_ << "has lagrange origin " << data->startpoint_ << "with xi-coordinates "
        << xi << "in element " << *fittingele << endl;

  LINALG::Matrix<numnode,1> shapeFcn(true);      // shape function
  LINALG::Matrix<3,numnode> shapeFcnDeriv(true); // shape function derivatives w.r.t xyz
  LINALG::Matrix<nsd,nsd> xji(true);             // invers of jacobian

  double deltaT = 0; // pseudo time-step size


  // data for the final back-tracking
  LINALG::Matrix<nsd,1> vel(true);                                                                // velocity data
  vector<LINALG::Matrix<nsd,1> >   veln(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));          // velocity at t^n
  vector<LINALG::Matrix<nsd,nsd> > velnDeriv1(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true));  // first derivation of velocity data

  LINALG::Matrix<1,1> pres(true);                                                                 // pressure data
  vector<LINALG::Matrix<1,nsd> > presnDeriv1(oldVectors_.size(),LINALG::Matrix<1,nsd>(true));     // first derivation of pressure data
  LINALG::Matrix<nsd,1> transportVeln(true);                                                      // transport velocity at Lagrangian origin (x_Lagr(t^n))

  int numele; // number of elements
  vector<const DRT::Element*> nodeeles;

  if ((data->startGid_.size() != 1) and
      (strcmp(backTrackingType,static_cast<const char*>("failing")) == 0))
    dserror("back-tracking shall be done only for one node here!");

  DRT::Node* node = discret_->gNode(data->startGid_[0]); // current node
  if (strcmp(backTrackingType,"failing")==0)
  {

    // get all surrounding elements around node for alternative algo
    dserror("why is the lagrangian origin a node? -> alternative algo?");
    addPBCelements(node,nodeeles);
    numele=nodeeles.size();
  }
  else // standard case
    numele=1;

  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd,numnode> nodevel(true);
  LINALG::Matrix<numnode,1> nodepre(true);
  // node velocities of the element nodes for the data that should be changed
  vector<LINALG::Matrix<nsd,numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  vector<LINALG::Matrix<numnode,1> > nodepresdata(oldVectors_.size(),LINALG::Matrix<numnode,1>(true));

  vector<LINALG::Matrix<nsd,1> > velValues(oldVectors_.size(),LINALG::Matrix<nsd,1>(true)); // velocity of the data that should be changed
  vector<double> presValues(oldVectors_.size(),0); // pressures of the data that should be changed

  for (int iele=0;iele<numele;iele++) // loop over elements containing the startpoint (usually one)
  {
    for (size_t index=0;index<oldVectors_.size();index++)
    {
      nodeveldata[index].Clear();
      nodepresdata[index].Clear();
    }

    DRT::Element* ele = NULL; // current element
    if (numele>1)
    { dserror(" numele > 1 for lagrangian origin -> alternative algo?");
//      ele = (DRT::Element*)nodeeles[iele];
//      LINALG::Matrix<nsd,1> coords(node->X());
//      LINALG::Matrix<nsd,1> vel(true); // dummy
//      bool elefound = false;
//      callXToXiCoords(ele,coords,xi,elefound);
//
//      if (!elefound) // possibly slave node looked for element of master node or vice versa
//      {
//        // get pbcnode
//        bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
//        DRT::Node* pbcnode = NULL;
//        findPBCNode(node,pbcnode,pbcnodefound);
//
//        // get local coordinates
//        LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
//        callXToXiCoords(ele,pbccoords,xi,elefound);
//
//        if (!elefound) // now something is really wrong...
//          dserror("element of a row node not on same processor as node?! BUG!");
//      }
    }
    else
      ele = fittingele;



    bool compute_deriv = true;

    pointdataXFEM<numnode,DISTYPE>(
        ele,
        xi,
        xji,
        shapeFcn,
        shapeFcnDeriv,
        compute_deriv
    );

    //-------------------------------------------------------
    // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
    std::vector<int> lm;

    for(int inode=0; inode< numnode; inode++)
    {
      DRT::Node* node = ele->Nodes()[inode];
      std::vector<int> dofs;
      dofset_old_->Dof(*node, data->nds_[inode], dofs );

      int size = dofs.size();

      for (int j=0; j< size; ++j)
      {
        lm.push_back(dofs[j]);
        //cout << "lm " << lm[j] << endl;
      }

    }

    // all vectors are based on the same map

    //-------------------------------------------------------


    if (iele==0)
    {
      extractNodalValuesFromVector<numnode,DISTYPE>(nodevel,nodepre, veln_,lm);
    }
    for (size_t index=0;index<oldVectors_.size();index++)
      extractNodalValuesFromVector<numnode,DISTYPE>(nodeveldata[index],nodepresdata[index],oldVectors_[index],lm);


    cout << "element " << ele->Id() << endl;
    cout << " nodevel " << nodevel << endl;
    cout << " nodepre " << nodepre << endl;

    if (iele==0) // compute transportvel just once!
    {
      // interpolate velocity and pressure values at starting point
      transportVeln.Multiply(nodevel, shapeFcn);

      cout << "transportVeln\t" <<  transportVeln << endl;

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
      {
        veln[index].Multiply(nodeveldata[index],shapeFcn);
        cout << "veln of vector " << index << " " << veln[index] << endl;
      }
      velnDeriv1[index].MultiplyNT(1.0,nodeveldata[index],shapeFcnDeriv,1.0);

      if (index==0)
        cout << *ele << "with nodevels " << nodeveldata[index] << ", shapefcnderiv " <<
        shapeFcnDeriv << " and summed up velnderiv currently is " << velnDeriv1[index] << endl;

      presnDeriv1[index].MultiplyTT(1.0,nodepresdata[index],shapeFcnDeriv,1.0);
    } // end loop over vectors to be read from
  } // end loop over elements containing the point (usually one)


  cout << "\n!!!!!!!!  COMPUTED DATA !!!!!!!!!!" << endl;


  for (size_t index=0;index<oldVectors_.size();index++)
  {
    velnDeriv1[index].Scale(1.0/numele);
    presnDeriv1[index].Scale(1.0/numele);

    vel.Multiply(1.0-Theta(data),velnDeriv1[index],transportVeln);    // v = (1-theta)*Dv^n/Dx*v^n
    vel.Multiply(Theta(data),data->velDeriv_[index],data->vel_,1.0);  // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
    vel.Update(1.0,veln[index],deltaT);                               // v = v_n + dt*(theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n)
    velValues[index]=vel;

    pres.Multiply(1.0-Theta(data),presnDeriv1[index],transportVeln);   // p = (1-theta)*Dp^n/Dx*v^n
    pres.Multiply(Theta(data),data->presDeriv_[index],data->vel_,1.0); // p = theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n
    pres.MultiplyTN(1.0,nodepresdata[index],shapeFcn,deltaT);          // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
    presValues[index] = pres(0);

    cout << "\n!!!!!!!!  COMPUTED DATA !!!!!!!!!!" << endl;
    cout << "velocity entry in vector \t" << index << "\tn " << vel << endl;
    cout << "pressure entry in vector \t" << index << "\t " << pres(0) << endl;
  } // loop over vectors to be set

  data->startOwner_ = vector<int>(1,myrank_);
  data->velValues_ = velValues;
  data->presValues_ = presValues;
  data->state_ = TimeIntData::doneStd_;



} // end backTracking



/*------------------------------------------------------------------------------------------------*
 * rewrite data for new computation                                              winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::newIteration_prepare(
    vector<RCP<Epetra_Vector> > newRowVectors
)
{
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
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
void XFEM::XFLUID_SemiLagrange::newIteration_nodalData(
    vector<RCP<Epetra_Vector> > newRowVectors
)
{
  const int nsd = 3;

  // data about column vectors required
  const Epetra_Map& newdofcolmap = *discret_->DofColMap();
  map<XFEM::DofKey,XFEM::DofGID> newNodalDofColDistrib;
  newdofman_->fillNodalDofColDistributionMap(newNodalDofColDistrib);

  vector<RCP<Epetra_Vector> > newColVectors;

  for (size_t index=0;index<newRowVectors.size();index++)
  {
    RCP<Epetra_Vector> tmpColVector = rcp(new Epetra_Vector(newdofcolmap,true));
    newColVectors.push_back(tmpColVector);
    LINALG::Export(*newRowVectors[index],*newColVectors[index]);
  }

  // computed data
  vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1(static_cast<int>(oldVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  vector<LINALG::Matrix<1,nsd> > presnpDeriv1(static_cast<int>(oldVectors_.size()),LINALG::Matrix<1,nsd>(true));

  vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1Tmp(static_cast<int>(oldVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  vector<LINALG::Matrix<1,nsd> > presnpDeriv1Tmp(static_cast<int>(oldVectors_.size()),LINALG::Matrix<1,nsd>(true));

  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    DRT::Node& node = data->node_;

    vector<const DRT::Element*> eles;
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
        dserror("xfem assembly type not yet implemented in time integration");
      };

      for (size_t i=0;i<newColVectors.size();i++)
      {
        velnpDeriv1[i]+=velnpDeriv1Tmp[i];
        presnpDeriv1[i]+=presnpDeriv1Tmp[i];
      }
    } // end loop over elements around node

    // set transport velocity at this node
    const int gid = node.Id();
    const set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(gid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
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
    //    cout << "after setting transportvel is " << data->vel_ << ", velderiv is " << velnpDeriv1[0]
    //         << " and presderiv is " << presnpDeriv1[0] << endl;
  }
}



/*------------------------------------------------------------------------------------------------*
 * reinitialize data for new computation                                         winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::reinitializeData()
{
  cout << "in SemiLagrange::reinitializeData" << endl;
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
        timeIntData_->push_back(TimeIntData(
            *currnode,
            LINALG::Matrix<nsd,1>(true),
            vector<LINALG::Matrix<nsd,nsd> >(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),
            vector<LINALG::Matrix<1,nsd> >(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)),
            dummyStartpoint,
//            (*phinp_)[lnodeid],
            1,
            0,
            vector<int>(1,-1),
            vector<int>(1,-1),
            INFINITY,
            TimeIntData::predictor_));
      else // other side than last FSI, but same side as old solution at last time step
      {
        for (vector<TimeIntData>::iterator data=timeIntData_->begin();
            data!=timeIntData_->end(); data++)
        {
          const int nodeid = currnode->Id();

          // 1) delete data
          if (data->node_.Id()==nodeid)
            timeIntData_->erase(data);

          // 2) reset value of old solution
          // get nodal velocities and pressures with help of the field set of node
          const std::set<XFEM::FieldEnr>& fieldEnrSet(newdofman_->getNodeDofSet(nodeid));
          for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
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
              cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << endl;
              dserror("unknown enrichment type");
              break;
            }
            } // end switch enrichment
          } // end loop over fieldenr
        } // end loop over nodes
      }
    }
  } // end loop over processor nodes

  startpoints();

  // test loop if all initial startpoints have been computed
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
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
void XFEM::XFLUID_SemiLagrange::computeNodalGradient(
    vector<RCP<Epetra_Vector> >& newColVectors,
    const Epetra_Map& newdofcolmap,
    map<XFEM::DofKey,XFEM::DofGID>& newNodalDofColDistrib,
    const DRT::Element* ele,
    DRT::Node* node,
    vector<LINALG::Matrix<3,3> >& velnpDeriv1,
    vector<LINALG::Matrix<1,3> >& presnpDeriv1
) const
{ dserror("fix computeNodalGradient");
//  const int nsd = 3;
//
//  for (size_t i=0;i<newColVectors.size();i++)
//  {
//    velnpDeriv1[i].Clear();
//    presnpDeriv1[i].Clear();
//  }
//
//  // shape fcn data
//  LINALG::Matrix<nsd,1> xi(true);
//  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);
//
//  // nodal data
//  LINALG::Matrix<nsd,2*numnode> nodevel(true);
//  LINALG::Matrix<1,2*numnode> nodepres(true);
//  LINALG::Matrix<nsd,1> nodecoords(true);
//
//  // dummies for function call
//  LINALG::Matrix<nsd,nsd> jacobiDummy(true);
//  LINALG::Matrix<numnode,1> shapeFcnDummy(true);
//  LINALG::Matrix<2*numnode,1> enrShapeFcnDummy(true);
//
//#ifdef COMBUST_NORMAL_ENRICHMENT
//  LINALG::Matrix<1,numnode> nodevelenr(true);
//  ApproxFuncNormalVector<2,2*numnode> shp;
//#else
//  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
//  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
//#endif
//
//  { // get local coordinates
//    bool elefound = false;
//    LINALG::Matrix<nsd,1> coords(node->X());
//    callXToXiCoords(ele,coords,xi,elefound);
//
//    if (!elefound) // possibly slave node looked for element of master node or vice versa
//    {
//      // get pbcnode
//      bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
//      DRT::Node* pbcnode = NULL;
//      findPBCNode(node,pbcnode,pbcnodefound);
//
//      // get local coordinates
//      LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
//      callXToXiCoords(ele,pbccoords,xi,elefound);
//
//      if (!elefound) // now something is really wrong...
//        dserror("element of a row node not on same processor as node?! BUG!");
//    }
//  }
//
//#ifdef COMBUST_NORMAL_ENRICHMENT
//  pointdataXFEMNormal<numnode,DISTYPE>(
//      ele,
//#ifdef COLLAPSE_FLAME_NORMAL
//      coords,
//#endif
//      xi,
//      jacobiDummy,
//      shapeFcnDummy,
//      enrShapeFcnDummy,
//      enrShapeXYPresDeriv1,
//      shp,
//      true);
//#else
//  pointdataXFEM<numnode,DISTYPE>(
//      (DRT::Element*)ele,
//      xi,
//      jacobiDummy,
//      shapeFcnDummy,
//      enrShapeFcnVel,
//      enrShapeFcnDummy,
//      enrShapeXYVelDeriv1,
//      enrShapeXYPresDeriv1,
//      true);
//#endif
//
//  //cout << "shapefcnvel is " << enrShapeFcnVel << ", velderiv is " << enrShapeXYVelDeriv1 << " and presderiv is " << enrShapeXYPresDeriv1 << endl;
//  for (size_t i=0;i<newColVectors.size();i++)
//  {
//#ifdef COMBUST_NORMAL_ENRICHMENT
//    elementsNodalData<numnode>(
//        ele,
//        newColVectors[i],
//        dofman_,
//        newdofcolmap,
//        newNodalDofColDistrib,
//        nodevel,
//        nodevelenr,
//        nodepres);
//
//    const int* nodeids = currele->NodeIds();
//    size_t velncounter = 0;
//
//    // vderxy = enr_derxy(j,k)*evelnp(i,k);
//    for (int inode = 0; inode < numnode; ++inode)
//    {
//      // standard shape functions are identical for all vector components
//      // e.g. shp.velx.dx.s == shp.vely.dx.s == shp.velz.dx.s
//      velnpDeriv1[i](0,0) += nodevel(0,inode)*shp.velx.dx.s(inode);
//      velnpDeriv1[i](0,1) += nodevel(0,inode)*shp.velx.dy.s(inode);
//      velnpDeriv1[i](0,2) += nodevel(0,inode)*shp.velx.dz.s(inode);
//
//      velnpDeriv1[i](1,0) += nodevel(1,inode)*shp.vely.dx.s(inode);
//      velnpDeriv1[i](1,1) += nodevel(1,inode)*shp.vely.dy.s(inode);
//      velnpDeriv1[i](1,2) += nodevel(1,inode)*shp.vely.dz.s(inode);
//
//      velnpDeriv1[i](2,0) += nodevel(2,inode)*shp.velz.dx.s(inode);
//      velnpDeriv1[i](2,1) += nodevel(2,inode)*shp.velz.dy.s(inode);
//      velnpDeriv1[i](2,2) += nodevel(2,inode)*shp.velz.dz.s(inode);
//
//      const int gid = nodeids[inode];
//      const std::set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(gid);
//
//      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
//          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
//      {
//        if (enrfield->getField() == XFEM::PHYSICS::Veln)
//        {
//          velnpDeriv1[i](0,0) += nodevelenr(0,velncounter)*shp.velx.dx.n(velncounter);
//          velnpDeriv1[i](0,1) += nodevelenr(0,velncounter)*shp.velx.dy.n(velncounter);
//          velnpDeriv1[i](0,2) += nodevelenr(0,velncounter)*shp.velx.dz.n(velncounter);
//
//          velnpDeriv1[i](1,0) += nodevelenr(0,velncounter)*shp.vely.dx.n(velncounter);
//          velnpDeriv1[i](1,1) += nodevelenr(0,velncounter)*shp.vely.dy.n(velncounter);
//          velnpDeriv1[i](1,2) += nodevelenr(0,velncounter)*shp.vely.dz.n(velncounter);
//
//          velnpDeriv1[i](2,0) += nodevelenr(0,velncounter)*shp.velz.dx.n(velncounter);
//          velnpDeriv1[i](2,1) += nodevelenr(0,velncounter)*shp.velz.dy.n(velncounter);
//          velnpDeriv1[i](2,2) += nodevelenr(0,velncounter)*shp.velz.dz.n(velncounter);
//
//          velncounter += 1;
//        }
//      }
//    }
//
//#else
//    elementsNodalData<numnode>(
//        (DRT::Element*)ele,
//        newColVectors[i],
//        newdofman_,
//        newdofcolmap,
//        newNodalDofColDistrib,
//        nodevel,
//        nodepres);
//
//    velnpDeriv1[i].MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1,1.0);
//    presnpDeriv1[i].MultiplyNT(1.0,nodepres,enrShapeXYPresDeriv1,1.0);
//#endif
//  }
} // end function compute nodal gradient



/*------------------------------------------------------------------------------------------------*
 * get the time integration factor theta fitting to the computation type         winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
double XFEM::XFLUID_SemiLagrange::Theta(TimeIntData* data) const
{
  double theta = -1.0;
  switch (data->type_)
  {
  case TimeIntData::predictor_: theta = 0.0; break;
  case TimeIntData::standard_ : theta = theta_default_; break;
  default: dserror("type not implemented");
  }
  if (theta < 0.0) dserror("something wrong");
  return theta;
} // end function theta



/*------------------------------------------------------------------------------------------------*
 * check if newton iteration searching for the Lagrangian origin has finished    winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_SemiLagrange::globalNewtonFinished(
    int counter
) const
{
  if (counter == newton_max_iter_*numproc_)
    return true; // maximal number of iterations reached
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if ((data->state_==TimeIntData::currSL_) or
        (data->state_==TimeIntData::nextSL_))
    {
      return false; // one node requires more data

    }
  }
  return true; // if no more node requires data, we are done
}



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export alternative algo data to neighbour proc                                winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::exportAlternativAlgoData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  vector<vector<TimeIntData> > dataVec(numproc_);

  // fill vectors with the data
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::failedSL_)
    {
      if (data->startOwner_.size()!=1)
      {
        vector<int> gids = data->startGid_;
        vector<int> owners = data->startOwner_;

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
    for (vector<TimeIntData>::iterator data=dataVec[dest].begin();
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

    for (vector<TimeIntData>::iterator data=dataVec[dest].begin();
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

    vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      LINALG::Matrix<nsd,1> vel;
      vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      vector<int> startGid;
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
          vector<int>(1,myrank_),
          (TimeIntData::type)newtype)); // startOwner is current proc
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors
} // end exportAlternativAlgoData



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::exportIterData(
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

    vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
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
    for (vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::nextSL_)
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
      }
    }

    dataSend.StartPacking();

    for (vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::nextSL_)
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
      }
    }

    clearState(TimeIntData::nextSL_);

    vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      LINALG::Matrix<nsd,1> vel;
      vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      int searchedProcs;
      int iter;
      vector<int> startGid;
      vector<int> startOwner;
      int newtype;

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
          (TimeIntData::type)newtype));
    } // end loop over number of points to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end if procfinished == false
} // end exportIterData
#endif // parallel


