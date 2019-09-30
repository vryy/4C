/*----------------------------------------------------------------------*/
/*! \file

\brief provides the SemiLagrangean class

\level 3

\maintainer  Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_cutwizard.H"

#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_xfem/xfem_dofset.H"

#include "xfluid_timeInt_std_SemiLagrange.H"



//#define DEBUG_SEMILAGRANGE

/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrange Back-Tracking algorithm constructor                                 schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_SemiLagrange::XFLUID_SemiLagrange(
    XFEM::XFLUID_TIMEINT_BASE& timeInt,  /// time integration base class object
    const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
        reconstr_method,                      /// reconstruction map for nodes and its dofsets
    INPAR::XFEM::XFluidTimeInt& timeIntType,  /// type of time integration
    const Teuchos::RCP<Epetra_Vector> veln,   /// velocity at time t^n in col map
    const double& dt,                         /// time step size
    const double& theta,                      /// OST theta
    bool initialize                           /// is initialization?
    )
    : XFLUID_STD(timeInt, reconstr_method, timeIntType, veln, dt, initialize),
      theta_default_(theta),
      relTolIncr_(1.0e-10),
      relTolRes_(1.0e-10)
{
  return;
}  // end constructor



/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrangean Back-Tracking main algorithm                                      schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::compute(
    std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn  // row
)
{
  const int nsd = 3;  // 3 dimensions for a 3d fluid element
  handleVectors(newRowVectorsn);

  // REMARK: in case of a new FGI iteration we have values at new position
  // they are used to compute nodal gradients used for non-predictor case
  // at the end of newIteration_prepare newVectors will be cleared to fill it with new information
  // from SL-algo
  newIteration_prepare(newVectors_);

  switch (FGIType_)
  {
    case FRS1FGI1_:
    {
      if (myrank_ == 0)
        IO::cout << "\nXFLUID_SemiLagrange::compute: case FRS1FGI1_ ..." << IO::flush;
      resetState(TimeIntData::basicStd_, TimeIntData::currSL_);
      break;
    }
    case FRSNot1_:
    {
      if (myrank_ == 0)
        IO::cout << "\nXFLUID_SemiLagrange::compute: case FRSNot1_ ..." << IO::flush;
      resetState(TimeIntData::doneStd_, TimeIntData::currSL_);
      break;
    }
    case FRS1FGINot1_:
    {
      if (myrank_ == 0)
        IO::cout << "\nXFLUID_SemiLagrange::compute: case FRS1FGINot1_ ..." << IO::flush;
      reinitializeData();
      resetState(TimeIntData::basicStd_, TimeIntData::currSL_);
      resetState(TimeIntData::doneStd_, TimeIntData::currSL_);
      break;
    }
    default:
      dserror("not implemented");
      break;
  }  // end switch


#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n----------------------------------------------------------------------------------"
              "------- ";
  IO::cout << "\nReconstruct data with SEMILAGRANGEAN algorithm for " << timeIntData_->size()
           << " dofsets ";
  IO::cout << "\n----------------------------------------------------------------------------------"
              "------- "
           << IO::endl;
#endif

  /*----------------------------------------------------*
   * first part: get the correct origin for the node    *
   * in a lagrangian point of view using a newton loop  *
   *----------------------------------------------------*/
  int counter = 0;  // loop counter to avoid infinite loops

  // loop over nodes which still don't have and may get a good startvalue
  while (true)
  {
    counter += 1;

    // counter limit because maximal max_iter Newton iterations with maximal
    // numproc processor changes per iteration (avoids infinite loop)
    if (!globalNewtonFinished(counter))
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n==============================================";
      IO::cout << "\n CONTINUE GLOBAL NEWTON (" << counter << ") on proc " << myrank_;
      IO::cout << "\n==============================================" << IO::endl;
#endif

      // loop over all nodes (their std-dofsets) that have been chosen for SEMI-Lagrangean
      // reconstruction
      for (std::vector<TimeIntData>::iterator data = timeIntData_->begin();
           data != timeIntData_->end(); data++)
      {
#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t * STD-SL algorithm for node " << data->node_.Id();
#endif

        //------------------------------------
        // find element the initial startpoint lies in, if not done yet
        //------------------------------------
        if (data->initial_eid_ == -1)  // initial point not found yet
        {
          bool initial_elefound =
              false;  // true if the element for a point was found on the processor
          DRT::Element* initial_ele = NULL;  // pointer to the element where start point lies in
          LINALG::Matrix<nsd, 1> initial_xi(
              true);  // local transformed coordinates of x w.r.t found ele

          // set the element pointer where the initial point lies in!
          elementSearch(initial_ele, data->initialpoint_, initial_xi, initial_elefound);

          if (!initial_elefound)
          {
            if (data->searchedProcs_ < numproc_)
            {
              // set state to nextSL to proceed with these data on the next proc
              data->state_ = TimeIntData::nextSL_;
              data->searchedProcs_ += 1;  // increment counter that the element the point lies in
                                          // has not been found on this processor
              data->initial_eid_ = -1;
            }
            else  // all procs searched -> initial point not in domain
            {
              data->state_ = TimeIntData::failedSL_;
              dserror(
                  "<<< WARNING! Initial point for node %d for finding the Lagrangean origin not in "
                  "domain! >>>",
                  data->node_.Id());
            }
          }
          else
          {
            // check if the initial point lies in the fluid domain
            //----------------------------------------------
            // get dofset w.r.t to old interface position
            bool step_np = false;
            std::vector<int> nds_curr;
            getNodalDofSet(
                initial_ele, data->initialpoint_, nds_curr, data->last_valid_vc_, step_np);

            if (nds_curr.size() == 0) dserror("no valid nds-vector for initial point found");
            if (data->last_valid_vc_ != NULL and  // not an uncut element
                data->last_valid_vc_->Position() != GEO::CUT::Point::outside)
            {
              dserror("initial point does not lie in the fluid");
            }

            data->initial_eid_ = initial_ele->Id();
            data->startpoint_ =
                data->initialpoint_;  // start with the initial point as startpoint approximation
            data->initial_ele_owner_ = myrank_;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t -> Initial point found in element " << data->initial_eid_;
#endif
          }
        }

#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t\t -> start with start point approximation: " << data->startpoint_;
#endif

        if (data->state_ ==
            TimeIntData::currSL_)  // do not proceed when nextSL is set for current data
        {
          //------------------------------------
          // find element the initial startpoint lies in
          // if found, then get the element information
          //------------------------------------

          // Initialization
          bool elefound = false;     // true if the element for a point was found on the processor
          DRT::Element* ele = NULL;  // pointer to the element where start point lies in
          LINALG::Matrix<nsd, 1> xi(true);   // local transformed coordinates of x w.r.t found ele
          LINALG::Matrix<nsd, 1> vel(true);  // velocity of the start point approximation

          // search for an element where the current startpoint lies in
          elementSearch(ele, data->startpoint_, xi, elefound);

          //------------------------------------
          // if element is found on this proc, the newton iteration to find a better startpoint can
          // start
          //------------------------------------
          if (elefound)
          {
#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t\t ... start point approximation found in element: " << ele->Id();
#endif
            //----------------------------------------------
            DRT::Element* initial_ele = NULL;

            if (!discret_->HaveGlobalElement(data->initial_eid_))
            {
              // ChangedSide check for intersections with all sides in the boundary-dis
              // this is not so efficient but should not be called so often
              // the check itself does not need information about the background elements
            }
            else
            {
              initial_ele = discret_->gElement(data->initial_eid_);
            }

            data->changedside_ =
                ChangedSide(ele, data->startpoint_, false, initial_ele, data->initialpoint_, false);

            //----------------------------------------------
            // get dofset w.r.t to old interface position
            bool step_np = false;
            std::vector<int> nds_curr;
            getNodalDofSet(ele, data->startpoint_, nds_curr, data->last_valid_vc_, step_np);


            //=========================================================
            // how to continue if no side changing comparison possible...
            //=========================================================
            if (data->changedside_)  // how to continue in newton loop, or stop the newton loop
            {
              if (!continueForChangingSide(&*data, ele, nds_curr))
                continue;  // continue with next TimintData
            }
            else  // point did not change the side
            {
              // this set is a valid fluid set on the right side of the interface
              data->last_valid_nds_ = nds_curr;
              data->last_valid_ele_ = ele->Id();
              data->nds_ = nds_curr;
            }

            if (data->changedside_ == false and nds_curr.size() == 0)
            {
              // std::cout << "node" << data->node_.Id() << std::endl;
              // std::cout << "initial point " << data->initialpoint_ << std::endl;
              // std::cout << "current startpoint " << data->startpoint_ << std::endl;
              // std::cout << "initial element Id " << data->initial_eid_ << std::endl;
              // std::cout << "initial element owner " << data->initial_ele_owner_ << std::endl;

              dserror("point did not change the side, but nds = empty?!.");
            }

            //-------------------------------------------------
            // Newton loop just for sensible points
            //-------------------------------------------------

            //----------------------------------------------
            // compute the velocity at startpoint
            LINALG::Matrix<nsd, nsd> vel_deriv(true);  // dummy matrix for velocity derivatives
            double pres = 0.0;                         // dummy variable for pressure
            LINALG::Matrix<1, nsd> pres_deriv(true);   // dummy matrix for the pressure derivatives

            getGPValues(
                ele, xi, nds_curr, *dofset_old_, vel, vel_deriv, pres, pres_deriv, veln_, false);

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t\t ... computed velocity at start point approximation: " << vel;
#endif

            //----------------------------------------------------------------------------------------
            // call the Newton loop to get the right Lagrangean origin
            // REMARK: if newton loop is converged then return the element,
            //         local coordinates and velocity at lagrangean origin
            NewtonLoop(ele, &*data, xi, vel, elefound);
            //----------------------------------------------------------------------------------------


            // if iteration can go on (that is when startpoint is on
            // correct interface side and iter < max_iter)
            if ((data->counter_ < newton_max_iter_) and (data->state_ == TimeIntData::currSL_))
            {
              // if element is not found in a newton step, look at another processor and so add
              // all according data to the vectors which will be sent to the next processor
              if (!elefound)
              {
                // here, a point has not been found on this proc for a second time
                data->searchedProcs_ = 2;
                data->state_ = TimeIntData::nextSL_;
              }
              else
              {
                //----------------------------------------------------------------------------------------
                // newton iteration converged to a good startpoint and so the data can be used to go
                // on
                if (data->accepted_)
                {
                  callBackTracking(ele, &*data, xi, "standard");
                }
                else  // converged Lagrangean origin, however not accepted
                {
                  //----------------------------------------------------------------------
                  // a Lagrangean origin has been found but it does not lie in the fluid
                  //----------------------------------------------------------------------

                  LINALG::Matrix<nsd, 1> proj_x(true);

                  FindNearestSurfPoint(data->startpoint_, proj_x, data->last_valid_vc_, "idispn");

                  proj_x = data->initialpoint_;

                  elementSearch(ele, proj_x, xi, elefound);

                  if (elefound)
                  {
                    // check if the projected point still remains in the same element
                    if (ele->Id() == data->last_valid_ele_)
                    {
                      // TODO: check if the point lies on Boundary of the last valid vc

                      // set the new startpoint
                      data->startpoint_ = proj_x;

                      callBackTracking(ele, &*data, xi, "standard");
                    }
                    else
                    {
                      dserror(
                          "projection of startpoint lies in another element compared to the point "
                          "to be projected");
                      data->state_ = TimeIntData::failedSL_;
                    }
                  }
                  else
                    dserror(
                        "element where the projection point lies in not available on this proc");
                }
                //----------------------------------------------------------------------------------------
              }
            }  // end if
            // maximum number of iterations reached or converged origin on wrong interface side
            else if (data->counter_ == newton_max_iter_ or (!data->accepted_))
            {
              // do not use the lagrangian origin since this case is strange and potential dangerous
              data->state_ = TimeIntData::failedSL_;

#ifdef DEBUG_SEMILAGRANGE
              IO::cout << " <<< WARNING: newton iteration to find start value did not converge! >>>"
                       << IO::endl;
#endif
            }  // not converged in max_iter
          }    // if(elefound)
          // if element is not found, look at another processor and so add all
          // according data to the vectors which will be sent to the next processor
          else  // (!elefound)
          {
            if (data->searchedProcs_ < numproc_)
            {
              data->state_ = TimeIntData::nextSL_;
              data->searchedProcs_ += 1;
            }
            else  // all procs searched -> point not in domain
            {
              data->state_ = TimeIntData::failedSL_;
              IO::cout << " <<< WARNING! Lagrangian start point not in domain! >>>" << IO::endl;
            }
          }  // end if elefound
        }
      }  // end loop over all nodes stored in TiminitData marked for Semilagrangean algorithm
    }    // !globalNewtonFinished(counter)
    else
    {
      // reset the state to failed
      resetState(TimeIntData::currSL_, TimeIntData::failedSL_);
    }

    //===================================================================
    //                     PARALLEL COMMUNICATION

    // export nodes and according data for which the startpoint isn't still found (next_ vector) to
    // next proc
    bool procDone = globalNewtonFinished();

#ifdef DEBUG_SEMILAGRANGE
    if (procDone)
    {
      IO::cout << "\n==============================================";
      IO::cout << "\n FINISHED GLOBAL NEWTON on proc " << myrank_;
      IO::cout << "\n==============================================" << IO::endl;
    }
#endif

    exportIterData(procDone);

    // convergencecheck: procfinished == 1 just if all procs have finished
    if (procDone)
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n-------------------------------------------------";
      IO::cout << "\n\t\t\t !!!!!!!!!! procDone!!!!!!!!";
      IO::cout << "\n-------------------------------------------------" << IO::endl;
#endif
      break;
    }

    //===================================================================

  }  // end while loop over searched nodes


  /*-----------------------------------------------------------------------------*
   * second part: get sensible startvalues for nodes where the algorithm failed, *
   * using another algorithm, and combine the "Done" and the "Failed" - vectors  *
   *-----------------------------------------------------------------------------*/
  if (FGIType_ == FRSNot1_)  // failed nodes stay equal after the first computation
    clearState(TimeIntData::failedSL_);
  else
  {
    //===================================================================
    //                     PARALLEL COMMUNICATION
    exportAlternativAlgoData();  // export data of failed nodes
    //===================================================================

    getDataForNotConvergedNodes();  // compute final data for failed nodes
  }



  /*-----------------------------------------------------------*
   * third part: set the computed values into the state vector *
   *-----------------------------------------------------------*/
  // send the computed startvalues for every node which needs
  // new start data to the processor where the node is
  exportFinalData();

  // now every proc has the whole data for the nodes and so the data can be set to the right place
  // now
  setFinalData();

#ifdef DEBUG
  if (counter > 8 * numproc_)  // too much loops shouldnt be if all this works
    std::cout << "WARNING: semiLagrangeExtrapolation seems to run an infinite loop!" << std::endl;
#endif
}  // end semiLagrangeExtrapolation



/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm                   schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::NewtonLoop(DRT::Element*& ele,  /// pointer to element
    TimeIntData* data,                                          /// current data
    LINALG::Matrix<3, 1>& xi,                                   /// local coordinates of point
    LINALG::Matrix<3, 1>& vel,                                  /// velocity at current point
    bool& elefound                                              /// is element found ?
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t -> XFLUID_SemiLagrange::NewtonLoop" << IO::endl;
#endif

  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // Initialization
  LINALG::Matrix<nsd, 1> residuum(true);  // residuum of the newton iteration
  LINALG::Matrix<nsd, 1> incr(true);      // increment of the newton system

  // coordinates of endpoint of Lagrangian characteristics
  LINALG::Matrix<nsd, 1> origNodeCoords(true);
  for (int i = 0; i < nsd; i++) origNodeCoords(i) = data->node_.X()[i] + data->dispnp_(i);

  //-------------------------------------------------------
  // initialize residual (Theta = 0 at predictor step)
  residuum.Clear();

  // data->vel_ = vel^(n+1) for FGI>1, vel = vel^n
  residuum.Update((1.0 - Theta(data)), vel, Theta(data), data->vel_);  // dt*v(data->startpoint_)
  residuum.Update(1.0, data->startpoint_, -1.0, origNodeCoords,
      dt_);  // R = data->startpoint_ - data->node_ + dt*v(data->startpoint_)

  //==========================================================
  // (re-)start the Newton-loop on this processor
  //==========================================================
  while (data->counter_ < newton_max_iter_)  // newton loop
  {
#ifdef DEBUG_SEMILAGRANGE
    IO::cout << "\n\t\t\t NewtonLoop(" << data->counter_ << "): residuum " << residuum << IO::endl;
#endif

    data->counter_ += 1;

    //-------------------------------------
    // compute a new Newton iteration
    //-------------------------------------
    NewtonIter(ele, data, xi, residuum, incr, elefound);


    //=========================================================
    // continue on this proc if the new startpoint approximation is also on this processor
    //=========================================================
    if (elefound)  // element of data->startpoint_ at this processor
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n\t\t\t\t ... elefound " << ele->Id();
#endif

      DRT::Element* initial_ele = NULL;

      if (!discret_->HaveGlobalElement(data->initial_eid_))
      {
        // TODO: modify the ChangedSide check for intersections with all sides in the boundary-dis
        // this is not so efficient but should not be called not so often
        // the check itself does not need information about the background elements
        //              dserror("element where initial point lies in not available on proc %d, no
        //              ChangedSide comparison possible", myrank_);
      }
      else
      {
        initial_ele = discret_->gElement(data->initial_eid_);
      }

      data->changedside_ =
          ChangedSide(ele, data->startpoint_, false, initial_ele, data->initialpoint_, false);

      bool step_np = false;  // new timestep or old timestep
      std::vector<int> nds_curr;
      getNodalDofSet(ele, data->startpoint_, nds_curr, data->last_valid_vc_, step_np);

      //=========================================================
      // how to continue if no side changing comparison possible...
      //=========================================================
      if (data->changedside_)  // how to continue in newton loop, or stop the newton loop
      {
        if (!continueForChangingSide(data, ele, nds_curr)) break;
      }
      else  // point did not change the sie
      {
        // this set is a valid fluid set on the right side of the interface
        data->last_valid_nds_ = nds_curr;
        data->last_valid_ele_ = ele->Id();
        data->nds_ = nds_curr;
      }

      //=========================================================
      // continue if the semi-lagrangean origin approximation did not change the side w.r.t the
      // initial start point
      //=========================================================

      //-------------------------------------------------------
      // compute the velocity at startpoint
      LINALG::Matrix<nsd, nsd> vel_deriv(true);  // dummy matrix
      double pres = 0.0;                         // dummy variable for pressure
      LINALG::Matrix<1, nsd> pres_deriv(true);   // dummy matrix for the pressure derivatives

      getGPValues(ele, xi, nds_curr, *dofset_old_, vel, vel_deriv, pres, pres_deriv, veln_, false);


#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n\t\t\t ... computed velocity at start point approximation: " << vel;
#endif

      //-------------------------------------------------------
      // reset residual
      residuum.Clear();
      residuum.Update(
          (1.0 - Theta(data)), vel, Theta(data), data->vel_);  // dt*v(data->startpoint_)
      residuum.Update(1.0, data->startpoint_, -1.0, origNodeCoords,
          dt_);  // R = data->startpoint_ - data->movNode_ + dt*v(data->startpoint_)

      //-------------------------------------------------------
      // convergence criterion
      if (data->startpoint_.Norm2() > 1e-3)
      {
        if (incr.Norm2() / data->startpoint_.Norm2() < relTolIncr_ &&
            residuum.Norm2() / data->startpoint_.Norm2() < relTolRes_)
        {
          if (data->changedside_ == false)
          {
            data->accepted_ = true;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!";
            IO::cout << "\n\t  LAGRANGEAN ORIGIN ACCEPTED";
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }
          else
          {
            data->accepted_ = false;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!";
            IO::cout << "\n\t  LAGRANGEAN ORIGIN NOT (!!!) ACCEPTED";
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }

          break;
        }
      }
      else
      {
        if (incr.Norm2() < relTolIncr_ && residuum.Norm2() < relTolRes_)
        {
          if (data->changedside_ == false)
          {
            data->accepted_ = true;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!";
            IO::cout << "\n\t  LAGRANGEAN ORIGIN ACCEPTED";
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }
          else
          {
            data->accepted_ = false;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!";
            IO::cout << "\n\t  LAGRANGEAN ORIGIN NOT (!!!) ACCEPTED";
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }

          break;
        }
      }

    }  // end if elefound is true
    //=========================================================
    // stop Newton loop on this proc since the new startpoint approximation is not on this processor
    //=========================================================
    else  // element of data->startpoint_ not at this processor
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout
          << "\t <<< !!! element not found on this proc -> stop Newton loop on this proc !!! >>>"
          << IO::endl;
#endif
      break;  // stop newton loop on this proc
    }
  }  // end while Newton loop



#ifdef DEBUG_SEMILAGRANGE
  // did newton iteration converge?
  if (data->counter_ == newton_max_iter_)
  {
    IO::cout << "\t <<< WARNING: newton iteration for finding start value not converged for point "
                "!!! >>>"
             << IO::endl;
  }
#endif


}  // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * One Newton iteration of the Semi-Lagrangian Back-Tracking algorithm               schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::NewtonIter(DRT::Element*& ele,  /// pointer to element to be updated
    TimeIntData* data,                                          /// current data to be updated
    LINALG::Matrix<3, 1>& xi,        /// local coordinates w.r.t ele to be updated
    LINALG::Matrix<3, 1>& residuum,  /// residual for semilagrangean backtracking
    LINALG::Matrix<3, 1>& incr,      /// computed increment for lagrangean origin to be updated
    bool& elefound                   /// element found ?
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t\t ... new iteration";
#endif

  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // Initialization
  LINALG::Matrix<nsd, 1> vel_dummy(true);    // dummy matrix for the velocity
  LINALG::Matrix<nsd, nsd> vel_deriv(true);  // matrix for the velocity derivatives
  double pres_dummy = 0.0;                   // dummy variable for pressure
  LINALG::Matrix<1, nsd> pres_deriv(true);   // dummy matrix for the pressure derivatives

  LINALG::Matrix<nsd, nsd> sysmat(true);  // matrix for the newton system

  // compute the velocity derivatives at startpoint
  getGPValues(
      ele, xi, data->nds_, *dofset_old_, vel_dummy, vel_deriv, pres_dummy, pres_deriv, veln_, true);

  // build sysmat
  // JAC = I + dt(1-theta)*velDerivXY
  sysmat.Update((1.0 - Theta(data)) * dt_, vel_deriv);  // dt*(1-theta)dN/dx

  for (int i = 0; i < nsd; i++) sysmat(i, i) += 1.0;  // I + dt*velDerivXY

  // invers system Matrix built
  sysmat.Invert();


  // solve Newton iteration
  incr.Clear();
  incr.Multiply(-1.0, sysmat, residuum);  // incr = -Systemmatrix^-1 * residuum

  // update iteration
  for (int i = 0; i < nsd; i++) data->startpoint_(i) += incr(i);

#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t\t ... new approximate startvalue is " << data->startpoint_(0) << " "
           << data->startpoint_(1) << " " << data->startpoint_(2) << IO::endl;
#endif

  // find the element the new approximation lies in
  elementSearch(ele, data->startpoint_, xi, elefound);


  return;
}  // end function NewtonIter


/*------------------------------------------------------------------------------------------------*
 * check if newton iteration searching for the Lagrangian origin has finished        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_SemiLagrange::globalNewtonFinished(int counter) const
{
  if (counter == newton_max_iter_ * numproc_) return true;  // maximal number of iterations reached
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if ((data->state_ == TimeIntData::currSL_) or (data->state_ == TimeIntData::nextSL_))
    {
      return false;  // one node requires more data
    }
  }
  return true;  // if no more node requires data, we are done
}


/*------------------------------------------------------------------------------------------------*
 * Decide how or if to continue when the startpoint approximation changed the side  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_SemiLagrange::continueForChangingSide(
    TimeIntData* data,  ///< current data to be updated
    DRT::Element* ele,  ///< pointer to element the current point lies in
    std::vector<int>&
        nds_curr  ///< nds-vector of current volumecell the current startpoint approximation lies in
)
{
#if (1)
  //--------------------------------------------------------------------------------------
  // ALTERNATIVE: CONTINUE NEWTON-ALGO when startvalue changed side during newton
  // maybe the newton turns back to the right interface side


  if (nds_curr == data->last_valid_nds_ and ele->Id() == data->last_valid_ele_)
  {
    // the new newton step is within the same element and has the same nds-vector(same cell set)
    // but changed the side, then we are at the tip of a thin structure -> failed
#ifdef DEBUG_SEMILAGRANGE
    IO::cout << "\n "
                "----------------------------------------------------------------------------------"
                "-------------";
    IO::cout << "\n <<< Startpoint approximation moved within one fld-vc, but the trace intersects "
                "the side >>>";
    IO::cout << "\n                          CHANGED SIDE ";
    IO::cout
        << "\n Newton stopped! We are at the tip of a thin structure! -> leave newton loop >>>";
    IO::cout << "\n "
                "----------------------------------------------------------------------------------"
                "-------------"
             << IO::endl;
#endif
    data->state_ = TimeIntData::failedSL_;

    return false;
  }
  else if (nds_curr != data->last_valid_nds_ and ele->Id() == data->last_valid_ele_)
  {
    // the new newton step is within the same element but has a different nds-vector
    // we are within the structure or changed the side completely
    // it can be that the newton iterations goes back on the right side
    // -> continue within this element using the last valid nds-vector

    nds_curr = data->last_valid_nds_;

    return true;
  }
  else if (ele->Id() != data->last_valid_ele_)
  {
    // within the newton the element and the side have changed
#ifdef DEBUG_SEMILAGRANGE
    IO::cout << " <<< Newton for lagrangian origin can not be continued, iteration changed the "
                "side and the element! -> leave newton loop >>>"
             << IO::endl;
#endif
    data->state_ = TimeIntData::failedSL_;

    return false;
  }
  else
    dserror("case not possible");


#else
  //--------------------------------------------------------------------------------------
  // ALTERNATIVE: STOP NEWTON-ALGO when startvalue changed side during newton
#ifdef DEBUG_SEMILAGRANGE
  std::cout << "!!! CHANGED SIDE within Newton loop !!!! -> leave newton loop" << std::endl;
#endif
  data->state_ = TimeIntData::failedSL_;
  return false;  // leave newton loop if point is on wrong domain side
  //--------------------------------------------------------------------------------------
#endif

  return false;
}


/*------------------------------------------------------------------------------------------------*
 * Computing final data where Semi-Lagrangian approach failed                        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::getDataForNotConvergedNodes()
{
  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // remark: all data have to be sent to the processor where
  //         the startpoint lies before calling this function
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->state_ == TimeIntData::failedSL_)
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "P " << myrank_ << " WARNING: failedSL -> alternative algo!" << IO::endl;
      IO::cout << "P " << myrank_ << " node " << data->node_.Id() << IO::endl;
      IO::cout << "P " << myrank_ << " use initial point: " << data->initialpoint_ << IO::endl;
#endif

      // Initialization
      DRT::Element* ele = NULL;  // pointer to the element where pseudo-Lagrangian origin lies in
      LINALG::Matrix<nsd, 1> xi(true);   // local coordinates of pseudo-Lagrangian origin
      LINALG::Matrix<nsd, 1> vel(true);  // velocity at pseudo-Lagrangian origin
      bool elefound = false;  // true if an element for a point was found on the processor

      // search for an element where the current startpoint lies in
      elementSearch(ele, data->initialpoint_, xi, elefound);

      // if found, give out all data at the startpoint
      if (elefound)
      {
        bool step_np = false;  // data w.r.t old interface position
        getNodalDofSet(ele, data->initialpoint_, data->nds_, data->last_valid_vc_, step_np);
      }
      else  // possibly slave node looked for element of master node or vice versa
      {
        dserror("element not found");
      }

      //-------------------------------------------------------------------
      // call the back Tracking computation based on the initial point
      // which is a rough approximation of the lagrangian origin
      callBackTracking(ele, &*data, xi, static_cast<const char*>("failing"));
      //-------------------------------------------------------------------

    }  // if(failedSL_)
  }    // end loop over nodes

  return;
}  // end getDataForNotConvergedNodes


/*------------------------------------------------------------------------------------------------*
 * rewrite data for new computation                                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::newIteration_prepare(
    std::vector<Teuchos::RCP<Epetra_Vector>> newRowVectors)
{
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    data->searchedProcs_ = 1;
    data->counter_ = 0;
    data->velValues_.clear();
    data->presValues_.clear();
  }

  newIteration_nodalData(newRowVectors);  // data at t^n+1 not used in predictor
  newRowVectors.clear();                  // no more needed
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                          schott 04/13 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::newIteration_nodalData(
    std::vector<Teuchos::RCP<Epetra_Vector>> newRowVectors)
{
  const int nsd = 3;

  std::vector<Teuchos::RCP<Epetra_Vector>> newColVectors;

  for (size_t index = 0; index < newRowVectors.size(); index++)
  {
    const Epetra_Map* newdofcolmap = discret_->DofColMap();

    Teuchos::RCP<Epetra_Vector> tmpColVector = Teuchos::rcp(new Epetra_Vector(*newdofcolmap, true));
    newColVectors.push_back(tmpColVector);
    LINALG::Export(*newRowVectors[index], *newColVectors[index]);
  }

  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->type_ == TimeIntData::predictor_)
      continue;  // no info at new interface position required

    if (data->type_ != TimeIntData::predictor_)
      dserror(
          "this function is used for the first time here, check the compute gradient "
          "functionality, check also the parallel export!");

    DRT::Node* node = &data->node_;

    //----------------------------------------------------------
    // Reconstruct nodal gradients for all vectors
    //----------------------------------------------------------

    // node velocities of the element nodes for the fields that should be changed
    std::vector<LINALG::Matrix<nsd, nsd>> avg_nodevelgraddata(
        newColVectors.size(), LINALG::Matrix<nsd, nsd>(true));

    // node pressures of the element nodes for the data that should be changed
    std::vector<LINALG::Matrix<1, nsd>> avg_nodepresgraddata(
        newColVectors.size(), LINALG::Matrix<1, nsd>(true));


    // determine the elements used for the nodal gradient computation
    // determine the corresponding nodal dofset vectors used for averaging the nodal gradients
    std::vector<DRT::Element*> eles_avg;
    std::vector<std::vector<int>> eles_avg_nds;


    GEO::CUT::Node* n = wizard_new_->GetNode(node->Id());

    if (n != NULL)
    {
      // get the nodal dofset w.r.t the Lagrangean origin
      const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cellset =
          n->GetNodalDofSet(0)->VolumeCellComposite();  // always the standard dofset

      // get for each adjacent element contained in the nodal dofset the first vc, its parent
      // element is used for the reconstruction REMARK: adjacent elements for that no elementhandle
      // exists in the cut won't be found here
      for (std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>::const_iterator cellset_it =
               cellset.begin();
           cellset_it != cellset.end(); cellset_it++)
      {
        // the first vc representing the set
        GEO::CUT::VolumeCell* vc = *((*cellset_it).begin());
        int peid = vc->ParentElement()->GetParentId();

        if (!discret_->HaveGlobalElement(peid))
          dserror("element %d for averaging not on proc %d", peid, myrank_);

        // get the element
        DRT::Element* e = discret_->gElement(peid);
        // get the vc's nds vector
        const std::vector<int> e_nds = vc->NodalDofSet();

        // add the element and the nds vector
        eles_avg.push_back(e);
        eles_avg_nds.push_back(e_nds);
      }
    }

    // check all surrounding elements
    int numele = node->NumElement();
    DRT::Element** eles = node->Elements();

    // add surrounding std uncut elements for that no elementhandle is available
    for (int i = 0; i < numele; i++)
    {
      GEO::CUT::ElementHandle* eh = wizard_new_->GetElement(eles[i]);

      if (eh != NULL)
        continue;  // element and the right nds-vec should have been found using the for-loop before

      // if we are here, then the element is a standard uncut element
      // and it is ensured that it has not been added to eles_avg yet
      std::vector<int> std_nds(eles[i]->NumNode(), 0);

      eles_avg.push_back(eles[i]);
      eles_avg_nds.push_back(std_nds);
    }

    if (eles_avg.size() == 0) dserror("there is no element for averaging");

    // compute the nodal gradients for velocity/acc and pressure component
    computeNodalGradient(newColVectors, node, eles_avg, eles_avg_nds, *dofset_new_,
        avg_nodevelgraddata, avg_nodepresgraddata);


    data->velDeriv_ = avg_nodevelgraddata;
    data->presDeriv_ = avg_nodepresgraddata;

    //----------------------------------------------------------
    // get the transport velocity at t^(n+1) at the node
    //----------------------------------------------------------

    //-------------------------------------------------------
    // get node location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
    std::vector<int> lm;
    std::vector<int> dofs;

    dofset_new_->Dof(dofs, node, 0);  // dofs for standard dofset

    for (int j = 0; j < 4; ++j)
    {
      lm.push_back(dofs[j]);
    }

    //-------------------------------------------------------
    // the first vector contains the velocity information
    LINALG::Matrix<3, 1> nodevel(true);
    LINALG::Matrix<1, 1> nodepre(true);
    extractNodalValuesFromVector<1>(nodevel, nodepre, newColVectors[0], lm);

    data->vel_ = nodevel;

    LINALG::Matrix<3, 1> nodedispnp(true);
    if (dispnp_ != Teuchos::null)  // is alefluid
    {
      //------------------------------------------------------- add ale disp
      // get node location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)

      LINALG::Matrix<1, 1> nodepredummy(true);
      extractNodalValuesFromVector<1>(nodedispnp, nodepredummy, dispnp_, lm);
    }

    data->dispnp_ = nodedispnp;
  }
  //----------------------------------------------------------------------------------
}



/*------------------------------------------------------------------------------------------------*
 * reinitialize data for new computation                                         winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::reinitializeData()
{
  dserror("adapt implementation of this function");
  dserror("adapt, how to get nds_np?");

  //  int nds_np = -1;
  //
  //  std::cout << "in SemiLagrange::reinitializeData" << std::endl;
  //  const int nsd = 3; // dimension
  //  LINALG::Matrix<nsd,1> dummyStartpoint; // dummy startpoint for comparison
  //  for (int i=0;i<nsd;i++) dummyStartpoint(i) = 777.777;
  //
  //  // fill curr_ structure with the data for the nodes which changed interface side
  //  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor
  //  nodes
  //  {
  //    DRT::Node* currnode = discret_->lColNode(lnodeid);
  //    // node on current processor which changed interface side
  //    if ((currnode->Owner() == myrank_) &&
  //        (interfaceSideCompare((*phinp_)[lnodeid],(*phinpi_)[lnodeid])==false))
  //    {
  //      if (interfaceSideCompare((*phinp_)[lnodeid],(*phin_)[lnodeid]) == false) // real new side
  //        timeIntData_->push_back(TimeIntData(
  //            *currnode,
  //            nds_np,
  //            LINALG::Matrix<nsd,1>(true),
  //            std::vector<LINALG::Matrix<nsd,nsd>
  //            >(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),
  //            std::vector<LINALG::Matrix<1,nsd> >(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)),
  //            dummyStartpoint,
  ////            (*phinp_)[lnodeid],
  //            1,
  //            0,
  //            std::vector<int>(1,-1),
  //            std::vector<int>(1,-1),
  //            INFINITY,
  //            TimeIntData::predictor_));
  //      else // other side than last FSI, but same side as old solution at last time step
  //      {
  //        for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
  //            data!=timeIntData_->end(); data++)
  //        {
  //          const int nodeid = currnode->Id();
  //
  //          // 1) delete data
  //          if (data->node_.Id()==nodeid)
  //            timeIntData_->erase(data);
  //
  //          // 2) reset value of old solution
  //          // get nodal velocities and pressures with help of the field set of node
  //          const std::set<XFEM::FieldEnr>& fieldEnrSet(newdofman_->getNodeDofSet(nodeid));
  //          for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
  //              fieldenr != fieldEnrSet.end();++fieldenr)
  //          {
  //            const DofKey dofkey(nodeid, *fieldenr);
  //            const int newdofpos = newNodalDofRowDistrib_.find(dofkey)->second;
  //            const int olddofpos = oldNodalDofColDistrib_.find(dofkey)->second;
  //            switch (fieldenr->getEnrichment().Type())
  //            {
  //            case XFEM::Enrichment::typeJump :
  //            case XFEM::Enrichment::typeKink : break; // just standard dofs
  //            case XFEM::Enrichment::typeStandard :
  //            case XFEM::Enrichment::typeVoid :
  //            {
  //              for (size_t index=0;index<newVectors_.size();index++) // reset standard dofs due
  //              to old solution
  //                (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] =
  //                    (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
  //              break;
  //            }
  //            case XFEM::Enrichment::typeUndefined : break;
  //            default :
  //            {
  //              std::cout <<
  //              fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) <<
  //              std::endl; dserror("unknown enrichment type"); break;
  //            }
  //            } // end switch enrichment
  //          } // end loop over fieldenr
  //        } // end loop over nodes
  //      }
  //    }
  //  } // end loop over processor nodes
  //
  //  startpoints();
  //
  //  // test loop if all initial startpoints have been computed
  //  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
  //      data!=timeIntData_->end(); data++)
  //  {
  //    if (data->startpoint_==dummyStartpoint)
  //      dserror("WARNING! No enriched node on one interface side found!\nThis "
  //          "indicates that the whole area is at one side of the interface!");
  //  } // end loop over nodes
}  // end function reinitializeData


/*------------------------------------------------------------------------------------------------*
 * call back-tracking of data at final Lagrangian origin of a point                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::callBackTracking(DRT::Element*& ele,  /// pointer to element
    TimeIntData* data,                                                /// data
    LINALG::Matrix<3, 1>& xi,                                         /// local coordinates
    const char* backTrackingType                                      /// type of backTracking
)
{
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      backTracking<numnode, DRT::Element::hex8>(ele, data, xi, backTrackingType);
    }
    break;
    case DRT::Element::hex20:
    {
      const int numnode =
          DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
      backTracking<numnode, DRT::Element::hex20>(ele, data, xi, backTrackingType);
    }
    break;
    default:
      dserror("xfem assembly type not yet implemented in time integration");
      break;
  };
}  // end backTracking


/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                       schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
template <const int numnode, DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::backTracking(DRT::Element*& fittingele,  /// pointer to element
    TimeIntData* data,                                                   /// data
    LINALG::Matrix<3, 1>& xi,                                            /// local coordinates
    const char* backTrackingType                                         /// type of backTracking
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n==============================================";
  IO::cout << "\n BACK-TRACKING on proc " << myrank_;
  IO::cout << "\n==============================================" << IO::endl;
#endif


  const int nsd = 3;  // dimension

  if ((strcmp(backTrackingType, static_cast<const char*>("standard")) != 0) and
      (strcmp(backTrackingType, static_cast<const char*>("failing")) != 0))
    dserror("backTrackingType not implemented");

#ifdef DEBUG_SEMILAGRANGE
  if (strcmp(backTrackingType, static_cast<const char*>("standard")) == 0)
  {
    std::cout << "\n--------------------------------------------------\n"
              << "\nnode: " << data->node_ << "\ncomputed LAGRANGEAN ORIGIN  (startpoint) "
              << data->startpoint_ << "with xi-coord. " << xi << "in element " << *fittingele
              << "\n--------------------------------------------------" << std::endl;
  }
  if (strcmp(backTrackingType, static_cast<const char*>("failing")) == 0)
  {
    std::cout << "\n--------------------------------------------------\n"
              << "\nnode: " << data->node_
              << "\nused <<<PSEUDO>>> LAGRANGEAN ORIGIN (initialpoint) " << data->initialpoint_
              << "with xi-coord. " << xi << "in element " << *fittingele
              << "\n--------------------------------------------------" << std::endl;
  }
#endif

  //---------------------------------------------------------------------------------
  // Initialization

  LINALG::Matrix<3, 1> lagrangeanOrigin(
      true);  // the applied Lagrangean origin (the real computed or an approximated)

  if (strcmp(backTrackingType, static_cast<const char*>("standard")) == 0)
  {
    lagrangeanOrigin = data->startpoint_;  // use the computed start point approximation or the
                                           // projected start point
  }
  else if (strcmp(backTrackingType, static_cast<const char*>("failing")) == 0)
  {
    lagrangeanOrigin = data->initialpoint_;  // use the initial guess for the Lagrangean origin

    IO::cout << "\n\tWARNING: Proc " << myrank_
             << ": SEMI-LAGRANGEAN algorithm: used the initial-guess instead of the real "
                "Lagrangean origin for node "
             << data->node_.Id() << IO::endl;
  }
  else
    dserror("backTrackingType not implemented");


  LINALG::Matrix<numnode, 1> shapeFcn(true);       // shape function
  LINALG::Matrix<3, numnode> shapeFcnDeriv(true);  // shape function derivatives w.r.t xyz
  LINALG::Matrix<nsd, nsd> xji(true);              // invers of jacobian

  double deltaT = 0;  // pseudo time-step size, used when the initial point is used instead of the
                      // computed lagrangean startpoint

  // data for the final back-tracking
  LINALG::Matrix<nsd, 1> vel(true);  // velocity data
  std::vector<LINALG::Matrix<nsd, nsd>> velnDeriv1(
      oldVectors_.size(), LINALG::Matrix<nsd, nsd>(true));  // first derivation of velocity data

  LINALG::Matrix<1, 1> pres(true);  // pressure data
  std::vector<LINALG::Matrix<1, nsd>> presnDeriv1(
      oldVectors_.size(), LINALG::Matrix<1, nsd>(true));  // first derivation of pressure data

  std::vector<LINALG::Matrix<nsd, 1>> veln(
      oldVectors_.size(), LINALG::Matrix<nsd, 1>(true));  // velocity at t^n
  LINALG::Matrix<nsd, 1> transportVeln(
      true);  // transport velocity at Lagrangian origin (x_Lagr(t^n))


  //---------------------------------------------------------------------------------
  // fill velocity and pressure data at nodes of element ...

  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd, numnode> nodevel(true);
  LINALG::Matrix<numnode, 1> nodepre(true);

  // node velocities of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<nsd, numnode>> nodeveldata(
      oldVectors_.size(), LINALG::Matrix<nsd, numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<numnode, 1>> nodepresdata(
      oldVectors_.size(), LINALG::Matrix<numnode, 1>(true));

  // velocity of the data that shall be changed
  std::vector<LINALG::Matrix<nsd, 1>> velValues(oldVectors_.size(), LINALG::Matrix<nsd, 1>(true));
  // pressures of the data that shall be changed
  std::vector<double> presValues(oldVectors_.size(), 0);


  for (size_t index = 0; index < oldVectors_.size(); index++)
  {
    nodeveldata[index].Clear();
    nodepresdata[index].Clear();
  }

  DRT::Element* ele = fittingele;  // current element

  //---------------------------------------------------------------------------------
  // get shape functions and derivatives at local coordinates

  bool compute_deriv = true;

  evalShapeAndDeriv<numnode, DISTYPE>(ele, xi, xji, shapeFcn, shapeFcnDeriv, compute_deriv);

  //-------------------------------------------------------
  // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
  std::vector<int> lm;

  for (int inode = 0; inode < numnode; inode++)
  {
    DRT::Node* node = ele->Nodes()[inode];
    std::vector<int> dofs;
    dofset_old_->Dof(dofs, node, data->nds_[inode]);


    for (int j = 0; j < 4; ++j)
    {
      lm.push_back(dofs[j]);
    }
  }

  //-------------------------------------------------------
  // all vectors are based on the same map

  extractNodalValuesFromVector<numnode>(nodevel, nodepre, veln_, lm);

  for (size_t index = 0; index < oldVectors_.size(); index++)
    extractNodalValuesFromVector<numnode>(
        nodeveldata[index], nodepresdata[index], oldVectors_[index], lm);



  //----------------------------------------------------------------------------------

  // node velocities of the element nodes for the fields that should be changed
  std::vector<std::vector<LINALG::Matrix<nsd, nsd>>> avg_nodevelgraddata;
  avg_nodevelgraddata.reserve(numnode);

  // node pressures of the element nodes for the data that should be changed
  std::vector<std::vector<LINALG::Matrix<1, nsd>>> avg_nodepresgraddata;
  avg_nodepresgraddata.reserve(numnode);

  for (size_t i = 0; i < numnode; i++)
  {
    std::vector<LINALG::Matrix<nsd, nsd>> tmp_vec(
        oldVectors_.size(), LINALG::Matrix<nsd, nsd>(true));
    avg_nodevelgraddata.push_back(tmp_vec);

    std::vector<LINALG::Matrix<1, nsd>> tmp_vec2(oldVectors_.size(), LINALG::Matrix<1, nsd>(true));
    avg_nodepresgraddata.push_back(tmp_vec2);
  }

  // Reconstruct nodal gradients
  for (int inode = 0; inode < numnode; inode++)
  {
    DRT::Node* node = (ele->Nodes())[inode];

    // determine the elements used for the nodal gradient computation
    // determine the corresponding nodal dofset vectors used for averaging the nodal gradients
    std::vector<DRT::Element*> eles_avg;
    std::vector<std::vector<int>> eles_avg_nds;


    GEO::CUT::Node* n = wizard_old_->GetNode(node->Id());

    if (n != NULL)
    {
      // get the nodal dofset w.r.t the Lagrangean origin
      const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cellset =
          n->GetNodalDofSet(data->nds_[inode])->VolumeCellComposite();


      // get for each adjacent element contained in the nodal dofset the first vc, its parent
      // element is used for the reconstruction REMARK: adjacent elements for that no elementhandle
      // exists in the cut won't be found here
      for (std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>::const_iterator cellset_it =
               cellset.begin();
           cellset_it != cellset.end(); cellset_it++)
      {
        // the first vc representing the set
        GEO::CUT::VolumeCell* vc = *((*cellset_it).begin());
        int peid = vc->ParentElement()->GetParentId();

        if (!discret_->HaveGlobalElement(peid))
          dserror("element %d for averaging not on proc %d", peid, myrank_);

        // get the element
        DRT::Element* e = discret_->gElement(peid);
        // get the vc's nds vector
        const std::vector<int> e_nds = vc->NodalDofSet();

        // add the element and the nds vector
        eles_avg.push_back(e);
        eles_avg_nds.push_back(e_nds);
      }
    }

    // check all surrounding elements
    int numele = node->NumElement();
    DRT::Element** eles = node->Elements();

    // add surrounding std uncut elements for that no elementhandle is available
    for (int i = 0; i < numele; i++)
    {
      GEO::CUT::ElementHandle* eh = wizard_old_->GetElement(eles[i]);

      if (eh != NULL)
        continue;  // element and the right nds-vec should have been found using the for-loop before

      // if we are here, then the element is a standard uncut element
      // and it is ensured that it has not been added to eles_avg yet
      std::vector<int> std_nds(eles[i]->NumNode(), 0);

      eles_avg.push_back(eles[i]);
      eles_avg_nds.push_back(std_nds);
    }

    if (eles_avg.size() == 0) dserror("there is no element for averaging");

    // compute the nodal gradients for velocity/acc and pressure component
    computeNodalGradient(oldVectors_, node, eles_avg, eles_avg_nds, *dofset_old_,
        avg_nodevelgraddata[inode], avg_nodepresgraddata[inode]);
  }
  //----------------------------------------------------------------------------------


  //---------------------------------------------------------------------------------
  // interpolate velocity and pressure values at starting point
  transportVeln.Multiply(nodevel, shapeFcn);

#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\t transportVeln\t" << transportVeln << IO::endl;
#endif

  //---------------------------------------------------------------------------------
  // computing pseudo time-step deltaT
  // remark: if x is the Lagrange-origin of node, deltaT = dt with respect to small errors.
  // if it is not, deltaT estimates the time x needs to move to node)
  if (data->type_ == TimeIntData::predictor_)
  {
    LINALG::Matrix<nsd, 1> diff(data->node_.X());
    for (int i = 0; i < nsd; ++i) diff(i) += data->dispnp_(i);
    diff -= lagrangeanOrigin;  // diff = x_Node - x_Appr

    double numerator = transportVeln.Dot(diff);             // numerator = v^T*(x_Node-x_Appr)
    double denominator = transportVeln.Dot(transportVeln);  // denominator = v^T*v

    if (denominator > 1e-15) deltaT = numerator / denominator;  // else deltaT = 0 as initialized

#ifdef DEBUG_SEMILAGRANGE
    IO::cout << " \t recomputed modified pseudo time-step size: " << deltaT << IO::endl;
#endif
  }
  else
    deltaT = dt_;


  // interpolate velocity and pressure gradients for all fields at starting point and get final
  // values
  for (size_t index = 0; index < oldVectors_.size(); index++)
  {
    veln[index].Multiply(nodeveldata[index], shapeFcn);

#if (0)  // use the non-averaged gradients interpolated within the element
    velnDeriv1[index].MultiplyNT(1.0, nodeveldata[index], shapeFcnDeriv, 1.0);
    presnDeriv1[index].MultiplyTT(1.0, nodepresdata[index], shapeFcnDeriv, 1.0);
#else  // use the averaged nodal gradients
    velnDeriv1[index].Clear();

    for (int i = 0; i < numnode; i++)
    {
      velnDeriv1[index].Update(shapeFcn(i), avg_nodevelgraddata[i][index], 1.0);
      presnDeriv1[index].Update(shapeFcn(i), avg_nodepresgraddata[i][index], 1.0);
    }

#endif
  }  // end loop over vectors to be read from


  for (size_t index = 0; index < oldVectors_.size(); index++)
  {
    vel.Multiply(1.0 - Theta(data), velnDeriv1[index], transportVeln);  // v = (1-theta)*Dv^n/Dx*v^n
    vel.Multiply(Theta(data), data->velDeriv_[index], data->vel_,
        1.0);  // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
    vel.Update(
        1.0, veln[index], deltaT);  // v = v_n + dt*(theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n)
    velValues[index] = vel;

    pres.Multiply(
        1.0 - Theta(data), presnDeriv1[index], transportVeln);  // p = (1-theta)*Dp^n/Dx*v^n
    pres.Multiply(Theta(data), data->presDeriv_[index], data->vel_,
        1.0);  // p = theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n
    pres.MultiplyTN(1.0, nodepresdata[index], shapeFcn,
        deltaT);  // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
    presValues[index] = pres(0);

#ifdef DEBUG_SEMILAGRANGE
    IO::cout << "\n***********************************************";
    IO::cout << "\n           RECONSTRUCTED VALUES for node " << (data->node_).Id();
    IO::cout << "\nvelocity entry in vector \t" << index << "\t " << vel;
    IO::cout << "pressure entry in vector \t" << index << "\t " << pres(0);
    IO::cout << "\n***********************************************" << IO::endl;
#endif
  }  // loop over vectors to be set

  data->velValues_ = velValues;
  data->presValues_ = presValues;
  data->state_ = TimeIntData::doneStd_;


  return;
}  // end backTracking


/*------------------------------------------------------------------------------------------------*
 * determine point's dofset in element ele w.r.t old or new interface position       schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::getNodalDofSet(DRT::Element* ele,  /// pointer to element
    LINALG::Matrix<3, 1>& x,                                       /// global coordinates of point
    std::vector<int>& nds,      /// determine the points dofset w.r.t old/new interface position
    GEO::CUT::VolumeCell*& vc,  /// valid fluid volumecell the point x lies in
    bool step_np                /// computation w.r.t old or new interface position?
)
{
  nds.clear();


#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t ... getNodalDofSet";
#endif


  Teuchos::RCP<GEO::CutWizard> wizard = step_np ? wizard_new_ : wizard_old_;

  GEO::CUT::ElementHandle* e = wizard->GetElement(ele);

  bool inside_structure = false;

  if (e != NULL)  // element in cut involved
  {
    GEO::CUT::plain_volumecell_set cells;
    e->GetVolumeCells(cells);

    if (cells.size() == 0)
      dserror("GEO::CUT::Element %d does not contain any volume cell", ele->Id());

    for (GEO::CUT::plain_volumecell_set::iterator cell_it = cells.begin(); cell_it != cells.end();
         cell_it++)
    {
      GEO::CUT::VolumeCell* cell = *cell_it;
      //      if(cell->Contains(x))
      // cell contains the point inside or on one of its boundaries and the cell is an outside
      // (fluid) cell
      if (((cell->IsThisPointInside(x) == "inside") or
              (cell->IsThisPointInside(x) == "onBoundary")) and
          cell->Position() == GEO::CUT::Point::outside)
      {
#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t\t\t -> Position of point w.r.t volumecell is "
                 << cell->IsThisPointInside(x) << " \t cell pos = " << cell->Position() << IO::endl;
#endif
        nds = cell->NodalDofSet();

        if ((int)nds.size() != ele->NumNode())
          dserror(
              " size of nds-vector != number of element nodes, why does the vc does not have nodal "
              "dofsets?!");

        vc = cell;

#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "nds-vector ";
        for (int i = 0; i < (int)nds.size(); i++)
        {
          IO::cout << " " << nds[i];
        }
        IO::cout << "\n";
#endif

        return;
      }
      // point lies within the structure or on Boundary
      else if ((cell->IsThisPointInside(x) == "inside" or
                   cell->IsThisPointInside(x) == "onBoundary") and
               (cell->Position() == GEO::CUT::Point::inside))
      {
#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t\t\t -> Position of point w.r.t volumecell is "
                 << cell->IsThisPointInside(x) << " \t cell pos = " << cell->Position() << IO::endl;
#endif
        // do not return before all the other vcs have been tested, maybe a fluid-cell with
        // onBoundary can be found
        inside_structure = true;
      }
      else
      {
        // the point does not lie inside this vc !
        //#ifdef DEBUG_SEMILAGRANGE
        //        IO::cout << "\n\t\t\t -> Position of cell " << cell->Position() << " and
        //        IsThisPointInside "<< cell->IsThisPointInside(x) << IO::endl;
        //#endif
      }
    }

    if (!inside_structure and (int)(nds.size()) != ele->NumNode())
    {
      GEO::CUT::plain_volumecell_set cells;
      e->GetVolumeCells(cells);

      IO::cout << "point: " << x << IO::endl;
      for (GEO::CUT::plain_volumecell_set::iterator cell_it = cells.begin(); cell_it != cells.end();
           cell_it++)
      {
        GEO::CUT::VolumeCell* cell = *cell_it;
        IO::cout << "vc-pos: " << cell->Position() << IO::endl;
      }

      dserror(
          "no valid nds-vector could be determined and point does not lie within the structure!");
    }

    // return if the structural volume cell is the only one which was found
    if (inside_structure)
    {
      nds.clear();
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n\t\t\t -> Position of point inside structure and not onBoundary of other "
                  "fluid-vcs -> reset nds to empty vector"
               << IO::endl;
#endif

      return;
    }

    IO::cout << "error: coordinates of point x " << x << " number of volumecells: " << cells.size()
             << IO::endl;
    dserror(
        "there is no volume cell in element %d which contains point with coordinates (%f,%f,%f) -> "
        "void element???",
        ele->Id(), x(0), x(1), x(2));
  }
  else  // standard element, all its nodes have dofset 0
  {
    int numnode = ele->NumNode();

    for (int inode = 0; inode < numnode; inode++)
    {
      nds.push_back(0);
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 * compute gradients at nodes for that SL-reconstruction is called                   schott 04/13 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::computeNodalGradient(
    const std::vector<Teuchos::RCP<Epetra_Vector>>&
        colVectors,                    ///< all vectors for that we reconstruct the their gradients
    DRT::Node* node,                   ///< node at which we reconstruct the gradients
    std::vector<DRT::Element*>& eles,  ///< elements around node used for the reconstruction
    std::vector<std::vector<int>>& ele_nds,  ///< corresonding elements nodal dofset information
    XFEM::XFEMDofSet& dofset,                ///< XFEM dofset
    std::vector<LINALG::Matrix<3, 3>>&
        velDeriv_avg,  ///< velocity/acc component derivatives for several vectors
    std::vector<LINALG::Matrix<1, 3>>&
        preDeriv_avg  ///< pressure-component derivatives for several vectors
    ) const
{
  const int nsd = 3;  // dimension

  for (size_t vec_index = 0; vec_index < colVectors.size(); vec_index++)
  {
    velDeriv_avg[vec_index].Clear();
    preDeriv_avg[vec_index].Clear();
  }

  int numele = (int)eles.size();

  if (numele != (int)ele_nds.size()) dserror("number of nds-vector != number of elements!");


  for (int iele = 0; iele < numele; iele++)
  {
    DRT::Element* e = eles[iele];

    // xi coordinates of node w.r.t this element
    LINALG::Matrix<nsd, 1> tmp_xi(true);
    bool indomain = false;  // dummy variable

    // xyz coordinates of the node
    LINALG::Matrix<nsd, 1> x_node(node->X());

    // get the local coordinates of the node w.r.t current element
    // Comment to the configuration:
    // Here we use the fact that x_node is a node position in reference configuration,
    // which means we get the same tmp_xi if we use the element in reference config, as
    // if callXToXi would be done in (n+1) config and x_node^(n+1)
    callXToXiCoords(e, x_node, tmp_xi, "reference", indomain);

    for (size_t tmp_index = 0; tmp_index < colVectors.size(); tmp_index++)
    {
      LINALG::Matrix<nsd, 1> vel(true);
      LINALG::Matrix<nsd, nsd> vel_deriv(true);
      double pres = 0.0;
      LINALG::Matrix<1, nsd> pres_deriv(true);

      getGPValues(e, tmp_xi, ele_nds[iele], dofset, vel, vel_deriv, pres, pres_deriv,
          colVectors[tmp_index], true);

      // add the current gradients
      velDeriv_avg[tmp_index].Update(1.0, vel_deriv, 1.0);
      preDeriv_avg[tmp_index].Update(1.0, pres_deriv, 1.0);
    }
  }  // end loop ele

  for (size_t vec_index = 0; vec_index < colVectors.size(); vec_index++)
  {
    velDeriv_avg[vec_index].Scale(1. / numele);
    preDeriv_avg[vec_index].Scale(1. / numele);
  }

}  // end function compute nodal gradient



/*------------------------------------------------------------------------------------------------*
 * get the time integration factor theta fitting to the computation type             schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
double XFEM::XFLUID_SemiLagrange::Theta(TimeIntData* data) const
{
  double theta = -1.0;

  switch (data->type_)
  {
    case TimeIntData::predictor_:
      theta = 0.0;
      break;
    case TimeIntData::standard_:
      theta = theta_default_;
      break;
    default:
      dserror("type not implemented");
      break;
  }

  if (theta < 0.0) dserror("something wrong");

  return theta;
}  // end function theta



/*------------------------------------------------------------------------------------------------*
 * export alternative algo data to neighbour proc                                    schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::exportAlternativAlgoData()
{
  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  std::vector<std::vector<TimeIntData>> dataVec(numproc_);

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->state_ == TimeIntData::failedSL_)
    {
      dataVec[data->initial_ele_owner_].push_back(*data);
    }
  }

  clearState(TimeIntData::failedSL_);
  timeIntData_->insert(timeIntData_->end(), dataVec[myrank_].begin(), dataVec[myrank_].end());

  dataVec[myrank_].clear();  // clear the set data from the vector

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest
  // higher neighbour...)
  for (int dest = (myrank_ + 1) % numproc_; dest != myrank_;
       dest = (dest + 1) % numproc_)  // dest is the target processor
  {
    // Initialization of sending
    DRT::PackBuffer dataSend;  // vector including all data that has to be send to dest proc

    // Initialization
    int source = myrank_ - (dest - myrank_);  // source proc (sends (dest-myrank_) far and gets from
                                              // (dest-myrank_) earlier)
    if (source < 0)
      source += numproc_;
    else if (source >= numproc_)
      source -= numproc_;

    // pack data to be sent
    for (std::vector<TimeIntData>::iterator data = dataVec[dest].begin();
         data != dataVec[dest].end(); data++)
    {
      if (data->state_ == TimeIntData::failedSL_)
      {
        packNode(dataSend, data->node_);
        DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend, data->vel_);
        DRT::ParObject::AddtoPack(dataSend, data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->dispnp_);
        DRT::ParObject::AddtoPack(dataSend, data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend, (int)data->type_);
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data = dataVec[dest].begin();
         data != dataVec[dest].end(); data++)
    {
      if (data->state_ == TimeIntData::failedSL_)
      {
        packNode(dataSend, data->node_);
        DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend, data->vel_);
        DRT::ParObject::AddtoPack(dataSend, data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->dispnp_);
        DRT::ParObject::AddtoPack(dataSend, data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend, (int)data->type_);
      }
    }

    // clear the no more needed data
    dataVec[dest].clear();

    std::vector<char> dataRecv;
    sendData(dataSend, dest, source, dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0, (double*)coords, 0);
      int nds_np;
      LINALG::Matrix<nsd, 1> vel;
      std::vector<LINALG::Matrix<nsd, nsd>> velDeriv;
      std::vector<LINALG::Matrix<1, nsd>> presDeriv;
      LINALG::Matrix<nsd, 1> dispnp;
      LINALG::Matrix<nsd, 1> initialpoint;
      int initial_eid;
      int initial_ele_owner;
      int newtype;

      unpackNode(posinData, dataRecv, node);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, nds_np);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, vel);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, velDeriv);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, presDeriv);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, dispnp);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, initialpoint);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, initial_eid);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, initial_ele_owner);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, newtype);

      timeIntData_->push_back(TimeIntData(node, nds_np, vel, velDeriv, presDeriv, dispnp,
          initialpoint, initial_eid, initial_ele_owner,
          (TimeIntData::type)newtype));  // startOwner is current proc
    }                                    // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  }  // end loop over processors
}  // end exportAlternativAlgoData



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                                   schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::exportIterData(bool& procDone)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t=============================";
  IO::cout << "\n\t  export Iteration Data  ";
  IO::cout << "\n\t=============================" << IO::endl;
#endif

  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // Initialization
  int dest = myrank_ + 1;  // destination proc (the "next" one)
  if (myrank_ == (numproc_ - 1)) dest = 0;

  int source = myrank_ - 1;  // source proc (the "last" one)
  if (myrank_ == 0) source = numproc_ - 1;


  /*-------------------------------------------*
   * first part: send procfinished in order to *
   * check whether all procs have finished     *
   *-------------------------------------------*/
  for (int iproc = 0; iproc < numproc_ - 1; iproc++)
  {
    DRT::PackBuffer dataSend;

    DRT::ParObject::AddtoPack(dataSend, static_cast<int>(procDone));
    dataSend.StartPacking();
    DRT::ParObject::AddtoPack(dataSend, static_cast<int>(procDone));

    std::vector<char> dataRecv;
    sendData(dataSend, dest, source, dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    size_t posinData = 0;
    int allProcsDone;

    // unpack received data
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, allProcsDone);

    if (allProcsDone == 0) procDone = 0;

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
    for (std::vector<TimeIntData>::iterator data = timeIntData_->begin();
         data != timeIntData_->end(); data++)
    {
      if (data->state_ == TimeIntData::nextSL_)
      {
        packNode(dataSend, data->node_);
        DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend, data->vel_);
        DRT::ParObject::AddtoPack(dataSend, data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->dispnp_);
        DRT::ParObject::AddtoPack(dataSend, data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend, data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend, data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend, data->counter_);
        DRT::ParObject::AddtoPack(dataSend, (int)data->type_);
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data = timeIntData_->begin();
         data != timeIntData_->end(); data++)
    {
      if (data->state_ == TimeIntData::nextSL_)
      {
        packNode(dataSend, data->node_);
        DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend, data->vel_);
        DRT::ParObject::AddtoPack(dataSend, data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend, data->dispnp_);
        DRT::ParObject::AddtoPack(dataSend, data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend, data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend, data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend, data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend, data->counter_);
        DRT::ParObject::AddtoPack(dataSend, (int)data->type_);
      }
    }

    clearState(TimeIntData::nextSL_);

    std::vector<char> dataRecv;
    sendData(dataSend, dest, source, dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0, (double*)coords, 0);
      int nds_np;
      LINALG::Matrix<nsd, 1> vel;
      std::vector<LINALG::Matrix<nsd, nsd>> velDeriv;
      std::vector<LINALG::Matrix<1, nsd>> presDeriv;
      LINALG::Matrix<nsd, 1> dispnp;
      LINALG::Matrix<nsd, 1> initialpoint;
      int initial_eid;
      int initial_ele_owner;
      LINALG::Matrix<nsd, 1> startpoint;
      int searchedProcs;
      int iter;
      int newtype;

      unpackNode(posinData, dataRecv, node);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, nds_np);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, vel);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, velDeriv);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, presDeriv);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, dispnp);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, initialpoint);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, initial_eid);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, initial_ele_owner);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, startpoint);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, searchedProcs);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, iter);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, newtype);

      timeIntData_->push_back(
          TimeIntData(node, nds_np, vel, velDeriv, presDeriv, dispnp, initialpoint, initial_eid,
              initial_ele_owner, startpoint, searchedProcs, iter, (TimeIntData::type)newtype));
    }  // end loop over number of points to get

    // processors wait for each other
    discret_->Comm().Barrier();
  }  // end if procfinished == false
}  // end exportIterData
