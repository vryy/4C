/*!------------------------------------------------------------------------------------------------*
\file startvalues.cpp

\brief provides the Startvalues class

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "startvalues.H"
#include "dof_management_element.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_combust/combust_defines.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"
#include "../linalg/linalg_utils.H"
#include <iostream>



/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrange Back-Tracking algorithm constructor                             winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
XFEM::Startvalues::Startvalues(
    const RCP<DRT::Discretization> discret,
    const RCP<DofManager> olddofman,
    const RCP<DofManager> dofman,
    vector<RCP<Epetra_Vector> > oldVectors,
    vector<RCP<Epetra_Vector> > newVectors,
    const RCP<Epetra_Vector> veln,
    const RCP<COMBUST::FlameFront> flamefront,
    const RCP<InterfaceHandle> interfacehandle_old,
    const Epetra_Map& olddofcolmap,
    const map<DofKey<onNode>, DofGID>& oldNodalDofColDistrib,
    const Epetra_Map& newdofrowmap,
    const map<DofKey<onNode>, DofGID>& newNodalDofRowDistrib,
    const double& dt,
    const double& theta,
    const double& flamespeed
) :
  veln_(veln),
  flamefront_(flamefront),
  oldVectors_(oldVectors),
  newVectors_(newVectors),
  discret_(discret),
  olddofman_(olddofman),
  dofman_(dofman),
  olddofcolmap_(olddofcolmap),
  oldNodalDofColDistrib_(oldNodalDofColDistrib),
  newdofrowmap_(newdofrowmap),
  newNodalDofRowDistrib_(newNodalDofRowDistrib),
  exporter_(discret_->Comm()),
  myrank_(discret_->Comm().MyPID()),
  numproc_(discret_->Comm().NumProc()),
  max_iter_(20),
  dt_(dt),
  theta_default_(theta),
  flamespeed_(flamespeed)
{
  // remark: in flamefront the phi vectors always have to fit to the current
  // discretization so no mapping from node col id to phi dof col id is needed
  phinpi_ = flamefront->Phin(); // last phi vector in column map
  phinpip_ = flamefront->Phinp(); // current phi vector in column map

  const int nsd = 3;

  basic_ = rcp(new vector<StartpointData>);
  curr_ = rcp(new vector<StartpointData>);
  next_ = rcp(new vector<StartpointData>);
  done_ = rcp(new vector<StartpointData>);
  failed_ = rcp(new vector<StartpointData>);



/*------------------------*
 * Initialization         *
 *------------------------*/
  for (int leleid=0; leleid<discret_->NumMyColElements(); leleid++)  // loop over processor nodes
  {
    DRT::Element* iele = discret_->lColElement(leleid);
    if (interfacehandle_old->ElementBisected(iele) || interfacehandle_old->ElementTouchedPlus(iele) || interfacehandle_old->ElementTouchedMinus(iele))
    {
      const int* nodeGids = iele->NodeIds(); // node gids
      for (int inode=0;inode<iele->NumNode();inode++)
      {
        if(olddofcolmap_.MyGID(nodeGids[inode]));
          oldEnrNodes_.insert(nodeGids[inode]);
      }
    }
  }

  // fill curr_ structure with the data for the nodes which changed interface side
  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
  {
    DRT::Node* currnode = discret_->lColNode(lnodeid);

    // node on current processor which changed interface side
    if ((currnode->Owner() == myrank_) &&
        (interfaceSideCompareCombust((*phinpip_)[lnodeid],(*phinpi_)[lnodeid]) == false))
    {
      basic_->push_back(StartpointData(
        *currnode,
        LINALG::Matrix<nsd,1>(true),
        vector<LINALG::Matrix<nsd,nsd> >(newVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),
        vector<LINALG::Matrix<1,nsd> >(newVectors_.size(),LINALG::Matrix<1,nsd>(true)),
        LINALG::Matrix<nsd,1>(true),
        0,
        (*phinpip_)[lnodeid],
        1,
        0,
        vector<int>(1,-1),
        vector<int>(1,-1),
        INFINITY));
    }
  } // end loop over processor nodes

  if (myrank_==0)
    cout << "---  computing new reference solution(s)" << endl;

  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrangian Back-Tracking main algorithm                                  winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::semiLagrangeBackTracking(
    vector<RCP<Epetra_Vector> > newRowVectors,
    bool predictor_step)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  if (predictor_step)
  {
    predictor_step_ = true;
    theta_curr_ = 0.0; // explicit scheme if predictor step

    startpoints(); // compute startpoint predictor for semi-Lagrangian origin

    *curr_ = *basic_;
    basic_->clear();
  }
  else
  {
    predictor_step_ = false;
    theta_curr_ = theta_default_; // standard theta if no predictor step

    newIteration_prepare(newRowVectors);

    *curr_=*done_;
    done_->clear();
  }



/*----------------------------------------------------*
 * first part: get the correct origin for the node    *
 * in a lagrangian point of view using a newton loop  *
 *----------------------------------------------------*/
#ifdef PARALLEL
  int counter = 0; // loop counter to avoid infinite loops
  bool procfinished = false; // true if movNodes is empty

  // loop over nodes which still don't have and may get a good startvalue
  while (true)
  {
    counter += 1;

    // counter limit because maximal max_iter newton iterations with maximal
    // numproc processor changes per iteration (avoids infinite loop)
    if (curr_->size()>0 && counter<max_iter_*numproc_)
    {
      procfinished = false; // processor has still nodes to examine
#endif
      // loop over all nodes which changed interface side
      // remark: negative loop so that deleted elements don't influence other elements position
      for (size_t nodeid=0;nodeid<curr_->size();nodeid++)
      {
        StartpointData& curr = (*curr_)[nodeid]; // current analysed node

//        cout << "on proc " << myrank_ << " iteration starts for\n" <<
//            curr.movNode_ << " with changed " << curr.changed_ << endl;

        // Initialization
        DRT::Element* ele = NULL; // pointer to the element where start point lies in
        LINALG::Matrix<nsd,1> xi(true); // local transformed coordinates of x
        LINALG::Matrix<nsd,1> vel(true); // velocity of the start point approximation
        double phin = 0.0; // phi-value of the start point approximation
        bool elefound = false;             // true if an element for a point was found on the processor
        bool ihSide = true;             // true if point lies on correct side
        bool stdBackTracking = true; // true if standard back-tracking shall be done
        double deltaT1; // at t^n + deltaT1 the lagrangian point crosses the interface

        // search for an element where the current startpoint lies in
        // if found, give out all data at the startpoint
        callElementSearch(ele,curr.startpoint_,xi,vel,phin,elefound); // template data is dummy data here

        // if element is not found, look at another processor and so add all
        // according data to the vectors which will be sent to the next processor
        if (!elefound)
        {
          if (curr.searchedProcs_ < numproc_)
          {
            StartpointData next(curr);
            next.searchedProcs_ += 1;
            next_->push_back(next);
          }
          else // all procs searched -> point not in domain
          {
//            cout << "velderiv size is " << curr.velDeriv_.size() << endl;
            failed_->push_back(StartpointData(curr));
            cout << "WARNING! Lagrangian start point not in domain!" << endl;
          }
        } // end if elefound
        else  // if element is found, the newton iteration to find a better startpoint can start
        {
          if (interfaceSideCompareCombust(curr.phiValue_,phin) == false) // Lagrangian origin and original node on different interface sides
          {
            ihSide = false;
          }
          else  // Newton loop just for sensible points
          {
            ihSide = true;
            NewtonLoop(ele,curr,xi,vel,phin,elefound,ihSide,stdBackTracking,deltaT1);
          }

          // if iteration can go on (that is when startpoint is on
          // correct interface side and iter < max_iter)
          if ((curr.iter_<max_iter_) and (ihSide==true))
          {
            // if element is not found in a newton step, look at another processor and so add
            // all according data to the vectors which will be sent to the next processor
            if (!elefound)
            {
              curr.searchedProcs_ = 2;
              next_->push_back(curr);
            }
            else // newton iteration converged to a good startpoint and so the data can be used to go on
            {
              const char* backTrackingType;
              if (stdBackTracking) backTrackingType = static_cast<const char*>("standard");
              else if (!stdBackTracking) backTrackingType = static_cast<const char*>("alternative");
              callBackTracking(
                  ele,curr,xi,backTrackingType,deltaT1);
            }
          } // end if ihSide

          // if startpoint lies on wrong side at any time, this algorithm has to stop
          //for the searched node which is done by setting changed = false
          if (ihSide==false)
            curr.changed_ = false;

          if (curr.iter_ == max_iter_) // maximum number of iterations reached
          {
            cout << "WARNING: start value set somehow sensible\n"
                    "  but newton iteration to find it didnt converge!" << endl;
            const char* backTrackingType;
            if (stdBackTracking) backTrackingType = static_cast<const char*>("standard");
            else if (!stdBackTracking) backTrackingType = static_cast<const char*>("alternative");
            callBackTracking(
                ele,curr,xi,backTrackingType,deltaT1);
          }
        } // end if over nodes which changed interface

//        cout << "on proc " << myrank_ << " after " << curr.iter_ <<
//            " iterations the startpoint is " << curr.startpoint_ << endl;
      } // end loop over all nodes with changed interface side
    } // end if movenodes is empty or max. iter reached
    else // maximum iteration number reached or no startpoints left to look for on this proc
      procfinished = true;

    // all nodes which got changed = false in this step need another algorithm to be set
    for (size_t nodeid=0;nodeid<curr_->size();nodeid++)
    {
      StartpointData& curr = (*curr_)[nodeid];

      if (!(*curr_)[nodeid].changed_)
      {
        failed_->push_back(StartpointData(curr));
      }
    }

#ifdef PARALLEL
    // export nodes and according data for which the startpoint isn't still found (next_ vector) to next proc
    exportIterData(procfinished);

    // convergencecheck: procfinished == 1 just if all procs have finished
    if (procfinished)
      break;

    *curr_ = *next_;
    next_->clear();
  } // end while loop over searched nodes
#endif



/*-----------------------------------------------------------------------------*
 * second part: get sensible startvalues for nodes where the algorithm failed, *
 * using another algorithm, and combine the "Done" and the "Failed" - vectors  *
 *-----------------------------------------------------------------------------*/

  // nodes which are still in curr-vector, ran into an infinite loop
  // and so have to be set otherwise (this should not happen!)
  if (!procfinished)
  {
    for (size_t inode=0;inode<curr_->size();inode++)
    {
      failed_->push_back(StartpointData(
          (*curr_)[inode].movNode_,
          (*curr_)[inode].vel_,
          (*curr_)[inode].velDeriv_,
          (*curr_)[inode].presDeriv_,
          (*curr_)[inode].startpoint_,
          (*curr_)[inode].phiValue_,
          (*curr_)[inode].startGid_,
          (*curr_)[inode].startOwner_));
      cout << "WARNING! Node ran into infinite loop in startvalues.cpp!" << endl;
    }
  }

  (*curr_).clear(); // no more needed

#ifdef PARALLEL
  exportAlternativAlgoData(); // export data of failed nodes
#endif

  getDataForNotConvergedNodes(); // compute final data for failed nodes
  (*failed_).clear(); // no more needed



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
  if (counter > 8*numproc_) // too much loops shouldnt be if all this works
    cout << "WARNING: semiLagrangeExtrapolation seems to run an infinite loop!" << endl;
#endif
} // end semiLagrangeExtrapolation



/*------------------------------------------------------------------------------------------------*
 * Compute startvalues for the interface-changing nodes                          winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::startpoints()
{
  //Initialization
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  const double TOL = 1.0e-3;

  // loop over processors
  for (int procid=0; procid<numproc_; procid++)
  {
    // loop over nodes which changed interface side
    for (size_t nodeid=0;nodeid<basic_->size();nodeid++)
    {
      StartpointData& basic = (*basic_)[nodeid]; // current analysed node (modified in basic_ vector also)
      bool currChanged = false; // startvalue changed on current proc?
      LINALG::Matrix<nsd,1> newNodeCoords(basic.movNode_.X()); // coords of endpoint of Lagrangian characteristics

      // loop over intersected elements on processor
      for (std::set<int>::const_iterator enrnode = oldEnrNodes_.begin();
          enrnode != oldEnrNodes_.end();
          enrnode++)
      {
        DRT::Node* lnodeold = discret_->gNode(*enrnode);//elenodeids[elenode]);  // node near interface

        LINALG::Matrix<nsd,1> oldNodeCoords(lnodeold->X());  // coords of potential startpoint

        // just look for points on the same interface side and on the same processor in row map
        if ((interfaceSideCompareCombust(basic.phiValue_,(*phinpi_)[lnodeold->LID()])) and
            (lnodeold->Owner()==myrank_))
        {
          LINALG::Matrix<nsd,1> diff;  // vector from old point at time n to new point at time n+1
          diff.Update(1.0,newNodeCoords,-1.0,oldNodeCoords);

          if (diff.Norm2() < 2*basic.dMin_) // possible new nearest node
          {
            // get nodal velocity of a node
            LINALG::Matrix<nsd,1> nodevel(true);  // velocity of "old" node at time n
            const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(*enrnode));
            for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
                fieldenr != fieldenrset.end();++fieldenr)
            {
              const DofKey<onNode> olddofkey(*enrnode,*fieldenr);
              const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;

              if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
              {
                if (fieldenr->getField() == XFEM::PHYSICS::Velx)
                  nodevel(0) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                if (fieldenr->getField() == XFEM::PHYSICS::Vely)
                  nodevel(1) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                if (fieldenr->getField() == XFEM::PHYSICS::Velz)
                  nodevel(2) = (*veln_)[olddofcolmap_.LID(olddofpos)];
              }
            } // end loop over fieldenr

            LINALG::Matrix<1,1> arc(true); // cosinus of angle between dist and vel(x_n+1)*diff.Norm2()
            arc.MultiplyTN(1.0/nodevel.Norm2(),diff,nodevel);
            double dist = diff.Norm2() + (diff.Norm2()-arc.Norm2());

            if (dist-basic.dMin_+TOL*dist < 0) // new nearest node (cosinus shall be near 1!)
            {
              basic.startpoint_.Update(1.0,newNodeCoords,-dt_,nodevel);
              basic.dMin_ = dist;
              basic.startGid_.clear();
              basic.startGid_.push_back(*enrnode);
              basic.startOwner_.clear();
              basic.startOwner_.push_back(myrank_);
              currChanged = true;
            } // end if new best startvalue
            else if ((dist-basic.dMin_ > -TOL*dist) and
                (dist-basic.dMin_ < TOL*dist)) // handles special case that two nodes are very similar near
            {
              basic.startpoint_.Update(0.5,newNodeCoords,-0.5*dt_,nodevel,0.5); // midpoint is 0.5*(x-dt*vel+old_startpoint)
              basic.dMin_ = (dist+basic.dMin_)/2.0;
              basic.startGid_.push_back(*enrnode);
              basic.startOwner_.push_back(myrank_);
            }
          } // end if possibly new best startvalue
        } // end if just points on the correct side
      } // end loop over intersected elements
      if (currChanged == true) basic.changed_ = true;
    } // end loop over nodes which changed interface side

#ifdef PARALLEL
    exportStartData();
#endif

  } // end loop over processors

  // test loop over all nodes which changed interface side
  for (size_t nodeid=0; nodeid<basic_->size(); nodeid++)
  {
    if (!(*basic_)[nodeid].changed_)
    {
      dserror("WARNING! No enriched node on one interface side found!\nThis "
          "indicates that the whole area is at one side of the interface!");
      break;
    }
  } // end loop over nodes
} // end startValuesFinder



/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::NewtonLoop(
    DRT::Element*& ele,
    StartpointData& curr,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& elefound,
    bool& ihSide,
    bool& stdBackTracking,
    double& deltaT1
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  stdBackTracking = true; // standard Newton loop -> standard back tracking

  // Initialization
  LINALG::Matrix<nsd,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<nsd,1> incr(true);                 // increment of the newton system

  const double relTolIncr = 1.0e-10;;  // tolerance for the increment
  const double relTolRes = 1.0e-10;    // tolerance for the residual

  LINALG::Matrix<nsd,1> origNodeCoords(true); // coordinates of endpoint of Lagrangian characteristics
  for (int i=0;i<nsd;i++)
    origNodeCoords(i) = curr.movNode_.X()[i];

  // initialize residual (Theta = 0 at predictor step)
  residuum.Clear();
  residuum.Update((1.0-theta_curr_),vel,theta_curr_,curr.vel_); // dt*v(curr.startpoint_)
  residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords,dt_);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)

  while(curr.iter_ < max_iter_)  // newton loop
  {
    curr.iter_ += 1;

    switch (ele->Shape())
    {
    case DRT::Element::hex8:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      NewtonIter<numnode,DRT::Element::hex8>(ele,curr,xi,vel,residuum,incr,phi,elefound);
    }
    break;
    case DRT::Element::hex20:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
      NewtonIter<numnode,DRT::Element::hex20>(ele,curr,xi,vel,residuum,incr,phi,elefound);
    }
    break;
    default:
      dserror("xfem assembly type not yet implemented in time integration");
    };

    if (elefound) // element of curr.startpoint_ at this processor
    {
      if (interfaceSideCompareCombust(curr.phiValue_,phi) == false)
      {
        ihSide = false;
        break; // leave newton loop if element is on wrong domain side
      }
      else
      {
        ihSide = true;

        // reset residual
        residuum.Clear();
        residuum.Update((1.0-theta_curr_),vel,theta_curr_,curr.vel_); // dt*v(curr.startpoint_)
        residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords,dt_);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)

        // convergence criterion
        if (curr.startpoint_.Norm2()>1e-3)
        {
          if (incr.Norm2()/curr.startpoint_.Norm2() < relTolIncr && residuum.Norm2()/curr.startpoint_.Norm2() < relTolRes)
            break;
        }
        else
        {
          if (incr.Norm2() < relTolIncr && residuum.Norm2() < relTolRes)
            break;
        }
      }
    }
    else // element of curr.startpoint_ not at this processor
      break; // stop newton loop on this proc
  }
#ifdef DEBUG
  // did newton iteration converge?
  if(curr.iter_ == max_iter_){cout << "WARNING: newton iteration for finding start value not converged for point\n" << endl;}
//    cout << "after " << curr.iter_ << " iterations the endpoint is\n" << xAppr << endl;
#endif
} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::Startvalues::NewtonIter(
    DRT::Element*& ele,
    StartpointData& curr,
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

#ifdef COMBUST_NORMAL_ENRICHMENT
  LINALG::Matrix<1,numnode> nodevelenr(true);
  LINALG::Matrix<nsd,nsd> vderxy(true);
  ApproxFuncNormalVector<2,2*numnode> shp;
#else
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true); // dummy
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
#endif

  // build systemMatrix
#ifdef COMBUST_NORMAL_ENRICHMENT
  pointdataXFEMNormal<numnode,DISTYPE>(
      ele,
#ifdef COLLAPSE_FLAME_NORMAL
      curr.startpoint_,
#endif
      xi,
      xji,
      shapeFcn,
      enrShapeFcnPres,
      enrShapeXYPresDeriv1,
      shp,
      false
  );

  elementsNodalData<numnode>(
      ele,
      veln_,
      olddofman_,
      olddofcolmap_,
      oldNodalDofColDistrib_,
      nodevel,
      nodevelenr,
      nodepres); // nodal data of the element nodes

  vderxy.Clear();
  { // build sysmat
    const int* nodeids = ele->NodeIds();
    size_t velncounter = 0;

    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    for (size_t inode = 0; inode < numnode; ++inode) // loop over element nodes
    {
      // standard shape functions are identical for all vector components
      // e.g. shp.velx.dx.s == shp.vely.dx.s == shp.velz.dx.s
      vderxy(0,0) += nodevel(0,inode)*shp.velx.dx.s(inode);
      vderxy(0,1) += nodevel(0,inode)*shp.velx.dy.s(inode);
      vderxy(0,2) += nodevel(0,inode)*shp.velx.dz.s(inode);

      vderxy(1,0) += nodevel(1,inode)*shp.vely.dx.s(inode);
      vderxy(1,1) += nodevel(1,inode)*shp.vely.dy.s(inode);
      vderxy(1,2) += nodevel(1,inode)*shp.vely.dz.s(inode);

      vderxy(2,0) += nodevel(2,inode)*shp.velz.dx.s(inode);
      vderxy(2,1) += nodevel(2,inode)*shp.velz.dy.s(inode);
      vderxy(2,2) += nodevel(2,inode)*shp.velz.dz.s(inode);

      const int gid = nodeids[inode];
      const std::set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(gid);

      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        if (enrfield->getField() == XFEM::PHYSICS::Veln)
        {
          vderxy(0,0) += nodevelenr(0,velncounter)*shp.velx.dx.n(velncounter);
          vderxy(0,1) += nodevelenr(0,velncounter)*shp.velx.dy.n(velncounter);
          vderxy(0,2) += nodevelenr(0,velncounter)*shp.velx.dz.n(velncounter);

          vderxy(1,0) += nodevelenr(0,velncounter)*shp.vely.dx.n(velncounter);
          vderxy(1,1) += nodevelenr(0,velncounter)*shp.vely.dy.n(velncounter);
          vderxy(1,2) += nodevelenr(0,velncounter)*shp.vely.dz.n(velncounter);

          vderxy(2,0) += nodevelenr(0,velncounter)*shp.velz.dx.n(velncounter);
          vderxy(2,1) += nodevelenr(0,velncounter)*shp.velz.dy.n(velncounter);
          vderxy(2,2) += nodevelenr(0,velncounter)*shp.velz.dz.n(velncounter);

          velncounter += 1;
        }
      }
    } // end loop over element nodes
  } // sysmat built

  sysmat.Update((1.0-theta_curr_)*dt_,vderxy);
#else
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

  sysmat.MultiplyNT((1.0-theta_curr_)*dt_,nodevel,enrShapeXYVelDeriv1); // (1-theta) * dt * v_nodes * dN/dx
#endif
  for (int i=0;i<nsd;i++)
    sysmat(i,i) += 1.0; // I + dt*velDerivXY
  sysmat.Invert();
  // invers system Matrix built

  //solve Newton iteration
  incr.Clear();
  incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

  // update iteration
  for (int i=0;i<nsd;i++)
    curr.startpoint_(i) += incr(i);
    //cout << "in newton loop: approximate startvalue is " << curr.startpoint_(0) << " " << curr.startpoint_(1) << " " << curr.startpoint_(2) << endl;

  //=============== update residuum================
  callElementSearch(ele, curr.startpoint_, xi, vel, phi, elefound);
} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * Computing final data where semi-Lagrangian approach failed                    winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::getDataForNotConvergedNodes(
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // remark: all data has to be sent to the processor where
  //         the startpoint lies before calling this function
  for (size_t inode=0;inode<failed_->size();inode++)
  {
    StartpointData failed = (*failed_)[inode]; // data of node where sl approach failed

    if (failed.startGid_.size() != 1)
      dserror("data for alternative nodes shall be computed for one node here");

    DRT::Node* node = discret_->gNode(failed.startGid_[0]); // failed node
    DRT::Element* nodesElement = node->Elements()[0]; // element of failed node
    failed.startpoint_=LINALG::Matrix<nsd,1>(node->X());

    LINALG::Matrix<nsd,1> x(node->X()); // coordinates of failed node
    LINALG::Matrix<nsd,1> xi(true); // local coordinates of failed node
    LINALG::Matrix<nsd,1> vel(true); // velocity at pseudo-Lagrangian origin
    double phi = 0;
    bool elefound = false;

/*------------------------------------------*
 * element data at pseudo-Lagrangian origin *
 * remark: an element must be found since   *
 *         the intersected node must be on  *
 *         the current processor            *
 *------------------------------------------*/
    callElementSearch(nodesElement,x,xi,vel,phi,elefound);

    if (!elefound)
      dserror("element of a row node of a proc should be in procs colmap!");
    else
    {
/*---------------------------------*
 * correction step:                *
 * vel = dt*velDeriv1*vel + vel    *
 * pres = pres + dt*presDeriv1*vel *
 * with pseudo time-step dt        *
 *---------------------------------*/
      callBackTracking(
          nodesElement,failed,xi,static_cast<const char*>("failing"),0.0);
    } // end if elefound
  } // end loop over nodes
} // end getDataForNotConvergedNodes



/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::callBackTracking(
    DRT::Element*& ele,
    StartpointData& data,
    LINALG::Matrix<3,1>& xi,
    const char* backTrackingType,
    double deltaT1
)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    backTracking<numnode,DRT::Element::hex8>(ele,data,xi,backTrackingType,deltaT1);
  }
  break;
  case DRT::Element::hex20:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    backTracking<numnode,DRT::Element::hex20>(ele,data,xi,backTrackingType,deltaT1);
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
void XFEM::Startvalues::backTracking(
    DRT::Element*& fittingele,
    StartpointData& data,
    LINALG::Matrix<3,1>& xi,
    const char* backTrackingType,
    double deltaT1
)
{
//  cout << "backtracking for node" << data.movNode_ << endl;
  const int nsd = 3;

  if ((strcmp(backTrackingType,static_cast<const char*>("standard"))!=0) and
      (strcmp(backTrackingType,static_cast<const char*>("alternative"))!=0) and
      (strcmp(backTrackingType,static_cast<const char*>("failing"))!=0))
    dserror("backTrackingType not implemented");

//  cout << data.movNode_ << "has lagrange origin " << data.startpoint_ << "with xi-coordinates "
//      << xi << "in element " << *fittingele << endl;
  LINALG::Matrix<nsd,nsd> xji(true); // invers of jacobian
  LINALG::Matrix<numnode,1> shapeFcn(true); // shape function
  double deltaT = 0; // pseudo time-step size

  // enriched shape functions and there derivatives in local coordinates (N * \Psi)
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);

  // data for the final back-tracking
  LINALG::Matrix<nsd,1> vel(true); // velocity data
  vector<LINALG::Matrix<nsd,1> > veln(oldVectors_.size(),LINALG::Matrix<nsd,1>(true)); // velocity at t^n
  vector<LINALG::Matrix<nsd,nsd> > velnDeriv1(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)); // first derivation of velocity data
  LINALG::Matrix<1,1> pres(true); // pressure data
  vector<LINALG::Matrix<1,nsd> > presnDeriv1(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)); // first derivation of pressure data
  LINALG::Matrix<nsd,1> transportVeln(true); // transport velocity at Lagrangian origin (x_Lagr(t^n))

  int numele;
  DRT::Element** nodeeles = NULL;

  if ((data.startGid_.size() != 1) and
      (strcmp(backTrackingType,static_cast<const char*>("failing")) == 0))
    dserror("back-tracking shall be done only for one node here!");

  DRT::Node* node = discret_->gNode(data.startGid_[0]);
  if (strcmp(backTrackingType,"failing")==0)
  {
    numele=node->NumElement();
    nodeeles = node->Elements();
  }
  else
  {
    numele=1;
  }

  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd,2*numnode> nodevel(true);
  // node velocities of the element nodes for the data that should be changed
  vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,2*numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size(),LINALG::Matrix<1,2*numnode>(true));
#ifdef COMBUST_NORMAL_ENRICHMENT
    vector<LINALG::Matrix<1,numnode> > nodevelenrdata(oldVectors_.size(),LINALG::Matrix<1,numnode>(true));
    LINALG::Matrix<1,numnode> nodevelenr(true);
#endif
    vector<LINALG::Matrix<nsd,1> > velValues(oldVectors_.size(),LINALG::Matrix<nsd,1>(true)); // velocity of the data that should be changed
    vector<double> presValues(oldVectors_.size(),0); // pressures of the data that should be changed

    LINALG::Matrix<1,1> phinpi(true); // phivalue (currently not required)

  for (int iele=0;iele<numele;iele++) // loop over elements containing the startpoint (usually one)
  {
    for (size_t index=0;index<oldVectors_.size();index++)
    {
      nodeveldata[index].Clear();
      nodepresdata[index].Clear();
#ifdef COMBUST_NORMAL_ENRICHMENT
      nodevelenrdata[index].Clear();
#endif
    }

    DRT::Element* ele = NULL;
    if (numele>1)
    {
      ele = nodeeles[iele];
      LINALG::Matrix<nsd,1> coords(node->X());
      LINALG::Matrix<nsd,1> vel(true); // dummy
      bool elefound = false;
      callElementSearch(ele,coords,xi,vel,phinpi(0),elefound);

      if ((!elefound) or
          (ele->Id()!=nodeeles[iele]->Id())) // node of element should lie in element...
      {
        dserror("node of an element is not lying in the element! BUG?!");
      }
    }
    else
    {
      ele = fittingele;
    }

    // evaluate data for the given point
#ifdef COMBUST_NORMAL_ENRICHMENT
#ifdef COLLAPSE_FLAME_NORMAL
    LINALG::Matrix<nsd,1> normal(node->X());
    normal(2) = 0.0;
    normal.Scale(-1.0/normal.Norm2());
#endif
    ApproxFuncNormalVector<2,2*numnode> shp;
    pointdataXFEMNormal<numnode,DISTYPE>(
        ele,
#ifdef COLLAPSE_FLAME_NORMAL
        normal,
#endif
        xi,
        xji,
        shapeFcn,
        enrShapeFcnPres,
        enrShapeXYPresDeriv1,
        shp,
        false
    );
#else
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
#endif

    const int* elenodeids = ele->NodeIds();

    int dofcounterVelx = 0;
    int dofcounterVely = 0;
    int dofcounterVelz = 0;
    int dofcounterPres = 0;
#ifdef COMBUST_NORMAL_ENRICHMENT
    int dofcounterVeln = 0;
#endif

    for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
    {
      // get nodal velocities and pressures with help of the field set of node
      const std::set<XFEM::FieldEnr>& fieldEnrSet(olddofman_->getNodeDofSet(elenodeids[nodeid]));
      for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
          fieldenr != fieldEnrSet.end();++fieldenr)
      {
        const DofKey<onNode> olddofkey(elenodeids[nodeid], *fieldenr);
        const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
        switch (fieldenr->getEnrichment().Type())
        {
        case XFEM::Enrichment::typeStandard :
        case XFEM::Enrichment::typeJump :
        case XFEM::Enrichment::typeVoidFSI :
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
#ifdef COMBUST_NORMAL_ENRICHMENT
          else if (fieldenr->getField() == XFEM::PHYSICS::Veln)
          {
            nodevelenr(0,dofcounterVeln) =  (*veln_)[olddofcolmap_.LID(olddofpos)];
            for (size_t index=0;index<oldVectors_.size();index++)
              nodevelenrdata[index](0,dofcounterVeln) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            dofcounterVeln++;
          }
#endif
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          {
            for (size_t index=0;index<oldVectors_.size();index++)
              nodepresdata[index](0,dofcounterPres) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
            dofcounterPres++;
          }
          else
          {
            cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << endl;
            dserror("not implemented physical field!");
          }
          break;
        }
        case XFEM::Enrichment::typeUndefined :
        default :
        {
          cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << endl;
          dserror("unknown enrichment type");
          break;
        }
        } // end switch enrichment
      } // end loop over fieldenr
    } // end loop over element nodes

    if (iele==0) // compute transportvel just once!
    {
#ifdef COMBUST_NORMAL_ENRICHMENT
      size_t velncounter = 0;
      for (int inode=0; inode<numnode; ++inode)
      {
        // standard shape functions are identical for all vector components
        // shp.velx.d0.s == shp.vely.d0.s == shp.velz.d0.s
        transportVeln(0) += nodevel(0,inode)*shp.velx.d0.s(inode);
        transportVeln(1) += nodevel(1,inode)*shp.vely.d0.s(inode);
        transportVeln(2) += nodevel(2,inode)*shp.velz.d0.s(inode);

        const set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(elenodeids[inode]);
        for (set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            transportVeln(0) += nodevelenr(0,velncounter)*shp.velx.d0.n(velncounter);
            transportVeln(1) += nodevelenr(0,velncounter)*shp.vely.d0.n(velncounter);
            transportVeln(2) += nodevelenr(0,velncounter)*shp.velz.d0.n(velncounter);

            velncounter += 1;
          }
        }
      } // end loop over elenodes
#else
      // interpolate velocity and pressure values at starting point
      transportVeln.Multiply(nodevel, enrShapeFcnVel);
#endif

      // computing pseudo time-step deltaT
      // remark: if x is the Lagrange-origin of node, deltaT = dt with respect to small errors.
      // if its not, deltaT estimates the time x needs to move to node)
      if (predictor_step_)
      {
        LINALG::Matrix<nsd,1> diff(data.movNode_.X());
        diff -= data.startpoint_; // diff = x_Node - x_Appr

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
#ifdef COMBUST_NORMAL_ENRICHMENT
      size_t velncounter = 0;
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (int inode = 0; inode < numnode; ++inode)
      {
        // standard shape functions are identical for all vector components
        // shp.velx.d0.s == shp.vely.d0.s == shp.velz.d0.s
        if (iele==0)
        {
          veln[index](0) += nodevel(0,inode)*shp.velx.d0.s(inode);
          veln[index](1) += nodevel(1,inode)*shp.vely.d0.s(inode);
          veln[index](2) += nodevel(2,inode)*shp.velz.d0.s(inode);
        }

        velnDeriv1[index](0,0) += nodeveldata[index](0,inode)*shp.velx.dx.s(inode);
        velnDeriv1[index](0,1) += nodeveldata[index](0,inode)*shp.velx.dy.s(inode);
        velnDeriv1[index](0,2) += nodeveldata[index](0,inode)*shp.velx.dz.s(inode);

        velnDeriv1[index](1,0) += nodeveldata[index](1,inode)*shp.vely.dx.s(inode);
        velnDeriv1[index](1,1) += nodeveldata[index](1,inode)*shp.vely.dy.s(inode);
        velnDeriv1[index](1,2) += nodeveldata[index](1,inode)*shp.vely.dz.s(inode);

        velnDeriv1[index](2,0) += nodeveldata[index](2,inode)*shp.velz.dx.s(inode);
        velnDeriv1[index](2,1) += nodeveldata[index](2,inode)*shp.velz.dy.s(inode);
        velnDeriv1[index](2,2) += nodeveldata[index](2,inode)*shp.velz.dz.s(inode);

        const std::set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(elenodeids[inode]);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            if (iele==0)
            {
              veln[index](0) += nodevelenr(0,velncounter)*shp.velx.d0.n(velncounter);
              veln[index](1) += nodevelenr(0,velncounter)*shp.vely.d0.n(velncounter);
              veln[index](2) += nodevelenr(0,velncounter)*shp.velz.d0.n(velncounter);
            }


            velnDeriv1[index](0,0) += nodeveldata[index](0,velncounter)*shp.velx.dx.n(velncounter);
            velnDeriv1[index](0,1) += nodeveldata[index](0,velncounter)*shp.velx.dy.n(velncounter);
            velnDeriv1[index](0,2) += nodeveldata[index](0,velncounter)*shp.velx.dz.n(velncounter);

            velnDeriv1[index](1,0) += nodeveldata[index](0,velncounter)*shp.vely.dx.n(velncounter);
            velnDeriv1[index](1,1) += nodeveldata[index](0,velncounter)*shp.vely.dy.n(velncounter);
            velnDeriv1[index](1,2) += nodeveldata[index](0,velncounter)*shp.vely.dz.n(velncounter);

            velnDeriv1[index](2,0) += nodeveldata[index](0,velncounter)*shp.velz.dx.n(velncounter);
            velnDeriv1[index](2,1) += nodeveldata[index](0,velncounter)*shp.velz.dy.n(velncounter);
            velnDeriv1[index](2,2) += nodeveldata[index](0,velncounter)*shp.velz.dz.n(velncounter);

            velncounter += 1;
          } // end if veln field
        } // end loop over fieldenrset
      } // end loop over element nodes
      if (index==0)
        cout << *ele << "with nodevels " << nodeveldata[0] << ", nodeenrvals " << nodevelenrdata[0]
            << " and summed up velnderiv currently is " << velnDeriv1[0] << endl;
#else
      if (iele==0)
        veln[index].Multiply(nodeveldata[index],enrShapeFcnVel);
      velnDeriv1[index].MultiplyNT(1.0,nodeveldata[index],enrShapeXYVelDeriv1,1.0);
//      if (index==0)
//        cout << *ele << "with nodevels " << nodeveldata[index] << ", shapefcnderiv " <<
//            enrShapeXYVelDeriv1 << " and summed up velnderiv currently is " << velnDeriv1[index] << endl;
#endif
      presnDeriv1[index].MultiplyNT(1.0,nodepresdata[index],enrShapeXYPresDeriv1,1.0);
    } // end loop over vectors to be read from
  } // end loop over elements containing the point (usually one)

  for (size_t index=0;index<oldVectors_.size();index++)
  {
    velnDeriv1[index].Scale(1.0/numele);
    presnDeriv1[index].Scale(1.0/numele);

    if ((strcmp(backTrackingType,"standard")==0) or
        (strcmp(backTrackingType,"failing"))==0)
    {
      vel.Multiply(1.0-theta_curr_,velnDeriv1[index],transportVeln); // v = (1-theta)*Dv^n/Dx*v^n
      vel.Multiply(theta_curr_,data.velDeriv_[index],data.vel_,1.0); // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
      vel.Update(1.0,veln[index],deltaT); // v = v_n + dt*(theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n)
      velValues[index]=vel;

      pres.Multiply(1.0-theta_curr_,presnDeriv1[index],transportVeln); // p = (1-theta)*Dp^n/Dx*v^n
      pres.Multiply(theta_curr_,data.presDeriv_[index],data.vel_,1.0); // p = theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n
      pres.Multiply(1.0,nodepresdata[index],enrShapeFcnPres,deltaT); // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
      presValues[index] = pres(0);
    }
    else
    {
      dserror("out of use!");
    }
  } // loop over vectors to be set

  done_->push_back(StartpointData(
      data.movNode_,
      data.startpoint_,
      data.phiValue_,
      data.startGid_,
      vector<int>(1,myrank_),
      velValues,
      presValues));
} // end backTracking



/*------------------------------------------------------------------------------------------------*
 * setting the final data in Epetra Vector for a node                            winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::setFinalData(
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  map<int,int> usedStartpoints;
  int numStartpoints;
  double newValue = 0.0;

  for (size_t inode=0;inode<done_->size();inode++)
  {
    vector<LINALG::Matrix<nsd,1> > velValues((*done_)[inode].velValues_); // velocities of the node
    vector<double> presValues((*done_)[inode].presValues_); // pressures of the node

    const int gnodeid = (*done_)[inode].movNode_.Id(); // global node id

    map<int,int>::iterator currstartpoint = usedStartpoints.find(gnodeid);
    if (currstartpoint==usedStartpoints.end()) // standard case and "standard alternative" case
    {
      usedStartpoints.insert(pair<int,int>(gnodeid,1));
      numStartpoints = 1;
    }
    else
    {
      currstartpoint->second +=1;
      numStartpoints = currstartpoint->second;
    }

    // set nodal velocities and pressures with help of the field set of node
    const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gnodeid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey<onNode> newdofkey(gnodeid, *fieldenr);
      const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;

      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
/*
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[0](0) << endl;
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[0](1) << endl;
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[0](2) << endl;
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
          cout << (*newVectors_[0])[newdofrowmap_.LID(newdofpos)] << " becomes " << presValues[0] << endl;
*/

        for (size_t index=0;index<newVectors_.size();index++)
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
            newValue = ((numStartpoints-1.0)/numStartpoints)*(*newVectors_[index])[newdofrowmap_.LID(newdofpos)]
                       + velValues[index](0)/numStartpoints;
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
            newValue = ((numStartpoints-1.0)/numStartpoints)*(*newVectors_[index])[newdofrowmap_.LID(newdofpos)]
                       + velValues[index](1)/numStartpoints;
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
            newValue = ((numStartpoints-1.0)/numStartpoints)*(*newVectors_[index])[newdofrowmap_.LID(newdofpos)]
                       + velValues[index](2)/numStartpoints;
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
            newValue = ((numStartpoints-1.0)/numStartpoints)*(*newVectors_[index])[newdofrowmap_.LID(newdofpos)]
                       + presValues[index]/numStartpoints;

          (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = newValue;
        }
      }
    } // end loop over fieldenr
  } // end loop over nodes
} // end setFinalData



/*------------------------------------------------------------------------------------------------*
 * evaluate element and data for a given point                                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::callElementSearch(
    DRT::Element*& ele,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& elefound
) const
{
  elefound = false; // becomes true in the algo if element is found
  DRT::Element* currele = NULL;

  int startid; // local row element id
  if (ele==NULL) startid = 0; // start with first local row element
  else           startid = -1; // pseudo-id so that id+1 will be 0

  //loop over elements
  for (int ieleid = startid;ieleid<discret_->NumMyRowElements();ieleid++)
  {
    // if ele != NULL and so initialized,
    // first it should be checked if it is fitting
    if (ieleid == -1)
    {
      currele = ele;
      ele = NULL; // ele will be set if an element is found finally
    }
    else
    {
      currele = discret_->lRowElement(ieleid);
    }

    switch (currele->Shape())
    {
    case DRT::Element::hex8:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      elementAndLocalCoords<numnode,DRT::Element::hex8>(currele,x,xi,vel,phi,elefound);
    }
    break;
    case DRT::Element::hex20:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
      elementAndLocalCoords<numnode,DRT::Element::hex20>(currele,x,xi,vel,phi,elefound);
    }
    break;
    default:
      dserror("xfem assembly type not yet implemented in time integration");
    };

    if (elefound)
    {
      ele = currele;
      break;
    }
  }
} // end function findElementAndLocalCoords



/*------------------------------------------------------------------------------------------------*
 * evaluate element and data for a given point                                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::Startvalues::elementAndLocalCoords(
    DRT::Element*& ele,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& elefound
) const
{
  elefound = false; // becomes true in the algo if element is found
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization of what every following Newton loop will need
  LINALG::Matrix<nsd,1> incr(true);        // increment of newton iteration
  LINALG::Matrix<nsd,1> residuum(true);    //residuum of newton iteration

  LINALG::Matrix<numnode,1> shapeFcn(true);               // shape functions
  LINALG::Matrix<nsd,numnode> shapeXiDeriv1(true);        // derivation of shape functions

  LINALG::Matrix<nsd,nsd> xjm(true);        //jacobian at start vector
  LINALG::Matrix<nsd,nsd> xji(true);        // invers of jacobian

  int iter = 0;   // iteration counter

  const int max_iter = 10;   // maximum number of iterations
  const double relTolRes = 1e-10;   // relative tolerance for residual
  const double relTolIncr = 1e-10;   // relaltive tolerance for increment
  const double XiTol = 1e-3; // max. distance from point to cell in xi-coordinates

  if (ele == NULL) // take first row ele
    ele = discret_->lRowElement(0);
  else
  {
    if (numnode < 3)
      dserror("No 2D or 3D element has less than 3 nodes... Bug!");
  }

  // Initialization of nodal coordinates
  LINALG::Matrix<nsd,numnode> nodecoords(true);  // node coordinates of the element
  for (int elenodeid=0;elenodeid<ele->NumNode();elenodeid++)
  {
    DRT::Node* currnode = discret_->gNode(ele->NodeIds()[elenodeid]);
    for (int i=0;i<nsd;i++)
      nodecoords(i,elenodeid) = currnode->X()[i];
  }

  iter = 0;   // iteration counter

  // evaluate shape functions at xi
  DRT::UTILS::shape_function_3D(shapeFcn, xi(0),xi(1),xi(2),DISTYPE);
  // evaluate derivative of shape functions at xi
  DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1, xi(0),xi(1),xi(2),DISTYPE);

  xjm.MultiplyNT(shapeXiDeriv1,nodecoords);   // jacobian
  xji.Invert(xjm);       // jacobian inverted

  residuum.Clear();
  residuum.Multiply(nodecoords, shapeFcn); //x_i
  residuum -= x;               // x_i - x
  residuum.Scale(-1.0);        // negative residuum -> RHS for Newton iteration

  // Newton loop to find local coordinates
  while(iter < max_iter)
  {
    iter += 1;

    //solve Newton iteration
    incr.Clear();
    incr.MultiplyTN(xji,residuum); // J^(-T)*residuum

    // update iteration
    for (int i=0;i<nsd;i++)
      xi(i) += incr(i);

    //=============== update residuum================
    //new shape functions and derivatives
    DRT::UTILS::shape_function_3D(shapeFcn, xi(0),xi(1),xi(2),DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1, xi(0),xi(1),xi(2),DISTYPE);

    // reset residual
    residuum.Clear();
    residuum.Multiply(nodecoords, shapeFcn); // x_i
    residuum -= x;               // x_i - x
    residuum.Scale(-1.0);        // negative residuum -> RHS for Newton iteration (new step)

    // convergence criterion
    if (incr.Norm2()/xi.Norm2() < relTolIncr && residuum.Norm2()/xi.Norm2() < relTolRes)
      break;

    // stop if point is too far away from element
    if (xi.Norm2()>3)
      break;

    // update jacobian and inverted jacobian
    xjm.MultiplyNT(shapeXiDeriv1,nodecoords);
    xji.Invert(xjm);
  } // end Newton loop

  // did newton iteration converge?
  if(iter == max_iter)
    cout << "WARNING: Newton iteration for finding local coordinates didn't converge!" << endl;

  // point near the element? (very high tolerance because of computational error)
  if (xi(0)<=1.0+XiTol && xi(1)<=1.0+XiTol && xi(2)<=1.0+XiTol
      && xi(0)>=-1.0-XiTol && xi(1)>=-1.0-XiTol && xi(2)>=-1.0-XiTol)
    elefound = true;

  if (elefound)
  {
    vel.Clear();
    const int* elenodeids = ele->NodeIds();  // nodeids of element

    LINALG::Matrix<numnode,1> nodephi(true); // nodal phivalues
    for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
      nodephi(nodeid,0) = (*phinpi_)[discret_->gNode(elenodeids[nodeid])->LID()];

    // get phivalue of point
    phi = nodephi.Dot(shapeFcn);

    LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true); // dummy
    LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true); // dummy

    // initialize nodal vectors
    LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
    LINALG::Matrix<1,2*numnode> nodepres(true); // node pressures, just for function call


    // evaluate data for the given point
#ifdef COMBUST_NORMAL_ENRICHMENT

    ApproxFuncNormalVector<2,2*numnode> shp;
    LINALG::Matrix<1,numnode> nodevelenr(true); // nodal enrichment values of normal enrichment

    pointdataXFEMNormal<numnode,DISTYPE>(
        ele,
#ifdef COLLAPSE_FLAME_NORMAL
        x,
#endif
        xi,
        xji,
        shapeFcn,
        enrShapeFcnPres,
        enrShapeXYPresDeriv1,
        shp,
        false
    );

    elementsNodalData<numnode>(
        ele,
        veln_,
        olddofman_,
        olddofcolmap_,
        oldNodalDofColDistrib_,
        nodevel,
        nodevelenr,
        nodepres); // nodal data of the element

    size_t velncounter = 0;
    for (size_t inode=0; inode<numnode; ++inode)
    {
      // standard shape functions are identical for all vector components
      // shp.velx.d0.s == shp.vely.d0.s == shp.velz.d0.s
      vel(0) += nodevel(0,inode)*shp.velx.d0.s(inode);
      vel(1) += nodevel(1,inode)*shp.vely.d0.s(inode);
      vel(2) += nodevel(2,inode)*shp.velz.d0.s(inode);

      const int gid = elenodeids[inode];
      const set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(gid);

      for (set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        if (enrfield->getField() == XFEM::PHYSICS::Veln)
        {
          vel(0) += nodevelenr(0,velncounter)*shp.velx.d0.n(velncounter);
          vel(1) += nodevelenr(0,velncounter)*shp.vely.d0.n(velncounter);
          vel(2) += nodevelenr(0,velncounter)*shp.velz.d0.n(velncounter);

          velncounter += 1;
        }
      }
    }
#else
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

    elementsNodalData<numnode>(
        ele,
        veln_,
        olddofman_,
        olddofcolmap_,
        oldNodalDofColDistrib_,
        nodevel,
        nodepres); // nodal data of the element

    // interpolate velocity and pressure values at starting point
    vel.Multiply(nodevel, enrShapeFcnVel);
//    cout << "element in which startpoint lies in is:\n" << *ele;
//    cout << "nodal velocities  are " << nodevel;
//    cout << "startpoint approximation: " << x;
//    cout << "in local xi-coordinates: " << xi;
//    cout << "nodal shape function values are " << shapeFcn;
//    cout << "nodal enr function values for vel in local coords are " << enrShapeFcnVel;
//    cout << "interpolated velocity at this startpoint is " << vel;
#endif
  } // end if fitting ele found
} // end function findElementAndLocalCoords



/*------------------------------------------------------------------------------------------------*
 * extract nodal pressures and velocities for element nodes                      winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode>
void XFEM::Startvalues::elementsNodalData(
    DRT::Element*& element,
    const RCP<Epetra_Vector> field, // Epetra_Vector fitting to field and dofDistribution
    const RCP<DofManager> dofman,
    const Epetra_Map& dofMap, // DofMap fitting to field and dofDistribution
    const map<DofKey<onNode>, DofGID>& dofDistribution, // dofDistribution fitting to Epetra_Vector and field
    LINALG::Matrix<3,2*numnode>& nodevel,
#ifdef COMBUST_NORMAL_ENRICHMENT
    LINALG::Matrix<1,numnode>& nodevelenr,
#endif
    LINALG::Matrix<1,2*numnode>& nodepres
) const
{
  nodevel.Clear();
#ifdef COMBUST_NORMAL_ENRICHMENT
    nodevelenr.Clear();
#endif
  nodepres.Clear();

  const int* elenodeids = element->NodeIds();  // nodegids of element nodes

  int dofcounterVelx = 0;
  int dofcounterVely = 0;
  int dofcounterVelz = 0;
  int dofcounterPres = 0;
#ifdef COMBUST_NORMAL_ENRICHMENT
  int dofcounterVeln = 0;
#endif

  for (int nodeid=0;nodeid<element->NumNode();nodeid++) // loop over element nodes
  {
    // get nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldEnrSet(dofman->getNodeDofSet(elenodeids[nodeid]));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
        fieldenr != fieldEnrSet.end();++fieldenr)
    {
      const DofKey<onNode> dofkey(elenodeids[nodeid], *fieldenr);
      const int dofpos = dofDistribution.find(dofkey)->second;
      switch (fieldenr->getEnrichment().Type())
      {
      case XFEM::Enrichment::typeStandard :
      case XFEM::Enrichment::typeJump :
      case XFEM::Enrichment::typeVoidFSI :
      case XFEM::Enrichment::typeVoid :
      case XFEM::Enrichment::typeKink :
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          nodevel(0,dofcounterVelx) = (*field)[dofMap.LID(dofpos)];
          dofcounterVelx++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          nodevel(1,dofcounterVely) = (*field)[dofMap.LID(dofpos)];
          dofcounterVely++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          nodevel(2,dofcounterVelz) = (*field)[dofMap.LID(dofpos)];
          dofcounterVelz++;
        }
#ifdef COMBUST_NORMAL_ENRICHMENT
        else if (fieldenr->getField() == XFEM::PHYSICS::Veln)
        {
          nodevelenr(0,dofcounterVeln) = (*field)[dofMap.LID(dofpos)];
          dofcounterVeln++;
        }
#endif
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          nodepres(0,dofcounterPres) = (*field)[dofMap.LID(dofpos)];
          dofcounterPres++;
        }
        else
        {
          cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << endl;
          dserror("not implemented physical field!");
        }
        break;
      }
      case XFEM::Enrichment::typeUndefined :
      default :
      {
        cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << endl;
        dserror("unknown enrichment type");
        break;
      }
      } // end switch enrichment
    } // end loop over fieldenr
  } // end loop over element nodes
}



#ifndef COMBUST_NORMAL_ENRICHMENT
/*------------------------------------------------------------------------------------------------*
 * compute data for an arbitrary point lying in a given element                  winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::Startvalues::pointdataXFEM(
    DRT::Element*& element,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,3>& xji,
    LINALG::Matrix<numnode,1>& shapeFcn,
    LINALG::Matrix<2*numnode,1>& enrShapeFcnVel,
    LINALG::Matrix<2*numnode,1>& enrShapeFcnPres,
    LINALG::Matrix<3,2*numnode>& enrShapeXYVelDeriv1,
    LINALG::Matrix<3,2*numnode>& enrShapeXYPresDeriv1,
    bool currentTimeStep
) const
{
  const int* elenodeids = element->NodeIds();  // nodeids of element
  const int nsd = 3;

  RCP<XFEM::DofManager> dofman;
  RCP<Epetra_Vector> phi;
  if (currentTimeStep)
  {
    dofman = dofman_;
    phi = phinpip_;
  }
  else
  {
    dofman = olddofman_;
    phi = phinpi_;
  }

  // clear data that should be filled
  shapeFcn.Clear();
  enrShapeFcnVel.Clear();
  enrShapeXYVelDeriv1.Clear();
  enrShapeFcnPres.Clear();
  enrShapeXYPresDeriv1.Clear();

  LINALG::Matrix<nsd,numnode> nodecoords(true); // node coordinates of the element
  for (size_t nodeid=0;nodeid<numnode;nodeid++)
  {
    DRT::Node* currnode = discret_->gNode(elenodeids[nodeid]);
    for (int i=0;i<nsd;i++)
      nodecoords(i,nodeid) = currnode->X()[i];
  }

  LINALG::Matrix<nsd,numnode> shapeXiDeriv1(true);
  // evaluate shape functions at xi
  DRT::UTILS::shape_function_3D(shapeFcn, xi(0),xi(1),xi(2),DISTYPE);
  // evaluate derivative of shape functions at xi
  DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1, xi(0),xi(1),xi(2),DISTYPE);

  LINALG::Matrix<nsd,nsd> xjm(true);
  xjm.MultiplyNT(shapeXiDeriv1,nodecoords);   // jacobian J = (dx/dxi)^T
  xji.Invert(xjm);       // jacobian inverted J^(-1) = dxi/dx

  LINALG::Matrix<nsd,numnode> shapeXYDeriv1(true); // first derivation of global shape functions
  shapeXYDeriv1.Multiply(xji,shapeXiDeriv1); // (dN/dx)^T = (dN/dxi)^T * J^(-T)
//  cout << "deriv1 of global shape fcn is " << shapeXYDeriv1 << endl;

  // second  derivative of (enriched) shape function in local and global coordinates
  LINALG::Matrix<2*nsd,numnode> shapeXYDeriv2(true); // just needed for function call
  LINALG::Matrix<2*nsd,2*numnode> enrShapeXYDeriv2(true); // just needed for function call so just one for vel and pres

  // nodal phivalues
  LINALG::Matrix<numnode,1> nodephi(true);
  for (int nodeid=0;nodeid<element->NumNode();nodeid++) // loop over element nodes
    nodephi(nodeid,0) = (*phi)[discret_->gNode(elenodeids[nodeid])->LID()];

  // create an element dof manager
  const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty; // ansatz map needed for eledofman
  Teuchos::RCP<XFEM::ElementDofManager> eleDofManager = rcp(new XFEM::ElementDofManager(*element,element_ansatz_empty,*dofman));

  // the enrichment functions may depend on the point
  // therefore the computation of the enrichment functions is called here
  // the gauss point is contained in shapeFcn!
  const XFEM::ElementEnrichmentValues enrvals(
      *element,*eleDofManager,nodephi,
      shapeFcn,shapeXYDeriv1,shapeXYDeriv2);

  // enriched shape functions and derivatives for nodal parameters (dofs)
  enrvals.ComputeModifiedEnrichedNodalShapefunction(
      XFEM::PHYSICS::Velx,shapeFcn,shapeXYDeriv1,shapeXYDeriv2,
      enrShapeFcnVel,enrShapeXYVelDeriv1,enrShapeXYDeriv2); // enrichment assumed to be equal for the 3 dimensions
  enrvals.ComputeModifiedEnrichedNodalShapefunction(
      XFEM::PHYSICS::Pres,shapeFcn,shapeXYDeriv1,shapeXYDeriv2,
      enrShapeFcnPres,enrShapeXYPresDeriv1,enrShapeXYDeriv2);
}
#endif



#ifdef COMBUST_NORMAL_ENRICHMENT
/*------------------------------------------------------------------------------------------------*
 * compute data for an arbitrary point lying in a given element                  winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::Startvalues::pointdataXFEMnormal(
    DRT::Element*& element,
#ifdef COLLAPSE_FLAME_NORMAL
    LINALG::Matrix<3,1> xyz,
#endif
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,3>& xji,
    LINALG::Matrix<numnode,1>& funct,
    LINALG::Matrix<2*numnode,1>& enrShapeFcnPres,
    LINALG::Matrix<3,2*numnode>& enrShapeXYPresDeriv1,
    ApproxFuncNormalVector<2,2*numnode>& shp,
    bool currentTimeStep
) const
{
  dserror("normal enrichment not working! probably some error in the shape functions (wrong 1st global derivative in cut elements!");
  const int* elenodeids = element->NodeIds();  // nodeids of element

  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

  RCP<XFEM::DofManager> dofman;
  RCP<Epetra_Vector> phi;
  if (currentTimeStep)
  {
    dofman = dofman_;
    phi = phinpip_;
  }
  else
  {
    dofman = olddofman_;
    phi = phinpi_;
  }

  // clear data that should be filled
  funct.Clear();
  shp.Clear();

  LINALG::Matrix<numnode,1> nodephi(true); // nodal phivalues
  LINALG::Matrix<nsd,numnode> nodecoords(true); // node coordinates of the element
  for (size_t nodeid=0;nodeid<numnode;nodeid++)
  {
    DRT::Node* currnode = discret_->gNode(elenodeids[nodeid]);
    nodephi(nodeid,0) = (*phi)[currnode->LID()];
    for (int i=0;i<nsd;i++)
      nodecoords(i,nodeid) = currnode->X()[i];
  }

  LINALG::Matrix<nsd,numnode> shapeXiDeriv1(true);
  // evaluate shape functions at xi
  DRT::UTILS::shape_function_3D(funct, xi(0),xi(1),xi(2),DISTYPE);
  // evaluate derivative of shape functions at xi
  DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1, xi(0),xi(1),xi(2),DISTYPE);

  LINALG::Matrix<nsd,nsd> xjm(true);
  xjm.MultiplyNT(shapeXiDeriv1,nodecoords);   // jacobian J = (dx/dxi)^T
  xji.Invert(xjm);       // jacobian inverted J^(-1) = dxi/dx

  LINALG::Matrix<nsd,numnode> derxy(true); // first derivation of global shape functions
  derxy.Multiply(xji,shapeXiDeriv1); // (dN/dx)^T = (dN/dxi)^T * J^(-T)
//  cout << "deriv1 of global shape fcn is " << shapeXYDeriv1 << endl;

  // second  derivative of (enriched) shape function in local and global coordinates
  LINALG::Matrix<2*nsd,numnode> shapeXYDeriv2(true); // just needed for function call
  LINALG::Matrix<2*nsd,2*numnode> enrShapeXYDeriv2(true); // just needed for function call so just one for vel and pres

  // create an element dof manager
  const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty; // ansatz map needed for eledofman
  Teuchos::RCP<XFEM::ElementDofManager> eleDofManager = rcp(new XFEM::ElementDofManager(*element,element_ansatz_empty,*dofman));

  // the enrichment functions may depend on the point
  // therefore the computation of the enrichment functions is called here
  // the gauss point is contained in shapeFcn!
  const XFEM::ElementEnrichmentValues enrvals(
      *element,*eleDofManager,nodephi,
      funct,derxy,shapeXYDeriv2);

  // get pointer to vector holding SMOOTHED G-function values at the fluid nodes
  Teuchos::RCP<Epetra_MultiVector> gradphiVec;
//  bool oldTimeStep = true;
//  if (oldTimeStep)
//    gradphiVec = flamefront_->GradPhiOld();
//  else
  gradphiVec = flamefront_->GradPhi();

  if (gradphiVec == Teuchos::null)
    dserror("No gradient of phi computed!");

  std::vector<double> gradphi; ///< smoothed G-function gradient at n+1

  DRT::UTILS::ExtractMyNodeBasedValues(element, gradphi,*gradphiVec);

  LINALG::Matrix<nsd,numnode> egradphi;
  unsigned ipos;
  for (size_t iparam=0; iparam<numnode; ++iparam)
  {
    ipos = iparam*3;
    egradphi(0,iparam) = gradphi[ipos];
    egradphi(1,iparam) = gradphi[ipos+1];
    egradphi(2,iparam) = gradphi[ipos+2];
  }

#ifdef COLLAPSE_FLAME_NORMAL
  LINALG::Matrix<nsd,1> normal(true);
  normal(0) = xyz(0);
  normal(1) = xyz(1);
  normal.Scale(-1.0/normal.Norm2());
#endif
  // enriched shape functions and derivatives for nodal parameters (dofs)
  enrvals.ComputeNormalShapeFunction(
      funct,
      derxy,
      shapeXYDeriv2,
      egradphi,
#ifdef COLLAPSE_FLAME_NORMAL
      normal,
      xyz,
#endif
      shp
  );

  for (size_t iparam = 0; iparam < numnode; ++iparam)
  {
    shp.velx.d0.s(iparam) = funct(iparam);
    shp.vely.d0.s(iparam) = funct(iparam);
    shp.velz.d0.s(iparam) = funct(iparam);

    shp.velx.dx.s(iparam) = derxy(0,iparam);
    shp.vely.dx.s(iparam) = derxy(0,iparam);
    shp.velz.dx.s(iparam) = derxy(0,iparam);

    shp.velx.dy.s(iparam) = derxy(1,iparam);
    shp.vely.dy.s(iparam) = derxy(1,iparam);
    shp.velz.dy.s(iparam) = derxy(1,iparam);

    shp.velx.dz.s(iparam) = derxy(2,iparam);
    shp.vely.dz.s(iparam) = derxy(2,iparam);
    shp.velz.dz.s(iparam) = derxy(2,iparam);
  }

  enrvals.ComputeModifiedEnrichedNodalShapefunction(
      XFEM::PHYSICS::Pres,funct,derxy,shapeXYDeriv2,
      enrShapeFcnPres,enrShapeXYPresDeriv1,enrShapeXYDeriv2);
}
#endif



/*------------------------------------------------------------------------------------------------*
 * rewrite data for new iteration of slbt                                        winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::newIteration_prepare(
    vector<RCP<Epetra_Vector> > newRowVectors
)
{
  for (size_t index=0;index<done_->size();index++)
  {
    StartpointData& currData = (*done_)[index];
    currData.movNode_ = *discret_->gNode(currData.movNode_.Id());
    currData.changed_ = 1; // initialize as true!!!
    currData.searchedProcs_ = 1;
    currData.iter_ = 0;
    currData.velValues_.clear();
    currData.presValues_.clear();
  }

  newIteration_nodalData(newRowVectors); // data at t^n+1 not used in predictor
  newRowVectors.clear(); // no more needed
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::newIteration_nodalData(
    vector<RCP<Epetra_Vector> > newRowVectors
)
{
  const int nsd = 3;

  // data about column vectors required
  const Epetra_Map& newdofcolmap = *discret_->DofColMap();
  map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID> newNodalDofColDistrib;
  dofman_->fillNodalDofColDistributionMap(newNodalDofColDistrib);

  vector<RCP<Epetra_Vector> > newColVectors;

  for (size_t index=0;index<newRowVectors.size();index++)
  {
    RCP<Epetra_Vector> tmpColVector = rcp(new Epetra_Vector(newdofcolmap,true));
    newColVectors.push_back(tmpColVector);
    LINALG::Export(*newRowVectors[index],*newColVectors[index]);
  }

  // computed data
  vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1(static_cast<int>(newVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  vector<LINALG::Matrix<1,nsd> > presnpDeriv1(static_cast<int>(newVectors_.size()),LINALG::Matrix<1,nsd>(true));

  vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1Tmp(static_cast<int>(newVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  vector<LINALG::Matrix<1,nsd> > presnpDeriv1Tmp(static_cast<int>(newVectors_.size()),LINALG::Matrix<1,nsd>(true));

  for (size_t index=0;index<done_->size();index++)
  {
    StartpointData& currData = (*done_)[index];
    DRT::Node& currnode = currData.movNode_;
    LINALG::Matrix<nsd,1> coords(currnode.X());

    for (size_t i=0;i<newColVectors.size();i++)
    {
      velnpDeriv1[i].Clear();
      presnpDeriv1[i].Clear();
    }

    for (int iele=0;iele<currnode.NumElement();iele++)
    {
      DRT::Element* currele = currnode.Elements()[iele];
// compute velnpderiv and presnpderiv

      switch (currele->Shape())
      {
      case DRT::Element::hex8:
      {
        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
        computeNodalGradient<numnode,DRT::Element::hex8>(newColVectors,newdofcolmap,newNodalDofColDistrib,currele,coords,velnpDeriv1Tmp,presnpDeriv1Tmp);
      }
      break;
      case DRT::Element::hex20:
      {
        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
        computeNodalGradient<numnode,DRT::Element::hex20>(newColVectors,newdofcolmap,newNodalDofColDistrib,currele,coords,velnpDeriv1Tmp,presnpDeriv1Tmp);
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
    const int gid = currnode.Id();
    const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey<onNode> newdofkey(gid, *fieldenr);

      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
            const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
            currData.vel_(0) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
            const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
            currData.vel_(1) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
            const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
            currData.vel_(2) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
        }
      }
    } // end loop over fieldenr
//    cout << "final velnpderiv is " << velnpDeriv1[0] << " and presnpderiv is " << presnpDeriv1[0] << endl;

    for (size_t i=0;i<newColVectors.size();i++)
    {
      velnpDeriv1[i].Scale(1.0/currnode.NumElement());
      presnpDeriv1[i].Scale(1.0/currnode.NumElement());
    }
//    cout << "scaled velnpderiv is " << velnpDeriv1[0] << " and presnpderiv is " << presnpDeriv1[0] << endl;

    currData.velDeriv_ = velnpDeriv1;
    currData.presDeriv_ = presnpDeriv1;
//    cout << "after setting transportvel is " << currData.vel_ << ", velderiv is " << velnpDeriv1[0]
//         << " and presderiv is " << presnpDeriv1[0] << endl;
  }
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::Startvalues::computeNodalGradient(
    vector<RCP<Epetra_Vector> >& newColVectors,
    const Epetra_Map& newdofcolmap,
    map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID>& newNodalDofColDistrib,
    DRT::Element*& ele,
    LINALG::Matrix<3,1> coords,
    vector<LINALG::Matrix<3,3> >& velnpDeriv1,
    vector<LINALG::Matrix<1,3> >& presnpDeriv1
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
  LINALG::Matrix<nsd,1> vel(true);

  // dummies for function call
  LINALG::Matrix<nsd,nsd> xji(true);
  LINALG::Matrix<numnode,1> shapeFcn(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);

#ifdef COMBUST_NORMAL_ENRICHMENT
  LINALG::Matrix<1,numnode> nodevelenr(true);
  ApproxFuncNormalVector<2,2*numnode> shp;
#else
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
#endif

  double phi = 0.0;
  bool elefound = false;

  callElementSearch(ele,coords,xi,vel,phi,elefound);
//cout << "xi coordinates are " << xi << endl;

  if (!elefound)
    dserror("element of a row node not on same processor as node?! BUG!");

#ifdef COMBUST_NORMAL_ENRICHMENT
  pointdataXFEMNormal<numnode,DISTYPE>(
      ele,
#ifdef COLLAPSE_FLAME_NORMAL
      coords,
#endif
      xi,
      xji,
      shapeFcn,
      enrShapeFcnPres,
      enrShapeXYPresDeriv1,
      shp,
      true);
#else
  pointdataXFEM<numnode,DISTYPE>(
      ele,
      xi,
      xji,
      shapeFcn,
      enrShapeFcnVel,
      enrShapeFcnPres,
      enrShapeXYVelDeriv1,
      enrShapeXYPresDeriv1,
      true);
#endif
//cout << "shapefcnvel is " << enrShapeFcnVel << ", velderiv is " << enrShapeXYVelDeriv1 << " and presderiv is " << enrShapeXYPresDeriv1 << endl;
  for (size_t i=0;i<newColVectors.size();i++)
  {
#ifdef COMBUST_NORMAL_ENRICHMENT
    elementsNodalData<numnode>(
        ele,
        newColVectors[i],
        dofman_,
        newdofcolmap,
        newNodalDofColDistrib,
        nodevel,
        nodevelenr,
        nodepres);

    const int* nodeids = currele->NodeIds();
    size_t velncounter = 0;

    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    for (int inode = 0; inode < numnode; ++inode)
    {
      // standard shape functions are identical for all vector components
      // e.g. shp.velx.dx.s == shp.vely.dx.s == shp.velz.dx.s
      velnpDeriv1[i](0,0) += nodevel(0,inode)*shp.velx.dx.s(inode);
      velnpDeriv1[i](0,1) += nodevel(0,inode)*shp.velx.dy.s(inode);
      velnpDeriv1[i](0,2) += nodevel(0,inode)*shp.velx.dz.s(inode);

      velnpDeriv1[i](1,0) += nodevel(1,inode)*shp.vely.dx.s(inode);
      velnpDeriv1[i](1,1) += nodevel(1,inode)*shp.vely.dy.s(inode);
      velnpDeriv1[i](1,2) += nodevel(1,inode)*shp.vely.dz.s(inode);

      velnpDeriv1[i](2,0) += nodevel(2,inode)*shp.velz.dx.s(inode);
      velnpDeriv1[i](2,1) += nodevel(2,inode)*shp.velz.dy.s(inode);
      velnpDeriv1[i](2,2) += nodevel(2,inode)*shp.velz.dz.s(inode);

      const int gid = nodeids[inode];
      const std::set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(gid);

      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        if (enrfield->getField() == XFEM::PHYSICS::Veln)
        {
          velnpDeriv1[i](0,0) += nodevelenr(0,velncounter)*shp.velx.dx.n(velncounter);
          velnpDeriv1[i](0,1) += nodevelenr(0,velncounter)*shp.velx.dy.n(velncounter);
          velnpDeriv1[i](0,2) += nodevelenr(0,velncounter)*shp.velx.dz.n(velncounter);

          velnpDeriv1[i](1,0) += nodevelenr(0,velncounter)*shp.vely.dx.n(velncounter);
          velnpDeriv1[i](1,1) += nodevelenr(0,velncounter)*shp.vely.dy.n(velncounter);
          velnpDeriv1[i](1,2) += nodevelenr(0,velncounter)*shp.vely.dz.n(velncounter);

          velnpDeriv1[i](2,0) += nodevelenr(0,velncounter)*shp.velz.dx.n(velncounter);
          velnpDeriv1[i](2,1) += nodevelenr(0,velncounter)*shp.velz.dy.n(velncounter);
          velnpDeriv1[i](2,2) += nodevelenr(0,velncounter)*shp.velz.dz.n(velncounter);

          velncounter += 1;
        }
      }
    }

#else
    elementsNodalData<numnode>(
        ele,
        newColVectors[i],
        dofman_,
        newdofcolmap,
        newNodalDofColDistrib,
        nodevel,
        nodepres);
//if (i==0) cout << "nodevels are " << nodevel << " and nodepres are " << nodepres << endl;
    velnpDeriv1[i].MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1,1.0);
    presnpDeriv1[i].MultiplyNT(1.0,nodepres,enrShapeXYPresDeriv1,1.0);
//if (i==0) cout << "curr velnpderiv is " << velnpDeriv1[i] << " and presnpderiv is " << presnpDeriv1[i] << endl;
#endif
  }
}



/*------------------------------------------------------------------------------------------------*
 * compare interface side of two points in combustion                            winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::Startvalues::interfaceSideCompareCombust(
    double phi1,
    double phi2
) const
{
  if (interfaceSideCombust(phi1) == interfaceSideCombust(phi2)) return true;
  else return false;
}



/*------------------------------------------------------------------------------------------------*
 * identify interface side of a point in combustion                              winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
int XFEM::Startvalues::interfaceSideCombust(
    double phi
) const
{
  if (phi >= 0) return 1;
  else return -1;
}



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export start data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::exportStartData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // destination proc (the "next" one)
  int dest = myrank_+1;
  if(myrank_ == (numproc_-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank_-1;
  if(myrank_ == 0)
    source = numproc_-1;

  DRT::PackBuffer dataSend;

  // packing the data
  DRT::ParObject::AddtoPack(dataSend,(int)basic_->size());
  for (size_t ipoint=0;ipoint<basic_->size();ipoint++)
  {
    StartpointData& basic = (*basic_)[ipoint];

    DRT::ParObject::AddtoPack(dataSend,basic.movNode_.Id());
    DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(basic.movNode_.X()));
    DRT::ParObject::AddtoPack(dataSend,basic.movNode_.Owner());
    DRT::ParObject::AddtoPack(dataSend,basic.vel_);
    DRT::ParObject::AddtoPack(dataSend,basic.velDeriv_);
    DRT::ParObject::AddtoPack(dataSend,basic.presDeriv_);
    DRT::ParObject::AddtoPack(dataSend,basic.startpoint_);
    DRT::ParObject::AddtoPack(dataSend,basic.changed_);
    DRT::ParObject::AddtoPack(dataSend,basic.phiValue_);
    DRT::ParObject::AddtoPack(dataSend,basic.searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,basic.iter_);
    DRT::ParObject::AddtoPack(dataSend,basic.startGid_);
    DRT::ParObject::AddtoPack(dataSend,basic.startOwner_);
    DRT::ParObject::AddtoPack(dataSend,basic.dMin_);
  }

  dataSend.StartPacking();

  DRT::ParObject::AddtoPack(dataSend,(int)basic_->size());
  for (size_t ipoint=0;ipoint<basic_->size();ipoint++)
  {
    StartpointData& basic = (*basic_)[ipoint];

    DRT::ParObject::AddtoPack(dataSend,basic.movNode_.Id());
    DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(basic.movNode_.X()));
    DRT::ParObject::AddtoPack(dataSend,basic.movNode_.Owner());
    DRT::ParObject::AddtoPack(dataSend,basic.vel_);
    DRT::ParObject::AddtoPack(dataSend,basic.velDeriv_);
    DRT::ParObject::AddtoPack(dataSend,basic.presDeriv_);
    DRT::ParObject::AddtoPack(dataSend,basic.startpoint_);
    DRT::ParObject::AddtoPack(dataSend,basic.changed_);
    DRT::ParObject::AddtoPack(dataSend,basic.phiValue_);
    DRT::ParObject::AddtoPack(dataSend,basic.searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,basic.iter_);
    DRT::ParObject::AddtoPack(dataSend,basic.startGid_);
    DRT::ParObject::AddtoPack(dataSend,basic.startOwner_);
    DRT::ParObject::AddtoPack(dataSend,basic.dMin_);
  }

  vector<int> lengthSend(1,0);
  lengthSend[0] = dataSend().size();
  int size_one = 1;

#ifdef DEBUG
  cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
#endif

  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter_.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
  // ... and receive length
  vector<int> lengthRecv(1,0);
  exporter_.Receive(source, length_tag, lengthRecv, size_one);
  exporter_.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter_.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);
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
  unsigned numberOfNodes = 0;

  // clear vector that should be filled
  basic_->clear();

  // unpack received data
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
  for (size_t inode=0;inode<numberOfNodes;inode++)
  {
    int gid;
    LINALG::Matrix<nsd,1> coords;
    int owner;
    LINALG::Matrix<nsd,1> vel;
    vector<LINALG::Matrix<nsd,nsd> > velDeriv;
    vector<LINALG::Matrix<1,nsd> > presDeriv;
    LINALG::Matrix<nsd,1> startpoint;
    int changed;
    double phiValue;
    int searchedProcs;
    int iter;
    vector<int> startGid;
    vector<int> startOwner;
    double dMin;

    DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,changed);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,iter);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,dMin);

    double coordinates[nsd];
    for (int dim=0;dim<nsd;dim++) coordinates[dim] = coords(dim);

    DRT::Node movNode(gid,coordinates,owner);

    basic_->push_back(StartpointData(
        movNode,
        vel,
        velDeriv,
        presDeriv,
        startpoint,
        changed,
        phiValue,
        searchedProcs,
        iter,
        startGid,
        startOwner,
        dMin));
  }

  // processors wait for each other
  discret_->Comm().Barrier();

} // end exportStartData



void XFEM::Startvalues::exportAlternativAlgoData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  vector<vector<StartpointData> > failedVec(numproc_);

  // fill vectors with the data
  for (size_t inode=0;inode<failed_->size();inode++)
  {
    StartpointData failed = (*failed_)[inode];
    for (size_t i=0;i<(*failed_)[inode].startOwner_.size();i++)
    {
      failed.startGid_.clear();
      failed.startGid_.push_back((*failed_)[inode].startGid_[i]);
      failed.startOwner_.clear();
      failed.startOwner_.push_back((*failed_)[inode].startOwner_[i]);
      failedVec[failed.startOwner_[0]].push_back(failed);
    }
  }

  (*failed_).clear(); // clear final failed vector
  (*failed_) = failedVec[myrank_]; // set entries of own proc
  failedVec[myrank_].clear(); // clear the set data from the vector

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    DRT::PackBuffer dataSend; // vector including all data that has to be send to dest proc
    vector<int> lengthSend(1,0);
    int size_one = 1;

    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;

    // pack data to be sent
    DRT::ParObject::AddtoPack(dataSend,(int)failedVec[dest].size());
    for (size_t inode=0;inode<failedVec[dest].size();inode++)
    {
      StartpointData failed = failedVec[dest][inode];

      DRT::ParObject::AddtoPack(dataSend,failed.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(failed.movNode_.X()));
      DRT::ParObject::AddtoPack(dataSend,failed.movNode_.Owner());
      DRT::ParObject::AddtoPack(dataSend,failed.vel_);
      DRT::ParObject::AddtoPack(dataSend,failed.velDeriv_);
      DRT::ParObject::AddtoPack(dataSend,failed.presDeriv_);
      DRT::ParObject::AddtoPack(dataSend,failed.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,failed.phiValue_);
      DRT::ParObject::AddtoPack(dataSend,failed.startGid_);
    }

    dataSend.StartPacking();

    DRT::ParObject::AddtoPack(dataSend,(int)failedVec[dest].size());
    for (size_t inode=0;inode<failedVec[dest].size();inode++)
    {
      StartpointData failed = failedVec[dest][inode];

      DRT::ParObject::AddtoPack(dataSend,failed.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(failed.movNode_.X()));
      DRT::ParObject::AddtoPack(dataSend,failed.movNode_.Owner());
      DRT::ParObject::AddtoPack(dataSend,failed.vel_);
      DRT::ParObject::AddtoPack(dataSend,failed.velDeriv_);
      DRT::ParObject::AddtoPack(dataSend,failed.presDeriv_);
      DRT::ParObject::AddtoPack(dataSend,failed.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,failed.phiValue_);
      DRT::ParObject::AddtoPack(dataSend,failed.startGid_);
    }

    // clear the no more needed data
    failedVec[dest].clear();

    lengthSend[0] = dataSend().size();

#ifdef DEBUG
    cout << "--- sending " << lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
# endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter_.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    vector<int> lengthRecv(1,0);
    exporter_.Receive(source, length_tag, lengthRecv, size_one);
    exporter_.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter_.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter_.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter_.Wait(req_data);

#ifdef DEBUG
    cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << endl;
#endif

    // pointer to current position of group of cells in global string (counts bytes)
    size_t posinData = 0;

    // unpack received data
    unsigned numberOfNodes;
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
    for (size_t inode=0;inode<numberOfNodes;inode++)
    {
      int gid;
      LINALG::Matrix<nsd,1> coords;
      int owner;
      LINALG::Matrix<nsd,1> vel;
      vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      vector<int> startGid;

      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);

      double coordinates[nsd];
      for (int dim=0;dim<nsd;dim++)
        coordinates[dim] = coords(dim);

      DRT::Node movNode(gid,coordinates,owner);

      failed_->push_back(StartpointData(
          movNode,
          vel,
          velDeriv,
          presDeriv,
          startpoint,
          phiValue,
          startGid,
          vector<int>(1,myrank_))); // startOwner is current proc
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors
} // end exportAlternativAlgoData



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::exportIterData(
  bool& procfinished
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  vector<int> lengthSend(1,0);
  int size_one = 1;

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

    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procfinished));
    dataSend.StartPacking();
    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procfinished));

    lengthSend[0] = dataSend().size();

    // send length of the data to be received ...
    int length_tag = 0;
    MPI_Request req_length_data;
    exporter_.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    vector<int> lengthRecv(1,0);
    exporter_.Receive(source, length_tag, lengthRecv, size_one);
    exporter_.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter_.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter_.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter_.Wait(req_data);

    // pointer to current position of group of cells in global string (counts bytes)
    size_t posinData = 0;
    int procfinishedNew;

    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,procfinishedNew);

    if (procfinishedNew==0)
    {
      procfinished = 0;
    }

    // processors wait for each other
    discret_->Comm().Barrier();
  }



/*--------------------------------------*
 * second part: if not all procs have   *
 * finished send data to neighbour proc *
 *--------------------------------------*/
  if (!procfinished)
  {
    DRT::PackBuffer dataSend;

    // packing the data
    DRT::ParObject::AddtoPack(dataSend,(int)next_->size());
    for (size_t ipoint=0;ipoint<next_->size();ipoint++)
    {
      StartpointData next = (*next_)[ipoint];

      DRT::ParObject::AddtoPack(dataSend,next.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(next.movNode_.X()));
      DRT::ParObject::AddtoPack(dataSend,next.movNode_.Owner());
      DRT::ParObject::AddtoPack(dataSend,next.vel_);
      DRT::ParObject::AddtoPack(dataSend,next.velDeriv_);
      DRT::ParObject::AddtoPack(dataSend,next.presDeriv_);
      DRT::ParObject::AddtoPack(dataSend,next.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,next.changed_);
      DRT::ParObject::AddtoPack(dataSend,next.phiValue_);
      DRT::ParObject::AddtoPack(dataSend,next.searchedProcs_);
      DRT::ParObject::AddtoPack(dataSend,next.iter_);
      DRT::ParObject::AddtoPack(dataSend,next.startGid_);
      DRT::ParObject::AddtoPack(dataSend,next.startOwner_);
    }

    dataSend.StartPacking();

    DRT::ParObject::AddtoPack(dataSend,(int)next_->size());
    for (size_t ipoint=0;ipoint<next_->size();ipoint++)
    {
      StartpointData next = (*next_)[ipoint];

      DRT::ParObject::AddtoPack(dataSend,next.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(next.movNode_.X()));
      DRT::ParObject::AddtoPack(dataSend,next.movNode_.Owner());
      DRT::ParObject::AddtoPack(dataSend,next.vel_);
      DRT::ParObject::AddtoPack(dataSend,next.velDeriv_);
      DRT::ParObject::AddtoPack(dataSend,next.presDeriv_);
      DRT::ParObject::AddtoPack(dataSend,next.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,next.changed_);
      DRT::ParObject::AddtoPack(dataSend,next.phiValue_);
      DRT::ParObject::AddtoPack(dataSend,next.searchedProcs_);
      DRT::ParObject::AddtoPack(dataSend,next.iter_);
      DRT::ParObject::AddtoPack(dataSend,next.startGid_);
      DRT::ParObject::AddtoPack(dataSend,next.startOwner_);
    }

    lengthSend[0] = dataSend().size();

#ifdef DEBUG
    cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
#endif

    // send length of the data to be received ...
    int length_tag = 0;
    MPI_Request req_length_data;
    exporter_.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    vector<int> lengthRecv(1,0);
    exporter_.Receive(source, length_tag, lengthRecv, size_one);
    exporter_.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter_.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter_.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter_.Wait(req_data);

#ifdef DEBUG
    cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << endl;
#endif


    // pointer to current position of group of cells in global string (counts bytes)
    size_t posinData = 0;

    // number of points sent
    unsigned numberOfNodes = 0;

    // clear next_-vector so it can be filled with new values
    next_->clear();

    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);

    for (size_t inode=0;inode<numberOfNodes;inode++) // loop over number of nodes to get
    {
      int gid;
      LINALG::Matrix<nsd,1> coords;
      int owner;
      LINALG::Matrix<nsd,1> vel;
      vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      int changed;
      double phiValue;
      int searchedProcs;
      int iter;
      vector<int> startGid;
      vector<int> startOwner;

      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,changed);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,iter);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);

      double coordinates[nsd];
      for (int dim=0;dim<nsd;dim++)
        coordinates[dim] = coords(dim);

      DRT::Node movNode(gid,coordinates,owner);

      next_->push_back(StartpointData(
          movNode,
          vel,
          velDeriv,
          presDeriv,
          startpoint,
          changed,
          phiValue,
          searchedProcs,
          iter,
          startGid,
          startOwner));
    } // end loop over number of points to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end if procfinished == false
} // end exportIterData



/*------------------------------------------------------------------------------------------------*
 * export final data to node's proc                                              winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::exportFinalData()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  vector<vector<StartpointData> > doneVec(numproc_);

  // fill vectors with the data
  for (size_t inode=0;inode<done_->size();inode++)
  {
    doneVec[(*done_)[inode].movNode_.Owner()].push_back((*done_)[inode]);
  }

  (*done_).clear(); // clear final vectors
  (*done_) = doneVec[myrank_]; // set final data of own processor
  doneVec[myrank_].clear(); // clear data about current proc

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    vector<int> lengthSend(1,0);
    int size_one = 1;

    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;

    DRT::PackBuffer dataSend;

    // pack data to be sent
    DRT::ParObject::AddtoPack(dataSend,(int)doneVec[dest].size());
    for (size_t inode=0;inode<doneVec[dest].size();inode++) // loop over number of nodes to be sent
    {
      StartpointData done = doneVec[dest][inode];

      DRT::ParObject::AddtoPack(dataSend,done.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,done.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,done.phiValue_);
      DRT::ParObject::AddtoPack(dataSend,done.startGid_);
      DRT::ParObject::AddtoPack(dataSend,done.startOwner_);
      DRT::ParObject::AddtoPack(dataSend,done.velValues_);
      DRT::ParObject::AddtoPack(dataSend,done.presValues_);
    }

    dataSend.StartPacking();

    DRT::ParObject::AddtoPack(dataSend,(int)doneVec[dest].size());
    for (size_t inode=0;inode<doneVec[dest].size();inode++) // loop over number of nodes to be sent
    {
      StartpointData done = doneVec[dest][inode];

      DRT::ParObject::AddtoPack(dataSend,done.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,done.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,done.phiValue_);
      DRT::ParObject::AddtoPack(dataSend,done.startGid_);
      DRT::ParObject::AddtoPack(dataSend,done.startOwner_);
      DRT::ParObject::AddtoPack(dataSend,done.velValues_);
      DRT::ParObject::AddtoPack(dataSend,done.presValues_);
    }

    // clear the no more needed data
    doneVec[dest].clear();

    lengthSend[0] = dataSend().size();

#ifdef DEBUG
    cout << "--- sending " << lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
# endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter_.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
    // ... and receive length
    vector<int> lengthRecv(1,0);
    exporter_.Receive(source, length_tag, lengthRecv, size_one);
    exporter_.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter_.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter_.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter_.Wait(req_data);

#ifdef DEBUG
    cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << endl;
#endif

    // pointer to current position of group of cells in global string (counts bytes)
    size_t posinData = 0;

    // unpack received data
    unsigned numberOfNodes;
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);

    for (size_t inode=0;inode<numberOfNodes;inode++) // loop over number of nodes to get
    {
      int gid;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      vector<int> startGid;
      vector<int> startOwner;
      vector<LINALG::Matrix<nsd,1> > velValues;
      vector<double> presValues;

      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velValues);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presValues);

      done_->push_back(StartpointData(
          *discret_->gNode(gid),
          startpoint,
          phiValue,
          startGid,
          startOwner,
          velValues,
          presValues));
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors
} // end exportfinalData

#endif // parallel

#endif // CCADISCRET



