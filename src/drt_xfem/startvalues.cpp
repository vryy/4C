/*!-----------------------------------------------------------------------------------------------*
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
    const RCP<XFEM::Enrichmentvalues> enrichvals,
    const RCP<InterfaceHandle> ih_npi,
    const Epetra_Map& olddofcolmap,
    const map<DofKey<onNode>, DofGID>& oldNodalDofColDistrib,
    const Epetra_Map& newdofrowmap,
    const map<DofKey<onNode>, DofGID>& newNodalDofRowDistrib,
    const double& dt,
    const double& theta,
    const double& flamespeed
) :
  veln_(veln),
  oldVectors_(oldVectors),
  newVectors_(newVectors),
  discret_(discret),
  olddofman_(olddofman),
  dofman_(dofman),
  oldJumpsAndKinks_(enrichvals->oldJumpAndKinkValues()),
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
    if (ih_npi->ElementBisected(iele) || ih_npi->ElementTouchedPlus(iele) || ih_npi->ElementTouchedMinus(iele))
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
        -1,
        -1,
        INFINITY));
    }
  } // end loop over processor nodes
  cout << " Computing new reference solution(s) for " << basic_->size() << " nodes..." << endl;
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrangian Back-Tracking main algorithm                                  winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::semiLagrangeBackTracking(
    vector<RCP<Epetra_Vector> > newRowVectors,
    bool predictor_step)
{
  // TODO restart testen
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element


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
    theta_curr_ = 0.5;//theta_default_; // standard theta if no predictor step

    prepareNewIteration(newRowVectors);

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
//    cout << "on proc " << myrank_ << " at beginning next size is " << next_->size() <<
//        ", curr size is " << curr_->size() << " and failed size is " << failed_->size() << endl;
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

        // if no startpoint found, interface is too far away from point
        // so that lagrangian point of view is senseless
        if (curr.changed_)
        {
          // Initialization
          DRT::Element* fittingele = NULL; // pointer to the element where start point lies in
          LINALG::Matrix<nsd,1> xi(true); // local transformed coordinates of x
          LINALG::Matrix<nsd,1> vel(true); // velocity of the start point approximation
          double phin = 0; // phi-value of the start point approximation
          bool elefound = false;             // true if an element for a point was found on the processor
          bool ihSide = true;             // true if point lies on correct side
          bool stdBackTracking = true; // true if standard back-tracking shall be done
          double deltaT1; // at t^n + deltaT1 the lagrangian point crosses the interface

          // search for an element where the current startpoint lies in
          // if found, give out all data at the startpoint
          elementAndLocalCoords(fittingele,curr.startpoint_,xi,vel,phin,elefound);

          // if element is not found, look at another processor and so add all
          // according data to the vectors which will be sent to the next processor
          if (!elefound)
          {
            if (curr.searchedProcs_ < numproc_)
            {
              StartpointData next(
                  curr.movNode_,
                  curr.vel_,
                  curr.velDeriv_,
                  curr.presDeriv_,
                  curr.startpoint_,
                  curr.changed_,
                  curr.phiValue_,
                  curr.searchedProcs_+1,
                  curr.iter_,
                  curr.startGid_,
                  curr.startOwner_);
              next_->push_back(next);
            }
            else // all procs searched -> point not in domain
            {
//            cout << "velderiv size is " << curr.velDeriv_.size() << endl;
              failed_->push_back(StartpointData(
                  curr.movNode_,
                  curr.vel_,
                  curr.velDeriv_,
                  curr.presDeriv_,
                  curr.startpoint_,
                  curr.phiValue_,
                  curr.startGid_,
                  curr.startOwner_));
              cout << "WARNING! Lagrangian start point not in domain!" << endl;
            }
          } // end if elefound

          // if element is found, the newton iteration to find a better startpoint can start
          else
          {
            if (interfaceSideCompareCombust(curr.phiValue_,phin) == false) // Lagrangian origin and original node on different interface sides
            {
//              if (!predictor_step_)
//              {
//                ihSide = true; // TODO change nomenclature... passt inhaltlich aktuell so, aber der name ist nicht treffend...
//                alternativeNewtonLoop(
//                    fittingele,
//                    curr,
//                    xi,
//                    vel,
//                    phin,
//                    elefound,
//                    ihSide,
//                    stdBackTracking,
//                    deltaT1);
//              }
//              else
                ihSide = false;
            }
            else  // Newton loop just for sensible points
            {
              ihSide = true;
              NewtonLoop(fittingele,curr,xi,vel,phin,elefound,ihSide,stdBackTracking,deltaT1);
            }

            // if iteration can go on (that is when startpoint is on
            // correct interface side and iter < max_iter)
            if (ihSide && curr.iter_<max_iter_)
            {

              // if element is not found in a newton step, look at another processor and so add
              // all according data to the vectors which will be sent to the next processor
              if (!elefound)
              {
                next_->push_back(StartpointData(
                    curr.movNode_,
                    curr.vel_,
                    curr.velDeriv_,
                    curr.presDeriv_,
                    curr.startpoint_,
                    curr.changed_,
                    curr.phiValue_,
                    2,
                    curr.iter_,
                    curr.startGid_,
                    curr.startOwner_));
              }

              // newton iteration converged to a good startpoint and so the data can be used to go on
              else
              {
                backTracking(
                    fittingele,curr,xi,stdBackTracking,deltaT1);
              }
            } // end if ihSide

            // if startpoint lies on wrong side at any time, this algorithm has to stop
            //for the searched node which is done by setting changed = false
            if (!ihSide)
              curr.changed_ = false;

            if (curr.iter_ == max_iter_) // maximum number of iterations reached
            {
              cout << "WARNING: start value set somehow sensible\n"
                      "  but newton iteration to find it didnt converge!" << endl;
              backTracking(
                  fittingele,curr,xi,stdBackTracking,deltaT1);
            }
          } // end if over nodes which changed interface
        } // end if startpoint changed
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
        failed_->push_back(StartpointData(
            curr.movNode_,
            curr.vel_,
            curr.velDeriv_,
            curr.presDeriv_,
            curr.startpoint_,
            curr.phiValue_,
            curr.startGid_,
            curr.startOwner_));
      }
    }

#ifdef PARALLEL
    // export nodes and according data for which the startpoint isn't still found (next_ vector) to next proc
//    cout << "on proc " << myrank_ << " before exporting next size is " << next_->size() <<
//        ", curr size is " << curr_->size() << " and failed size is " << failed_->size() << endl;
    exportIterData(procfinished);
//    cout << "on proc " << myrank_ << " after exporting next size is " << next_->size() <<
//        ", curr size is " << curr_->size() << " and failed size is " << failed_->size() << endl;

    // convergencecheck: procfinished == 1 just if all procs have finished
    if (procfinished)
      break;

    *curr_ = *next_;
    next_->clear();
  } // end while loop over searched nodes
#endif
//cout << "here after newton loops with failed size " << failed_->size() << endl;


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
//  cout << "failed size before exporting on proc " << myrank_ << " is " << failed_->size() << endl;
#ifdef PARALLEL
  exportAlternativAlgoData(); // export data of failed nodes
#endif
//  cout << "failed size after exporting on proc " << myrank_ << " is " << failed_->size() << endl;

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
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

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

        // just look for points on the same interface side
        if (interfaceSideCompareCombust(basic.phiValue_,(*phinpi_)[lnodeold->LID()]))
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
            double dist = diff.Norm2() + 2*(diff.Norm2()-arc.Norm2());
//            if (basic.movNode_.Id()==5650) cout << "diff is " << diff.Norm2() << " and angle is " << arc.Norm2() << endl;
            if (dist<basic.dMin_) // new nearest node (cosinus shall be near 1!)
            {
//if (basic.movNode_.Id()==5650) cout << "for " << basic.movNode_ << " new startnode is " << *lnodeold << endl;
              basic.startpoint_.Update(1.0,newNodeCoords,-dt_,nodevel);
              basic.dMin_ = dist;
              basic.startGid_ = *enrnode;
              basic.startOwner_ = myrank_;
              currChanged = true;
          } // end if new best startvalue
          }
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
  {//cout << (*basic_)[nodeid].movNode_ << " has startpoint " << (*basic_)[nodeid].startpoint_ << endl;
    if (!(*basic_)[nodeid].changed_)
    {
      dserror("WARNING! No point on one interface side found!\nThis indicates "
          "that the whole area is at one side of the interface!");
      break;
    }
  } // end loop over nodes
} // end startValuesFinder



/*------------------------------------------------------------------------------------------------*
 * Computing final data where semi-Lagrangian approach failed                    winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::getDataForNotConvergedNodes(
)
{
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

  // remark: all data has to be sent to the processor where
  //         the startpoint lies before calling this function
  for (size_t inode=0;inode<failed_->size();inode++)
  {
    StartpointData failed = (*failed_)[inode]; // data of node where sl approach failed

    DRT::Node* node = discret_->gNode(failed.startGid_); // failed node
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
    elementAndLocalCoords(
        nodesElement,x,xi,vel,phi,elefound);
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
      backTracking(
          nodesElement,failed,xi,true,0.0);
    } // end if elefound
  } // end loop over nodes
} // end getDataForNotConvergedNodes



/*------------------------------------------------------------------------------------------------*
 * alternativ Newton loop including least squares approach to determine time                      *
 * when particle moved over the interface side                                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::alternativeNewtonLoop(
    DRT::Element*& fittingele,
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

  dserror("currently out of use!");

  const size_t numnode = 8;  // 8 nodes for hex8
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
//#define version1
//#define version2
#define version3
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8; // only hex8 implemented
  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

  stdBackTracking = false; // alternative Newton loop -> alternative back tracking

  // Initialization
  LINALG::Matrix<numnode,1> shapeFcn(true);          // shape functions
  LINALG::Matrix<nsd,nsd> xji(true);        // invers of jacobian
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true); // dummy
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true); // dummy
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true); // dummy

  LINALG::Matrix<nsd,1> deltaV(true); // v^n - v^n+1_i
  LINALG::Matrix<nsd,1> deltaX(true); // x - x_{lagr,i} - deltaT*v^n+1
  double normDeltaV;
  LINALG::Matrix<1,nsd> deltaT1Deriv(true);
  LINALG::Matrix<nsd,nsd> velDeriv(true);
  LINALG::Matrix<1,1> tmpProduct(true);

  LINALG::Matrix<nsd,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<nsd,nsd> sysmat(true); // matrix for the newton system
  LINALG::Matrix<nsd,1> incr(true);                 // increment of the newton system

  LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
  LINALG::Matrix<1,2*numnode> nodepres(true); // node pressures, just for function call

#ifdef version2
  LINALG::Matrix<nsd,nsd> deltaXDeriv(true);
  LINALG::Matrix<nsd,nsd> tmpMat;
  LINALG::Matrix<nsd,1> tmpVec;
#endif



  const double relTolIncr = 1.0e-10;;  // tolerance for the increment
  const double relTolRes = 1.0e-10;    // tolerance for the residual

  LINALG::Matrix<nsd,1> origNodeCoords(true); // coordinates of endpoint of Lagrangian characteristics
  for (size_t i=0;i<nsd;i++)
    origNodeCoords(i) = curr.movNode_.X()[i];

  // initialize residual and required data
  cout << "x^n: " << origNodeCoords << "v(x^n): " << curr.vel_ << "x^n+1: "
      << curr.startpoint_ << "v(x^n+1): " << vel << endl;  deltaX.Update(1.0,origNodeCoords,-1.0,curr.startpoint_);
  deltaX.Update(-dt_,curr.vel_,1.0); cout << "deltaX is " << deltaX << endl;

  deltaV.Update(1.0,vel,-1.0,curr.vel_); cout << "deltaV is " << deltaV << endl;
  normDeltaV = deltaV.Norm2();

  LINALG::Matrix<nsd,1> normal(true); normal(0) = 1.0;
  tmpProduct.MultiplyTN(1.0/curr.vel_.Norm2(),normal,curr.vel_);

#ifdef version3
  // compute time when Lagrangian point crosses the interface
  {
    // get material from first (arbitrary!) element adjacent to this node
    const Teuchos::RCP<MAT::Material> matlistptr = fittingele->Material();
    dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());

    // density burnt domain
    Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
    dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
    const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
    const double rhoplus = matplus->Density();

    // density unburnt domain
    Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
    dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
    const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
    const double rhominus = matminus->Density();

    if (curr.phiValue_>=0)
    {
      deltaT1 = dt_-curr.phiValue_/(rhominus/rhoplus*flamespeed_*tmpProduct(0));
      cout << "deltaT1 is " << deltaT1 << endl;
    }
    else
    {
      deltaT1 = dt_-curr.phiValue_/(flamespeed_*tmpProduct(0));
      cout << "deltaT1 is " << deltaT1 << endl;
    }
  }
#endif

#ifdef version2
  deltaT1 = deltaV.Dot(deltaX)/(normDeltaV*normDeltaV); cout << "deltaT1 is " << deltaT1 << endl;

  velDeriv.MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1);
  for (size_t i=0;i<nsd;i++)
  {
    for (size_t j=0;j<nsd;j++)
    {
      if (i==j)
        tmpMat(i,j) = 1.0 + deltaT1*velDeriv(i,j);
      else
        tmpMat(i,j) = deltaT1*velDeriv(i,j);
    }
  }
  residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords); // dt*v(curr.startpoint_)
  residuum.Update(deltaT1,vel,dt_-deltaT1,curr.vel_,1.0);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)
  residuum.MultiplyTN(tmpMat,residuum);
#else
  residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords); // dt*v(curr.startpoint_)
  residuum.Update(deltaT1,vel,dt_-deltaT1,curr.vel_,1.0);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)
#endif

  while(curr.iter_ < max_iter_)  // newton loop
  {
    curr.iter_ += 1;

    { // build systemMatrix
      pointdataXFEM<nsd,numnode,DISTYPE>(
          fittingele,
          xi,
          xji,
          shapeFcn,
          enrShapeFcnVel,
          enrShapeFcnPres,
          enrShapeXYVelDeriv1,
          enrShapeXYPresDeriv1,
          phinpi_,
          olddofman_
      );

      elementsNodalData<nsd,numnode>(
          fittingele,
          veln_,
          olddofman_,
          olddofcolmap_,
          oldNodalDofColDistrib_,
          nodevel,
          nodepres); // nodal data of the element nodes

#ifdef version1
      velDeriv.MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1);
      deltaT1Deriv.UpdateT(-1.0,deltaV);
      deltaT1Deriv.MultiplyTN(1.0,deltaX,velDeriv,1.0);
      deltaT1Deriv.MultiplyTN(-2.0*deltaT1,deltaV,velDeriv);
      deltaT1Deriv.Scale(normDeltaV*normDeltaV);
      tmpProduct.Multiply(deltaT1Deriv,deltaV);
      for (size_t i=0;i<nsd;i++)
      {
        for (size_t j=0;j<nsd;j++)
        {
          if (i==j)
            sysmat(i,j) = 1.0 + tmpProduct(0); // I + dt*v_nodes*dN/dx
          else
            sysmat(i,j) = tmpProduct(0);
        }
      }
      sysmat.Update(deltaT1,velDeriv,1.0);
#endif
#ifdef version2

      for (size_t i=0;i<nsd;i++)
      {
        for (size_t j=0;j<nsd;j++)
        {
          if (i==j)
            deltaXDeriv(i,j) = 1.0 - deltaT1 * velDeriv(i,j); // I + dt*v_nodes*dN/dx
          else
            deltaXDeriv(i,j) = - deltaT1 * velDeriv(i,j);
        }
      }
      sysmat.Update(deltaXDeriv);
      sysmat.Multiply(deltaT1,velDeriv,deltaXDeriv,1.0);
#endif
#ifdef version3
      sysmat.Update(deltaT1,velDeriv);
      for (size_t i=0;i<nsd;i++)
        sysmat(i,i) += 1.0;
#endif
      sysmat.Invert(); cout << "residuum is " << residuum << " and sysmat is " << sysmat << endl;
    } // invers system Matrix built

    //solve Newton iteration
    incr.Clear();
    incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

    // update iteration
    for (size_t i=0;i<nsd;i++)
      curr.startpoint_(i) += incr(i);
//      cout << "in newton loop approximate startvalue is " << curr.startpoint_ << endl;

    //=============== update residuum================
    elementAndLocalCoords(fittingele, curr.startpoint_, xi, vel, phi, elefound);

    if (elefound) // element of curr.startpoint_ at this processor
    {
      if (interfaceSideCompareCombust(curr.phiValue_,phi) == true)
      {
#ifdef DEBUG
        cout << "phiOrigValue " << curr.phiValue_ << endl;
        cout << "phiApprValue " << phi << " und phiApprSign" << interfaceSideCombust(phi) << endl;
        cout << "phivalues have different sign\nlagragian origin lies on other side as target" << endl;
#endif
        ihSide = true;
        NewtonLoop(fittingele,curr,xi,vel,phi,elefound,ihSide,stdBackTracking,deltaT1);
        break; // leave newton loop if element is on wrong domain side
      }
      else
      {
        ihSide = true;

        // reset residual and required data
        cout << "x^n: " << origNodeCoords << "v(x^n): " << curr.vel_ << "x^n+1: "
            << curr.startpoint_ << "v(x^n+1): " << vel << endl;
        deltaX.Update(1.0,origNodeCoords,-1.0,curr.startpoint_);
        deltaX.Update(-dt_,curr.vel_,1.0); cout << "deltaX is " << deltaX << endl;

        deltaV.Update(1.0,vel,-1.0,curr.vel_); cout << "deltaV is " << deltaV << endl;
        normDeltaV = deltaV.Norm2();

        deltaT1 = dt_-curr.phiValue_/(flamespeed_*tmpProduct(0));cout << "deltaT1 is " << deltaT1 << endl;

#ifdef version2
//        velDeriv.MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1);
//        for (size_t i=0;i<nsd;i++)
//        {
//          for (size_t j=0;j<nsd;j++)
//          {
//            if (i==j)
//              tmpMat(i,j) = 1.0 + deltaT1*velDeriv(i,j);
//            else
//              tmpMat(i,j) = deltaT1*velDeriv(i,j);
//          }
//        }
//        residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords); // dt*v(curr.startpoint_)
//        residuum.Update(deltaT1,vel,dt_-deltaT1,curr.vel_,1.0);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)
//        residuum.MultiplyTN(tmpMat,residuum);
#else
        residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords); // dt*v(curr.startpoint_)
        residuum.Update(deltaT1,vel,dt_-deltaT1,curr.vel_,1.0);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)
#endif

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
}


/*------------------------------------------------------------------------------------------------*
 * evaluate element and data for a given point                                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::elementAndLocalCoords(
    DRT::Element*& fittingele,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    bool& elefound
) const
{
  elefound = false; // becomes true in the algo if element is found
  const size_t numnode = 8;  // 8 nodes for hex8
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8; // only hex8 implemented

  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

  // Initialization of what every following Newton loop will need
  DRT::Element* currele = NULL;                 // pointer to the currently viewed element
  LINALG::Matrix<nsd,numnode> nodecoords(true);  // node coordinates of the element

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
  //loop over elements
  for (int ieleid=(((fittingele==NULL) ? 0 : -1));ieleid<discret_->NumMyRowElements();ieleid++)
  {
    // if fittingele != NULL and so initialized,
    // first it should be checked if it is fitting
    if (ieleid == -1)
    {
      currele = fittingele;
      fittingele = NULL;
    }
    else
      currele = discret_->lRowElement(ieleid);

    // Initialization
    for (int elenodeid=0;elenodeid<currele->NumNode();elenodeid++)
    {
      DRT::Node* currnode = discret_->gNode(currele->NodeIds()[elenodeid]);
      for (size_t i=0;i<nsd;i++)
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
      for (size_t i=0;i<nsd;i++)
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
    {
      fittingele = currele;
      elefound = true;
      break;
    } // end if near the element
  } // end loop over processor elements

  if (fittingele != NULL) // fittingele found
  {
    LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
    LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true); // dummy
    LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true); // dummy
    LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true); // dummy

    // evaluate data for the given point
    pointdataXFEM<nsd,numnode,DISTYPE>(
        fittingele,
        xi,
        xji,
        shapeFcn,
        enrShapeFcnVel,
        enrShapeFcnPres,
        enrShapeXYVelDeriv1,
        enrShapeXYPresDeriv1,
        phinpi_,
        olddofman_
    );

    const int* elenodeids = fittingele->NodeIds();  // nodeids of element

    LINALG::Matrix<numnode,1> nodephi(true); // nodal phivalues
    for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
      nodephi(nodeid,0) = (*phinpi_)[discret_->gNode(elenodeids[nodeid])->LID()];

    // get phivalue of point
    phi = nodephi.Dot(shapeFcn);

    // initialize nodal vectors
    LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
    LINALG::Matrix<1,2*numnode> nodepres(true); // node pressures, just for function call

    elementsNodalData<nsd,numnode>(
        fittingele,
        veln_,
        olddofman_,
        olddofcolmap_,
        oldNodalDofColDistrib_,
        nodevel,
        nodepres); // nodal data of the element

    // interpolate velocity and pressure values at starting point
    vel.Multiply(nodevel, enrShapeFcnVel);

//    cout << "element in which startpoint lies in is:\n" << *fittingele;
//    cout << "nodal velocities  are " << nodevel;
//    cout << "startpoint approximation: " << x;
//    cout << "in local xi-coordinates: " << xi;
//    cout << "nodal shape function values are " << shapeFcn;
//    cout << "nodal enr function values for vel in local coords are " << enrShapeFcnVel;
//    cout << "interpolated velocity at this startpoint is " << vel;
  } // end if fitting ele found
} // end function findElementAndLocalCoords



/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::NewtonLoop(
    DRT::Element*& fittingele,
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
  const size_t numnode = 8;  // 8 nodes for hex8
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8; // only hex8 implemented
  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

  stdBackTracking = true; // standard Newton loop -> standard back tracking

  // Initialization
  LINALG::Matrix<numnode,1> shapeFcn(true);          // shape functions
  LINALG::Matrix<nsd,nsd> xji(true);        // invers of jacobian
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true); // dummy
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true); // dummy
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true); // dummy

  LINALG::Matrix<nsd,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<nsd,nsd> sysmat(true); // matrix for the newton system
  LINALG::Matrix<nsd,1> incr(true);                 // increment of the newton system

  LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
  LINALG::Matrix<1,2*numnode> nodepres(true); // node pressures, just for function call

  const double relTolIncr = 1.0e-10;;  // tolerance for the increment
  const double relTolRes = 1.0e-10;    // tolerance for the residual

  LINALG::Matrix<nsd,1> origNodeCoords(true); // coordinates of endpoint of Lagrangian characteristics
  for (size_t i=0;i<nsd;i++)
    origNodeCoords(i) = curr.movNode_.X()[i];

  // initialize residual (Theta = 0 at predictor step)
  residuum.Clear();
  residuum.Update((1.0-theta_curr_),vel,theta_curr_,curr.vel_); // dt*v(curr.startpoint_)
  residuum.Update(1.0,curr.startpoint_,-1.0,origNodeCoords,dt_);  // R = curr.startpoint_ - curr.movNode_ + dt*v(curr.startpoint_)

  while(curr.iter_ < max_iter_)  // newton loop
  {
    curr.iter_ += 1;

    { // build systemMatrix
      pointdataXFEM<nsd,numnode,DISTYPE>(
          fittingele,
          xi,
          xji,
          shapeFcn,
          enrShapeFcnVel,
          enrShapeFcnPres,
          enrShapeXYVelDeriv1,
          enrShapeXYPresDeriv1,
          phinpi_,
          olddofman_
      );

      elementsNodalData<nsd,numnode>(
          fittingele,
          veln_,
          olddofman_,
          olddofcolmap_,
          oldNodalDofColDistrib_,
          nodevel,
          nodepres); // nodal data of the element nodes

      sysmat.MultiplyNT((1.0-theta_curr_)*dt_,nodevel,enrShapeXYVelDeriv1); // (1-theta) * dt * v_nodes * dN/dx
      for (size_t i=0;i<nsd;i++)
        sysmat(i,i) += 1.0; // I + dt*v_nodes*dN/dx
      sysmat.Invert();
    } // invers system Matrix built

    //solve Newton iteration
    incr.Clear();
    incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

    // update iteration
    for (size_t i=0;i<nsd;i++)
      curr.startpoint_(i) += incr(i);
//      cout << "in newton loop approximate startvalue is " << curr.startpoint_ << endl;

    //=============== update residuum================
    elementAndLocalCoords(fittingele, curr.startpoint_, xi, vel, phi, elefound);

    if (elefound) // element of curr.startpoint_ at this processor
    {
      if (interfaceSideCompareCombust(curr.phiValue_,phi) == false)
      {
#ifdef DEBUG
        cout << "phiOrigValue " << curr.phiValue_ << endl;
        cout << "phiApprValue " << phi << " und phiApprSign" << interfaceSideCombust(phi) << endl;
        cout << "phivalues have different sign\nlagragian origin lies on other side as target" << endl;
#endif
//        if (predictor_step_)
          ihSide = false;
//        else
//        {
//          ihSide = true;
//          alternativeNewtonLoop(fittingele,curr,xi,vel,phi,elefound,ihSide,stdBackTracking,deltaT1);
//        }
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
 * back-tracking of data at final Lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::backTracking(
    DRT::Element*& fittingele,
    StartpointData& data,
    LINALG::Matrix<3,1>& xi,
    bool stdBackTracking,
    double deltaT1
)
{
  const int numnode = 8;
  const int nsd = 3;
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8;

  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

//  cout << node << "has lagrange origin " << x << "\nwith xi-coordinates "
//      << xi << "\n in element " << *fittingele << endl;
  const int* elenodeids = fittingele->NodeIds();  // nodeids of element

  LINALG::Matrix<nsd,nsd> xji(true); // invers of jacobian
  LINALG::Matrix<numnode,1> shapeFcn(true); // shape function

  // enriched shape functions and there derivatives in local coordinates (N * \Psi)
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);

  // data for the final back-tracking
  LINALG::Matrix<nsd,1> vel(true); // velocity data
  LINALG::Matrix<nsd,nsd> velnDeriv1(true); // first derivation of velocity data
  LINALG::Matrix<nsd,nsd> velnpDeriv1(true); // first derivation of velocity data
  LINALG::Matrix<1,1> pres(true); // pressure data
  LINALG::Matrix<1,nsd> presnDeriv1(true); // first derivation of pressure data
  LINALG::Matrix<1,nsd> presnpDeriv1(true); // first derivation of pressure data
  LINALG::Matrix<nsd,1> transportVeln(true); // transport velocity at Lagrangian origin (x_Lagr(t^n))

  pointdataXFEM<nsd,numnode,DISTYPE>(
      fittingele,
      xi,
      xji,
      shapeFcn,
      enrShapeFcnVel,
      enrShapeFcnPres,
      enrShapeXYVelDeriv1,
      enrShapeXYPresDeriv1,
      phinpi_,
      olddofman_
  );

  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd,2*numnode> nodevel(true);
  // node velocities of the element nodes for the data that should be changed
  vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,2*numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size(),LINALG::Matrix<1,2*numnode>(true));
  // nodal phivalues
  LINALG::Matrix<numnode,1> nodephi(true);
  for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
    nodephi(nodeid,0) = (*phinpi_)[discret_->gNode(elenodeids[nodeid])->LID()];

  int dofcounterVelx = 0;
  int dofcounterVely = 0;
  int dofcounterVelz = 0;
  int dofcounterPres = 0;

  for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
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
          nodevel(0,dofcounterVelx) = (*veln_)[olddofcolmap_.LID(olddofpos)];
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](0,dofcounterVelx) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVelx++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          nodevel(1,dofcounterVely) = (*veln_)[olddofcolmap_.LID(olddofpos)];
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](1,dofcounterVely) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVely++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
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

  transportVeln.Multiply(nodevel, enrShapeFcnVel);
  LINALG::Matrix<1,1> phinpi;
  phinpi.MultiplyTN(shapeFcn,nodephi);

  // computing pseudo time-step deltaT
  // remark: if x is the Lagrange-origin of node, deltaT = dt with respect to small errors.
  // if its not, deltaT estimates the time x needs to move to node)
  double deltaT = 0; // pseudo time-step
  {
    LINALG::Matrix<nsd,1> diff(data.movNode_.X());
    diff -= data.startpoint_; // diff = x_Node - x_Appr

    double numerator = transportVeln.Dot(diff); // numerator = v^T*(x_Node-x_Appr)
    double denominator = transportVeln.Dot(transportVeln); // denominator = v^T*v

    if (denominator>1e-15) deltaT = numerator/denominator; // else deltaT = 0 as initialized
  }

  vector<LINALG::Matrix<nsd,1> > velValues(oldVectors_.size(),LINALG::Matrix<nsd,1>(true)); // velocity of the data that should be changed
  vector<double> presValues(oldVectors_.size(),0); // pressures of the data that should be changed

  // interpolate velocity and pressure gradients for all fields at starting point and get final values
  for (size_t index=0;index<oldVectors_.size();index++)
  {
    vel.Clear();
    velnDeriv1.Clear();
    pres.Clear();
    presnDeriv1.Clear();

    velnDeriv1.MultiplyNT(nodeveldata[index],enrShapeXYVelDeriv1);
    presnDeriv1.MultiplyNT(nodepresdata[index],enrShapeXYPresDeriv1);
//    if (index==0)
//        cout << "backtracking data: velnderiv1 is " << velnDeriv1 << ", transportveln is "
//             << transportVeln << ",velnpderiv1 is " << data.velDeriv_[index] <<
//             ", transportvelnp is " << data.vel_ << "presnderiv1 is " << presnDeriv1 <<
//             " and presnpderiv1 is " << data.presDeriv_[index] << endl;
    if (stdBackTracking)
    {
      vel.Multiply(1.0-theta_curr_,velnDeriv1,transportVeln); // v = (1-theta)*Dv^n/Dx*v^n
      vel.Multiply(theta_curr_,data.velDeriv_[index],data.vel_,1.0); // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
      vel.Multiply(1.0,nodeveldata[index],enrShapeFcnVel,deltaT); // v = v_n + dt*(theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n)
      velValues[index]=vel;
//      if (index == 0) cout << "after std back-tracking vel is " << vel << endl;

      pres.Multiply(1.0-theta_curr_,presnDeriv1,transportVeln); // p = (1-theta)*Dp^n/Dx*v^n
      pres.Multiply(theta_curr_,data.presDeriv_[index],data.vel_,1.0); // p = theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n
      pres.Multiply(1.0,nodepresdata[index],enrShapeFcnPres,deltaT); // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
      presValues[index] = pres(0);
//      if (index == 0) cout << "after std back-tracking pres is " << pres(0) << endl;
    }
    else
    {
      vel.Multiply(deltaT1,velnDeriv1,transportVeln);
      vel.Multiply(deltaT-deltaT1,data.velDeriv_[index],data.vel_,1.0);
      vel.Multiply(1.0,nodeveldata[index],enrShapeFcnVel,1.0);
      LINALG::Matrix<nsd,1> vel_jump(true); vel_jump(0)=2.0; vel+=vel_jump; // TODO verallgemeinern
      velValues[index]=vel;
//      if (index == 0) cout << "after alter back-tracking vel is " << vel << endl;

      pres.Multiply(deltaT1,presnDeriv1,transportVeln);
      pres.Multiply(deltaT-deltaT1,data.presDeriv_[index],data.vel_,1.0);
      pres.Multiply(1.0,nodepresdata[index],enrShapeFcnPres,1.0);
      presValues[index] = pres(0)+8.0;//presJump+phinpi*presKink; // TODO verallgemeinern
      pres.Multiply(nodepresdata[index],enrShapeFcnPres);
//      if (index == 0) cout << "after alter back-tracking pres is " << pres(0) <<
//          ". it was " << pres << endl;
    }
  }

  done_->push_back(StartpointData(
      data.movNode_,
      data.startpoint_,
      data.phiValue_,
      data.startGid_,
      myrank_,
      velValues,
      presValues));
} // end getFinalStartvalues



/*------------------------------------------------------------------------------------------------*
 * setting the final data in Epetra Vector for a node                            winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::setFinalData(
) const
{
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

  for (size_t inode=0;inode<done_->size();inode++)
  {
    vector<LINALG::Matrix<nsd,1> > velValues((*done_)[inode].velValues_); // velocities of the node
    vector<double> presValues((*done_)[inode].presValues_); // pressures of the node

    const int gnodeid = (*done_)[inode].movNode_.Id(); // global node id

    // set nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gnodeid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey<onNode> newdofkey(gnodeid, *fieldenr);
      const int newdofpos = newNodalDofRowDistrib_.find(newdofkey)->second;

      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          for (size_t index=0;index<newVectors_.size();index++)
          {
//            if (index == 0)
//              cout << (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[index](0) << endl;
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = velValues[index](0);
          }
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          for (size_t index=0;index<newVectors_.size();index++)
          {
//            if (index == 0)
//              cout << (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[index](1) << endl;
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = velValues[index](1);
          }
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          for (size_t index=0;index<newVectors_.size();index++)
          {
//            if (index == 0)
//              cout << (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] << " becomes " << velValues[index](2) << endl;
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = velValues[index](2);
          }
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          for (size_t index=0;index<newVectors_.size();index++)
          {
//            if (index == 0)
//              cout << (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] << " becomes " << presValues[index] << endl;
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = presValues[index];
          }
        }
      }
    } // end loop over fieldenr
  } // end loop over nodes
} // end setFinalData



/*------------------------------------------------------------------------------------------------*
 * compare interface side of two points in combustion                            winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::Startvalues::interfaceSideCompareCombust(
    double phi1,
    double phi2
) const
{
  if (interfaceSideCombust(phi1) == interfaceSideCombust(phi2)) return true;
  else return false; // TODO how should an interface be handled which is completely on cell boundarys
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



/*------------------------------------------------------------------------------------------------*
 * extract nodal pressures and velocities for element nodes                      winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
template<size_t nsd, size_t numnode>
void XFEM::Startvalues::elementsNodalData(
    DRT::Element*& element,
    const RCP<Epetra_Vector> field, // Epetra_Vector fitting to field and dofDistribution
    const RCP<DofManager> dofman,
    const Epetra_Map& dofMap, // DofMap fitting to field and dofDistribution
    const map<DofKey<onNode>, DofGID>& dofDistribution, // dofDistribution fitting to Epetra_Vector and field
    LINALG::Matrix<nsd,2*numnode>& nodevel,
    LINALG::Matrix<1,2*numnode>& nodepres
) const
{
  nodevel.Clear();
  nodepres.Clear();

  const int* elenodeids = element->NodeIds();  // nodegids of element nodes

  int dofcounterVelx = 0;
  int dofcounterVely = 0;
  int dofcounterVelz = 0;
  int dofcounterPres = 0;

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



/*------------------------------------------------------------------------------------------------*
 * compute data for an arbitrary point lying in a given element                  winklmaier 10/10 *
 *------------------------------------------------------------------------------------------------*/
template<size_t nsd, size_t numnode, DRT::Element::DiscretizationType DISTYPE>
void XFEM::Startvalues::pointdataXFEM(
    DRT::Element*& element,
    LINALG::Matrix<nsd,1>& xi,
    LINALG::Matrix<nsd,nsd>& xji,
    LINALG::Matrix<numnode,1>& shapeFcn,
    LINALG::Matrix<2*numnode,1>& enrShapeFcnVel,
    LINALG::Matrix<2*numnode,1>& enrShapeFcnPres,
    LINALG::Matrix<nsd,2*numnode>& enrShapeXYVelDeriv1,
    LINALG::Matrix<nsd,2*numnode>& enrShapeXYPresDeriv1,
    RCP<Epetra_Vector> phi,
    RCP<XFEM::DofManager> dofman
) const
{
  const int* elenodeids = element->NodeIds();  // nodeids of element

  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

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
    for (size_t i=0;i<nsd;i++)
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



/*------------------------------------------------------------------------------------------------*
 * rewrite data for new iteration of slbt                                        winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::prepareNewIteration(
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

  nodalDataAtNewIteration(newRowVectors); // data at t^n+1 not used in predictor
  newRowVectors.clear(); // no more needed
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::nodalDataAtNewIteration(
    vector<RCP<Epetra_Vector> > newRowVectors
)
{
  const int nsd = 3;
  const int numnode = 8;
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8;
  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");

  // data about column vectors required
  const Epetra_Map& newdofcolmap = *discret_->DofColMap();
  map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID> newNodalDofColDistrib;
  dofman_->fillNodalDofColDistributionMap(newNodalDofColDistrib);

//cout << *newRowVectors[0] << endl;
//cout << endl<< endl << endl << *newRowVectors[1] << endl;
  vector<RCP<Epetra_Vector> > newColVectors;//(static_cast<int>(newVectors_.size()),rcp(new Epetra_Vector(newdofcolmap)));

  for (size_t index=0;index<newRowVectors.size();index++)
  {
    RCP<Epetra_Vector> tmpColVector = rcp(new Epetra_Vector(newdofcolmap,true));
    newColVectors.push_back(tmpColVector);
    LINALG::Export(*newRowVectors[index],*newColVectors[index]);
  }

  // computed data
  vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1(static_cast<int>(newVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
  vector<LINALG::Matrix<1,nsd> > presnpDeriv1(static_cast<int>(newVectors_.size()),LINALG::Matrix<1,nsd>(true));

  // shape fcn data
  LINALG::Matrix<nsd,1> xi(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);

  // nodal data
  LINALG::Matrix<nsd,2*numnode> nodevel(true);
  LINALG::Matrix<1,2*numnode> nodepres(true);
  LINALG::Matrix<nsd,1> nodecoords(true);
  LINALG::Matrix<nsd,1> vel(true);

  // dummies for function call
  LINALG::Matrix<nsd,nsd> xji(true);
  LINALG::Matrix<numnode,1> shapeFcn(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);
  double phi = 0.0;
  bool elefound = false;

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
      elementAndLocalCoords(currele,coords,xi,vel,phi,elefound);
//cout << "xi coordinates are " << xi << endl;
      pointdataXFEM<nsd,numnode,DISTYPE>(
          currele,
          xi,
          xji,
          shapeFcn,
          enrShapeFcnVel,
          enrShapeFcnPres,
          enrShapeXYVelDeriv1,
          enrShapeXYPresDeriv1,
          phinpip_,
          dofman_);
//cout << "shapefcnvel is " << enrShapeFcnVel << ", velderiv is " << enrShapeXYVelDeriv1 << " and presderiv is " << enrShapeXYPresDeriv1 << endl;
      for (size_t i=0;i<newColVectors.size();i++)
      {
        elementsNodalData<nsd,numnode>(
            currele,
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
        if (iele==0 && i==0) // velocity
          currData.vel_.Multiply(nodevel,enrShapeFcnVel);
      }
    }
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



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export start data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::exportStartData()
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
    int startGid;
    int startOwner;
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
    for (size_t dim=0;dim<nsd;dim++) coordinates[dim] = coords(dim);

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
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  vector<vector<StartpointData> > failedVec(numproc_);

  // fill vectors with the data
  for (size_t inode=0;inode<failed_->size();inode++)
  {
    failedVec[(*failed_)[inode].startOwner_].push_back((*failed_)[inode]);
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
      int startGid;

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
      for (size_t dim=0;dim<nsd;dim++)
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
          myrank_)); // startOwner is current proc
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
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

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
      int startGid;
      int startOwner;

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
      for (size_t dim=0;dim<nsd;dim++)
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
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

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
      int startGid;
      int startOwner;
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



