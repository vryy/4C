/*!
\file dof_distribution_switcher.cpp

\brief provides the dofmanager classes

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
#ifdef CCADISCRET

#include "dof_management.H"
#include "dof_management_element.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include <iostream>
#include "startvalues.H"


XFEM::Startvalues::Startvalues(
    const RCP<DRT::Discretization> discret,
    const RCP<DofManager> olddofman,
    const RCP<DofManager> dofman,
    vector<RCP<Epetra_Vector> > oldVectors,
    vector<RCP<Epetra_Vector> > newVectors,
    const RCP<Epetra_Vector> veln,
    const RCP<COMBUST::FlameFront> flamefront,
    const RCP<InterfaceHandle> ihold,
    Epetra_Map olddofcolmap,
    map<DofKey<onNode>, DofGID> oldNodalDofColDistrib,
    Epetra_Map newdofrowmap,
    map<DofKey<onNode>, DofGID> nodalDofRowDistrib
//    map<int,set<DofKey<onNode> > > oldEnrNodalDofSet,
//    map<DofKey<onNode>,DofGID> oldEnrNodalDofDistrib,
//    map<int,double> oldEnrValues
) :
  veln_(veln),
  oldVectors_(oldVectors),
  newVectors_(newVectors),
  ih_old_(ihold),
  discret_(discret),
  olddofman_(olddofman),
  dofman_(dofman),
  olddofcolmap_(olddofcolmap),
  oldNodalDofColDistrib_(oldNodalDofColDistrib),
  newdofrowmap_(newdofrowmap),
  nodalDofRowDistrib_(nodalDofRowDistrib),
  exporter_(discret_->Comm()),
  myrank_(discret_->Comm().MyPID()),
  numproc_(discret_->Comm().NumProc())
//  oldEnrNodalDofSet_(oldEnrNodalDofSet),
//  oldEnrNodalDofDistrib_(oldEnrNodalDofDistrib),
//  oldEnrValues_(oldEnrValues)
{
  
#ifdef PARALLEL
  //-------------------------------
  // prepare parallel communication
  //-------------------------------
  // get communicator, e.g. from fluid discretization
//  numproc_ = comm_->NumProc();
//  myrank_ = comm_->MyPID();
//  exporter_ = rcp(new DRT::Exporter(*comm_));
//  exporter_(*comm_);
#else
  const int numproc = 1;
#endif
  phiOldCol_ = flamefront->Phin();
  phiNewCol_ = flamefront->Phinp();
  // pointer to phivector at time n with rowmap
//  phiOldRow_ = rcp(new Epetra_Vector(*discret_->NodeRowMap()));
//  LINALG::Export(*phiOldCol_,*phiOldRow_);
//  // pointer to phivector at time n+1 with rowmap
//  phiNewRow_ = rcp(new Epetra_Vector(*discret_->NodeRowMap()));
//  LINALG::Export(*flamefront->Phinp(),*phiNewRow_);
  
//  cout << "phiold is " << *phiOldRow_;
//  cout << "phioldCol is " << *phiOldCol_;
  
  discret_->Comm().Barrier();
  // initialize empty structure vectors
  curr_.clear();
  next_.clear();
  done_.clear();
  failed_.clear();
  
  // fill curr_ structure with the data for the nodes which changed interface side
  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
  {
    // node which changed interface side (determined by sign of phi at time n and n+1)
    if ((discret_->lColNode(lnodeid)->Owner() == myrank_) &&
        ((((*phiNewCol_)[lnodeid]>=0) ? 1 : -1) != (((*phiOldCol_)[lnodeid]>=0) ? 1 : -1)))
    {
      StartpointData curr(*discret_->lColNode(lnodeid),LINALG::Matrix<3,1>(true),0,
          ((*phiNewCol_)[lnodeid]>=0) ? 1 : -1,1,0,-1,-1,INFINITY);
      curr_.push_back(curr);
    }
  }
  // TODO make a "side-compair function" so just this has to be changed for other applications
  
  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
  {
    if (discret_->lColNode(lnodeid)->Owner() == myrank_)
    {
      const int gid = discret_->lColNode(lnodeid)->Id();
      
      const set<XFEM::FieldEnr>& oldfieldenrset(olddofman_->getNodeDofSet(gid));
      for (set<XFEM::FieldEnr>::const_iterator oldfieldenr = oldfieldenrset.begin();
          oldfieldenr != oldfieldenrset.end();
          oldfieldenr++)
      {
        if (oldfieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
        {
          // not standard dof exists -> node enriched
          oldEnrNodes_.insert(gid); // TODO one set for every interface side would be better
          break;
        }
      }
    }
  } // end loop over processor nodes
  
  init_ = curr_;
  
//  cout << "number of enriched nodes at old timestep on proc " << myrank_ <<
//      " is " << oldEnrNodes_.size() << endl;
//  for (set<int>::const_iterator i=oldEnrNodes_.begin();
//      i != oldEnrNodes_.end();i++)
//    cout << "node " << *i << " was enriched at old timestep" << endl;
//  discret_->Comm().Barrier();
//  cout << "phinew is " << *phiNewRow_;
  
  return;
}

/*------------------------------------------------------------------------------------------------*
 * semi-lagrangian back-tracing method                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::semiLagrangeExtrapolation(
    const Teuchos::RCP<COMBUST::FlameFront>& flamefront_,
    const double& dt
)
{
  //TODO until now no attention to boundary conditions
  // TODO global startwerte genauer machen über lokales verfahren:
  //      phi = phi + dt * d(phi)/dt liefert O(dt^2) startnäherung statt O(dt)
  // TODO global/lokal trennen, lokal templaten
  // TODO testen ob mehrere fgiters funktionieren
  // TODO entweder enrichments manuell rein oder geeignetes non-enriched element wählen
  // TODO restart testen
  // TODO testen, ob fgiter>1 funzt
  
//  cout << "old phi on proc " << myrank_ << ": " << *phiOldRow_;
//  cout << "new phi on proc " << myrank_ << ": " << *phiNewRow_;
  
/*------------------------*
 * Initialization         *
 *------------------------*/
  
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element

//#ifdef DEBUG
  cout << " Computing new startdata for " << curr_.size() << " nodes..." << endl;
//#endif
  
  
/*---------------------------------------*
 * first part: get a sensible and good  *
 * start value for an interface changing *
 * node in a lagrangian point of view    *
 *---------------------------------------*/
//  cout << "before function call on proc " << myrank_ << endl;
//  for (size_t i=0;i<curr_.size();i++)
//    cout << "on proc " << myrank_ << " startpoint for "
//        << curr_[i].movNode_ << " is searched" << endl;
  
  startValuesStartValue(dt);
  
//  discret_->Comm().Barrier();
//  cout << "after calling startValuesStartValue - function on proc " <<
//      myrank_ << " size is " << curr_.size() << endl;
//  for (size_t i=0;i<curr_.size();i++)
//    cout << "on proc " << myrank_ << " startpoint for "
//        << curr_[i].movNode_ << " is " << curr_[i].startpoint_;

  
/*----------------------------------------------------*
 * second part: get the correct origin for the node    *
 * in a lagrangian point of view using a newton loop  *
 *----------------------------------------------------*/
  
#ifdef PARALLEL
  int counter = 0; // loop counter to avoid infinite loops
  bool procfinished; // true if movNodes is empty
  
  // loop over nodes which still dont have and may get a good startvalue
  while (true)
  {
    next_.clear();
//    cout << curr_.size() << " nodes on proc " << myrank_ << endl;
    counter += 1;
    
    // counter limit because maximal 10 newton iterations with maximal
    // numproc processor changes per iteration (avoids infinite loop)
    if (curr_.size()>0 && counter<10*numproc_)
    {
      procfinished = false; // processor has still nodes to examine
#endif
      
      // loop over all nodes which changed interface side
      for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
      {
        StartpointData curr = curr_[nodeid];
        
//        cout << "on proc " << myrank_ <<
//            " iteration starts for\n" << curr.movNode_ << endl;
        
        // if no startpoint found, interface is too far away from point
        // so that lagrangian point of view is senseless
        if (curr.changed_)
        {
          // Initialization
          DRT::Element* fittingele = NULL;  // pointer to the element where start point lies in
          LINALG::Matrix<nsd,1> xi(true);     // local transformed coordinates of x
          LINALG::Matrix<nsd,1> vel(true);    // velocity of the start point approximation
          double phin = 0;   // phi-value of the start point approximation
          bool elefound = false;             // true if an element for a point was found on the processor
          bool ihSide = false;             // true if point lies on correct side
          
          // search for an element where the current startpoint lies in
          // if found, give out all data at the startpoint
          elementAndLocalCoords(fittingele,curr.startpoint_,xi,vel,phin,elefound);
          
          // if element is not found, look at another processor
          // and so add all according data to the vectors
          // which will be sent to the next processor
          if (!elefound)
          {
            if (curr.searchedProcs_ < numproc_)
            {
              StartpointData next(
                  curr.movNode_,curr.startpoint_,curr.changed_,curr.phiSign_,
                  curr.searchedProcs_+1,curr.iter_,curr.startGid_,curr.startOwner_);
              next_.push_back(next);
            }
            else // all procs searched -> point not in domain
            {
              StartpointData failed(curr.movNode_,curr.startGid_,curr.startOwner_);
              failed_.push_back(failed);
              cout << "WARNING! Lagrangian start point not in domain!" << endl;
            }
          }
          
          // if element is found, the newton iteration 
          // to find a better start can start
          else
          {
            NewtonLoop(
                fittingele,curr.movNode_,dt,curr.startpoint_,xi,vel,
                phin,curr.phiSign_,elefound,ihSide,curr.iter_);
            
            // if iteration can go on (that is when startpoint is on
            // correct interface side and iter < 10)
            if (ihSide && curr.iter_<20) // TODO define maxiter here...
            {
              
              // if element is not found in a newton step, look at
              // another processor and so add all according data to
              // the vectors which will be sent to the next processor
              if (!elefound)
              {
                StartpointData next(
                    curr.movNode_,curr.startpoint_,curr.changed_,curr.phiSign_,
                    2,curr.iter_,curr.startGid_,curr.startOwner_);
                next_.push_back(next);
              }
              
              // newton iteration converged to a good startpoint
              // and so the data can be used to go on
              else
              {
                backTracking(
                    fittingele,curr.movNode_,curr.startpoint_,xi);
              }
            } // end if ihSide
            
            // if startpoint lies on wrong side at any time,
            // this algorithm has to stop for the searched node
            // which is done by setting changed = false
            if (!ihSide)
              curr_[nodeid].changed_ = false;
            
            if (curr.iter_==20) // maximum number of iterations reached//TODO sollte maxiter sein und übergeben werden
            {
              cout << "WARNING: start value set somehow sensible\n"
                      "  but newton iteration to find it didnt converge!" << endl;
              backTracking(
                  fittingele,curr.movNode_,curr.startpoint_,xi);
            }
          } // end if over nodes which changed interface
        } // end if startpoint changed
//        cout << "on proc " << myrank_ << " after " << curr.iter_ <<
//            " iterations the startpoint is " << curr.startpoint_ << endl;
      } // end loop over all nodes with changed interface side
    } // end if movenodes is empty or max. iter reached
    
    // maximum iteration number reached or
    // no startpoints left to look for on this proc
    else
      procfinished = true;
    
    // all nodes which got changed = false in this step
    // need another algorithm to be set
    for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
    {
      if (!curr_[nodeid].changed_)
      {
        StartpointData failed(
            curr_[nodeid].movNode_,
            curr_[nodeid].startGid_,
            curr_[nodeid].startOwner_
        );
        failed_.push_back(failed);
        cout << "WARNING! " << curr_[nodeid].movNode_ <<
            "\nhas no sensible startpoint in an Lagrangian point of view" << endl;
      }
    }
    
//    cout << "before calling iter-export function on proc " << myrank_ << endl;
//    discret_->Comm().Barrier();
    
#ifdef PARALLEL
    // export nodes and according data for which
    // the startpoint isnt still found to next proc
    exportIterData(procfinished);
    
//    cout << "after calling iter-export function on proc " << myrank_ << endl;
//    discret_->Comm().Barrier();
    
    // convergencecheck: procfinished == 1 just if all procs have finished
    // then the loop can stop because no more nodes to look for are left
    if (procfinished)
      break;
    
//    cout << "on proc " << myrank_ << " after newton loop " << counter << endl;
  } // end while loop over searched nodes
#endif
  
//  cout << "after newton interation on proc " << myrank_ << endl;
//  discret_->Comm().Barrier();
  
  
/*------------------------------------------------------*
 * third part: get sensible startvalues for nodes      *
 * where the algorithm failed, using another algorithm, *
 * and combine the "Done" and the "Failed" - vectors    *
 *------------------------------------------------------*/
  
  // nodes which are still in curr-vector, ran into an
  // infinite loop and so have to be set otherwise
  if (!procfinished)
  {
    for (size_t inode=0;inode<curr_.size();inode++)
    {
      StartpointData failed(
          curr_[inode].movNode_,
          curr_[inode].startGid_,
          curr_[inode].startOwner_);
      failed_.push_back(failed);
    }
  }
  
  curr_.clear(); // no more needed
  
//  cout << "after final filling of failed on proc " << myrank_ << endl;
//  
//  cout << "\n\nnow every node should be either in failed or in done vector!" << endl;
//  for (size_t i=0;i<done_.size();i++)
//    cout << "on proc " << myrank_ << " done node is " <<
//        done_[i].movNode_ << "\nwith velocity " << done_[i].velValues_[3]
//        << " and pressure " << done_[i].presValues_[3] << endl;
//  for (size_t i=0;i<failed_.size();i++)
//    cout << "on proc " << myrank_ << " failed node is "
//        << failed_[i].movNode_ << "\nwith startGid " << failed_[i].startGid_ <<
//         " and startOwner " << failed_[i].startOwner_ << endl;
  
//  // TODO remove if rest of algorithm works again...
//  failed_ = curr_;
  
  exportAlternativAlgoData();
  getDataForNotConvergedNodes();
  
//  cout << "after handling failed nodes on proc " << myrank_ << endl;
//  
//  failed_.clear(); // no more needed
//  cout << "\n\nnow every node should be in done vector!" << endl;
//  for (size_t i=0;i<done_.size();i++)
//    cout << "on proc " << myrank_ << " done node is " << done_[i].movNode_ <<
//        "\nwith velocity " << done_[i].velValues_[3] << " and pressure " <<
//        done_[i].presValues_[3] << endl;
//  discret_->Comm().Barrier();
  
  
  
/*-----------------------------------------------------------*
 * fifth part: set the computed values into the state vector *
 *-----------------------------------------------------------*/
  
#ifdef PARALLEL
  // send the computed startvalues for every node which needs
  // new start data to the processor where the node is
//  cout << "before final exporting on proc " << myrank_ << endl;
//  discret_->Comm().Barrier();
  exportFinalData();
//  cout << "after final exporting on proc " << myrank_ << endl;
//  discret_->Comm().Barrier();
#endif
  
//  cout << "\n\nnow every node should be in done vector on its own processor!" << endl;
//  for (size_t i=0;i<done_.size();i++)
//    cout << "on proc " << myrank_ << " done node is " << done_[i].movNode_ <<
//        "\nwith velocity " << done_[i].velValues_[3] << " and pressure " <<
//        done_[i].presValues_[3] << endl;
  
  // now every proc has the whole data for the nodes
  // and so the data can be set to the right place now
  setFinalStartData();
//  cout << "setting done on proc " << myrank_ << endl;
#ifdef DEBUG
  // TODO @Martin this does not compile in the DEBUG version "numproc not defined n this range"
  //if (counter > 8*numproc) // too much loops shouldnt be if all this works
  //  cout << "WARNING: semiLagrangeExtrapolation seems to run an infinite loop!" << endl;
#endif
  
  curr_ = init_; // for next vector TODO to be changed when all is done with one vector of epetra vectors
  
} // end semiLagrangeExtrapolation
//TODO test if seriell version works



/*------------------------------------------------------------------------------------------------*
 * finding startvalues for the new semi-lagrangian startvalues                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::startValuesStartValue(
    const double& dt
)
{
  //Initialization
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
//  vector<DRT::Element*> elevector_; // elements intersected by old interface at time n
//  
//  // fill element vector with intersected elements
//  for (int leleid=0;leleid<discret_->NumMyRowElements();leleid++)
//  {
//    DRT::Element* ele = discret_->lRowElement(leleid);
//    if (ih_old_->ElementIntersected(ele->Id()))
//      elevector_.push_back(ele);
//  }
  
  // loop over processors
  for (int procid=0; procid<numproc_; procid++)
  {
    // loop over nodes which changed interface side
    for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
    {
      StartpointData curr = curr_[nodeid];
      bool currChanged = false; // startvalue for startvalue changed on this proc?
      LINALG::Matrix<nsd,1> newNodeCoords(curr.movNode_.X()); // coords of endpoint
      
      // loop over intersected elements on processor
      for (std::set<int>::const_iterator enrnode = oldEnrNodes_.begin();
          enrnode != oldEnrNodes_.end();
          enrnode++)
      {
//        const int* elenodeids = elevector_[interele]->NodeIds();  // nodegids of element
//        
//        // loop over nodes of intersected element
//        for (int elenode=0;elenode<elevector_[interele]->NumNode();elenode++)
//        {
          DRT::Node* lnodeold = discret_->gNode(*enrnode);//elenodeids[elenode]);  // node near interface
          LINALG::Matrix<nsd,1> oldNodeCoords(lnodeold->X());  // coords of potential startpoint
          
          // just look for points on the same interface side
          if (curr.phiSign_ == (((*phiOldCol_)[lnodeold->LID()]>=0) ? 1 : -1))
          {
            LINALG::Matrix<nsd,1> diff;  // vector from old point at time n to new point at time n+1
            diff.Update(1.0,newNodeCoords,-1.0,oldNodeCoords);
            
            if (diff.Norm2()<curr.dMin_)
            {
              // get nodevelocity of a node
              LINALG::Matrix<nsd,1> nodevel(true);  // velocity of "old" node at time n
              const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(*enrnode));//elenodeids[elenode]));
              for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
                  fieldenr != fieldenrset.end();++fieldenr)
              {
                const DofKey<onNode> olddofkey(*enrnode,*fieldenr);//elenodeids[elenode], *fieldenr);
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
              
              // set new data
              curr.startpoint_.Update(1.0,newNodeCoords,-dt,nodevel);
              curr.dMin_ = diff.Norm2();
              curr.startGid_ = *enrnode;//elenodeids[elenode];
              curr.startOwner_ = discret_->gNode(*enrnode)->Owner();//elenodeids[elenode])->Owner();
              currChanged = true;
            } // end if new best startvalue
          } // end if just points on the correct side
//        } // end loop over nodes of intersected element
      } // end loop over intersected elements
      if (currChanged == true)
      {
        curr.changed_ = true;
        curr_[nodeid] = curr;
      }
    } // end loop over nodes which changed interface side
    
#ifdef PARALLEL
//    discret_->Comm().Barrier();
//    cout << "before exporting on proc " << myrank_ << endl;
//    discret_->Comm().Barrier();
//    for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
//      cout << "on proc " << myrank_ << curr_[nodeid].movNode_ <<
//          " got startvalue " << curr_[nodeid].startpoint_ << endl;
    
    exportStartData();
//    discret_->Comm().Barrier();
//    cout << "after exporting on proc " << myrank_ << endl;
//    for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
//      cout << "on proc " << myrank_ << curr_[nodeid].movNode_ <<
//          " has startvalue " << curr_[nodeid].startpoint_ << endl;
//    discret_->Comm().Barrier();
#endif
    
  } // end loop over processors
  
  // test loop over all nodes which changed interface side
  for (size_t nodeid=0; nodeid<curr_.size(); nodeid++)
  {
    if (!curr_[nodeid].changed_)
    {
      cout << "WARNING! No point on one interface side found!\nThis indicates "
          "that the whole area is at one side of the interface!" << endl;
      break;
    }
  } // end loop over nodes
//  cout << "on proc " << myrank_ << " size is " << curr_.size() << endl;
} // end startValuesStartValueFinder
// TODO if it turns out that this search for a startnode is computational
// expensive (i don't think so) then something else like a search tree has
// to be build and could be used then


/*------------------------------------------------------------------------------------------------*
 * setting data where semi-lagrangian approach failed                            winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::getDataForNotConvergedNodes(
)
{
  //Initialization
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  
  // now data is on proc where the startnode is -> final data can be computed here
  for (size_t inode=0;inode<failed_.size();inode++)
  {
    StartpointData failed = failed_[inode];
    
    DRT::Node* node = discret_->gNode(failed.startGid_);
    DRT::Element* nodesElement = node->Elements()[0];
    
    LINALG::Matrix<nsd,1> x(node->X());
    LINALG::Matrix<nsd,1> xi(true);
    LINALG::Matrix<nsd,1> vel(true);
    double phi = 0;
    bool elefound = false;
    
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
      backTracking(nodesElement,failed.movNode_,x,xi);
    } // end if elefound
  } // end loop over nodes
} // end getDataForNotConvergedNodes



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
  DRT::Element* currele;                 // pointer to the currently viewed element
  LINALG::Matrix<nsd,numnode> nodecoords;  // node coordinates of the element
  
  LINALG::Matrix<nsd,1> incr;        // increment of newton iteration
  LINALG::Matrix<nsd,1> residuum;    //residuum of newton iteration
  
  LINALG::Matrix<numnode,1> shapeFcn;               // shape functions
  LINALG::Matrix<nsd,numnode> shapeXiDeriv1;        // derivation of shape functions
  
  LINALG::Matrix<nsd,nsd> xjm;        //jacobian at start vector
  LINALG::Matrix<nsd,nsd> xji;        // invers of jacobian
  
  int iter;   // iteration counter
  
  const int max_iter = 10;   // maximum number of iterations
  const double relTolRes = 1e-010;   // relative tolerance for residual
  const double relTolIncr = 1e-010;   // relaltive tolerance for increment
  const double XiTol = 1e-6; // max. distance from point to cell in xi-coordinates
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
    
    // get the element coordinates for the point x with help of a
    // nonlinear solution in 3D element coordinates via Newton iteration
    
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
    if (xi(0)<=1+XiTol && xi(1)<=1+XiTol && xi(2)<=1+XiTol
        && xi(0)>=-1-XiTol && xi(1)>=-1-XiTol && xi(2)>=-1-XiTol)
    {
      fittingele = currele;
      elefound = true;
      break;
    } // end if near the element
  } // end loop over processor elements
  
  if (fittingele != NULL) // fittingele found
  {
    const int* elenodeids = fittingele->NodeIds();  // nodegids of element nodes
//    cout << "element is " << *fittingele << endl;
    
    // shape function in local coordinates
    DRT::UTILS::shape_function_3D(shapeFcn,xi(0),xi(1),xi(2),DISTYPE);
    LINALG::Matrix<nsd,numnode> shapeXiDeriv1;
    LINALG::Matrix<2*nsd,numnode> shapeXiDeriv2;
    
    // nodal phivalues
    LINALG::Matrix<numnode,1> nodephi;
    for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
      nodephi(nodeid,0) = (*phiOldCol_)[discret_->gNode(elenodeids[nodeid])->LID()];
    
    // get phivalue of point
    phi = nodephi.Dot(shapeFcn);
//    cout << "node phivalues " << nodephi << "interpolated phivalue " << phi;
    
    // initialize nodal vectors
    LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
    
    // enriched shape functions in local coordinates (N * \Psi)
    LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
    
    // create an element dof manager
    Teuchos::RCP<XFEM::ElementDofManager> eleDofManager;
    
    // create an empty element ansatz map
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty;
    eleDofManager = rcp(new XFEM::ElementDofManager(*fittingele,element_ansatz_empty,*olddofman_));
    
    // the enrichment functions may depend on the point
    // therefore the computation of the enrichment functions is called here
    // the gauss point is contained in shapeFcn!
    const XFEM::ElementEnrichmentValues enrvals(
        *fittingele,*eleDofManager,nodephi,
        shapeFcn,shapeXiDeriv1,shapeXiDeriv2);
    
    // enriched shape functions and derivatives for nodal parameters (dofs)
    enrvals.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx,shapeFcn,enrShapeFcnVel);
//      cout << "complete shapefunction for pressure is " << enrShapeFcnPres;
    
    int dofcounterVelx = 0;
    int dofcounterVely = 0;
    int dofcounterVelz = 0;
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
        case XFEM::Enrichment::typeVoidFSI : // TODO changes anything if void enrichment is used?
        case XFEM::Enrichment::typeVoid :
        case XFEM::Enrichment::typeKink :
        {
          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
          {
            nodevel(0,dofcounterVelx) = (*veln_)[olddofcolmap_.LID(olddofpos)];
            dofcounterVelx++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
          {
            nodevel(1,dofcounterVely) = (*veln_)[olddofcolmap_.LID(olddofpos)];
            dofcounterVely++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
          {
            nodevel(2,dofcounterVelz) = (*veln_)[olddofcolmap_.LID(olddofpos)];
            dofcounterVelz++;
          }
          else if (fieldenr->getField() == XFEM::PHYSICS::Pres); // not needed here
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
    
    // interpolate velocity and pressure values at starting point
    vel.Multiply(nodevel, enrShapeFcnVel);
    
//    cout << "element in which startpoint lies in is:\n" << *fittingele;
//    cout << "nodal velocities  are " << nodevel;
//    cout << "number of enriched nodes is " << dofcounterPres-numnode << endl;
//    cout << "startpoint approximation: " << x;
//    cout << "in local xi-coordinates: " << xi;
//    cout << "nodal shape function values are " << shapeFcn;
//    cout << "nodal enr function values for vel in local coords are " << enrShapeFcnVel;
//    cout << "interpolated velocity at this startpoint is " << vel;
  } // end if fitting ele found
} // end function findElementAndLocalCoords
// TODO if it turns out that this search for a element is computational
// too expensive (i don't think so) then something else like a search
// tree has to be build and could be used then


/*------------------------------------------------------------------------------------------------*
 * main Newton loop for semi-lagrangian algorithm                                winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::NewtonLoop(
    DRT::Element*& fittingele,
    DRT::Node& xNode,
    const double& dt,
    LINALG::Matrix<3,1>& xAppr,
    LINALG::Matrix<3,1>& xi,
    LINALG::Matrix<3,1>& vel,
    double& phi,
    int& phiSign,
    bool& elefound,
    bool& ihSide,
    int& iter
)
{
  
  const size_t numnode = 8;  // 8 nodes for hex8
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
    
  //cout << "Startwert\n" << xAppr << endl;
  //cout << "phiOrigSign " << phiSign << endl;
  //cout << "phiApprValue " << phi(0,0) << " und phiApprSign" << ((phi(0,0)>=0) ? 1 : -1) << endl;
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8; // only hex8 implemented
  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");
  
  if (phiSign != ((phi>=0) ? 1 : -1))
  {
    
#ifdef DEBUG
    cout << "phiOrigSign " << phiSign << endl;
    cout << "phiApprValue " << phi << " und phiApprSign" << ((phi>=0) ? 1 : -1) << endl;
    cout << "phivalues have different sign\nlagragian origin lies on other side as target" << endl;
#endif
    
    ihSide = false;
  }
  else  // Newton loop just for sensible points
  {
    ihSide = true;
    
    // Initialization
    LINALG::Matrix<numnode,1> shapeFcn;          // shape functions
    LINALG::Matrix<nsd,numnode> shapeXiDeriv1;        // derivation of shape functions
    
    LINALG::Matrix<nsd,nsd> xjm;        //jacobian at start vector
    LINALG::Matrix<nsd,nsd> xji;        // invers of jacobian
    
    LINALG::Matrix<nsd,numnode> nodevels;    // matrix with node-velocities for the startpoints element
    LINALG::Matrix<nsd,1> residuum;          // residuum of the newton iteration
    LINALG::Matrix<nsd,nsd> systemMatrix;      // matrix for the newton system
    LINALG::Matrix<nsd,1> incr;              // increment of the newton system
    
    const int max_iter = 10; // maximum number of iterations
    const double relTolIncr = 1.0e-10;;  // tolerance for the increment
    const double relTolRes = 1.0e-10;    // tolerance for the residual
    
    LINALG::Matrix<nsd,numnode> nodecoords; // node coordinates of the element
    for (int elenodeid=0;elenodeid<fittingele->NumNode();elenodeid++)
    {
      DRT::Node* currnode = discret_->gNode(fittingele->NodeIds()[elenodeid]);
      for (size_t i=0;i<nsd;i++)
        nodecoords(i,elenodeid) = currnode->X()[i];
    }
    
    LINALG::Matrix<nsd,1> origNodeCoords;
    for (size_t i=0;i<nsd;i++)
      origNodeCoords(i) = xNode.X()[i];
    
    // evaluate shape functions at xi
    DRT::UTILS::shape_function_3D(shapeFcn, xi(0),xi(1),xi(2),DISTYPE);
    // evaluate derivative of shape functions at xi
    DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1, xi(0),xi(1),xi(2),DISTYPE);
    
    xjm.MultiplyNT(shapeXiDeriv1,nodecoords);   // jacobian J = (dx/dxi)^T
    xji.Invert(xjm);       // jacobian inverted J^(-1) = dxi/dx
    
    LINALG::Matrix<nsd,numnode> shapeXYDeriv1; // first derivation of global shape functions
    shapeXYDeriv1.Multiply(xji,shapeXiDeriv1); // (dN/dx)^T = (dN/dxi)^T * J^(-T)
//    cout << "1) deriv1 of global shape fcn is " << shapeXYDeriv1 << endl;

    // initialize residual
    residuum.Clear();
    residuum.Update(dt, vel); // dt*v(xAppr)
    residuum.Update(1.0,xAppr,-1.0,origNodeCoords,1.0);  // R = xAppr - xNode + dt*v(xAppr)
    
    while(iter < max_iter)  // newton loop
    {
      iter += 1;
      
      { // build Matrix
        const int* elenodeids = fittingele->NodeIds();  // nodeids of element
        systemMatrix.Clear();
        
        // initialize nodal vectors
        LINALG::Matrix<nsd,2*numnode> nodevel(true); // node velocities of the element nodes
        
        // second derivative of shape function in local and global coordinates
        LINALG::Matrix<2*nsd,numnode> shapeXiDeriv2; // just needed for function call
        LINALG::Matrix<2*nsd,numnode> shapeXYDeriv2(true);
        
        // enriched shape functions and there derivatives in local coordinates (N * \Psi)
        LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
        LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
        LINALG::Matrix<2*nsd,2*numnode> enrShapeXYDeriv2(true); // just needed for function call so just one for vel and pres
        
        // nodal phivalues
        LINALG::Matrix<numnode,1> nodephi;
        for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
          nodephi(nodeid,0) = (*phiOldCol_)[discret_->gNode(elenodeids[nodeid])->LID()];
          
        // create an element dof manager
        Teuchos::RCP<XFEM::ElementDofManager> eleDofManager;
        const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty; // ansatz map needed for eledofman
        eleDofManager = rcp(new XFEM::ElementDofManager(*fittingele,element_ansatz_empty,*olddofman_));
        
        // the enrichment functions may depend on the point
        // therefore the computation of the enrichment functions is called here
        // the gauss point is contained in shapeFcn!
        const XFEM::ElementEnrichmentValues enrvals(
            *fittingele,*eleDofManager,nodephi,
            shapeFcn,shapeXYDeriv1,shapeXYDeriv2);
//        cout << "2) deriv1 of global shape fcn is " << shapeXYDeriv1 << endl;
          
        // enriched shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,shapeFcn,shapeXYDeriv1,shapeXYDeriv2,
            enrShapeFcnVel,enrShapeXYVelDeriv1,enrShapeXYDeriv2); // enrichment assumed to be equal for the 3 dimensions
//        cout << "enr shapefunction for vel is " << enrShapeFcnVel << endl;
//        cout << "enr shapefunction for deriv of vel is " << enrShapeXYVelDeriv1 << endl;
          
        int dofcounterVelx = 0;
        int dofcounterVely = 0;
        int dofcounterVelz = 0;
        for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
        {
          // get nodal velocities and pressures with help of the field set of node
          const std::set<XFEM::FieldEnr>& fieldEnrSet(dofman_->getNodeDofSet(elenodeids[nodeid]));
          for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
              fieldenr != fieldEnrSet.end();++fieldenr)
          {
            const DofKey<onNode> olddofkey(elenodeids[nodeid], *fieldenr);
            const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
            switch (fieldenr->getEnrichment().Type())
            {
            case XFEM::Enrichment::typeStandard :
            case XFEM::Enrichment::typeJump :
            case XFEM::Enrichment::typeVoidFSI : // TODO changes anything if void enrichment is used?
            case XFEM::Enrichment::typeVoid :
            case XFEM::Enrichment::typeKink :
            {
              if (fieldenr->getField() == XFEM::PHYSICS::Velx)
              {
                nodevel(0,dofcounterVelx) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                dofcounterVelx++;
              }
              else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
              {
                nodevel(1,dofcounterVely) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                dofcounterVely++;
              }
              else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
              {
                nodevel(2,dofcounterVelz) = (*veln_)[olddofcolmap_.LID(olddofpos)];
                dofcounterVelz++;
              }
              else if (fieldenr->getField() == XFEM::PHYSICS::Pres); // not needed here
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
//        cout << "nodal velocities of this element are " << nodevel;
        
        systemMatrix.MultiplyNT(dt,nodevel,enrShapeXYVelDeriv1); // v_nodes * dN/dx
        for (size_t i=0;i<nsd;i++)
          systemMatrix(i,i) += 1.0; // I + dt*v_nodes*dN/dx
        systemMatrix.Invert();
      } // system Matrix built
      
      //solve Newton iteration
      incr.Clear();
      incr.Multiply(-1.0,systemMatrix,residuum); // incr = -Systemmatrix^-1 * residuum
      
      // update iteration
      for (size_t i=0;i<nsd;i++)
        xAppr(i) += incr(i);
      
      //=============== update residuum================
      elementAndLocalCoords(fittingele, xAppr, xi, vel, phi, elefound);
      
      if (elefound) // element of xAppr at this processor
      {
//        cout << " in ele " << *fittingele << " \nwith xAppr "
//            << xAppr << "and xi " << xi;
        if (phiSign != ((phi>=0) ? 1 : -1))
        {
          
#ifdef DEBUG
          cout << "phiOrigSign " << phiSign << endl;
          cout << "phiApprValue " << phi << " und phiApprSign" << ((phi>=0) ? 1 : -1) << endl;
          cout << "phivalues have different sign\nlagragian origin lies on other side as target" << endl;
#endif
          
          ihSide = false;
          break; // leave newton loop if element is on wrong domain side
        }
        else
        {
          ihSide = true;
          
          // reset residual
          residuum.Clear();
          residuum.Update(dt, vel); // dt*v(xAppr)
          residuum.Update(1.0,xAppr,-1.0,origNodeCoords,+1.0);  // R = xAppr - xNode + dt*v(xAppr)
//          cout << " and tols are " << incr.Norm2()/xAppr.Norm2()
//              << " and " << residuum.Norm2()/xAppr.Norm2() << endl;
          
          // convergence criterion
          if (xAppr.Norm2()>1e-3)
          {
            if (incr.Norm2()/xAppr.Norm2() < relTolIncr && residuum.Norm2()/xAppr.Norm2() < relTolRes)
              break;
          }
          else
          {
            if (incr.Norm2() < relTolIncr && residuum.Norm2() < relTolRes)
              break;
          }
          
          // update jacobian and inverted jacobian
          xjm.MultiplyNT(shapeXiDeriv1,nodecoords);
          xji.Invert(xjm);
        }
      }
      else // element of xAppr not at this processor
        break; // stop newton loop on this proc
    }
    
#ifdef DEBUG
    // did newton iteration converge?
    if(iter == max_iter){cout << "WARNING: newton iteration for finding start value not converged for point\n" << endl;}
//    cout << "after " << iter << " iterations the endpoint is\n" << xAppr << endl;
#endif
    
  } // end if point in correct domain
}



/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::backTracking(
    DRT::Element*& fittingele,
    DRT::Node& node,
    LINALG::Matrix<3,1>& x,
    LINALG::Matrix<3,1>& xi
)
{
  const DRT::Element::DiscretizationType DISTYPE = DRT::Element::hex8; // only hex8 implemented
  if (DISTYPE != DRT::Element::hex8)
    dserror("element type not implemented until now!");
  
//  cout << node << "has lagrange origin " << x << "\nwith xi-coordinates "
//      << xi << "\n in element " << *fittingele << endl;
  const int numnode = 8; // 8 nodes for hex8
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  const int* elenodeids = fittingele->NodeIds();  // nodegids of element nodes
  
  LINALG::Matrix<numnode,1> shapeFcn; // shape function
  DRT::UTILS::shape_function_3D(shapeFcn,xi(0),xi(1),xi(2),DISTYPE);
  // first derivative of shape function in local coordinates
  LINALG::Matrix<nsd,numnode> shapeXiDeriv1;
  DRT::UTILS::shape_function_3D_deriv1(shapeXiDeriv1,xi(0),xi(1),xi(2),DISTYPE);
//  cout << "deriv1 of local shape fcn is " << shapeXiDeriv1 << endl;
  
  LINALG::Matrix<nsd,nsd> xji; // invers of jacobian
  { // Computing the invers jacobian
    LINALG::Matrix<nsd,nsd> xjm; //jacobian dx/dxi at start vector
    LINALG::Matrix<nsd,numnode> nodecoords; // node coordinates of the element
    
    // evaluate nodecoords
    for (int elenodeid=0;elenodeid<fittingele->NumNode();elenodeid++)
    {
      DRT::Node* currnode = discret_->gNode(fittingele->NodeIds()[elenodeid]);
      for (size_t i=0;i<nsd;i++)
        nodecoords(i,elenodeid) = currnode->X()[i];
    }
    
    xjm.MultiplyNT(shapeXiDeriv1,nodecoords);   // jacobian J = (dx/dxi)^T
    xji.Invert(xjm);       // jacobian inverted J^(-1) = dxi/dx
  }
  
  LINALG::Matrix<nsd,numnode> shapeXYDeriv1; // first derivation of global shape functions
  shapeXYDeriv1.Multiply(xji,shapeXiDeriv1); // (dN/dx)^T = (dN/dxi)^T * J^(-T)
  
  // second derivative of shape function in local and global coordinates
  LINALG::Matrix<2*nsd,numnode> shapeXiDeriv2; // just needed for function call
  LINALG::Matrix<2*nsd,numnode> shapeXYDeriv2(true);
  
  // enriched shape functions and there derivatives in local coordinates (N * \Psi)
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);
  LINALG::Matrix<2*nsd,2*numnode> enrShapeXYDeriv2(true); // just needed for function call so just one for vel and pres
  
  // nodal phivalues
  LINALG::Matrix<numnode,1> nodephi;
  for (int nodeid=0;nodeid<fittingele->NumNode();nodeid++) // loop over element nodes
    nodephi(nodeid,0) = (*phiOldCol_)[discret_->gNode(elenodeids[nodeid])->LID()];
  
  // create an element dof manager
  Teuchos::RCP<XFEM::ElementDofManager> eleDofManager;
  const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz_empty; // element ansatz map needed for eledofman
  eleDofManager = rcp(new XFEM::ElementDofManager(*fittingele,element_ansatz_empty,*olddofman_));
  
  // the enrichment functions may depend on the gauss point
  // therefore the computation of the enrichment functions is called here
  // the gauss point is contained in shapeFcn!
  const XFEM::ElementEnrichmentValues enrvals(
      *fittingele,*eleDofManager,nodephi,
      shapeFcn,shapeXYDeriv1,shapeXYDeriv2);
//    cout << "deriv1 of global shape fcn is " << shapeXYDeriv1 << endl;

  // enriched shape functions and derivatives for nodal parameters (dofs)
  enrvals.ComputeModifiedEnrichedNodalShapefunction(
      XFEM::PHYSICS::Velx,shapeFcn,shapeXYDeriv1,shapeXYDeriv2,
      enrShapeFcnVel,enrShapeXYVelDeriv1,enrShapeXYDeriv2); // enrichment assumed to be equal for the 3 dimensions
  enrvals.ComputeModifiedEnrichedNodalShapefunction(
      XFEM::PHYSICS::Pres,shapeFcn,shapeXYDeriv1,shapeXYDeriv2,
      enrShapeFcnPres,enrShapeXYPresDeriv1,enrShapeXYDeriv2);
//    cout << "enr shapefunction for pressure is " << enrShapeFcnPres;
//    cout << "enr shapefunction for deriv of pres is " << enrShapeXYPresDeriv1 << endl;
  
  
  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd,2*numnode> nodevel(true);
  // node velocities of the element nodes for the data that should be changed
  vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size());
  // node pressures of the element nodes for the data that should be changed
  vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size());
  
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
      case XFEM::Enrichment::typeVoidFSI : // TODO changes anything if void enrichment is used?
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
  
  LINALG::Matrix<3,1> transportVel; // transport velocity
  transportVel.Multiply(nodevel, enrShapeFcnVel);
  
  // computing pseudo time-step deltaT
  // remark: if x is the Lagrange-origin of node, deltaT = dt with respect to small errors.
  // if its not, deltaT estimates the time x needs to move to node)
  double deltaT; // pseudo time-step
  {
    LINALG::Matrix<nsd,1> diff(node.X());
    diff -= x; // diff = x_Node - x_Appr
    
    double numerator = transportVel.Dot(diff); // numerator = v^T*(x_Node-x_Appr)
    double denominator = transportVel.Dot(transportVel); // denominator = v^T*v
    
    if (denominator<1e-15)
      deltaT = 0;
    else
      deltaT = numerator/denominator;
  }
  
  vector<LINALG::Matrix<3,1> > velValues(oldVectors_.size()); // velocity of the data that should be changed
  vector<double> presValues(oldVectors_.size()); // pressures of the data that should be changed
  
  // interpolate velocity and pressure gradients for all fields at starting point and get final values
  for (size_t index=0;index<oldVectors_.size();index++)
  {
    LINALG::Matrix<nsd,1> vel; // velocity data
    vel.Multiply(nodeveldata[index],enrShapeFcnVel);
    LINALG::Matrix<nsd,nsd> velDeriv1; // first derivation of velocity data
    velDeriv1.MultiplyNT(nodeveldata[index],enrShapeXYVelDeriv1);
    vel.Multiply(deltaT,velDeriv1,transportVel,1.0); // vel = dt*velDeriv1*transportVel + vel
    velValues[index]=vel;
    
    LINALG::Matrix<1,1> pres; // pressure data
    pres.Multiply(nodepresdata[index],enrShapeFcnPres);
//    if (index == 3)
//      cout << "for " << node << " pres was " << pres;
    LINALG::Matrix<1,nsd> presDeriv1; // first derivation of pressure data
    presDeriv1.MultiplyNT(nodepresdata[index],enrShapeXYPresDeriv1);
    pres.Multiply(deltaT,presDeriv1,transportVel,1.0); // pres = dt*presDeriv1*transportVel + pres
    presValues[index] = pres(0); 
//    if (index == 3)
//      cout << " and became " << pres << endl;
  }
  
//  cout << "first derivation of velocity is " << velDeriv1;
//  cout << "first derivation of pressure is " << presDeriv1;
//  cout << "nodal velocities are " << nodeveldata[3];
//  cout << "nodal pressures are " << nodepresdata[3];
  
  StartpointData done(node,velValues,presValues);
  done_.push_back(done);
//              cout << "backtracking was done on proc " << myrank_ <<
//                  " with final pres " << pres;
} // end getFinalStartvalues



/*------------------------------------------------------------------------------------------------*
 * setting the final data in Epetra Vector for a node                            winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::setFinalStartData(
) const
{
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  
  for (size_t inode=0;inode<done_.size();inode++)
  {
    vector<LINALG::Matrix<nsd,1> > velValues(done_[inode].velValues_); // velocities of the node
    vector<double> presValues(done_[inode].presValues_); // pressures of the node
    
    const int gnodeid = done_[inode].movNode_.Id(); // global node id
    
    // set nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gnodeid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey<onNode> newdofkey(gnodeid, *fieldenr);
      const int newdofpos = nodalDofRowDistrib_.find(newdofkey)->second;
      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
//          cout << (*velnp_)[newdofrowmap_.LID(newdofpos)] << " becomes " << nodeVel(0) << endl;
          for (size_t index=0;index<newVectors_.size();index++)
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = velValues[index](0);
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
//          cout << (*velnp_)[newdofrowmap_.LID(newdofpos)] << " becomes " << nodeVel(1) << endl;
          for (size_t index=0;index<newVectors_.size();index++)
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = velValues[index](1);
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
//          cout << (*velnp_)[newdofrowmap_.LID(newdofpos)] << " becomes " << nodeVel(2) << endl;
          for (size_t index=0;index<newVectors_.size();index++)
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = velValues[index](2);
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
//          cout << (*velnp_)[newdofrowmap_.LID(newdofpos)] << " becomes " << done_[inode].pres_ << endl;
          for (size_t index=0;index<newVectors_.size();index++)
            (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] = presValues[index];
        }
      }
    } // end loop over fieldenr
//    cout << *discret_->gNode(gnodeid) << "\n gets velocity\n"
//    << nodeVel << "and pressure " << pres[inode] << endl;
  } // end loop over nodes
} // end setFinalStartData



void XFEM::Startvalues::setEnrichmentValues(
)//TODO statt DofGID lieber int...
{
  // Initialization
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
//  std::set<int> oldEnrNodes; // set with all gids of enriched nodes on proc
  std::set<int> newEnrNodes;
//  double halfJump = 0.0;
  // TODO jump is also unknown for some fields in 3D (general everywhere...)
  
//  cout << "sizes on proc " << myrank_ << " are: " << oldEnrNodalDofSet_.size() << ", "
//     << oldEnrNodalDofDistrib_.size() << " and " << allOldEnrValues[0].size() << endl;
  
  
/*---------------------------------------------------------------*
 * fill enr_ structure with data about nodes which changed      *
 * interface side and identify all enriched nodes for later use  *
 *---------------------------------------------------------------*/
  
  // 
  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
  {
    if (discret_->lColNode(lnodeid)->Owner() == myrank_)
    {
    const DRT::Node* lnode = discret_->lColNode(lnodeid);
    const int gid = lnode->Id();
    
    const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey<onNode> newdofkey(gid, *fieldenr);
//      const int newdofpos = nodalDofColDistrib_.find(newdofkey)->second;
      
      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
        map<DofKey<onNode>, DofGID>::const_iterator olddof = oldNodalDofColDistrib_.find(newdofkey);
        
        // olddof found and on same side -> enrichment existed before and can be used
        if ((olddof != oldNodalDofColDistrib_.end()) &&
            ((((*phiNewCol_)[lnodeid]>=0) ? 1 : -1) == (((*phiOldCol_)[lnodeid]>=0) ? 1 : -1)))
        {
//          if (fieldenr->getField() == XFEM::PHYSICS::Velx)
//            halfJump = 0.5;
//          else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
//            halfJump = 0.5;
//          
//          const double phiQuotient = (*phiNewRow_)[lnodeid]/(*phiOldRow_)[lnodeid]; // phi_n+1(x)/phi_n(x)
//          
//          for (size_t i=0;i<valueVector.size();i++)
//          {
//            cout << "normal enrichment value was " << (*valueVector[i])[newdofrowmap_.LID(newdofpos)];
//            const double currOldEnr = (*valueVector[i])[newdofrowmap_.LID(newdofpos)];
//            (*valueVector[i])[newdofrowmap_.LID(newdofpos)] = 
//                halfJump + (currOldEnr-halfJump)*phiQuotient;
//            cout << " and became " << (*valueVector[i])[newdofrowmap_.LID(newdofpos)] << endl;
//            // new_enr = jump_size/2 + (old_enr-jump_size/2)*phi_new/phi_old
//          }
//          halfJump = 0.0;
        // TODO uncomment if correct formula for new enrichment value is found
        } // end if dof found
        else // olddof not found or enrichment changed side -> new value needed
        {
          StartpointData enr(*lnode,(((*phiNewCol_)[lnodeid]>=0) ? 1 : -1),INFINITY,0.0,map<DofKey<onNode>,vector<double> >());
          enr_.push_back(enr);
          break;
          // all other enrichments for this node will be handled later
          // here the node should be added just one time, so loop should break
        } // end if dof not found
      } // end if jump enrichment
    } // end loop over fieldenr
    
    
//    const set<XFEM::FieldEnr>& oldfieldenrset(olddofman_->getNodeDofSet(gid));
//    for (set<XFEM::FieldEnr>::const_iterator oldfieldenr = oldfieldenrset.begin();
//        oldfieldenr != oldfieldenrset.end();
//        oldfieldenr++)
//    {
//      if (oldfieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
//      {
//        // not standard dof exists -> node enriched
//        oldEnrNodes.insert(gid); // TODO one set for every interface side would be better
//        break;
//      }
//    }
    
    
//    const set<XFEM::FieldEnr>& newfieldenrset(dofman_->getNodeDofSet(gid));
//    for (set<XFEM::FieldEnr>::const_iterator newfieldenr = newfieldenrset.begin();
//        newfieldenr != newfieldenrset.end();
//        newfieldenr++)
//    {
//      if (newfieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
//      {
//        // not standard dof exists -> node enriched
//        newEnrNodes.insert(gid); // TODO will be used later...
//        break;
//      }
//    }
    }
  } // end loop over processor nodes
  
//  discret_->Comm().Barrier();
//  cout << "number of points which need new enrichment values on proc " <<
//      myrank_ << " is " << enr_.size() << endl;
//  for (size_t i=0;i<enr_.size();i++)
//    cout << enr_[i].movNode_ << endl;
//  discret_->Comm().Barrier();
//  cout << "number of enriched nodes at new timestep on proc " << myrank_ <<
//      " is " << newEnrNodes.size() << endl;
//  for (set<int>::const_iterator i=newEnrNodes.begin();
//      i != newEnrNodes.end();i++)
//    cout << "node " << *i << " is enriched at new timestep" << endl;
//  discret_->Comm().Barrier();

/*--------------------------------------------------*
 * get the nearest node for the enriched nodes      *
 * which need a completely new enrichment value     *
 *--------------------------------------------------*/
  
  
//  cout << "\n\n\n" << endl;
//  for (map<int,set<XFEM::DofKey<onNode> > >::const_iterator i=oldEnrNodalDofSet_.begin();
//      i != oldEnrNodalDofSet_.end();i++)
//  {
//    cout << "node with gid " << i->first << " has dofkeys" << endl;
//    for(set<XFEM::DofKey<onNode> >::const_iterator j=i->second.begin();
//        j != i->second.end();j++)
//      cout << *j << endl;
//  }
//  cout << "\n\n\n" << endl;
//  for (map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID>::const_iterator i=oldEnrNodalDofDistrib_.begin();
//      i != oldEnrNodalDofDistrib_.end();
//      i++)
//    cout << i->first << ", " << i->second << endl;
//  cout << "\n\n\n" << endl;
//  for (map<int,double>::const_iterator i=allOldEnrValues[3].begin();
//      i != allOldEnrValues[3].end();i++)
//    cout << i->first << "   " << i->second << endl;
//  cout << "\n\n\n" << endl;
//  discret_->Comm().Barrier();
  
  // loop over processors
  for (int procid=0; procid<numproc_; procid++)
  {
    // loop over nodes with new enrichment
    for (size_t nodeid=0;nodeid<enr_.size();nodeid++)
    {
      StartpointData enr = enr_[nodeid];
      bool currChanged = false; // new nearest node on this proc?
      LINALG::Matrix<nsd,1> newNodeCoords(enr.movNode_.X()); // coords of endpoint
      const int movNodeGid = enr.movNode_.Id(); // global id of enriched node
      
//      cout << enr.movNode_ << " searches for a new nearest enriched node on proc "
//          << myrank_ << " and has data:\nphisign " << enr.phiSign_ << ", dist " <<
//          enr.dMin_ << " and phivalue " << enr.phiValue_<< endl;
      
      // loop over enriched nodes
      for (std::set<int>::const_iterator inode=oldEnrNodes_.begin();
          inode != oldEnrNodes_.end();
          inode++)
      {
         DRT::Node* lnodeold = discret_->gNode(*inode);  // enriched node
         
//         cout << *lnodeold << "is tested now!" << endl;
         
         // just look at enriched nodes on same interface side
        if (enr.phiSign_ == (((*phiOldCol_)[lnodeold->LID()]>=0) ? 1 : -1))
        { // TODO enrichments with phi = 0 can be used by both sides and
          // can be set by both sides because of the enrichments characteristics
          LINALG::Matrix<nsd,1> oldNodeCoords(lnodeold->X());  // coords of potential startpoint
          
          LINALG::Matrix<nsd,1> diff;  // vector from old point at time n to new point at time n+1
          diff.Update(1.0,newNodeCoords,-1.0,oldNodeCoords);
          
          if (diff.Norm2()<enr.dMin_)
          {
            enr.enrValues_.clear();
            
            // set new data
            currChanged = true;
            enr.dMin_ = diff.Norm2();
            enr.phiValue_ = (*phiOldCol_)[lnodeold->LID()];
            
            const set<FieldEnr> fieldenrset = olddofman_->getNodeDofSet(*inode);
            // loop over enrichments for this node
            // remark: every field should have same number of enrichments
            //         since they all should base on same map
            for (set<FieldEnr>::const_iterator fieldenr=fieldenrset.begin();
                fieldenr != fieldenrset.end();
                fieldenr++)
            {
              if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
              {
                const DofKey<onNode> olddofkey(*inode, *fieldenr);
                const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
                
                // change gid of dofkeys of enrichment values to the moving node gid
                // so it can be found later by the according value that has to be set
                const DofKey<onNode> newdofkey(movNodeGid,*fieldenr);
//                cout << "here with dof " << olddof->first << " and pos " << olddof->second << endl;
                
                vector<double> currEnrValues;
                for (size_t i=0;i<oldVectors_.size();i++)
                {
                  currEnrValues.push_back((*oldVectors_[i])[olddofcolmap_.LID(olddofpos)]);
//                  cout << "new size is " << currEnrValues.size() << endl;
                }
                enr.enrValues_.insert(make_pair(newdofkey,currEnrValues));
              } // end if dofkey found
            } // end loop over fieldenrset
//            cout << "nearest, at old time enriched node is " << *lnodeold << "with enrvalues\n";
//            for (map<DofKey<onNode>,vector<double> >::const_iterator i=enr.enrValues_.begin();
//                i!=enr.enrValues_.end();
//                i++)
//              cout << i->second[3] << "    ";
//            cout << endl;
          } // end if new nearest enriched point
        } // end loop over nodes of intersected element
      } // end loop over enriched nodes
      if (currChanged == true)
      {
        enr_[nodeid] = enr;
//        cout << "new vel enr values are:" << endl;
//        for (map<DofKey<onNode>,vector<double> >::const_iterator enrdata=enr.enrValues_.begin();
//            enrdata != enr.enrValues_.end();enrdata++)
//          cout << enrdata->second[3] << endl;
      }
    } // end loop over nodes with new enrichments
    
#ifdef PARALLEL
//    discret_->Comm().Barrier();
//    cout << "before exporting on proc " << myrank_ << endl;
//    for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
//      cout << "on proc " << myrank_ << curr_[nodeid].movNode_ <<
//          " got startvalue " << curr_[nodeid].startpoint_ << endl;
//    discret_->Comm().Barrier();
//    cout << "here before exporting on proc " << myrank_ << endl;
    
    exportEnrichmentData();
    
//    discret_->Comm().Barrier();
//    cout << "after exporting on proc " << myrank_ << endl;
//    for (size_t nodeid=0;nodeid<curr_.size();nodeid++)
//      cout << "on proc " << myrank_ << curr_[nodeid].movNode_ <<
//          " has startvalue " << curr_[nodeid].startpoint_ << endl;
//    cout << "here after exporting on proc " << myrank_ << endl;
//    discret_->Comm().Barrier();
#endif
    
  } // end loop over processors
  
//  discret_->Comm().Barrier();
//  cout << "every node is on its own processor now and has a nearest enriched node" << endl;
//  cout << "number of points with new enrichment values on proc " << myrank_ <<
//      " is " << enr_.size() << endl;
  
  // now every node is on its own processor again and has a nearest enriched node
  for (size_t inode=0;inode<enr_.size();inode++) // loop over newly enriched nodes
  {
    StartpointData enr = enr_[inode];
    const int gid = enr.movNode_.Id();
    
    const set<XFEM::FieldEnr>& fieldenrset(dofman_->getNodeDofSet(gid));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
        fieldenr != fieldenrset.end();++fieldenr)
    {
      const DofKey<onNode> newdofkey(gid, *fieldenr);
      const int newdofpos = nodalDofRowDistrib_.find(newdofkey)->second;
      
      // now nodes with new enrichments can be set
      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
//        const double phiQuotient = (*phiNewRow_)[newdofrowmap_.LID(gid)]/enr.phiValue_; // phi_n+1(x)/phi_n(x_near)
//        
//        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
//          halfJump = 0.5;
//        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
//          halfJump = 0.5;
        
        // find dofkey in enrichment values map
        std::map<DofKey<onNode>,vector<double> >::const_iterator currEnrValues = enr.enrValues_.find(newdofkey);
        if (currEnrValues != enr.enrValues_.end())
        {
          for (size_t i=0;i<newVectors_.size();i++)
          {
//            cout << "enrichment value was " << (*newVectors_[i])[newdofrowmap_.LID(newdofpos)];
            const double currOldEnr = currEnrValues->second[i];
            (*newVectors_[i])[newdofrowmap_.LID(newdofpos)] = currOldEnr;
//                halfJump + (currOldEnr-halfJump)*phiQuotient;
//            cout << " and became " << (*newVectors_[i])[newdofrowmap_.LID(newdofpos)] << endl;
            // new_enr = jump_size/2 + (old_enr-jump_size/2)*phi_new/phi_old
          }
        }
        else
        {
          cout << newdofkey << "wasn't found!" << endl;
//          for (std::map<DofKey<onNode>,vector<double> >::const_iterator i=enr.enrValues_.begin();
//              i != enr.enrValues_.end();i++)
//            cout << i->first << endl;
          dserror("bug! It seems that no enrichments at old timestep are available but they should be!");
        }
//        halfJump = 0.0;
      } // end if jump enrichment
    } // end loop over fieldenr
  } // end loop over newly enriched nodes
} // end function setEnrichmentValues



# ifdef PARALLEL
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
  
  // vector including all data that has to be send to the next proc
  vector<char> dataSend;
  
  // packing the data
  DRT::ParObject::AddtoPack(dataSend,curr_.size());
  for (size_t ipoint=0;ipoint<curr_.size();ipoint++)
  {
    StartpointData curr = curr_[ipoint];
    
    DRT::ParObject::AddtoPack(dataSend,curr.movNode_.Id());
    DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(curr.movNode_.X()));
    DRT::ParObject::AddtoPack(dataSend,curr.movNode_.Owner());
    DRT::ParObject::AddtoPack(dataSend,curr.startpoint_);
    DRT::ParObject::AddtoPack(dataSend,curr.changed_);
    DRT::ParObject::AddtoPack(dataSend,curr.phiSign_);
    DRT::ParObject::AddtoPack(dataSend,curr.searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend,curr.iter_);
    DRT::ParObject::AddtoPack(dataSend,curr.startGid_);
    DRT::ParObject::AddtoPack(dataSend,curr.startOwner_);
    DRT::ParObject::AddtoPack(dataSend,curr.dMin_);
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
  // ... and receive length
  vector<int> lengthRecv(1,0);
  exporter_.Receive(source, length_tag, lengthRecv, size_one);
  exporter_.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
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
  curr_.clear();
  
  // unpack received data
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
  for (size_t inode=0;inode<numberOfNodes;inode++)
  {
    int gid;
    LINALG::Matrix<nsd,1> coords;
    int owner;
    LINALG::Matrix<nsd,1> startpoint;
    int changed;
    int phiSign;
    int searchedProcs;
    int iter;
    int startGid;
    int startOwner;
    double dMin;
    
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,changed);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiSign);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,iter);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,dMin);
    
    double coordinates[nsd];
    for (size_t dim=0;dim<nsd;dim++)
      coordinates[dim] = coords(dim);
    
    DRT::Node movNode(gid,coordinates,owner);
    
    StartpointData curr(
        movNode,startpoint,changed,phiSign,
        searchedProcs,iter,startGid,startOwner,dMin);
    curr_.push_back(curr);
  }
  
  // check if all sizes fit
//  if (gids.size()!=numberOfNodes ||
//      coords.size()!=numberOfNodes ||
//      owners.size()!=numberOfNodes ||
//      curr_.startpoints_.size()!=numberOfNodes ||
//      dMin.size()!=numberOfNodes ||
//      curr_.changed_.size()!=numberOfNodes ||
//      curr_.phiSign_.size()!=numberOfNodes ||
//      curr_.startGid_.size()!=numberOfNodes ||
//      curr_.startOwner_.size()!=numberOfNodes)
//    dserror( "sending of starting point data failed!");
} // end exportStartData



void XFEM::Startvalues::exportAlternativAlgoData()
{
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  
  // array of vectors which stores data for
  // every processor in one vector
  vector<vector<StartpointData> > failedVec(numproc_);
  
  // fill vectors with the data
  for (size_t inode=0;inode<failed_.size();inode++)
  {
    failedVec[failed_[inode].startOwner_].push_back(failed_[inode]);
  }
  
  // clear final vectors
  failed_.clear();
  
  // set values for own proc
  failed_ = failedVec[myrank_];
  
  // clear the set data from the vector
  failedVec[myrank_].clear();
  
  
  // send data to the processor where the point lies
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    vector<char> dataSend; // vector including all data that has to be send to dest proc
    vector<int> lengthSend(1,0);
    int size_one = 1;
    
    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;
    
    dataSend.clear();
//    cout << "before sending on proc " << myrank_ << " sizes are "
//        << gid[dest].size() << ", " << movGid[dest].size() <<
//        ", " << movCoords[dest].size() << " and " <<
//        movOwner[dest].size() << endl;
//    comm_->Barrier();
    
    // pack data to be sent
    DRT::ParObject::AddtoPack(dataSend,failedVec[dest].size());
    for (size_t inode=0;inode<failedVec[dest].size();inode++)
    {
      StartpointData failed = failedVec[dest][inode];
      
      DRT::ParObject::AddtoPack(dataSend,failed.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(failed.movNode_.X()));
      DRT::ParObject::AddtoPack(dataSend,failed.movNode_.Owner());
      DRT::ParObject::AddtoPack(dataSend,failed.startGid_);
    }
    
    
    // clear the no more needed data
    failedVec[dest].clear();
    
    lengthSend[0] = dataSend.size();
    
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
    exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
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
    size_t numberOfNodes;
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
    for (size_t inode=0;inode<numberOfNodes;inode++)
    {
      int gid;
      LINALG::Matrix<nsd,1> coords;
      int owner;
      int startGid;
      
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      
      double coordinates[nsd];
      for (size_t dim=0;dim<nsd;dim++)
        coordinates[dim] = coords(dim);
      
      DRT::Node movNode(gid,coordinates,owner);
      
      StartpointData failed(movNode,startGid,myrank_); // startOwner is current proc
      failed_.push_back(failed);
    }
//    cout << "after sending on proc " << myrank_ << " sizes are "
//        << gidNew.size() << ", " << movGidNew.size() <<
//        ", " << movCoordsNew.size() << " and " <<
//        movOwnerNew.size() << endl;
    
    // check if sizes fit
//    if (gidNew.size()!=movGidNew.size() ||
//        gidNew.size()!=movCoordsNew.size() ||
//        gidNew.size()!=movOwnerNew.size())
//      dserror("sending of alternative starting point data failed!");
    
  } // end loop over processors
} // end exportfinalData



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                               winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::exportIterData(
  bool& procfinished
)
{
  const size_t nsd = 3; // 3 dimensions for a 3d fluid element
  
  // Initialization
  int dest = myrank_+1; // destination proc (the "next" one)
  if(myrank_ == (numproc_-1))
    dest = 0;
  
  int source = myrank_-1; // source proc (the "last" one)
  if(myrank_ == 0)
    source = numproc_-1;
  
  vector<char> dataSend; // vector including all data that has to be send to the next proc
  vector<int> lengthSend(1,0);
  int size_one = 1;
  
  
  
/*-------------------------*
 * convergence check first *
 *-------------------------*/
  for (int iproc=0;iproc<numproc_-1;iproc++)
  {
    // packing the data
    DRT::ParObject::AddtoPack(dataSend,procfinished);
    lengthSend[0] = dataSend.size();
    
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
    exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    vector<char> dataRecv(lengthRecv[0]);
    exporter_.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter_.Wait(req_data);
    
    // unpacking the data
    size_t posinData = 0;
    bool procfinishedNew;
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,procfinishedNew);
    
    // setting converged
    if (!procfinishedNew)
      procfinished = false;
  }
  
  if (!procfinished) // just send data if its needed
  {
    dataSend.clear();
    
    // packing the data
    DRT::ParObject::AddtoPack(dataSend,next_.size());
    for (size_t ipoint=0;ipoint<next_.size();ipoint++)
    {
      StartpointData next = next_[ipoint];
      
      DRT::ParObject::AddtoPack(dataSend,next.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(next.movNode_.X()));
      DRT::ParObject::AddtoPack(dataSend,next.movNode_.Owner());
      DRT::ParObject::AddtoPack(dataSend,next.startpoint_);
      DRT::ParObject::AddtoPack(dataSend,next.changed_);
      DRT::ParObject::AddtoPack(dataSend,next.phiSign_);
      DRT::ParObject::AddtoPack(dataSend,next.searchedProcs_);
      DRT::ParObject::AddtoPack(dataSend,next.iter_);
      DRT::ParObject::AddtoPack(dataSend,next.startGid_);
      DRT::ParObject::AddtoPack(dataSend,next.startOwner_);
      
//      cout << "data on proc " << myrank_ << " before exporting:\n" << next.movNode_ <<
//          " with startpoint " << next.startpoint_ << "changed " << next.changed_ <<
//          ", phisign " << next.phiSign_ << ", searchedprocs " << next.searchedProcs_ <<
//          ", iter " << next.iter_ << ", startgid " << next.startGid_ <<
//          " and startowner " << next.startOwner_ << endl;
    }
    
    lengthSend[0] = dataSend.size();
    
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
    exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
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
    size_t numberOfNodes = 0;
    
    // clear curr_ - vector so it can be filled with new values
    curr_.clear();
    
    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
    for (size_t inode=0;inode<numberOfNodes;inode++)
    {
      int gid;
      LINALG::Matrix<nsd,1> coords;
      int owner;
      LINALG::Matrix<nsd,1> startpoint;
      int changed;
      int phiSign;
      int searchedProcs;
      int iter;
      int startGid;
      int startOwner;
      
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,changed);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiSign);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,iter);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
      
      double coordinates[nsd];
      for (size_t dim=0;dim<nsd;dim++)
        coordinates[dim] = coords(dim);
      
      DRT::Node movNode(gid,coordinates,owner);
      
      StartpointData curr(
          movNode,startpoint,changed,phiSign,
          searchedProcs,iter,startGid,startOwner);
      curr_.push_back(curr);
//      cout << "data on proc " << myrank_ << " after exporting:\n" << curr.movNode_ <<
//          " with startpoint " << curr.startpoint_ << "changed " << curr.changed_ <<
//          ", phisign " << curr.phiSign_ << ", searchedprocs " << curr.searchedProcs_ <<
//          ", iter " << curr.iter_ << ", startgid " << curr.startGid_ <<
//          " and startowner " << curr.startOwner_ << endl;
    }
    
    // check if all sizes fit
//    if (gids.size()!=numberOfNodes ||
//        coords.size()!=numberOfNodes ||
//        owners.size()!=numberOfNodes ||
//        curr_.startpoints_.size()!=numberOfNodes ||
//        curr_.changed_.size()!=numberOfNodes ||
//        curr_.phiSign_.size()!=numberOfNodes ||
//        curr_.searchedProcs_.size()!=numberOfNodes ||
//        curr_.iter_.size()!=numberOfNodes ||
//        curr_.startGid_.size()!=numberOfNodes ||
//        curr_.startpoints_.size()!=numberOfNodes)
//      dserror("sending of starting point data failed!");
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
  for (size_t inode=0;inode<done_.size();inode++)
  {
    doneVec[done_[inode].movNode_.Owner()].push_back(done_[inode]);
  }
  
  // clear final vectors
  done_.clear();
  
  // set gid, vel and pres for own processor
  done_ = doneVec[myrank_];
  
  // clear data about current proc
  doneVec[myrank_].clear();
  
  // send data to the processor where the point lies
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    vector<char> dataSend; // vector including all data that has to be send to dest proc
    vector<int> lengthSend(1,0);
    int size_one = 1;
    
    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;
    
    dataSend.clear();
    
    // pack data to be sent
    DRT::ParObject::AddtoPack(dataSend,doneVec[dest].size());
    for (size_t inode=0;inode<doneVec[dest].size();inode++)
    {
      StartpointData done = doneVec[dest][inode];
      
      DRT::ParObject::AddtoPack(dataSend,done.movNode_.Id());
      DRT::ParObject::AddtoPack(dataSend,done.velValues_);
      DRT::ParObject::AddtoPack(dataSend,done.presValues_);
    }
    
    // clear the no more needed data
    doneVec[dest].clear();
    
    lengthSend[0] = dataSend.size();
    
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
    exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
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
    size_t numberOfNodes;
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
    for (size_t inode=0;inode<numberOfNodes;inode++)
    {
      int gid;
      vector<LINALG::Matrix<nsd,1> > velValues;
      vector<double> presValues;
      
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velValues);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presValues);
      
      StartpointData done(*discret_->gNode(gid),velValues,presValues); // startOwner is current proc
      done_.push_back(done);
    }
    
//    // check if sizes fit
//    if (velNew.size()!=presNew.size() ||
//        velNew.size()!=gidNew.size())
//      dserror("sending of final starting point data failed!");
  } // end loop over processors
} // end exportfinalData



/*------------------------------------------------------------------------------------------------*
 * export enrichment data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Startvalues::exportEnrichmentData()
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
  
  
//  cout << "data that should be sent" << endl;
  
  // packing the data
  DRT::ParObject::AddtoPack(dataSend,enr_.size());
  for (size_t ipoint=0;ipoint<enr_.size();ipoint++)
  {
    StartpointData enr = enr_[ipoint];
    
//    cout << enr.movNode_ << "\nwith phisign " << enr.phiSign_ << ", dist " <<
//        enr.dMin_ << ", phivalue " << enr.phiValue_ << " and enrvalues:" << endl;
//    for (std::map<DofKey<onNode>,vector<double> >::const_iterator j = enr.enrValues_.begin();
//        j != enr.enrValues_.end();j++)
//      cout << j->first << " with values " << j->second[0] << " and " << j->second[2] << endl;
    
    vector<char> data;
    DRT::ParObject::AddtoPack(dataSend,enr.movNode_.Id());
    DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(enr.movNode_.X()));
    DRT::ParObject::AddtoPack(dataSend,enr.movNode_.Owner());
    DRT::ParObject::AddtoPack(dataSend,enr.phiSign_);
    DRT::ParObject::AddtoPack(dataSend,enr.dMin_);
    DRT::ParObject::AddtoPack(dataSend,enr.phiValue_);
    DRT::ParObject::AddtoPack(dataSend,enr.enrValues_.size());

    // add map to pack
    for (std::map<DofKey<onNode>,vector<double> >::const_iterator currEnrVal = enr.enrValues_.begin();
        currEnrVal != enr.enrValues_.end();
        currEnrVal++)
    {
      vector<char> data;
      currEnrVal->first.Pack(data);
      
      DRT::ParObject::AddtoPack(dataSend,data);
      DRT::ParObject::AddtoPack(dataSend,currEnrVal->second);
    }
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
  // ... and receive length
  vector<int> lengthRecv(1,0);
  exporter_.Receive(source, length_tag, lengthRecv, size_one);
  exporter_.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter_.ISend(myrank_, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
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
  enr_.clear();
  
  // unpack received data
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,numberOfNodes);
  for (size_t inode=0;inode<numberOfNodes;inode++)
  {
    int gid;
    LINALG::Matrix<nsd,1> coords;
    int owner;
    int phiSign;
    double dMin;
    double phiValue;
    size_t numEnr;
    map<DofKey<onNode>,vector<double> > enrValues;
//    cout << "here1 on proc " << myrank_ << endl;
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,gid);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,owner);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiSign);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,dMin);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,numEnr);
//    cout << "here2 on proc " << myrank_ << endl;
    // extract map-data
    for (size_t currEnrVal = 0; currEnrVal < numEnr; currEnrVal++)
    {
      vector<char> data;
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,data);
      DofKey<onNode> dofkey(data);
      
      vector<double> enrVals;
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,enrVals);
      
      enrValues.insert(make_pair(dofkey,enrVals));
    }
//    cout << "here4 on proc " << myrank_ << endl;
    double coordinates[nsd];
    for (size_t dim=0;dim<nsd;dim++)
      coordinates[dim] = coords(dim);
    
    DRT::Node movNode(gid,coordinates,owner);
    
    StartpointData enr(movNode,phiSign,dMin,phiValue,enrValues);
    enr_.push_back(enr);
  }
  
//  discret_->Comm().Barrier();
//  cout << "received data:" << endl;
//  for (size_t ipoint=0;ipoint<enr_.size();ipoint++)
//  {
//    StartpointData enr = enr_[ipoint];
//    cout << enr.movNode_ << "\nwith phisign " << enr.phiSign_ << ", dist " <<
//        enr.dMin_ << ", phivalue " << enr.phiValue_ << " and enrvalues:" << endl;
//    for (std::map<DofKey<onNode>,vector<double> >::const_iterator j = enr.enrValues_.begin();
//        j != enr.enrValues_.end();j++)
//      cout << j->first << " with values " << j->second[0] << " and " << j->second[2] << endl;
//  }
  
  // check if all sizes fit
//  if (gids.size()!=numberOfNodes ||
//      coords.size()!=numberOfNodes ||
//      owners.size()!=numberOfNodes ||
//      curr_.startpoints_.size()!=numberOfNodes ||
//      dMin.size()!=numberOfNodes ||
//      curr_.changed_.size()!=numberOfNodes ||
//      curr_.phiSign_.size()!=numberOfNodes ||
//      curr_.startGid_.size()!=numberOfNodes ||
//      curr_.startOwner_.size()!=numberOfNodes)
//    dserror( "sending of starting point data failed!");
} // end exportEnrichmentData
  //TODO any barriers needed in exportfiles? (with more than 2 procs)

#endif // parallel

#endif // CCADISCRET



