/*!----------------------------------------------------------------------
\file drt_periodicbc.cpp

\brief  update dofrowmap and dofsets for periodic boundary conditions

        o effects distribution of master and slave nodes. Master and
          slave are owned by one proc afterwards.
        o effects ghosting of nodes. Master and slave nodes are ghosted
          in pairs.
        o passes list of coupled nodes to the dofset
        

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE


#include "drt_periodicbc.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
PeriodicBoundaryConditions::PeriodicBoundaryConditions
  (RefCountPtr<DRT::Discretization> actdis)
  : discret_(actdis)
{
  // get periodic surface boundary conditions
  discret_->GetCondition("SurfacePeriodic",mysurfpbcs_);

  // set number of pairs of periodic boundary conditions
  numpbcpairs_ = mysurfpbcs_.size()/2;

  // create map that will be connecting master to slave nodes owned by
  // this proc
  //       master node -> list of his slave node(s)
  allcoupledrownodes_=rcp(new map<int,vector<int> >);
  
  // create map that will be connecting master to slave nodes owned or
  // ghosted by this proc
  //       master node -> list of his slave node(s)
  allcoupledcolnodes_=rcp(new map<int,vector<int> >);
  
  return;

}// PeriodicBoundaryConditions(RefCountPtr<DRT::Discretization> actdis)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Proceed all pairs of periodic boundary conditions and create         |
 | complete coupling map                                     gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void PeriodicBoundaryConditions::UpdateDofsForPeriodicBoundaryConditions()
{


  // map from global masternodeids (on this proc) to global slavenodeids
  // for a single condition
  map<int,int> midtosid;

  // pointers to master and slave condition
  DRT::Condition* mastercond;
  DRT::Condition* slavecond;

  // global master node Ids and global slave node Ids
  vector <int>* masternodeids;
  vector <int>* slavenodeids;
  
  //----------------------------------------------------------------------
  //          LOOP PAIRS OF PERIODIC BOUNDARY CONDITIONS
  //----------------------------------------------------------------------

  // loop pairs of periodic boundary conditions
  for (int pbcid=0;pbcid<numpbcpairs_;++pbcid)
  {
    //--------------------------------------------------
    // get master and slave condition pair
    {
      for (unsigned numcond=0;numcond<mysurfpbcs_.size();++numcond)
      {
        vector<int>* myid
          = mysurfpbcs_[numcond]->Get<vector<int> >("Id of periodic boundary condition");
        if (myid[0][0] == pbcid)
        {
          vector<int>* mymasterslavetoggle
            = mysurfpbcs_[numcond]->Get<vector<int> >("Is slave periodic boundary condition");
          
          if(mymasterslavetoggle[0][0]==0)
          {
            mastercond =mysurfpbcs_[numcond];
          }
          else if (mymasterslavetoggle[0][0]==1)
          {
            slavecond =mysurfpbcs_[numcond];
          }
          else
          {
            dserror("pbc is neither master nor slave");
          }
        }
      }
    }

    //--------------------------------------------------
    // vector specifying the plane of this pair
    //
    //     |                           |
    //     |                           |
    //     |      parallel planes      |
    //     |-------------------------->|
    //     |                           |
    //     |                           |
    //     |                           |
    //   slave                      master
    //
    //
    vector <int>* dofsforpbcplane  =
      mastercond->Get<vector<int> >("degrees of freedom for the pbc plane");
        

    //--------------------------------------------------
    // get global master node Ids and global slave node Ids
    masternodeids = mastercond->Get<vector<int> >("Node Ids");
    slavenodeids  = slavecond ->Get<vector<int> >("Node Ids");


    //----------------------------------------------------------------------
    //      CONSTRUCT NODE MATCHING BETWEEN MASTER AND SLAVE NODES
    //                        FOR THIS CONDITION
    //----------------------------------------------------------------------
    
    // clear map from global masternodeids (on this proc) to global
    // slavenodeids --- it belongs to this pair of periodic boundary
    // conditions!!!
    midtosid.clear();

    // get map master on this proc -> slave on some proc
    CreateNodeCouplingForSinglePBC(
      midtosid,
      *masternodeids,
      *slavenodeids ,
      *dofsforpbcplane);

    // Add the connectivity from this condition to the connectivity
    // of all previously processed periodic boundary conditions.
    // Redistribute the nodes (rownodes+ghosting)
    // Assign the same degrees of freedom to coupled nodes
    AddConnectivityRedistributeAndCreateDofCoupling(
      midtosid,pbcid);
      
  } // end loop pairs of periodic boundary conditions

  // eventually call METIS to optimally distribute the nodes --- up to
  // now, a periodic boundary condition might remove all nodes from a
  // proc ...
  
  return;
  
}// UpdateDofsForPeriodicBoundaryConditions()



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Couple nodes for specific pair of periodic boundary conditions       |
 |                                                  (public) gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void PeriodicBoundaryConditions::CreateNodeCouplingForSinglePBC(
  map<int,int>        &midtosid,
  const vector <int>   masternodeids,
  const vector <int>   slavenodeids,
  const vector <int>   dofsforpbcplane
  )
{


  // these are just parameter definitions for the octree search algorithm
  double tol            = 1E-12;
  int    maxnodeperleaf = 25;

  //----------------------------------------------------------------------
  //                   BUILD PROCESSOR LOCAL OCTREE
  //----------------------------------------------------------------------
  // build processor local octree
  NodeMatchingOctree nodematchingoctree(
    discret_      ,
    masternodeids ,
    maxnodeperleaf,
    tol
    );
      

  //----------------------------------------------------------------------
  //  SEARCH CLOSEST NODES IN OCTREES ON ALL PROCESSORS
  //----------------------------------------------------------------------

  // create connectivity for this condition in this direction
  {

    // create map from gid masternode -> gid corresponding slavenode
    nodematchingoctree.CreateGlobalNodeMatching(
      slavenodeids   ,
      dofsforpbcplane,
      midtosid
      );
  }

  return;
} // PeriodicBoundaryConditions::CreateNodeCouplingForSinglePBC


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Add the connectivity from this condition to the connectivity         |
 | of all previously processed periodic boundary conditions.            |
 | Redistribute the nodes and assign the dofs to the                    |
 | current distribution of nodes                    (public) gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void PeriodicBoundaryConditions::AddConnectivityRedistributeAndCreateDofCoupling(
  map<int,int>       &midtosid,
  const int           pbcid
  )
{

  // the "inverse" mapping of allcoupled(row/col)nodes
  //       slave node -> his master node (list of size 1)
  RefCountPtr<map<int,vector<int> > > inversenodecoupling;
  inversenodecoupling=rcp(new map<int,vector<int> >);

  
  // rcp to the constructed rowmap
  RefCountPtr<Epetra_Map> newrownodemap;
  
  //----------------------------------------------------------------------
  //  ADD THE CONNECTIVITY FROM THIS CONDITION TO THE CONNECTIVITY OF
  //                           ALL CONDITIONS
  //----------------------------------------------------------------------
  {
    map<int,int>::iterator iter;
    int  masterid;
    int  slaveid;
    bool alreadyinlist;
        
    for( iter = midtosid.begin(); iter != midtosid.end(); ++iter )
    {
      // get id of masternode and slavenode
      masterid  = iter->first;
      slaveid   = iter->second;

      // is masterid already in allcoupledrownodes?
      {
        alreadyinlist=false;

        map<int,vector<int> >::iterator found;

        found = allcoupledrownodes_->find(masterid);
        if(found != allcoupledrownodes_->end())
        {
          // masterid is already in the list --- i.e., the master is the
          // master of a previous condition. Simply append the slave id here
          alreadyinlist=true;
          found->second.push_back(slaveid);
        }
        // masterid is not in the list yet. -> new entry
        if (alreadyinlist==false)
        {
          (*allcoupledrownodes_)[masterid].push_back(slaveid);
        } // end if not in map
      }
    } // end insert entries of midtosid into the allcoupledrownodes map

    //---------------------------------------------------------------
    //     COMPLETE MATCHING FOR NODES WITH MORE THAN ONE PBC
    //---------------------------------------------------------------
    // complete matching (we are able to do this because of the
    // communication step which was done in the last iteration
    // (first step is no problem at all ...)
    // this might be expensive --- but since I do not expect to have
    // many periodic boundary conditions (more than three) and the
    // number of nodes concerned is very small (edge and cornernodes,
    // for example in a 32x32x32 channel flow we talk about 33 nodes)
    // I think we could afford these intermediate redistributions
    // 
    {
      map<int,vector<int> >::iterator curr;
      vector<int >         ::iterator alreadyin;

      bool toggle = false;
      for( curr = allcoupledrownodes_->begin(); curr != allcoupledrownodes_->end(); ++curr )
      {
        // we have a node which is coupled to more then one nodes
        if(curr->second.size()>1)
        {
          map<int,vector<int> >::iterator found;

          // the list of matching nodes may be incomplete
          for (unsigned rr=0;rr<curr->second.size();++rr)
          {
            found = allcoupledrownodes_->find(curr->second[rr]);
                
            if(found != allcoupledrownodes_->end())
            {
              bool infoused=false;
              // close the connectivity using the slave node which was the
              // masternode of the previous condition
              for (unsigned mm = 0;mm<found->second.size();++mm)
              {
                toggle = false;
                // add only if not already in list
                for( alreadyin   = curr->second.begin();
                     alreadyin  != curr->second.end()  ;
                     ++alreadyin )
                {
                  if(*alreadyin == found->second[mm])
                  {
                    toggle  = true;
                    infoused= true;
                    break; // iter alreadyin
                  }
                }

                if(toggle!=true)
                {
                  curr->second.push_back(found->second[mm]);
                }
              }
              // the information was used and the entry is not needed anymore
              if(infoused=true)
              {
                found->second.clear();
              }
            } // end if we have a further connectivity information...
          } // end loop all coupled nodes to this one
        } // end if we have a node which is coupled to more then one nodes
      } // end for curr

      // erase all redundant (now empty) connections
      for( curr = allcoupledrownodes_->begin(); curr != allcoupledrownodes_->end();)
      {
        if(curr->second.empty())
        {
          allcoupledrownodes_->erase(curr++);
        }
        else
        {
          ++curr;
        }
      }
    } // end complete matching
  }
      
  //----------------------------------------------------------------------
  // REDISTRIBUTE ALL MASTER/SLAVE PAIRS OF THIS CONDITION TO THE
  //                      SAME PROCESSOR
  //----------------------------------------------------------------------
  {
    // a list of all nodes on this proc
    vector<int> nodesonthisproc(discret_->NodeRowMap()->NumMyElements());
      
    // get all node gids of nodes on this proc
    discret_->NodeRowMap()->MyGlobalElements(&nodesonthisproc[0]);

    // remove all node gids of slave nodes on this proc
    {
      vector<int>::iterator curr;   
      for( curr = nodesonthisproc.begin(); curr != nodesonthisproc.end(); )
      {
        // get the node if it's available on the processor
        if(discret_->HaveGlobalNode(*curr))
        {
          DRT::Node* actnode = discret_->gNode(*curr);

          // get all periodic boundary conditions on this node
          vector<DRT::Condition*> thiscond;
          actnode->GetCondition("SurfacePeriodic",thiscond);

          // loop them and check, whether this is a pbc slave node of the
          // current condition
          vector<int>* myid;
            
          bool erased=false;
            
          for (unsigned numcond=0;numcond<thiscond.size();++numcond)
          {
            myid  = thiscond[numcond]->Get<vector<int> >("Id of periodic boundary condition");

            if (myid[0][0] == pbcid)
            {
              vector<int>* mymasterslavetoggle
                = thiscond[numcond]->Get<vector<int> >("Is slave periodic boundary condition");
          
              if(mymasterslavetoggle[0][0]==1)
              {
                // erase the coupled nodes from the map --- they are redundant
                allcoupledrownodes_->erase(*curr);
                // erase id from vector --- be careful, curr is increased by return value!
                curr=nodesonthisproc.erase(curr);
                erased=true;
                break;
              } // end is slave?
            }// end if (myid[0][0] == pbcid)
          } // end loop this conditions

          if(erased==false)
          {
            ++curr;
          }
            
        } // is the node available on the proc?
        else
        {
          dserror("RowNode not available on proc");
        }
      }// iteration over all nodes associated with this proc
    }
        
    // append slavenodes to this list of nodes on this proc
    {
      map<int,int>::iterator iter;
      int slaveid;
      for( iter = midtosid.begin(); iter != midtosid.end(); ++iter )
      {
        // get id of slavenodes
        slaveid  = iter->second;     
        nodesonthisproc.push_back(slaveid);
      } // end append slavenodes of this proc
    }

    //--------------------------------------------------
    // build noderowmap for new distribution of nodes
    newrownodemap = rcp(new Epetra_Map(discret_->NumGlobalNodes(),
                                       nodesonthisproc.size(),
                                       &nodesonthisproc[0],
                                       0,
                                       discret_->Comm()));

    // create nodal graph of problem, according to old RowNodeMap
    RefCountPtr<Epetra_CrsGraph> oldnodegraph = discret_->BuildNodeGraph();

    // export the graph to newrownodemap
    Epetra_CrsGraph nodegraph(Copy,*newrownodemap,108,false);

    {
      Epetra_Export exporter(*discret_->NodeRowMap(),*newrownodemap);
      int err = nodegraph.Export(*oldnodegraph,exporter,Add);
      if (err<0) dserror("Graph export returned err=%d",err);
    }
    nodegraph.FillComplete();
    nodegraph.OptimizeStorage();
    
    // build nodecolmap for new distribution of nodes
    const Epetra_BlockMap cntmp = nodegraph.ColMap();

    RefCountPtr<Epetra_Map> newcolnodemap;

    newcolnodemap = rcp(new Epetra_Map(-1,
                                       cntmp.NumMyElements(),
                                       cntmp.MyGlobalElements(),
                                       0,
                                       discret_->Comm()));

    // create the inverse map --- slavenode -> masternode
    inversenodecoupling->clear();
        
    for(map<int,vector<int> >::iterator curr = allcoupledrownodes_->begin();
        curr != allcoupledrownodes_->end();
        ++curr )
    {
      for(unsigned rr=0;rr<curr->second.size();++rr)
      {
        (*inversenodecoupling)[curr->second[rr]].push_back(curr->first);
      }
    }
    //----------------------------------------------------------------------
    //       GHOSTED NODES NEED INFORMATION ON THEIR COUPLED NODES
    //----------------------------------------------------------------------
    *allcoupledcolnodes_=(*allcoupledrownodes_);
#ifdef PARALLEL
    {
      // create an exporter 
      DRT::Exporter exportconnectivity(*newrownodemap,
                                       *newcolnodemap,
                                       discret_->Comm());
          
      // export information on all master->slave couplings (with multiple
      // couplings)
      exportconnectivity.Export(*allcoupledcolnodes_);

      // export the inverse slave->master matching without multiple couplings
      exportconnectivity.Export(*inversenodecoupling);
    }
#endif

    // to assign the degrees of freedom, we have to make sure that coupled
    // nodes are only ghosted in pairs --- so modify the colmap
    {
      // mycolnodes contains all nodes which will be stored on this proc
      // according to the colmap constructed
      vector<int> mycolnodes(newcolnodemap->NumMyElements());
      newcolnodemap->MyGlobalElements(&(mycolnodes[0]));

      // determine all ghosted slave nodes in this vector which do not have
      // a ghosted master on this proc --- we have to fetch it to be able
      // to assign the dofs
      for (map<int,vector<int> >::iterator curr=inversenodecoupling->begin();
           curr!=inversenodecoupling->end();
           ++curr)
      {
        if(curr->second.empty())
        {
          dserror("inverse slave-master matching incomplete");
        }
        int mymaster=curr->second[0];
        if(newcolnodemap->LID(mymaster)<0)
        {
          // was master already added to the list of (ghosted) nodes?
          vector<int>::iterator found;
          found=find(mycolnodes.begin(),mycolnodes.end(),mymaster);
          // no, it's not inside
          if(found==mycolnodes.end())
          {
            mycolnodes.push_back(mymaster);
          }
        }
      }

      // now reconstruct the extended colmap 
      newcolnodemap = rcp(new Epetra_Map(-1,
                                         mycolnodes.size(),
                                         &mycolnodes[0],
                                         0,
                                         discret_->Comm()));

      *allcoupledcolnodes_=(*allcoupledrownodes_);
#ifdef PARALLEL
      // the new master-ghost nodes need their information about
      // connectivity
      {
        // create an exporter 
        DRT::Exporter exportconnectivity(*newrownodemap,
                                         *newcolnodemap,
                                         discret_->Comm());
        // export information on all slave->master couplings (with multiple
        // couplings)
        exportconnectivity.Export(*allcoupledcolnodes_);
      }
#endif
    }


    // Create a new DofSet special for Periodic boundary conditions with
    // this type of node coupling

    // create a new dofset specialisation for periodic boundary conditions

    discret_->ReplaceDofSet(rcp(new DRT::PBCDofSet(allcoupledcolnodes_)));
//    discret_->SetPeriodicBoundaryConditionCouples();

    //--------------------------------------------------
    // redistribute the nodes
    //
    // this contains a call to FillComplete and assigns the same
    // degree of freedom to the matching nodes
    discret_->Redistribute(*newrownodemap,*newcolnodemap);

    // throw away old nodegraph
    oldnodegraph = null;
  }

  return;
}// PeriodicBoundaryConditions::AddConnectivityRedistributeAndCreateDofCoupling


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor dtor  (public)                                 gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
PeriodicBoundaryConditions::~PeriodicBoundaryConditions()
{
  return;
}// ~PeriodicBoundaryConditions()

#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */
