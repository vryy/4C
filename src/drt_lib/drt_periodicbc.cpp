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

  if(mysurfpbcs_.empty())
  {
    discret_->GetCondition("LinePeriodic",mysurfpbcs_);
  }

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


  if(numpbcpairs_>0)
  {
    // -------------------------------------------------------------------
    // create timers and time monitor
    // -------------------------------------------------------------------
    timepbctot_         = TimeMonitor::getNewTimer("0) pbc routine total"                                 );
    timepbcmidtosid_    = TimeMonitor::getNewTimer("1)   +create midtosid maps"                           );
    timepbcmidoct_      = TimeMonitor::getNewTimer("2)      +build local octrees"                         );
    timepbcmidmatch_    = TimeMonitor::getNewTimer("3)      +search closest nodes in octrees on all procs");
    timepbcaddcon_      = TimeMonitor::getNewTimer("4)   +add connectivity to previous conditions"        );
    timepbcreddis_      = TimeMonitor::getNewTimer("5)   +Redistribute the nodes"                         );
    timepbcmakeghostmap_= TimeMonitor::getNewTimer("6)      +build rowmap and temporary colmap"           );
    timepbcghost_       = TimeMonitor::getNewTimer("7)      +repair ghosting"                             );
    timepbcrenumdofs_   = TimeMonitor::getNewTimer("8)      +call discret->Redistribute"                  );
  }



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
  if(numpbcpairs_>0)
  {
    // time measurement --- start TimeMonitor tm0
    tm0_ref_        = rcp(new TimeMonitor(*timepbctot_ ));

    // map from global masternodeids (on this proc) to global slavenodeids
    // for a single condition
    map<int,int> midtosid;

    // pointers to master and slave condition
    DRT::Condition* mastercond=NULL;
    DRT::Condition* slavecond =NULL;

    // global master node Ids and global slave node Ids
    const vector <int>* masternodeids;
    const vector <int>* slavenodeids;

    //----------------------------------------------------------------------
    //          LOOP PAIRS OF PERIODIC BOUNDARY CONDITIONS
    //----------------------------------------------------------------------


    // loop pairs of periodic boundary conditions
    for (int pbcid=0;pbcid<numpbcpairs_;++pbcid)
    {
      if (discret_->Comm().MyPID() == 0)
      {
        cout << "pbc pair " << pbcid << ": ";
        fflush(stdout);
      }

      // time measurement --- start TimeMonitor tm1
      tm1_ref_        = rcp(new TimeMonitor(*timepbcmidtosid_ ));

      //--------------------------------------------------
      // get master and slave condition pair
      {
        for (unsigned numcond=0;numcond<mysurfpbcs_.size();++numcond)
        {
          const vector<int>* myid
            = mysurfpbcs_[numcond]->Get<vector<int> >("Id of periodic boundary condition");
          if (myid[0][0] == pbcid)
          {
            const vector<int>* mymasterslavetoggle
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
      const vector <int>* dofsforpbcplane  =
        mastercond->Get<vector<int> >("degrees of freedom for the pbc plane");


      //--------------------------------------------------
      // get global master node Ids and global slave node Ids
      masternodeids = mastercond->Nodes();
      slavenodeids  = slavecond ->Nodes();


      //----------------------------------------------------------------------
      //      CONSTRUCT NODE MATCHING BETWEEN MASTER AND SLAVE NODES
      //                        FOR THIS CONDITION
      //----------------------------------------------------------------------

      // clear map from global masternodeids (on this proc) to global
      // slavenodeids --- it belongs to this pair of periodic boundary
      // conditions!!!
      midtosid.clear();

      if (discret_->Comm().MyPID() == 0)
      {
        cout << " creating midtosid-map ... ";
        fflush(stdout);
      }

      // get map master on this proc -> slave on some proc
      CreateNodeCouplingForSinglePBC(
        midtosid,
        *masternodeids,
        *slavenodeids ,
        *dofsforpbcplane);
      // time measurement --- this causes the TimeMonitor tm1 to stop here
      tm1_ref_ = null;

      if (discret_->Comm().MyPID() == 0)
      {
        cout << "adding connectivity to previous pbcs ... ";
        fflush(stdout);
      }

      // time measurement --- start TimeMonitor tm4
      tm4_ref_        = rcp(new TimeMonitor(*timepbcaddcon_ ));


      //----------------------------------------------------------------------
      //      ADD CONNECTIVITY TO CONNECTIVITY OF ALL PREVIOUS PBCS
      //----------------------------------------------------------------------
      // Add the connectivity from this condition to the connectivity
      // of all previously processed periodic boundary conditions.
      // Redistribute the nodes (rownodes+ghosting)
      // Assign the same degrees of freedom to coupled nodes
      AddConnectivity(midtosid,pbcid);

      // time measurement --- this causes the TimeMonitor tm4 to stop here
      tm4_ref_ = null;

      if (discret_->Comm().MyPID() == 0)
      {
        cout << " done\n";
        fflush(stdout);
      }
    } // end loop pairs of periodic boundary conditions


    //----------------------------------------------------------------------
    //         REDISTRIBUTE ACCORDING TO THE GENERATED CONNECTIVITY
    //----------------------------------------------------------------------

    // time measurement --- start TimeMonitor tm5
    tm5_ref_        = rcp(new TimeMonitor(*timepbcreddis_ ));


    if (discret_->Comm().MyPID() == 0)
    {
      cout << "Redistributing ... ";
      fflush(stdout);
    }

    RedistributeAndCreateDofCoupling();

    if (discret_->Comm().MyPID() == 0)
    {
      cout << " done\n";
      fflush(stdout);
    }

    // time measurement --- this causes the TimeMonitor tm5 to stop here
    tm5_ref_ = null;

    // eventually call METIS to optimally distribute the nodes --- up to
    // now, a periodic boundary condition might remove all nodes from a
    // proc ...


    // time measurement --- this causes the TimeMonitor tm0 to stop here
    //                                                (call of destructor)
    tm0_ref_ = null;


    if(discret_->Comm().MyPID()==0)
    {
      cout<<endl<<endl;
    }
    TimeMonitor::summarize();

    if(discret_->Comm().MyPID()==0)
    {
      cout<<endl<<endl;
    }
  } // end if numpbcpairs_>0
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
  double tol            = 1E-9;
  int    maxnodeperleaf = 25;

  //----------------------------------------------------------------------
  //                   BUILD PROCESSOR LOCAL OCTREE
  //----------------------------------------------------------------------
  // time measurement --- start TimeMonitor tm2
  tm2_ref_        = rcp(new TimeMonitor(*timepbcmidoct_ ));

  // build processor local octree
  DRT::Utils::NodeMatchingOctree nodematchingoctree(
    *discret_     ,
    masternodeids ,
    maxnodeperleaf,
    tol
    );
  // time measurement --- this causes the TimeMonitor tm2 to stop here
  tm2_ref_ = null;

  //----------------------------------------------------------------------
  //  SEARCH CLOSEST NODES IN OCTREES ON ALL PROCESSORS
  //----------------------------------------------------------------------

  // time measurement --- start TimeMonitor tm3
  tm3_ref_        = rcp(new TimeMonitor(*timepbcmidmatch_ ));
  // create connectivity for this condition in this direction
  {

    // create map from gid masternode -> gid corresponding slavenode
    nodematchingoctree.CreateGlobalNodeMatching(
      slavenodeids   ,
      dofsforpbcplane,
      midtosid
      );
  }

  // time measurement --- this causes the TimeMonitor tm3 to stop here
  tm3_ref_ = null;

  return;
} // PeriodicBoundaryConditions::CreateNodeCouplingForSinglePBC


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Add the connectivity from this condition to the connectivity         |
 | of all previously processed periodic boundary conditions.            |
 |                                                  (public) gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void PeriodicBoundaryConditions::AddConnectivity(
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
    // complete matching --- we are able to do this because of the
    // communication step (first step is no problem at all ...)

    {
      if(pbcid>0)
      {

        // 1) each proc generates a list of his multiple coupled masters
        // 2) the list is communicated in a round robin pattern to all the
        //    other procs.
        // 3) the proc checks the package from each proc and inserts missing
        //    links into the multiple coupling


        //--------------------------------------------------------------------
        // -> 1) create a list of multiple master
        // Communicate multiple couplings for completion...
        vector < vector <int> > multiplecouplings;
        for( iter = midtosid.begin(); iter != midtosid.end(); ++iter )
        {
          // get id of masternode and the node itself
          masterid  = iter->first;
          DRT::Node* actnode = discret_->gNode(masterid);

          // get all periodic boundary conditions on this node
          vector<DRT::Condition*> thiscond;
          actnode->GetCondition("SurfacePeriodic",thiscond);

          if(thiscond.empty())
          {
            actnode->GetCondition("LinePeriodic",thiscond);
          }

          // loop them and check, whether this is a pbc pure master node
          // for all previous conditions
          unsigned ntimesmaster = 0;
          for (unsigned numcond=0;numcond<thiscond.size();++numcond)
          {
            const vector<int>* mymasterslavetoggle
              = thiscond[numcond]->Get<vector<int> >("Is slave periodic boundary condition");

            if(mymasterslavetoggle[0][0]==0)
            {
              ++ntimesmaster;
            } // end is slave?
          } // end loop this conditions

          if(ntimesmaster==thiscond.size())
          {
            // yes, we have such a pure master node
            vector <int> thiscoupling;
            thiscoupling.push_back(masterid);
            for(vector<int>::iterator rr=(*allcoupledrownodes_)[masterid].begin();
                rr!=(*allcoupledrownodes_)[masterid].end();
                ++rr)
            {
              thiscoupling.push_back(*rr);
            }

            // add it to the list of multiple coupled masters on this proc
            multiplecouplings.push_back(thiscoupling);
          }
        }

        //--------------------------------------------------------------------
        // -> 2) round robin loop

#ifdef PARALLEL
        int myrank  =discret_->Comm().MyPID();
#endif
        int numprocs=discret_->Comm().NumProc();

        vector<char> sblock;
        vector<char> rblock;


#ifdef PARALLEL
        // create an exporter for point to point comunication
        DRT::Exporter exporter(discret_->Comm());
#endif

        for (int np=0;np<numprocs;np++)
        {
          // pack multiple couplings
          sblock.clear();
          for(unsigned rr=0;rr<multiplecouplings.size();++rr)
          {
            DRT::ParObject::AddtoPack(sblock,multiplecouplings[rr]);
          }
#ifdef PARALLEL
          MPI_Request request;
          int         tag    =myrank;

          int         frompid=myrank;
          int         topid  =(myrank+1)%numprocs;

          int         length=sblock.size();

          exporter.ISend(frompid,topid,
                         &(sblock[0]),sblock.size(),
                         tag,request);

          rblock.clear();

          // receive from predecessor
          frompid=(myrank+numprocs-1)%numprocs;
          exporter.ReceiveAny(frompid,tag,rblock,length);

          if(tag!=(myrank+numprocs-1)%numprocs)
          {
            dserror("received wrong message (ReceiveAny)");
          }

          exporter.Wait(request);

          {
            // for safety
            exporter.Comm().Barrier();
          }
#else
          // dummy communication
          rblock.clear();
          rblock=sblock;
#endif

          //--------------------------------------------------
          // Unpack received block.
          multiplecouplings.clear();

          int index = 0;
          while (index < (int)rblock.size())
          {
            vector<int> onecoup;
            DRT::ParObject::ExtractfromPack(index,rblock,onecoup);
            multiplecouplings.push_back(onecoup);
          }

          //--------------------------------------------------
          // -> 3) Try to complete the matchings

          for(unsigned rr=0;rr<multiplecouplings.size();++rr)
          {
            for(unsigned mm=1;mm<multiplecouplings[rr].size();++mm)
            {
              int possiblemaster=(multiplecouplings[rr])[mm];

              map<int,vector<int> >::iterator found
                = allcoupledrownodes_->find(possiblemaster);

              if(found != allcoupledrownodes_->end())
              {
                // close the connectivity using the slave node which was the
                // masternode of the previous condition
                for (unsigned mm = 0;mm<found->second.size();++mm)
                {
                  bool dontdoit=false;
                  for(unsigned i = 1 ;i<multiplecouplings[rr].size();++i)
                  {
                    if(found->second[mm]==(multiplecouplings[rr])[i])
                    {
                      dontdoit=true;
                    }
                  }
                  if(!dontdoit)
                  {
                    multiplecouplings[rr].push_back(found->second[mm]);
                  }
                }
                allcoupledrownodes_->erase(found);
              } // end if we have a further connectivity information...
            }
          }
        }

        // add this information to the map of all coupled nodes
        for(unsigned rr=0;rr<multiplecouplings.size();++rr)
        {
          int multimaster=(multiplecouplings[rr])[0];

          for(unsigned mm=((*allcoupledrownodes_)[multimaster].size()+1);
              mm<multiplecouplings[rr].size();
              ++mm)
          {
            int thisslave = (multiplecouplings[rr])[mm];
            (*allcoupledrownodes_)[multimaster].push_back(thisslave);
          }
        }
      }
    } // end complete matching
  }

  return;
}// PeriodicBoundaryConditions::AddConnectivity


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Redistribute the nodes and assign the dofs to the                    |
 | current distribution of nodes                    (public) gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void PeriodicBoundaryConditions::RedistributeAndCreateDofCoupling(
  )
{
  // the "inverse" mapping of allcoupled(row/col)nodes
  //       slave node -> his master node (list of size 1)
  RefCountPtr<map<int,vector<int> > > inversenodecoupling;
  inversenodecoupling=rcp(new map<int,vector<int> >);


  // rcp to the constructed rowmap
  RefCountPtr<Epetra_Map> newrownodemap;

  {
    // time measurement --- start TimeMonitor tm6
    tm6_ref_        = rcp(new TimeMonitor(*timepbcmakeghostmap_ ));

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

          if(thiscond.empty())
          {
            actnode->GetCondition("LinePeriodic",thiscond);
          }

          // loop them and check, whether this is a pbc slave node of the
          // current condition
          bool erased=false;

          for (unsigned numcond=0;numcond<thiscond.size();++numcond)
          {
            const vector<int>* mymasterslavetoggle
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
      for(map<int,vector<int> >::iterator curr = allcoupledrownodes_->begin();
          curr != allcoupledrownodes_->end();
          ++curr )
      {
        for (vector<int>::iterator iter=curr->second.begin();iter!=curr->second.end();++iter)
        {
          int slaveid  = *iter;

          nodesonthisproc.push_back(slaveid);
        }
      }
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

    // time measurement --- this causes the TimeMonitor tm6 to stop here
    tm6_ref_ = null;


    // time measurement --- start TimeMonitor tm7
    tm7_ref_        = rcp(new TimeMonitor(*timepbcghost_ ));

    //----------------------------------------------------------------------
    //       GHOSTED NODES NEED INFORMATION ON THEIR COUPLED NODES
    //----------------------------------------------------------------------

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

    // time measurement --- this causes the TimeMonitor tm7 to stop here
    tm7_ref_ = null;


    // time measurement --- start TimeMonitor tm8
    tm8_ref_        = rcp(new TimeMonitor(*timepbcrenumdofs_ ));

    // Create a new DofSet special for Periodic boundary conditions with
    // this type of node coupling

    // create a new dofset specialisation for periodic boundary conditions

    discret_->ReplaceDofSet(rcp(new DRT::PBCDofSet(allcoupledcolnodes_)));

    //--------------------------------------------------
    // redistribute the nodes
    //
    // this contains a call to FillComplete and assigns the same
    // degree of freedom to the matching nodes

    discret_->Redistribute(*newrownodemap,*newcolnodemap);

    // time measurement --- this causes the TimeMonitor tm8 to stop here
    tm8_ref_ = null;

    // throw away old nodegraph
    oldnodegraph = null;

  }

  return;
}// PeriodicBoundaryConditions::RedistributeAndCreateDofCoupling





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

#endif /* CCADISCRET       */
