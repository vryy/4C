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
 | complete coupling map                         (public)    gammi 05/07|
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

    if(discret_->Comm().MyPID()==0)
    {
      cout << "Generate new dofset";
      cout<<endl<<endl;
    }

    // fetch all slaves to the proc of the master
    PutAllSlavesToMastersProc();


    if(discret_->Comm().NumProc()>1)
    {
      if(discret_->Comm().MyPID()==0)
      {
        cout << "\n---------------------------------------------\n";
        cout << "Call METIS \n";
        cout<<endl;
      }

      // eventually call METIS to optimally distribute the nodes --- up to
      // now, a periodic boundary condition might remove all nodes from a
      // proc ...
      BalanceLoadUsingMetis();
    }
    // time measurement --- this causes the TimeMonitor tm0 to stop here
    //                                              (call of destructor)
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

    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      const Epetra_Map* noderowmap = discret_->NodeRowMap();

      int mypid   =discret_->Comm().MyPID()  ;
      int numprocs=discret_->Comm().NumProc();

      int countslave=0;
      for(map<int,vector<int> >::iterator iter=allcoupledcolnodes_->begin();
	  iter!=allcoupledcolnodes_->end();++iter)
      {
	for(vector<int>::iterator viter=iter->second.begin();
	    viter!=iter->second.end();++viter)
	{
	  ++countslave;
	}
      }

      vector<int> my_n_nodes   (numprocs,0);
      vector<int>    n_nodes   (numprocs,0);
      vector<int> my_n_master  (numprocs,0);
      vector<int>    n_master  (numprocs,0);
      vector<int> my_n_slave   (numprocs,0);
      vector<int>    n_slave   (numprocs,0);
      vector<int> my_n_elements(numprocs,0);
      vector<int>    n_elements(numprocs,0);
      vector<int> my_n_ghostele(numprocs,0);
      vector<int>    n_ghostele(numprocs,0);
      vector<int> my_n_dof     (numprocs,0);
      vector<int>    n_dof     (numprocs,0);

      my_n_nodes   [mypid]=noderowmap->NumMyElements ();
      my_n_master  [mypid]=allcoupledcolnodes_->size();
      my_n_slave   [mypid]=countslave;
      my_n_elements[mypid]=discret_->NumMyColElements();
      my_n_ghostele[mypid]=discret_->NumMyColElements()-discret_->NumMyRowElements();
      my_n_dof     [mypid]=dofrowmap->NumMyElements();

      discret_->Comm().SumAll(&my_n_nodes[0]   ,&n_nodes[0]   ,numprocs);
      discret_->Comm().SumAll(&my_n_master[0]  ,&n_master[0]  ,numprocs);
      discret_->Comm().SumAll(&my_n_slave[0]   ,&n_slave[0]   ,numprocs);
      discret_->Comm().SumAll(&my_n_elements[0],&n_elements[0],numprocs);
      discret_->Comm().SumAll(&my_n_ghostele[0],&n_ghostele[0],numprocs);
      discret_->Comm().SumAll(&my_n_dof[0]     ,&n_dof[0]     ,numprocs);

      if(discret_->Comm().MyPID()==0)
      {
	printf("   +-----+---------------+--------------+-------------+-----------------+----------------+-----------------+\n");
	printf("   | PID |    n_nodes    |   n_master   |   n_slave   |    n_elements   |   n_ghostele   |      n_dof      |\n");
	printf("   +-----+---------------+--------------+-------------+-----------------+----------------+-----------------+\n");
	for(int npid=0;npid<numprocs;++npid)
	{
	  printf("   | %3d | %13d | %12d | %11d | %15d | %14d | %15d |\n",npid,n_nodes[npid],n_master[npid],n_slave[npid],n_elements[npid],n_ghostele[npid],n_dof[npid]);
	  printf("   +-----+---------------+--------------+-------------+-----------------+----------------+-----------------+\n");
	}
	cout << endl <<endl;
      }
    }
  } // end if numpbcpairs_>0
  return;

}// UpdateDofsForPeriodicBoundaryConditions()


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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |                                                                      |
 | generate master->slave connectivities                                |
 | o allcoupledrownodes_                                                |
 | o allcoupledcolnodes_ (including ghosted master/slave nodes)         |
 |                                                                      |
 | send slave nodes to master proc.                                     |
 |                                                                      |
 | Generate a new dofset in which slaves do not have their own dofs     |
 | anymore; they just point to the master's dofs                        |
 |                                                           gammi 11/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void PeriodicBoundaryConditions::PutAllSlavesToMastersProc()
{
  // clear old data
  allcoupledrownodes_=rcp(new map<int,vector<int> >);
  allcoupledcolnodes_=rcp(new map<int,vector<int> >);

  // map from global masternodeids (on this proc) to global slavenodeids
  // for a single condition
  map<int,vector<int> > midtosid;

  // pointers to master and slave condition
  DRT::Condition* mastercond=NULL;
  DRT::Condition* slavecond =NULL;

  // global master node Ids and global slave node Ids
  vector <int> masternodeids;
  vector <int> slavenodeids;

  //----------------------------------------------------------------------
  //                     LOOP PERIODIC DIRECTIONS
  //----------------------------------------------------------------------

  vector<string> planes;
  planes.push_back("xy");
  planes.push_back("xz");
  planes.push_back("yz");
  planes.push_back("xyz");

  // the id of the plane --- will be counted in the loop....
  int num=0;

  // loop over periodic directions/planes
  for(vector<string>::iterator thisplane=planes.begin();
      thisplane!=planes.end();
      ++thisplane)
  {
    // loop over all three layers (we allow three layers since
    // the code should be able to deal with up to cubic splines
    // which couple three layers of nodes)
    for(int nlayer=0;nlayer<3;++nlayer)
    {
      // master and slave sets for this periodic direction
      std::set<int> masterset;
      std::set<int> slaveset;
      // possible angles of rotation for slave plane for each pbc pair
      vector<double> rotangles(numpbcpairs_);

      // absolute node matching tolerance for octree
      double abs_tol=0.0;

      // a toggle to indicate whether tolerance for octree was already set,
      // if so check if all values are equal
      bool tol_set=false;

      //----------------------------------------------------
      // in the following, we loop all periodic boundary
      // conditions which have the prescribed periodic
      // direction.
      // For every Master condition, we add the nodes into
      // the set of all masternodeids for periodic boundary
      // conditions with this homogeneous direction.
      // The same is done for the slave conditions.

      // loop pairs of periodic boundary conditions
      for (int pbcid=0;pbcid<numpbcpairs_;++pbcid)
      {

        //--------------------------------------------------
        // get master and slave condition pair with id pbcid

        for (unsigned numcond=0;numcond<mysurfpbcs_.size();++numcond)
        {
          const vector<int>* myid
            = mysurfpbcs_[numcond]->Get<vector<int> >("Id of periodic boundary condition");
          const vector<int>* mylayer
            = mysurfpbcs_[numcond]->Get<vector<int> >("Layer of periodic boundary condition");
          // yes, I am the condition with id pbcid and in the desired layer
          if (myid[0][0] == pbcid && (mylayer[0][0]+1) == nlayer)
          {
            const string* mymasterslavetoggle
              = mysurfpbcs_[numcond]->Get<string>("Is slave periodic boundary condition");

            if(*mymasterslavetoggle=="Master")
            {
              mastercond =mysurfpbcs_[numcond];

              //--------------------------------------------------
              // check whether this periodic boundary condition belongs
              // to thisplane

              const string* dofsforpbcplanename
                =
                mastercond->Get<string>("degrees of freedom for the pbc plane");

              if(*dofsforpbcplanename == *thisplane)
              {
                // add all master nodes to masterset

                //--------------------------------------------------
                // get global master node Ids
                const vector <int>* masteridstoadd;

                masteridstoadd = mastercond->Nodes();

                for(vector<int>::const_iterator idtoadd =(*masteridstoadd).begin();
                    idtoadd!=(*masteridstoadd).end();
                    ++idtoadd)
                {
                  masterset.insert(*idtoadd);
                }

                // check for angle of rotation (has to be zero for master plane)
                const double angle = mastercond->GetDouble("Angle of rotation");
                if (abs(angle) > EPS13)
                  dserror("Angle is not zero for master plane: %f",angle);
              }
            }
            else if (*mymasterslavetoggle=="Slave")
            {
              slavecond =mysurfpbcs_[numcond];

              //--------------------------------------------------
              // check whether this periodic boundary condition belongs
              // to thisplane
              const string* dofsforpbcplanename = slavecond->Get<string>("degrees of freedom for the pbc plane");

              if(*dofsforpbcplanename == *thisplane)
              {
                // add all slave nodes to slaveset

                //--------------------------------------------------
                // get global slave node Ids
                const vector <int>* slaveidstoadd;

                slaveidstoadd = slavecond->Nodes();

                for(vector<int>::const_iterator idtoadd =(*slaveidstoadd).begin();
                    idtoadd!=(*slaveidstoadd).end();
                    ++idtoadd)
                {
                  slaveset.insert(*idtoadd);
                }

                // check for angle of rotation of slave plane and store it
                const double angle = slavecond->GetDouble("Angle of rotation");
                if (abs(angle)> EPS13)
                {
                  if ((*thisplane != "xz")&&(*thisplane != "yz"))
                    dserror("Rotation of slave plane only implemented for xz and yz planes");
                  else
                  {
                    rotangles[pbcid] = angle*PI/180.0;  //convert from DEG to RAD!
                    if (pbcid > 0)
                    {
                      if (rotangles[pbcid] != rotangles[pbcid-1])
                        dserror("Angle has to be the same for all pairs in pbc");
                    }
                  }
                }
              }
            }
            else
            {
              dserror("pbc is neither master nor slave");
            }

            // set tolerance for octree
            const double tol = (mysurfpbcs_[numcond])->GetDouble("Tolerance for nodematching in octree");
            
            if(!tol_set)
            {
              abs_tol = tol;

              tol_set=true;
            }
            else
            {
              if(fabs(abs_tol-tol)>1e-5)
              {
                dserror("none matching tolerances %12.5e neq %12.5e for nodmatching octree. All values in direction %s have to match\n",abs_tol,tol,(*thisplane).c_str());
              }
            }
          } // end if i am the right condition in the right layer
        } // end loop over conditions
      } // end loop pairs of periodic boundary conditions


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
      
      // we transform the three strings "xy", "xz", "yz" into integer
      // values dofsforpbcplanename
      vector<int> dofsforpbcplane(2);

      // this is a char-operation:
      //
      //       x -> 0
      //       y -> 1
      //       z -> 2
      //
      // it is based on the fact that the letters x, y and z are
      // consecutive in the ASCII table --- 'x' is the ASCII
      // calue of x ....

      if (*thisplane == "xyz")
      {
        // nodes in exact the same position are coupled
        dofsforpbcplane.clear();
      }
      else
      {
        dofsforpbcplane[0] = thisplane->c_str()[0] - 'x';
        dofsforpbcplane[1] = thisplane->c_str()[1] - 'x';
      }
      //--------------------------------------------------
      // we just write the sets into vectors
      (masternodeids).clear();
      (slavenodeids ).clear();

      for(std::set<int>::iterator appendednode = masterset.begin();
          appendednode != masterset.end();
          ++appendednode)
      {
        masternodeids.push_back(*appendednode);
      }

      for(std::set<int>::iterator appendednode = slaveset.begin();
          appendednode != slaveset.end();
          ++appendednode)
      {
        slavenodeids.push_back(*appendednode);
      }

      //----------------------------------------------------------------------
      //      CONSTRUCT NODE MATCHING BETWEEN MASTER AND SLAVE NODES
      //                        FOR THIS DIRECTION
      //----------------------------------------------------------------------

      // clear map from global masternodeids (on this proc) to global
      // slavenodeids --- it belongs to this master slave pair!!!
      midtosid.clear();

      if (discret_->Comm().MyPID() == 0)
      {
        cout << " creating layer " << nlayer << " of midtosid-map in " << *thisplane << " direction ... ";
        fflush(stdout);
      }

      // get map master on this proc -> slave on some proc
      CreateNodeCouplingForSinglePBC(
        midtosid,
        masternodeids,
        slavenodeids ,
        dofsforpbcplane,
        rotangles[0],
        abs_tol);
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
      AddConnectivity(midtosid,num);

      // time measurement --- this causes the TimeMonitor tm4 to stop here
      tm4_ref_ = null;

      if (discret_->Comm().MyPID() == 0)
      {
        cout << " done\n";
        fflush(stdout);
      }

      ++num;
    } // end loop over layers
  } // end loop over planes

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

  return;
} // PutAllSlavesToMastersProc()

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Couple nodes for specific pair of periodic boundary conditions       |
 |                                                           gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void PeriodicBoundaryConditions::CreateNodeCouplingForSinglePBC(
  map<int,vector<int> > &midtosid,
  const vector <int>     masternodeids,
  const vector <int>     slavenodeids,
  const vector <int>     dofsforpbcplane,
  const double           rotangle,
  const double           abstol
  )
{


  // these are just parameter definitions for the octree search algorithm
  double tol            = abstol;
  int    maxnodeperleaf = 250;

  //----------------------------------------------------------------------
  //                   BUILD PROCESSOR LOCAL OCTREE
  //----------------------------------------------------------------------
  // time measurement --- start TimeMonitor tm2
  tm2_ref_        = rcp(new TimeMonitor(*timepbcmidoct_ ));

  // build processor local octree
  DRT::UTILS::NodeMatchingOctree nodematchingoctree(
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
        rotangle,
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
 |                                                           gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void PeriodicBoundaryConditions::AddConnectivity(
  map<int,vector<int> > &midtosid,
  const int              pbcid
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
    map<int,vector<int> >::iterator iter;
    int  masterid;
    int  slaveid;
    bool alreadyinlist;

    for( iter = midtosid.begin(); iter != midtosid.end(); ++iter )
    {
      // get id of masternode and slavenode
      masterid  = iter->first;

      vector<int>::iterator i;
      for(i=(iter->second).begin();i!=(iter->second).end();++i)
      {
        slaveid = *i;
        if (slaveid == masterid)
          dserror("Node %d is master AND slave node of periodic boundary condition", masterid);

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
            const string* mymasterslavetoggle
              = thiscond[numcond]->Get<string>("Is slave periodic boundary condition");

            if(*mymasterslavetoggle=="Master")
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
 | current distribution of nodes                             gammi 05/07|
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
            const string* mymasterslavetoggle
              = thiscond[numcond]->Get<string>("Is slave periodic boundary condition");

            if(*mymasterslavetoggle=="Slave")
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

    discret_->ReplaceDofSet(rcp(new PBCDofSet(allcoupledcolnodes_)));

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
 | o adjust weights of slavenodes                                       |
 |   they need a small weight since they do not contribute any dofs     |
 |   to the linear system                                               |
 |                                                                      |
 | o compute connectivity                                               |
 |   iterate all elements on this proc including ghosted ones. Include  |
 |   connections between master and slave nodes                         |
 |                                                                      |
 | o set weights of edges between master/slave pairs to a high value    |
 |   in order to keep both on the same proc when redistributing         |
 |                                                                      |
 | o gather all data to proc 1, do partitioning using METIS             |
 |                                                                      |
 | o redistribute nodes without assigning dofs                          |
 |                                                                      |
 | o repair master/slave distribution, finally assign dofs              |
 |                                                           gammi 11/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void PeriodicBoundaryConditions::BalanceLoadUsingMetis()
{

  if (discret_->Comm().NumProc()>1)
  {
    const Epetra_Map* noderowmap = discret_->NodeRowMap();

    // weights for graph partition
    Epetra_Vector weights(*noderowmap,false);
    weights.PutScalar(10.0);

    // ----------------------------------------
    // loop masternodes to adjust weights of slavenodes
    // they need a small weight since they do not contribute any dofs
    // to the linear system
    {
      map<int,vector<int> >::iterator masterslavepair;

      for(masterslavepair =allcoupledcolnodes_->begin();
          masterslavepair!=allcoupledcolnodes_->end()  ;
          ++masterslavepair)
      {
        // get masternode
        DRT::Node*  master = discret_->gNode(masterslavepair->first);

        if(master->Owner()!=discret_->Comm().MyPID())
        {
          continue;
        }

        // loop slavenodes associated with master
        for(vector<int>::iterator iter=masterslavepair->second.begin();
            iter!=masterslavepair->second.end();++iter)
        {
          double initval=1.0;
          int gid       =*iter;

          weights.ReplaceGlobalValues(1,&initval,&gid);
        }
      }
    }

    // allocate graph
    RefCountPtr<Epetra_CrsGraph> nodegraph = rcp(new Epetra_CrsGraph(Copy,*noderowmap,108,false));

    // -------------------------------------------------------------
    // iterate all elements on this proc including ghosted ones
    // compute connectivity

    // standard part without master<->slave coupling
    // Note:
    // if a proc stores the appropiate ghosted elements, the resulting
    // graph will be the correct and complete graph of the distributed
    // discretization even if nodes are not ghosted.

    for (int nele=0;nele<discret_->NumMyColElements();++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lColElement(nele);

      // get its nodes and nodeids
      const int  nnode   = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      for (int row=0; row<nnode; ++row)
      {
        const int rownode = nodeids[row];

        // insert into line of graph only when this proc owns the node
        if (!noderowmap->MyGID(rownode)) continue;

        // insert all neighbours from element in the graph
        for (int col=0; col<nnode; ++col)
        {
          int colnode = nodeids[col];
          int err = nodegraph->InsertGlobalIndices(rownode,1,&colnode);
          if (err<0) dserror("nodegraph->InsertGlobalIndices returned err=%d",err);
        }
      }
    }

    // -------------------------------------------------------------
    // additional coupling between master and slave
    // we do not only connect master and slave nodes but if a master/slave
    // is connected to a master/slave, we connect the corresponding slaves/master
    // as well

    for (int nele=0;nele<discret_->NumMyColElements();++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lColElement(nele);

      // get its nodes and nodeids
      const int  nnode   = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      for (int row=0; row<nnode; ++row)
      {
        const int rownode = nodeids[row];

        // insert into line of graph only when this proc owns the node
        if (!noderowmap->MyGID(rownode)) continue;

        map<int,vector<int> >::iterator masterslavepair = allcoupledcolnodes_->find(rownode);
        if(masterslavepair!=allcoupledcolnodes_->end())
        {
          // get all masternodes of this element
          for (int col=0; col<nnode; ++col)
          {
            int colnode = nodeids[col];

            map<int,vector<int> >::iterator othermasterslavepair = allcoupledcolnodes_->find(colnode);
            if(othermasterslavepair!=allcoupledcolnodes_->end())
            {
              // add connection to all slaves

              for(vector<int>::iterator iter=othermasterslavepair->second.begin();
                  iter!=othermasterslavepair->second.end();++iter)
              {
                int othermastersslaveindex = *iter;
                int masterindex            = rownode;
                int err = nodegraph->InsertGlobalIndices(rownode,1,&othermastersslaveindex);
                if (err<0) dserror("nodegraph->InsertGlobalIndices returned err=%d",err);

                if (noderowmap->MyGID(*iter))
                {
                  err = nodegraph->InsertGlobalIndices(*iter,1,&masterindex);
                  if (err<0) dserror("nodegraph->InsertGlobalIndices returned err=%d",err);
                }
              }
            }
          }
        }
      }
    }

    // finalize construction of initial graph
    int err = nodegraph->FillComplete();
    if (err) dserror("graph->FillComplete returned %d",err);

    const int myrank   = nodegraph->Comm().MyPID();
    const int numproc  = nodegraph->Comm().NumProc();

    if (numproc>1)
    {
      // proc that will do the serial partitioning
      // the graph is collapsed to this proc
      // Normally this would be proc 0 but 0 always has so much to do.... ;-)
      int workrank=1;

      // get rowmap of the graph
      const Epetra_BlockMap& tmp = nodegraph->RowMap();
      Epetra_Map rowmap(tmp.NumGlobalElements(),tmp.NumMyElements(),
                        tmp.MyGlobalElements(),0,nodegraph->Comm());

      // -------------------------------------------------------------
      // build a target map that stores everything on proc workrank
      // We have arbirtary gids here and we do not tell metis about
      // them. So we have to keep rowrecv until the redistributed map is
      // build.

      // rowrecv is a fully redundant vector (size of number of nodes)
      vector<int> rowrecv(rowmap.NumGlobalElements());

      // after AllreduceEMap rowrecv contains
      //
      // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
      // * | | .... | | * | | .... | | * ..........  * | | .... | | *
      // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
      //   gids stored     gids stored                  gids stored
      //  on first proc  on second proc                 on last proc
      //
      // the ordering of the gids on the procs is arbitrary (as copied
      // from the map)
      LINALG::AllreduceEMap(rowrecv, rowmap);

      // construct an epetra map from the list of gids
      Epetra_Map tmap(rowmap.NumGlobalElements(),
                      // if ..........    then ............... else
                      (myrank == workrank) ? (int)rowrecv.size() : 0,
                      &rowrecv[0],
                      0,
                      rowmap.Comm());

      // export the graph to tmap
      Epetra_CrsGraph tgraph(Copy,tmap,108,false);
      Epetra_Export exporter(rowmap,tmap);
      {
        int err = tgraph.Export(*nodegraph,exporter,Add);
        if (err<0) dserror("Graph export returned err=%d",err);
      }
      tgraph.FillComplete();
      tgraph.OptimizeStorage();

      // export the weights to tmap
      Epetra_Vector tweights(tmap,false);
      err = tweights.Export(weights,exporter,Insert);
      if (err<0) dserror("Vector export returned err=%d",err);

      // metis requests indexes. So we need a reverse lookup from gids
      // to indexes.
      map<int,int> idxmap;
      // xadj points from index i to the index of the
      // first adjacent node
      vector<int> xadj  (rowmap.NumGlobalElements()+1);
      // a list of adjacent nodes, adressed using xadj
      vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound

      // This is a vector of size n that upon successful completion stores the partition vector of the graph
      vector<int> part(tmap.NumMyElements());

      // construct reverse lookup for all procs
      for (unsigned i=0; i<rowrecv.size(); ++i)
      {
        idxmap[rowrecv[i]] = i;
      }

      if (myrank==workrank)
      {
        // ----------------------------------------

        // rowrecv(i)       rowrecv(i+1)                      node gids
        //     ^                 ^
        //     |                 |
        //     | idxmap          | idxmap
        //     |                 |
        //     v                 v
        //     i                i+1                       equivalent indices
        //     |                 |
        //     | xadj            | xadj
        //     |                 |
        //     v                 v
        //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
        //    | | | | | | | | | | | ............... | | |      adjncy
        //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
        //
        //    |       i's       |    (i+1)'s
        //    |    neighbours   |   neighbours           (numbered by equivalent indices)
        //

        int count=0;
        xadj[0] = 0;
        for (int row=0; row<tgraph.NumMyRows(); ++row)
        {
          int grid = tgraph.RowMap().GID(row);
          int numindices;
          int* lindices;
          int err = tgraph.ExtractMyRowView(row,numindices,lindices);
          if (err) dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d",err);

          for (int col=0; col<numindices; ++col)
          {
            int gcid = tgraph.ColMap().GID(lindices[col]);
            if (gcid==grid) continue;
            adjncy[count] = idxmap[gcid];
            ++count;
          }
          xadj[row+1] = count;
        }
      }

      // broadcast xadj
      tmap.Comm().Broadcast(&xadj[0],xadj.size(),workrank);

      // broadcast adjacence (required for edge weights)
      int adjncysize = (int)adjncy.size();
      tmap.Comm().Broadcast(&adjncysize,1,workrank);
      adjncy.resize(adjncysize);
      tmap.Comm().Broadcast(&adjncy[0],adjncysize,workrank);

      // -------------------------------------------------------------
      // set a fully redundant vector of weights for edges
      vector<int> ladjwgt(adjncy.size(),0);
      vector<int>  adjwgt(adjncy.size(),0);

      for(vector<int>::iterator iter =ladjwgt.begin();
          iter!=ladjwgt.end();
          ++iter)
      {
        *iter=0;
      }

      // loop all master nodes on this proc
      map<int,vector<int> >::iterator masterslavepair;

      for(masterslavepair =allcoupledcolnodes_->begin();
          masterslavepair!=allcoupledcolnodes_->end()  ;
          ++masterslavepair)
      {
        // get masternode
        DRT::Node*  master = discret_->gNode(masterslavepair->first);

        if(master->Owner()!=myrank)
        {
          continue;
        }

        map<int,int>::iterator paul=idxmap.find(master->Id());
        if (paul == idxmap.end())
        {
          dserror("master not in reverse lookup");
        }

        // inverse lookup
        int masterindex=idxmap[master->Id()];

        // loop slavenodes
        for(vector<int>::iterator iter=masterslavepair->second.begin();
            iter!=masterslavepair->second.end();++iter)
        {
          DRT::Node*  slave = discret_->gNode(*iter);

          if(slave->Owner()!=myrank)
          {
            dserror("own master but not slave\n");
          }

          int slaveindex=idxmap[slave->Id()];

          map<int,int>::iterator foo=idxmap.find(slave->Id());
          if (foo == idxmap.end())
          {
            dserror("slave not in reverse lookup");
          }

          // -------------------------------------------------------------
          // connections between master and slavenodes are very strong
          // we do not want to partition between master and slave nodes
          for(int j = xadj[masterindex];j<xadj[masterindex+1];++j)
          {
            if(adjncy[j] == slaveindex)
            {
              ladjwgt[j]=100;
            }
          }

          for(int j = xadj[slaveindex];j<xadj[slaveindex+1];++j)
          {
            if(adjncy[j] == masterindex)
            {
              ladjwgt[j]=100;
            }
          }
        }
      }

      // do communication to aquire edge weight information from all procs
      tmap.Comm().SumAll(&ladjwgt[0], &adjwgt[0], adjwgt.size());

      // the standard edge weight is one
      for(vector<int>::iterator iter =adjwgt.begin();
          iter!=adjwgt.end();
          ++iter)
      {
        if(*iter==0)
        *iter=1;
      }

      // the reverse lookup is not required anymore
      idxmap.clear();

      // -------------------------------------------------------------
      // do partitioning using metis on workrank
      if (myrank==workrank)
      {
        // the vertex weights
        vector<int> vwgt(tweights.MyLength());
        for (int i=0; i<tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

        // 0 No weights (vwgts and adjwgt are NULL)
        // 1 Weights on the edges only (vwgts = NULL)
        // 2 Weights on the vertices only (adjwgt = NULL)
        // 3 Weights both on vertices and edges.
        int wgtflag=3;
        // 0 C-style numbering is assumed that starts from 0
        // 1 Fortran-style numbering is assumed that starts from 1
        int numflag=0;
        // The number of parts to partition the graph.
        int npart=numproc;
        // This is an array of 5 integers that is used to pass parameters for the various phases of the algorithm.
        // If options[0]=0 then default values are used. If options[0]=1, then the remaining four elements of
        // options are interpreted as follows:
        // options[1]    Determines matching type. Possible values are:
        //               1 Random Matching (RM)
        //               2 Heavy-Edge Matching (HEM)
        //               3 Sorted Heavy-Edge Matching (SHEM) (Default)
        //               Experiments has shown that both HEM and SHEM perform quite well.
        // options[2]    Determines the algorithm used during initial partitioning. Possible values are:
        //               1 Region Growing (Default)
        // options[3]    Determines the algorithm used for re%Gï¬%@nement. Possible values are:
        //               1 Early-Exit Boundary FM re%Gï¬%@nement (Default)
        // options[4]    Used for debugging purposes. Always set it to 0 (Default).
        int options[5] = { 0,3,1,1,0 };
        // Upon successful completion, this variable stores the number of edges that are cut by the partition.
        int edgecut=0;
        // The number of vertices in the graph.
        int nummyele = tmap.NumMyElements();

        cout << "proc " <<  myrank << " repartition graph using metis\n";
        if (numproc<8) // better for smaller no. of partitions
        {
#ifdef PARALLEL
          METIS_PartGraphRecursive(&nummyele,
                                   &xadj[0],
                                   &adjncy[0],
                                   &vwgt[0],
                                   &adjwgt[0],
                                   &wgtflag,
                                   &numflag,
                                   &npart,
                                   options,
                                   &edgecut,
                                   &part[0]);

          cout << "METIS_PartGraphRecursive produced edgecut of " << edgecut << "\n";
          fflush(stdout);
#endif
        }
        else
        {
#ifdef PARALLEL
          METIS_PartGraphKway(&nummyele,
                              &xadj[0],
                              &adjncy[0],
                              &vwgt[0],
                              &adjwgt[0],
                              &wgtflag,
                              &numflag,
                              &npart,
                              options,
                              &edgecut,
                              &part[0]);
#endif
        }


      } // if (myrank==workrank)

      // broadcast partitioning result
      int size = tmap.NumMyElements();
      tmap.Comm().Broadcast(&size,1,workrank);
      part.resize(size);
      tmap.Comm().Broadcast(&part[0],size,workrank);

      // loop part and count no. of nodes belonging to me
      // (we reuse part to save on memory)
      int count=0;
      for (int i=0; i<size; ++i)
      if (part[i]==myrank)
      {
        part[count] = rowrecv[i];
        ++count;
      }

      // rowrecv is done
      rowrecv.clear();

      // create map with new layout
      Epetra_Map newmap(size,count,&part[0],0,nodegraph->Comm());

      // create the new graph and export to it
      RefCountPtr<Epetra_CrsGraph> newnodegraph;

      newnodegraph = rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
      Epetra_Export exporter2(nodegraph->RowMap(),newmap);
      err = newnodegraph->Export(*nodegraph,exporter2,Add);
      if (err<0) dserror("Graph export returned err=%d",err);
      newnodegraph->FillComplete();
      newnodegraph->OptimizeStorage();

      // the rowmap will become the new distribution of nodes
      const Epetra_BlockMap rntmp = newnodegraph->RowMap();
      Epetra_Map newnoderowmap(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,discret_->Comm());

      // the column map will become the new ghosted distribution of nodes
      const Epetra_BlockMap Mcntmp = newnodegraph->ColMap();
      Epetra_Map newnodecolmap(-1,Mcntmp.NumMyElements(),Mcntmp.MyGlobalElements(),0,discret_->Comm());
      // do the redistribution without assigning dofs
      discret_->Redistribute(newnoderowmap,newnodecolmap,false,true,true);


      if(discret_->Comm().MyPID()==0)
      {
        cout << "---------------------------------------------\n";
        cout << "Repair Master->Slave connection, generate final dofset";
        cout<<endl<<endl;
      }

      // assign the new dofs, make absolutely sure that we always
      // have all slaves to a master
      // the finite edge weights are not a 100% warranty for that...
      PutAllSlavesToMastersProc();
    }
  }

  return;
}// BalanceLoadUsingMetis

#endif /* CCADISCRET       */
