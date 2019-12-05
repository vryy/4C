/*---------------------------------------------------------------------*/
/*! \file

\brief  update dofrowmap and dofsets for periodic boundary conditions

        o effects distribution of master and slave nodes. Master and
          slave are owned by one proc afterwards.
        o effects ghosting of nodes. Master and slave nodes are ghosted
          in pairs.
        o passes list of coupled nodes to the dofset

\level 0


\maintainer  Martin Kronbichler

*/
/*---------------------------------------------------------------------*/


#include "drt_periodicbc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_matchingoctree.H"
#include "../drt_lib/drt_dofset_pbc.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
PeriodicBoundaryConditions::PeriodicBoundaryConditions(
    Teuchos::RCP<DRT::Discretization> actdis, bool verbose)
    : discret_(actdis), verbose_(verbose), pbcdofset_(Teuchos::null)
{
  // get periodic surface boundary conditions
  discret_->GetCondition("SurfacePeriodic", mysurfpbcs_);

  if (mysurfpbcs_.empty())
  {
    discret_->GetCondition("LinePeriodic", mysurfpbcs_);
  }

  // set number of pairs of periodic boundary conditions
  numpbcpairs_ = mysurfpbcs_.size() / 2;

  // create map that will be connecting master to slave nodes owned by
  // this proc
  //       master node -> list of his slave node(s)
  allcoupledrownodes_ = Teuchos::rcp(new std::map<int, std::vector<int>>);

  // create map that will be connecting master to slave nodes owned or
  // ghosted by this proc
  //       master node -> list of his slave node(s)
  allcoupledcolnodes_ = Teuchos::rcp(new std::map<int, std::vector<int>>);


  if (numpbcpairs_ > 0)
  {
    // -------------------------------------------------------------------
    // create timers and time monitor
    // -------------------------------------------------------------------
    timepbctot_ = Teuchos::TimeMonitor::getNewTimer("0) pbc routine total");
    timepbcmidtosid_ = Teuchos::TimeMonitor::getNewTimer("1)   +create midtosid maps");
    timepbcmidoct_ = Teuchos::TimeMonitor::getNewTimer("2)      +build local octrees");
    timepbcmidmatch_ =
        Teuchos::TimeMonitor::getNewTimer("3)      +search closest nodes in octrees on all procs");
    timepbcaddcon_ =
        Teuchos::TimeMonitor::getNewTimer("4)   +add connectivity to previous conditions");
    timepbcreddis_ = Teuchos::TimeMonitor::getNewTimer("5)   +Redistribute the nodes");
    timepbcmakeghostmap_ =
        Teuchos::TimeMonitor::getNewTimer("6)      +build rowmap and temporary colmap");
    timepbcghost_ = Teuchos::TimeMonitor::getNewTimer("7)      +repair ghosting");
    timepbcrenumdofs_ = Teuchos::TimeMonitor::getNewTimer("8)      +call discret->Redistribute");
  }

  return;

}  // PeriodicBoundaryConditions(Teuchos::RCP<DRT::Discretization> actdis)


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
  if (numpbcpairs_ > 0)
  {
    // time measurement --- start TimeMonitor tm0
    tm0_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbctot_));

    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << "Generate new dofset for discretization " << discret_->Name();
      std::cout << std::endl << std::endl;
    }

    // fetch all slaves to the proc of the master
    PutAllSlavesToMastersProc();


    if (discret_->Comm().NumProc() > 1)
    {
      if (discret_->Comm().MyPID() == 0 && verbose_)
      {
        std::cout << "\n---------------------------------------------\n";
        std::cout << "Call METIS \n";
        std::cout << std::endl;
      }

      // eventually call METIS to optimally distribute the nodes --- up to
      // now, a periodic boundary condition might remove all nodes from a
      // proc ...
      BalanceLoadUsingMetis();
    }
    // time measurement --- this causes the TimeMonitor tm0 to stop here
    //                                              (call of destructor)
    tm0_ref_ = Teuchos::null;

    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << std::endl << std::endl;
    }

    if (verbose_)
    {
      Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm =
          COMM_UTILS::toTeuchosComm<int>(discret_->Comm());
      Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
    }

    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << std::endl << std::endl;
    }

    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      const Epetra_Map* noderowmap = discret_->NodeRowMap();

      int mypid = discret_->Comm().MyPID();
      int numprocs = discret_->Comm().NumProc();

      int countslave = 0;
      for (std::map<int, std::vector<int>>::iterator iter = allcoupledcolnodes_->begin();
           iter != allcoupledcolnodes_->end(); ++iter)
      {
        for (std::vector<int>::iterator viter = iter->second.begin(); viter != iter->second.end();
             ++viter)
        {
          ++countslave;
        }
      }

      std::vector<int> my_n_nodes(numprocs, 0);
      std::vector<int> n_nodes(numprocs, 0);
      std::vector<int> my_n_master(numprocs, 0);
      std::vector<int> n_master(numprocs, 0);
      std::vector<int> my_n_slave(numprocs, 0);
      std::vector<int> n_slave(numprocs, 0);
      std::vector<int> my_n_elements(numprocs, 0);
      std::vector<int> n_elements(numprocs, 0);
      std::vector<int> my_n_ghostele(numprocs, 0);
      std::vector<int> n_ghostele(numprocs, 0);
      std::vector<int> my_n_dof(numprocs, 0);
      std::vector<int> n_dof(numprocs, 0);

      my_n_nodes[mypid] = noderowmap->NumMyElements();
      my_n_master[mypid] = allcoupledcolnodes_->size();
      my_n_slave[mypid] = countslave;
      my_n_elements[mypid] = discret_->NumMyColElements();
      my_n_ghostele[mypid] = discret_->NumMyColElements() - discret_->NumMyRowElements();
      my_n_dof[mypid] = dofrowmap->NumMyElements();

      discret_->Comm().SumAll(&my_n_nodes[0], &n_nodes[0], numprocs);
      discret_->Comm().SumAll(&my_n_master[0], &n_master[0], numprocs);
      discret_->Comm().SumAll(&my_n_slave[0], &n_slave[0], numprocs);
      discret_->Comm().SumAll(&my_n_elements[0], &n_elements[0], numprocs);
      discret_->Comm().SumAll(&my_n_ghostele[0], &n_ghostele[0], numprocs);
      discret_->Comm().SumAll(&my_n_dof[0], &n_dof[0], numprocs);

      if (discret_->Comm().MyPID() == 0 && verbose_)
      {
        printf(
            "   "
            "+-----+---------------+--------------+-------------+-----------------+----------------"
            "+-----------------+\n");
        printf(
            "   | PID |    n_nodes    |   n_master   |   n_slave   |    n_elements   |   "
            "n_ghostele   |      n_dof      |\n");
        printf(
            "   "
            "+-----+---------------+--------------+-------------+-----------------+----------------"
            "+-----------------+\n");
        for (int npid = 0; npid < numprocs; ++npid)
        {
          printf("   | %3d | %13d | %12d | %11d | %15d | %14d | %15d |\n", npid, n_nodes[npid],
              n_master[npid], n_slave[npid], n_elements[npid], n_ghostele[npid], n_dof[npid]);
          printf(
              "   "
              "+-----+---------------+--------------+-------------+-----------------+--------------"
              "--+-----------------+\n");
        }
        std::cout << std::endl << std::endl;
      }
    }
  }  // end if numpbcpairs_>0
  return;

}  // UpdateDofsForPeriodicBoundaryConditions()


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
}  // ~PeriodicBoundaryConditions()


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
  if (numpbcpairs_ > 0)
  {
    // clear old data
    allcoupledrownodes_ = Teuchos::rcp(new std::map<int, std::vector<int>>);
    allcoupledcolnodes_ = Teuchos::rcp(new std::map<int, std::vector<int>>);

    // map from global masternodeids (on this proc) to global slavenodeids
    // for a single condition
    std::map<int, std::vector<int>> midtosid;

    // pointers to master and slave condition
    DRT::Condition* mastercond = NULL;
    DRT::Condition* slavecond = NULL;

    // global master node Ids and global slave node Ids
    std::vector<int> masternodeids;
    std::vector<int> slavenodeids;

    //----------------------------------------------------------------------
    //                     LOOP PERIODIC DIRECTIONS
    //----------------------------------------------------------------------

    std::vector<std::string> planes;
    planes.push_back("xy");
    planes.push_back("xz");
    planes.push_back("yz");
    planes.push_back("xyz");

    // the id of the plane --- will be counted in the loop....
    int num = 0;

    // loop over periodic directions/planes
    for (std::vector<std::string>::iterator thisplane = planes.begin(); thisplane != planes.end();
         ++thisplane)
    {
      // loop over all three layers (we allow three layers since
      // the code should be able to deal with up to cubic splines
      // which couple three layers of nodes)
      for (int nlayer = 0; nlayer < 3; ++nlayer)
      {
        //----------------------------------------------------
        // in the following, we loop all periodic boundary
        // conditions which have the prescribed periodic
        // direction.
        // For every Master condition, we add the nodes into
        // the set of all masternodeids for periodic boundary
        // conditions with this homogeneous direction.
        // The same is done for the slave conditions.
        // loop pairs of periodic boundary conditions
        for (int pbcid = 0; pbcid < numpbcpairs_; ++pbcid)
        {
          // master and slave sets for this periodic direction
          std::set<int> masterset;
          std::set<int> slaveset;
          // possible angles of rotation for slave plane for each pbc pair
          std::vector<double> rotangles(numpbcpairs_);

          // absolute node matching tolerance for octree
          double abs_tol = 0.0;

          // a toggle to indicate whether tolerance for octree was already set,
          // if so check if all values are equal
          bool tol_set = false;

          //--------------------------------------------------
          // get master and slave condition pair with id pbcid

          for (unsigned numcond = 0; numcond < mysurfpbcs_.size(); ++numcond)
          {
            const std::vector<int>* myid =
                mysurfpbcs_[numcond]->Get<std::vector<int>>("Id of periodic boundary condition");
            const std::vector<int>* mylayer =
                mysurfpbcs_[numcond]->Get<std::vector<int>>("Layer of periodic boundary condition");
            // yes, I am the condition with id pbcid and in the desired layer

            if (myid[0][0] == pbcid && (mylayer[0][0] + 1) == nlayer)
            {
              const std::string* mymasterslavetoggle =
                  mysurfpbcs_[numcond]->Get<std::string>("Is slave periodic boundary condition");

              if (*mymasterslavetoggle == "Master")
              {
                mastercond = mysurfpbcs_[numcond];

                //--------------------------------------------------
                // check whether this periodic boundary condition belongs
                // to thisplane

                const std::string* dofsforpbcplanename =
                    mastercond->Get<std::string>("degrees of freedom for the pbc plane");

                if (*dofsforpbcplanename == *thisplane)
                {
                  // add all master nodes to masterset

                  //--------------------------------------------------
                  // get global master node Ids
                  const std::vector<int>* masteridstoadd;

                  masteridstoadd = mastercond->Nodes();

                  for (std::vector<int>::const_iterator idtoadd = (*masteridstoadd).begin();
                       idtoadd != (*masteridstoadd).end(); ++idtoadd)
                  {
                    // we only add row nodes to the set
                    if (discret_->HaveGlobalNode(*idtoadd))
                      if (discret_->gNode(*idtoadd)->Owner() == discret_->Comm().MyPID())
                        masterset.insert(*idtoadd);
                  }

                  // check for angle of rotation (has to be zero for master plane)
                  const double angle = mastercond->GetDouble("Angle of rotation");
                  if (abs(angle) > EPS13) dserror("Angle is not zero for master plane: %f", angle);
                }
              }
              else if (*mymasterslavetoggle == "Slave")
              {
                slavecond = mysurfpbcs_[numcond];

                //--------------------------------------------------
                // check whether this periodic boundary condition belongs
                // to thisplane
                const std::string* dofsforpbcplanename =
                    slavecond->Get<std::string>("degrees of freedom for the pbc plane");

                if (*dofsforpbcplanename == *thisplane)
                {
                  // add all slave nodes to slaveset

                  //--------------------------------------------------
                  // get global slave node Ids
                  const std::vector<int>* slaveidstoadd;

                  slaveidstoadd = slavecond->Nodes();

                  for (std::vector<int>::const_iterator idtoadd = (*slaveidstoadd).begin();
                       idtoadd != (*slaveidstoadd).end(); ++idtoadd)
                  {
                    // we only add row nodes to the set
                    if (discret_->HaveGlobalNode(*idtoadd))
                      if (discret_->gNode(*idtoadd)->Owner() == discret_->Comm().MyPID())
                        slaveset.insert(*idtoadd);
                  }

                  // check for angle of rotation of slave plane and store it
                  const double angle = slavecond->GetDouble("Angle of rotation");
                  if (abs(angle) > EPS13)
                  {
                    if ((*thisplane != "xz") && (*thisplane != "yz"))
                      dserror("Rotation of slave plane only implemented for xz and yz planes");
                    else
                    {
                      rotangles[pbcid] = angle * PI / 180.0;  // convert from DEG to RAD!
                      if (pbcid > 0)
                      {
                        if (rotangles[pbcid] != rotangles[pbcid - 1])
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
              const double tol =
                  (mysurfpbcs_[numcond])->GetDouble("Tolerance for nodematching in octree");

              if (!tol_set)
              {
                abs_tol = tol;

                tol_set = true;
              }
              else
              {
                if (fabs(abs_tol - tol) > 1e-5)
                {
                  dserror(
                      "none matching tolerances %12.5e neq %12.5e for nodmatching octree. All "
                      "values in direction %s have to match\n",
                      abs_tol, tol, (*thisplane).c_str());
                }
              }
            }  // end if i am the right condition in the right layer
          }    // end loop over conditions



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

          // we transform the three std::strings "xy", "xz", "yz" into integer
          // values dofsforpbcplanename
          std::vector<int> dofsforpbcplane(2);

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
          (slavenodeids).clear();

          for (std::set<int>::iterator appendednode = masterset.begin();
               appendednode != masterset.end(); ++appendednode)
          {
            masternodeids.push_back(*appendednode);
          }

          for (std::set<int>::iterator appendednode = slaveset.begin();
               appendednode != slaveset.end(); ++appendednode)
          {
            slavenodeids.push_back(*appendednode);
          }

          //----------------------------------------------------------------------
          //      CONSTRUCT NODE MATCHING BETWEEN MASTER AND SLAVE NODES
          //                        FOR THIS DIRECTION
          //----------------------------------------------------------------------

          // time measurement --- start TimeMonitor tm1
          tm1_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcmidtosid_));

          // clear map from global masternodeids (on this proc) to global
          // slavenodeids --- it belongs to this master slave pair!!!
          midtosid.clear();

          if (discret_->Comm().MyPID() == 0 && verbose_ && pbcid == numpbcpairs_ - 1)
          {
            std::cout << " creating layer " << nlayer << " of midtosid-map in " << *thisplane
                      << " direction ... ";
            fflush(stdout);
          }

          // get map master on this proc -> slave on some proc
          CreateNodeCouplingForSinglePBC(
              midtosid, masternodeids, slavenodeids, dofsforpbcplane, rotangles[0], abs_tol);
          // time measurement --- this causes the TimeMonitor tm1 to stop here
          tm1_ref_ = Teuchos::null;

          if (discret_->Comm().NumProc() == 1)
          {
            if (masternodeids.size() != midtosid.size())
            {
              // before throwing dserror, print helpful information to screen
              for (size_t i = 0; i < masternodeids.size(); i++)
              {
                int mid = masternodeids[i];
                bool found = false;
                std::map<int, std::vector<int>>::iterator curr;
                for (curr = midtosid.begin(); curr != midtosid.end(); ++curr)
                {
                  if (curr->first == mid)
                  {
                    found = true;
                    break;
                  }
                }
                if (not found)
                {
                  const double* x = discret_->gNode(mid)->X();
                  std::cout << "\nmaster node not found in midtosid list: " << mid
                            << "  coord: x=" << x[0] << " y=" << x[1] << " z=" << x[2];
                }
              }
              // now it is time for the dserror
              dserror("have %d masters in midtosid list, %d expected\n", midtosid.size(),
                  masternodeids.size());
            }
          }

          if (discret_->Comm().MyPID() == 0 && verbose_ && pbcid == numpbcpairs_ - 1)
          {
            std::cout << "adding connectivity to previous pbcs ... ";
            fflush(stdout);
          }

          // time measurement --- start TimeMonitor tm4
          tm4_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcaddcon_));

          //----------------------------------------------------------------------
          //      ADD CONNECTIVITY TO CONNECTIVITY OF ALL PREVIOUS PBCS
          //----------------------------------------------------------------------
          // Add the connectivity from this condition to the connectivity
          // of all previously processed periodic boundary conditions.
          // Redistribute the nodes (rownodes+ghosting)
          // Assign the same degrees of freedom to coupled nodes
          AddConnectivity(midtosid, num);

          // time measurement --- this causes the TimeMonitor tm4 to stop here
          tm4_ref_ = Teuchos::null;

          if (discret_->Comm().MyPID() == 0 && verbose_ && pbcid == numpbcpairs_ - 1)
          {
            std::cout << "done.\n";
            fflush(stdout);
          }
        }  // end loop pairs of periodic boundary conditions
        ++num;
      }  // end loop over layers
    }    // end loop over planes

    //----------------------------------------------------------------------
    //         REDISTRIBUTE ACCORDING TO THE GENERATED CONNECTIVITY
    //----------------------------------------------------------------------

    // time measurement --- start TimeMonitor tm5
    tm5_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcreddis_));


    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << "Redistributing: \n";
      fflush(stdout);
    }

    RedistributeAndCreateDofCoupling();

    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << "... done\n";
      fflush(stdout);
    }

    // time measurement --- this causes the TimeMonitor tm5 to stop here
    tm5_ref_ = Teuchos::null;
  }  // if (numpbcpairs_ > 2)
  return;
}  // PutAllSlavesToMastersProc()

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
    std::map<int, std::vector<int>>& midtosid, const std::vector<int> masternodeids,
    const std::vector<int> slavenodeids, const std::vector<int> dofsforpbcplane,
    const double rotangle, const double abstol)
{
  // these are just parameter definitions for the octree search algorithm
  double tol = abstol;
  int maxnodeperleaf = 250;

  //----------------------------------------------------------------------
  //                   BUILD PROCESSOR LOCAL OCTREE
  //----------------------------------------------------------------------
  // time measurement --- start TimeMonitor tm2
  tm2_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcmidoct_));

  // build processor local octree
  DRT::UTILS::NodeMatchingOctree nodematchingoctree = DRT::UTILS::NodeMatchingOctree();

  nodematchingoctree.Init(*discret_, masternodeids, maxnodeperleaf, tol);
  nodematchingoctree.Setup();
  // time measurement --- this causes the TimeMonitor tm2 to stop here
  tm2_ref_ = Teuchos::null;

  //----------------------------------------------------------------------
  //  SEARCH CLOSEST NODES IN OCTREES ON ALL PROCESSORS
  //----------------------------------------------------------------------

  // time measurement --- start TimeMonitor tm3
  tm3_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcmidmatch_));
  // create connectivity for this condition in this direction
  {
    // create map from gid masternode -> gid corresponding slavenode
    nodematchingoctree.CreateGlobalEntityMatching(
        slavenodeids, dofsforpbcplane, rotangle, midtosid);
  }

  // time measurement --- this causes the TimeMonitor tm3 to stop here
  tm3_ref_ = Teuchos::null;

  return;
}  // PeriodicBoundaryConditions::CreateNodeCouplingForSinglePBC


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
    std::map<int, std::vector<int>>& midtosid, const int pbcid)
{
  // the "inverse" mapping of allcoupled(row/col)nodes
  //       slave node -> his master node (list of size 1)
  Teuchos::RCP<std::map<int, std::vector<int>>> inversenodecoupling;
  inversenodecoupling = Teuchos::rcp(new std::map<int, std::vector<int>>);

  // Teuchos::rcp to the constructed rowmap
  Teuchos::RCP<Epetra_Map> newrownodemap;

  //----------------------------------------------------------------------
  //  ADD THE CONNECTIVITY FROM THIS CONDITION TO THE CONNECTIVITY OF
  //                           ALL CONDITIONS
  //----------------------------------------------------------------------
  {
    std::map<int, std::vector<int>>::iterator iter;
    int masterid;
    int slaveid;
    bool alreadyinlist;

    for (iter = midtosid.begin(); iter != midtosid.end(); ++iter)
    {
      // get id of masternode and slavenode
      masterid = iter->first;

      std::vector<int>::iterator i;
      for (i = (iter->second).begin(); i != (iter->second).end(); ++i)
      {
        slaveid = *i;
        if (slaveid == masterid)
          dserror("Node %d is master AND slave node of periodic boundary condition", masterid);

        // is masterid already in allcoupledrownodes?
        {
          alreadyinlist = false;

          std::map<int, std::vector<int>>::iterator found;

          found = allcoupledrownodes_->find(masterid);
          if (found != allcoupledrownodes_->end())
          {
            // masterid is already in the list --- i.e., the master is the
            // master of a previous condition. Simply append the slave id here
            alreadyinlist = true;
            found->second.push_back(slaveid);
          }
          // masterid is not in the list yet. -> new entry
          if (alreadyinlist == false)
          {
            (*allcoupledrownodes_)[masterid].push_back(slaveid);
          }  // end if not in map
        }
      }
    }  // end insert entries of midtosid into the allcoupledrownodes map

    //---------------------------------------------------------------
    //     COMPLETE MATCHING FOR NODES WITH MORE THAN ONE PBC
    //---------------------------------------------------------------
    // complete matching --- we are able to do this because of the
    // communication step (first step is no problem at all ...)

    {
      if (pbcid > 0)
      {
        // 1) each proc generates a list of his multiple coupled masters
        // 2) the list is communicated in a round robin pattern to all the
        //    other procs.
        // 3) the proc checks the package from each proc and inserts missing
        //    links into the multiple coupling


        //--------------------------------------------------------------------
        // -> 1) create a list of multiple master
        // Communicate multiple couplings for completion...
        std::map<int, std::vector<int>> multiplecouplings;
        for (iter = midtosid.begin(); iter != midtosid.end(); ++iter)
        {
          // get id of masternode and the node itself
          masterid = iter->first;
          DRT::Node* actnode = discret_->gNode(masterid);

          // get all periodic boundary conditions on this node
          std::vector<DRT::Condition*> thiscond;
          actnode->GetCondition("SurfacePeriodic", thiscond);

          if (thiscond.empty())
          {
            actnode->GetCondition("LinePeriodic", thiscond);
          }

          // loop them and check, whether this is a pbc pure master node
          // for all previous conditions
          unsigned ntimesmaster = 0;
          for (unsigned numcond = 0; numcond < thiscond.size(); ++numcond)
          {
            const std::string* mymasterslavetoggle =
                thiscond[numcond]->Get<std::string>("Is slave periodic boundary condition");

            if (*mymasterslavetoggle == "Master")
            {
              ++ntimesmaster;
            }  // end is slave?
          }    // end loop this conditions

          if (ntimesmaster == thiscond.size())
          {
            // yes, we have such a pure master node
            std::vector<int> thiscoupling;
            for (std::vector<int>::iterator rr = (*allcoupledrownodes_)[masterid].begin();
                 rr != (*allcoupledrownodes_)[masterid].end(); ++rr)
            {
              thiscoupling.push_back(*rr);
            }

            // add it to the list of multiple coupled masters on this proc
            multiplecouplings[masterid] = thiscoupling;
          }
        }

        //--------------------------------------------------------------------
        // -> 2) round robin loop

        const int numproc = discret_->Comm().NumProc();
        const int myrank = discret_->Comm().MyPID();            // me
        const int torank = (myrank + 1) % numproc;              // to
        const int fromrank = (myrank + numproc - 1) % numproc;  // from

        DRT::Exporter exporter(discret_->Comm());


        for (int irobin = 0; irobin < numproc; ++irobin)
        {
          std::vector<char> sdata;
          std::vector<char> rdata;

          // ---- pack data for sending -----
          {
            DRT::PackBuffer data;

            std::vector<int> mids;
            for (std::map<int, std::vector<int>>::const_iterator iter = multiplecouplings.begin();
                 iter != multiplecouplings.end(); ++iter)
              mids.push_back(iter->first);

            DRT::ParObject::AddtoPack(data, mids);
            for (std::map<int, std::vector<int>>::const_iterator iter = multiplecouplings.begin();
                 iter != multiplecouplings.end(); ++iter)
              DRT::ParObject::AddtoPack(data, iter->second);

            data.StartPacking();

            DRT::ParObject::AddtoPack(data, mids);
            for (std::map<int, std::vector<int>>::const_iterator iter = multiplecouplings.begin();
                 iter != multiplecouplings.end(); ++iter)
              DRT::ParObject::AddtoPack(data, iter->second);

            std::swap(sdata, data());
          }

          // ---- send ----
          MPI_Request request;
          exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1337, request);

          // ---- receive ----
          int length = rdata.size();
          int tag = -1;
          int from = -1;
          exporter.ReceiveAny(from, tag, rdata, length);
          if (tag != 1337 or from != fromrank)
            dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
                myrank, from, myrank);

          // ---- unpack ----
          {
            multiplecouplings.clear();
            size_t pos = 0;
            std::vector<int> mids;
            DRT::ParObject::ExtractfromPack(pos, rdata, mids);

            for (std::vector<int>::const_iterator iter = mids.begin(); iter != mids.end(); ++iter)
            {
              std::vector<int> slvs;
              DRT::ParObject::ExtractfromPack(pos, rdata, slvs);
              multiplecouplings[*iter] = slvs;
            }
          }

          // wait for all communication to finish
          exporter.Wait(request);
          discret_->Comm().Barrier();  // I feel better this way ;-)

          //--------------------------------------------------
          // -> 3) Try to complete the matchings

          for (std::map<int, std::vector<int>>::iterator mciter = multiplecouplings.begin();
               mciter != multiplecouplings.end(); ++mciter)
          {
            size_t len = mciter->second.size();
            for (size_t mm = 0; mm < len; ++mm)  // this cannot be done through an iterator since we
                                                 // append to the vector while looping over it.
            {
              int possiblemaster = (mciter->second)[mm];

              std::map<int, std::vector<int>>::iterator found =
                  allcoupledrownodes_->find(possiblemaster);

              if (found != allcoupledrownodes_->end())
              {
                // close the connectivity using the slave node which was the
                // masternode of the previous condition
                for (std::vector<int>::iterator fsiter = found->second.begin();
                     fsiter != found->second.end(); ++fsiter)
                {
                  bool doit = true;
                  for (std::vector<int>::const_iterator innersiter = mciter->second.begin();
                       innersiter != mciter->second.end(); ++innersiter)
                    if (*fsiter == *innersiter) doit = false;

                  if (doit) mciter->second.push_back(*fsiter);
                }
                allcoupledrownodes_->erase(found);
              }  // end if we have a further connectivity information...
            }
          }
        }

        // add this information to the map of all coupled nodes
        for (std::map<int, std::vector<int>>::iterator mciter = multiplecouplings.begin();
             mciter != multiplecouplings.end(); ++mciter)
        {
          std::map<int, std::vector<int>>::iterator found =
              allcoupledrownodes_->find(mciter->first);

          if (found != allcoupledrownodes_->end())
          {
            for (size_t i = found->second.size(); i < mciter->second.size(); ++i)
              found->second.push_back(mciter->second[i]);
          }
        }
      }
    }  // end complete matching
  }

  return;
}  // PeriodicBoundaryConditions::AddConnectivity


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

void PeriodicBoundaryConditions::RedistributeAndCreateDofCoupling()
{
  // the "inverse" mapping of allcoupled(row/col)nodes
  //       slave node -> his master node (list of size 1)
  Teuchos::RCP<std::map<int, std::vector<int>>> inversenodecoupling;
  inversenodecoupling = Teuchos::rcp(new std::map<int, std::vector<int>>);

  // Teuchos::rcp to the constructed rowmap
  Teuchos::RCP<Epetra_Map> newrownodemap;

  {
    // time measurement --- start TimeMonitor tm6
    tm6_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcmakeghostmap_));

    // make sure we have a filled discretisation at this place
    // dofs are not required yet, they are assigned after redistribution
    // accessing the noderowmap requires a 'completed' discretization
    if (!discret_->Filled())
    {
      discret_->FillComplete(false, false, false);
    }

    // a list of all nodes on this proc
    std::vector<int> nodesonthisproc(discret_->NodeRowMap()->NumMyElements());

    // get all node gids of nodes on this proc
    discret_->NodeRowMap()->MyGlobalElements(&nodesonthisproc[0]);

    std::set<int> nodeset;

    for (std::vector<int>::const_iterator rr = nodesonthisproc.begin(); rr != nodesonthisproc.end();
         ++rr)
    {
      nodeset.insert(*rr);
    }

    // -----------------------------------------------
    // remove all node gids of slave nodes on this proc

    // get all periodic boundary conditions on this node
    std::vector<DRT::Condition*> thiscond;

    std::vector<DRT::Condition*> linecond;
    discret_->GetCondition("LinePeriodic", linecond);

    for (std::vector<DRT::Condition*>::iterator cond = linecond.begin(); cond != linecond.end();
         ++cond)
    {
      thiscond.push_back(*cond);
    }
    std::vector<DRT::Condition*> surfcond;
    discret_->GetCondition("SurfacePeriodic", surfcond);
    for (std::vector<DRT::Condition*>::iterator cond = surfcond.begin(); cond != surfcond.end();
         ++cond)
    {
      thiscond.push_back(*cond);
    }

    int myerase = 0;
    int numerase = 0;

    int mycerase = 0;
    int numcerase = 0;

    for (unsigned numcond = 0; numcond < thiscond.size(); ++numcond)
    {
      const std::string* mymasterslavetoggle =
          thiscond[numcond]->Get<std::string>("Is slave periodic boundary condition");

      if (*mymasterslavetoggle == "Slave")
      {
        const std::vector<int>* slaveidstodel;

        slaveidstodel = thiscond[numcond]->Nodes();

        for (std::vector<int>::const_iterator idtodel = (*slaveidstodel).begin();
             idtodel != (*slaveidstodel).end(); ++idtodel)
        {
          if (discret_->HaveGlobalNode(*idtodel))
          {
            // erase the coupled nodes from the map --- they are redundant
            allcoupledrownodes_->erase(*idtodel);

            DRT::Node* actnode = discret_->gNode(*idtodel);

            // check for row nodesactnodes ??????????????????
            if (actnode->Owner() != discret_->Comm().MyPID())
            {
              continue;
            }


            ++mycerase;

            std::set<int>::iterator curr = nodeset.find(*idtodel);
            if (curr != nodeset.end())
            {
              // erase id from vector
              nodeset.erase(curr);
              ++myerase;
            }
          }
        }
      }
    }

    discret_->Comm().SumAll(&myerase, &numerase, 1);
    discret_->Comm().SumAll(&mycerase, &numcerase, 1);
    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << " Erased " << numerase << " slaves from nodeset.\n";
      std::cout << " Erased " << numcerase << " from the map of all ";
      std::cout << "coupled rownodes.\n";
      std::cout << "        (we want a master->slaves map, so all entries ";
      std::cout << "slave->... are deleted)\n";
    }

    nodesonthisproc.clear();

    for (std::set<int>::iterator rr = nodeset.begin(); rr != nodeset.end(); ++rr)
    {
      nodesonthisproc.push_back(*rr);
    }

    int mynumappend = 0;
    int numappend = 0;

    // append slavenodes to this list of nodes on this proc
    {
      for (std::map<int, std::vector<int>>::iterator curr = allcoupledrownodes_->begin();
           curr != allcoupledrownodes_->end(); ++curr)
      {
        for (std::vector<int>::iterator iter = curr->second.begin(); iter != curr->second.end();
             ++iter)
        {
          int slaveid = *iter;

          nodesonthisproc.push_back(slaveid);
          ++mynumappend;
        }
      }
    }

    discret_->Comm().SumAll(&mynumappend, &numappend, 1);
    if (discret_->Comm().MyPID() == 0 && verbose_)
    {
      std::cout << " Appended " << numappend << " ids which belong to slave ";
      std::cout << "nodes that are coupled to a master\n";
      std::cout << "        node. They will be fetched to the master's ";
      std::cout << "procs, an their total \n";
      std::cout << "        number has to equal the number of slaves ";
      std::cout << "erased from the nodeset.\n";
    }



    {
      int mymax = 0;
      int max = 0;
      int mymin = 1000;
      int min = 0;

      int myallcouplednodes = allcoupledrownodes_->size();
      int allcouplednodes = 0;

      for (std::map<int, std::vector<int>>::iterator curr = allcoupledrownodes_->begin();
           curr != allcoupledrownodes_->end(); ++curr)
      {
        if ((int)curr->second.size() > mymax)
        {
          mymax = curr->second.size();
        }
        if ((int)curr->second.size() < mymin)
        {
          mymin = curr->second.size();
        }
      }

      discret_->Comm().SumAll(&myallcouplednodes, &allcouplednodes, 1);
      discret_->Comm().MaxAll(&mymax, &max, 1);
      discret_->Comm().MinAll(&mymin, &min, 1);

      if (discret_->Comm().MyPID() == 0 && verbose_)
      {
        std::cout << " The layout is generated: " << allcouplednodes
                  << " masters are coupled to at least " << min << " and up to " << max
                  << " slaves,\n";
        std::cout << "        all master/slave couples are sent to the same proc.\n";
      }
    }


    {
      int myn = (int)nodesonthisproc.size();
      int gn = 0;

      discret_->Comm().SumAll(&myn, &gn, 1);

      if (gn != discret_->NumGlobalNodes())
      {
        dserror(
            "Unmatching numbers of nodes before and after call Redistribution. Nodemap constructor "
            "will crash.\n");
      }
    }

    //--------------------------------------------------
    // build noderowmap for new distribution of nodes
    newrownodemap = Teuchos::rcp(new Epetra_Map(discret_->NumGlobalNodes(), nodesonthisproc.size(),
        &nodesonthisproc[0], 0, discret_->Comm()));

    // create nodal graph of problem, according to old RowNodeMap
    Teuchos::RCP<Epetra_CrsGraph> oldnodegraph = discret_->BuildNodeGraph();

    // export the graph to newrownodemap
    Epetra_CrsGraph nodegraph(Copy, *newrownodemap, 108, false);

    {
      Epetra_Export exporter(*discret_->NodeRowMap(), *newrownodemap);
      int err = nodegraph.Export(*oldnodegraph, exporter, Add);
      if (err < 0) dserror("Graph export returned err=%d", err);
    }
    nodegraph.FillComplete();
    nodegraph.OptimizeStorage();

    // build nodecolmap for new distribution of nodes
    const Epetra_BlockMap cntmp = nodegraph.ColMap();

    Teuchos::RCP<Epetra_Map> newcolnodemap;

    newcolnodemap = Teuchos::rcp(
        new Epetra_Map(-1, cntmp.NumMyElements(), cntmp.MyGlobalElements(), 0, discret_->Comm()));

    // time measurement --- this causes the TimeMonitor tm6 to stop here
    tm6_ref_ = Teuchos::null;


    // time measurement --- start TimeMonitor tm7
    tm7_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcghost_));

    //----------------------------------------------------------------------
    //       GHOSTED NODES NEED INFORMATION ON THEIR COUPLED NODES
    //----------------------------------------------------------------------

    // create the inverse map --- slavenode -> masternode
    inversenodecoupling->clear();

    for (std::map<int, std::vector<int>>::iterator curr = allcoupledrownodes_->begin();
         curr != allcoupledrownodes_->end(); ++curr)
    {
      for (unsigned rr = 0; rr < curr->second.size(); ++rr)
      {
        (*inversenodecoupling)[curr->second[rr]].push_back(curr->first);
      }
    }

    *allcoupledcolnodes_ = (*allcoupledrownodes_);
#ifdef PARALLEL
    {
      // create an exporter
      DRT::Exporter exportconnectivity(*newrownodemap, *newcolnodemap, discret_->Comm());

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
      std::vector<int> mycolnodes(newcolnodemap->NumMyElements());
      newcolnodemap->MyGlobalElements(&(mycolnodes[0]));

      // determine all ghosted slave nodes in this vector which do not have
      // a ghosted master on this proc --- we have to fetch it to be able
      // to assign the dofs
      for (std::map<int, std::vector<int>>::iterator curr = inversenodecoupling->begin();
           curr != inversenodecoupling->end(); ++curr)
      {
        if (curr->second.empty())
        {
          dserror("inverse slave-master matching incomplete");
        }
        int mymaster = curr->second[0];
        if (newcolnodemap->LID(mymaster) < 0)
        {
          // was master already added to the list of (ghosted) nodes?
          std::vector<int>::iterator found;
          found = find(mycolnodes.begin(), mycolnodes.end(), mymaster);
          // no, it's not inside
          if (found == mycolnodes.end())
          {
            mycolnodes.push_back(mymaster);
          }
        }
      }

      // We need to do this again in order to get the new master-slave pairs
      // that might have been added in the previous loop over the inversenodecoupling
#ifdef PARALLEL
      {
        // now reconstruct the extended colmap
        newcolnodemap = Teuchos::rcp(
            new Epetra_Map(-1, mycolnodes.size(), &mycolnodes[0], 0, discret_->Comm()));

        *allcoupledcolnodes_ = (*allcoupledrownodes_);

        // create an exporter
        DRT::Exporter exportconnectivity(*newrownodemap, *newcolnodemap, discret_->Comm());

        // export information on all master->slave couplings (with multiple
        // couplings)
        exportconnectivity.Export(*allcoupledcolnodes_);
      }
#endif

      // determine all ghosted master nodes in this vector which do not have
      // all their slaves ghosted on this proc --- we have to fetch them to be able
      // to assign the dofs
      for (std::map<int, std::vector<int>>::iterator curr = allcoupledcolnodes_->begin();
           curr != allcoupledcolnodes_->end(); ++curr)
      {
        if (curr->second.empty())
        {
          dserror("master-slave matching incomplete");
        }
        for (size_t i = 0; i < curr->second.size(); ++i)
        {
          if (newcolnodemap->LID(curr->second[i]) < 0)
          {
            // was slave already added to the list of (ghosted) nodes?
            std::vector<int>::iterator found;
            found = find(mycolnodes.begin(), mycolnodes.end(), curr->second[i]);
            // no, it's not inside
            if (found == mycolnodes.end())
            {
              mycolnodes.push_back(curr->second[i]);
            }
          }
        }
      }

      // now reconstruct the extended colmap
      newcolnodemap =
          Teuchos::rcp(new Epetra_Map(-1, mycolnodes.size(), &mycolnodes[0], 0, discret_->Comm()));

      *allcoupledcolnodes_ = (*allcoupledrownodes_);
#ifdef PARALLEL
      // the new master-ghost nodes need their information about
      // connectivity
      {
        // create an exporter
        DRT::Exporter exportconnectivity(*newrownodemap, *newcolnodemap, discret_->Comm());
        // export information on all slave->master couplings (with multiple
        // couplings)
        exportconnectivity.Export(*allcoupledcolnodes_);
      }
#endif
    }

    // time measurement --- this causes the TimeMonitor tm7 to stop here
    tm7_ref_ = Teuchos::null;


    // time measurement --- start TimeMonitor tm8
    tm8_ref_ = Teuchos::rcp(new Teuchos::TimeMonitor(*timepbcrenumdofs_));

    // check whether we have already passed a PBCDofSet to the discretization
    // If we did not the regular DofSet is replaced with a PBCDofSet. This will
    // lead to a new offset of the DofGIDs and therefore make exporting of vectors
    // impossible.
    // If we already passed a PBCDofSet to the discretization we simply update the
    // DofSet and therefore maintain a correct DofGIDs-offset.
    if (pbcdofset_ == Teuchos::null)
    {
      // create a new dofset specialisation for periodic boundary conditions
      // the 'true' flag makes sure that the pbc dofset replaces the old
      // dofset also in the static_dofsets_.
      pbcdofset_ = Teuchos::rcp(new DRT::PBCDofSet(allcoupledcolnodes_));
      discret_->ReplaceDofSet(0, pbcdofset_, true);
    }
    else
    {
      // the discretization already has a pbc dofset, we merely need to update it
      // (a replace dofset is also not needed since we are working on pointer)
      pbcdofset_->SetCoupledNodes(allcoupledcolnodes_);
      pbcdofset_->Reset();
    }

    //--------------------------------------------------
    // redistribute the nodes
    //
    // this contains a call to FillComplete and assigns the same
    // degree of freedom to the matching nodes

    discret_->Redistribute(*newrownodemap, *newcolnodemap);

    // time measurement --- this causes the TimeMonitor tm8 to stop here
    tm8_ref_ = Teuchos::null;

    // throw away old nodegraph
    oldnodegraph = Teuchos::null;
  }

  return;
}  // PeriodicBoundaryConditions::RedistributeAndCreateDofCoupling

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
  if (discret_->Comm().NumProc() > 1)
  {
    const Epetra_Map* noderowmap = discret_->NodeRowMap();

    // weights for graph partition
    Epetra_Vector weights(*noderowmap, false);
    weights.PutScalar(10.0);

    // apply weight of special elements
    for (int inode = 0; inode < noderowmap->NumMyElements(); ++inode)
    {
      int gid = noderowmap->GID(inode);
      DRT::Node* node = discret_->gNode(gid);
      if (!node) dserror("cant find node");
      double weight = 0.0;
      DRT::Element** surrele = node->Elements();

      // returns 10.0 for standard elements
      for (int k = 0; k < node->NumElement(); ++k)
        weight = std::max(weight, surrele[k]->EvaluationCost());

      weights.ReplaceMyValue(inode, 0, weight);
    }
    // ----------------------------------------
    // loop masternodes to adjust weights of slavenodes
    // they need a small weight since they do not contribute any dofs
    // to the linear system
    {
      std::map<int, std::vector<int>>::iterator masterslavepair;

      for (masterslavepair = allcoupledcolnodes_->begin();
           masterslavepair != allcoupledcolnodes_->end(); ++masterslavepair)
      {
        // get masternode
        DRT::Node* master = discret_->gNode(masterslavepair->first);

        if (master->Owner() != discret_->Comm().MyPID())
        {
          continue;
        }

        // loop slavenodes associated with master
        for (std::vector<int>::iterator iter = masterslavepair->second.begin();
             iter != masterslavepair->second.end(); ++iter)
        {
          double initval = 1.0;
          int gid = *iter;

          weights.ReplaceGlobalValues(1, &initval, &gid);
        }
      }
    }

    // allocate graph
    Teuchos::RCP<Epetra_CrsGraph> nodegraph =
        Teuchos::rcp(new Epetra_CrsGraph(Copy, *noderowmap, 108, false));

    // -------------------------------------------------------------
    // iterate all elements on this proc including ghosted ones
    // compute connectivity

    // standard part without master<->slave coupling
    // Note:
    // if a proc stores the appropiate ghosted elements, the resulting
    // graph will be the correct and complete graph of the distributed
    // discretization even if nodes are not ghosted.

    for (int nele = 0; nele < discret_->NumMyColElements(); ++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lColElement(nele);

      // get its nodes and nodeids
      const int nnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      for (int row = 0; row < nnode; ++row)
      {
        const int rownode = nodeids[row];

        // insert into line of graph only when this proc owns the node
        if (!noderowmap->MyGID(rownode)) continue;

        // insert all neighbours from element in the graph
        for (int col = 0; col < nnode; ++col)
        {
          int colnode = nodeids[col];
          int err = nodegraph->InsertGlobalIndices(rownode, 1, &colnode);
          if (err < 0) dserror("nodegraph->InsertGlobalIndices returned err=%d", err);
        }
      }
    }

    // -------------------------------------------------------------
    // additional coupling between master and slave
    // we do not only connect master and slave nodes but if a master/slave
    // is connected to a master/slave, we connect the corresponding slaves/master
    // as well

    for (int nele = 0; nele < discret_->NumMyColElements(); ++nele)
    {
      // get the element
      DRT::Element* ele = discret_->lColElement(nele);

      // get its nodes and nodeids
      const int nnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      for (int row = 0; row < nnode; ++row)
      {
        const int rownode = nodeids[row];

        // insert into line of graph only when this proc owns the node
        if (!noderowmap->MyGID(rownode)) continue;

        std::map<int, std::vector<int>>::iterator masterslavepair =
            allcoupledcolnodes_->find(rownode);
        if (masterslavepair != allcoupledcolnodes_->end())
        {
          // get all masternodes of this element
          for (int col = 0; col < nnode; ++col)
          {
            int colnode = nodeids[col];

            std::map<int, std::vector<int>>::iterator othermasterslavepair =
                allcoupledcolnodes_->find(colnode);
            if (othermasterslavepair != allcoupledcolnodes_->end())
            {
              // add connection to all slaves

              for (std::vector<int>::iterator iter = othermasterslavepair->second.begin();
                   iter != othermasterslavepair->second.end(); ++iter)
              {
                int othermastersslaveindex = *iter;
                int masterindex = rownode;
                int err = nodegraph->InsertGlobalIndices(rownode, 1, &othermastersslaveindex);
                if (err < 0) dserror("nodegraph->InsertGlobalIndices returned err=%d", err);

                if (noderowmap->MyGID(*iter))
                {
                  err = nodegraph->InsertGlobalIndices(*iter, 1, &masterindex);
                  if (err < 0) dserror("nodegraph->InsertGlobalIndices returned err=%d", err);
                }
              }
            }
          }
        }
      }
    }

    // finalize construction of initial graph
    int err = nodegraph->FillComplete();
    if (err) dserror("graph->FillComplete returned %d", err);

    const int myrank = nodegraph->Comm().MyPID();
    const int numproc = nodegraph->Comm().NumProc();

    if (numproc > 1)
    {
      // proc that will do the serial partitioning
      // the graph is collapsed to this proc
      // Normally this would be proc 0 but 0 always has so much to do.... ;-)
      int workrank = 1;

      // get rowmap of the graph
      const Epetra_BlockMap& tmp = nodegraph->RowMap();
      Epetra_Map rowmap(tmp.NumGlobalElements(), tmp.NumMyElements(), tmp.MyGlobalElements(), 0,
          nodegraph->Comm());

      // -------------------------------------------------------------
      // build a target map that stores everything on proc workrank
      // We have arbirtary gids here and we do not tell metis about
      // them. So we have to keep rowrecv until the redistributed map is
      // build.

      // rowrecv is a fully redundant vector (size of number of nodes)
      std::vector<int> rowrecv(rowmap.NumGlobalElements());

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
          (myrank == workrank) ? (int)rowrecv.size() : 0, &rowrecv[0], 0, rowmap.Comm());

      // export the graph to tmap
      Epetra_CrsGraph tgraph(Copy, tmap, 108, false);
      Epetra_Export exporter(rowmap, tmap);
      {
        int err = tgraph.Export(*nodegraph, exporter, Add);
        if (err < 0) dserror("Graph export returned err=%d", err);
      }
      tgraph.FillComplete();
      tgraph.OptimizeStorage();

      // export the weights to tmap
      Epetra_Vector tweights(tmap, false);
      err = tweights.Export(weights, exporter, Insert);
      if (err < 0) dserror("Vector export returned err=%d", err);

      // metis requests indexes. So we need a reverse lookup from gids
      // to indexes.
      std::map<int, int> idxmap;
      // xadj points from index i to the index of the
      // first adjacent node
      std::vector<int> xadj(rowmap.NumGlobalElements() + 1);
      // a list of adjacent nodes, adressed using xadj
      std::vector<int> adjncy(tgraph.NumGlobalNonzeros());  // the size is an upper bound

      // This is a vector of size n that upon successful completion stores the partition vector of
      // the graph
      std::vector<int> part(tmap.NumMyElements());

      // construct reverse lookup for all procs
      for (unsigned i = 0; i < rowrecv.size(); ++i)
      {
        idxmap[rowrecv[i]] = i;
      }

      if (myrank == workrank)
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

        int count = 0;
        xadj[0] = 0;
        for (int row = 0; row < tgraph.NumMyRows(); ++row)
        {
          int grid = tgraph.RowMap().GID(row);
          int numindices;
          int* lindices;
          int err = tgraph.ExtractMyRowView(row, numindices, lindices);
          if (err) dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d", err);

          for (int col = 0; col < numindices; ++col)
          {
            int gcid = tgraph.ColMap().GID(lindices[col]);
            if (gcid == grid) continue;
            adjncy[count] = idxmap[gcid];
            ++count;
          }
          xadj[row + 1] = count;
        }
      }

      // broadcast xadj
      tmap.Comm().Broadcast(&xadj[0], xadj.size(), workrank);

      // broadcast adjacence (required for edge weights)
      int adjncysize = (int)adjncy.size();
      tmap.Comm().Broadcast(&adjncysize, 1, workrank);
      adjncy.resize(adjncysize);
      tmap.Comm().Broadcast(&adjncy[0], adjncysize, workrank);

      // -------------------------------------------------------------
      // set a fully redundant vector of weights for edges
      std::vector<int> ladjwgt(adjncy.size(), 0);
      std::vector<int> adjwgt(adjncy.size(), 0);

      for (std::vector<int>::iterator iter = ladjwgt.begin(); iter != ladjwgt.end(); ++iter)
      {
        *iter = 0;
      }

      // loop all master nodes on this proc
      std::map<int, std::vector<int>>::iterator masterslavepair;

      for (masterslavepair = allcoupledcolnodes_->begin();
           masterslavepair != allcoupledcolnodes_->end(); ++masterslavepair)
      {
        // get masternode
        DRT::Node* master = discret_->gNode(masterslavepair->first);

        if (master->Owner() != myrank)
        {
          continue;
        }

        std::map<int, int>::iterator paul = idxmap.find(master->Id());
        if (paul == idxmap.end())
        {
          dserror("master not in reverse lookup");
        }

        // inverse lookup
        int masterindex = idxmap[master->Id()];

        // loop slavenodes
        for (std::vector<int>::iterator iter = masterslavepair->second.begin();
             iter != masterslavepair->second.end(); ++iter)
        {
          DRT::Node* slave = discret_->gNode(*iter);

          if (slave->Owner() != myrank)
          {
            dserror("own master but not slave\n");
          }

          int slaveindex = idxmap[slave->Id()];

          std::map<int, int>::iterator foo = idxmap.find(slave->Id());
          if (foo == idxmap.end())
          {
            dserror("slave not in reverse lookup");
          }

          // -------------------------------------------------------------
          // connections between master and slavenodes are very strong
          // we do not want to partition between master and slave nodes
          for (int j = xadj[masterindex]; j < xadj[masterindex + 1]; ++j)
          {
            if (adjncy[j] == slaveindex)
            {
              ladjwgt[j] = 100;
            }
          }

          for (int j = xadj[slaveindex]; j < xadj[slaveindex + 1]; ++j)
          {
            if (adjncy[j] == masterindex)
            {
              ladjwgt[j] = 100;
            }
          }
        }
      }

      // do communication to aquire edge weight information from all procs
      tmap.Comm().SumAll(&ladjwgt[0], &adjwgt[0], adjwgt.size());

      // the standard edge weight is one
      for (std::vector<int>::iterator iter = adjwgt.begin(); iter != adjwgt.end(); ++iter)
      {
        if (*iter == 0) *iter = 1;
      }

      // the reverse lookup is not required anymore
      idxmap.clear();

      // -------------------------------------------------------------
      // do partitioning using metis on workrank
      if (myrank == workrank)
      {
        // the vertex weights
        std::vector<int> vwgt(tweights.MyLength());
        for (int i = 0; i < tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

        // 0 No weights (vwgts and adjwgt are NULL)
        // 1 Weights on the edges only (vwgts = NULL)
        // 2 Weights on the vertices only (adjwgt = NULL)
        // 3 Weights both on vertices and edges.
        int wgtflag = 3;
        // 0 C-style numbering is assumed that starts from 0
        // 1 Fortran-style numbering is assumed that starts from 1
        int numflag = 0;
        // The number of parts to partition the graph.
        int npart = numproc;
        // This is an array of 5 integers that is used to pass parameters for the various phases of
        // the algorithm. If options[0]=0 then default values are used. If options[0]=1, then the
        // remaining four elements of options are interpreted as follows: options[1]    Determines
        // matching type. Possible values are:
        //               1 Random Matching (RM)
        //               2 Heavy-Edge Matching (HEM)
        //               3 Sorted Heavy-Edge Matching (SHEM) (Default)
        //               Experiments has shown that both HEM and SHEM perform quite well.
        // options[2]    Determines the algorithm used during initial partitioning. Possible values
        // are:
        //               1 Region Growing (Default)
        // options[3]    Determines the algorithm used for re%Gï¬%@nement. Possible values
        // are:
        //               1 Early-Exit Boundary FM re%Gï¬%@nement (Default)
        // options[4]    Used for debugging purposes. Always set it to 0 (Default).
        int options[5] = {0, 3, 1, 1, 0};
        // Upon successful completion, this variable stores the number of edges that are cut by the
        // partition.
        int edgecut = 0;
        // The number of vertices in the graph.
        int nummyele = tmap.NumMyElements();

        if (verbose_)
        {
          std::cout << "proc " << myrank << " repartition graph using metis\n";
        }
        if (numproc < 8)  // better for smaller no. of partitions
        {
#ifdef PARALLEL
          METIS_PartGraphRecursive(&nummyele, &xadj[0], &adjncy[0], &vwgt[0], &adjwgt[0], &wgtflag,
              &numflag, &npart, options, &edgecut, &part[0]);

          if (verbose_)
          {
            std::cout << "METIS_PartGraphRecursive produced edgecut of " << edgecut << "\n";
            fflush(stdout);
          }
#endif
        }
        else
        {
#ifdef PARALLEL
          METIS_PartGraphKway(&nummyele, &xadj[0], &adjncy[0], &vwgt[0], &adjwgt[0], &wgtflag,
              &numflag, &npart, options, &edgecut, &part[0]);
#endif
        }


      }  // if (myrank==workrank)

      // broadcast partitioning result
      int size = tmap.NumMyElements();
      tmap.Comm().Broadcast(&size, 1, workrank);
      part.resize(size);
      tmap.Comm().Broadcast(&part[0], size, workrank);

      // loop part and count no. of nodes belonging to me
      // (we reuse part to save on memory)
      int count = 0;
      for (int i = 0; i < size; ++i)
        if (part[i] == myrank)
        {
          part[count] = rowrecv[i];
          ++count;
        }

      // rowrecv is done
      rowrecv.clear();

      // create map with new layout
      Epetra_Map newmap(size, count, &part[0], 0, nodegraph->Comm());

      // create the new graph and export to it
      Teuchos::RCP<Epetra_CrsGraph> newnodegraph;

      newnodegraph = Teuchos::rcp(new Epetra_CrsGraph(Copy, newmap, 108, false));
      Epetra_Export exporter2(nodegraph->RowMap(), newmap);
      err = newnodegraph->Export(*nodegraph, exporter2, Add);
      if (err < 0) dserror("Graph export returned err=%d", err);
      newnodegraph->FillComplete();
      newnodegraph->OptimizeStorage();

      // the rowmap will become the new distribution of nodes
      const Epetra_BlockMap rntmp = newnodegraph->RowMap();
      Epetra_Map newnoderowmap(
          -1, rntmp.NumMyElements(), rntmp.MyGlobalElements(), 0, discret_->Comm());

      // the column map will become the new ghosted distribution of nodes
      const Epetra_BlockMap Mcntmp = newnodegraph->ColMap();
      Epetra_Map newnodecolmap(
          -1, Mcntmp.NumMyElements(), Mcntmp.MyGlobalElements(), 0, discret_->Comm());
      // do the redistribution without assigning dofs
      discret_->Redistribute(newnoderowmap, newnodecolmap, false, true, true);


      if (discret_->Comm().MyPID() == 0 && verbose_)
      {
        std::cout << "---------------------------------------------\n";
        std::cout << "Repair Master->Slave connection, generate final dofset";
        std::cout << std::endl << std::endl;
      }

      // assign the new dofs, make absolutely sure that we always
      // have all slaves to a master
      // the finite edge weights are not a 100% warranty for that...
      PutAllSlavesToMastersProc();
    }
  }

  return;
}  // BalanceLoadUsingMetis
