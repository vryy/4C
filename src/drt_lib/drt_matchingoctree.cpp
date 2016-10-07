/*!----------------------------------------------------------------------
\file drt_matchingoctree.cpp

\brief Search closest node in given set of nodes using an octree search

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/


#include "drt_matchingoctree.H"
#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::MatchingOctree::MatchingOctree() :
  discret_(NULL),
  tol_(-1.0),
  masterentityids_(NULL),
  maxtreenodesperleaf_(-1),
  issetup_(false),
  isinit_(false)
{
  return;
}// MatchingOctree::MatchingOctree

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::UTILS::MatchingOctree::Init(
    const DRT::Discretization&       actdis,
    const std::vector <int> &        masternodeids,
    const int                        maxnodeperleaf,
    const double                     tol )
{
  SetIsSetup(false);

  discret_ = &actdis;
  masterentityids_ = &masternodeids;
  maxtreenodesperleaf_ = maxnodeperleaf;
  tol_ = tol;

  SetIsInit(true);
  return 0;
}// MatchingOctree::Init

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::UTILS::MatchingOctree::Setup()
{
  CheckIsInit();

  const int nummygids = masterentityids_->size();

  // extract all masternodes on this proc from the list masternodeids
  std::vector <int> masternodesonthisproc;

  // construct octree if proc has masternodes
  if(nummygids>0)
  {
    // create initial bounding box for all nodes
    //
    //                 +-            -+
    //                 |  xmin  xmax  |
    //                 |  ymin  ymax  |
    //                 |  zmin  zmax  |
    //                 +-            -+
    //
    Epetra_SerialDenseMatrix initialboundingbox(3,2);
    double pointcoord[3];
    CalcPointCoordinate(discret_,masterentityids_->at(0),&pointcoord[0]);
    for (int dim=0;dim<3;dim++)
    {
      initialboundingbox(dim,0)=pointcoord[dim]-tol_;
      initialboundingbox(dim,1)=pointcoord[dim]+tol_;

      // store coordinates of one point in master plane (later on, one
      // coordinate of the masternode will be substituted by the coordinate
      // of the master plane)
      masterplanecoords_.push_back(pointcoord[dim]);
    }

    for(unsigned locn=0;locn<(unsigned)nummygids;locn++)
    {
      // check if entity is on this proc
      if(not CheckHaveEntity(discret_, masterentityids_->at(locn)))
        dserror("MatchingOctree can only be constructed with entities,\n"
            "which are either owned, or ghosted by calling proc.");

      masternodesonthisproc.push_back(masterentityids_->at(locn));

      CalcPointCoordinate(discret_,masternodesonthisproc[locn],&pointcoord[0]);

      for (int dim=0;dim<3;dim++)
      {
        initialboundingbox(dim,0)=std::min(initialboundingbox(dim,0),
            pointcoord[dim]-tol_);
        initialboundingbox(dim,1)=std::max(initialboundingbox(dim,1),
            pointcoord[dim]+tol_);
      }
    }

    // create octree root --- initial layer is 0
    // all other layers are generated down here by recursive calls
    int initlayer = 0;

    octreeroot_ = CreateOctreeElement(
        masternodesonthisproc,
        initialboundingbox,
        initlayer);
  }

  SetIsSetup(true);
  return 0;
}// MatchingOctree::Setup

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::MatchingOctree::~MatchingOctree()
{
  return;
}// MatchingOctree::~MatchingOctree

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::MatchingOctree::SearchClosestEntityOnThisProc(
  const std::vector<double>& x,
  int           & idofclosestpoint,
  double        & distofclosestpoint,
  bool          searchsecond
  )
{
  // flag
  bool nodeisinbox=false;

  nodeisinbox = octreeroot_->IsPointInBoundingBox(x);

  if (nodeisinbox==true)
  {
    // the node is inside the bounding box. So maybe the closest one is
    // here on this proc -> search for it

    if(octreeroot_==Teuchos::null)
    {
      dserror("No root for octree on proc");
    }

    Teuchos::RCP<OctreeElement> octreeele = octreeroot_;

    while(octreeele->IsLeaf()==false)
    {
      octreeele = octreeele->ReturnChildContainingPoint(x);

      if(octreeele==Teuchos::null)
      {
        dserror("Child is nullpointer");
      }
    }

    // now get closest point in leaf
    octreeele->SearchClosestNodeInLeaf(
        x,
        idofclosestpoint,
        distofclosestpoint,
        tol_,
        searchsecond);
  }

  return nodeisinbox;
}// MatchingOctree::SearchClosestNodeOnThisProc

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MatchingOctree::CreateGlobalEntityMatching(
  const std::vector<int>&          slavenodeids,
  const std::vector<int>&          dofsforpbcplane,
  const double                     rotangle,
  std::map<int,std::vector<int> >& midtosid
  )
{
  CheckIsInit();
  CheckIsSetup();

  int myrank  =discret_->Comm().MyPID();
  int numprocs=discret_->Comm().NumProc();

  // map from global masternodeids to distances to their global slave
  // counterpart
  std::map<int,double > diststom;

  // 1) each proc generates a list of his slavenodes
  // 2) the list is communicated in a round robin pattern to all the
  //    other procs.
  //
  // 3) the proc checks the package from each proc and calcs the min
  //    distance on each --- the result is kept if distance is smaller
  //    than on the preceding processors

  //--------------------------------------------------------------------
  // -> 1) create a list of slave nodes on this proc. Pack it.
  std::vector<char> sblockofnodes;
  std::vector<char> rblockofnodes;

  sblockofnodes.clear();
  rblockofnodes.clear();

  DRT::PackBuffer data;

  for(int globn=0;globn<(int)slavenodeids.size();globn++)
  {
    // is this slavenode on this proc?
    if(CheckHaveEntity(discret_,slavenodeids[globn]))
    {
      // Add node to list of nodes which will be sent to the next proc
      PackEntity(data,discret_,slavenodeids[globn]);
    } // end if slavenode on proc
  } // end loop globn

  data.StartPacking();

  for(int globn=0;globn<(int)slavenodeids.size();globn++)
  {
    // is this slavenode on this proc?
    if(CheckHaveEntity(discret_,slavenodeids[globn]))
    {
      // Add node to list of nodes which will be sent to the next proc
      PackEntity(data,discret_,slavenodeids[globn]);

    } // end if slavenode on proc
  } // end loop globn

  swap( sblockofnodes, data() );

  //--------------------------------------------------------------------
  // -> 2) round robin loop
#ifdef PARALLEL
  // create an exporter for point to point comunication
  DRT::Exporter exporter(discret_->Comm());
#endif

  for (int np=0;np<numprocs;np++)
  {
    //--------------------------------------------------
    // Send block to next proc. Receive a block from the last proc
    if(np >0) // in the first step, we keep all nodes on this proc
    {
#ifdef PARALLEL
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblockofnodes.size();

      exporter.ISend(frompid,topid,
                     &(sblockofnodes[0]),sblockofnodes.size(),
                     tag,request);

      // make sure that you do not think you received something if
      // you didn't
      if(rblockofnodes.empty()==false)
      {
        dserror("rblockofnodes not empty");
      }

      rblockofnodes.clear();

      // receive from predecessor
      frompid=(myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblockofnodes,length);

      if(tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);

      {
        // for safety
        exporter.Comm().Barrier();
      }
#endif
    }
    else
    {
      // dummy communication
      rblockofnodes = sblockofnodes;
    }

    //--------------------------------------------------
    // Unpack block.
    std::vector<char>::size_type index = 0;
    while (index < rblockofnodes.size())
    {
      // extract node data from blockofnodes
      std::vector<char> data;
      DRT::ParObject::ExtractfromPack(index,rblockofnodes,data);

      // allocate an "empty node". Fill it with info from
      // extracted node data
      Teuchos::RCP<DRT::ParObject> o = Teuchos::rcp(DRT::UTILS::Factory(data));

      // check type of ParObject, and return gid
      const int id = CheckValidEntityType(o);

      //----------------------------------------------------------------
      // there is nothing to do if there are no master nodes on this
      // proc
      if(masterplanecoords_.empty()!=true)
      {
        double pointcoord[3];
        CalcPointCoordinate(o.getRawPtr(),&pointcoord[0]);

        // get its coordinates
        std::vector <double> x(3);

        if (abs(rotangle) < EPS13)
        {
          for (int dim=0;dim<3;dim++)
          {
            x[dim] = pointcoord[dim];
          }
        }
        else
        {
          // if there is a rotationally symmetric periodic boundary condition:
          // rotate slave plane for making it parallel to the master plane
          x[0] = pointcoord[0]*cos(rotangle) + pointcoord[1]*sin(rotangle);
          x[1] = pointcoord[0]*(-sin(rotangle)) + pointcoord[1]*cos(rotangle);
          x[2] = pointcoord[2];
        }

        // Substitute the coordinate normal to the master plane by the
        // coordinate of the masterplane
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

        // get direction for parallel translation
        if(!dofsforpbcplane.empty())
        {
          int dir=-1;

          for (int dim=0;dim<3;dim++)
          {
            if(dofsforpbcplane[0]==dim || dofsforpbcplane[1]==dim)
            {
              // direction dim is in plane
                continue;
            }
            else
            {
              dir=dim;
            }
          }

          if (dir<0)
          {
            dserror("Unable to get direction orthogonal to plane");
          }

          // substitute x value
          x[dir] = masterplanecoords_[dir];
        }
        //--------------------------------------------------------
        // 3) now search for closest master point on this proc
        int        idofclosestpoint;
        double     distofclosestpoint;

        bool       nodeisinbox;

        nodeisinbox=this->SearchClosestEntityOnThisProc(
          x,
          idofclosestpoint,
          distofclosestpoint);

        // If x is not in the bounding box on this proc, its probably not
        // matching a point in the box. We do nothing.
        if(nodeisinbox==true)
        {

          std::map<int,std::vector<int> >::iterator found
      = midtosid.find(idofclosestpoint);

          if( found != midtosid.end() )
          {
            // if this is true, we already have an entry in the list and
            // have to check whether this is a better value
            if(diststom[idofclosestpoint] > distofclosestpoint)
            {
        (midtosid[idofclosestpoint]).clear();
        (midtosid[idofclosestpoint]).push_back(id);
              diststom[idofclosestpoint] = distofclosestpoint;
            }
      else if(diststom[idofclosestpoint]<distofclosestpoint+1e-9
        &&
        diststom[idofclosestpoint]>distofclosestpoint-1e-9
        )
      {
        (midtosid[idofclosestpoint]).push_back(id);
      }
          }
          else
          {
            // this is the first estimate for a closest point
      (midtosid[idofclosestpoint]).clear();
      (midtosid[idofclosestpoint]).push_back(id);
            diststom[idofclosestpoint] = distofclosestpoint;

          }
        } // end nodeisinbox==true

      } // end if (masterplanecoords_.empty()!=true)
    }



    //----------------------------------------------------------------
    // prepare to send nodes to next proc (keep list).

    // the received nodes will be sent to the next proc
    sblockofnodes = rblockofnodes;

    // we need a new receive buffer
    rblockofnodes.clear();

#ifdef PARALLEL
    {
      // for safety
      exporter.Comm().Barrier();
    }
#endif
  } // end loop np

  return;
} // MatchingOctree::CreateGlobalNodeMatching

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MatchingOctree::FindMatch(
                                   const DRT::Discretization& slavedis,
                                   const std::vector<int>& slavenodeids,
                                   std::map<int,std::pair<int,double> >& coupling)
{
  CheckIsInit();
  CheckIsSetup();

  int numprocs = discret_->Comm().NumProc();

  if (slavedis.Comm().NumProc()!=numprocs)
    dserror("compared discretizations must live on same procs");

  //------
  // For XFEM based FSI problems, explicit node-matching procedure is not necessary
  // FSI (boundary) discretization is created from soliddis by considering the nodes
  // that belongs to FSICoupling conditions. All the nodes in the FSI interface
  // is contained in soliddis and they have same node numbers
  //------
  if ( DRT::Problem::Instance()->ProblemType() == prb_fsi_xfem or
       DRT::Problem::Instance()->ProblemType() == prb_fsi_crack )
  {
    for( unsigned node = 0 ; node < slavenodeids.size(); node++ )
    {
      const int slaid = slavenodeids[node];
      coupling[slaid] = std::make_pair( slaid, 0.0 );
    }

    return;
  }

  // 1) each proc generates a list of his slavenodes
  //
  // 2) the list is communicated in a round robin pattern to all the
  //    other procs.
  //
  // 3) the proc checks the package from each proc and calcs the min
  //    distance on each --- the result is kept if distance is smaller
  //    than on the preceding processors

  //--------------------------------------------------------------------
  // -> 1) create a list of slave nodes on this proc. Pack it.
  std::vector<char> sblockofnodes;
  std::vector<char> rblockofnodes;

  DRT::PackBuffer data;

  for (unsigned globn=0; globn<slavenodeids.size(); ++globn)
  {
    // is this slavenode on this proc?
    if (CheckHaveEntity(&slavedis, slavenodeids[globn]))
    {
        // Add node to list of nodes which will be sent to the next proc
        PackEntity(data,&slavedis,slavenodeids[globn]);
    }
  }

  data.StartPacking();

  for (unsigned globn=0; globn<slavenodeids.size(); ++globn)
  {
    // is this slavenode on this proc?
    if (CheckHaveEntity(&slavedis, slavenodeids[globn]))
    {
        // Add node to list of nodes which will be sent to the next proc
        PackEntity(data,&slavedis,slavenodeids[globn]);
    }
  }

  swap( sblockofnodes, data() );

  //--------------------------------------------------------------------
  // -> 2) round robin loop

  // create an exporter for point to point comunication
  // We do all communication with the communicator of the original
  // discretization.
  DRT::Exporter exporter(discret_->Comm());

  for (int np=0; np<numprocs; np++)
  {
    //--------------------------------------------------
    // Send block to next proc. Receive a block from the last proc
    if (np > 0) // in the first step, we keep all nodes on this proc
    {
#ifdef PARALLEL
      int myrank   = discret_->Comm().MyPID();
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblockofnodes.size();

      exporter.ISend(frompid,topid,
                     &(sblockofnodes[0]),sblockofnodes.size(),
                     tag,request);

      // make sure that you do not think you received something if
      // you didn't
      if (not rblockofnodes.empty())
      {
        dserror("rblockofnodes not empty");
      }

      rblockofnodes.clear();

      // receive from predecessor
      frompid = (myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblockofnodes,length);

      if (tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);
#else
      dserror("How did you get here? Go away!");
#endif
    }
    else
    {
      // no need to communicate
      swap(rblockofnodes,sblockofnodes);
    }

    //--------------------------------------------------
    // Unpack block.
    std::vector<char>::size_type index = 0;
    while (index < rblockofnodes.size())
    {
      // extract node data from blockofnodes
      std::vector<char> data;
      DRT::ParObject::ExtractfromPack(index,rblockofnodes,data);

      // allocate an "empty node". Fill it with info from
      // extracted node data
      Teuchos::RCP<DRT::ParObject> o = Teuchos::rcp(DRT::UTILS::Factory(data));

      // cast ParObject to specific type and return id
      const int id = CheckValidEntityType(o);

      //----------------------------------------------------------------
      // there is nothing to do if there are no master nodes on this
      // proc
      if (not masterplanecoords_.empty())
      {
        double pointcoord[3];
        CalcPointCoordinate(o.getRawPtr(),&pointcoord[0]);

        // get its coordinates
        std::vector<double> x(pointcoord, pointcoord+3);

        //--------------------------------------------------------
        // 3) now search for closest master point on this proc
        int    gid;
        double dist;

        // If x is not in the bounding box on this proc, its probably not
        // matching a point in the box. We do nothing.
        if (SearchClosestEntityOnThisProc(x, gid, dist))
        {
          std::map<int,std::pair<int,double> >::iterator found = coupling.find(gid);

          // search for second point with same distance, if found gid is already in coupling
          if (found != coupling.end())
          {
            if (SearchClosestEntityOnThisProc(x, gid, dist, true))
            {
              found = coupling.find(gid);
            }
          }

          // we are interested in the closest match
          if (found==coupling.end() or coupling[gid].second > dist)
          {
            coupling[gid] = std::make_pair(id, dist);
          }
        }
      }
    }

    //----------------------------------------------------------------
    // prepare to send nodes to next proc (keep list).

    // the received nodes will be sent to the next proc
    swap(sblockofnodes,rblockofnodes);

    // we need a new receive buffer
    rblockofnodes.clear();
  }
  return;
} // MatchingOctree::FindMatch

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::MatchingOctree::FillSlaveToMasterGIDMapping(
    const DRT::Discretization& slavedis,
    const std::vector<int>& slavenodeids,
    std::map<int,std::vector<double> >& coupling)
{
  CheckIsInit();
  CheckIsSetup();

  int numprocs = discret_->Comm().NumProc();

  if (slavedis.Comm().NumProc()!=numprocs)
    dserror("compared discretizations must live on same procs");

  // 1) each proc generates a list of his slavenodes
  //
  // 2) the list is communicated in a round robin pattern to all the
  //    other procs.
  //
  // 3) the proc checks the package from each proc and calcs the min
  //    distance on each --- the result is kept if distance is smaller
  //    than on the preceding processors

  //--------------------------------------------------------------------
  // -> 1) create a list of slave nodes on this proc. Pack it.
  std::vector<char> sblockofnodes;
  std::vector<char> rblockofnodes;

  DRT::PackBuffer data;

  for (unsigned globn=0; globn<slavenodeids.size(); ++globn)
  {
    PackEntity(data,&slavedis,slavenodeids[globn]);
  }

  data.StartPacking();

  for (unsigned globn=0; globn<slavenodeids.size(); ++globn)
  {
    PackEntity(data,&slavedis,slavenodeids[globn]);
  }

  swap( sblockofnodes, data() );

  //--------------------------------------------------------------------
  // -> 2) round robin loop

  // create an exporter for point to point comunication
  // We do all communication with the communicator of the original
  // discretization.
  DRT::Exporter exporter(discret_->Comm());

  for (int np=0; np<numprocs; np++)
  {
    //--------------------------------------------------
    // Send block to next proc. Receive a block from the last proc
    if (np > 0) // in the first step, we keep all nodes on this proc
    {
#ifdef PARALLEL
      int myrank   = discret_->Comm().MyPID();
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblockofnodes.size();

      exporter.ISend(frompid,topid,
                     &(sblockofnodes[0]),sblockofnodes.size(),
                     tag,request);

      // make sure that you do not think you received something if
      // you didn't
      if (not rblockofnodes.empty())
      {
        dserror("rblockofnodes not empty");
      }

      rblockofnodes.clear();

      // receive from predecessor
      frompid = (myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblockofnodes,length);

      if (tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);
#else
      dserror("How did you get here? Go away!");
#endif
    }
    else
    {
      // no need to communicate
      swap(rblockofnodes,sblockofnodes);
    }

    //--------------------------------------------------
    // Unpack block.
    std::vector<char>::size_type index = 0;
    while (index < rblockofnodes.size())
    {
      // extract node data from blockofnodes
      std::vector<char> data;
      UnPackEntity(index,rblockofnodes,data);

      // allocate an "empty node". Fill it with info from
      // extracted node data
      Teuchos::RCP<DRT::ParObject> o = Teuchos::rcp(DRT::UTILS::Factory(data));

      const int id = CheckValidEntityType(o);

      //----------------------------------------------------------------
      // there is nothing to do if there are no master nodes on this
      // proc
      if (not masterplanecoords_.empty())
      {
        double pointcoord[3];
        CalcPointCoordinate(o.getRawPtr(),&pointcoord[0]);

        // get its coordinates
        std::vector<double> x(pointcoord, pointcoord+3);

        //--------------------------------------------------------
        // 3) now search for closest master point on this proc
        int    gid;
        double dist;

        // If x is not in the bounding box on this proc, its probably not
        // matching a point in the box. We do nothing.
        if (SearchClosestEntityOnThisProc(x, gid, dist))
        {
          std::map<int,std::vector<double> >::iterator found
              = coupling.find(id);

          // search for second point with same distance,
          // if found gid is already in coupling
          if (found != coupling.end())
          {
            if (SearchClosestEntityOnThisProc(x, gid, dist, true))
            {
              found = coupling.find(id);
            }
          }

          // we are interested in the closest match
          if (found==coupling.end() or (coupling[id])[1] > dist)
          {
            if(dist <= tol_)
            {
              bool isrownode = CheckEntityOwner(discret_,gid);
              std::vector<double> myvec(3); // initialize vector
              myvec[0]=(double)gid; // save master gid in vector
              myvec[1]=dist;        // save distance in vector
              myvec[2]=(double)isrownode; // save master row col info in vector
              coupling[id] = myvec; // copy vector to map
            }
          }
        }
      }
    }

    //----------------------------------------------------------------
    // prepare to send nodes to next proc (keep list).

    // the received nodes will be sent to the next proc
    swap(sblockofnodes,rblockofnodes);

    // we need a new receive buffer
    rblockofnodes.clear();
  }
  return;
} // MatchingOctree::FillSlaveToMasterGIDMapping


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::NodeMatchingOctree::NodeMatchingOctree() :
  MatchingOctree()
{
  return;
} // NodeMatchingOctree::NodeMatchingOctree

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::NodeMatchingOctree::~NodeMatchingOctree()
{
  return;
}// NodeMatchingOctree::~NodeMatchingOctree

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::NodeMatchingOctree::CalcPointCoordinate(
    const DRT::Discretization* dis,
    const int id,
    double* coord)
{
  DRT::Node* actnode = dis->gNode(id);

  const int dim = 3;

  for(int idim=0;idim<dim;idim++)
    coord[idim]=actnode->X()[idim];

  return;
} // NodeMatchingOctree::CalcPointCoordinate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//! calc unique coordinate of entity
void  DRT::UTILS::NodeMatchingOctree::CalcPointCoordinate(
    DRT::ParObject* entity,
    double* coord)
{
  DRT::Node* actnode = dynamic_cast<DRT::Node*>(entity);
  if (actnode==NULL)
    dserror("dynamic_cast failed");

  const int dim = 3;

  for(int idim=0;idim<dim;idim++)
    coord[idim]=actnode->X()[idim];

  return;
} // NodeMatchingOctree::CalcPointCoordinate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::NodeMatchingOctree::CheckHaveEntity(
    const DRT::Discretization* dis,
    const int id)
{
  return dis->HaveGlobalNode(id);
} // NodeMatchingOctree::CheckHaveEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::NodeMatchingOctree::CheckEntityOwner(
    const DRT::Discretization* dis,
    const int id)
{
  return (dis->gNode(id)->Owner() == dis->Comm().MyPID());
} // NodeMatchingOctree::CheckEntityOwner

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::NodeMatchingOctree::PackEntity(
    PackBuffer& data,
    const DRT::Discretization* dis,
    const int id)
{
  // get the slavenode
  DRT::Node* actnode = dis->gNode(id);
  // Add node to list of nodes which will be sent to the next proc
  DRT::ParObject::AddtoPack(data,actnode);

  return;
} // NodeMatchingOctree::PackEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::NodeMatchingOctree::UnPackEntity(
    std::vector<char>::size_type& index,
    std::vector<char>& rblockofnodes,
    std::vector<char>& data)
{
  DRT::ParObject::ExtractfromPack(index,rblockofnodes,data);
  return;
} // NodeMatchingOctree::UnPackEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::UTILS::NodeMatchingOctree::CheckValidEntityType(Teuchos::RCP<DRT::ParObject> o)
{
  // cast ParObject to Node
  DRT::Node* actnode = dynamic_cast<DRT::Node*>(o.get());
  if (actnode==NULL)
    dserror("unpack of invalid data");

  return actnode->Id();
} // NodeMatchingOctree::CheckValidEntityType

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::OctreeElement>
DRT::UTILS::NodeMatchingOctree::CreateOctreeElement(
    std::vector <int> &              nodeidstoadd,
    Epetra_SerialDenseMatrix&        boundingboxtoadd,
    int                              layer)
{
  Teuchos::RCP<DRT::UTILS::OctreeElement> newtreeelement =
      Teuchos::rcp(new OctreeNodalElement());

  newtreeelement->Init(
      *discret_,
      nodeidstoadd,
      boundingboxtoadd,
      layer,
      maxtreenodesperleaf_,
      tol_);

  newtreeelement->Setup();

  return newtreeelement;
}// NodeMatchingOctree::CreateOctreeElement


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ElementMatchingOctree::ElementMatchingOctree() :
  MatchingOctree()
{
  return;
} // ElementMatchingOctree::ElementMatchingOctree

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ElementMatchingOctree::~ElementMatchingOctree()
{
  return;
}// ElementMatchingOctree::~ElementMatchingOctree

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ElementMatchingOctree::CalcPointCoordinate(
    const DRT::Discretization* dis,
    const int id,
    double* coord)
{
  DRT::Element* actele = dis->gElement(id);

  const int numnode = actele->NumNode();
  const int dim = 3;

  for(int idim=0;idim<dim;idim++)
    coord[idim]=0.0;

  for(int node=0;node<numnode;node++)
    for(int idim=0;idim<dim;idim++)
      coord[idim]+=(actele->Nodes())[node]->X()[idim];

  return;
} // ElementMatchingOctree::CalcPointCoordinate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void  DRT::UTILS::ElementMatchingOctree::CalcPointCoordinate(
    DRT::ParObject* entity,
    double* coord)
{
  DRT::Element* actele = dynamic_cast<DRT::Element*>(entity);
  if (actele==NULL)
    dserror("dynamic_cast failed");

  DRT::Node** nodes = actele->Nodes();
  if(nodes==NULL)
    dserror("could not get pointer to nodes");

  const int numnode = actele->NumNode();
  const int dim = 3;

  for(int idim=0;idim<dim;idim++)
    coord[idim]=0.0;

  for(int node=0;node<numnode;node++)
    for(int idim=0;idim<dim;idim++)
      coord[idim]+=(nodes[node]->X())[idim];

  return;
} // ElementMatchingOctree::CalcPointCoordinate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::ElementMatchingOctree::CheckHaveEntity(
    const DRT::Discretization* dis,
    const int id)
{
  return dis->HaveGlobalElement(id);
} // ElementMatchingOctree::CheckHaveEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::ElementMatchingOctree::CheckEntityOwner(
    const DRT::Discretization* dis,
    const int id)
{
  return (dis->gElement(id)->Owner() == dis->Comm().MyPID());
} // ElementMatchingOctree::CheckHaveEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ElementMatchingOctree::PackEntity(
    PackBuffer& data,
    const DRT::Discretization* dis,
    const int id)
{
  // get the slavenode
  DRT::Element* actele = dis->gElement(id);
  DRT::Node** nodes = actele->Nodes();
  // Add node to list of nodes which will be sent to the next proc
  DRT::ParObject::AddtoPack(data,actele->NumNode());
  DRT::ParObject::AddtoPack(data,actele);
  for(int node=0;node<actele->NumNode();node++)
    DRT::ParObject::AddtoPack(data,nodes[node]);

  return;
} // ElementMatchingOctree::PackEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::ElementMatchingOctree::UnPackEntity(
    std::vector<char>::size_type& index,
    std::vector<char>& rblockofnodes,
    std::vector<char>& data)
{
  nodes_.clear();
  int numnode = DRT::ParObject::ExtractInt(index,rblockofnodes);
  DRT::ParObject::ExtractfromPack(index,rblockofnodes,data);

  for(int node=0;node<numnode;node++)
  {
    std::vector<char> nodedata;
    DRT::ParObject::ExtractfromPack(index,rblockofnodes,nodedata);
    Teuchos::RCP<DRT::ParObject> o = Teuchos::rcp(DRT::UTILS::Factory(nodedata));
    Teuchos::RCP<DRT::Node> actnode = Teuchos::rcp_dynamic_cast<DRT::Node >(o);
    if(actnode==Teuchos::null)
      dserror("cast from ParObject to Node failed");
    nodes_.insert(std::pair<int,Teuchos::RCP<DRT::Node> >(actnode->Id(),actnode));
  }

  return;
} // ElementMatchingOctree::UnPackEntity

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::UTILS::ElementMatchingOctree::CheckValidEntityType(Teuchos::RCP<DRT::ParObject> o)
{
  // cast ParObject to element
  DRT::Element* actele = dynamic_cast<DRT::Element*>(o.get());
  if (actele==NULL)
    dserror("unpack of invalid data");

  // set nodal pointers for this element
  actele->BuildNodalPointers(nodes_);

  return actele->Id();
} // ElementMatchingOctree::CheckValidEntityType

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::OctreeElement>
DRT::UTILS::ElementMatchingOctree::CreateOctreeElement(
    std::vector <int> &              nodeidstoadd,
    Epetra_SerialDenseMatrix&        boundingboxtoadd,
    int                              layer)
{
  Teuchos::RCP<DRT::UTILS::OctreeElement> newtreeelement =
      Teuchos::rcp(new OctreeElementElement());

  newtreeelement->Init(
      *discret_,
      nodeidstoadd,
      boundingboxtoadd,
      layer,
      maxtreenodesperleaf_,
      tol_);

  newtreeelement->Setup();

  return newtreeelement;
}// ElementMatchingOctree::CreateOctreeElement


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::OctreeNodalElement::OctreeNodalElement() :
  OctreeElement()
{
  return;
} //OctreeElement()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::OctreeNodalElement::CalcPointCoordinate(
    const DRT::Discretization* dis,
    const int id,
    double* coord)
{
  DRT::Node* actnode = dis->gNode(id);

  const int dim = 3;

  for(int idim=0;idim<dim;idim++)
    coord[idim]=actnode->X()[idim];

  return;
} // OctreeNodalElement::CalcPointCoordinate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::OctreeElement>
DRT::UTILS::OctreeNodalElement::CreateOctreeElement(
    std::vector <int> &              nodeidstoadd,
    Epetra_SerialDenseMatrix&        boundingboxtoadd,
    int                              layer)
{
  Teuchos::RCP<DRT::UTILS::OctreeElement> newtreeelement =
      Teuchos::rcp(new OctreeNodalElement());

  newtreeelement->Init(
      *discret_,
      nodeidstoadd,
      boundingboxtoadd,
      layer,
      maxtreenodesperleaf_,
      tol_);

  newtreeelement->Setup();

  return newtreeelement;
}// OctreeNodalElement::CreateOctreeElement


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::OctreeElementElement::OctreeElementElement() :
  OctreeElement()
{
  return;
}// OctreeElementElement::OctreeElementElement

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::OctreeElementElement::CalcPointCoordinate(
    const DRT::Discretization* dis,
    const int id,
    double* coord)
{
  DRT::Element* actele = dis->gElement(id);

  const int numnode = actele->NumNode();
  const int dim = 3;

  for(int idim=0;idim<dim;idim++)
    coord[idim]=0.0;

  for(int node=0;node<numnode;node++)
    for(int idim=0;idim<dim;idim++)
      coord[idim]+=(actele->Nodes())[node]->X()[idim];

  return;
} // OctreeElementElement::CalcPointCoordinate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::OctreeElement>
DRT::UTILS::OctreeElementElement::CreateOctreeElement(
    std::vector <int> &              nodeidstoadd,
    Epetra_SerialDenseMatrix&        boundingboxtoadd,
    int                              layer)
{
  Teuchos::RCP<DRT::UTILS::OctreeElement> newtreeelement =
      Teuchos::rcp(new OctreeElementElement());

  newtreeelement->Init(
      *discret_,
      nodeidstoadd,
      boundingboxtoadd,
      layer,
      maxtreenodesperleaf_,
      tol_);

  newtreeelement->Setup();

  return newtreeelement;
}// OctreeElementElement::CreateOctreeElement


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::OctreeElement::OctreeElement() :
  discret_(NULL),
  layer_(-1),
  maxtreenodesperleaf_(-1),
  tol_(-1.0),
  issetup_(false),
  isinit_(false)
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::UTILS::OctreeElement::Init(
    const DRT::Discretization&       actdis,
    std::vector <int>&               nodeidstoadd,
    const Epetra_SerialDenseMatrix&  boundingboxtoadd,
    const int                        layer,
    const int                        maxnodeperleaf,
    const double                     tol)
{
  SetIsSetup(false);

  discret_=&actdis;
  boundingbox_=boundingboxtoadd;
  nodeids_=nodeidstoadd;
  layer_=layer;
  maxtreenodesperleaf_=maxnodeperleaf;
  tol_=tol;

  SetIsInit(true);
  return 0;
} // OctreeElement::Init

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::UTILS::OctreeElement::Setup()
{
  CheckIsInit();

  if(layer_>200)
  {
    dserror("max.depth of octree: 200. Can't append further children\n");
  }

  int numnodestoadd = nodeids_.size();
  // if number of slavenodes on this proc is too large split the element
  if (numnodestoadd>maxtreenodesperleaf_)
  {
    // mean coordinate value in direction of the longest edge
    double mean[3];
    double pointcoord[3];

    for(int dim=0;dim<3;dim++)
    {
      mean[dim]=0.0;
      pointcoord[dim]=0.0;
    }

    // calculate mean coordinate for all directions
    for(int locn=0;locn<numnodestoadd;locn++)
    {
      CalcPointCoordinate(discret_,nodeids_.at(locn),&pointcoord[0]);

      for(int dim=0;dim<3;dim++)
      {
        mean[dim]+=pointcoord[dim];
      }
    }

    for(int dim=0;dim<3;dim++)
    {
      mean[dim]=mean[dim]/numnodestoadd;
    }

    // direction specifies which side will be cut (the one with the largest
    // value of the mean distance to the "lower" boundary)
    int direction=0;


    // the maximum distance of the mean value from the boundary of the box
    double maxdist   =0;
    double wheretocut=0;

    // get the coordinate with the maximum distance
    for(int dim=0;dim<3;dim++)
    {
      double thisdist= std::min(mean[dim]-boundingbox_(dim,0),
                            boundingbox_(dim,1)-mean[dim]);

      if(maxdist<thisdist)
      {
        maxdist   = thisdist;
        wheretocut= mean[dim];
        direction = dim;
      }
    }

    // Why choose the coordinate with the maximum distance from both edges
    // and not simply the longest edge?
    //
    // Here is what happens if you do so:
    //
    //
    //  +-------------+
    //  |             |
    //  X             |
    //  |             |     (*)
    //  X             |
    //  |             |
    //  |             |
    //  +-------------+
    //   longest edge
    //
    // This leads to two children:
    //
    //
    // +--------------+
    // |              |
    // |X             |
    // |              |
    // |X             |
    // |              |
    // |              |
    // +--------------+
    //
    // and
    //
    // +-+
    // | |
    // |X|
    // | |
    // |X|
    // | |
    // | |
    // +-+
    //
    // This means an endless loop from the problem (*) ...

    // create bounding boxes for the children
    Epetra_SerialDenseMatrix   childboundingbox1(3,2);
    Epetra_SerialDenseMatrix   childboundingbox2(3,2);

    // initialise the new bounding boxes with the old values
    for (int dim=0;dim<3;dim++)
    {
      childboundingbox1(dim,0)=boundingbox_(dim,0);
      childboundingbox1(dim,1)=boundingbox_(dim,1);

      childboundingbox2(dim,0)=boundingbox_(dim,0);
      childboundingbox2(dim,1)=boundingbox_(dim,1);
    }

    // replace the boundaries of component direction to divide cells
    // create initial bounding box for all nodes
    //
    // example direction = 1 , e. g. y-direction:
    //
    // mean = ymin + maxdist
    //
    //   +-            -+      +-               -+   +-               -+
    //   |  xmin  xmax  |      |  xmin  xmax     |   |  xmin     xmax  |
    //   |  ymin  ymax  | ---> |  ymin  mean+eps | + |  mean-eps ymax  |
    //   |  zmin  zmax  |      |  zmin  zmax     |   |  zmin     zmax  |
    //   +-            -+      +-               -+   +-               -+
    //
    //                         lower bounding box    upper bounding box
    //


    childboundingbox1(direction,1)=wheretocut+tol_;
    childboundingbox2(direction,0)=wheretocut-tol_;

    // distribute nodes to children
    std::vector <int> childnodeids1;
    std::vector <int> childnodeids2;
    for(int locn=0;locn<(int)nodeids_.size();locn++)
    {
      CalcPointCoordinate(discret_,nodeids_.at(locn),&pointcoord[0]);

      // node is in "lower" bounding box
      if(pointcoord[direction]<childboundingbox1(direction,1))
      {
        childnodeids1.push_back(nodeids_.at(locn));
      }
      // node is in "upper" bounding box
      if(pointcoord[direction]>childboundingbox2(direction,0))
      {
        childnodeids2.push_back(nodeids_.at(locn));
      }
    }

    // we do not need the full node id vector anymore --- it was distributed
    // to the children --> throw it away
    nodeids_.clear();

    // append children to parent
    octreechild1_ = CreateOctreeElement(childnodeids1,childboundingbox1,layer_+1);

    octreechild2_ = CreateOctreeElement(childnodeids2,childboundingbox2,layer_+1);

  } // end number of slavenodes on this proc is too large split the element
  else
  {
    if ((int)nodeids_.size() == 0)
    {
      dserror("Trying to create leaf with no nodes. Stop.");
    }
  }

  SetIsSetup(true);
  return 0;
} // OctreeElement::Setup

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::OctreeElement::SearchClosestNodeInLeaf(
  const std::vector <double> & x,
  int             & idofclosestpoint,
  double          & distofclosestpoint,
  const double    & elesize,
  bool              searchsecond
  )
{
  CheckIsInit();
  CheckIsSetup();

  double thisdist;
  std::vector <double> dx(3);
  double pointcoord[3];

  // the first node is the guess for the closest node
  CalcPointCoordinate(discret_,nodeids_.at(0),&pointcoord[0]);
  for (int dim=0;dim<3;dim++)
  {
    dx[dim]=pointcoord[dim]-x[dim];
  }

  distofclosestpoint = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  idofclosestpoint   = nodeids_.at(0);

  // now loop the others and check whether they are better
  for (int nn=1;nn<(int)nodeids_.size();nn++)
  {
    CalcPointCoordinate(discret_,nodeids_.at(nn),&pointcoord[0]);

    for (int dim=0;dim<3;dim++)
    {
      dx[dim]=pointcoord[dim]-x[dim];
    }
    thisdist = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if (thisdist < (distofclosestpoint - 1e-02*elesize))
    {
      distofclosestpoint = thisdist;
      idofclosestpoint = nodeids_.at(nn);
    }
    else
    {
      if ((abs(thisdist - distofclosestpoint) < 1e-02*elesize) & (searchsecond == true))
      {
        distofclosestpoint = thisdist;
        idofclosestpoint = nodeids_.at(nn);
      }
    }
  }

  return;
}// OctreeElement::SearchClosestNodeInLeaf

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::OctreeElement::IsPointInBoundingBox(
  const std::vector <double> &x
  )
{
  CheckIsInit();
  CheckIsSetup();

  bool nodeinboundingbox=true;

  for (int dim=0;dim<3;dim++)
  {
    // check whether node is outside of bounding box
    if((x[dim]<this->boundingbox_(dim,0))
       ||
       (x[dim]>this->boundingbox_(dim,1)))
    {
      nodeinboundingbox=false;
    }
  }

  return nodeinboundingbox;
} //OctreeElement::IsPointInBoundingBox

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::OctreeElement>
DRT::UTILS::OctreeElement::ReturnChildContainingPoint(
  const std::vector <double> &x
  )
{
  CheckIsInit();
  CheckIsSetup();

  Teuchos::RCP<OctreeElement> nextelement;

  if (this->octreechild1_ == Teuchos::null || this->octreechild2_ == Teuchos::null)
  {
    dserror("Asked leaf element for further children.");
  }

  if(this->octreechild1_->IsPointInBoundingBox(x) == true)
  {

    nextelement=this->octreechild1_;
  }
  else if(this->octreechild2_->IsPointInBoundingBox(x) == true)
  {
    nextelement=this->octreechild2_;
  }
  else
  {
    dserror("point in no bounding box of children, but in parent bounding box!");
  }

  return nextelement;
} //OctreeElement::ReturnChildContainingPoint

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::OctreeElement::IsLeaf()
{
  bool isleaf=true;

  if (this->nodeids_.size()==0)
  {
    isleaf=false;
  }

  return isleaf;
} //OctreeElement::IsLeaf

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::OctreeElement::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Leaf in Layer " << layer_ << " Nodes ";

  for (int nn=0; nn<(int)nodeids_.size(); ++nn)
  {
    os << nodeids_.at(nn) << " ";
  }
  os << std::endl;
  return;
} // OctreeElement::Print(ostream& os)

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::OctreeElement::~OctreeElement()
{
  return;
}// ~OctreeElement()
