/*!----------------------------------------------------------------------
\file drt_nodematchingoctree.cpp

\brief Search closest node in given set of nodes using an octree search

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/


#include "drt_nodematchingoctree.H"
#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
\brief Constructor (public)

<pre>

 Set up processor local octree                               gammi 04/07

<pre>

\param    Teuchos::RCP<DRT::Discretization> (i) discretisation
\param    const vector <int> &             (i) list of masternodeids
\param    int                              (i) parameter for octree
\param    double                           (i) tolerance for octree

\return void

 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

DRT::UTILS::NodeMatchingOctree::NodeMatchingOctree(
  const DRT::Discretization&       actdis,
  const std::vector <int> &        masternodeids,
  int                              maxnodeperleaf,
  double                           tol
  ):
  // call constructor for "nontrivial" objects
  discret_(actdis),
  tol_(tol)
{
  // extract all masternodes on this proc from the list masternodeids
  std::vector <int> masternodesonthisproc;

  for(int locn=0;locn<(int)masternodeids.size();locn++)
  {
    // if node is on this proc
    if(discret_.HaveGlobalNode(masternodeids[locn]))
    {
      // if node is not ghosted
      if (discret_.gNode(masternodeids[locn])->Owner() == discret_.Comm().MyPID())
      {
        // this masternode is on this proc and is not a ghosted one
        masternodesonthisproc.push_back(masternodeids[locn]);
      }
    }
  }

  // construct octree if proc has masternodes
  if(masternodesonthisproc.size()>0)
  {
    // create initial bounding box for all nodes
    //
    //                 +-            -+
    //                 |  xmin  xmax  |
    //                 |  ymin  ymax  |
    //                 |  zmin  zmax  |
    //                 +-            -+
    //
    Epetra_SerialDenseMatrix   initialboundingbox(3,2);

    DRT::Node* actnode = discret_.gNode(masternodesonthisproc[0]);
    for (int dim=0;dim<3;dim++)
    {
      initialboundingbox(dim,0)=actnode->X()[dim]-tol_;
      initialboundingbox(dim,1)=actnode->X()[dim]+tol_;

      // store coordinates of one point in master plane (later on, one
      // coordinate of the masternode will be substituted by the coordinate
      // of the master plane)
      masterplanecoords_.push_back(actnode->X()[dim]);
    }

    for(unsigned locn=0;locn<masternodesonthisproc.size();locn++)
    {
      DRT::Node* actnode = discret_.gNode(masternodesonthisproc[locn]);
      for (int dim=0;dim<3;dim++)
      {
        initialboundingbox(dim,0)=std::min(initialboundingbox(dim,0),
                                       actnode->X()[dim]-tol_);
        initialboundingbox(dim,1)=std::max(initialboundingbox(dim,1),
                                       actnode->X()[dim]+tol_);
      }
    }

    // create octree root --- initial layer is 0
    // all other layers are generated down here by recursive calls
    int initlayer = 0;

    octreeroot_ = Teuchos::rcp(new OctreeElement(discret_,
                                                       masternodesonthisproc,
                                                       initialboundingbox,
                                                       initlayer,
                                                       maxnodeperleaf,
                                                       tol_));
  }

  return;
} // NodeMatchingOctree::NodeMatchingOctree

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Search for closest (slave) nodes on all processors to given (master) |
 | nodeset (only in bounding box of masternodes)                        |
 |                                                           gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::UTILS::NodeMatchingOctree::CreateGlobalNodeMatching(
  const std::vector<int>    &     slavenodeids,
  const std::vector<int>    &     dofsforpbcplane,
  const double               rotangle,
  std::map<int,std::vector<int> >&     midtosid
  )
{
  int myrank  =discret_.Comm().MyPID();
  int numprocs=discret_.Comm().NumProc();

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
    if(discret_.HaveGlobalNode(slavenodeids[globn]))
    {
      // get the slavenode
      DRT::Node* actnode = discret_.gNode(slavenodeids[globn]);

      // take only nodes which are not ghosted
      if (actnode->Owner() == myrank)
      {
        // Add node to list of nodes which will be sent to the next proc
        actnode->Pack(data);
      }
    } // end if slavenode on proc
  } // end loop globn

  data.StartPacking();

  for(int globn=0;globn<(int)slavenodeids.size();globn++)
  {
    // is this slavenode on this proc?
    if(discret_.HaveGlobalNode(slavenodeids[globn]))
    {
      // get the slavenode
      DRT::Node* actnode = discret_.gNode(slavenodeids[globn]);

      // take only nodes which are not ghosted
      if (actnode->Owner() == myrank)
      {
        // Add node to list of nodes which will be sent to the next proc
        actnode->Pack(data);
      }
    } // end if slavenode on proc
  } // end loop globn

  swap( sblockofnodes, data() );

  //--------------------------------------------------------------------
  // -> 2) round robin loop
#ifdef PARALLEL
  // create an exporter for point to point comunication
  DRT::Exporter exporter(discret_.Comm());
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

      // cast ParObject to Node
      DRT::Node* actnode = dynamic_cast<DRT::Node*>(o.get());

      //----------------------------------------------------------------
      // there is nothing to do if there are no master nodes on this
      // proc
      if(masterplanecoords_.empty()!=true)
      {
        // get its coordinates
        std::vector <double> x(3);

        if (abs(rotangle) < EPS13)
        {
          for (int dim=0;dim<3;dim++)
          {
            x[dim] = actnode->X()[dim];
          }
        }
        else
        {
          // if there is a rotationally symmetric periodic boundary condition:
          // rotate slave plane for making it parallel to the master plane
          x[0] = actnode->X()[0]*cos(rotangle) + actnode->X()[1]*sin(rotangle);
          x[1] = actnode->X()[0]*(-sin(rotangle)) + actnode->X()[1]*cos(rotangle);
          x[2] = actnode->X()[2];
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

        nodeisinbox=this->SearchClosestNodeOnThisProc(
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
        (midtosid[idofclosestpoint]).push_back(actnode->Id());
              diststom[idofclosestpoint] = distofclosestpoint;
            }
      else if(diststom[idofclosestpoint]<distofclosestpoint+1e-9
        &&
        diststom[idofclosestpoint]>distofclosestpoint-1e-9
        )
      {
        (midtosid[idofclosestpoint]).push_back(actnode->Id());
      }
          }
          else
          {
            // this is the first estimate for a closest point
      (midtosid[idofclosestpoint]).clear();
      (midtosid[idofclosestpoint]).push_back(actnode->Id());
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
} // CreateGlobalNodeMatching



void DRT::UTILS::NodeMatchingOctree::FindMatch(const DRT::Discretization& slavedis,
                                   const std::vector<int>& slavenodeids,
                                   std::map<int,std::pair<int,double> >& coupling)
{
  int numprocs = discret_.Comm().NumProc();

  int slaverank = slavedis.Comm().MyPID();
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
    if (slavedis.HaveGlobalNode(slavenodeids[globn]))
    {
      // get the slavenode
      DRT::Node* actnode = slavedis.gNode(slavenodeids[globn]);

      // take only nodes which are not ghosted
      if (actnode->Owner() == slaverank)
      {
        // Add node to list of nodes which will be sent to the next proc
        DRT::ParObject::AddtoPack(data,actnode);
      }
    }
  }

  data.StartPacking();

  for (unsigned globn=0; globn<slavenodeids.size(); ++globn)
  {
    // is this slavenode on this proc?
    if (slavedis.HaveGlobalNode(slavenodeids[globn]))
    {
      // get the slavenode
      DRT::Node* actnode = slavedis.gNode(slavenodeids[globn]);

      // take only nodes which are not ghosted
      if (actnode->Owner() == slaverank)
      {
        // Add node to list of nodes which will be sent to the next proc
        DRT::ParObject::AddtoPack(data,actnode);
      }
    }
  }

  swap( sblockofnodes, data() );

  //--------------------------------------------------------------------
  // -> 2) round robin loop

  // create an exporter for point to point comunication
  // We do all communication with the communicator of the original
  // discretization.
  DRT::Exporter exporter(discret_.Comm());

  for (int np=0; np<numprocs; np++)
  {
    //--------------------------------------------------
    // Send block to next proc. Receive a block from the last proc
    if (np > 0) // in the first step, we keep all nodes on this proc
    {
#ifdef PARALLEL
      int myrank   = discret_.Comm().MyPID();
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

      // cast ParObject to Node
      DRT::Node* actnode = dynamic_cast<DRT::Node*>(o.get());
      if (actnode==NULL)
        dserror("unpack of invalid data");

      //----------------------------------------------------------------
      // there is nothing to do if there are no master nodes on this
      // proc
      if (not masterplanecoords_.empty())
      {
        // get its coordinates
        std::vector<double> x(actnode->X(), actnode->X()+3);

        //--------------------------------------------------------
        // 3) now search for closest master point on this proc
        int    gid;
        double dist;

        // If x is not in the bounding box on this proc, its probably not
        // matching a point in the box. We do nothing.
        if (SearchClosestNodeOnThisProc(x, gid, dist))
        {
          std::map<int,std::pair<int,double> >::iterator found = coupling.find(gid);

          // search for second point with same distance, if found gid is already in coupling
          if (found != coupling.end())
          {
            if (SearchClosestNodeOnThisProc(x, gid, dist, true))
            {
              found = coupling.find(gid);
            }
          }

          // we are interested in the closest match
          if (found==coupling.end() or coupling[gid].second > dist)
          {
            coupling[gid] = std::make_pair(actnode->Id(), dist);
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
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Search closest node to given node in local octree         gammi 05/07|
 | returns false if node is not in bounding box of local octree         |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
bool DRT::UTILS::NodeMatchingOctree::SearchClosestNodeOnThisProc(
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
    octreeele->SearchClosestNodeInLeaf(x,idofclosestpoint,distofclosestpoint,tol_,searchsecond);
  }

  return nodeisinbox;
}// NodeMatchingOctree::SearchClosestNodeOnThisProc

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor (public)                                       gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
DRT::UTILS::NodeMatchingOctree::~NodeMatchingOctree()
{
  return;
}// NodeMatchingOctree::~NodeMatchingOctree



//======================================================================



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Constructor (public)                                                 |
 | Create one element in octree                              gammi 05/07|
  *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

DRT::UTILS::OctreeElement::OctreeElement(
  const DRT::Discretization&       actdis,
  std::vector <int> &              nodeidstoadd,
  Epetra_SerialDenseMatrix&        boundingboxtoadd,
  int                              layer,
  int                              maxnodeperleaf,
  double                           tol
  ):
// call constructor for "nontrivial" objects
  discret_(actdis)
{

  boundingbox_ = boundingboxtoadd;
  layer_       = layer;

  if(layer_>200)
  {
    dserror("max.depth of octree: 200. Can't append further children\n");
  }

  int numnodestoadd = nodeidstoadd.size();
  // if number of slavenodes on this proc is too large split the element
  if (numnodestoadd>maxnodeperleaf)
  {
    // mean coordinate value in direction of the longest edge
    double mean[3];

    for(int dim=0;dim<3;dim++)
    {
      mean[dim]=0;
    }

    // calculate mean coordinate for all directions
    for(int locn=0;locn<numnodestoadd;locn++)
    {
      DRT::Node* actnode = discret_.gNode(nodeidstoadd[locn]);

      for(int dim=0;dim<3;dim++)
      {
        mean[dim]+=actnode->X()[dim];
      }
    }

    for(int dim=0;dim<3;dim++)
    {
      mean[dim]=mean[dim]/numnodestoadd;
    }

    // direction specifies which side will be cut (the one with the largest
    // value of the mean distance to the "lower" boundary)
    int    direction=0;


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


    childboundingbox1(direction,1)=wheretocut+tol;
    childboundingbox2(direction,0)=wheretocut-tol;

    // distribute nodes to children
    std::vector <int> childnodeids1;
    std::vector <int> childnodeids2;
    for(int locn=0;locn<(int)nodeidstoadd.size();locn++)
    {
      double coordinate=discret_.gNode(nodeidstoadd[locn])->X()[direction];

      // node is in "lower" bounding box
      if(coordinate<childboundingbox1(direction,1))
      {
        childnodeids1.push_back(nodeidstoadd[locn]);
      }
      // node is in "upper" bounding box
      if(coordinate>childboundingbox2(direction,0))
      {
        childnodeids2.push_back(nodeidstoadd[locn]);
      }
    }

    // we do not need the full node id vector anymore --- it was distributed
    // to the children --> throw it away
    nodeidstoadd.clear();

    // append children to parent
    octreechild1_ = Teuchos::rcp(new OctreeElement(discret_,
                                                         childnodeids1,
                                                         childboundingbox1,
                                                         layer_+1,
                                                         maxnodeperleaf,
                                                         tol));
    octreechild2_ = Teuchos::rcp(new OctreeElement(discret_,
                                                         childnodeids2,
                                                         childboundingbox2,
                                                         layer_+1,
                                                         maxnodeperleaf,
                                                         tol));
  } // end number of slavenodes on this proc is too large split the element
  else
  {
    if ((int)nodeidstoadd.size() == 0)
    {
      dserror("Trying to create leaf with no nodes. Stop.");
    }

    // we have a leave element
    nodeids_     = nodeidstoadd;
  }
  return;
} //OctreeElement()

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Check if a point is in the bounding box of an element (public)       |
 |                                                           gammi 05/07|
  *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

bool DRT::UTILS::OctreeElement::IsPointInBoundingBox(
  const std::vector <double> &x
  )
{
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

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Return a child containing the node                    (public)       |
 |                                                           gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

Teuchos::RCP<DRT::UTILS::OctreeElement> DRT::UTILS::OctreeElement::ReturnChildContainingPoint(
  const std::vector <double> &x
  )
{
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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Check if a point is in the bounding box of an element (public)       |
 |                                                           gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

bool DRT::UTILS::OctreeElement::IsLeaf()
{
  bool isleaf=true;

  if (this->nodeids_.size()==0)
  {
    isleaf=false;
  }

  return isleaf;
} //OctreeElement::IsLeaf


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Print some information on the octree element (public)                |
 |                                                           gammi 05/07|
  *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void DRT::UTILS::OctreeElement::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Leaf in Layer " << layer_ << " Nodes ";

  for (int nn=0; nn<(int)nodeids_.size(); ++nn)
  {
    os << nodeids_[nn] << " ";
  }
  os << std::endl;
  return;
} // OctreeElement::Print(ostream& os)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Return closest point in leaf                 (public)                |
 |                                                           gammi 05/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::UTILS::OctreeElement::SearchClosestNodeInLeaf(
  const std::vector <double> & x,
  int             & idofclosestpoint,
  double          & distofclosestpoint,
  const double    & elesize,
  bool              searchsecond
  )
{
  double          thisdist;
  std::vector <double> dx(3);

  // the first node is the guess for the closest node
  DRT::Node* actnode = discret_.gNode(this->nodeids_[0]);
  for (int dim=0;dim<3;dim++)
  {
    dx[dim]=actnode->X()[dim]-x[dim];
  }

  distofclosestpoint = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  idofclosestpoint   = this->nodeids_[0];

  // now loop the others and check whether they are better
  for (int nn=1;nn<(int)this->nodeids_.size();nn++)
  {
    actnode = discret_.gNode(this->nodeids_[nn]);

    for (int dim=0;dim<3;dim++)
    {
      dx[dim]=actnode->X()[dim]-x[dim];
    }
    thisdist = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

    if (thisdist < (distofclosestpoint - 1e-02*elesize))
    {
      distofclosestpoint = thisdist;
      idofclosestpoint = this->nodeids_[nn];
    }
    else
    {
      if ((abs(thisdist - distofclosestpoint) < 1e-02*elesize) & (searchsecond == true))
      {
        distofclosestpoint = thisdist;
        idofclosestpoint = this->nodeids_[nn];
      }
    }
  }


  return;
}// OctreeElement::SearchClosestNodeInLeaf

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 04/07|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
DRT::UTILS::OctreeElement::~OctreeElement()
{
  return;
}// ~OctreeElement()
