/*!----------------------------------------------------------------------
\file drt_dofset_transparent.cpp

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_dofset_transparent.H"
#include "../linalg/linalg_utils.H"

DRT::TransparentDofSet::TransparentDofSet(RCP<DRT::Discretization> sourcedis,bool parallel) :
DRT::DofSet(),
parallel_(parallel)
{
  sourcedis_.push_back(sourcedis);
  return;
}

DRT::TransparentDofSet::TransparentDofSet(vector<RCP<DRT::Discretization> > sourcedis,bool parallel) :
DRT::DofSet(),
sourcedis_(sourcedis),
parallel_(parallel)
{
  return;
}

int DRT::TransparentDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{

  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);

  if(!parallel_)
  {
    TransferDegreesOfFreedom(sourcedis_, dis, start);
  }
  else
  {
    ParallelTransferDegreesOfFreedom(sourcedis_, dis, start);
  }

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

/// Assign dof numbers for new discretization using dof numbering from source discretization.
void DRT::TransparentDofSet::TransferDegreesOfFreedom(
        const vector<RCP<DRT::Discretization> > sourcedis,
        const DRT::Discretization& newdis,
        const int start
        )
{
  for (int i=0; i<int(sourcedis.size()); i++)
  {
    if (!sourcedis[i]->DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!sourcedis[i]->NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!sourcedis[i]->ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");
  }

    if (!newdis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!newdis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!newdis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");

    // build dofrowmap based on dofrowvec which contains the dofs of each node from the (multiple) source dis;
    // dofs of the first node are followed by dofs of the second node and so on
    vector<int> dofrowvec;
    dofrowvec.reserve(dofrowmap_->NumMyElements());
    for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
    {
      const DRT::Node* newnode = newdis.lRowNode(inode);

      // loop over all source discretizations to fill each node with dofs from the different source discrets
      vector<int> dofs;
      for (int i=0; i<int(sourcedis.size()); i++)
      {
        const DRT::Node* sourcenode = sourcedis[i]->gNode(newnode->Id());
        vector<int> ldofs = sourcedis[i]->Dof(sourcenode);
        for(int l=0; l<int(ldofs.size()); l++)
        {
          dofs.push_back(ldofs[l]);
        }
      }

      const int newlid = newnode->LID();
      const int numdofs = (*numdfcolnodes_)[newlid];
      if (numdofs > 0)
      {
        if(numdofs > (int)dofs.size())
          dserror("numdofs %d > %d for node %d",numdofs,(int)dofs.size(),newnode->Id());
        (*idxcolnodes_)[newlid] = dofs[0];
        for (int idof = 0; idof < numdofs; ++idof)
        {
          dofrowvec.push_back(dofs[idof]);
        }
      }
    }

    dofrowmap_ = rcp(new Epetra_Map(-1, dofrowvec.size(), &dofrowvec[0], 0, newdis.Comm()));

    // build dofcolmap based on dofcolvec which contains the dofs of each node from the (multiple) source dis;
    // dofs of the first node are followed by dofs of the second node and so on
    vector<int> dofcolvec;
    dofcolvec.reserve(dofcolmap_->NumMyElements());
    for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
    {
      const DRT::Node* newnode = newdis.lColNode(inode);

      // loop over all source discretizations to fill each node with dofs from the different source discrets
      vector<int> dofs;
      for (int i=0; i<int(sourcedis.size()); i++)
      {
        const DRT::Node* sourcenode = sourcedis[i]->gNode(newnode->Id());

        const int lid = sourcenode->LID();
        if (lid==-1)
        {
          dserror("required node %d not on proc",newnode->Id());
        }
        vector<int> ldofs = sourcedis[i]->Dof(sourcenode);
        for(int l=0; l<int(ldofs.size()); l++)
          dofs.push_back(ldofs[l]);
      }

      const int newlid = newnode->LID();
      //const int newfirstidx = (*idxcolnodes_)[newlid];
      const int numdofs = (*numdfcolnodes_)[newlid];
      if (numdofs > 0)
      {
        if(numdofs > (int)dofs.size())
          dserror("numdofs %d > %d for node %d",numdofs,(int)dofs.size(),newnode->Id());
        (*idxcolnodes_)[newlid] = dofs[0];
        for (int idof = 0; idof < numdofs; ++idof)
        {
          dofcolvec.push_back(dofs[idof]);
        }
      }
    }

    dofcolmap_ = rcp(new Epetra_Map(-1, dofcolvec.size(), &dofcolvec[0], 0, newdis.Comm()));

}


/// Assign dof numbers for new discretization using dof numbering from source discretization.
void DRT::TransparentDofSet::ParallelTransferDegreesOfFreedom(
  const vector<RCP<DRT::Discretization> > sourcedis,
  const DRT::Discretization& newdis,
  const int start
  )
{
  for (int i=0; i<int(sourcedis.size()); i++)
  {
    if (!sourcedis[i]->DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
    if (!sourcedis[i]->NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
    if (!sourcedis[i]->ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");
  }

  if (!newdis.DofRowMap()->UniqueGIDs()) dserror("DofRowMap is not unique");
  if (!newdis.NodeRowMap()->UniqueGIDs()) dserror("NodeRowMap is not unique");
  if (!newdis.ElementRowMap()->UniqueGIDs()) dserror("ElementRowMap is not unique");

  // list all my rownode ids

  //
  // we need a mapping
  //
  // colnode gid -> vector<int> dofs of sourcenode
  //
  // problem: sourcenode not necessarily on this proc -> communicate
  //
  // the idea is to search for the sourcerownode on some proc and to get
  // this unique number
  //
  map<int,vector<int> >  gid_to_dofs;
  for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
  {
    const DRT::Node* newnode = newdis.lColNode(inode);
    int gid=newnode->Id();
    vector<int> emptyvec;
    gid_to_dofs.insert(pair<int, vector<int> >(gid,emptyvec));
  }

  {
#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(sourcedis[0]->Comm());

    // necessary variables
    MPI_Request request;
#endif

    // define send and receive blocks
    vector<char> sblock;
    vector<char> rblock;

    // get number of processors and the current processors id
    int numproc=sourcedis[0]->Comm().NumProc();
    int myrank=sourcedis[0]->Comm().MyPID();

    //----------------------------------------------------------------------
    // communication is done in a round robin loop
    //----------------------------------------------------------------------
    for (int np=0;np<numproc+1;++np)
    {
      // in the first step, we cannot receive anything
      if(np >0)
      {
#ifdef PARALLEL
        ReceiveBlock(numproc,myrank,rblock,exporter,request);
#else
        rblock=sblock;
#endif

        // Unpack info from the receive block from the last proc
        UnpackLocalSourceDofs(gid_to_dofs,rblock);
      }

      // in the last step, we keep everything on this proc
      if(np < numproc)
      {
        // -----------------------
        // do what we wanted to do
        SetSourceDofsAvailableOnThisProc(gid_to_dofs);

        // Pack info into block to send
        DRT::PackBuffer data;
        PackLocalSourceDofs(gid_to_dofs,data);
        data.StartPacking();
        PackLocalSourceDofs(gid_to_dofs,data);
        gid_to_dofs.clear();
        swap( sblock, data() );

#ifdef PARALLEL
        SendBlock(numproc,myrank,sblock,exporter,request);
#endif
      }
    }
  }

  set<int> slaveset;

  std::vector<DRT::Condition*> mypbcs;

  // get periodic surface boundary conditions
  sourcedis_[0]->GetCondition("SurfacePeriodic",mypbcs);

  if(mypbcs.empty())
  {
    sourcedis_[0]->GetCondition("LinePeriodic",mypbcs);
  }

  for (unsigned numcond=0;numcond<mypbcs.size();++numcond)
  {
    DRT::Condition* thiscond= mypbcs[numcond];

    // see whether we have a slave condition
    const string* mymasterslavetoggle
      = thiscond->Get<string>("Is slave periodic boundary condition");

    // test whether in all source discretizations the same periodic bc are slave or master
    if(int(sourcedis.size()) != 1)
      cout<<"Warning! ParallelTransferDegreesOfFreedom hasn't been tested for this application!!"<<endl;
    for (int i=1; i<int(sourcedis.size()); i++)
    {
      std::vector<DRT::Condition*> testpbcs;
      sourcedis_[i]->GetCondition("SurfacePeriodic",testpbcs);

      if(testpbcs.empty())
      {
        sourcedis_[i]->GetCondition("LinePeriodic",testpbcs);
      }

      const string* testmasterslavetoggle
        = testpbcs[numcond]->Get<string>("Is slave periodic boundary condition");

      if( *testmasterslavetoggle != *mymasterslavetoggle)
        dserror("the same periodic boundary conditions are different for different source discretizations");
    } // test end

    // normal prcedure with periodic bcs
    if(!(*mymasterslavetoggle=="Master"))
    {
      const vector <int>* pbcids;
      pbcids = (*thiscond).Nodes();

      for(vector<int>::const_iterator iter=pbcids->begin();iter!=pbcids->end();++iter)
      {
        slaveset.insert(*iter);
      }
    }
  }


  // build dofrowmap based on dofrowvec which contains the dofs of each node from the (multiple) source dis;
  // dofs of the first node are followed by dofs of the second node and so on
  vector<int> dofrowvec;
  dofrowvec.reserve(dofrowmap_->NumMyElements());
  for (int inode = 0; inode != newdis.NumMyRowNodes(); ++inode)
  {
    const DRT::Node* newnode = newdis.lRowNode(inode);

    const vector<int> dofs = gid_to_dofs[newnode->Id()];

    const int newlid = newnode->LID();
    const int numdofs = (*numdfcolnodes_)[newlid];

    if(numdofs!=(int)dofs.size())
    {
      printf("This is node %d  (%12.5e,%12.5e,%12.5e)\n",
             newnode->Id(),
             newnode->X()[0],
             newnode->X()[1],
             newnode->X()[2]);

      dserror("spooky, isn't it? dofs to overwrite %d != %d dofs.size() to set \n",numdofs,dofs.size());
    }

    if (numdofs > 0)
    {
      (*idxcolnodes_)[newlid] = dofs[0];

      // slave-dofs must not enter the dofrowvec (if master&slave are on different procs)
      std::set<int>::iterator curr=slaveset.find(newnode->Id());

      if(curr==slaveset.end())
      {
        for (int idof = 0; idof < numdofs; ++idof)
        {
          dofrowvec.push_back(dofs[idof]);
        }
      }
    }
  }

  dofrowmap_ = rcp(new Epetra_Map(-1, dofrowvec.size(), &dofrowvec[0], 0, newdis.Comm()));

  // build dofcolmap based on dofcolvec which contains the dofs of each node from the (multiple) source dis;
  // dofs of the first node are followed by dofs of the second node and so on
  vector<int> dofcolvec;
  dofcolvec.reserve(dofcolmap_->NumMyElements());
  for (int inode = 0; inode != newdis.NumMyColNodes(); ++inode)
  {
    const DRT::Node* newnode = newdis.lColNode(inode);

    const vector<int> dofs = gid_to_dofs[newnode->Id()];

    const int newlid = newnode->LID();
    //const int newfirstidx = (*idxcolnodes_)[newlid];
    const int numdofs = (*numdfcolnodes_)[newlid];
    if (numdofs > 0)
    {
      (*idxcolnodes_)[newlid] = dofs[0];
      if(numdofs!=(int)dofs.size())
      dserror("numdofs %d!=%d for node %d",numdofs,dofs.size(),newnode->Id());

      for (int idof = 0; idof < numdofs; ++idof)
      {
        dofcolvec.push_back(dofs[idof]);
      }
    }
  }

  dofcolmap_ = rcp(new Epetra_Map(-1, dofcolvec.size(), &dofcolvec[0], 0, newdis.Comm()));

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | write dof info into map                                    (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::TransparentDofSet::SetSourceDofsAvailableOnThisProc(
  map<int,vector<int> >   & gid_to_dofs
  )
{
  for(map<int,vector<int> >::iterator curr=gid_to_dofs.begin();
    curr!=gid_to_dofs.end();++curr)
  {

    // test whether all source nodes in all discrets are on this proc
    int lid = 0;
    for (int i=0; i<int(sourcedis_.size()); i++)
    {
      if( sourcedis_[i]->NodeRowMap()->LID(curr->first) == -1)
        lid = sourcedis_[i]->NodeRowMap()->LID(curr->first);
    }

    if (lid>-1)
    {
      curr->second.clear();
      // Gather dofs of all source discretizations
      for (int i=0; i<int(sourcedis_.size()); i++)
      {
        const DRT::Node* sourcenode = sourcedis_[i]->gNode(curr->first);

        const vector<int> dofs = sourcedis_[i]->Dof(sourcenode);

        for (vector<int>::const_iterator iter=dofs.begin();iter!=dofs.end();++iter)
        {
          curr->second.push_back(*iter);
        }
      }
    }
    else
    {
      int numproc=sourcedis_[0]->Comm().NumProc();
      if(numproc==1)
      {
        dserror("I have a one-processor problem but the node is not on the proc. sourcedis_->NodeRowMap() is probably currupted.");
      }
    }
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | pack all values into a send block                          (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::TransparentDofSet::PackLocalSourceDofs(
  map<int,vector<int> >   & gid_to_dofs  ,
  DRT::PackBuffer         & sblock
  )
{
  int size=gid_to_dofs.size();

  // add size  to sendblock
  DRT::ParObject::AddtoPack(sblock,size);

  for(map<int,vector<int> >::iterator curr=gid_to_dofs.begin();
    curr!=gid_to_dofs.end();++curr)
  {
    int         gid    = curr->first;
    vector<int> mydofs = curr->second;
    int numdofs = (int)mydofs.size();

    DRT::ParObject::AddtoPack(sblock,gid);
    DRT::ParObject::AddtoPack(sblock,numdofs);
    for(int ll=0;ll<numdofs;++ll)
    {
      DRT::ParObject::AddtoPack(sblock,mydofs[ll]);
    }
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | unpack all values contained in receive block               (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::TransparentDofSet::UnpackLocalSourceDofs(
  map<int,vector<int> >   & gid_to_dofs  ,
  vector<char>            & rblock
  )
{
  gid_to_dofs.clear();

  // position to extract
  vector<char>::size_type position = 0;

  // extract size
  int size=0;
  DRT::ParObject::ExtractfromPack(position,rblock,size);

  for(int rr=0;rr<size;++rr)
  {
    int         gid    = -1;
    vector<int> mydofs ;
    int numdofs = 0;

    DRT::ParObject::ExtractfromPack(position,rblock,gid);
    DRT::ParObject::ExtractfromPack(position,rblock,numdofs);
    
    for(int ll=0;ll<numdofs;++ll)
    {
      int thisdof=0;

      DRT::ParObject::ExtractfromPack(position,rblock,thisdof);
      mydofs.push_back(thisdof);
    }

    gid_to_dofs.insert(pair<int, vector<int> >(gid,mydofs));
  }

  rblock.clear();
  return;
}

#ifdef PARALLEL
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | receive a block in the round robin communication pattern   (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::TransparentDofSet::ReceiveBlock(
  int             numproc ,
  int             myrank  ,
  vector<char>   & rblock,
  DRT::Exporter  & exporter,
  MPI_Request    & request)
{

  // necessary variables

  int         length =-1;
  int         frompid=(myrank+numproc-1)%numproc;
  int         tag    =frompid;

  // make sure that you do not think you received something if
  // you didn't
  if(rblock.empty()==false)
  {
    dserror("rblock not empty");
  }

  // receive from predecessor
  exporter.ReceiveAny(frompid,tag,rblock,length);

  if(tag!=(myrank+numproc-1)%numproc)
  {
    dserror("received wrong message (ReceiveAny)");
  }

  exporter.Wait(request);

  // for safety
  exporter.Comm().Barrier();

  return;
} // ReceiveBlock
#endif



#ifdef PARALLEL
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | send a block in the round robin communication pattern      (private) |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void DRT::TransparentDofSet::SendBlock(
  int             numproc ,
  int             myrank  ,
  vector<char>  & sblock  ,
  DRT::Exporter & exporter,
  MPI_Request   & request )
{

  // Send block to next proc.
  int         tag    =myrank;
  int         frompid=myrank;
  int         topid  =(myrank+1)%numproc;

  exporter.ISend(frompid,topid,
                 &(sblock[0]),sblock.size(),
                 tag,request);


  // for safety
  exporter.Comm().Barrier();

  return;
} //SendBlock
#endif
#endif  // #ifdef CCADISCRET
