/*!----------------------------------------------------------------------
\file drt_transfer_turb_inflow.cpp

\brief Methods to transfer a turbulent inflow profile from a (usually
periodic boundary condition) separate domain to a name Dirichlet
boundary of the actual domain

\level 3

<pre>
\maintainer Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#include "drt_transfer_turb_inflow.H"
#include "../drt_lib/drt_matchingoctree.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 03/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::TransferTurbulentInflowCondition::TransferTurbulentInflowCondition(
  Teuchos::RCP<DRT::Discretization>  dis    ,
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps
    )
  : dis_(dis),
    dbcmaps_(dbcmaps),
    curve_(-1),
    numveldof_(3)
{
  active_=false;

  // vector of pointers to all node clouds i.e. conditions to couple
  std::vector<DRT::Condition*>     nodecloudstocouple;

  // get surfaces to couple
  dis_->GetCondition("TransferTurbulentInflow",nodecloudstocouple);

  if(not nodecloudstocouple.empty())
  {
    // activate
    active_=true;

    // master and slave sets to couple
    std::set<int> masterset;
    std::set<int> slaveset;

    // global master node Ids and global slave node Ids
    std::vector <int> masternodeids;
    std::vector <int> slavenodeids;

    // the (at the moment) one and only direction to couple
    int dir=-1;

    // loop all conditions and check whether they are of master or slave
    // type
    for(std::vector<DRT::Condition*>::iterator cond=nodecloudstocouple.begin();
        cond!=nodecloudstocouple.end();
        ++cond)
    {
      // get id, direction info and toggle
      int        id       =-1;
      int        direction=-1;
      ToggleType toggle   =none;

      GetData(id,direction,toggle,*cond);

      if(dir==-1)
      {
        dir=direction;
      }
      else
      {
        if(dir!=direction)
        {
          dserror("multiple directions are not supported yet");
        }
      }

      if(id!=0)
      {
        dserror("expecting only one group of coupling surfaces (up to now), its %d",id);
      }

      switch (toggle)
      {
      case master:
      {
        //--------------------------------------------------
        // get global master node Ids
        const std::vector <int>* masteridstoadd;

        masteridstoadd = (*cond)->Nodes();

        for(std::vector<int>::const_iterator idtoadd =(*masteridstoadd).begin();
            idtoadd!=(*masteridstoadd).end();
            ++idtoadd)
        {
          // we construct the local octree only with nodes owned by this proc
          if(dis_->HaveGlobalNode(*idtoadd))
            if(dis_->gNode(*idtoadd)->Owner() == dis_->Comm().MyPID())
              masterset.insert(*idtoadd);
        }

        break;
      }
      case slave:
      {
        //--------------------------------------------------
        // get global slave node Ids
        const std::vector <int>* slaveidstoadd;

        slaveidstoadd = (*cond)->Nodes();

        for(std::vector<int>::const_iterator idtoadd =(*slaveidstoadd).begin();
            idtoadd!=(*slaveidstoadd).end();
            ++idtoadd)
        {
          // we only try to match owned nodes of each proc
          if(dis_->HaveGlobalNode(*idtoadd))
            if(dis_->gNode(*idtoadd)->Owner() == dis_->Comm().MyPID())
              slaveset.insert(*idtoadd);
        }

        break;
      }
      default :
        dserror("toggle non master or slave");
      }

    }

    //--------------------------------------------------
    // just write sets into vectors
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

    // these are just parameter definitions for the octree search algorithm
    const double tol            = 1e-6;
    const int    maxnodeperleaf = 250;
    const double rotangle       = 0.0;

    std::vector<int> dofsforpbcplane(2);
    {
      int mm=0;
      for(int rr=0;rr<3;++rr)
      {
        if(rr!=dir)
        {
          dofsforpbcplane[mm]=rr;
          ++mm;
        }
      }
    }

    // build processor local octree
    DRT::UTILS::NodeMatchingOctree nodematchingoctree = DRT::UTILS::NodeMatchingOctree();
    nodematchingoctree.Init(*dis_,masternodeids,maxnodeperleaf,tol);
    nodematchingoctree.Setup();

    // create map from gid masternode -> gid corresponding slavenode
    nodematchingoctree.CreateGlobalEntityMatching(
        slavenodeids   ,
        dofsforpbcplane,
        rotangle,
        midtosid_
    );
    // sanity check
    for (std::map<int,std::vector<int> >::iterator pair=midtosid_.begin();
         pair!=midtosid_.end();
         ++pair)
    {
      if(pair->second.size()!=1)
      {
        dserror("expected one node to match, got %d out of %d",
        pair->second.size(),
        slavenodeids.size());
      }
    }
  }
  return;
} // TransferTurbulentInflowCondition::TransferTurbulentInflowCondition

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Perform transfer process (public)                        gammi 03/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::Transfer(
  const Teuchos::RCP<Epetra_Vector> veln ,
  Teuchos::RCP<Epetra_Vector>       velnp,
  const double                      time)
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();

  std::vector<int>             mymasters;
  std::vector<std::vector<double> > mymasters_vel(numveldof_);

  if(active_)
  {

      // initialization of time curve factor
      double curvefac    = 1.0;
      if (curve_ >= 0) // yes, we have a time curve
      {
        // time factor for the intermediate step
        if(time >= 0.0)
        {
          curvefac = DRT::Problem::Instance()->Funct(curve_).EvaluateTime(time);
        }
        else
        {
          // do not compute an "alternative" curvefac here since a negative time value
          // indicates an error.
          dserror("Negative time value: time = %f",time);
        }
      }
      else // we do not have a time curve --- time curve factor is constant equal 1
      {
        // nothing to do
      }

    // collect masters on this proc and associated velocities
    for (std::map<int,std::vector<int> >::iterator pair=midtosid_.begin();
         pair!=midtosid_.end();
         ++pair)
    {
      int gid=pair->first;

      if(dis_->HaveGlobalNode(gid))
      {
        mymasters.push_back(gid);

        DRT::Node* master=dis_->gNode(gid);

        std::vector<int> masterdofs=dis_->Dof(master);

        for(int rr=0;rr<3;++rr)
        {
          int lid = dofrowmap->LID(masterdofs[rr]);

          (mymasters_vel[rr]).push_back(((*veln)[lid])*curvefac);
        }
      }
      else
      {
  dserror("master %d in midtosid but not on proc. This was unexpected",gid);
      }
    }

#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(dis_->Comm());

    // necessary variables
    MPI_Request request;
#endif

    // define send and receive blocks
    std::vector<char> sblock;
    std::vector<char> rblock;

    // get number of processors and the current processors id
    int numproc=dis_->Comm().NumProc();

    //----------------------------------------------------------------------
    // communication is done in a round robin loop
    //----------------------------------------------------------------------
    for (int np=0;np<numproc+1;++np)
    {
      // in the first step, we cannot receive anything
      if(np >0)
      {
#ifdef PARALLEL
        ReceiveBlock(rblock,exporter,request);
#else
        rblock=sblock;
#endif

        // Unpack info from the receive block from the last proc
        UnpackLocalMasterValues(mymasters,mymasters_vel,rblock);
      }

      // in the last step, we keep everything on this proc
      if(np < numproc)
      {
        // -----------------------
        // do what we wanted to do
        SetValuesAvailableOnThisProc(mymasters,mymasters_vel,velnp);

        // Pack info into block to send
        DRT::PackBuffer data;
        PackLocalMasterValues(mymasters,mymasters_vel,data);
        data.StartPacking();
        PackLocalMasterValues(mymasters,mymasters_vel,data);
        swap( sblock, data() );

#ifdef PARALLEL
        SendBlock(sblock,exporter,request);
#endif
      }
    }
  }
  return;
} // Transfer

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor dtor  (public)                                 gammi 03/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::TransferTurbulentInflowCondition::~TransferTurbulentInflowCondition()
{
  return;
}// ~TransferTurbulentInflowCondition

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | get condition id etc (private)                            gammi 03/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::GetData(
  int                 & id       ,
  int                 & direction,
  ToggleType          & type     ,
  const DRT::Condition* cond     )
{

  const std::vector<int>* myid = cond->Get<std::vector<int> >("id");
  id=(*myid)[0];

  const std::string* mydirection = cond->Get<std::string>("transfer direction");
  if      (*mydirection == "x")
  {
    direction=0;
  }
  else if (*mydirection == "y")
  {
    direction=1;
  }
  else if (*mydirection == "z")
  {
    direction=2;
  }
  else
  {
    dserror("unknown direction");
  }

  const std::string* mytoggle = cond->Get<std::string>("toggle");
  if      (*mytoggle == "master")
  {
    type=master;
  }
  else if (*mytoggle == "slave")
  {
    type=slave;
  }
  else
  {
    dserror("expecting either master or slave");
  }

  // find out whether we will use a time curve
  if (curve_ == -1)
  {
    const std::vector<int>* curve = cond->Get<std::vector<int> >("curve");

    // set curve number
    if (curve) curve_ = (*curve)[0];
  }

  return;
}


#ifdef PARALLEL
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | receive a block in the round robin communication pattern   (private) |
 |                                                          gammi 03/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::ReceiveBlock(
    std::vector<char>   & rblock,
    DRT::Exporter  & exporter,
    MPI_Request    & request)
{
  // get number of processors and the current processors id
  int numproc=dis_->Comm().NumProc();
  int myrank =dis_->Comm().MyPID();

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
} // TransferTurbulentInflowCondition::ReceiveBlock
#endif



#ifdef PARALLEL
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | send a block in the round robin communication pattern      (private) |
 |                                                          gammi 03/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::SendBlock(
    std::vector<char>  & sblock  ,
    DRT::Exporter & exporter,
    MPI_Request   & request )
{
  // get number of processors and the current processors id
  int numproc=dis_->Comm().NumProc();
  int myrank =dis_->Comm().MyPID();

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
} // TransferTurbulentInflowCondition::SendBlock
#endif



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | unpack all master values contained in receive block        (private) |
 |                                                          gammi 03/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::UnpackLocalMasterValues(
    std::vector<int>             & mymasters    ,
    std::vector<std::vector<double> > & mymasters_vel,
    std::vector<char>            & rblock
  )
{
  mymasters.clear();

  midtosid_.clear();

  if((int)mymasters_vel.size()!=numveldof_)
  {
    dserror("expecting three spatial dimensions in mymasters_vel to unpack into");
  }

  for(int rr=0;rr<numveldof_;++rr)
  {
    (mymasters_vel[rr]).clear();
  }

  // position to extract
  std::vector<char>::size_type position = 0;

  // extract size
  int size=0;
  DRT::ParObject::ExtractfromPack(position,rblock,size);

  // extract master ids
  for(int i=0;i<size;++i)
  {
    int id;

    DRT::ParObject::ExtractfromPack(position,rblock,id);
    mymasters.push_back(id);

    std::map<int,std::vector<int> >::iterator iter=midtosid_.find(id);

    if(iter!=midtosid_.end())
    {
      iter->second.clear();
    }
    else
    {
      midtosid_[id].clear();
    }
  }

  // extract slave ids
  for(int rr=0;rr<size;++rr)
  {
    int slavesize;

    DRT::ParObject::ExtractfromPack(position,rblock,slavesize);

    for(int ll=0;ll<slavesize;++ll)
    {
      int sid;
      DRT::ParObject::ExtractfromPack(position,rblock,sid);

      std::map<int,std::vector<int> >::iterator iter=midtosid_.find(mymasters[rr]);

      if(iter!=midtosid_.end())
      {
  iter->second.push_back(sid);
      }
      else
      {
  dserror("master id %d was not in midtosid_",mymasters[rr]);
      }
    }

    if(midtosid_[mymasters[rr]].size()<1)
    {
      dserror("require at least one slave to master %d, got %d",mymasters[rr],midtosid_[mymasters[rr]].size());
    }
  }

  // extract values (first u, then v, then w)
  for(int mm=0;mm<numveldof_;++mm)
  {
    for(int rr=0;rr<size;++rr)
    {
      double value;

      DRT::ParObject::ExtractfromPack(position,rblock,value);

      (mymasters_vel[mm]).push_back(value);
    }
  }

  rblock.clear();
  return;
} // UnpackLocalMasterValues



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | pack all master values into a send block                   (private) |
 |                                                          gammi 03/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::PackLocalMasterValues(
    std::vector<int>             & mymasters    ,
    std::vector<std::vector<double> > & mymasters_vel,
    DRT::PackBuffer         & sblock
    )
{
  int size=mymasters.size();

  if(mymasters_vel.size()!=(unsigned)numveldof_)
  {
    dserror("expecting three spatial dimensions in mymasters_vel to pack");
  }

  for(int rr=0;rr<numveldof_;++rr)
  {
    if((int)(mymasters_vel[rr]).size()!=size)
    {
      dserror("at least one of the components of mymasters_vel has the wrong size");
    }
  }

  // add size  to sendblock
  DRT::ParObject::AddtoPack(sblock,size);

  // add master ids
  for(int rr=0;rr<size;++rr)
  {
    DRT::ParObject::AddtoPack(sblock,mymasters[rr]);
  }

  // add slave ids
  for(int rr=0;rr<size;++rr)
  {
    std::map<int,std::vector<int> >::iterator iter=midtosid_.find(mymasters[rr]);

    if(iter==midtosid_.end())
    {
      dserror("tried to pack slaves to master master %d, got none",mymasters[rr]);
    }
    else
    {
      std::vector<int> slaves=iter->second;

      int slavesize=(int)slaves.size();

      DRT::ParObject::AddtoPack(sblock,slavesize);
      for(int ll=0;ll<slavesize;++ll)
      {
  DRT::ParObject::AddtoPack(sblock,slaves[ll]);
      }
    }
  }

  // add values (first u, then v, then w)
  for(int mm=0;mm<numveldof_;++mm)
  {
    for(int rr=0;rr<size;++rr)
    {
      DRT::ParObject::AddtoPack(sblock,(mymasters_vel[mm])[rr]);
    }
  }

  return;
} // PackLocalMasterValues

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | for all values avaible on the processor, do the final setting of the |
 | value                          (private)                 gammi 03/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowCondition::SetValuesAvailableOnThisProc(
    std::vector<int>                 & mymasters,
    std::vector<std::vector<double> >     & mymasters_vel,
    Teuchos::RCP<Epetra_Vector>   velnp)
{
  const Teuchos::RCP<const Epetra_Map> activedbcdofs=dbcmaps_->CondMap();

  for(unsigned nn=0;nn<mymasters.size();++nn)
  {
    std::map<int,std::vector<int> >::iterator iter=midtosid_.find(mymasters[nn]);

    if(iter!=midtosid_.end())
    {
      std::vector<int> myslaves(iter->second);

      for(std::vector<int>::iterator sid=myslaves.begin();sid!=myslaves.end();++sid)
      {
        // is this slave id on this proc?
        if(dis_->NodeRowMap()->MyGID(*sid))
        {
          DRT::Node* slave=dis_->gNode(*sid);

          // get dofs
          std::vector<int> slavedofs=dis_->Dof(slave);

          for(int rr=0;rr<numveldof_;++rr)
          {
            int gid = slavedofs[rr];

            // only set if DBC is active, otherwise throw error
            if(activedbcdofs->MyGID(gid))
            {
              double value=(mymasters_vel[rr])[nn];

              velnp->ReplaceGlobalValues(1,&value,&gid);
            }
            else
            {
              int    id=slave->Id();

              double x=slave->X()[0];
              double y=slave->X()[1];
              double z=slave->X()[2];

              dserror("Dirichlet condition required on slave node (%12.5e,%12.5e,%12.5e), id %d, dof %d of transfer condition",x,y,z,id,rr);
            }
          }
        }
      }
    }
  }

  return;
} // SetValuesAvailableOnThisProc


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                        bk 09/14|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::TransferTurbulentInflowConditionXW::TransferTurbulentInflowConditionXW(
  Teuchos::RCP<DRT::Discretization>  dis    ,
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps
    )
  : TransferTurbulentInflowCondition(dis,dbcmaps)
{
  numveldof_=6;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor dtor  (public)                                    bk 09/14|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::TransferTurbulentInflowConditionXW::~TransferTurbulentInflowConditionXW()
{
  return;
}// ~TransferTurbulentInflowCondition (XW)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Perform transfer process (public)                           bk 09/14|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowConditionXW::Transfer(
  const Teuchos::RCP<Epetra_Vector> veln ,
  Teuchos::RCP<Epetra_Vector>       velnp,
  const double                      time)
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();

  std::vector<int>             mymasters;
  //there can be up to 6 velocity dofs per node (3 +3 virtual dofs)
  std::vector<std::vector<double> > mymasters_vel(numveldof_);

  if(active_)
  {

      // initialization of time curve factor
      double curvefac    = 1.0;
      if (curve_ >= 0) // yes, we have a time curve
      {
        // time factor for the intermediate step
        if(time >= 0.0)
        {
          curvefac = DRT::Problem::Instance()->Funct(curve_).EvaluateTime(time);
        }
        else
        {
          // do not compute an "alternative" curvefac here since a negative time value
          // indicates an error.
          dserror("Negative time value: time = %f",time);
        }
      }
      else // we do not have a time curve --- time curve factor is constant equal 1
      {
        // nothing to do
      }

    // collect masters on this proc and associated velocities
    for (std::map<int,std::vector<int> >::iterator pair=midtosid_.begin();
         pair!=midtosid_.end();
         ++pair)
    {
      int gid=pair->first;

      if(dis_->HaveGlobalNode(gid))
      {
        mymasters.push_back(gid);

        DRT::Node* master=dis_->gNode(gid);

        std::vector<int> masterdofs=dis_->Dof(master);

        for(int rr=0;rr<3;++rr)
        {
          int lid = dofrowmap->LID(masterdofs[rr]);

          (mymasters_vel[rr]).push_back(((*veln)[lid])*curvefac);
        }

        //in xwall, we have another virtual node right after this node
        if(dis_->NumDof(master)==8)
        {
          for(int rr=4;rr<7;++rr)
          {
            int lid = dofrowmap->LID(masterdofs[rr]);

            (mymasters_vel[rr-1]).push_back(((*veln)[lid])*curvefac);
          }
        }
        else
          for(int rr=4;rr<7;++rr)
          {
            (mymasters_vel[rr-1]).push_back(0.0);
          }
      }
      else
      {
  dserror("master %d in midtosid but not on proc. This was unexpected",gid);
      }
    }

#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(dis_->Comm());

    // necessary variables
    MPI_Request request;
#endif

    // define send and receive blocks
    std::vector<char> sblock;
    std::vector<char> rblock;

    // get number of processors and the current processors id
    int numproc=dis_->Comm().NumProc();

    //----------------------------------------------------------------------
    // communication is done in a round robin loop
    //----------------------------------------------------------------------
    for (int np=0;np<numproc+1;++np)
    {
      // in the first step, we cannot receive anything
      if(np >0)
      {
#ifdef PARALLEL
        ReceiveBlock(rblock,exporter,request);
#else
        rblock=sblock;
#endif

        // Unpack info from the receive block from the last proc
        UnpackLocalMasterValues(mymasters,mymasters_vel,rblock);
      }

      // in the last step, we keep everything on this proc
      if(np < numproc)
      {
        // -----------------------
        // do what we wanted to do
        SetValuesAvailableOnThisProc(mymasters,mymasters_vel,velnp);

        // Pack info into block to send
        DRT::PackBuffer data;
        PackLocalMasterValues(mymasters,mymasters_vel,data);
        data.StartPacking();
        PackLocalMasterValues(mymasters,mymasters_vel,data);
        swap( sblock, data() );

#ifdef PARALLEL
        SendBlock(sblock,exporter,request);
#endif
      }
    }
  }
  return;
} // Transfer (XW)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | for all values avaible on the processor, do the final setting of the |
 | value                          (private)                    bk 09/14 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowConditionXW::SetValuesAvailableOnThisProc(
    std::vector<int>                 & mymasters,
    std::vector<std::vector<double> >     & mymasters_vel,
    Teuchos::RCP<Epetra_Vector>   velnp)
{
  const Teuchos::RCP<const Epetra_Map> activedbcdofs=dbcmaps_->CondMap();

  for(unsigned nn=0;nn<mymasters.size();++nn)
  {
    std::map<int,std::vector<int> >::iterator iter=midtosid_.find(mymasters[nn]);

    if(iter!=midtosid_.end())
    {
      std::vector<int> myslaves(iter->second);

      for(std::vector<int>::iterator sid=myslaves.begin();sid!=myslaves.end();++sid)
      {
        // is this slave id on this proc?
        if(dis_->NodeRowMap()->MyGID(*sid))
        {
          DRT::Node* slave=dis_->gNode(*sid);

          // get dofs
          std::vector<int> slavedofs=dis_->Dof(slave);

          for(int rr=0;rr<3;++rr)
          {
            int gid = slavedofs[rr];

            // only set if DBC is active, otherwise throw error
            if(activedbcdofs->MyGID(gid))
            {
              double value=(mymasters_vel[rr])[nn];

              velnp->ReplaceGlobalValues(1,&value,&gid);
            }
            else
            {
              int    id=slave->Id();

              double x=slave->X()[0];
              double y=slave->X()[1];
              double z=slave->X()[2];

              dserror("Dirichlet condition required on slave node (%12.5e,%12.5e,%12.5e), id %d, dof %d of transfer condition",x,y,z,id,rr);
            }
          }

          //and treat xwall dofs
          if(dis_->NumDof(slave)==8)
            for(int rr=4;rr<7;++rr)
            {
              int gid = slavedofs[rr];

              // only set if DBC is active, otherwise throw error
              if(activedbcdofs->MyGID(gid))
              {
                double value=(mymasters_vel[rr-1])[nn];

                velnp->ReplaceGlobalValues(1,&value,&gid);
              }
              else
                dserror("xwall dofs don't have active dbc for transfer");
            }
        }
      }
    }
  }

  return;
} // SetValuesAvailableOnThisProc (XW)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                        bk 09/14|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::TransferTurbulentInflowConditionNodal::TransferTurbulentInflowConditionNodal(
  Teuchos::RCP<DRT::Discretization>  dis    ,
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps
    )
  : TransferTurbulentInflowCondition(dis,dbcmaps)
{
  numveldof_=1;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor dtor  (public)                                    bk 09/14|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::TransferTurbulentInflowConditionNodal::~TransferTurbulentInflowConditionNodal()
{
  return;
}// ~TransferTurbulentInflowCondition (Nodal)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Perform transfer process (public)                           bk 09/14|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowConditionNodal::Transfer(
  const Teuchos::RCP<Epetra_Vector> invec ,
  Teuchos::RCP<Epetra_Vector>       outvec,
  const double                      time)
{


  std::vector<int>             mymasters;
  std::vector<std::vector<double> > mymasters_vec(numveldof_);

  if(active_)
  {

    // collect masters on this proc and associated velocities
    for (std::map<int,std::vector<int> >::iterator pair=midtosid_.begin();
         pair!=midtosid_.end();
         ++pair)
    {
      int gid=pair->first;

      if(dis_->HaveGlobalNode(gid))
      {
        mymasters.push_back(gid);

        DRT::Node* master=dis_->gNode(gid);

        std::vector<int> masterdofs=dis_->Dof(master);

        // and the 7th value is filled with the wall shear stress
        int lnodeid=dis_->NodeRowMap()->LID(gid);
        (mymasters_vec[0]).push_back((*invec)[lnodeid]);
      }
      else
      {
  dserror("master %d in midtosid but not on proc. This was unexpected",gid);
      }
    }

#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(dis_->Comm());

    // necessary variables
    MPI_Request request;
#endif

    // define send and receive blocks
    std::vector<char> sblock;
    std::vector<char> rblock;

    // get number of processors and the current processors id
    int numproc=dis_->Comm().NumProc();

    //----------------------------------------------------------------------
    // communication is done in a round robin loop
    //----------------------------------------------------------------------
    for (int np=0;np<numproc+1;++np)
    {
      // in the first step, we cannot receive anything
      if(np >0)
      {
#ifdef PARALLEL
        ReceiveBlock(rblock,exporter,request);
#else
        rblock=sblock;
#endif

        // Unpack info from the receive block from the last proc
        UnpackLocalMasterValues(mymasters,mymasters_vec,rblock);
      }

      // in the last step, we keep everything on this proc
      if(np < numproc)
      {
        // -----------------------
        // do what we wanted to do
        SetValuesAvailableOnThisProc(mymasters,mymasters_vec,outvec);

        // Pack info into block to send
        DRT::PackBuffer data;
        PackLocalMasterValues(mymasters,mymasters_vec,data);
        data.StartPacking();
        PackLocalMasterValues(mymasters,mymasters_vec,data);
        swap( sblock, data() );

#ifdef PARALLEL
        SendBlock(sblock,exporter,request);
#endif
      }
    }
  }
  return;
} // Transfer (Nodal)

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | for all values avaible on the processor, do the final setting of the |
 | value                          (private)                    bk 09/14 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::TransferTurbulentInflowConditionNodal::SetValuesAvailableOnThisProc(
    std::vector<int>                 & mymasters,
    std::vector<std::vector<double> >     & mymasters_vec,
    Teuchos::RCP<Epetra_Vector>   outvec)
{

  for(unsigned nn=0;nn<mymasters.size();++nn)
  {
    std::map<int,std::vector<int> >::iterator iter=midtosid_.find(mymasters[nn]);

    if(iter!=midtosid_.end())
    {
      std::vector<int> myslaves(iter->second);

      for(std::vector<int>::iterator sid=myslaves.begin();sid!=myslaves.end();++sid)
      {
        // is this slave id on this proc?
        if(dis_->NodeRowMap()->MyGID(*sid))
          outvec->ReplaceGlobalValue(*sid,0,(mymasters_vec[0])[nn]);
      }
    }
  }

  return;
} // SetValuesAvailableOnThisProc (Nodal)


