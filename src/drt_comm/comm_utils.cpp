/*----------------------------------------------------------------------*/
/*!
\file comm_utils.cpp

\brief Helper class for everything that deals with communication

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include <Epetra_MpiComm.h>
#include <mpi.h>

#include <vector>
#include <sstream>
#include <string>

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#include "comm_utils.H"


/*----------------------------------------------------------------------*
 | create communicator                                      ghamm 02/12 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::CreateComm(int argc, char** argv)
{
  // for coupled simulations: color = 1 for BACI and color = 0 for other programs
  // so far: either nested parallelism within BACI or coupling with further
  // executables is possible
  // default values without nested parallelism
  int myrank = -1;
  int size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int color = 1;
  int ngroup = 1;
  NP_TYPE npType = no_nested_parallelism;

  // parse command line and separate configuration arguments
  std::vector<char*> conf(0);
  for(int i=1; i<argc; i++)
  {
    std::string temp=argv[i];
    if (temp.substr( 0, 1 ) == "-")
    {
      conf.push_back(argv[i]);
    }
  }

  // grouplayout will be filled accordingly to the user given input
  std::vector<int> grouplayout;
  bool ngroupisset = false;
  bool nptypeisset = false;
  for(int i=0; i<int(conf.size()); i++)
  {
    //----------------------------------------------------------------
    // determine number of groups and the proc distribution
    //----------------------------------------------------------------
    std::string numgroup(conf[i]);
    if (numgroup.substr( 0, 8 ) == "-ngroup=")
    {
      ngroupisset = true;
      ngroup = atoi( numgroup.substr( 8, std::string::npos ).c_str() );

      // read out argument after ngroup=
      std::string glayout;
      if(i+1 < int(conf.size()))
        glayout=conf[i+1];
      else
        glayout="dummy";

      // case with given group layout
      if (glayout.substr( 0, 9 ) == "-glayout=")
      {
        glayout = glayout.substr( 9, std::string::npos ).c_str();

        istringstream layout( glayout );
        int sumprocs=0;

        while (layout)
        {
          string s;
          if (!getline( layout, s, ',' )) break;
          grouplayout.push_back( atoi(s.c_str()) );
          sumprocs += atoi(s.c_str());
        }

        // final check whether a correct group layout is specified
        if(ngroup != int(grouplayout.size()) or size != sumprocs or ngroup<1)
        {
          if(myrank == (size-1)) // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
          {
            printf("\n\nNumber of procs (%d) and number of groups (%d) does not match given group layout! \n",size,ngroup);
            printf("Example mpirun -np 4 baci-release -ngroup=2 -glayout=1,3 \n");
            printf("Try again!\n");
          }
          MPI_Finalize();
          exit(1);
        }
      }
      // case without given group layout
      else
      {
        if(myrank == (size-1)) // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
        {
          printf("\n\n\nINFO: Group layout is not specified. Default is equal size of the groups.\n");
        }
        if ((size % ngroup) != 0 and myrank == (size-1))
        {
          printf("\n\nNumber of processors (%d) cannot be divided by the number of groups (%d)!",size,ngroup);
          printf("Try again!\n");
          MPI_Finalize();
          exit(1);
        }

        // equal size of the groups
        for(int k=0; k<ngroup; k++)
        {
          grouplayout.push_back(size/ngroup);
        }
      }

      // the color is specified: procs are distributed to the groups with increasing global rank
      color = -1;
      int gsum = 0;
      do
      {
        color++;
        gsum += grouplayout[color];
      }
      while(gsum <= myrank);

      cout << "Nested parallelism layout: Global rank: " << myrank <<" is in group: "<< color << endl;

    } // end if (numgroup.substr( 0, 8 ) == "-ngroup=")


    //----------------------------------------------------------------
    // nested parallelism type
    //----------------------------------------------------------------
    std::string nptype(conf[i]);
    if (nptype.substr( 0, 8 ) == "-nptype=")
    {
      nptypeisset = true;
      nptype = nptype.substr( 8, std::string::npos ).c_str();
      if(nptype == "copyDatFile")
        npType = copy_dat_file;
      else if(nptype == "everyGroupReadDatFile")
        npType = every_group_read_dat_file;
      else if(nptype == "separateDatFiles")
        npType = separate_dat_files;
      else
      {
        if (myrank == (size-1)) // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
        {
          printf("\n\nOnly copyDatFile, everyGroupReadDatFile and separateDatFiles is available for nptype=  \n\n");
          printf("Try again!\n");
        }
        MPI_Finalize();
        exit(1);
      }
    }

  } // end for(int i=0; i<int(conf.size()); i++)


  if( (int(conf.size()) > 1)  and  (ngroupisset == false or nptypeisset == false ) )
  {
    if (myrank == (size-1)) // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
    {
      printf("\n\nAt least -nptype= and -ngroup= must be specified for nested parallelism. -glayout is optional (behind -ngroup).  \n\n");
      printf("Try again!\n");
    }
    MPI_Finalize();
    exit(1);
  }


  // do the splitting of the communicator
  MPI_Comm  mpi_local_comm;
  MPI_Comm_split(MPI_COMM_WORLD,color,myrank,&mpi_local_comm);

  Teuchos::RCP<Epetra_Comm> lcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));

  // the global communicator is created
  MPI_Comm mpi_global_comm;

  if(ngroup == 1)
  {
    mpi_global_comm = mpi_local_comm;
    // TODO: INCA needs color = 1 and BACI needs color = 0, then the proceeding line can be removed
    color = 0;
  }
  else
  {
    // TODO: consider a second executable that is coupled to BACI in case of nested parallelism
    // TODO: the procs owned by another executable have to be removed from world_group, e.g. MPI_Group_excl
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm_create(MPI_COMM_WORLD, world_group, &mpi_global_comm);
    MPI_Group_free(&world_group);
  }

  Teuchos::RCP<Epetra_Comm> gcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_global_comm));

  // mapping of local proc ids to global proc ids
  std::map<int, int> lpidgpid;
  int localsize = lcomm->NumProc();
  for(int lpid=0; lpid<localsize; lpid++)
  {
    lpidgpid[lpid]=gcomm->MyPID()-lcomm->MyPID()+lpid;
  }

  // nested parallelism group is given to the global problem
  DRT::Problem::Instance()->NPGroup(color, ngroup, lpidgpid, lcomm, gcomm, npType);

  return;
}


/*----------------------------------------------------------------------*
 | constructor nested parallelism group                     ghamm 03/12 |
 *----------------------------------------------------------------------*/
COMM_UTILS::NestedParGroup::NestedParGroup(
  int groupId,
  int ngroup,
  std::map<int, int> lpidgpid,
  Teuchos::RCP<Epetra_Comm> lcomm,
  Teuchos::RCP<Epetra_Comm> gcomm,
  NP_TYPE npType
  ) :
  groupId_(groupId),
  ngroup_(ngroup),
  lpidgpid_(lpidgpid),
  lcomm_(lcomm),
  gcomm_(gcomm),
  npType_(npType)
{
  return;
}


/*----------------------------------------------------------------------*
 | local proc id  of global proc id is returned             ghamm 03/12 |
 *----------------------------------------------------------------------*/
int COMM_UTILS::NestedParGroup::LPID(int GPID)
{
  std::map<int, int>::iterator it = lpidgpid_.begin();
  while(it != lpidgpid_.end())
  {
      if(it->second == GPID)
        return it->first;
      ++it;
  }
  // if GPID is not part of the current group
  printf("\n\n\nERROR: GPID (%d) is not in this group (%d) \n\n\n\n", GPID, groupId_);
  MPI_Abort(rcp_dynamic_cast<Epetra_MpiComm>(gcomm_,true)->GetMpiComm(),EXIT_FAILURE);
  exit(1);

  return -1;
}


/*----------------------------------------------------------------------*
 | broadcast all discretizations from bcast to all other groups gee 03/12 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::BroadcastDiscretizations(const int bgroup)
{
  Teuchos::RCP<DRT::Problem> problem = DRT::Problem::Instance();
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = problem->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  
  int sbcaster = -1;
  int bcaster = -1;
  int numfield = 0;
  vector<int> numdis;
  // only proc 0 of group bgroup
  if (group->GroupId()==bgroup && lcomm->MyPID()==0) 
  {
    sbcaster = group->GPID(0);
    numfield = problem->NumFields();
    numdis.resize(numfield);
    for (int i=0; i<numfield; ++i) 
      numdis[i] = problem->NumDis(i);
  }
  gcomm->MaxAll(&sbcaster,&bcaster,1);
  gcomm->Broadcast(&numfield,1,bcaster);
  numdis.resize(numfield);
  gcomm->Broadcast(&numdis[0],numfield,bcaster);

  for (int i=0; i<numfield; ++i)
  {
    for (int j=0; j<numdis[i]; ++j)
    {
      RCP<DRT::Discretization> dis = Teuchos::null;
      string name;
      string distype;
      vector<char> data;
      if (gcomm->MyPID()==bcaster) 
      {
        DRT::Container cont;
        dis = problem->Dis(i,j);
        name = dis->Name();
        distype = problem->SpatialApproximation();
        cont.Add("disname",name);
        cont.Add("distype",distype);
        DRT::PackBuffer buffer;
        cont.Pack(buffer);
        buffer.StartPacking();
        cont.Pack(buffer);
        std::swap(data,buffer());
      }
      // create map to export from proc 0 of group bgroup to all
      {
        int snummyelements = 0;
        int rnummyelements = 1;
        int myelements = 0;
        if (gcomm->MyPID()==bcaster) snummyelements = 1;
        Epetra_Map source(-1,snummyelements,&myelements,0,*gcomm);
        Epetra_Map target(-1,rnummyelements,&myelements,0,*gcomm);
        DRT::Exporter exporter(source,target,*gcomm);
        map<int,vector<char> > smap;
        if (gcomm->MyPID()==bcaster) smap[0] = data;
        exporter.Export<char>(smap);
        data = smap[0];
      }
      DRT::Container cont;
      vector<char> singledata;
      vector<char>::size_type index = 0;
      cont.ExtractfromPack(index,data,singledata);
      cont.Unpack(singledata);
      const string* rname = cont.Get<string>("disname");
      name = *rname;
      const string* rdistype = cont.Get<string>("distype");
      distype = *rdistype;
      // allocate or get the discretization
      if (group->GroupId()==bgroup) dis = problem->Dis(i,j);
      else
      {
        if (distype=="Nurbs")
          dis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(name,lcomm));
        else
          dis = Teuchos::rcp(new DRT::Discretization(name,lcomm));
      }
      // copy the discretization to the other groups
      for (int k=0; k<group->NumGroups(); ++k)
      {
        if (k==bgroup) continue; // broadcasting group does not copy to itself
        int color = MPI_UNDEFINED;
        if (group->GroupId()==k || bgroup==group->GroupId())
          color = 1;
        MPI_Comm intercomm;
        Epetra_MpiComm* mpicomm = dynamic_cast<Epetra_MpiComm*>(gcomm.get());
        if (!mpicomm) dserror("dyncast failed");
        MPI_Comm_split(mpicomm->Comm(),color,gcomm->MyPID(),&intercomm);
        
        
        if (group->GroupId()==k || bgroup==group->GroupId())
        {
          Teuchos::RCP<Epetra_MpiComm> icomm = Teuchos::rcp(new Epetra_MpiComm(intercomm));
          NPDuplicateDiscretization(bgroup,k,group,dis,icomm);
          icomm = Teuchos::null;
          MPI_Comm_free(&intercomm);
        }
        fflush(stdout);
        gcomm->Barrier();
      }


    } // for (int j=0; j<numdis[i]; ++j)
  } // for (int i=0; i<numfield; ++i)

  fflush(stdout);
  gcomm->Barrier();
  printf("gpid %d numfield %d numdis %d\n",gcomm->MyPID(),numfield,numdis[0]);
  return;
}



/*----------------------------------------------------------------------*
 | distribute a discretization from one group to one other    gee 03/12 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::NPDuplicateDiscretization(
                                 const int sgroup, 
                                 const int rgroup,
                                 Teuchos::RCP<NestedParGroup> group,
                                 Teuchos::RCP<DRT::Discretization> dis,
                                 Teuchos::RCP<Epetra_MpiComm> icomm)
{
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
#if 0
  vector<int> sbcaster(2,-1);
  vector<int> bcaster(2,-1);
  vector<int> sgsize(2,-1);
  vector<int> gsize(2,-1);
  if (group->GroupId()==sgroup && lcomm->MyPID()==0)
  {
    sbcaster[0] = icomm->MyPID();
    sgsize[0] = group->GroupSize();
  }
  if (group->GroupId()==rgroup && lcomm->MyPID()==0)
  {
    sbcaster[1] = icomm->MyPID();
    sgsize[1] = group->GroupSize();
  }
  icomm->MaxAll(&sbcaster[0],&bcaster[0],2);
  icomm->MaxAll(&sgsize[0],&gsize[0],2);
  printf("proc %d sgroup %d rgroup %d bcaster %d %d gsize %d %d\n",icomm->MyPID(),sgroup,rgroup,bcaster[0],bcaster[1],gsize[0],gsize[1]);
#endif

  // create a common discretization that we then fill with all stuff from sender group
  string name = dis->Name();  
  string type = DRT::Problem::Instance()->SpatialApproximation();
  RCP<DRT::Discretization> commondis;
  if (type=="Nurbs")
    commondis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(name,icomm));
  else
    commondis = Teuchos::rcp(new DRT::Discretization(name,icomm));

  // sender group fills commondis with elements and nodes and conditions
  // note that conditions are fully redundant
  map<int,vector<char> >        condmap;
  map<int,RCP<DRT::Container> > condnamemap;
  vector<int>                   myrowelements;
  if (group->GroupId()==sgroup)
  {
    if (!dis->Filled()) dis->FillComplete(false,false,false);
    
    myrowelements.resize(dis->ElementRowMap()->NumMyElements(),-1);
    
    // loop all my elements
    for (int i=0; i<dis->ElementRowMap()->NumMyElements(); ++i)
    {
      myrowelements[i] = dis->ElementRowMap()->GID(i);
      DRT::Element* ele = dis->lRowElement(i);
      if (myrowelements[i] != ele->Id()) dserror("Element global id mismatch");
      RCP<DRT::Element> newele = rcp(ele->Clone());
      newele->SetOwner(icomm->MyPID());
      commondis->AddElement(newele);
    }
    // loop all my nodes
    for (int i=0; i<dis->NodeRowMap()->NumMyElements(); ++i)
    {
      DRT::Node* node = dis->lRowNode(i);
      RCP<DRT::Node> newnode = rcp(node->Clone());
      newnode->SetOwner(icomm->MyPID());
      commondis->AddNode(newnode);
    }
    // loop conditions
    vector<string> condnames;
    dis->GetConditionNames(condnames);
    for (int i=0; i<(int)condnames.size(); ++i)
    {
      vector<DRT::Condition*> conds;
      dis->GetCondition(condnames[i],conds);
      RCP<vector<char> > data = dis->PackCondition(condnames[i]);
      if (lcomm->MyPID()==0) 
      {
        condmap[i] = *data;
        condnamemap[i] = rcp(new DRT::Container());
        condnamemap[i]->Add("condname",condnames[i]);
      }
      commondis->UnPackCondition(data,condnames[i]);
    }
  } // if (group->GroupId()==sgroup)
  
  // conditions have to be fully redundant, need to copy to all procs in intercomm
  {
    int nummyelements = (int)condmap.size();
    vector<int> myelements(nummyelements,-1);
    for (int i=0; i<nummyelements; ++i) myelements[i] = i;
    Epetra_Map smap(-1,nummyelements,&myelements[0],0,*icomm);
    nummyelements = smap.NumGlobalElements();
    myelements.resize(nummyelements);  
    for (int i=0; i<nummyelements; ++i) myelements[i] = i;
    Epetra_Map rmap(-1,nummyelements,&myelements[0],0,*icomm);
    DRT::Exporter exporter(smap,rmap,*icomm);
    exporter.Export<char>(condmap);
    exporter.Export<DRT::Container>(condnamemap);
  }

  // all rgroup procs loop and add condition
  if (group->GroupId()==rgroup)
  {
    map<int,vector<char> >::iterator fool1 = condmap.begin();
    map<int,RCP<DRT::Container> >::iterator fool2 = condnamemap.begin();
    for ( ; fool1 != condmap.end(); ++fool1)
    {
      const string* name = fool2->second->Get<string>("condname");
      commondis->UnPackCondition(rcp(&(fool1->second),false),*name);
      ++fool2;
    }
  }
  
  // finalize the commondis
  {
    RCP<Epetra_Map> roweles = rcp(new Epetra_Map(-1,(int)myrowelements.size(),&myrowelements[0],-1,*icomm));
    RCP<Epetra_Map> coleles;
    RCP<Epetra_Map> rownodes;
    RCP<Epetra_Map> colnodes;
    DRT::UTILS::PartUsingParMetis(commondis,roweles,rownodes,colnodes,icomm,false);
    commondis->BuildElementRowColumn(*rownodes,*colnodes,roweles,coleles);
    commondis->ExportRowNodes(*rownodes);
    commondis->ExportRowElements(*roweles);
    commondis->ExportColumnNodes(*colnodes);
    commondis->ExportColumnElements(*coleles);
    commondis->FillComplete(false,false,false);
  }
  
  // move everything from sender group to receiver group
  // build a target rowmap for nodes and elements
  
  
  // receiving group build its own discretization
    

  icomm->Barrier();
  return;
}
/*----------------------------------------------------------------------*/
