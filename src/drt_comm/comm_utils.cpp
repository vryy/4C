/*----------------------------------------------------------------------*/
/*! \file

\brief Helper class for everything that deals with communication

\level 0

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include "comm_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_globalproblem_enums.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_utils_rebalancing.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>

#include <vector>
#include <sstream>
#include <string>
#include <iomanip>


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
  NestedParallelismType npType = no_nested_parallelism;

  // parse command line and separate configuration arguments
  std::vector<char*> conf(0);
  for (int i = 1; i < argc; i++)
  {
    std::string temp = argv[i];
    if (temp.substr(0, 1) == "-")
    {
      conf.push_back(argv[i]);
    }
  }

  // grouplayout will be filled accordingly to the user given input
  std::vector<int> grouplayout;
  bool ngroupisset = false;
  bool nptypeisset = false;
  for (int i = 0; i < int(conf.size()); i++)
  {
    // fill std::string with current argument
    std::string argument(conf[i]);
    //----------------------------------------------------------------
    // determine number of groups and the proc distribution
    //----------------------------------------------------------------
    if (argument.substr(0, 8) == "-ngroup=")
    {
      ngroupisset = true;
      ngroup = atoi(argument.substr(8, std::string::npos).c_str());

      // read out argument after ngroup=
      std::string glayout;
      if (i + 1 < int(conf.size()))
        glayout = conf[i + 1];
      else
        glayout = "dummy";

      // case with given group layout
      if (glayout.substr(0, 9) == "-glayout=")
      {
        glayout = glayout.substr(9, std::string::npos).c_str();

        std::istringstream layout(glayout);
        int sumprocs = 0;

        while (layout)
        {
          std::string s;
          if (!getline(layout, s, ',')) break;
          grouplayout.push_back(atoi(s.c_str()));
          sumprocs += atoi(s.c_str());
        }

        // final check whether a correct group layout is specified
        if (ngroup != int(grouplayout.size()) or size != sumprocs or ngroup < 1)
        {
          if (myrank == (size - 1))  // myrank == 0 is eventually not within BACI (i.e. coupling to
                                     // external codes)
          {
            printf(
                "\n\nNumber of procs (%d) and number of groups (%d) does not match given group "
                "layout! \n",
                size, ngroup);
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
        if (myrank ==
            (size -
                1))  // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
        {
          printf(
              "\n\n\nINFO: Group layout is not specified. Default is equal size of the groups.\n");
        }
        if ((size % ngroup) != 0)
        {
          if (myrank == (size - 1))
          {
            printf(
                "\n\nNumber of processors (%d) cannot be divided by the number of groups (%d)!\n",
                size, ngroup);
            printf("Try again!\n");
          }
          MPI_Finalize();
          exit(1);
        }

        // equal size of the groups
        for (int k = 0; k < ngroup; k++)
        {
          grouplayout.push_back(size / ngroup);
        }
      }

      // the color is specified: procs are distributed to the groups with increasing global rank
      color = -1;
      int gsum = 0;
      do
      {
        color++;
        gsum += grouplayout[color];
      } while (gsum <= myrank);

#ifdef DEBUG
      std::cout << "Nested parallelism layout: Global rank: " << myrank << " is in group: " << color
                << std::endl;
#endif

    }  // end if (argument.substr( 0, 8 ) == "-ngroup=")

    //----------------------------------------------------------------
    // nested parallelism type
    //----------------------------------------------------------------
    else if (argument.substr(0, 8) == "-nptype=")
    {
      nptypeisset = true;
      argument = argument.substr(8, std::string::npos).c_str();
      if (argument == "copyDatFile")
        npType = copy_dat_file;
      else if (argument == "everyGroupReadDatFile")
        npType = every_group_read_dat_file;
      else if (argument == "separateDatFiles")
        npType = separate_dat_files;
      else if (argument == "nestedMultiscale")
      {
        npType = separate_dat_files;
        // the color is specified: only two groups and group one (macro problem) is distributed over
        // all processors
        color = -1;
        if (myrank % (int)(size / grouplayout[0]) == 0 and
            myrank < (grouplayout[0] * (int)(size / grouplayout[0])))
          color = 0;
        else
          color = 1;
      }
      else if (argument.substr(0, 9) == "diffgroup")
      {
        npType = no_nested_parallelism;
        ngroup = 2;
        color = atoi(argument.substr(9, std::string::npos).c_str());
      }
      else
      {
        if (myrank ==
            (size -
                1))  // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
        {
          printf(
              "\n\nOnly copyDatFile, everyGroupReadDatFile and separateDatFiles is available for "
              "nptype=  \n\n");
          printf("Try again!\n");
        }
        MPI_Finalize();
        exit(1);
      }
    }

    //----------------------------------------------------------------
    // check for valid arguments that can be used in baci.cpp
    //----------------------------------------------------------------
    else if ((argument.substr(0, 9) != "-glayout=") and (argument.substr(0, 2) != "-v") and
             (argument.substr(0, 2) != "-h") and (argument.substr(0, 6) != "--help") and
             (argument.substr(0, 2) != "-p") and (argument.substr(0, 12) != "--parameters") and
             (argument.substr(0, 2) != "-d") and (argument.substr(0, 9) != "--datfile") and
             (argument.substr(0, 13) != "--interactive"))
    {
      printf(
          "\n\n You have specified an argument ( %s ) for BACI starting with a \"-\" that is not "
          "valid!\n",
          argument.c_str());
      printf("Please refer to ./baci-release --help and try again!\n");
      MPI_Finalize();
      exit(1);
    }

  }  // end for(int i=0; i<int(conf.size()); i++)


  if ((int(conf.size()) > 1) and (ngroupisset == false or nptypeisset == false))
  {
    if (myrank ==
        (size - 1))  // myrank == 0 is eventually not within BACI (i.e. coupling to external codes)
    {
      printf(
          "\n\nAt least -nptype= and -ngroup= must be specified for nested parallelism. -glayout "
          "is optional (behind -ngroup).  \n\n");
      printf("Try again!\n");
    }
    MPI_Finalize();
    exit(1);
  }


  // do the splitting of the communicator
  MPI_Comm mpi_local_comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, myrank, &mpi_local_comm);

  Teuchos::RCP<Epetra_Comm> lcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));

  // the global communicator is created
  Teuchos::RCP<Epetra_Comm> gcomm;

  if (ngroup == 1)
  {
    gcomm = lcomm;
    // TODO: INCA needs color = 1 and BACI needs color = 0, then the proceeding line can be removed
    color = 0;
  }
  else
  {
    // TODO: consider a second executable that is coupled to BACI in case of nested parallelism
    // TODO: the procs owned by another executable have to be removed from world_group, e.g.
    // MPI_Group_excl
    MPI_Comm mpi_global_comm;
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm_create(MPI_COMM_WORLD, world_group, &mpi_global_comm);
    MPI_Group_free(&world_group);

    gcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_global_comm));
  }

  // mapping of local proc ids to global proc ids
  std::map<int, int> lpidgpid;
  int localsize = lcomm->NumProc();
  for (int lpid = 0; lpid < localsize; lpid++)
  {
    lpidgpid[lpid] = gcomm->MyPID() - lcomm->MyPID() + lpid;
  }

  // nested parallelism group is given to the global problem
  DRT::Problem::Instance()->NPGroup(color, ngroup, lpidgpid, lcomm, gcomm, npType);

  // info for the nested parallelism user
  if (lcomm->MyPID() == 0 && ngroup > 1)
    printf("Nested parallelism layout: Group %d has %d processors.\n ", color, lcomm->NumProc());
  fflush(stdout);

  // for sync of output
  gcomm->Barrier();
  gcomm->Barrier();

  return;
}


/*----------------------------------------------------------------------*
 | constructor nested parallelism group                     ghamm 03/12 |
 *----------------------------------------------------------------------*/
COMM_UTILS::NestedParGroup::NestedParGroup(int groupId, int ngroup, std::map<int, int> lpidgpid,
    Teuchos::RCP<Epetra_Comm> lcomm, Teuchos::RCP<Epetra_Comm> gcomm, NestedParallelismType npType)
    : groupId_(groupId),
      ngroup_(ngroup),
      lpidgpid_(lpidgpid),
      lcomm_(lcomm),
      gcomm_(gcomm),
      subcomm_(Teuchos::null),
      npType_(npType)
{
  return;
}

/*----------------------------------------------------------------------*
 | copy constructor nested parallelism group                 kehl 08/16 |
 *----------------------------------------------------------------------*/
COMM_UTILS::NestedParGroup::NestedParGroup(NestedParGroup& npgroup)
{
  groupId_ = npgroup.GroupId();
  ngroup_ = npgroup.NumGroups();
  lcomm_ = npgroup.LocalComm();
  gcomm_ = npgroup.GlobalComm();
  subcomm_ = Teuchos::null;
  npType_ = npgroup.NpType();

  int size = npgroup.GroupSize();
  for (int i = 0; i < size; i++) lpidgpid_[i] = npgroup.GPID(i);

  return;
}


/*----------------------------------------------------------------------*
 | local proc id  of global proc id is returned             ghamm 03/12 |
 *----------------------------------------------------------------------*/
int COMM_UTILS::NestedParGroup::LPID(int GPID)
{
  std::map<int, int>::iterator it = lpidgpid_.begin();
  while (it != lpidgpid_.end())
  {
    if (it->second == GPID) return it->first;
    ++it;
  }
  // if GPID is not part of the current group
  printf("\n\n\nERROR: GPID (%d) is not in this group (%d) \n\n\n\n", GPID, groupId_);
  MPI_Abort(Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm_, true)->GetMpiComm(), EXIT_FAILURE);
  exit(1);

  return -1;
}


/*----------------------------------------------------------------------*
 | set sub communicator                                     ghamm 04/12 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::NestedParGroup::SetSubComm(Teuchos::RCP<Epetra_Comm> subcomm)
{
  subcomm_ = subcomm;
  return;
}

#if (0)
/*----------------------------------------------------------------------*
 | broadcast all discretizations from bcast to all other groups gee 03/12 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::BroadcastDiscretizations(const int bgroup)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = problem->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();

  int sbcaster = -1;
  int bcaster = -1;
  int numfield = 0;
  std::vector<int> numdis;
  // only proc 0 of group bgroup
  if (group->GroupId() == bgroup && lcomm->MyPID() == 0)
  {
    sbcaster = group->GPID(0);
    numfield = problem->NumFields();
    //    numdis.resize(numfield);
    //    for (int i=0; i<numfield; ++i)
    //      numdis[i] = 1; //obsolete! problem->NumDis(i);
  }
  gcomm->MaxAll(&sbcaster, &bcaster, 1);  // communicate proc null of bgroup
  gcomm->Broadcast(&numfield, 1, bcaster);
  //  numdis.resize(numfield);
  //  gcomm->Broadcast(&numdis[0],numfield,bcaster);

  const std::vector<std::string> disnames = DRT::Problem::Instance()->GetDisNames();

  for (int i = 0; i < numfield; ++i)
  {
    //    for (int j=0; j<numdis[i]; ++j)
    //    {
    Teuchos::RCP<DRT::Discretization> dis = Teuchos::null;
    std::string name;
    ShapeFunctionType distype;
    std::vector<char> data;
    if (gcomm->MyPID() == bcaster)
    {
      DRT::Container cont;
      dis = problem->GetDis(disnames[i]);
      name = dis->Name();
      distype = problem->SpatialApproximationType();
      cont.Add("disname", name);
      cont.Add("distype", (int)distype);
      DRT::PackBuffer buffer;
      cont.Pack(buffer);
      buffer.StartPacking();
      cont.Pack(buffer);
      std::swap(data, buffer());
    }
    // create map to export from proc 0 of group bgroup to all
    {
      int snummyelements = 0;
      int rnummyelements = 1;
      int myelements = 0;
      if (gcomm->MyPID() == bcaster) snummyelements = 1;
      Epetra_Map source(-1, snummyelements, &myelements, 0, *gcomm);
      Epetra_Map target(-1, rnummyelements, &myelements, 0, *gcomm);
      DRT::Exporter exporter(source, target, *gcomm);
      std::map<int, std::vector<char>> smap;
      if (gcomm->MyPID() == bcaster) smap[0] = data;
      exporter.Export<char>(smap);
      data = smap[0];
    }
    DRT::Container cont;
    std::vector<char> singledata;
    std::vector<char>::size_type index = 0;
    cont.ExtractfromPack(index, data, singledata);
    cont.Unpack(singledata);
    const std::string* rname = cont.Get<std::string>("disname");
    name = *rname;
    distype = (ShapeFunctionType)cont.GetInt("distype");
    // allocate or get the discretization
    if (group->GroupId() == bgroup)
    {
      dis = problem->GetDis(disnames[i]);
    }
    else
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          dis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(name, lcomm));
          break;
        }
        default:
        {
          dis = Teuchos::rcp(new DRT::Discretization(name, lcomm));
          break:
        }
      }
    }
    // copy the discretization to the other groups
    for (int k = 0; k < group->NumGroups(); ++k)
    {
      if (k == bgroup) continue;  // broadcasting group does not copy to itself
      int color = MPI_UNDEFINED;
      if (group->GroupId() == k || bgroup == group->GroupId()) color = 1;
      MPI_Comm intercomm;
      Epetra_MpiComm* mpicomm = dynamic_cast<Epetra_MpiComm*>(gcomm.get());
      if (!mpicomm) dserror("dyncast failed");
      MPI_Comm_split(mpicomm->Comm(), color, gcomm->MyPID(), &intercomm);


      if (group->GroupId() == k || bgroup == group->GroupId())
      {
        Teuchos::RCP<Epetra_MpiComm> icomm = Teuchos::rcp(new Epetra_MpiComm(intercomm));
        NPDuplicateDiscretization(bgroup, k, group, dis, icomm);
        icomm = Teuchos::null;
        MPI_Comm_free(&intercomm);
      }
      if (group->GroupId() == k)
      {
        const std::string disname = dis->Name();
        problem->AddDis(disname, dis);
      }
      gcomm->Barrier();
    }  // for (int k=0; k<group->NumGroups(); ++k)
       //    } // for (int j=0; j<numdis[i]; ++j)
  }    // for (int i=0; i<numfield; ++i)
  return;
}

#endif
/*----------------------------------------------------------------------*
 | broadcast all discretizations from group 0 to all other groups using
 | a ponzi scheme                                          biehler 05/13 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::BroadcastDiscretizations(int instance)
{
  // source of all discretizations is always group 0
  int bgroup = 0;
  // group to which discretizations are sent
  int tgroup = -1;

  DRT::Problem* problem = DRT::Problem::Instance(instance);
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = problem->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();

  // number of levels needed for the ponzi scheme
  int numlevels;
  // get next larger power of two
  numlevels = (int)(log(group->NumGroups()) / log(2) + 1);

  int sbcaster = -1;
  int bcaster = -1;
  int numfield = 0;
  std::vector<int> numdis;
  // only proc 0 of group bgroup
  if (group->GroupId() == bgroup && lcomm->MyPID() == 0)
  {
    sbcaster = group->GPID(0);
    numfield = problem->NumFields();
  }
  gcomm->MaxAll(&sbcaster, &bcaster, 1);  // communicate proc null of bgroup
  gcomm->Broadcast(&numfield, 1, bcaster);

  const std::vector<std::string> disnames = problem->GetDisNames();

  for (int i = 0; i < numfield; ++i)
  {
    Teuchos::RCP<DRT::Discretization> dis = Teuchos::null;
    std::string name;
    ShapeFunctionType distype;
    std::vector<char> data;
    if (gcomm->MyPID() == bcaster)
    {
      DRT::Container cont;
      dis = problem->GetDis(disnames[i]);
      name = dis->Name();
      distype = problem->SpatialApproximationType();
      cont.Add("disname", name);
      cont.Add("distype", (int)distype);
      DRT::PackBuffer buffer;
      cont.Pack(buffer);
      buffer.StartPacking();
      cont.Pack(buffer);
      std::swap(data, buffer());
    }
    // create map to export from proc 0 of group bgroup to all
    {
      int snummyelements = 0;
      int rnummyelements = 1;
      int myelements = 0;
      if (gcomm->MyPID() == bcaster) snummyelements = 1;
      Epetra_Map source(-1, snummyelements, &myelements, 0, *gcomm);
      Epetra_Map target(-1, rnummyelements, &myelements, 0, *gcomm);
      DRT::Exporter exporter(source, target, *gcomm);
      std::map<int, std::vector<char>> smap;
      if (gcomm->MyPID() == bcaster) smap[0] = data;
      exporter.Export<char>(smap);
      data = smap[0];
    }
    DRT::Container cont;
    std::vector<char> singledata;
    std::vector<char>::size_type index = 0;
    cont.ExtractfromPack(index, data, singledata);
    cont.Unpack(singledata);
    const std::string* rname = cont.Get<std::string>("disname");
    name = *rname;
    distype = (ShapeFunctionType)cont.GetInt("distype");
    // allocate or get the discretization
    if (group->GroupId() == bgroup)
    {
      dis = problem->GetDis(disnames[i]);
    }
    else
    {
      switch (distype)
      {
        case ShapeFunctionType::shapefunction_nurbs:
        {
          dis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(name, lcomm));
          break;
        }
        default:
        {
          dis = Teuchos::rcp(new DRT::Discretization(name, lcomm));
          break;
        }
      }
    }
    // copy the discretization to the other groups
    for (int k = numlevels; k > 0; --k)
    {
      // reset source and target group in every step
      tgroup = -1;
      bgroup = 0;

      int color = MPI_UNDEFINED;
      gcomm->Barrier();

      // color all sending groups
      if (group->GroupId() % (int)(pow(2.0, k)) == 0)
      {
        color = group->GroupId() / k;
        bgroup = group->GroupId();
      }
      // color all receivig groups
      if ((group->GroupId() % (int)(pow(2.0, k - 1))) == 0 &&
          !(group->GroupId() % (int)(pow(2.0, k)) == 0))
      {
        color = (group->GroupId() - pow(2.0, k - 1)) / k;
        tgroup = group->GroupId();
      }
      // std::cout << " GroupID " << group->GroupId() << "color " << color << std::endl;
      gcomm->Barrier();

      MPI_Comm intercomm;
      Epetra_MpiComm* mpicomm = dynamic_cast<Epetra_MpiComm*>(gcomm.get());
      if (!mpicomm) dserror("dyncast failed");
      MPI_Comm_split(mpicomm->Comm(), color, gcomm->MyPID(), &intercomm);


      if (bgroup == group->GroupId() || tgroup == group->GroupId())
      {
        Teuchos::RCP<Epetra_MpiComm> icomm = Teuchos::rcp(new Epetra_MpiComm(intercomm));
        // if group have same size use simpler broadcasting approach to ensure same
        // parallel distribution across all groups
        if (lcomm->NumProc() * 2 == icomm->NumProc())
          NPDuplicateDiscretizationEqualGroupSize(bgroup, tgroup, group, dis, icomm);
        else
          NPDuplicateDiscretization(bgroup, tgroup, group, dis, icomm);
        icomm = Teuchos::null;
        MPI_Comm_free(&intercomm);
      }
      if (group->GroupId() == tgroup)
      {
        const std::string disname = dis->Name();
        problem->AddDis(disname, dis);
      }
      gcomm->Barrier();
    }  // for (int k=numlevels; k>0; --k)
  }    // for (int i=0; i<numfield; ++i)
  return;
}



/*----------------------------------------------------------------------*
 | distribute a discretization from one group to one other    gee 03/12 |
 *----------------------------------------------------------------------*/
void COMM_UTILS::NPDuplicateDiscretization(const int sgroup, const int rgroup,
    Teuchos::RCP<NestedParGroup> group, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<Epetra_MpiComm> icomm)
{
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();

  std::vector<int> sbcaster(2, -1);
  std::vector<int> bcaster(2, -1);
  std::vector<int> sgsize(2, -1);
  std::vector<int> gsize(2, -1);
  if (group->GroupId() == sgroup && lcomm->MyPID() == 0)
  {
    sbcaster[0] = icomm->MyPID();
    sgsize[0] = group->GroupSize();
  }
  if (group->GroupId() == rgroup && lcomm->MyPID() == 0)
  {
    sbcaster[1] = icomm->MyPID();
    sgsize[1] = group->GroupSize();
  }
  icomm->MaxAll(&sbcaster[0], &bcaster[0], 2);
  icomm->MaxAll(&sgsize[0], &gsize[0], 2);
  // printf("proc %d sgroup %d rgroup %d bcaster %d %d gsize %d
  // %d\n",icomm->MyPID(),sgroup,rgroup,bcaster[0],bcaster[1],gsize[0],gsize[1]);

  // create a common discretization that we then fill with all stuff from sender group
  std::string name = dis->Name();
  ShapeFunctionType type = DRT::Problem::Instance()->SpatialApproximationType();
  Teuchos::RCP<DRT::Discretization> commondis;
  switch (type)
  {
    case ShapeFunctionType::shapefunction_nurbs:
    {
      commondis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(name, icomm));
      dserror("For Nurbs this method needs additional features!");
      break;
    }
    default:
    {
      commondis = Teuchos::rcp(new DRT::Discretization(name, icomm));
      break;
    }
  }

  // --------------------------------------
  // sender group fills commondis with elements and nodes and conditions
  // note that conditions are fully redundant
  std::map<int, std::vector<char>> condmap;
  std::map<int, Teuchos::RCP<DRT::Container>> condnamemap;
  std::vector<int> myrowelements;
  if (group->GroupId() == sgroup)
  {
    if (!dis->Filled()) dis->FillComplete(false, false, false);

    myrowelements.resize(dis->ElementRowMap()->NumMyElements(), -1);

    // loop all my elements
    for (int i = 0; i < dis->ElementRowMap()->NumMyElements(); ++i)
    {
      myrowelements[i] = dis->ElementRowMap()->GID(i);
      DRT::Element* ele = dis->lRowElement(i);
      if (myrowelements[i] != ele->Id()) dserror("Element global id mismatch");
      Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
      newele->SetOwner(icomm->MyPID());
      commondis->AddElement(newele);
    }
    // loop all my nodes
    for (int i = 0; i < dis->NodeRowMap()->NumMyElements(); ++i)
    {
      DRT::Node* node = dis->lRowNode(i);
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
      newnode->SetOwner(icomm->MyPID());
      commondis->AddNode(newnode);
    }
    // loop conditions
    std::vector<std::string> condnames;
    dis->GetConditionNames(condnames);
    for (int i = 0; i < (int)condnames.size(); ++i)
    {
      std::vector<DRT::Condition*> conds;
      dis->GetCondition(condnames[i], conds);
      Teuchos::RCP<std::vector<char>> data = dis->PackCondition(condnames[i]);
      if (lcomm->MyPID() == 0)
      {
        condmap[i] = *data;
        condnamemap[i] = Teuchos::rcp(new DRT::Container());
        condnamemap[i]->Add("condname", condnames[i]);
      }
      commondis->UnPackCondition(data, condnames[i]);
    }
  }  // if (group->GroupId()==sgroup)

  // --------------------------------------
  // conditions have to be fully redundant, need to copy to all procs in intercomm
  {
    int nummyelements = (int)condmap.size();
    std::vector<int> myelements(nummyelements, -1);
    for (int i = 0; i < nummyelements; ++i) myelements[i] = i;
    Epetra_Map smap(-1, nummyelements, &myelements[0], 0, *icomm);
    nummyelements = smap.NumGlobalElements();
    myelements.resize(nummyelements);
    for (int i = 0; i < nummyelements; ++i) myelements[i] = i;
    Epetra_Map rmap(-1, nummyelements, &myelements[0], 0, *icomm);
    DRT::Exporter exporter(smap, rmap, *icomm);
    exporter.Export<char>(condmap);
    exporter.Export<DRT::Container>(condnamemap);
  }

  // all rgroup procs loop and add condition to their dis and to the commondis
  if (group->GroupId() == rgroup)
  {
    std::map<int, std::vector<char>>::iterator fool1 = condmap.begin();
    std::map<int, Teuchos::RCP<DRT::Container>>::iterator fool2 = condnamemap.begin();
    for (; fool1 != condmap.end(); ++fool1)
    {
      const std::string* name = fool2->second->Get<std::string>("condname");
      commondis->UnPackCondition(Teuchos::rcp(&(fool1->second), false), *name);
      dis->UnPackCondition(Teuchos::rcp(&(fool1->second), false), *name);
      ++fool2;
    }
  }

  //-------------------------------------- finalize the commondis
  {
    Teuchos::RCP<Epetra_Map> roweles =
        Teuchos::rcp(new Epetra_Map(-1, (int)myrowelements.size(), &myrowelements[0], -1, *icomm));
    Teuchos::RCP<Epetra_Map> coleles;
    Teuchos::RCP<Epetra_Map> rownodes;
    Teuchos::RCP<Epetra_Map> colnodes;
    DRT::UTILS::REBALANCING::ComputeRebalancedNodeMaps(
        commondis, roweles, rownodes, colnodes, lcomm, false, lcomm->NumProc());
    commondis->BuildElementRowColumn(*rownodes, *colnodes, roweles, coleles);
    commondis->ExportRowNodes(*rownodes);
    commondis->ExportRowElements(*roweles);
    commondis->ExportColumnNodes(*colnodes);
    commondis->ExportColumnElements(*coleles);
    commondis->FillComplete(false, false, false);
  }

  //-------------------------------------- build a target rowmap for elements
  std::vector<int> targetgids;
  {
    const Epetra_Map* elerowmap = commondis->ElementRowMap();
    int nummyelements = elerowmap->NumMyElements();
    targetgids.resize(nummyelements);
    const int* myglobalelements = elerowmap->MyGlobalElements();
    // receiver group keeps its own gids
    if (group->GroupId() == rgroup)
    {
      for (int i = 0; i < nummyelements; ++i) targetgids[i] = myglobalelements[i];
    }
    else
      targetgids.resize(0);  // sender group wants to get rid of everything

    // each proc of sender group broadcast its gids
    // They are then distributed equally in the receiver group
    for (int proc = 0; proc < icomm->NumProc(); ++proc)
    {
      // check whether proc is a 'sender' at all
      int willsend = 0;
      if (icomm->MyPID() == proc && group->GroupId() == sgroup) willsend = 1;
      icomm->Broadcast(&willsend, 1, proc);
      if (!willsend) continue;  // proc is not a sender, goto next proc

      // proc send his nummyelements
      int nrecv = nummyelements;
      icomm->Broadcast(&nrecv, 1, proc);
      // allocate recv buffer
      std::vector<int> recvbuff(nrecv, 0);
      if (icomm->MyPID() == proc)
        for (int i = 0; i < nrecv; ++i) recvbuff[i] = myglobalelements[i];
      icomm->Broadcast(&recvbuff[0], nrecv, proc);
      // receivers share these gids equally
      if (group->GroupId() == rgroup)
      {
        int myshare = nrecv / gsize[1];
        int rest = nrecv % gsize[1];
        // printf("sendproc %d proc %d have %d gids nrecv %d myshare %d rest
        // %d\n",proc,icomm->MyPID(),(int)targetgids.size(),nrecv,myshare,rest); fflush(stdout);
        for (int i = 0; i < myshare; ++i)
          targetgids.push_back(recvbuff[lcomm->MyPID() * myshare + i]);
        if (lcomm->MyPID() == lcomm->NumProc() - 1)
          for (int i = 0; i < rest; ++i)
            targetgids.push_back(recvbuff[lcomm->NumProc() * myshare + i]);
        // printf("sendproc %d proc %d have %d gids\n",proc,icomm->MyPID(),(int)targetgids.size());
        // fflush(stdout);
      }
    }  // for (int proc=0; proc<icomm->NumProc(); ++proc)
  }
  Epetra_Map targetrowele(-1, (int)targetgids.size(), &targetgids[0], 0, *icomm);

  // -------------------------------------- build target row map for nodes
  {
    const Epetra_Map* noderowmap = commondis->NodeRowMap();
    int nummyelements = noderowmap->NumMyElements();
    targetgids.resize(nummyelements);
    const int* myglobalelements = noderowmap->MyGlobalElements();
    // receiver group keeps its own gids
    if (group->GroupId() == rgroup)
    {
      for (int i = 0; i < nummyelements; ++i) targetgids[i] = myglobalelements[i];
    }
    else
      targetgids.resize(0);  // sender group wants to get rid of everything

    // each proc of sender group broadcast its gids
    // They are then distributed equally in the receiver group
    for (int proc = 0; proc < icomm->NumProc(); ++proc)
    {
      // check whether proc is a 'sender' at all
      int willsend = 0;
      if (icomm->MyPID() == proc && group->GroupId() == sgroup) willsend = 1;
      icomm->Broadcast(&willsend, 1, proc);
      if (!willsend) continue;  // proc is not a sender, goto next proc

      // proc send his nummyelements
      int nrecv = nummyelements;
      icomm->Broadcast(&nrecv, 1, proc);
      // allocate recv buffer
      std::vector<int> recvbuff(nrecv, 0);
      if (icomm->MyPID() == proc)
        for (int i = 0; i < nrecv; ++i) recvbuff[i] = myglobalelements[i];
      icomm->Broadcast(&recvbuff[0], nrecv, proc);
      // receivers share these gids equally
      if (group->GroupId() == rgroup)
      {
        int myshare = nrecv / gsize[1];
        int rest = nrecv % gsize[1];
        // printf("sendproc %d proc %d have %d gids nrecv %d myshare %d rest
        // %d\n",proc,icomm->MyPID(),(int)targetgids.size(),nrecv,myshare,rest); fflush(stdout);
        for (int i = 0; i < myshare; ++i)
          targetgids.push_back(recvbuff[lcomm->MyPID() * myshare + i]);
        if (lcomm->MyPID() == lcomm->NumProc() - 1)
          for (int i = 0; i < rest; ++i)
            targetgids.push_back(recvbuff[lcomm->NumProc() * myshare + i]);
        // printf("sendproc %d proc %d have %d gids\n",proc,icomm->MyPID(),(int)targetgids.size());
        // fflush(stdout);
      }
    }  // for (int proc=0; proc<icomm->NumProc(); ++proc)
  }
  Epetra_Map targetrownode(-1, (int)targetgids.size(), &targetgids[0], 0, *icomm);

  // -------------- perform export of elements and nodes
  commondis->ExportRowElements(targetrowele);
  commondis->ExportRowNodes(targetrownode);

  //---------------loop elements and nodes and copy them to the rgroup discretization
  if (group->GroupId() == rgroup)
  {
    for (int i = 0; i < targetrowele.NumMyElements(); ++i)
    {
      DRT::Element* ele = commondis->gElement(targetrowele.GID(i));
      Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
      newele->SetOwner(lcomm->MyPID());
      dis->AddElement(newele);
    }
    for (int i = 0; i < targetrownode.NumMyElements(); ++i)
    {
      DRT::Node* node = commondis->gNode(targetrownode.GID(i));
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
      newnode->SetOwner(lcomm->MyPID());
      dis->AddNode(newnode);
    }
  }

  //------------------------------------------- complete the discretization on the rgroup
  if (group->GroupId() == rgroup)
  {
    Teuchos::RCP<Epetra_Map> roweles = Teuchos::rcp(new Epetra_Map(
        -1, targetrowele.NumMyElements(), targetrowele.MyGlobalElements(), -1, *lcomm));
    Teuchos::RCP<Epetra_Map> coleles;
    Teuchos::RCP<Epetra_Map> rownodes;
    Teuchos::RCP<Epetra_Map> colnodes;
    DRT::UTILS::REBALANCING::ComputeRebalancedNodeMaps(
        dis, roweles, rownodes, colnodes, lcomm, false, lcomm->NumProc());
    dis->BuildElementRowColumn(*rownodes, *colnodes, roweles, coleles);
    dis->ExportRowNodes(*rownodes);
    dis->ExportRowElements(*roweles);
    dis->ExportColumnNodes(*colnodes);
    dis->ExportColumnElements(*coleles);
    dis->FillComplete(true, true, true);
    dis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(dis)));
  }
  icomm->Barrier();
  return;
}



/*-----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void COMM_UTILS::NPDuplicateDiscretizationEqualGroupSize(const int sgroup, const int rgroup,
    Teuchos::RCP<NestedParGroup> group, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<Epetra_MpiComm> icomm)
{
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();

  std::vector<int> sbcaster(2, -1);
  std::vector<int> bcaster(2, -1);
  std::vector<int> sgsize(2, -1);
  std::vector<int> gsize(2, -1);
  if (group->GroupId() == sgroup && lcomm->MyPID() == 0)
  {
    sbcaster[0] = icomm->MyPID();
    sgsize[0] = group->GroupSize();
  }
  if (group->GroupId() == rgroup && lcomm->MyPID() == 0)
  {
    sbcaster[1] = icomm->MyPID();
    sgsize[1] = group->GroupSize();
  }
  icomm->MaxAll(&sbcaster[0], &bcaster[0], 2);
  icomm->MaxAll(&sgsize[0], &gsize[0], 2);

  // create a common discretization that we then fill with all stuff from sender group
  std::string name = dis->Name();
  ShapeFunctionType type = DRT::Problem::Instance()->SpatialApproximationType();
  Teuchos::RCP<DRT::Discretization> commondis;
  switch (type)
  {
    case ShapeFunctionType::shapefunction_nurbs:
    {
      commondis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(name, icomm));
      dserror("For Nurbs this method needs additional features!");
      break;
    }
    default:
    {
      commondis = Teuchos::rcp(new DRT::Discretization(name, icomm));
      break;
    }
  }


  // --------------------------------------
  // sender group fills commondis with elements and nodes and conditions
  // note that conditions are fully redundant
  std::map<int, std::vector<char>> condmap;
  std::map<int, Teuchos::RCP<DRT::Container>> condnamemap;
  std::vector<int> myrowelements;
  if (group->GroupId() == sgroup)
  {
    if (!dis->Filled()) dis->FillComplete(false, false, false);

    myrowelements.resize(dis->ElementRowMap()->NumMyElements(), -1);

    // loop all my elements
    for (int i = 0; i < dis->ElementRowMap()->NumMyElements(); ++i)
    {
      myrowelements[i] = dis->ElementRowMap()->GID(i);
      DRT::Element* ele = dis->lRowElement(i);
      if (myrowelements[i] != ele->Id()) dserror("Element global id mismatch");
      Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
      newele->SetOwner(icomm->MyPID());
      commondis->AddElement(newele);
    }
    // loop all my nodes
    for (int i = 0; i < dis->NodeRowMap()->NumMyElements(); ++i)
    {
      DRT::Node* node = dis->lRowNode(i);
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
      newnode->SetOwner(icomm->MyPID());
      commondis->AddNode(newnode);
    }
    // loop conditions
    std::vector<std::string> condnames;
    dis->GetConditionNames(condnames);
    for (int i = 0; i < (int)condnames.size(); ++i)
    {
      std::vector<DRT::Condition*> conds;
      dis->GetCondition(condnames[i], conds);
      Teuchos::RCP<std::vector<char>> data = dis->PackCondition(condnames[i]);
      if (lcomm->MyPID() == 0)
      {
        condmap[i] = *data;
        condnamemap[i] = Teuchos::rcp(new DRT::Container());
        condnamemap[i]->Add("condname", condnames[i]);
      }
      commondis->UnPackCondition(data, condnames[i]);
    }
  }  // if (group->GroupId()==sgroup)

  // --------------------------------------
  // conditions have to be fully redundant, need to copy to all procs in intercomm
  {
    int nummyelements = (int)condmap.size();
    std::vector<int> myelements(nummyelements, -1);
    for (int i = 0; i < nummyelements; ++i) myelements[i] = i;
    Epetra_Map smap(-1, nummyelements, &myelements[0], 0, *icomm);
    nummyelements = smap.NumGlobalElements();
    myelements.resize(nummyelements);
    for (int i = 0; i < nummyelements; ++i) myelements[i] = i;
    Epetra_Map rmap(-1, nummyelements, &myelements[0], 0, *icomm);
    DRT::Exporter exporter(smap, rmap, *icomm);
    exporter.Export<char>(condmap);
    exporter.Export<DRT::Container>(condnamemap);
  }
  // all rgroup procs loop and add condition to their dis and to the commondis
  if (group->GroupId() == rgroup)
  {
    std::map<int, std::vector<char>>::iterator fool1 = condmap.begin();
    std::map<int, Teuchos::RCP<DRT::Container>>::iterator fool2 = condnamemap.begin();
    for (; fool1 != condmap.end(); ++fool1)
    {
      const std::string* name = fool2->second->Get<std::string>("condname");
      commondis->UnPackCondition(Teuchos::rcp(&(fool1->second), false), *name);
      dis->UnPackCondition(Teuchos::rcp(&(fool1->second), false), *name);
      ++fool2;
    }
  }
  commondis->SetupGhosting(true, true, true);
  commondis->FillComplete(false, false, false);
  //-------------------------------------- build a target rowmap for elements
  std::vector<int> targetgids(0);
  MPI_Comm mpi_icomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(icomm)->GetMpiComm();
  {
    const Epetra_Map* elerowmap = commondis->ElementRowMap();
    int nummyelements = elerowmap->NumMyElements();
    int* myglobalelements = elerowmap->MyGlobalElements();
    // each proc of sender group broadcast its gids
    // They are then distributed equally in the receiver group
    int tag;
    for (int proc = 0; proc < icomm->NumProc(); ++proc)
    {
      if (icomm->MyPID() == proc && group->GroupId() == sgroup)
      {
        tag = 1337;
        // compare data
        int lengthSend = nummyelements;
        // first: send length of data
        MPI_Send(&lengthSend, 1, MPI_INT, (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag,
            mpi_icomm);
        // second: send data
        tag = 2674;
        MPI_Send(&myglobalelements[0], lengthSend, MPI_INT,
            (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag, mpi_icomm);
      }
      if (icomm->MyPID() == proc && group->GroupId() == rgroup)
      {
        // first: receive length of data
        int lengthRecv = 0;
        MPI_Status status;
        tag = 1337;
        MPI_Recv(&lengthRecv, 1, MPI_INT, (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag,
            mpi_icomm, &status);
        if (lengthRecv == 0) dserror("Length of data received from second run is zero.");
        // second: receive data
        tag = 2674;
        targetgids.resize(lengthRecv);
        MPI_Recv(&targetgids[0], lengthRecv, MPI_INT,
            (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag, mpi_icomm, &status);
      }
    }  // for (int proc=0; proc<icomm->NumProc(); ++proc)
  }
  Epetra_Map targetrowele(-1, (int)targetgids.size(), &targetgids[0], 0, *icomm);
  // -------------------------------------- build target row map for nodes
  {
    const Epetra_Map* noderowmap = commondis->NodeRowMap();
    int nummyelements = noderowmap->NumMyElements();
    // targetgids.resize(nummyelements);
    targetgids.resize(0);
    int* myglobalelements = noderowmap->MyGlobalElements();

    for (int proc = 0; proc < icomm->NumProc(); ++proc)
    {
      if (icomm->MyPID() == proc && group->GroupId() == sgroup)
      {
        int tag = 1337;
        // compare data
        int lengthSend = nummyelements;
        // first: send length of data
        tag = 1337;
        MPI_Send(&lengthSend, 1, MPI_INT, (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag,
            mpi_icomm);

        // second: send data
        tag = 2674;
        MPI_Send(&myglobalelements[0], lengthSend, MPI_INT,
            (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag, mpi_icomm);
      }
      if (icomm->MyPID() == proc && group->GroupId() == rgroup)
      {
        // first: receive length of data
        int lengthRecv = 0;

        MPI_Status status;
        int tag = 1337;
        MPI_Recv(&lengthRecv, 1, MPI_INT, (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag,
            mpi_icomm, &status);
        if (lengthRecv == 0) dserror("Length of data received from second run is zero.");
        // second: receive data
        tag = 2674;
        targetgids.resize(lengthRecv);
        MPI_Recv(&targetgids[0], lengthRecv, MPI_INT,
            (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag, mpi_icomm, &status);
      }
    }  // for (int proc=0; proc<icomm->NumProc(); ++proc)
  }
  Epetra_Map targetrownode(-1, (int)targetgids.size(), &targetgids[0], 0, *icomm);
  // -------------- perform export of elements and nodes
  commondis->ExportRowElements(targetrowele);
  commondis->ExportRowNodes(targetrownode);


  //---------------loop elements and nodes and copy them to the rgroup discretization
  if (group->GroupId() == rgroup)
  {
    for (int i = 0; i < targetrowele.NumMyElements(); ++i)
    {
      DRT::Element* ele = commondis->gElement(targetrowele.GID(i));
      Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
      newele->SetOwner(lcomm->MyPID());
      dis->AddElement(newele);
    }
    for (int i = 0; i < targetrownode.NumMyElements(); ++i)
    {
      DRT::Node* node = commondis->gNode(targetrownode.GID(i));
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
      newnode->SetOwner(lcomm->MyPID());
      dis->AddNode(newnode);
    }
  }
  //-------------------------------- complete the discretization on the rgroup
  if (group->GroupId() == rgroup)
  {
    dis->SetupGhosting(true, true, true);
    dis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(dis)));
  }
  icomm->Barrier();
  return;
}
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compare vectors from different parallel baci runs        ghamm 10/14 |
 |                                                                      |
 | You can add COMM_UTILS::CompareVectors([insert any Epetra_MultiVector|
 | here], name_of_vector); in your code which will lead to a comparison |
 | of the given vector for different executables and/or configurations  |
 | command for starting this feature:                                   |
 | mpirun -np 1 ./baci-release -nptype=diffgroup0 input.dat xxx_ser : -np 3 ./other-baci-release
 -nptype=diffgroup1 other-input.dat xxx_par | | | Do not forget to include the header ( #include
 "../drt_comm/         | | comm_utils.H" ) for the compare method, otherwise it won't compile.  | |
 | | Further nice options are that you can compare results from different | | executables used for
 running the same simulation                     | | Note: you need to add the CompareVectors method
 in both executables  | | at the same position in the code                                     |
 *----------------------------------------------------------------------*/
void COMM_UTILS::CompareVectors(
    Teuchos::RCP<const Epetra_MultiVector> vec, const char* name, double tol /*= 1.0e-14*/
)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = problem->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  MPI_Comm mpi_lcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(lcomm)->GetMpiComm();
  MPI_Comm mpi_gcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm)->GetMpiComm();

  int result = -1;
  MPI_Comm_compare(mpi_gcomm, mpi_lcomm, &result);
  if (result == 0)
  {
    IO::cout << "WARNING:: Vectors " << name
             << " cannot be compared because second baci run is missing" << IO::endl;
    return;
  }

  // do stupid conversion from Epetra_BlockMap to Epetra_Map
  const Epetra_BlockMap& vecblockmap = vec->Map();
  Teuchos::RCP<Epetra_Map> vecmap = Teuchos::rcp(new Epetra_Map(vecblockmap.NumGlobalElements(),
      vecblockmap.NumMyElements(), vecblockmap.MyGlobalElements(), 0, vec->Comm()));

  // gather data of vector to compare on gcomm proc 0 and last gcomm proc
  Teuchos::RCP<Epetra_Map> proc0map;
  if (lcomm->MyPID() == gcomm->MyPID())
    proc0map = LINALG::AllreduceOverlappingEMap(*vecmap, 0);
  else
    proc0map = LINALG::AllreduceOverlappingEMap(*vecmap, lcomm->NumProc() - 1);

  // export full vectors to the two desired processors
  Teuchos::RCP<Epetra_MultiVector> fullvec =
      Teuchos::rcp(new Epetra_MultiVector(*proc0map, vec->NumVectors(), true));
  LINALG::Export(*vec, *fullvec);

  const int myglobalrank = gcomm->MyPID();
  // last proc in gcomm sends its data to proc 0 which does the comparison
  if (myglobalrank == 0)
  {
    // compare names
    int lengthRecv = 0;
    std::vector<char> receivename;
    MPI_Status status;
    // first: receive length of name
    int tag = 1336;
    MPI_Recv(&lengthRecv, 1, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);
    if (lengthRecv == 0) dserror("Length of name received from second run is zero.");

    // second: receive name
    tag = 2672;
    receivename.resize(lengthRecv);
    MPI_Recv(&receivename[0], lengthRecv, MPI_CHAR, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // do comparison of names
    if (std::strcmp(name, &receivename[0]))
      dserror(
          "comparison of different vectors: group 0 (%s) and group 1 (%s)", name, &receivename[0]);

    // compare data
    lengthRecv = 0;
    std::vector<double> receivebuf;
    // first: receive length of data
    tag = 1337;
    MPI_Recv(&lengthRecv, 1, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);
    // also enable comparison of empty vectors
    if (lengthRecv == 0 && fullvec->MyLength() != lengthRecv)
      dserror("Length of data received from second run is incorrect.");

    // second: receive data
    tag = 2674;
    receivebuf.resize(lengthRecv);
    MPI_Recv(&receivebuf[0], lengthRecv, MPI_DOUBLE, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // start comparison
    int mylength = fullvec->MyLength() * vec->NumVectors();
    if (mylength != lengthRecv)
      dserror("length of received data (%i) does not match own data (%i)", lengthRecv, mylength);

    double maxdiff = 0.0;
    for (int i = 0; i < mylength; ++i)
    {
      double difference = std::abs(fullvec->Values()[i] - receivebuf[i]);
      if (difference > tol)
      {
        std::stringstream diff;
        diff << std::scientific << std::setprecision(16) << maxdiff;
        std::cout << "vectors " << name << " do not match, difference in row "
                  << fullvec->Map().GID(i) << " between entries is: " << diff.str().c_str()
                  << std::endl;
      }
      maxdiff = std::max(maxdiff, difference);
    }
    if (maxdiff > tol)
    {
      std::stringstream diff;
      diff << std::scientific << std::setprecision(16) << maxdiff;
      dserror("vectors %s do not match, maximum difference between entries is: %s", name,
          diff.str().c_str());
    }
    else
    {
      IO::cout << "compared vectors " << name << " of length: " << mylength
               << " which are identical." << IO::endl;
    }
  }
  else if (myglobalrank == gcomm->NumProc() - 1)
  {
    // compare names
    // include terminating \0 of char array
    int lengthSend = std::strlen(name) + 1;
    // first: send length of name
    int tag = 1336;
    MPI_Send(&lengthSend, 1, MPI_INT, 0, tag, mpi_gcomm);

    // second: send name
    tag = 2672;
    MPI_Send(const_cast<char*>(name), lengthSend, MPI_CHAR, 0, tag, mpi_gcomm);

    // compare data
    lengthSend = fullvec->MyLength() * vec->NumVectors();
    // first: send length of data
    tag = 1337;
    MPI_Send(&lengthSend, 1, MPI_INT, 0, tag, mpi_gcomm);

    // second: send data
    tag = 2674;
    MPI_Send(&fullvec->Values()[0], lengthSend, MPI_DOUBLE, 0, tag, mpi_gcomm);
  }

  // force all procs to stay here until proc 0 has checked the vectors
  gcomm->Barrier();

  return;
}


/*----------------------------------------------------------------------*
 | see comment of CompareVectors above                      ghamm 12/14 |
 | from LINALG::SparseOperator to CrsMatrix, just do:                   |
 | Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(yoursparseoperator)->EpetraMatrix()
 *----------------------------------------------------------------------*/
void COMM_UTILS::CompareSparseMatrices(
    Teuchos::RCP<Epetra_CrsMatrix> matrix, const char* name, double tol /*= 1.0e-14*/
)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = problem->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  MPI_Comm mpi_lcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(lcomm)->GetMpiComm();
  MPI_Comm mpi_gcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm)->GetMpiComm();
  const int myglobalrank = gcomm->MyPID();

  int result = -1;
  MPI_Comm_compare(mpi_gcomm, mpi_lcomm, &result);
  if (result == 0)
  {
    IO::cout << "WARNING:: Matrices " << name
             << " cannot be compared because second baci run is missing" << IO::endl;
    return;
  }

  const Epetra_Map& originalmap = matrix->RowMap();

  // gather data of vector to compare on gcomm proc 0 and last gcomm proc
  Teuchos::RCP<Epetra_Map> serialmap;
  if (lcomm->MyPID() == gcomm->MyPID())
    serialmap = LINALG::AllreduceOverlappingEMap(originalmap, 0);
  else
    serialmap = LINALG::AllreduceOverlappingEMap(originalmap, lcomm->NumProc() - 1);

  // export full matrices to the two desired processors
  Teuchos::RCP<Epetra_Import> serialimporter =
      Teuchos::rcp(new Epetra_Import(*serialmap, originalmap));
  Teuchos::RCP<Epetra_CrsMatrix> serialCrsMatrix =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *serialmap, 0));
  serialCrsMatrix->Import(*matrix, *serialimporter, Insert);
  serialCrsMatrix->FillComplete();

  // fill data of matrices to container which can be easily communicated via MPI
  std::vector<int> data_indices;
  data_indices.reserve(serialCrsMatrix->NumMyNonzeros() * 2);
  std::vector<double> data_values;
  data_values.reserve(serialCrsMatrix->NumMyNonzeros());
  if (myglobalrank == 0 || myglobalrank == gcomm->NumProc() - 1)
  {
    for (int i = 0; i < serialmap->NumMyElements(); ++i)
    {
      int rowgid = serialmap->GID(i);
      int NumEntries;
      double* Values;
      int* Indices;
      int err = serialCrsMatrix->ExtractMyRowView(i, NumEntries, Values, Indices);
      if (err != 0) dserror("ExtractMyRowView error: %d", err);

      for (int j = 0; j < NumEntries; ++j)
      {
        // store row and col gid in order to compare them on proc 0 and for detailed error output
        // information
        data_indices.push_back(rowgid);
        data_indices.push_back(Indices[j]);
        data_values.push_back(Values[j]);
      }
    }
  }

  // last proc in gcomm sends its data to proc 0 which does the comparison
  if (myglobalrank == 0)
  {
    // compare names
    int lengthRecv = 0;
    std::vector<char> receivename;
    MPI_Status status;
    // first: receive length of name
    int tag = 1336;
    MPI_Recv(&lengthRecv, 1, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);
    if (lengthRecv == 0) dserror("Length of name received from second run is zero.");

    // second: receive name
    tag = 2672;
    receivename.resize(lengthRecv);
    MPI_Recv(&receivename[0], lengthRecv, MPI_CHAR, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // do comparison of names
    if (std::strcmp(name, &receivename[0]))
      dserror(
          "comparison of different vectors: group 0 (%s) and group 1 (%s)", name, &receivename[0]);

    // compare data: indices
    lengthRecv = 0;
    std::vector<int> receivebuf_indices;
    // first: receive length of data
    tag = 1337;
    MPI_Recv(&lengthRecv, 1, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);
    // also enable comparison of empty matrices
    if (lengthRecv == 0 && (int)data_indices.size() != lengthRecv)
      dserror("Length of data received from second run is incorrect.");

    // second: receive data
    tag = 2674;
    receivebuf_indices.resize(lengthRecv);
    MPI_Recv(
        &receivebuf_indices[0], lengthRecv, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // start comparison
    int mylength = data_indices.size();
    if (mylength != lengthRecv)
      dserror("length of received data (%i) does not match own data (%i)", lengthRecv, mylength);

    double maxdiff = 0.0;
    for (int i = 0; i < mylength; ++i)
    {
      if (data_indices[i] != receivebuf_indices[i])
      {
        bool iscolindex = data_indices[i] % 2;
        dserror("%s index of matrix %s does not match: group 0 (%i) and group 1 (%i)",
            iscolindex == 0 ? "row" : "col", name, data_indices[i], receivebuf_indices[i]);
      }
    }
    IO::cout << "indices of compared matrices " << name << " of length: " << mylength
             << " are identical." << IO::endl;

    // compare data: values
    lengthRecv = 0;
    std::vector<double> receivebuf_values;
    // first: receive length of data
    tag = 1338;
    MPI_Recv(&lengthRecv, 1, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);
    // also enable comparison of empty matrices
    if (lengthRecv == 0 && (int)data_values.size() != lengthRecv)
      dserror("Length of data received from second run is incorrect.");

    // second: receive data
    tag = 2676;
    receivebuf_values.resize(lengthRecv);
    MPI_Recv(&receivebuf_values[0], lengthRecv, MPI_DOUBLE, gcomm->NumProc() - 1, tag, mpi_gcomm,
        &status);

    // start comparison
    mylength = data_values.size();
    if (mylength != lengthRecv)
      dserror("length of received data (%i) does not match own data (%i)", lengthRecv, mylength);

    maxdiff = 0.0;
    for (int i = 0; i < mylength; ++i)
    {
      double difference = std::abs(data_values[i] - receivebuf_values[i]);
      if (difference > tol)
      {
        std::stringstream diff;
        diff << std::scientific << std::setprecision(16) << maxdiff;
        std::cout << "matrices " << name << " do not match, difference in row "
                  << data_indices[2 * i] << " , col: " << data_indices[2 * i + 1]
                  << " between entries is: " << diff.str().c_str() << std::endl;
      }
      maxdiff = std::max(maxdiff, difference);
    }
    if (maxdiff > tol)
    {
      std::stringstream diff;
      diff << std::scientific << std::setprecision(16) << maxdiff;
      dserror("matrices %s do not match, maximum difference between entries is: %s in row", name,
          diff.str().c_str());
    }
    else
    {
      IO::cout << "values of compared matrices " << name << " of length: " << mylength
               << " are identical." << IO::endl;
    }
  }
  else if (myglobalrank == gcomm->NumProc() - 1)
  {
    // compare names
    // include terminating \0 of char array
    int lengthSend = std::strlen(name) + 1;
    // first: send length of name
    int tag = 1336;
    MPI_Send(&lengthSend, 1, MPI_INT, 0, tag, mpi_gcomm);

    // second: send name
    tag = 2672;
    MPI_Send(const_cast<char*>(name), lengthSend, MPI_CHAR, 0, tag, mpi_gcomm);

    // compare data: indices
    lengthSend = data_indices.size();
    // first: send length of data
    tag = 1337;
    MPI_Send(&lengthSend, 1, MPI_INT, 0, tag, mpi_gcomm);

    // second: send data
    tag = 2674;
    MPI_Send(&data_indices[0], lengthSend, MPI_INT, 0, tag, mpi_gcomm);

    // compare data: values
    lengthSend = data_values.size();
    // first: send length of data
    tag = 1338;
    MPI_Send(&lengthSend, 1, MPI_INT, 0, tag, mpi_gcomm);

    // second: send data
    tag = 2676;
    MPI_Send(&data_values[0], lengthSend, MPI_DOUBLE, 0, tag, mpi_gcomm);
  }

  // force all procs to stay here until proc 0 has checked the matrices
  gcomm->Barrier();

  return;
}
