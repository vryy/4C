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
#ifdef CCADISCRET

#include <Epetra_MpiComm.h>
#include <mpi.h>

#include <vector>
#include <sstream>
#include <string>

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include "comm_utils.H"
#include "../drt_lib/drt_globalproblem.H"


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


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
