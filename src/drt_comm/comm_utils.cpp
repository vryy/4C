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

  if (argc > 4)
  {
    std::vector<int> grouplayout;
    for(int i=4; i<argc; i++)
    {
      //----------------------------------------------------------------
      // determine number of groups and the proc distribution
      //----------------------------------------------------------------
      std::string numgroup(argv[i]);
      if (numgroup.substr( 0, 7 ) == "ngroup=")
      {
        ngroup = atoi( numgroup.substr( 7, std::string::npos ).c_str() );

        // read out argument after ngroup=
        std::string glayout;
        if(i+1 < argc)
          glayout=argv[i+1];
        else
          glayout="dummy";

        // case with given group layout
        if (glayout.substr( 0, 8 ) == "glayout=")
        {
          glayout = glayout.substr( 8, std::string::npos ).c_str();

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
          if(ngroup != int(grouplayout.size()) || size != sumprocs)
          {
            dserror("Number of procs (%d) and number of groups (%d) does not match given group layout!",size,ngroup);
          }
        }
        // case without given group layout
        else
        {
          if(myrank == 0)
          {
            printf("INFO: Group layout is not specified. Default is equal size of the groups.\n");
          }
          if((size % ngroup) != 0)
            dserror("Number of processors (%d) cannot be divided by the number of groups (%d)!",size,ngroup);

          // equal size of the groups
          for(int k=0; k<ngroup; k++)
          {
            grouplayout.push_back(size/ngroup);
          }
        }
      }

      //----------------------------------------------------------------
      // nested parallelism type
      //----------------------------------------------------------------
      std::string nptype(argv[i]);
      if (nptype.substr( 0, 7 ) == "nptype=")
      {
        nptype = nptype.substr( 7, std::string::npos ).c_str();
        if(nptype == "copyDatFile");
          //do nothing so far
        else if(nptype == "microMacro")
          dserror("microMacro is not yet available for nptype");
      }
      //----------------------------------------------------------------
      // further input and output files
      //----------------------------------------------------------------
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

    cout<<"Nested parallelism layout:" << endl <<"Global rank: "<<myrank<<" is in group: "<<color<<endl;
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

  // insert the new communicators into the global problem
  DRT::Problem::Instance()->SetCommunicators(color, ngroup, lcomm, gcomm);

  return;
}


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
