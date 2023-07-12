/*----------------------------------------------------------------------*/
/*! \file

\brief Helper class for everything that deals with communication

\level 0

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 01/12 |
 *----------------------------------------------------------------------*/
#include "comm_utils.H"
#include "io_pstream.H"
#include "linalg_utils_densematrix_communication.H"
#include "linalg_utils_sparse_algebra_manipulation.H"

#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_MpiComm.h>

#include <vector>
#include <sstream>
#include <string>
#include <iomanip>


/*----------------------------------------------------------------------*
 | create communicator                                      ghamm 02/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<COMM_UTILS::Communicators> COMM_UTILS::CreateComm(std::vector<std::string> argv)
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
  std::vector<std::string> conf(0);
  for (std::size_t i = 1; i < argv.size(); i++)
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
  for (std::size_t i = 0; i < conf.size(); i++)
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
      if (i + 1 < conf.size())
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

  // nested parallelism group is created
  Teuchos::RCP<COMM_UTILS::Communicators> communicators =
      Teuchos::rcp(new COMM_UTILS::Communicators(color, ngroup, lpidgpid, lcomm, gcomm, npType));

  // info for the nested parallelism user
  if (lcomm->MyPID() == 0 && ngroup > 1)
    printf("Nested parallelism layout: Group %d has %d processors.\n ", color, lcomm->NumProc());
  fflush(stdout);

  // for sync of output
  gcomm->Barrier();
  gcomm->Barrier();

  return communicators;
}

/*----------------------------------------------------------------------*
 | constructor communicators                                ghamm 03/12 |
 *----------------------------------------------------------------------*/
COMM_UTILS::Communicators::Communicators(int groupId, int ngroup, std::map<int, int> lpidgpid,
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
 | local proc id  of global proc id is returned             ghamm 03/12 |
 *----------------------------------------------------------------------*/
int COMM_UTILS::Communicators::LPID(int GPID)
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
void COMM_UTILS::Communicators::SetSubComm(Teuchos::RCP<Epetra_Comm> subcomm)
{
  subcomm_ = subcomm;
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool COMM_UTILS::AreDistributedVectorsIdentical(const COMM_UTILS::Communicators& communicators,
    Teuchos::RCP<const Epetra_MultiVector> vec, const char* name, double tol /*= 1.0e-14*/
)
{
  Teuchos::RCP<Epetra_Comm> lcomm = communicators.LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = communicators.GlobalComm();
  MPI_Comm mpi_lcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(lcomm)->GetMpiComm();
  MPI_Comm mpi_gcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm)->GetMpiComm();

  int result = -1;
  MPI_Comm_compare(mpi_gcomm, mpi_lcomm, &result);
  if (result == 0)
  {
    IO::cout << "WARNING:: Vectors " << name
             << " cannot be compared because second baci run is missing" << IO::endl;
    return false;
  }

  // do stupid conversion from Epetra_BlockMap to Epetra_Map
  const Epetra_BlockMap& vecblockmap = vec->Map();
  Teuchos::RCP<Epetra_Map> vecmap = Teuchos::rcp(new Epetra_Map(vecblockmap.NumGlobalElements(),
      vecblockmap.NumMyElements(), vecblockmap.MyGlobalElements(), 0, vec->Comm()));

  // gather data of vector to compare on gcomm proc 0 and last gcomm proc
  Teuchos::RCP<Epetra_Map> proc0map;
  if (lcomm->MyPID() == gcomm->MyPID())
    proc0map = CORE::LINALG::AllreduceOverlappingEMap(*vecmap, 0);
  else
    proc0map = CORE::LINALG::AllreduceOverlappingEMap(*vecmap, lcomm->NumProc() - 1);

  // export full vectors to the two desired processors
  Teuchos::RCP<Epetra_MultiVector> fullvec =
      Teuchos::rcp(new Epetra_MultiVector(*proc0map, vec->NumVectors(), true));
  CORE::LINALG::Export(*vec, *fullvec);

  const int myglobalrank = gcomm->MyPID();
  double maxdiff = 0.0;
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
    MPI_Recv(
        receivename.data(), lengthRecv, MPI_CHAR, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // do comparison of names
    if (std::strcmp(name, receivename.data()))
      dserror("comparison of different vectors: communicators 0 (%s) and communicators 1 (%s)",
          name, receivename.data());

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
    MPI_Recv(
        receivebuf.data(), lengthRecv, MPI_DOUBLE, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // start comparison
    int mylength = fullvec->MyLength() * vec->NumVectors();
    if (mylength != lengthRecv)
      dserror("length of received data (%i) does not match own data (%i)", lengthRecv, mylength);

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
    if (maxdiff <= tol)
    {
      IO::cout << "compared vectors " << name << " of length: " << mylength
               << " which are identical." << IO::endl;
      result = 1;
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
    MPI_Send(fullvec->Values(), lengthSend, MPI_DOUBLE, 0, tag, mpi_gcomm);
  }

  // force all procs to stay here until proc 0 has checked the vectors
  gcomm->Broadcast(&maxdiff, 1, 0);
  if (maxdiff > tol)
  {
    std::stringstream diff;
    diff << std::scientific << std::setprecision(16) << maxdiff;
    dserror("vectors %s do not match, maximum difference between entries is: %s", name,
        diff.str().c_str());
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool COMM_UTILS::AreDistributedSparseMatricesIdentical(
    const COMM_UTILS::Communicators& communicators, Teuchos::RCP<Epetra_CrsMatrix> matrix,
    const char* name, double tol /*= 1.0e-14*/
)
{
  Teuchos::RCP<Epetra_Comm> lcomm = communicators.LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = communicators.GlobalComm();
  MPI_Comm mpi_lcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(lcomm)->GetMpiComm();
  MPI_Comm mpi_gcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm)->GetMpiComm();
  const int myglobalrank = gcomm->MyPID();

  int result = -1;
  MPI_Comm_compare(mpi_gcomm, mpi_lcomm, &result);
  if (result == 0)
  {
    IO::cout << "WARNING:: Matrices " << name
             << " cannot be compared because second baci run is missing" << IO::endl;
    return false;
  }

  const Epetra_Map& originalmap = matrix->RowMap();

  // gather data of vector to compare on gcomm proc 0 and last gcomm proc
  Teuchos::RCP<Epetra_Map> serialmap;
  if (lcomm->MyPID() == gcomm->MyPID())
    serialmap = CORE::LINALG::AllreduceOverlappingEMap(originalmap, 0);
  else
    serialmap = CORE::LINALG::AllreduceOverlappingEMap(originalmap, lcomm->NumProc() - 1);

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
  double maxdiff = 0.0;
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
    MPI_Recv(
        receivename.data(), lengthRecv, MPI_CHAR, gcomm->NumProc() - 1, tag, mpi_gcomm, &status);

    // do comparison of names
    if (std::strcmp(name, receivename.data()))
      dserror("comparison of different vectors: communicators 0 (%s) and communicators 1 (%s)",
          name, receivename.data());

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
    MPI_Recv(receivebuf_indices.data(), lengthRecv, MPI_INT, gcomm->NumProc() - 1, tag, mpi_gcomm,
        &status);

    // start comparison
    int mylength = data_indices.size();
    if (mylength != lengthRecv)
      dserror("length of received data (%i) does not match own data (%i)", lengthRecv, mylength);

    for (int i = 0; i < mylength; ++i)
    {
      if (data_indices[i] != receivebuf_indices[i])
      {
        bool iscolindex = data_indices[i] % 2;
        dserror(
            "%s index of matrix %s does not match: communicators 0 (%i) and communicators 1 (%i)",
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
    MPI_Recv(receivebuf_values.data(), lengthRecv, MPI_DOUBLE, gcomm->NumProc() - 1, tag, mpi_gcomm,
        &status);

    // start comparison
    mylength = data_values.size();
    if (mylength != lengthRecv)
      dserror("length of received data (%i) does not match own data (%i)", lengthRecv, mylength);

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
    if (maxdiff <= tol)
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
    MPI_Send(data_indices.data(), lengthSend, MPI_INT, 0, tag, mpi_gcomm);

    // compare data: values
    lengthSend = data_values.size();
    // first: send length of data
    tag = 1338;
    MPI_Send(&lengthSend, 1, MPI_INT, 0, tag, mpi_gcomm);

    // second: send data
    tag = 2676;
    MPI_Send(data_values.data(), lengthSend, MPI_DOUBLE, 0, tag, mpi_gcomm);
  }

  // force all procs to stay here until proc 0 has checked the matrices
  gcomm->Broadcast(&maxdiff, 1, 0);
  if (maxdiff > tol)
  {
    std::stringstream diff;
    diff << std::scientific << std::setprecision(16) << maxdiff;
    dserror("matrices %s do not match, maximum difference between entries is: %s in row", name,
        diff.str().c_str());
  }

  return true;
}
