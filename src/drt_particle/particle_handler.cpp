/*----------------------------------------------------------------------*/
/*!

\brief Handler to control particle simulations

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#include "particle_handler.H"

#include "../drt_binstrategy/drt_meshfree_multibin.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <unordered_set>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PARTICLE::ParticleHandler::ParticleHandler()
    : binstrategy_(Teuchos::null), myrank_(-1), bincolmap_(Teuchos::null)
{
  // empty constructor
}

void PARTICLE::ParticleHandler::Init(
    int myrank, Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy)
{
  binstrategy_ = binstrategy;
  myrank_ = myrank;
}

void PARTICLE::ParticleHandler::Setup()
{
  // so far nothing to do
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PARTICLE::ParticleHandler::ParticleHandler(const Epetra_Comm& comm)
    : binstrategy_(Teuchos::rcp(new BINSTRATEGY::BinningStrategy(comm))),
      myrank_(comm.MyPID()),
      bincolmap_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleHandler::DistributeParticlesToBins(
    Teuchos::RCP<Epetra_Map> const& particlerowmap)
{
  std::list<Teuchos::RCP<DRT::Node>> homelessparticles;
  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = binstrategy_->BinDiscret()->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node, false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBinsRoundRobin(homelessparticles);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleHandler::RemoveAllParticles()
{
  // 1st) loop over bins and remove initial particle info
  const int numrowbin = BinStrategy()->BinDiscret()->NumMyColElements();
  for (int ibin = 0; ibin < numrowbin; ++ibin)
  {
    DRT::Element* actele = BinStrategy()->BinDiscret()->lColElement(ibin);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele)->DeleteNodes();
  }

  // 2nd) initial particles need to be removed from bindis_
  BinStrategy()->BinDiscret()->DeleteNodes();
}

/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::ParticleHandler::DistributeBinsToProcs()
{
  // create an initial equal distribution of row bins
  Teuchos::RCP<Epetra_Map> rowbins =
      binstrategy_->CreateLinearMapForNumbin(binstrategy_->BinDiscret()->Comm());

  const int maxband = 26;
  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *rowbins, maxband, false));

  // fill all local entries into the graph
  {
    for (int lid = 0; lid < rowbins->NumMyElements(); ++lid)
    {
      const int binId = rowbins->GID(lid);

      std::vector<int> neighbors;
      binstrategy_->GetNeighborBinIds(binId, neighbors);

      int err = graph->InsertGlobalIndices(binId, (int)neighbors.size(), &neighbors[0]);
      if (err < 0)
        dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err, binId);
    }
  }

  // finish graph
  graph->FillComplete();
  graph->OptimizeStorage();

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "PARTITION");

  Epetra_CrsGraph* balanced_graph = NULL;
  try
  {
    balanced_graph = Isorropia::Epetra::createBalancedCopy(*graph, paramlist);
  }
  catch (std::exception& exc)
  {
    std::cout << "Isorropia::createBalancedCopy threw "
              << "exception '" << exc.what() << "' on proc " << myrank_ << std::endl;
    dserror("Error within Isorropia (graph balancing)");
  }

  // obtain the row map
  Teuchos::RCP<Epetra_CrsGraph> rcp_balanced_graph = Teuchos::rcp(balanced_graph);
  rcp_balanced_graph->FillComplete();
  rcp_balanced_graph->OptimizeStorage();
  rowbins = Teuchos::rcp(new Epetra_Map(-1, rcp_balanced_graph->RowMap().NumMyElements(),
      rcp_balanced_graph->RowMap().MyGlobalElements(), 0, binstrategy_->BinDiscret()->Comm()));

  // fill bins into discret
  BinStrategy()->FillBinsIntoBinDiscretization(rowbins);

  // return binrowmap
  return rowbins;
}

/*----------------------------------------------------------------------*
| dynamic load balancing for bin distribution               ghamm 08/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_CrsGraph> PARTICLE::ParticleHandler::CreateGraph()
{
  const Epetra_Map* oldrowmap = binstrategy_->BinDiscret()->ElementRowMap();

  const int maxband = 26;
  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *oldrowmap, maxband, false));

  // fill all local entries into the graph
  {
    for (int lid = 0; lid < oldrowmap->NumMyElements(); ++lid)
    {
      const int binId = oldrowmap->GID(lid);

      std::vector<int> neighbors;
      binstrategy_->GetNeighborBinIds(binId, neighbors);

      int err = graph->InsertGlobalIndices(binId, (int)neighbors.size(), &neighbors[0]);
      if (err < 0)
        dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d", err, binId);
    }
  }

  return graph;
}

/*----------------------------------------------------------------------*
| fill particles into their correct bin on according proc   ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleHandler::FillParticlesIntoBinsRoundRobin(
    std::list<Teuchos::RCP<DRT::Node>>& homelessparticles)
{
  const int numproc = binstrategy_->BinDiscret()->Comm().NumProc();
  const int myrank = binstrategy_->BinDiscret()->Comm().MyPID();  // me
  const int torank = (myrank + 1) % numproc;                      // to
  const int fromrank = (myrank + numproc - 1) % numproc;          // from

  DRT::Exporter exporter(binstrategy_->BinDiscret()->Comm());

  for (int irobin = 0; irobin < numproc; ++irobin)
  {
    std::vector<char> sdata;
    std::vector<char> rdata;

    // ---- pack data for sending -----
    {
      DRT::PackBuffer data;
      for (std::list<Teuchos::RCP<DRT::Node>>::const_iterator currparticle =
               homelessparticles.begin();
           currparticle != homelessparticles.end(); ++currparticle)
      {
        //        cout << " Id:" << (*currparticle)->Id() << " was packed on proc: " << myrank_ <<
        //        endl;
        (*currparticle)->Pack(data);
      }
      data.StartPacking();
      for (std::list<Teuchos::RCP<DRT::Node>>::const_iterator currparticle =
               homelessparticles.begin();
           currparticle != homelessparticles.end(); ++currparticle)
      {
        (*currparticle)->Pack(data);
        binstrategy_->BinDiscret()->DeleteNode((*currparticle)->Id());
      }
      std::swap(sdata, data());
    }


    // ---- send ----
    MPI_Request request;
    exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);


    // ---- receive ----
    int length = rdata.size();
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234 or from != fromrank)
      dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank,
          from, myrank);


    // ---- unpack ----
    {
      // Put received nodes either into discretization or into list of homeless particles
      homelessparticles.clear();
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::vector<char> data;
        DRT::ParObject::ExtractfromPack(index, rdata, data);
        // this Teuchos::rcp holds the memory of the node
        Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
        Teuchos::RCP<DRT::Node> node = Teuchos::rcp_dynamic_cast<DRT::Node>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        // process received particle
        const double* currpos = node->X();
        PlaceNodeCorrectly(node, currpos, homelessparticles);
      }
    }


    // wait for all communication to finish
    exporter.Wait(request);
    binstrategy_->BinDiscret()->Comm().Barrier();  // I feel better this way ;-)
  }                                                // end for irobin

  if (homelessparticles.size())
  {
    std::cout << " There are " << homelessparticles.size()
              << " particles which have left the computational domain on rank " << myrank
              << std::endl;
    // erase everything that is left
    homelessparticles.clear();
  }

  return;
}

/*----------------------------------------------------------------------*
| fill particles into their correct bin on according proc   ghamm 03/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>> PARTICLE::ParticleHandler::FillParticlesIntoBinsRemoteIdList(
    std::list<Teuchos::RCP<DRT::Node>>& homelessparticles)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleHandler::FillParticlesIntoBinsRemoteIdList");
  const int numproc = binstrategy_->BinDiscret()->Comm().NumProc();
  Teuchos::RCP<std::list<int>> removedparticles = Teuchos::rcp(new std::list<int>(0));

  // parallel case
  // ---- find new host procs for particles -----
  const int fullsize = (int)homelessparticles.size();
  std::vector<int> targetbinIdlist;
  targetbinIdlist.reserve(fullsize);
  std::list<Teuchos::RCP<DRT::Node>>::const_iterator hlp;
  for (hlp = homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
  {
    const int binId = binstrategy_->ConvertPosToGid((*hlp)->X());
    targetbinIdlist.push_back(binId);
  }

  // get proc which will be the future host of homeless particles
  std::vector<int> pidlist(fullsize);
  {
    // only unique id lists are accepted in RemoteIDList
    // 1) make gid list unique
    std::set<int> unique_targetbinIdlist(targetbinIdlist.begin(), targetbinIdlist.end());
    std::vector<int> uniquevec_targetbinIdlist(
        unique_targetbinIdlist.begin(), unique_targetbinIdlist.end());
    const int uniquesize = (int)unique_targetbinIdlist.size();

    // 2) communication
    std::vector<int> unique_pidlist(uniquesize);
    int err = binstrategy_->BinDiscret()->ElementRowMap()->RemoteIDList(
        uniquesize, uniquevec_targetbinIdlist.data(), unique_pidlist.data(), NULL);
    if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d", err);

    // 3) build full pid list via lookup table
    std::map<int, int> lookuptable;
    for (int s = 0; s < uniquesize; ++s)
      lookuptable.insert(
          lookuptable.end(), std::pair<int, int>(uniquevec_targetbinIdlist[s], unique_pidlist[s]));
    for (int s = 0; s < fullsize; ++s) pidlist[s] = lookuptable[targetbinIdlist[s]];
  }

  // ---- pack data for sending -----
  std::map<int, std::vector<char>> sdata;
  std::vector<int> targetprocs(numproc, 0);
  int iter = 0;
  for (hlp = homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
  {
    Teuchos::RCP<DRT::Node> iterhomelessparticle = *hlp;

    // ---- pack data for sending -----
    const int targetproc = pidlist[iter];
    if (targetproc != -1)
    {
      DRT::PackBuffer data;
      iterhomelessparticle->Pack(data);
      data.StartPacking();
      iterhomelessparticle->Pack(data);
      binstrategy_->BinDiscret()->DeleteNode(iterhomelessparticle->Id());
      sdata[targetproc].insert(sdata[targetproc].end(), data().begin(), data().end());
      targetprocs[targetproc] = 1;
    }
    else
    {
      const int removeid = iterhomelessparticle->Id();
      removedparticles->push_back(removeid);
      binstrategy_->BinDiscret()->DeleteNode(removeid);
    }
    ++iter;
  }
  if (removedparticles->size() != 0)
    std::cout << " There are " << removedparticles->size()
              << " particles which have left the computational domain on rank " << myrank_
              << std::endl;
  homelessparticles.clear();

  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  binstrategy_->BinDiscret()->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // ---- send ----
  DRT::Exporter exporter(binstrategy_->BinDiscret()->Comm());
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end(); ++p)
  {
    exporter.ISend(
        myrank_, p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  ReceiveParticlesAndFillThemInBins(summedtargets[myrank_], exporter, homelessparticles);

  // wait for all communications to finish
  {
    for (int i = 0; i < length; ++i) exporter.Wait(request[i]);
  }

  binstrategy_->BinDiscret()->Comm().Barrier();  // I feel better this way ;-)

  return removedparticles;
}

/*-----------------------------------------------------------------------------------------*
| fill particles into their correct bin on according proc using ghosting   eichinger 02/17 |
 *-----------------------------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>> PARTICLE::ParticleHandler::FillParticlesIntoBinsUsingGhosting(
    std::list<Teuchos::RCP<DRT::Node>>& homelessparticles)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleHandler::FillParticlesIntoBinsUsingGhosting");

  const int numproc = binstrategy_->BinDiscret()->Comm().NumProc();
  Teuchos::RCP<std::list<int>> removedparticles = Teuchos::rcp(new std::list<int>(0));

  // parallel case
  // ---- find new host procs for particles -----
  std::list<Teuchos::RCP<DRT::Node>>::const_iterator hlp;
  std::map<int, std::list<Teuchos::RCP<DRT::Node>>> towhomwhat;
  int binId, binowner;
  for (hlp = homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
  {
    binId = binstrategy_->ConvertPosToGid((*hlp)->X());
    if (binId == -1)
    {
      binowner = -1;
    }
    else
    {
#ifdef DEBUG
      // safety check
      if (not binstrategy_->BinDiscret()->HaveGlobalElement(binId))
        dserror(
            "To transfer particles using ghosting you need to provide a one layer ghosting,"
            " that is not the case. Bin with gid %i not ghosted on rank %i ",
            binId, myrank_);
#endif
      binowner = binstrategy_->BinDiscret()->gElement(binId)->Owner();
    }
    towhomwhat[binowner].push_back((*hlp));
  }

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // ---- pack data for sending -----
  std::map<int, std::vector<char>> sdata;
  std::vector<int> targetprocs(numproc, 0);
  std::map<int, std::list<Teuchos::RCP<DRT::Node>>>::const_iterator p;
  for (p = towhomwhat.begin(); p != towhomwhat.end(); ++p)
  {
    if (p->first != -1)
    {
      std::list<Teuchos::RCP<DRT::Node>>::const_iterator iter;
      for (iter = p->second.begin(); iter != p->second.end(); ++iter)
      {
        DRT::PackBuffer data;
        (*iter)->Pack(data);
        data.StartPacking();
        (*iter)->Pack(data);
        binstrategy_->BinDiscret()->DeleteNode((*iter)->Id());
        sdata[p->first].insert(sdata[p->first].end(), data().begin(), data().end());
      }
      targetprocs[p->first] = 1;
    }
    else
    {
      std::list<Teuchos::RCP<DRT::Node>>::const_iterator iter;
      for (iter = p->second.begin(); iter != p->second.end(); ++iter)
      {
        const int removeid = (*iter)->Id();
        removedparticles->push_back(removeid);
        binstrategy_->BinDiscret()->DeleteNode(removeid);
      }
    }
  }

  // ---- send ----
  DRT::Exporter exporter(binstrategy_->BinDiscret()->Comm());
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end(); ++p)
  {
    exporter.ISend(
        myrank_, p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");


  if (removedparticles->size() != 0)
    std::cout << "There are " << removedparticles->size()
              << " particles which have "
                 "left the computational domain on rank "
              << myrank_ << std::endl;
  homelessparticles.clear();

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  binstrategy_->BinDiscret()->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // ---- receive -----
  ReceiveParticlesAndFillThemInBins(summedtargets[myrank_], exporter, homelessparticles);

  // wait for all communications to finish
  for (int i = 0; i < length; ++i) exporter.Wait(request[i]);

  // should be no time operation (if we have done everything correctly)
  binstrategy_->BinDiscret()->Comm().Barrier();

  return removedparticles;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void PARTICLE::ParticleHandler::ReceiveParticlesAndFillThemInBins(int const numrec,
    DRT::Exporter& exporter, std::list<Teuchos::RCP<DRT::Node>>& homelessparticles)
{
  // ---- receive ----
  for (int rec = 0; rec < numrec; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234) dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // ---- unpack ----
    {
      // Put received nodes into discretization
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::vector<char> data;
        DRT::ParObject::ExtractfromPack(index, rdata, data);
        // this Teuchos::rcp holds the memory of the node
        Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
        Teuchos::RCP<DRT::Node> node = Teuchos::rcp_dynamic_cast<DRT::Node>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        // process received particle
        const double* currpos = node->X();
        PlaceNodeCorrectly(node, currpos, homelessparticles);
        if (homelessparticles.size())
          dserror(
              "particle (id: %i) was sent to proc %i but corresponding bin (gid: %i) "
              " is missing",
              node->Id(), myrank_, binstrategy_->ConvertPosToGid(currpos));
      }
    }
  }
}

/*----------------------------------------------------------------------*
| node is placed into the correct row bin                   ghamm 09/12 |
 *----------------------------------------------------------------------*/
bool PARTICLE::ParticleHandler::PlaceNodeCorrectly(Teuchos::RCP<DRT::Node> node,
    const double* currpos, std::list<Teuchos::RCP<DRT::Node>>& homelessparticles)
{
  //  std::cout << "on proc: " << myrank_ << " node with ID: " << node->Id() << " and owner: " <<
  //  node->Owner() << " arrived in PlaceNodeCorrectly" << std::endl;
  const int binId = binstrategy_->ConvertPosToGid(currpos);

  // check whether the current node belongs into a bin on this proc
  const bool found = binstrategy_->BinDiscret()->HaveGlobalElement(binId);

  // either fill particle into correct bin on this proc or mark it as homeless
  if (found == true)
  {
    DRT::MESHFREE::MeshfreeMultiBin* currbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(binstrategy_->BinDiscret()->gElement(binId));
#ifdef DEBUG
    if (currbin == NULL)
      dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    // check whether it is a row bin
    if (currbin->Owner() == myrank_)  // row bin
    {
      //      std::cout << "on proc: " << myrank_ << " for node " << node->Id() << " a row bin was
      //      found" << std::endl;
      // node already exists (either row or ghost)
      if (binstrategy_->BinDiscret()->HaveGlobalNode(node->Id()) == true)
      {
        DRT::Node* existingnode = binstrategy_->BinDiscret()->gNode(node->Id());
        // existing node is a row node, this means that node is equal existingnode
        if (existingnode->Owner() == myrank_)
        {
          //          std::cout << "on proc: " << myrank_ << " existingnode row node " <<
          //          existingnode->Id() << " (ID from outside node: " << node->Id() << ") is added
          //          to element: " << currbin->Id() << std::endl;

          // assign node to the correct bin
          currbin->AddNode(existingnode);
        }
        else  // delete existing node, insert received node into discretization and add it to
              // correct bin
        {
          // delete existing node
          binstrategy_->BinDiscret()->DeleteNode(existingnode->Id());
          // update ownership
          node->SetOwner(myrank_);
          // add node
          binstrategy_->BinDiscret()->AddNode(node);
          // assign node to the correct bin
          currbin->AddNode(node.get());

          //        std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to
          //        the discretization and assigned to element: " << currbin->Id() << std::endl;
        }
      }
      else  // fill newly received node into discretization
      {
        // change owner of the node to this proc and add it to the discretization
        node->SetOwner(myrank_);
        binstrategy_->BinDiscret()->AddNode(node);
        //        std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to the
        //        discretization and assigned to element: " << currbin->Id() << std::endl;
        // assign node to the correct bin
        currbin->AddNode(node.get());
      }
      return true;
    }
    else  // ghost bin
    {
      homelessparticles.push_back(node);
      //      std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to
      //      homeless because of becoming a future ghost node" << std::endl;
      return false;
    }
  }
  else  // bin not found on this proc
  {
    homelessparticles.push_back(node);
    //    std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to homeless
    //    because bin is not on this proc " << std::endl;
    return false;
  }
}

/*----------------------------------------------------------------------*
| setup ghosting of bins and particles                      ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleHandler::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap)
{
  // 1st step: ghosting of bins
  {
    // gather bins of rowmap and all its neighbors (row + ghost)
    std::set<int> bins;
    for (int lid = 0; lid < binrowmap->NumMyElements(); ++lid)
    {
      int binId = binrowmap->GID(lid);
      std::vector<int> binvec;
      // get neighboring bins
      binstrategy_->GetNeighborAndOwnBinIds(binId, binvec);
      bins.insert(binvec.begin(), binvec.end());
    }  // end for lid

    // remove non-existing ghost bins from original bin set
    {
      // create copy of column bins
      std::set<int> ghostbins(bins);
      // find ghost bins and check for existence
      for (int lid = 0; lid < binrowmap->NumMyElements(); ++lid)
      {
        const int gid = binrowmap->GID(lid);
        std::set<int>::iterator iter = ghostbins.find(gid);
        if (iter != ghostbins.end()) ghostbins.erase(iter);
      }
      // only ghost bins remain
      std::vector<int> ghostbins_vec(ghostbins.begin(), ghostbins.end());
      const int size = (int)ghostbins.size();
      std::vector<int> pidlist(size);
      const int err = binrowmap->RemoteIDList(size, ghostbins_vec.data(), pidlist.data(), NULL);
      if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d", err);

      for (int i = 0; i < size; ++i)
      {
        if (pidlist[i] == -1)
        {
          std::set<int>::iterator iter = bins.find(ghostbins_vec[i]);
          if (iter == bins.end()) dserror("bin id is missing in bin set");
          // erase non-existing id
          bins.erase(iter);
        }
      }
    }

    // copy bingids to a vector and create bincolmap
    std::vector<int> bincolmap(bins.begin(), bins.end());
    bincolmap_ = Teuchos::rcp(new Epetra_Map(
        -1, (int)bincolmap.size(), &bincolmap[0], 0, binstrategy_->BinDiscret()->Comm()));

    if (bincolmap_->NumGlobalElements() == 1 && bincolmap_->Comm().NumProc() > 1)
      dserror("one bin cannot be run in parallel -> reduce CUTOFF_RADIUS");

    // make sure that all procs are either filled or unfilled
    binstrategy_->BinDiscret()->CheckFilledGlobally();

    // create ghosting for bins (each knowing its particle ids)
    binstrategy_->BinDiscret()->ExtendedGhosting(*bincolmap_, true, false, true, false);
  }


#ifdef DEBUG
  // check whether each proc has only particles that are within bins on this proc
  for (int k = 0; k < BinStrategy()->BinDiscret()->NumMyColElements(); k++)
  {
    int binid = BinStrategy()->BinDiscret()->lColElement(k)->Id();
    DRT::Node** particles = BinStrategy()->BinDiscret()->lColElement(k)->Nodes();

    for (int iparticle = 0; iparticle < BinStrategy()->BinDiscret()->lColElement(k)->NumNode();
         iparticle++)
    {
      int gidofbin = BinStrategy()->ConvertPosToGid(particles[iparticle]->X());
      if (gidofbin != binid)
        dserror("after ghosting: particle which should be in bin no. %i is in %i", gidofbin, binid);
    }
  }
#endif

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>> PARTICLE::ParticleHandler::TransferParticles(
    bool const fill_using_ghosting)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleHandler::TransferParticles");

  // set of homeless particles
  std::list<Teuchos::RCP<DRT::Node>> homelessparticles;

  std::vector<int> examinedbins(binstrategy_->BinDiscret()->NumMyRowElements(), 0);
  // first run over particles and then process whole bin in which particle is located
  // until all particles have been checked
  const int numrownodes = binstrategy_->BinDiscret()->NumMyRowNodes();
  for (int i = 0; i < numrownodes; ++i)
  {
    DRT::Node* currparticle = binstrategy_->BinDiscret()->lRowNode(i);

#ifdef DEBUG
    if (currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
#endif

    DRT::MESHFREE::MeshfreeMultiBin* currbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currparticle->Elements()[0]);
    // as checked above, there is only one element in currele array
    const int binId = currbin->Id();
    const int rlid = binstrategy_->BinDiscret()->ElementRowMap()->LID(binId);

    // if a bin has already been examined --> continue with next particle
    if (examinedbins[rlid]) continue;
    // else: bin is examined for the first time --> new entry in examinedbins_
    else
      examinedbins[rlid] = 1;

    DRT::Node** particles = currbin->Nodes();
    std::vector<int> tobemoved(0);
    for (int iparticle = 0; iparticle < currbin->NumNode(); ++iparticle)
    {
      // get current node
      DRT::Node* currnode = particles[iparticle];

      // transform to array
      std::vector<double> pos(3, 0.0);
      for (int dim = 0; dim < 3; ++dim) pos[dim] = currnode->X()[dim];

      const int gidofbin = binstrategy_->ConvertPosToGid(&pos[0]);
      // particle has left current bin
      if (gidofbin != binId)
      {
        // (looping over nodes and deleting at the same time is detrimental)
        tobemoved.push_back(currnode->Id());
        // find new bin for particle
        PlaceNodeCorrectly(Teuchos::rcp(currnode, false), &pos[0], homelessparticles);
      }
    }

    // finally remove nodes from their old bin
    for (size_t iter = 0; iter < tobemoved.size(); ++iter) currbin->DeleteNode(tobemoved[iter]);
  }

#ifdef DEBUG
  if (homelessparticles.size())
    std::cout << "There are " << homelessparticles.size() << " homeless particles on proc"
              << myrank_ << std::endl;
#endif

  // store particles that have left the computational domain
  Teuchos::RCP<std::list<int>> deletedparticles = Teuchos::rcp(new std::list<int>(0));

  //---------------------------------------------------------------------------
  // numproc == 1
  //---------------------------------------------------------------------------
  if (binstrategy_->BinDiscret()->Comm().NumProc() == 1)
  {
    if (homelessparticles.size())
    {
      std::cout << " There are " << homelessparticles.size()
                << " particles which have left the"
                   " computational domain on rank "
                << myrank_ << std::endl;
      std::list<Teuchos::RCP<DRT::Node>>::const_iterator hlp;
      for (hlp = homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
      {
        const int removeid = (*hlp)->Id();
        deletedparticles->push_back(removeid);
        binstrategy_->BinDiscret()->DeleteNode(removeid);
      }
      homelessparticles.clear();
    }
    return deletedparticles;
  }

  //---------------------------------------------------------------------------
  // numproc > 1
  //---------------------------------------------------------------------------
  // homeless particles are sent to their new processors where they are inserted into their correct
  // bin
  if (fill_using_ghosting)
  {
    deletedparticles = FillParticlesIntoBinsUsingGhosting(homelessparticles);
  }
  else
  {
    deletedparticles = FillParticlesIntoBinsRemoteIdList(homelessparticles);
  }

  return deletedparticles;
}

/*----------------------------------------------------------------------*
| bins are distributed to the processors based on an        ghamm 11/12 |
| underlying discretization                                             |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::ParticleHandler::DistributeBinsToProcsBasedOnUnderlyingDiscret(
    Teuchos::RCP<DRT::Discretization> underlyingdis, std::map<int, std::set<int>>& rowelesinbin)
{
  //--------------------------------------------------------------------
  // 1st step: exploiting bounding box idea for scatra elements and bins
  //--------------------------------------------------------------------

  binstrategy_->DistributeRowElementsToBinsUsingEleXAABB(underlyingdis, rowelesinbin);

  //--------------------------------------------------------------------
  // 2nd step: decide which proc will be owner of each bin
  //--------------------------------------------------------------------

  std::vector<int> rowbins;
  {
    // fill all bins which contain row elements to a vector on each proc
    const int numbin = rowelesinbin.size();
    std::vector<int> mybinswithroweles;
    mybinswithroweles.reserve(numbin);
    for (std::map<int, std::set<int>>::const_iterator it = rowelesinbin.begin();
         it != rowelesinbin.end(); ++it)
      mybinswithroweles.push_back(it->first);

    std::vector<int> maxrank_per_bin;

    // loop over all procs in which proc i determines which bins it will own finally
    for (int iproc = 0; iproc < underlyingdis->Comm().NumProc(); ++iproc)
    {
      // proc i determines number of bins with row elements and broadcasts
      int numbin_proc = 0;
      if (iproc == myrank_) numbin_proc = numbin;
      underlyingdis->Comm().Broadcast(&numbin_proc, 1, iproc);

      // allocate vector on each proc with number of row bins on proc i
      std::vector<int> binswithroweles_proci;

      // proc i broadcasts its bin ids
      if (iproc == myrank_)
        binswithroweles_proci = mybinswithroweles;
      else
        binswithroweles_proci.resize(numbin_proc, 0);
      underlyingdis->Comm().Broadcast(&binswithroweles_proci[0], numbin_proc, iproc);

      // each proc feeds its number of eles for the requested bin ids
      std::vector<int> numeles_per_bin(numbin_proc, 0);
      for (int i = 0; i < numbin_proc; ++i)
      {
        std::map<int, std::set<int>>::const_iterator it =
            rowelesinbin.find(binswithroweles_proci[i]);
        if (it != rowelesinbin.end()) numeles_per_bin[i] = (int)it->second.size();
      }

      // find maximum number of eles in each bin over all procs (init with -1)
      std::vector<int> maxnumeles_per_bin(numbin_proc, -1);
      underlyingdis->Comm().MaxAll(&numeles_per_bin[0], &maxnumeles_per_bin[0], numbin_proc);


#ifdef DEBUG
      // safety check - there is no empty bin allowed
      for (int i = 0; i < numbin_proc; ++i)
      {
        if (maxnumeles_per_bin[i] < 1)
          dserror("empty bin found on proc %d which should not be possible here", myrank_);
      }
#endif

      // it is possible that several procs have the same number of eles in a bin
      // therefore, only proc which has maximum number of eles in a bin writes its rank
      std::vector<int> myrank_per_bin(numbin_proc, -1);
      for (int i = 0; i < numbin_proc; ++i)
      {
        if (numeles_per_bin[i] == maxnumeles_per_bin[i]) myrank_per_bin[i] = myrank_;
      }

      // find maximum myrank for each bin over all procs (init with -1)
      std::vector<int> maxrank_per_bin_tmp(numbin_proc, -1);
      underlyingdis->Comm().MaxAll(&myrank_per_bin[0], &maxrank_per_bin_tmp[0], numbin_proc);

      // store max rank on proc i
      if (iproc == myrank_) maxrank_per_bin = maxrank_per_bin_tmp;
    }

    // assign bins to proc if it has the highest rank
    for (int i = 0; i < numbin; ++i)
    {
      if (myrank_ == maxrank_per_bin[i])
      {
        const int gid = mybinswithroweles[i];
        Teuchos::RCP<DRT::Element> bin =
            DRT::UTILS::Factory("MESHFREEMULTIBIN", "dummy", gid, myrank_);
        binstrategy_->BinDiscret()->AddElement(bin);
        rowbins.push_back(gid);
      }
    }
  }

  // return binrowmap (without having called FillComplete on bindis_ so far)
  return Teuchos::rcp(
      new Epetra_Map(-1, (int)rowbins.size(), &rowbins[0], 0, binstrategy_->BinDiscret()->Comm()));
}

/*-----------------------------------------------------------------------------*
 | build reduced bin col map based on boundary row bins       eichinger 01/17  |
 *-----------------------------------------------------------------------------*/
void PARTICLE::ParticleHandler::GetNeighbouringBinsOfParticleContainingBoundaryRowBins(
    std::set<int>& colbins) const
{
  colbins.clear();

  std::list<DRT::Element*> const boundaryrowbins = binstrategy_->BoundaryRowBins();

  if (boundaryrowbins.size() == 0)
    dserror("Boundary row bins unknown, call function DetermineBoundaryRowBins() first!");

  // loop over boundary row bins and add neighbors of filled row bins
  std::list<DRT::Element*>::const_iterator it;
  for (it = boundaryrowbins.begin(); it != boundaryrowbins.end(); ++it)
  {
    if ((*it)->NumNode() != 0)
    {
      std::vector<int> binvec;
      binvec.reserve(26);
      // get neighboring bins
      binstrategy_->GetNeighborBinIds((*it)->Id(), binvec);
      colbins.insert(binvec.begin(), binvec.end());
    }
  }
}
