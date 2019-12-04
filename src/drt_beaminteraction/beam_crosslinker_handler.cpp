/*----------------------------------------------------------------------*/
/*! \file

\brief Handler to control beam crosslinker simulations

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#include "../drt_binstrategy/drt_meshfree_multibin.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <unordered_set>
#include "beam_crosslinker_handler.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamCrosslinkerHandler::BeamCrosslinkerHandler()
    : binstrategy_(Teuchos::null), myrank_(-1), bincolmap_(Teuchos::null)
{
  // empty constructor
}

void BEAMINTERACTION::BeamCrosslinkerHandler::Init(
    int myrank, Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy)
{
  binstrategy_ = binstrategy;
  myrank_ = myrank;
}

void BEAMINTERACTION::BeamCrosslinkerHandler::Setup()
{
  // so far nothing to do
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamCrosslinkerHandler::DistributeLinkerToBins(
    Teuchos::RCP<Epetra_Map> const& linkerrowmap)
{
  std::list<Teuchos::RCP<DRT::Node>> homelesslinker;
  for (int lid = 0; lid < linkerrowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = binstrategy_->BinDiscret()->gNode(linkerrowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node, false), currpos, homelesslinker);
  }

  // start round robin loop to fill linker into their correct bins
  FillLinkerIntoBinsRoundRobin(homelesslinker);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamCrosslinkerHandler::RemoveAllLinker()
{
  // 1st) loop over bins and remove initial linker info
  const int numrowbin = BinStrategy()->BinDiscret()->NumMyColElements();
  for (int ibin = 0; ibin < numrowbin; ++ibin)
  {
    DRT::Element* actele = BinStrategy()->BinDiscret()->lColElement(ibin);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele)->DeleteNodes();
  }

  // 2nd) initial linker need to be removed from bindis_
  BinStrategy()->BinDiscret()->DeleteNodes();
}

/*----------------------------------------------------------------------*
| fill linker into their correct bin on according proc   ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamCrosslinkerHandler::FillLinkerIntoBinsRoundRobin(
    std::list<Teuchos::RCP<DRT::Node>>& homelesslinker)
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
      for (std::list<Teuchos::RCP<DRT::Node>>::const_iterator currlinker = homelesslinker.begin();
           currlinker != homelesslinker.end(); ++currlinker)
      {
        //        cout << " Id:" << (*currlinker)->Id() << " was packed on proc: " << myrank_ <<
        //        endl;
        (*currlinker)->Pack(data);
      }
      data.StartPacking();
      for (std::list<Teuchos::RCP<DRT::Node>>::const_iterator currlinker = homelesslinker.begin();
           currlinker != homelesslinker.end(); ++currlinker)
      {
        (*currlinker)->Pack(data);
        binstrategy_->BinDiscret()->DeleteNode((*currlinker)->Id());
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
      // Put received nodes either into discretization or into list of homeless linker
      homelesslinker.clear();
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::vector<char> data;
        DRT::ParObject::ExtractfromPack(index, rdata, data);
        // this Teuchos::rcp holds the memory of the node
        Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
        Teuchos::RCP<DRT::Node> node = Teuchos::rcp_dynamic_cast<DRT::Node>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        // process received linker
        const double* currpos = node->X();
        PlaceNodeCorrectly(node, currpos, homelesslinker);
      }
    }


    // wait for all communication to finish
    exporter.Wait(request);
    binstrategy_->BinDiscret()->Comm().Barrier();  // I feel better this way ;-)
  }                                                // end for irobin

  if (homelesslinker.size())
  {
    std::cout << " There are " << homelesslinker.size()
              << " linker which have left the computational domain on rank " << myrank << std::endl;
    // erase everything that is left
    homelesslinker.clear();
  }

  return;
}

/*----------------------------------------------------------------------*
| fill linker into their correct bin on according proc   ghamm 03/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>>
BEAMINTERACTION::BeamCrosslinkerHandler::FillLinkerIntoBinsRemoteIdList(
    std::list<Teuchos::RCP<DRT::Node>>& homelesslinker)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::BeamCrosslinkerHandler::FillLinkerIntoBinsRemoteIdList");
  const int numproc = binstrategy_->BinDiscret()->Comm().NumProc();
  Teuchos::RCP<std::list<int>> removedlinker = Teuchos::rcp(new std::list<int>(0));

  // parallel case
  // ---- find new host procs for linker -----
  const int fullsize = (int)homelesslinker.size();
  std::vector<int> targetbinIdlist;
  targetbinIdlist.reserve(fullsize);
  std::list<Teuchos::RCP<DRT::Node>>::const_iterator hlp;
  for (hlp = homelesslinker.begin(); hlp != homelesslinker.end(); ++hlp)
  {
    const int binId = binstrategy_->ConvertPosToGid((*hlp)->X());
    targetbinIdlist.push_back(binId);
  }

  // get proc which will be the future host of homeless linker
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
  for (hlp = homelesslinker.begin(); hlp != homelesslinker.end(); ++hlp)
  {
    Teuchos::RCP<DRT::Node> iterhomelesslinker = *hlp;

    // ---- pack data for sending -----
    const int targetproc = pidlist[iter];
    if (targetproc != -1)
    {
      DRT::PackBuffer data;
      iterhomelesslinker->Pack(data);
      data.StartPacking();
      iterhomelesslinker->Pack(data);
      binstrategy_->BinDiscret()->DeleteNode(iterhomelesslinker->Id());
      sdata[targetproc].insert(sdata[targetproc].end(), data().begin(), data().end());
      targetprocs[targetproc] = 1;
    }
    else
    {
      const int removeid = iterhomelesslinker->Id();
      removedlinker->push_back(removeid);
      binstrategy_->BinDiscret()->DeleteNode(removeid);
    }
    ++iter;
  }
  if (removedlinker->size() != 0)
    std::cout << " There are " << removedlinker->size()
              << " linker which have left the computational domain on rank " << myrank_
              << std::endl;
  homelesslinker.clear();

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

  ReceiveLinkerAndFillThemInBins(summedtargets[myrank_], exporter, homelesslinker);

  // wait for all communications to finish
  {
    for (int i = 0; i < length; ++i) exporter.Wait(request[i]);
  }

  binstrategy_->BinDiscret()->Comm().Barrier();  // I feel better this way ;-)

  return removedlinker;
}

/*-----------------------------------------------------------------------------------------*
| fill linker into their correct bin on according proc using ghosting   eichinger 02/17 |
 *-----------------------------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>>
BEAMINTERACTION::BeamCrosslinkerHandler::FillLinkerIntoBinsUsingGhosting(
    std::list<Teuchos::RCP<DRT::Node>>& homelesslinker)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::BeamCrosslinkerHandler::FillLinkerIntoBinsUsingGhosting");

  const int numproc = binstrategy_->BinDiscret()->Comm().NumProc();
  Teuchos::RCP<std::list<int>> removedlinker = Teuchos::rcp(new std::list<int>(0));

  // parallel case
  // ---- find new host procs for linker -----
  std::list<Teuchos::RCP<DRT::Node>>::const_iterator hlp;
  std::map<int, std::list<Teuchos::RCP<DRT::Node>>> towhomwhat;
  int binId, binowner;
  for (hlp = homelesslinker.begin(); hlp != homelesslinker.end(); ++hlp)
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
            "To transfer linker using ghosting you need to provide a one layer ghosting,"
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
        removedlinker->push_back(removeid);
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


  if (removedlinker->size() != 0)
    std::cout << "There are " << removedlinker->size()
              << " linker which have "
                 "left the computational domain on rank "
              << myrank_ << std::endl;
  homelesslinker.clear();

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  binstrategy_->BinDiscret()->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // ---- receive -----
  ReceiveLinkerAndFillThemInBins(summedtargets[myrank_], exporter, homelesslinker);

  // wait for all communications to finish
  for (int i = 0; i < length; ++i) exporter.Wait(request[i]);

  // should be no time operation (if we have done everything correctly)
  binstrategy_->BinDiscret()->Comm().Barrier();

  return removedlinker;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamCrosslinkerHandler::ReceiveLinkerAndFillThemInBins(
    int const numrec, DRT::Exporter& exporter, std::list<Teuchos::RCP<DRT::Node>>& homelesslinker)
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

        // process received linker
        const double* currpos = node->X();
        PlaceNodeCorrectly(node, currpos, homelesslinker);
        if (homelesslinker.size())
          dserror(
              "linker (id: %i) was sent to proc %i but corresponding bin (gid: %i) "
              " is missing",
              node->Id(), myrank_, binstrategy_->ConvertPosToGid(currpos));
      }
    }
  }
}

/*----------------------------------------------------------------------*
| node is placed into the correct row bin                   ghamm 09/12 |
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::BeamCrosslinkerHandler::PlaceNodeCorrectly(Teuchos::RCP<DRT::Node> node,
    const double* currpos, std::list<Teuchos::RCP<DRT::Node>>& homelesslinker)
{
  //  std::cout << "on proc: " << myrank_ << " node with ID: " << node->Id() << " and owner: " <<
  //  node->Owner() << " arrived in PlaceNodeCorrectly" << std::endl;
  const int binId = binstrategy_->ConvertPosToGid(currpos);

  // check whether the current node belongs into a bin on this proc
  const bool found = binstrategy_->BinDiscret()->HaveGlobalElement(binId);

  // either fill linker into correct bin on this proc or mark it as homeless
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
      homelesslinker.push_back(node);
      //      std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to
      //      homeless because of becoming a future ghost node" << std::endl;
      return false;
    }
  }
  else  // bin not found on this proc
  {
    homelesslinker.push_back(node);
    //    std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to homeless
    //    because bin is not on this proc " << std::endl;
    return false;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>> BEAMINTERACTION::BeamCrosslinkerHandler::TransferLinker(
    bool const fill_using_ghosting)
{
  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::BeamCrosslinkerHandler::TransferLinker");

  // set of homeless linker
  std::list<Teuchos::RCP<DRT::Node>> homelesslinker;

  std::vector<int> examinedbins(binstrategy_->BinDiscret()->NumMyRowElements(), 0);
  // first run over linker and then process whole bin in which linker is located
  // until all linker have been checked
  const int numrownodes = binstrategy_->BinDiscret()->NumMyRowNodes();
  for (int i = 0; i < numrownodes; ++i)
  {
    DRT::Node* currlinker = binstrategy_->BinDiscret()->lRowNode(i);

#ifdef DEBUG
    if (currlinker->NumElement() != 1) dserror("ERROR: A linker is assigned to more than one bin!");
#endif

    DRT::MESHFREE::MeshfreeMultiBin* currbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currlinker->Elements()[0]);
    // as checked above, there is only one element in currele array
    const int binId = currbin->Id();
    const int rlid = binstrategy_->BinDiscret()->ElementRowMap()->LID(binId);

    // if a bin has already been examined --> continue with next linker
    if (examinedbins[rlid]) continue;
    // else: bin is examined for the first time --> new entry in examinedbins_
    else
      examinedbins[rlid] = 1;

    DRT::Node** linker = currbin->Nodes();
    std::vector<int> tobemoved(0);
    for (int ilinker = 0; ilinker < currbin->NumNode(); ++ilinker)
    {
      // get current node
      DRT::Node* currnode = linker[ilinker];

      // transform to array
      std::vector<double> pos(3, 0.0);
      for (int dim = 0; dim < 3; ++dim) pos[dim] = currnode->X()[dim];

      const int gidofbin = binstrategy_->ConvertPosToGid(&pos[0]);
      // linker has left current bin
      if (gidofbin != binId)
      {
        // (looping over nodes and deleting at the same time is detrimental)
        tobemoved.push_back(currnode->Id());
        // find new bin for linker
        PlaceNodeCorrectly(Teuchos::rcp(currnode, false), &pos[0], homelesslinker);
      }
    }

    // finally remove nodes from their old bin
    for (size_t iter = 0; iter < tobemoved.size(); ++iter) currbin->DeleteNode(tobemoved[iter]);
  }

#ifdef DEBUG
  if (homelesslinker.size())
    std::cout << "There are " << homelesslinker.size() << " homeless linker on proc" << myrank_
              << std::endl;
#endif

  // store linker that have left the computational domain
  Teuchos::RCP<std::list<int>> deletedlinker = Teuchos::rcp(new std::list<int>(0));

  //---------------------------------------------------------------------------
  // numproc == 1
  //---------------------------------------------------------------------------
  if (binstrategy_->BinDiscret()->Comm().NumProc() == 1)
  {
    if (homelesslinker.size())
    {
      std::cout << " There are " << homelesslinker.size()
                << " linker which have left the"
                   " computational domain on rank "
                << myrank_ << std::endl;
      std::list<Teuchos::RCP<DRT::Node>>::const_iterator hlp;
      for (hlp = homelesslinker.begin(); hlp != homelesslinker.end(); ++hlp)
      {
        const int removeid = (*hlp)->Id();
        deletedlinker->push_back(removeid);
        binstrategy_->BinDiscret()->DeleteNode(removeid);
      }
      homelesslinker.clear();
    }
    return deletedlinker;
  }

  //---------------------------------------------------------------------------
  // numproc > 1
  //---------------------------------------------------------------------------
  // homeless linker are sent to their new processors where they are inserted into their correct
  // bin
  if (fill_using_ghosting)
  {
    deletedlinker = FillLinkerIntoBinsUsingGhosting(homelesslinker);
  }
  else
  {
    deletedlinker = FillLinkerIntoBinsRemoteIdList(homelesslinker);
  }

  return deletedlinker;
}

/*-----------------------------------------------------------------------------*
 | build reduced bin col map based on boundary row bins       eichinger 01/17  |
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamCrosslinkerHandler::GetNeighbouringBinsOfLinkerContainingBoundaryRowBins(
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
