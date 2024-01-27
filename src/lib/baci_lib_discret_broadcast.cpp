/*---------------------------------------------------------------------*/
/*! \file

\brief Duplicate and broadcast full discretizations

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_comm_exporter.H"
#include "baci_comm_utils.H"
#include "baci_global_data.H"
#include "baci_global_data_enums.H"
#include "baci_io.H"
#include "baci_lib_discret.H"
#include "baci_nurbs_discret.H"
#include "baci_rebalance.H"

#include <Epetra_FECrsGraph.h>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | broadcast all discretizations from group 0 to all other groups using
 | a ponzi scheme                                          biehler 05/13 |
 *----------------------------------------------------------------------*/
void DRT::BroadcastDiscretizations(GLOBAL::Problem& problem)
{
  // source of all discretizations is always group 0
  int bgroup = 0;
  // group to which discretizations are sent
  int tgroup = -1;

  Teuchos::RCP<CORE::COMM::Communicators> group = problem.GetCommunicators();
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
    numfield = problem.NumFields();
  }
  gcomm->MaxAll(&sbcaster, &bcaster, 1);  // communicate proc null of bgroup
  gcomm->Broadcast(&numfield, 1, bcaster);

  const std::vector<std::string> disnames = problem.GetDisNames();

  for (int i = 0; i < numfield; ++i)
  {
    Teuchos::RCP<DRT::Discretization> dis = Teuchos::null;
    std::string name;
    CORE::FE::ShapeFunctionType distype;
    std::vector<char> data;
    if (gcomm->MyPID() == bcaster)
    {
      DRT::Container cont;
      dis = problem.GetDis(disnames[i]);
      name = dis->Name();
      distype = problem.SpatialApproximationType();
      cont.Add("disname", name);
      cont.Add("distype", (int)distype);
      CORE::COMM::PackBuffer buffer;
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
      CORE::COMM::Exporter exporter(source, target, *gcomm);
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
    distype = (CORE::FE::ShapeFunctionType)cont.GetInt("distype");
    // allocate or get the discretization
    if (group->GroupId() == bgroup)
    {
      dis = problem.GetDis(disnames[i]);
    }
    else
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
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
          NPDuplicateDiscretizationEqualGroupSize(bgroup, tgroup, group, dis, distype, icomm);
        else
          NPDuplicateDiscretization(bgroup, tgroup, group, dis, distype, icomm);
        icomm = Teuchos::null;
        MPI_Comm_free(&intercomm);
      }
      if (group->GroupId() == tgroup)
      {
        const std::string disname = dis->Name();
        problem.AddDis(disname, dis);
      }
      gcomm->Barrier();
    }  // for (int k=numlevels; k>0; --k)
  }    // for (int i=0; i<numfield; ++i)
  return;
}

/*----------------------------------------------------------------------*
 | distribute a discretization from one group to one other    gee 03/12 |
 *----------------------------------------------------------------------*/
void DRT::NPDuplicateDiscretization(const int sgroup, const int rgroup,
    Teuchos::RCP<CORE::COMM::Communicators> group, Teuchos::RCP<DRT::Discretization> dis,
    CORE::FE::ShapeFunctionType distype, Teuchos::RCP<Epetra_MpiComm> icomm)
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
  icomm->MaxAll(sbcaster.data(), bcaster.data(), 2);
  icomm->MaxAll(sgsize.data(), gsize.data(), 2);
  // printf("proc %d sgroup %d rgroup %d bcaster %d %d gsize %d
  // %d\n",icomm->MyPID(),sgroup,rgroup,bcaster[0],bcaster[1],gsize[0],gsize[1]);

  // create a common discretization that we then fill with all stuff from sender group
  std::string name = dis->Name();
  Teuchos::RCP<DRT::Discretization> commondis;
  switch (distype)
  {
    case CORE::FE::ShapeFunctionType::nurbs:
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
    Epetra_Map smap(-1, nummyelements, myelements.data(), 0, *icomm);
    nummyelements = smap.NumGlobalElements();
    myelements.resize(nummyelements);
    for (int i = 0; i < nummyelements; ++i) myelements[i] = i;
    Epetra_Map rmap(-1, nummyelements, myelements.data(), 0, *icomm);
    CORE::COMM::Exporter exporter(smap, rmap, *icomm);
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
    Teuchos::RCP<Epetra_Map> roweles = Teuchos::rcp(
        new Epetra_Map(-1, (int)myrowelements.size(), myrowelements.data(), -1, *icomm));

    Teuchos::RCP<const Epetra_CrsGraph> nodegraph = CORE::REBALANCE::BuildGraph(commondis, roweles);

    Teuchos::ParameterList rebalanceParams;
    rebalanceParams.set<std::string>("num parts", std::to_string(lcomm->NumProc()));

    const auto& [rownodes, colnodes] =
        CORE::REBALANCE::RebalanceNodeMaps(nodegraph, rebalanceParams);

    commondis->Redistribute(*rownodes, *colnodes, false, false, false);
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
      icomm->Broadcast(recvbuff.data(), nrecv, proc);
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
  Epetra_Map targetrowele(-1, (int)targetgids.size(), targetgids.data(), 0, *icomm);

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
      icomm->Broadcast(recvbuff.data(), nrecv, proc);
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
  Epetra_Map targetrownode(-1, (int)targetgids.size(), targetgids.data(), 0, *icomm);

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

    Teuchos::RCP<const Epetra_CrsGraph> nodegraph = CORE::REBALANCE::BuildGraph(dis, roweles);

    Teuchos::ParameterList rebalanceParams;
    rebalanceParams.set<std::string>("num parts", std::to_string(lcomm->NumProc()));

    const auto& [rownodes, colnodes] =
        CORE::REBALANCE::RebalanceNodeMaps(nodegraph, rebalanceParams);

    dis->Redistribute(*rownodes, *colnodes, false, false, false);
    dis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(dis)));
  }
  icomm->Barrier();
}

/*-----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::NPDuplicateDiscretizationEqualGroupSize(const int sgroup, const int rgroup,
    Teuchos::RCP<CORE::COMM::Communicators> group, Teuchos::RCP<DRT::Discretization> dis,
    CORE::FE::ShapeFunctionType distype, Teuchos::RCP<Epetra_MpiComm> icomm)
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
  icomm->MaxAll(sbcaster.data(), bcaster.data(), 2);
  icomm->MaxAll(sgsize.data(), gsize.data(), 2);

  // create a common discretization that we then fill with all stuff from sender group
  std::string name = dis->Name();
  Teuchos::RCP<DRT::Discretization> commondis;
  switch (distype)
  {
    case CORE::FE::ShapeFunctionType::nurbs:
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
    Epetra_Map smap(-1, nummyelements, myelements.data(), 0, *icomm);
    nummyelements = smap.NumGlobalElements();
    myelements.resize(nummyelements);
    for (int i = 0; i < nummyelements; ++i) myelements[i] = i;
    Epetra_Map rmap(-1, nummyelements, myelements.data(), 0, *icomm);
    CORE::COMM::Exporter exporter(smap, rmap, *icomm);
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
        MPI_Send(myglobalelements, lengthSend, MPI_INT,
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
        MPI_Recv(targetgids.data(), lengthRecv, MPI_INT,
            (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag, mpi_icomm, &status);
      }
    }  // for (int proc=0; proc<icomm->NumProc(); ++proc)
  }
  Epetra_Map targetrowele(-1, (int)targetgids.size(), targetgids.data(), 0, *icomm);
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
        MPI_Send(myglobalelements, lengthSend, MPI_INT,
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
        MPI_Recv(targetgids.data(), lengthRecv, MPI_INT,
            (proc + icomm->NumProc() / 2) % icomm->NumProc(), tag, mpi_icomm, &status);
      }
    }  // for (int proc=0; proc<icomm->NumProc(); ++proc)
  }
  Epetra_Map targetrownode(-1, (int)targetgids.size(), targetgids.data(), 0, *icomm);
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

BACI_NAMESPACE_CLOSE
