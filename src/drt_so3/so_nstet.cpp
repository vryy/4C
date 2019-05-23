/*!----------------------------------------------------------------------**##

\brief NStet element

\level 3

\maintainer Christoph Meier
*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

#include "so_nstet.H"
#include "so_surface.H"
#include "so_line.H"

DRT::ELEMENTS::NStetType DRT::ELEMENTS::NStetType::instance_;

DRT::ELEMENTS::NStetType& DRT::ELEMENTS::NStetType::Instance() { return instance_; }


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
DRT::ParObject* DRT::ELEMENTS::NStetType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::NStet* object = new DRT::ELEMENTS::NStet(-1, -1);
  object->Unpack(data);
  return object;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStetType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "NSTET4")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::NStet(id, owner));
    return ele;
  }
  return Teuchos::null;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStetType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::NStet(id, owner));
  return ele;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void DRT::ELEMENTS::NStetType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
  np = 0;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void DRT::ELEMENTS::NStetType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void DRT::ELEMENTS::NStetType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["NSTET4"];

  defs["TET4"].AddIntVector("TET4", 4).AddNamedInt("MAT").AddNamedString("KINEM");
}


/*-----------------------------------------------------------------------
 |  ctor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::NStet(int id, int owner)
    : DRT::Element(id, owner), material_(0), V_(-1.0), nxyz_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::NStet(const DRT::ELEMENTS::NStet& old)
    : DRT::Element(old), material_(old.material_), V_(old.V_), nxyz_(old.nxyz_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::~NStet() { return; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // material_
  AddtoPack(data, material_);
  // stresstype_
  AddtoPack(data, stresstype_);
  // V_
  AddtoPack(data, V_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // material_
  ExtractfromPack(position, data, material_);
  // stresstype_
  stresstype_ = static_cast<StressType>(ExtractInt(position, data));
  // V_
  ExtractfromPack(position, data, V_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::so_nstet_expol(
    LINALG::Matrix<1, 6>& stresses, LINALG::Matrix<4, 6>& nodalstresses)
{
  LINALG::Matrix<4, 1> expol;
  expol(0, 0) = 1.0;
  expol(1, 0) = 1.0;
  expol(2, 0) = 1.0;
  expol(3, 0) = 1.0;
  nodalstresses.Multiply(expol, stresses);
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Print(std::ostream& os) const
{
  os << "NStet ";
  Element::Print(os);
  return;
}


/*====================================================================*/
/* 4-node tetrahedra node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (ksi1, ksi2, ksi3, ksi4) of nodes
 * of a common tetrahedron [-1,1]x[-1,1]x[-1,1]
 *  4-node hexahedron: node 0,1,...,3
 *
 * -----------------------
 *- this is the numbering used in GiD & EXODUS!!
 *      3-
 *      |\ ---
 *      |  \    ---
 *      |    \      ---
 *      |      \        -2
 *      |        \       /\
 *      |          \   /   \
 *      |            X      \
 *      |          /   \     \
 *      |        /       \    \
 *      |      /           \   \
 *      |    /               \  \
 *      |  /                   \ \
 *      |/                       \\
 *      0--------------------------1
 */
/*====================================================================*/

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  gee 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NStet::Volumes()
{
  dserror("volume not impl. yet");
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}


/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             gee 05/08|
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NStet::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gee 05/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::NStet::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      DRT::UTILS::buildLines, this);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitElementsandMaps(std::map<int, DRT::ELEMENTS::NStet*>& elecids,
    std::map<int, DRT::Node*>& noderids, const int myrank, const int numproc,
    DRT::Discretization& dis)
{
  const int numele = dis.NumMyColElements();

  std::vector<int> ctmp;
  std::vector<int> rtmp;

  for (int i = 0; i < numele; ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::NStet* actele = dynamic_cast<DRT::ELEMENTS::NStet*>(dis.lColElement(i));
    if (!actele) dserror("cast to NStet* failed");

    // init the element
    actele->InitElement();

    // register element in list of column nstet elements
    elecids[actele->Id()] = actele;
    ctmp.push_back(actele->Id());
    if (actele->Owner() == myrank) rtmp.push_back(actele->Id());

    // compute a map of all row nodes adjacent to a NStet element
    for (int j = 0; j < actele->NumNode(); ++j)
    {
      DRT::Node* node = actele->Nodes()[j];
      if (myrank == node->Owner()) noderids[node->Id()] = node;
    }
  }  // i

  elecmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)ctmp.size(), &ctmp[0], 0, dis.Comm()));
  elermap_ = Teuchos::rcp(new Epetra_Map(-1, (int)rtmp.size(), &rtmp[0], 0, dis.Comm()));

  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitAdjacency(std::map<int, DRT::ELEMENTS::NStet*>& elecids,
    std::map<int, DRT::Node*>& noderids, std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& adjele,
    std::map<int, std::map<int, DRT::Node*>>& adjnode, std::map<int, std::vector<int>>& adjlm,
    DRT::Discretization& dis)
{
  std::map<int, DRT::Node*>::iterator node;
  for (node = noderids.begin(); node != noderids.end(); ++node)
  {
    DRT::Node* nodeL = node->second;
    const int nodeidL = nodeL->Id();

    // list of adjacent elements
    std::vector<DRT::ELEMENTS::NStet*> myadjele(0);
    for (int j = 0; j < nodeL->NumElement(); ++j)
    {
      const int eleid = node->second->Elements()[j]->Id();
      std::map<int, DRT::ELEMENTS::NStet*>::iterator ele = elecids_.find(eleid);
      if (ele == elecids_.end()) continue;
      myadjele.push_back(ele->second);
    }
    adjele[nodeidL] = myadjele;

    // patch of all nodes adjacent to adjacent elements
    std::map<int, DRT::Node*> nodepatch;
    for (int j = 0; j < (int)myadjele.size(); ++j)
    {
      DRT::Node** nodes = myadjele[j]->Nodes();
      for (int k = 0; k < myadjele[j]->NumNode(); ++k) nodepatch[nodes[k]->Id()] = nodes[k];
    }
    adjnode_[nodeidL] = nodepatch;

    // lm and lmowner arrays
    const int numnodepatch = (int)nodepatch.size();
    const int ndofperpatch = numnodepatch * 3;

    // location and ownership vector of nodal patch
    std::vector<int> lm(ndofperpatch);
    std::map<int, DRT::Node*>::iterator pnode;
    int count = 0;
    for (pnode = nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
    {
      const std::vector<int>& dofs = dis.Dof(pnode->second);
      for (int j = 0; j < (int)dofs.size(); ++j) lm[count++] = dofs[j];
    }
    adjlm[nodeidL] = lm;
  }  // for (node=noderids.begin(); node != noderids.end(); ++node)
  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISnode(std::map<int, int>& misnodesmap,
    std::map<int, DRT::Node*>& rnodes, const int myrank, const int numproc,
    DRT::Discretization& dis)
{
  std::map<int, DRT::Node*>::iterator node;

  std::vector<int> misnodes(0);      // indicator which nodes are mis
  std::vector<int> deletednodes(0);  // all nodes that where deleted (for communication)

  for (int proc = 0; proc < numproc; ++proc)
  {
    if (proc == myrank)
    {
      // count how many elements are adjacent to a node and sort this information
      std::vector<std::pair<int, int>> count;
      for (node = rnodes.begin(); node != rnodes.end(); node++)
        count.push_back(
            std::pair<int, int>((int)adjele_[node->second->Id()].size(), node->second->Id()));
      // sort in ascending order
      sort(count.begin(), count.end());
      int n = count.size() - 1;
      while (n > -1)
      {
        if (count[n].first < MIS_MIN_PATCHSIZE)
          break;  // no patches with MIS_MIN_PATCHSIZE elements or less
        const int actnodeid = count[n].second;
        misnodes.push_back(actnodeid);
#ifdef PRINT_NSTET
        std::cout << "Proc " << myrank << " MIS    " << actnodeid << " numele " << count[n].first
                  << std::endl;
#endif
        // delete all neighboring nodes from its patch from map
        std::map<int, DRT::Node*>& nodepatch = adjnode_[actnodeid];
        std::map<int, DRT::Node*>::iterator neighbor;
        for (neighbor = nodepatch.begin(); neighbor != nodepatch.end(); ++neighbor)
        {
          deletednodes.push_back(neighbor->first);
          rnodes.erase(neighbor->first);
        }

        // rebuild the hierarchy of nodes
        count.clear();
        for (node = rnodes.begin(); node != rnodes.end(); node++)
          count.push_back(
              std::pair<int, int>((int)adjele_[node->second->Id()].size(), node->second->Id()));
        sort(count.begin(), count.end());
        n = count.size() - 1;
#if 0
        for (int i=n; i>=0; --i)
          std::cout << "Proc " << myrank << " n " << n << " still on stack node " << count[i].second << " numele " << count[i].first << std::endl;
#endif
      }
    }  // if (proc==mypid)
    fflush(stdout);


    // broadcast mis nodes of proc
    int size = (int)deletednodes.size();
    dis.Comm().Broadcast(&size, 1, proc);
    if (proc != myrank) deletednodes.resize(size);
    dis.Comm().Broadcast(&deletednodes[0], size, proc);

    // all other procs have to remove nodes adjacent to mis nodes from
    // their potential list
    if (myrank != proc)
    {
      for (int i = 0; i < size; ++i)
      {
        node = rnodes.find(deletednodes[i]);
        if (node != rnodes.end()) rnodes.erase(node);
      }
    }

    dis.Comm().Barrier();
    deletednodes.clear();
  }  // for (proc=0; proc<numproc; ++proc)

  // convert the misnodes vector to a map because its easier to search
  for (int i = 0; i < (int)misnodes.size(); ++i) misnodesmap[misnodes[i]] = misnodes[i];
  misnodes.clear();

#ifdef PRINT_NSTET
  // look for left over nodes in rnodes just to know
  for (node = rnodes.begin(); node != rnodes.end(); node++)
    std::cout << "Proc " << myrank << " Not MIS and NOT ADJ " << *(node->second) << std::endl;
#endif

  fflush(stdout);
  dis.Comm().Barrier();
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISpatchesGreedyI(std::map<int, int>& misnodesmap,
    std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
    std::map<int, std::vector<int>>& pstab_cid_mis,
    std::map<int, std::vector<double>>& pstab_cid_mis_weight,
    std::map<int, DRT::ELEMENTS::NStet*>& elecids,
    std::map<int, DRT::ELEMENTS::NStet*>& elecids_full, std::map<int, DRT::Node*>& noderids,
    const int myrank, const int numproc, DRT::Discretization& dis)
{
  //----------------------------------------------------------------------
  // assign mis nodes all surrounding elements (greedy phase 1)
  for (int proc = 0; proc < numproc; ++proc)
  {
    std::vector<int> sendeles(0);
    std::vector<int> sendelemis(0);           // gid of mis node belonging to that element
    std::vector<double> sendelemisweight(0);  // weight of mis node belonging to that element
    std::map<int, int>::iterator mis;
    if (proc == myrank)
    {
      for (mis = misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        std::vector<DRT::ELEMENTS::NStet*> eles(0);
        DRT::Node* misnode = noderids[mis->first];
#ifdef PRINT_NSTET
        std::cout << "Proc " << myrank << " MIS node " << mis->first << " patchsize "
                  << misnode->NumElement() << std::endl;
#endif
        for (int i = 0; i < misnode->NumElement(); ++i)
        {
          if (elecids.find(misnode->Elements()[i]->Id()) == elecids.end()) continue;
          eles.push_back((DRT::ELEMENTS::NStet*)misnode->Elements()[i]);
          elecids.erase(misnode->Elements()[i]->Id());
        }
        pstab_adjele[mis->first] = eles;
      }


      // make communication vector of all elements taken
      // make map cele -> MIS node
      for (mis = misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        std::vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];
        for (int i = 0; i < (int)adjele.size(); ++i)
        {
          pstab_cid_mis[adjele[i]->Id()].push_back(mis->first);
          pstab_cid_mis_weight[adjele[i]->Id()].push_back(1.0);
          sendeles.push_back(adjele[i]->Id());
          sendelemis.push_back(mis->first);
          sendelemisweight.push_back(1.0);
        }
      }

    }  // if (proc==myrank)

    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size, 1, proc);
    if (proc != myrank)
    {
      sendeles.resize(size);
      sendelemis.resize(size);
      sendelemisweight.resize(size);
    }
    dis.Comm().Broadcast(&sendeles[0], size, proc);
    dis.Comm().Broadcast(&sendelemis[0], size, proc);
    dis.Comm().Broadcast(&sendelemisweight[0], size, proc);

    // all other procs remove the already taken elements from their list
    if (myrank != proc)
    {
      // delete already taken elements from my column element map
      // look whether I have any of the communicated elements in my column map
      // If so, put the corresponding MIS node in my map
      std::map<int, DRT::ELEMENTS::NStet*>::iterator fool;
      for (int i = 0; i < size; ++i)
      {
        elecids.erase(sendeles[i]);
        fool = elecids_full.find(sendeles[i]);
        if (fool == elecids_full.end()) continue;
        pstab_cid_mis[sendeles[i]].push_back(sendelemis[i]);
        pstab_cid_mis_weight[sendeles[i]].push_back(sendelemisweight[i]);
      }

    }  // if (myrank != proc)

    sendeles.clear();
    sendelemis.clear();
    dis.Comm().Barrier();
  }  // for (int proc=0; proc<numproc; ++proc)
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISpatchesGreedyII(std::map<int, int>& misnodesmap,
    std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
    std::map<int, std::vector<double>>& pstab_adjele_weight,
    std::map<int, std::vector<int>>& pstab_cid_mis,
    std::map<int, std::vector<double>>& pstab_cid_mis_weight,
    std::map<int, DRT::ELEMENTS::NStet*>& elecids,
    std::map<int, DRT::ELEMENTS::NStet*>& elecids_full, std::map<int, DRT::Node*>& noderids,
    const int myrank, const int numproc, DRT::Discretization& dis)
{
  for (int proc = 0; proc < numproc; ++proc)
  {
    std::vector<int> sendeles(0);
    std::vector<int> sendelemis(0);  // gid of mis node belonging to that element
    std::vector<double> sendelemisweight(0);
    std::map<int, int>::iterator mis;
    if (proc == myrank)
    {
      for (mis = misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        std::vector<DRT::ELEMENTS::NStet*> eles(0);
        std::vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];

        // all elements already contained in this patch have weight 1 on this
        // mis node
        std::vector<double> weights;
        for (int i = 0; i < (int)adjele.size(); ++i) weights.push_back(1.0);
        pstab_adjele_weight[mis->first] = weights;

        // build a nodal patch
        std::map<int, DRT::Node*> nodepatch;
        for (int i = 0; i < (int)adjele.size(); ++i)
          for (int j = 0; j < adjele[i]->NumNode(); ++j)
            nodepatch[adjele[i]->Nodes()[j]->Id()] = adjele[i]->Nodes()[j];

        std::map<int, DRT::Node*>::iterator fool;
        for (fool = nodepatch.begin(); fool != nodepatch.end(); ++fool)
        {
          DRT::Node* node = fool->second;
          for (int k = 0; k < node->NumElement(); ++k)
          {
            // a candidate to be added to this patch
            DRT::Element* ele = node->Elements()[k];

            // check whether element already taken
            if (elecids.find(ele->Id()) == elecids.end()) continue;

            // check whether the element shares at least three nodes with the patch.
            // (because if it shares 3 out of 4 nodes, it shares a face with the patch)
            int numshare = 0;
            for (int l = 0; l < ele->NumNode(); ++l)
              if (nodepatch.find(ele->Nodes()[l]->Id()) != nodepatch.end()) numshare++;
            if (numshare < 3) continue;  // do not take elements that do not share a face

            // yes, we add this element to the patch
            pstab_adjele[mis->first].push_back((DRT::ELEMENTS::NStet*)ele);
            pstab_adjele_weight[mis->first].push_back(1.0);
            elecids.erase(ele->Id());
            sendeles.push_back(ele->Id());
            sendelemis.push_back(mis->first);
            sendelemisweight.push_back(1.0);
            pstab_cid_mis[ele->Id()].push_back(mis->first);
            pstab_cid_mis_weight[ele->Id()].push_back(1.0);
#ifdef PRINT_NSTET
            std::cout << "Proc " << myrank << " leftover NStet " << ele->Id()
                      << " found full weight MIS node " << mis->second << std::endl;
#endif
          }  // k

        }  // for (fool=nodepatch.begin(); fool != nodepatch.end(); ++fool)

      }  // for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)

    }  // if (proc==myrank)

    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size, 1, proc);
    if (proc != myrank)
    {
      sendeles.resize(size);
      sendelemis.resize(size);
      sendelemisweight.resize(size);
    }
    dis.Comm().Broadcast(&sendeles[0], size, proc);
    dis.Comm().Broadcast(&sendelemis[0], size, proc);
    dis.Comm().Broadcast(&sendelemisweight[0], size, proc);

    // all other procs remove the already taken elements from their list
    // of not-yet taken elements
    if (myrank != proc)
    {
      // delete already taken elements from my column element map
      for (int i = 0; i < size; ++i) elecids.erase(sendeles[i]);

      // look whether I have any of the communicated elements in my column map
      // If so, put the corresponding MIS node in my map
      std::map<int, DRT::ELEMENTS::NStet*>::iterator fool;
      for (int i = 0; i < (int)sendeles.size(); ++i)
      {
        fool = elecids_full.find(sendeles[i]);
        if (fool == elecids_full.end()) continue;
        pstab_cid_mis[sendeles[i]].push_back(sendelemis[i]);
        pstab_cid_mis_weight[sendeles[i]].push_back(sendelemisweight[i]);
      }
    }  // if (myrank != proc)

    sendeles.clear();
    sendelemis.clear();
    dis.Comm().Barrier();
    fflush(stdout);

  }  // for (int proc=0; proc<numproc; ++proc)
  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISpatchesGreedyIII(std::map<int, int>& misnodesmap,
    std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
    std::map<int, std::vector<double>>& pstab_adjele_weight,
    std::map<int, std::vector<int>>& pstab_cid_mis,
    std::map<int, std::vector<double>>& pstab_cid_mis_weight,
    std::map<int, DRT::ELEMENTS::NStet*>& elecids,
    std::map<int, DRT::ELEMENTS::NStet*>& elecids_full, std::map<int, DRT::Node*>& noderids,
    const int myrank, const int numproc, DRT::Discretization& dis)
{
  // assign patches weights of elements that are already taken
  // and no longer present in elecids
  std::map<int, int>::iterator mis;
  for (mis = misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
  {
    std::vector<DRT::ELEMENTS::NStet*> eles(0);
    std::vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];

    // all elements in this patch have weight 1 on this mis node
    std::vector<double> weights;
    for (int i = 0; i < (int)adjele.size(); ++i) weights.push_back(1.0);
    pstab_adjele_weight[mis->first] = weights;
  }

  fflush(stdout);
  dis.Comm().Barrier();

  // key is leftover elecid
  // std::vector<int> is patches this element is (multiple times) connected to
  std::map<int, std::vector<int>> eleonpatches;
  // key is leftover elecids
  // value is number of patches (mutliple times) interested in this element
  std::map<int, int> elenumpatches;

  std::map<int, DRT::ELEMENTS::NStet*>::iterator fool;
#ifdef PRINT_NSTET
  for (int proc = 0; proc < numproc; ++proc)
  {
    if (proc == myrank)
    {
#endif
      for (fool = elecids.begin(); fool != elecids.end(); ++fool)
      {
        // leftover element and its nodes
        DRT::ELEMENTS::NStet* ele = fool->second;
        const int numnode = ele->NumNode();
        DRT::Node** nodes = ele->Nodes();
#ifdef PRINT_NSTET
        std::cout << "Proc " << myrank << " LEFTOVER " << *ele << " " << std::endl;
        fflush(stdout);
#endif
        std::vector<std::vector<int>> nodeonpatch(numnode);
        std::map<int, double> patchmap;

        // find out on which patches that I have these nodes are
        std::map<int, std::vector<DRT::ELEMENTS::NStet*>>::iterator patches;
        for (patches = pstab_adjele.begin(); patches != pstab_adjele.end(); ++patches)
        {
          std::vector<DRT::ELEMENTS::NStet*>& eles = patches->second;

          // build nodal patch because it does not exist yet here
          std::map<int, DRT::Node*> nodalpatch;
          for (int i = 0; i < (int)eles.size(); ++i)
            for (int j = 0; j < eles[i]->NumNode(); ++j)
              nodalpatch[eles[i]->Nodes()[j]->Id()] = eles[i]->Nodes()[j];

          // check all nodes of the element in question against this nodal patch
          for (int i = 0; i < numnode; ++i)
            if (nodalpatch.find(nodes[i]->Id()) != nodalpatch.end())
              nodeonpatch[i].push_back(patches->first);

        }  // for (patches=pstab_adjele.begin(); patches != pstab_adjele.end(); ++patches)

        for (int i = 0; i < numnode; ++i)
        {
          for (int j = 0; j < (int)nodeonpatch[i].size(); ++j)
            eleonpatches[ele->Id()].push_back(nodeonpatch[i][j]);
          elenumpatches[ele->Id()] = (int)eleonpatches[ele->Id()].size();
        }

      }  // for (fool=elecids.begin(); fool != elecids.end(); ++fool)
#ifdef PRINT_NSTET
    }
    fflush(stdout);
    dis.Comm().Barrier();
  }
#endif

  // do communication of how many patches have interest in the leftover elements
  std::map<int, int> recvelenumpatches;
  std::map<int, int>::iterator recvelenum;
  std::map<int, std::vector<int>> recveleonpatches;
  std::map<int, std::vector<int>>::iterator recveleon;
  for (int proc = 0; proc < numproc; ++proc)
  {
    std::vector<int> sendeles(0);
    std::vector<int> sendnums(0);

    if (proc == myrank)
    {
      std::map<int, std::vector<int>>::iterator ele;
      for (ele = eleonpatches.begin(); ele != eleonpatches.end(); ++ele)
      {
        sendnums.push_back(elenumpatches[ele->first]);
        sendeles.push_back(ele->first);                    // element id
        sendeles.push_back((int)ele->second.size());       // number of patches
        for (int i = 0; i < (int)ele->second.size(); ++i)  // patches
          sendeles.push_back(ele->second[i]);
      }

    }  // if (proc==myrank)

    int sizenum = sendnums.size();
    dis.Comm().Broadcast(&sizenum, 1, proc);
    if (proc != myrank) sendnums.resize(sizenum);
    dis.Comm().Broadcast(&sendnums[0], sizenum, proc);

    int sizeele = sendeles.size();
    dis.Comm().Broadcast(&sizeele, 1, proc);
    if (proc != myrank) sendeles.resize(sizeele);
    dis.Comm().Broadcast(&sendeles[0], sizeele, proc);

    if (proc != myrank)
    {
      int* i = &sendeles[0];
      while (i < &sendeles[sizeele - 1] + 1)
      {
        const int eleid = *i;
        i++;
        const int num = *i;
        i++;
        int* j = i;
        i += num;
        fool = elecids.find(eleid);
        if (fool == elecids.end()) continue;  // I have no interest in this element
        recvelenumpatches[eleid] += num;
        for (int k = 0; k < num; ++k) recveleonpatches[eleid].push_back(j[k]);
      }
    }  // if (proc != myrank)

    dis.Comm().Barrier();
    sendeles.clear();
    sendnums.clear();

  }  // for (int proc=0; proc<numproc; ++proc)

  // now all receiving and sending is done, so add up
  for (recvelenum = recvelenumpatches.begin(); recvelenum != recvelenumpatches.end(); ++recvelenum)
    elenumpatches[recvelenum->first] += recvelenum->second;
  for (recveleon = recveleonpatches.begin(); recveleon != recveleonpatches.end(); ++recveleon)
  {
    const int eleid = recveleon->first;
    std::vector<int>& ps = recveleon->second;
    std::map<int, std::vector<int>>::iterator ele = eleonpatches.find(eleid);
    if (ele == eleonpatches.end())
      eleonpatches[eleid] = ps;
    else
      for (int i = 0; i < (int)ps.size(); ++i) eleonpatches[eleid].push_back(ps[i]);
  }



#ifdef PRINT_NSTET
  for (int proc = 0; proc < numproc; ++proc)  // serialization for printing only
  {
    if (proc == myrank)
    {
#endif
      std::map<int, std::vector<int>>::iterator ele;
      for (ele = eleonpatches.begin(); ele != eleonpatches.end(); ++ele)
      {
        fool = elecids.find(ele->first);
        if (fool == elecids.end()) continue;  // I have no interest in this element

        const int eleid = ele->first;
        const int num = elenumpatches[eleid];
        const int num2 = (int)ele->second.size();
        if (num != num2) dserror("num %d != num2 %d", num, num2);
        if (!num) dserror("Element %d not connected to any patch", eleid);

        const double weight = 1.0 / num;

        std::vector<int>& patches = ele->second;
        std::map<int, double> patchweights;  // key: patch id value: total weight the patch gets
        for (int i = 0; i < (int)patches.size(); ++i) patchweights[patches[i]] = 0;
        for (int i = 0; i < (int)patches.size(); ++i) patchweights[patches[i]] += weight;

        // Add element and its weight to all patches its connected to
        std::map<int, double>::iterator p;
        for (p = patchweights.begin(); p != patchweights.end(); ++p)
        {
          pstab_adjele_weight[p->first].push_back(p->second);
          pstab_cid_mis[eleid].push_back(p->first);
          pstab_cid_mis_weight[eleid].push_back(p->second);
          std::map<int, std::vector<DRT::ELEMENTS::NStet*>>::iterator patch =
              pstab_adjele.find(p->first);
          if (patch == pstab_adjele.end()) continue;  // not a patch I own
          patch->second.push_back(fool->second);
        }
        elecids.erase(fool);

#ifdef PRINT_NSTET
        std::cout << "Proc " << myrank << " NStet " << eleid << " numpatch "
                  << elenumpatches[ele->first] << " on";
        for (p = patchweights.begin(); p != patchweights.end(); ++p)
          printf(" p %3d w %10.5e ", p->first, p->second);
        std::cout << std::endl;
        fflush(stdout);
#endif

      }  // for (ele=eleonpatches.begin(); ele != eleonpatches.end(); ++ele)
#ifdef PRINT_NSTET
    }  // if (proc==myrank)
    fflush(stdout);
    dis.Comm().Barrier();
  }  // for (int proc=0; proc<numproc; ++proc)
#endif

  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISAdjacency(
    std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
    std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& adjele,
    std::map<int, std::map<int, DRT::Node*>>& pstab_adjnode,
    std::map<int, std::map<int, DRT::Node*>>& adjnode, std::map<int, std::vector<int>>& pstab_adjlm,
    std::map<int, std::vector<int>>& adjlm, const int myrank, const int numproc,
    DRT::Discretization& dis)
{
  std::map<int, std::vector<DRT::ELEMENTS::NStet*>>::iterator mis;
  for (mis = pstab_adjele.begin(); mis != pstab_adjele.end(); ++mis)
  {
    int id = mis->first;
#ifdef PRINT_NSTET
    int mispatchsize = (int)mis->second.size();
    int patchsize = (int)adjele[id].size();
    std::cout << "Proc " << myrank << " MIS " << id << " mispatchsize " << mispatchsize
              << " patchsize " << patchsize << std::endl;
#endif
    // adjnode
    std::map<int, DRT::Node*> nodepatch;
    std::vector<DRT::ELEMENTS::NStet*>& ele = mis->second;
    for (int j = 0; j < (int)ele.size(); ++j)
    {
      DRT::Node** nodes = ele[j]->Nodes();
      for (int k = 0; k < ele[j]->NumNode(); ++k) nodepatch[nodes[k]->Id()] = nodes[k];
    }
    pstab_adjnode[id] = nodepatch;

    // lm array
    std::vector<int> lm;
    std::map<int, DRT::Node*>::iterator pnode;
    for (pnode = nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
    {
      const std::vector<int>& dofs = dis.Dof(pnode->second);
      for (int j = 0; j < (int)dofs.size(); ++j) lm.push_back(dofs[j]);
    }
    pstab_adjlm[id] = lm;
  }  // for (mis=pstab_adjele.begin(); mis != pstab_adjele.end(); ++mis)
  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::ELEMENTS::NStetType::InitMISStressMap(
    std::map<int, std::vector<int>>& pstab_cid_mis, DRT::Discretization& dis)
{
  std::map<int, std::vector<int>>::iterator fool;
  std::map<int, int> ngidmap;
  std::vector<int> ngid;
  for (fool = pstab_cid_mis.begin(); fool != pstab_cid_mis.end(); ++fool)
  {
    std::vector<int>& mis = fool->second;
    for (int i = 0; i < (int)mis.size(); ++i) ngidmap[mis[i]] = mis[i];
  }

  std::map<int, int>::iterator fool2;
  for (fool2 = ngidmap.begin(); fool2 != ngidmap.end(); ++fool2) ngid.push_back(fool2->first);

  return (Teuchos::rcp(new Epetra_Map(-1, (int)ngid.size(), &ngid[0], 0, dis.Comm())));
}



/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NStetType::Initialize(DRT::Discretization& dis)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::NStetType::Initialize");

  const int myrank = dis.Comm().MyPID();
  const int numproc = dis.Comm().NumProc();

  //----------------------------------------------------------------------
  // init elements, make maps of column elements and row nodes
  InitElementsandMaps(elecids_, noderids_, myrank, numproc, dis);

  //----------------------------------------------------------------------
  // compute adjacency for each row node
  // make patch of adjacent elements
  // make patch of adjacent nodes (including center node itself)
  // make lm for nodal patch
  InitAdjacency(elecids_, noderids_, adjele_, adjnode_, adjlm_, dis);

  //----------------------------------------------------------------------
  // build parallel maximum independent set of nodes (MIS-nodes)
  // this is done in pseudo-serial, as a true parallel MIS algorithm ispretty difficult
#ifndef PUSOSOLBERG
  std::map<int, int> misnodesmap;
  std::map<int, DRT::Node*> rnodes = noderids_;  // working copy of row nodes
  InitMISnode(misnodesmap, rnodes, myrank, numproc, dis);


  //----------------------------------------------------------------------
  // each MIS node is associated with a patch of column elements of which
  // it takes the full integration area
  // Additionally, there are leftover column elements (not adjacent to any MIS node)
  // they are taken by a greedy algorithm in a second phase

  //----------------------------------------------------------------------
  std::map<int, DRT::ELEMENTS::NStet*> elecids = elecids_;
  InitMISpatchesGreedyI(misnodesmap, pstab_adjele_, pstab_cid_mis_, pstab_cid_mis_weight_, elecids,
      elecids_, noderids_, myrank, numproc, dis);


  //----------------------------------------------------------------------
  // assign leftover elements that are well connected to a patch to that patch
  // this is a distance 2 patch search
  int rn = 1;
  int count = 1;
  while (rn)
  {
#ifdef PRINT_NSTET
    if (!myrank) std::cout << "Greedy II round " << count << std::endl;
    fflush(stdout);
    dis.Comm().Barrier();
#endif

    InitMISpatchesGreedyII(misnodesmap, pstab_adjele_, pstab_adjele_weight_, pstab_cid_mis_,
        pstab_cid_mis_weight_, elecids, elecids_, noderids_, myrank, numproc, dis);

    int sn = (int)elecids.size();
    dis.Comm().SumAll(&sn, &rn, 1);
    if (count == NSTET_NUMGREEDY) break;
    count++;
  }


  //----------------------------------------------------------------------
  // split remaining elements among patches that its nodes belong to
  InitMISpatchesGreedyIII(misnodesmap, pstab_adjele_, pstab_adjele_weight_, pstab_cid_mis_,
      pstab_cid_mis_weight_, elecids, elecids_, noderids_, myrank, numproc, dis);

  //----------------------------------------------------------------------
  // test whether all column elements on all procs have been assigned a patch
  std::map<int, DRT::ELEMENTS::NStet*>::iterator ele;
  if ((int)elecids.size() != 0)
  {
    for (ele = elecids.begin(); ele != elecids.end(); ++ele)
    {
      std::cout << "Proc " << myrank << " leftover NStet with no MIS node " << *ele->second
                << std::endl;
    }
    dserror("Proc %d has the above column elements leftover", myrank);
  }

  // test whether all column elements on this proc know their MIS node_pos
  for (ele = elecids_.begin(); ele != elecids_.end(); ++ele)
  {
    std::map<int, std::vector<int>>::iterator fool = pstab_cid_mis_.find(ele->first);
    if (fool == pstab_cid_mis_.end())
    {
      std::cout << "This element did not find its MIS node:\n" << *ele->second << std::endl;
      dserror("Element %d did not find its MIS node", ele->first);
    }
  }
#endif

  //----------------------------------------------------------------------
  // have to build adjnode and adjlm arrays for the patches
  InitMISAdjacency(
      pstab_adjele_, adjele_, pstab_adjnode_, adjnode_, pstab_adjlm_, adjlm_, myrank, numproc, dis);


  //----------------------------------------------------------------------
  // create an overlapping map that contains stress data of MIS nodes on all procs
  // that will need it for stress output
  pstab_misstressout_ = InitMISStressMap(pstab_cid_mis_, dis);

  return 0;
}
