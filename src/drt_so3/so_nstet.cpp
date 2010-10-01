/*!----------------------------------------------------------------------**##
\file so_nstet.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>
#include "so_nstet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::NStetType DRT::ELEMENTS::NStetType::instance_;


DRT::ParObject* DRT::ELEMENTS::NStetType::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::NStet* object = new DRT::ELEMENTS::NStet(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStetType::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="NSTET4" )
  {
    RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NStet(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::NStetType::Create( const int id, const int owner )
{
  RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::NStet(id,owner));
  return ele;
}


void DRT::ELEMENTS::NStetType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::NStetType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::NStetType::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["NSTET4"];

  defs["TET4"]
    .AddIntVector("TET4",4)
    .AddNamedInt("MAT")
    ;
}


/*-----------------------------------------------------------------------
 |  ctor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::NStet(int id, int owner) :
DRT::Element(id,owner),
material_(0),
V_(-1.0),
nxyz_(),
FisNew_(false),
F_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::NStet(const DRT::ELEMENTS::NStet& old) :
DRT::Element(old),
material_(old.material_),
V_(old.V_),
nxyz_(old.nxyz_),
FisNew_(old.FisNew_),
F_(old.F_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              gee 05/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::NStet::~NStet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // material_
  AddtoPack(data,material_);
  // stresstype_
  AddtoPack(data,stresstype_);
  // V_
  AddtoPack(data,V_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // material_
  ExtractfromPack(position,data,material_);
  // stresstype_
  ExtractfromPack(position,data,stresstype_);
  // V_
  ExtractfromPack(position,data,V_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::so_nstet_expol(LINALG::Matrix<1,6>& stresses,
                                        LINALG::Matrix<4,6>& nodalstresses)
{
  LINALG::Matrix<4,1> expol;
  expol(0,0)=1.0;
  expol(1,0)=1.0;
  expol(2,0)=1.0;
  expol(3,0)=1.0;
  nodalstresses.Multiply(expol,stresses);
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::Print(ostream& os) const
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
vector<RCP<DRT::Element> > DRT::ELEMENTS::NStet::Volumes()
{
  dserror("volume not impl. yet");
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             gee 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NStet::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gee 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::NStet::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitElementsandMaps(
                           map<int,DRT::ELEMENTS::NStet*>& elecids,
                           map<int,DRT::Node*>&            noderids,
                           const int                       myrank,
                           const int                       numproc,
                           DRT::Discretization&            dis)
{
  const int numele = dis.NumMyColElements();
  for (int i=0; i<numele; ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::NStet* actele = 
                    dynamic_cast<DRT::ELEMENTS::NStet*>(dis.lColElement(i));
    if (!actele) dserror("cast to NStet* failed");

    // init the element
    actele->InitElement();

    // register element in list of column nstet elements
    elecids[actele->Id()] = actele;

    // compute a map of all row nodes adjacent to a NStet element
    for (int j=0; j<actele->NumNode(); ++j) 
    {
      DRT::Node* node = actele->Nodes()[j];
      if (myrank == node->Owner())
        noderids[node->Id()] = node;
    }
  }
  
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitAdjacency(
                     map<int,DRT::ELEMENTS::NStet*>&          elecids,
                     map<int,DRT::Node*>&                     noderids,
                     map<int,vector<DRT::ELEMENTS::NStet*> >& adjele,
                     map<int,map<int,DRT::Node*> >&           adjnode,
                     map<int,vector<int> >&                   adjlm,
                     DRT::Discretization&                     dis)
{
  std::map<int,DRT::Node*>::iterator node;
  for (node=noderids.begin(); node != noderids.end(); ++node)
  {
    DRT::Node* nodeL  = node->second;
    const int nodeidL = nodeL->Id();

    // list of adjacent elements
    vector<DRT::ELEMENTS::NStet*> myadjele(0);
    for (int j=0; j<nodeL->NumElement(); ++j)
    {
      const int eleid = node->second->Elements()[j]->Id();
      std::map<int,DRT::ELEMENTS::NStet*>::iterator ele = elecids_.find(eleid);
      if (ele==elecids_.end()) continue;
      myadjele.push_back(ele->second);
    }
    adjele[nodeidL] = myadjele;

    // patch of all nodes adjacent to adjacent elements
    map<int,DRT::Node*> nodepatch;
    for (int j=0; j<(int)myadjele.size(); ++j) 
    {
      DRT::Node** nodes = myadjele[j]->Nodes();
      for (int k=0; k<myadjele[j]->NumNode(); ++k)
        nodepatch[nodes[k]->Id()] = nodes[k];
    }
    adjnode_[nodeidL] = nodepatch;

    // lm and lmowner arrays
    const int numnodepatch = (int)nodepatch.size();
    const int ndofperpatch = numnodepatch*3;

    // location and ownership vector of nodal patch
    vector<int> lm(ndofperpatch);
    std::map<int,DRT::Node*>::iterator pnode;
    int count=0;
    for (pnode=nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
    {
      const vector<int>& dofs = dis.Dof(pnode->second);
      for (int j=0; j<(int)dofs.size(); ++j)
        lm[count++]        = dofs[j];
    }
    adjlm[nodeidL] = lm;
  } // for (node=noderids.begin(); node != noderids.end(); ++node)
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISnode(
                   map<int,int>&                   misnodesmap,
                   map<int,DRT::Node*>&            rnodes,
                   const int                       myrank,
                   const int                       numproc,
                   DRT::Discretization&            dis)
{
  map<int,DRT::Node*>::iterator node;
  
  vector<int> misnodes(0); // indicator which nodes are mis
  vector<int> smisnodes(0);// same but for communication

  for (int proc=0; proc<numproc; ++proc)
  {
    if (proc==myrank)
    {
      for (node=rnodes.begin(); node != rnodes.end(); )
      {
        const int actnodeid = node->second->Id();
        // make this node a mis node 
        misnodes.push_back(actnodeid);
        
        cout << myrank << " MIS    " << *(node->second) << endl;
        
        // delete all neighboring nodes from its patch from map
        map<int,DRT::Node*>& nodepatch = adjnode_[actnodeid];
        std::map<int,DRT::Node*>::iterator neighbor;
        for (neighbor=nodepatch.begin(); neighbor != nodepatch.end(); ++neighbor)
        {
          if (neighbor->first == actnodeid) continue; // do not delete myself
          rnodes.erase(neighbor->first);
        }
        rnodes.erase(node++);
      }
    } // if (proc==mypid)

    // broadcast mis nodes of proc
    int size = misnodes.size();
    dis.Comm().Broadcast(&size,1,proc);
    if (proc==myrank) smisnodes = misnodes;
    else              smisnodes.resize(size);
    dis.Comm().Broadcast(&smisnodes[0],size,proc);

    // all other procs have to remove nodes adjacent to mis nodes from 
    // their potential list
    if (myrank!=proc)
    {
      for (node=rnodes.begin(); node != rnodes.end();)
      {
        const int actnodeid = node->second->Id();
        map<int,DRT::Node*>& nodepatch = adjnode_[actnodeid];
        std::map<int,DRT::Node*>::iterator neighbor;
        bool foundit = false;
        for (neighbor=nodepatch.begin(); neighbor != nodepatch.end(); ++neighbor)
        {
          const int neighborid = neighbor->first;
          for (int i=0; i<size; ++i)
            if (neighborid == smisnodes[i])
            {
              foundit = true;
              break;
            }
          if (foundit) break;
        }
        if (foundit)
        {
          //cout << myrank << " deloff " << *(node->second) << endl;
          rnodes.erase(node++);
        }
        else node++;
      }
    }
    dis.Comm().Barrier();
    smisnodes.clear();
  } // for (proc=0; proc<numproc; ++proc)

  // convert the misnodes vector to a map because its easier to search
  for (int i=0; i<(int)misnodes.size(); ++i)
    misnodesmap[misnodes[i]] = misnodes[i];
  misnodes.clear();


  // look for left over nodes in rnodes just to know
  for (node=rnodes.begin(); node != rnodes.end(); node++)
    cout << myrank << " Not MIS and NOT ADJ " << *(node->second) << endl;
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISpatchesGreedyI(
                           map<int,int>&                            misnodesmap,
                           map<int,vector<DRT::ELEMENTS::NStet*> >& pstab_adjele,
                           map<int,vector<int> >&                   pstab_cid_mis,
                           map<int,DRT::ELEMENTS::NStet*>&          elecids,
                           map<int,DRT::ELEMENTS::NStet*>&          elecids_full,
                           map<int,DRT::Node*>&                     noderids,
                           const int                                myrank,
                           const int                                numproc,
                           DRT::Discretization&                     dis)
{
  //----------------------------------------------------------------------
  // assign mis nodes all surrounding row elements (greedy phase 1)
  for (int proc=0; proc<numproc; ++proc)
  {
    vector<int> sendeles(0);
    vector<int> sendelemis(0); // gid of mis node belonging to that element
    map<int,int>::iterator mis;
    if (proc==myrank)
    {
      for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        vector<DRT::ELEMENTS::NStet*> eles(0);
        DRT::Node* misnode = noderids[mis->first];
        for (int i=0; i<misnode->NumElement(); ++i)
        {
          if (elecids.find(misnode->Elements()[i]->Id()) == elecids.end()) continue;
          eles.push_back((DRT::ELEMENTS::NStet*)misnode->Elements()[i]);
          elecids.erase(misnode->Elements()[i]->Id());
        }
        pstab_adjele[mis->first] = eles;
      } 
      
      // Unassign mis nodes with patches that are too small
      for (mis=misnodesmap.begin(); mis != misnodesmap.end();)
      {
        vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];
        const int patchsize = (int)adjele.size();
        cout << "Proc " << myrank << " MIS node " << mis->first << " patchsize " << patchsize << endl;
        if (patchsize >= MIS_MIN_PATCHSIZE) mis++;
        else
        {
          cout << "Proc " << myrank << " erased MIS node " << mis->first << " patchsize " << patchsize << endl;
          for (int i=0; i<patchsize; ++i) elecids[adjele[i]->Id()] = adjele[i];
          pstab_adjele.erase(mis->first);
          misnodesmap.erase(mis++);
        }
      } 
      
      // make communication vector of all elements taken
      // make map cele -> MIS node
      for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];
        for (int i=0; i<(int)adjele.size(); ++i) 
        {
          sendeles.push_back(adjele[i]->Id());
          pstab_cid_mis[adjele[i]->Id()].push_back(mis->first);
          sendelemis.push_back(mis->first);
        }
      }
      
    } // if (proc==myrank)
    
    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size,1,proc);
    if (proc != myrank) 
    {
      sendeles.resize(size);
      sendelemis.resize(size);
    }
    dis.Comm().Broadcast(&sendeles[0],size,proc);
    dis.Comm().Broadcast(&sendelemis[0],size,proc);
    
    // all other procs remove the already taken elements from their list
    if (myrank != proc)
    {
      // delete already taken elements from my column element map
      for (int i=0; i<size; ++i) elecids.erase(sendeles[i]);
      
      // look whether I have any of the communicated elements in my column map
      // If so, put the corresponding MIS node in my map
      map<int,DRT::ELEMENTS::NStet*>::iterator fool;
      for (int i=0; i<(int)sendeles.size(); ++i)
      {
        fool = elecids_full.find(sendeles[i]);
        if (fool == elecids_full.end()) continue;
        if (fool->first != sendeles[i]) dserror("gid of element mismatch");
        pstab_cid_mis[sendeles[i]].push_back(sendelemis[i]);
      }
    } // if (myrank != proc)
    
    sendeles.clear();
    sendelemis.clear();
    dis.Comm().Barrier();
  } // for (int proc=0; proc<numproc; ++proc)
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISpatchesGreedyII(
                           map<int,int>&                            misnodesmap,
                           map<int,vector<DRT::ELEMENTS::NStet*> >& pstab_adjele,
                           map<int,vector<double> >&                pstab_adjele_weight,
                           map<int,vector<int> >&                   pstab_cid_mis,
                           map<int,DRT::ELEMENTS::NStet*>&          elecids,
                           map<int,DRT::ELEMENTS::NStet*>&          elecids_full,
                           map<int,DRT::Node*>&                     noderids,
                           const int                                myrank,
                           const int                                numproc,
                           DRT::Discretization&                     dis)
{
  for (int proc=0; proc<numproc; ++proc)
  {
    vector<int> sendeles(0);
    vector<int> sendelemis(0); // gid of mis node belonging to that element
    map<int,int>::iterator mis;
    if (proc==myrank)
    {
      for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
      {
        vector<DRT::ELEMENTS::NStet*> eles(0);
        vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];
        
        // all elements already contained in this patch have weight 1 on this
        // mis node
        vector<double> weights;
        for (int i=0; i<(int)adjele.size(); ++i)
          weights.push_back(1.0);
        pstab_adjele_weight[mis->first] = weights;
        
        // build a nodal patch
        map<int,DRT::Node*> nodepatch;
        for (int i=0; i<(int)adjele.size(); ++i)
          for (int j=0; j<adjele[i]->NumNode(); ++j)
            nodepatch[adjele[i]->Nodes()[j]->Id()] = adjele[i]->Nodes()[j];
        
        map<int,DRT::Node*>::iterator fool;
        for (fool=nodepatch.begin(); fool != nodepatch.end(); ++fool)
        {
          DRT::Node* node = fool->second;
          //if (node->Id() == mis->first) continue;
          for (int k=0; k<node->NumElement(); ++k)
          {
            // a candidate to be added to this patch
            DRT::Element* ele = node->Elements()[k];
            
            // check whether element already taken
            if (elecids.find(ele->Id()) == elecids.end()) continue;
            
            // check whether this element is completely off-processor
            // that is, none of its nodes belongs to me
            bool iownanode = false;
            for (int l=0; l<ele->NumNode(); ++l)
              if (dis.NodeRowMap()->LID(ele->Nodes()[l]->Id())!=-1)
              {
                iownanode = true;
                break;
              }
            if (!iownanode) 
            {
              // cout << "Proc " << myrank << " Yes, it does happen 1\n"; fflush(stdout);
              continue; // don't take this element, its off-processor
            }
            
            // check whether the element shares at least three nodes with the patch.
            // (because if it shares 3 out of 4 nodes, it shares a face with the patch)
            int numshare = 0;
            for (int l=0; l<ele->NumNode(); ++l)
              if (nodepatch.find(ele->Nodes()[l]->Id()) != nodepatch.end())
                numshare++;
            if (numshare<3) continue;
            
            // yes, we add this element to the patch
            pstab_adjele[mis->first].push_back((DRT::ELEMENTS::NStet*)ele);
            pstab_adjele_weight[mis->first].push_back(1.0);
            elecids.erase(ele->Id());
            sendeles.push_back(ele->Id());
            sendelemis.push_back(mis->second);
            pstab_cid_mis[ele->Id()].push_back(mis->first);
            
            cout << "Proc " << myrank << " leftover NStet " << ele->Id() << 
                    " found full weight MIS node " << mis->second << endl;

          } // k

        } // for (fool=nodepatch.begin(); fool != nodepatch.end(); ++fool)
        
      } // for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)

    } // if (proc==myrank)

    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size,1,proc);
    if (proc != myrank) 
    {
      sendeles.resize(size);
      sendelemis.resize(size);
    }
    dis.Comm().Broadcast(&sendeles[0],size,proc);
    dis.Comm().Broadcast(&sendelemis[0],size,proc);

    // all other procs remove the already taken elements from their list
    // of not-yet taken elements
    if (myrank != proc)
    {
      // delete already taken elements from my column element map
      for (int i=0; i<size; ++i) elecids.erase(sendeles[i]);
      
      // look whether I have any of the communicated elements in my column map
      // If so, put the corresponding MIS node in my map
      map<int,DRT::ELEMENTS::NStet*>::iterator fool;
      for (int i=0; i<(int)sendeles.size(); ++i)
      {
        fool = elecids_.find(sendeles[i]);
        if (fool == elecids_full.end()) continue;
        pstab_cid_mis[sendeles[i]].push_back(sendelemis[i]);
      }
    } // if (myrank != proc)
    
    sendeles.clear();
    sendelemis.clear();
    dis.Comm().Barrier();

  } // for (int proc=0; proc<numproc; ++proc)
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISpatchesGreedyIII(
                           map<int,int>&                            misnodesmap,
                           map<int,vector<DRT::ELEMENTS::NStet*> >& pstab_adjele,
                           map<int,vector<double> >&                pstab_adjele_weight,
                           map<int,vector<int> >&                   pstab_cid_mis,
                           map<int,DRT::ELEMENTS::NStet*>&          elecids,
                           map<int,DRT::ELEMENTS::NStet*>&          elecids_full,
                           map<int,DRT::Node*>&                     noderids,
                           const int                                myrank,
                           const int                                numproc,
                           DRT::Discretization&                     dis)
{
  for (int proc=0; proc<numproc; ++proc)
  {
    vector<int> sendeles(0);
    vector<int> sendelemis(0); // gid of mis node belonging to that element
    map<int,DRT::ELEMENTS::NStet*>::iterator fool;

    // assign patches weights of elements that are already taken 
    // and no longer present in elecids
    map<int,int>::iterator mis;
    for (mis=misnodesmap.begin(); mis != misnodesmap.end(); ++mis)
    {
      vector<DRT::ELEMENTS::NStet*> eles(0);
      vector<DRT::ELEMENTS::NStet*>& adjele = pstab_adjele[mis->first];
      
      // all elements in this patch have weight 1 on this mis node
      vector<double> weights;
      for (int i=0; i<(int)adjele.size(); ++i)
        weights.push_back(1.0);
      pstab_adjele_weight[mis->first] = weights;
    }

    if (proc==myrank)
    {
      for (fool=elecids.begin(); fool != elecids.end();)
      {
        // leftover element and its nodes
        DRT::ELEMENTS::NStet* ele = fool->second;
        const int numnode = ele->NumNode();
        DRT::Node** nodes = ele->Nodes();
        
        cout << "Proc " << myrank <<
                " leftover NStet with no MIS node " << *ele << " " << endl; fflush(stdout);
        
        vector<vector<int> > nodeonpatch(numnode);
        
        // find out on which patch these nodes are
        map<int,vector<DRT::ELEMENTS::NStet*> >::iterator patches;
        for (patches=pstab_adjele.begin(); patches != pstab_adjele.end(); ++patches)
        {
          vector<DRT::ELEMENTS::NStet*>& eles = patches->second;

          // build nodal patch
          map<int,DRT::Node*> nodalpatch;
          for (int i=0; i<(int)eles.size(); ++i)
            for (int j=0; j<eles[i]->NumNode(); ++j)
              nodalpatch[eles[i]->Nodes()[j]->Id()] = eles[i]->Nodes()[j];
          
          // check all row nodes against this nodal patch
          for (int i=0; i<numnode; ++i)
          {
            // only row nodes get weights
            if (nodes[i]->Owner() != myrank) continue;
            if (nodalpatch.find(nodes[i]->Id()) != nodalpatch.end())
              nodeonpatch[i].push_back(patches->first);
          }
        } // patches
        
#if 1
        for (int i=0; i<numnode; ++i)
        {
          cout << "Proc " << myrank << " Node " << nodes[i]->Id() << " on ";
          if ((int)nodeonpatch[i].size() == 0) cout << "NO PATCH!";
          for (int j=0; j<(int)nodeonpatch[i].size(); ++j)
            cout << nodeonpatch[i][j] << " ";
          cout << endl;
        }
#endif
        
        // build a map of the patches the nodes are on
        // count the splitting portions
        int count=0;
        map<int,double> patchmap;
        for (int i=0; i<numnode; ++i)
        {
          for (int j=0; j<(int)nodeonpatch[i].size(); ++j)
          {
            patches = pstab_adjele.find(nodeonpatch[i][j]);
            if (patches==pstab_adjele.end()) dserror("???");
            patchmap[nodeonpatch[i][j]] = 0.0;
            count++;
          }
        }
        if (!count) 
        {
          fool++;
          continue; // this element has no business on this proc
        }
        
        // distribute the weights among the patches
        for (int i=0; i<numnode; ++i)
        {
          // each node gets 1/4 of the element
          // this 1/4 is split among the present patches of that node
          double weight = 0.25;
          if (nodeonpatch[i].size()) weight /= (double)nodeonpatch[i].size();
          else weight = 0.0;
          for (int j=0; j<(int)nodeonpatch[i].size(); ++j)
            patchmap[nodeonpatch[i][j]] += weight;
        }

        map<int,double>::iterator p;
#if 1
        for (p=patchmap.begin(); p != patchmap.end(); ++p)
          cout << "Proc " << myrank << " patch " << p->first << " weight " << p->second << endl; fflush(stdout);
#endif

        // go through the patch map and add ele to patches
        // Also add the weight
        for (p=patchmap.begin(); p != patchmap.end(); ++p)
        {
          pstab_adjele[p->first].push_back(ele);
          pstab_adjele_weight[p->first].push_back(p->second);
          pstab_cid_mis[ele->Id()].push_back(p->first);
        }
        
        elecids.erase(fool++);

        
      } // for (fool=elecids.begin(); fool != elecids.end();)


    } // if (proc==myrank)

#if 0
    int size = (int)sendeles.size();
    dis.Comm().Broadcast(&size,1,proc);

    if (proc != myrank) 
    {
      sendeles.resize(size);
      sendelemis.resize(size);
    }
    dis.Comm().Broadcast(&sendeles[0],size,proc);
    dis.Comm().Broadcast(&sendelemis[0],size,proc);

    // all other procs remove the already taken elements from their list
    // of not-yet taken elements
    if (myrank != proc)
    {
      // delete already taken elements from my column element map
      for (int i=0; i<size; ++i) elecids.erase(sendeles[i]);
      
      // look whether I have any of the communicated elements in my column map
      // If so, put the corresponding MIS node in my map
      for (int i=0; i<(int)sendeles.size(); ++i)
      {
        fool = elecids_.find(sendeles[i]);
        if (fool == elecids_full.end()) continue;
        pstab_cid_mis[sendeles[i]].push_back(sendelemis[i]);
      }
    } // if (myrank != proc)
    
    sendeles.clear();
    sendelemis.clear();
    dis.Comm().Barrier();
#endif


  } // for (int proc=0; proc<numproc; ++proc)
  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::InitMISAdjacency(
                        map<int,vector<DRT::ELEMENTS::NStet*> >& pstab_adjele,
                        map<int,vector<DRT::ELEMENTS::NStet*> >& adjele,
                        map<int,map<int,DRT::Node*> >&           pstab_adjnode,
                        map<int,bool>&                           pstab_ident_patch,
                        map<int,vector<int> >&                   pstab_adjlm,
                        const int                                myrank,
                        const int                                numproc,
                        DRT::Discretization&                     dis)
{
  map<int,vector<DRT::ELEMENTS::NStet*> >::iterator mis;
  for (mis=pstab_adjele.begin(); mis != pstab_adjele.end(); ++mis)
  {
    int id = mis->first;
    int mispatchsize = (int)mis->second.size();
    int patchsize = (int)adjele[id].size();
    cout << "Proc " << myrank 
         << " MIS " << id
         << " mispatchsize " << mispatchsize
         << " patchsize " << patchsize << endl;
    if (mispatchsize==patchsize) 
    {
      mis->second.clear();
      pstab_ident_patch[id] = true;
    }
    else // patch not identical
    {
      // adjnode
      pstab_ident_patch[id] = false;
      map<int,DRT::Node*> nodepatch;
      vector<DRT::ELEMENTS::NStet*>& ele = mis->second;
      for (int j=0; j<(int)ele.size(); ++j) 
      {
        DRT::Node** nodes = ele[j]->Nodes();
        for (int k=0; k<ele[j]->NumNode(); ++k)
          nodepatch[nodes[k]->Id()] = nodes[k];
      }
      pstab_adjnode[id] = nodepatch;
      
      // lm array
      vector<int> lm;
      std::map<int,DRT::Node*>::iterator pnode;
      for (pnode=nodepatch.begin(); pnode != nodepatch.end(); ++pnode)
      {
        const vector<int>& dofs = dis.Dof(pnode->second);
        for (int j=0; j<(int)dofs.size(); ++j)
          lm.push_back(dofs[j]);
      }
      pstab_adjlm_[id] = lm;
    } // else
  } // for (mis=pstab_adjele.begin(); mis != pstab_adjele.end(); ++mis)
  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 09/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::ELEMENTS::NStetType::InitMISStressMap(
                                          map<int,vector<int> >& pstab_cid_mis, 
                                          DRT::Discretization& dis)
{
  map<int,vector<int> >::iterator fool;
  map<int,int> ngidmap;
  vector<int> ngid;
  for (fool=pstab_cid_mis.begin(); fool != pstab_cid_mis.end(); ++fool)
  {
    vector<int>& mis = fool->second;
    for (int i=0; i<(int)mis.size(); ++i)
      ngidmap[mis[i]] = mis[i];
  }
  
  map<int,int>::iterator fool2;
  for (fool2=ngidmap.begin(); fool2 != ngidmap.end(); ++fool2)
    ngid.push_back(fool2->first);
    
  return(rcp(new Epetra_Map(-1,(int)ngid.size(),&ngid[0],0,dis.Comm())));
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

  InitElementsandMaps(elecids_,noderids_,myrank,numproc,dis);
  
  //----------------------------------------------------------------------
  // compute adjacency for each row node
  // make patch of adjacent elements
  // make patch of adjacent nodes (including center node itself)
  // make lm for nodal patch

  InitAdjacency(elecids_,noderids_,adjele_,adjnode_,adjlm_,dis);
  
  //----------------------------------------------------------------------
  // build parallel maximum independent set of nodes (MIS-nodes)
  // this is done in pseudo-serial, as a true parallel MIS algorithm ispretty difficult
  
  map<int,int> misnodesmap;
  map<int,DRT::Node*> rnodes = noderids_; // working copy of row nodes
  InitMISnode(misnodesmap,rnodes,myrank,numproc,dis);
  

  //----------------------------------------------------------------------
  // each MIS node is associated with a patch of column elements of which
  // it takes the full integration area
  // Additionally, there are leftover column elements (not adjacent to any MIS node)
  // they are taken by a greedy algorithm in a second phase

  //----------------------------------------------------------------------

  map<int,DRT::ELEMENTS::NStet*> elecids = elecids_;
  InitMISpatchesGreedyI(misnodesmap,pstab_adjele_,pstab_cid_mis_,
                        elecids,elecids_,noderids_,myrank,numproc,dis);
  
 


  //----------------------------------------------------------------------
  // assign all leftover elements that are well connected to a patch to that patch
  // this is a distance 2 patch search

  //InitMISpatchesGreedyII(misnodesmap,pstab_adjele_,pstab_adjele_weight_,pstab_cid_mis_,
  //                       elecids,elecids_,noderids_,myrank,numproc,dis);

  
  //----------------------------------------------------------------------
  // split remaining elements among patches that its nodes belong to
  InitMISpatchesGreedyIII(misnodesmap,pstab_adjele_,pstab_adjele_weight_,pstab_cid_mis_,
                          elecids,elecids_,noderids_,myrank,numproc,dis);

  //----------------------------------------------------------------------
  // test whether all column elements on all procs have been assigned a patch
  map<int,DRT::ELEMENTS::NStet*>::iterator ele;
  if ((int)elecids.size() != 0)
  {
    for (ele=elecids.begin(); ele != elecids.end(); ++ele)
    {
      cout << "Proc " << myrank <<
              " leftover NStet with no MIS node " << *ele->second << endl;
    }
    dserror("Proc %d has the above column elements leftover",myrank);
  }
  
  // test whether all column elements on this proc know their MIS node_pos
  for (ele=elecids_.begin(); ele != elecids_.end(); ++ele)
  {
    map<int,vector<int> >::iterator fool = pstab_cid_mis_.find(ele->first);
    if (fool==pstab_cid_mis_.end())
    {
      cout << "This element did not find its MIS node:\n" << *ele->second << endl;
      dserror("Element %d did not find its MIS node",ele->first);
    }
  }

  //----------------------------------------------------------------------
  // have to build adjnode and adjlm arrays for the patches
  InitMISAdjacency(pstab_adjele_,adjele_,pstab_adjnode_,pstab_ident_patch_,pstab_adjlm_,
                   myrank,numproc,dis);


  //----------------------------------------------------------------------
  // create an overlapping map that contains stress data of MIS nodes on all procs
  // that will need it for stress output
  pstab_misstressout_ = InitMISStressMap(pstab_cid_mis_,dis);
  
  return 0;
}



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
