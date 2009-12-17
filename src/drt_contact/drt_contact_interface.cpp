/*!----------------------------------------------------------------------
\file drt_contact_interface.cpp
\brief One contact interface

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef PARALLEL
#include "Epetra_SerialComm.h"
#endif
#include "drt_contact_interface.H"
#include "drt_contact_integrator.H"
#include "drt_contact_coupling2d.H"
#include "drt_contact_coupling3d.H"
#include "drt_cdofset.H"
#include "contactdefines.H"
#include "drt_contact_binarytree.H"
#include "drt_contact_binarytree_self.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Interface::Interface(const int id, const Epetra_Comm& comm,
                              const int dim,
                              const Teuchos::ParameterList& icontact,
                              bool selfcontact) :
id_(id),
comm_(comm),
dim_(dim),
icontact_(icontact),
shapefcn_(Interface::Undefined),
selfcontact_(selfcontact)
{
  RCP<Epetra_Comm> com = rcp(Comm().Clone());
  if (Dim()!=2 && Dim()!=3) dserror("ERROR: Contact problem must be 2D or 3D");
  procmap_.clear();
  idiscret_ = rcp(new DRT::Discretization((string)"Contact Interface",com));
  contactsegs_.Reshape(0, 0);

  // overwrite shape function type
  INPAR::CONTACT::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::CONTACT::ShapeFcn>(IParams(),"SHAPEFCN");
  if (shapefcn == INPAR::CONTACT::shape_dual)
   shapefcn_ = Interface::DualFunctions;
  else if (shapefcn == INPAR::CONTACT::shape_standard)
    shapefcn_ = Interface::StandardFunctions;

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::Interface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "Contact Interface Id " << id_ << endl;
    os << "Contact Interface Discretization:" << endl;
  }
  os << Discret();
  return;
}

/*----------------------------------------------------------------------*
 |  finalize construction of interface (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FillComplete()
{
  // we'd like to call idiscret_.FillComplete(true,false,false) but this
  // will assign all nodes new degrees of freedom which we don't want.
  // We would like to use the degrees of freedom that were stored in the
  // contact nodes. To do so, we have to create and set our own
  // version of a DofSet class before we call FillComplete on the
  // interface discretization.
  // Our special dofset class will not assign new dofs but will assign the
  // dofs stored in the nodes.
  {
    RCP<CONTACT::CDofSet> cdofset = rcp(new CONTACT::CDofSet());
    Discret().ReplaceDofSet(cdofset);
    // do not assign dofs yet, we'll do this below after
    // shuffling around of nodes and elements (saves time)
    Discret().FillComplete(false,false,false);
  }

#ifdef CONTACTBOUNDMOD
  // only applicable for 2D problems up to now
  if (Dim()==3) dserror("ERROR: Boundary node modification not yet impl. for 3D");

  // detect boundary nodes on slave side
  // ---------------------------------------------------------------------
  // Within our dual Lagrange multiplier framework, results can be
  // improved by a modified treatment of these "boundary nodes". Their
  // status is changed to MASTER and consequently they will NOT carry
  // Lagrange multipliers later on. In order to sustain the partition
  // of unity property of the dual shape functions on the adjacent slave
  // elements, the dual shape functions of the adjacent nodes will be
  // modified later on! This way, the Mortar operator entries of these
  // "boundary nodes" are transfered to the neighboring slave nodes!
  // ---------------------------------------------------------------------
  for (int i=0; i<(Discret().NodeRowMap())->NumMyElements();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lRowNode(i));

    // candidates are slave nodes with only 1 adjacent CElement
    if (node->IsSlave() && node->NumElement()==1)
    {
      //case1: linear shape functions, boundary nodes already found
      if ((node->Elements()[0])->NumNode() == 2)
      {
        node->SetBound()=true;
        node->SetSlave()=false;
      }
      //case2: quad. shape functions, middle nodes must be sorted out
      else if (node->Id() != (node->Elements()[0])->NodeIds()[2])
      {
        node->SetBound()=true;
        node->SetSlave()=false;
      }
    }
  }
#endif // #ifdef CONTACTBOUNDMOD

  // later we will export node and element column map to FULL overlap,
  // thus store the standard column maps first
  // get standard nodal column map (overlap=1)
  oldnodecolmap_ = rcp (new Epetra_Map(*(Discret().NodeColMap())));
  // get standard element column map (overlap=1)
  oldelecolmap_ = rcp (new Epetra_Map(*(Discret().ElementColMap())));

  // create interface local communicator
  // find all procs that have business on this interface (own nodes/elements)
  // build a Epetra_Comm that contains only those procs
  // this intra-communicator will be used to handle most stuff on this
  // interface so the interface will not block all other procs
  {
#ifdef PARALLEL
    vector<int> lin(Comm().NumProc());
    vector<int> gin(Comm().NumProc());
    for (int i=0; i<Comm().NumProc(); ++i)
      lin[i] = 0;

    // check ownership of any elements / nodes
    const Epetra_Map* noderowmap = Discret().NodeRowMap();
    const Epetra_Map* elerowmap  = Discret().ElementRowMap();

    if (noderowmap->NumMyElements() || elerowmap->NumMyElements())
      lin[Comm().MyPID()] = 1;

    Comm().MaxAll(&lin[0],&gin[0],Comm().NumProc());
    lin.clear();

    // build global -> local communicator PID map
    // we need this when calling Broadcast() on lComm later
    int counter = 0;
    for (int i=0; i<Comm().NumProc(); ++i)
    {
      if (gin[i])
        procmap_[i]=counter++;
      else
        procmap_[i]=-1;
    }

    // typecast the Epetra_Comm to Epetra_MpiComm
    RCP<Epetra_Comm> copycomm = rcp(Comm().Clone());
    Epetra_MpiComm* epetrampicomm = dynamic_cast<Epetra_MpiComm*>(copycomm.get());
    if (!epetrampicomm)
      dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");

    // split the communicator into participating and none-participating procs
    int color;
    int key = Comm().MyPID();
    // I am taking part in the new comm if I have any ownership
    if (gin[Comm().MyPID()])
      color = 0;
    // I am not taking part in the new comm
    else
      color = MPI_UNDEFINED;

    // tidy up
    gin.clear();

    // create the local communicator
    MPI_Comm  mpi_global_comm = epetrampicomm->GetMpiComm();
    RCP<MPI_Comm> mpi_local_comm  = rcp(new MPI_Comm());
    MPI_Comm_split(mpi_global_comm,color,key,mpi_local_comm.get());

    // create the new Epetra_MpiComm
    if (*mpi_local_comm.get() == MPI_COMM_NULL)
      lcomm_ = null;
    else
      lcomm_ = rcp(new Epetra_MpiComm(*mpi_local_comm.get()));

#else  // the easy serial case
    RCP<Epetra_Comm> copycomm = rcp(Comm().Clone());
    Epetra_SerialComm* serialcomm = dynamic_cast<Epetra_SerialComm*>(copycomm.get());
    if (!serialcomm)
      dserror("ERROR: casting Epetra_Comm -> Epetra_SerialComm failed");
    lcomm_ = rcp(new Epetra_SerialComm(*serialcomm));
#endif // #ifdef PARALLEL
  }

  // to ease our search algorithms we'll afford the luxury to ghost all nodes
  // on all processors. To do so, we'll take the nodal row map and export it
  // to full overlap. Then we export the discretization to full overlap
  // column map. This way, also contact elements will be fully ghosted on all
  // processors.
  // Note that we'll do ghosting only on procs that do own or ghost any of the
  // nodes in the natural distribution of idiscret_!
  {
    const Epetra_Map* noderowmap = Discret().NodeRowMap();
    const Epetra_Map* elerowmap  = Discret().ElementRowMap();

    // fill my own row node ids
    vector<int> sdata(noderowmap->NumMyElements());
    for (int i=0; i<noderowmap->NumMyElements(); ++i)
      sdata[i] = noderowmap->GID(i);

    // build tprocs and numproc containing processors participating
    // in this interface
    vector<int> stproc(0);
    // a processor participates in the interface, if it owns any of
    // the nodes or elements
    if (noderowmap->NumMyElements() || elerowmap->NumMyElements())
      stproc.push_back(Comm().MyPID());
    vector<int> rtproc(0);
    vector<int> allproc(Comm().NumProc());
    for (int i=0; i<Comm().NumProc(); ++i) allproc[i] = i;
    LINALG::Gather<int>(stproc,rtproc,Comm().NumProc(),&allproc[0],Comm());
    vector<int> rdata;

    // gather all gids of nodes redundantly
    LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],Comm());

    // build completely overlapping map (on participating processors)
    RCP<Epetra_Map> newnodecolmap = rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,Comm()));
    sdata.clear();
    stproc.clear();
    rdata.clear();
    allproc.clear();
    // rtproc still in use

    // do the same business for elements
    sdata.resize(elerowmap->NumMyElements());
    rdata.resize(0);
    for (int i=0; i<elerowmap->NumMyElements(); ++i)
      sdata[i] = elerowmap->GID(i);
    // gather of element gids redundantly on processors that have business on the interface
    LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],Comm());
    rtproc.clear();

    // build complete overlapping map of elements (on participating processors)
    RCP<Epetra_Map> newelecolmap = rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,Comm()));
    sdata.clear();
    rdata.clear();

    // redistribute the discretization of the interface according to the
    // new column layout
    Discret().ExportColumnNodes(*newnodecolmap);
    Discret().ExportColumnElements(*newelecolmap);

    // make sure discretization is complete
    Discret().FillComplete(true,false,false);
  }

  // need row and column maps of slave and master nodes / elements / dofs
  // separately so we can easily adress them
  UpdateMasterSlaveSets();

  // get out of here if not participating in interface
  if (!lComm()) return;

  // warning
#ifdef CONTACTGMSHCTN
  if (Dim()==3 && Comm().MyPID()==0)
  {
    cout << "\n******************************************************************\n";
    cout << "GMSH output of all contact tree nodes in 3D needs a lot of memory!\n";
    cout << "******************************************************************\n";
  }
#endif //CONTACTGMSHCTN

  // binary tree search
  if (SearchAlg()==INPAR::CONTACT::search_binarytree)
  {
    //*****SELF CONTACT*****
    if (SelfContact())
    {
      // set state in interface to intialize all kinds of quantities
      RCP<Epetra_Vector> zero =rcp(new Epetra_Vector(*idiscret_->DofRowMap()));
      SetState("displacement",zero);

      // create fully overlapping map of all contact elements
      RCP<Epetra_Map> elefullmap = rcp(new Epetra_Map(*idiscret_->ElementColMap()));

      // create binary tree object for self contact search
      binarytreeself_ = rcp(new CONTACT::BinaryTreeSelf(Discret(),elefullmap,Dim(),SearchParam()));

    }
    //*****TWO BODY CONTACT*****
    else
    {
      // create binary tree object for contact search and setup tree
  	  binarytree_ = rcp(new CONTACT::BinaryTree(Discret(),selecolmap_,melefullmap_,Dim(),SearchParam()));

  	  // initialize active contact nodes via binarytree
  	  // binarytree_->SearchContactInit(binarytree_->Sroot(), binarytree_->Mroot());
    }
  }
  // no binary tree search
  else
  {
    if (SelfContact()) dserror("ERROR: Binarytree search needed for self contact");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                 popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::UpdateMasterSlaveSets()
{
  //********************************************************************
  // NODES
  //********************************************************************
  // need row and column maps of slave and master nodes separately so we
  // can easily adress them
  {
    const Epetra_Map* noderowmap = Discret().NodeRowMap();
    const Epetra_Map* nodecolmap = Discret().NodeColMap();

    vector<int> sc;          // slave column map
    vector<int> sr;          // slave row map
    vector<int> scfull;      // slave full map
    vector<int> mc;          // master column map
    vector<int> mr;          // master row map
    vector<int> mcfull;      // master full map
    vector<int> scb;         // slave column map + boundary nodes
    vector<int> scfullb;     // slave full map + boundary nodes
    vector<int> mcfullb;     // master full map - boundary nodes

    for (int i=0; i<nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      bool isslave = dynamic_cast<CONTACT::CNode*>(Discret().gNode(gid))->IsSlave();
      bool isonbound = dynamic_cast<CONTACT::CNode*>(Discret().gNode(gid))->IsOnBound();
      if (oldnodecolmap_->MyGID(gid))
      {
        if (isslave || isonbound) scb.push_back(gid);
        if (isslave) sc.push_back(gid);
        else         mc.push_back(gid);
      }
      if (isslave || isonbound) scfullb.push_back(gid);
      else                      mcfullb.push_back(gid);
      if (isslave) scfull.push_back(gid);
      else         mcfull.push_back(gid);
      if (!noderowmap->MyGID(gid)) continue;
      if (isslave) sr.push_back(gid);
      else         mr.push_back(gid);
    }

    snoderowmap_ = rcp(new Epetra_Map(-1,(int)sr.size(),&sr[0],0,Comm()));
    snodefullmap_ = rcp(new Epetra_Map(-1,(int)scfull.size(),&scfull[0],0,Comm()));
    snodecolmap_ = rcp(new Epetra_Map(-1,(int)sc.size(),&sc[0],0,Comm()));
    mnoderowmap_ = rcp(new Epetra_Map(-1,(int)mr.size(),&mr[0],0,Comm()));
    mnodefullmap_ = rcp(new Epetra_Map(-1,(int)mcfull.size(),&mcfull[0],0,Comm()));
    mnodecolmap_ = rcp(new Epetra_Map(-1,(int)mc.size(),&mc[0],0,Comm()));

    snodecolmapbound_ = rcp(new Epetra_Map(-1,(int)scb.size(),&scb[0],0,Comm()));
    snodefullmapbound_ = rcp(new Epetra_Map(-1,(int)scfullb.size(),&scfullb[0],0,Comm()));
    mnodefullmapnobound_ = rcp(new Epetra_Map(-1,(int)mcfullb.size(),&mcfullb[0],0,Comm()));
  }

  //********************************************************************
  // ELEMENTS
  //********************************************************************
  // do the same business for elements
  // (get row and column maps of slave and master elements seperately)
  {
    const Epetra_Map* elerowmap = Discret().ElementRowMap();
    const Epetra_Map* elecolmap = Discret().ElementColMap();

    vector<int> sc;          // slave column map
    vector<int> sr;          // slave row map
    vector<int> scfull;      // slave full map
    vector<int> mc;          // master column map
    vector<int> mr;          // master row map
    vector<int> mcfull;      // master full map

    for (int i=0; i<elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      bool isslave = dynamic_cast<CONTACT::CElement*>(Discret().gElement(gid))->IsSlave();
      if (oldelecolmap_->MyGID(gid))
      {
        if (isslave) sc.push_back(gid);
        else         mc.push_back(gid);
      }
      if (isslave) scfull.push_back(gid);
      else         mcfull.push_back(gid);
      if (!elerowmap->MyGID(gid)) continue;
      if (isslave) sr.push_back(gid);
      else         mr.push_back(gid);
    }

    selerowmap_ = rcp(new Epetra_Map(-1,(int)sr.size(),&sr[0],0,Comm()));
    selefullmap_ = rcp(new Epetra_Map(-1,(int)scfull.size(),&scfull[0],0,Comm()));
    selecolmap_ = rcp(new Epetra_Map(-1,(int)sc.size(),&sc[0],0,Comm()));
    melerowmap_ = rcp(new Epetra_Map(-1,(int)mr.size(),&mr[0],0,Comm()));
    melefullmap_ = rcp(new Epetra_Map(-1,(int)mcfull.size(),&mcfull[0],0,Comm()));
    melecolmap_ = rcp(new Epetra_Map(-1,(int)mc.size(),&mc[0],0,Comm()));
  }

  //********************************************************************
  // DOFS
  //********************************************************************
  // do the same business for dofs
  // (get row and column maps of slave and master dofs seperately)
  {
    const Epetra_Map* noderowmap = Discret().NodeRowMap();
    const Epetra_Map* nodecolmap = Discret().NodeColMap();

    vector<int> sc;          // slave column map
    vector<int> sr;          // slave row map
    vector<int> scfull;      // slave full map
    vector<int> mc;          // master column map
    vector<int> mr;          // master row map
    vector<int> mcfull;      // master full map

    for (int i=0; i<nodecolmap->NumMyElements();++i)
    {
      int gid = nodecolmap->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);
      bool isslave = cnode->IsSlave();

      if (oldnodecolmap_->MyGID(gid))
      {
        if (isslave)
          for (int j=0;j<cnode->NumDof();++j)
            sc.push_back(cnode->Dofs()[j]);
        else
          for (int j=0;j<cnode->NumDof();++j)
            mc.push_back(cnode->Dofs()[j]);
      }

      if (isslave)
        for (int j=0;j<cnode->NumDof();++j)
          scfull.push_back(cnode->Dofs()[j]);
      else
        for (int j=0;j<cnode->NumDof();++j)
          mcfull.push_back(cnode->Dofs()[j]);

      if (!noderowmap->MyGID(gid)) continue;

      if (isslave)
        for (int j=0;j<cnode->NumDof();++j)
          sr.push_back(cnode->Dofs()[j]);
      else
        for (int j=0;j<cnode->NumDof();++j)
          mr.push_back(cnode->Dofs()[j]);
    }

    sdofrowmap_ = rcp(new Epetra_Map(-1,(int)sr.size(),&sr[0],0,Comm()));
    sdoffullmap_ = rcp(new Epetra_Map(-1,(int)scfull.size(),&scfull[0],0,Comm()));
    sdofcolmap_ = rcp(new Epetra_Map(-1,(int)sc.size(),&sc[0],0,Comm()));
    mdofrowmap_ = rcp(new Epetra_Map(-1,(int)mr.size(),&mr[0],0,Comm()));
    mdoffullmap_ = rcp(new Epetra_Map(-1,(int)mcfull.size(),&mcfull[0],0,Comm()));
    mdofcolmap_ = rcp(new Epetra_Map(-1,(int)mc.size(),&mc[0],0,Comm()));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  initialize / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

    //reset nodal normal and tangents and jumps
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
      //node->jump()[j]=0.0;
    }

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((node->GetDerivN()).size());++j)
      (node->GetDerivN())[j].clear();
    (node->GetDerivN()).resize(0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
      (node->GetDerivTxi())[j].clear();
    (node->GetDerivTxi()).resize(0);
    for (int j=0;j<(int)((node->GetDerivTeta()).size());++j)
      (node->GetDerivTeta())[j].clear();
    (node->GetDerivTeta()).resize(0);

    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();
    for (int j=0;j<(int)((node->GetM()).size());++j)
      (node->GetM())[j].clear();
    for (int j=0;j<(int)((node->GetMmod()).size());++j)
      (node->GetMmod())[j].clear();

    (node->GetD()).resize(0);
    (node->GetM()).resize(0);
    (node->GetMmod()).resize(0);

    // reset derivative map of Mortar matrices
    (node->GetDerivD()).clear();
    (node->GetDerivM()).clear();

    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((node->GetDerivZ()).size()); ++j)
      (node->GetDerivZ())[j].clear();
    (node->GetDerivZ()).resize(0);

    // reset derivative map of jump
    for (int j=0; j<(int)((node->GetDerivJump()).size()); ++j)
      (node->GetDerivJump())[j].clear();
    (node->GetDerivJump()).resize(0);

    // reset nodal weighted gap and derivative
    node->Getg() = 1.0e12;
    (node->GetDerivG()).clear();

    // reset SNodes and Mnodes
    node->GetSNodes().clear();
    node->GetMNodes().clear();

    // reset feasible projection status
    node->HasProj() = false;
  }

  // loop over all elements to reset contact candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  // reset matrix containing interface contact segments (gmsh)
  CSegs().Shape(0,0);

  return;
}

/*----------------------------------------------------------------------*
 |  set current and old deformation state                      popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::SetState(const string& statename, const RCP<Epetra_Vector> vec)
{
  // ***WARNING:*** This is not possible here, as idiscret_->SetState()
  // needs all procs around, not only the interface local ones!
  // get out of here if not participating in interface
  // if (!lComm())
  //   return;

  if (statename=="displacement")
  {
    // set displacements in interface discretization
    idiscret_->SetState(statename, vec);

    // Get vec to full overlap
    // RCP<const Epetra_Vector> global = idiscret_->GetState(statename);

    // alternative method to get vec to full overlap
    Epetra_Vector global(*idiscret_->DofColMap(),false);
    LINALG::Export(*vec,global);

    // loop over all nodes to set current displacement
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));
      const int numdof = node->NumDof();
      vector<double> mydisp(numdof);
      vector<int> lm(numdof);

      for (int j=0;j<numdof;++j)
        lm[j]=node->Dofs()[j];

      DRT::UTILS::ExtractMyValues(global,mydisp,lm);

      // add mydisp[2]=0 for 2D problems
      if (mydisp.size()<3)
        mydisp.resize(3);

      // set current configuration and displacement
      for (int j=0;j<3;++j)
      {
        node->u()[j]=mydisp[j];
        node->xspatial()[j]=node->X()[j]+mydisp[j];
      }
    }

    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->Area()=element->ComputeArea();
    }
  }

  if (statename=="olddisplacement")
  {
    // set displacements in interface discretization
    idiscret_->SetState(statename, vec);

    // Get vec to full overlap
    // RCP<const Epetra_Vector> global = idiscret_->GetState(statename);

    // alternative method to get vec to full overlap
    Epetra_Vector global(*idiscret_->DofColMap(),false);
    LINALG::Export(*vec,global);

    // loop over all nodes to set current displacement
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));
      const int numdof = node->NumDof();
      vector<double> myolddisp(numdof);
      vector<int> lm(numdof);

      for (int j=0;j<numdof;++j)
        lm[j]=node->Dofs()[j];

      DRT::UTILS::ExtractMyValues(global,myolddisp,lm);

      // add mydisp[2]=0 for 2D problems
      if (myolddisp.size()<3)
        myolddisp.resize(3);

      // set current configuration and displacement
      for (int j=0;j<3;++j)
      {
        node->uold()[j]=myolddisp[j];
      }
    }

    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->Area()=element->ComputeArea();
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Evaluate()
{
  // interface needs to be complete
  if (!Filled() && Comm().MyPID()==0)
      dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm()) return;

#ifdef CONTACTFDNORMAL
  // FD check of normal derivatives
  FDCheckNormalDeriv();
#endif // #ifdef CONTACTFDNORMAL

  // loop over proc's slave nodes of the interface
  // use standard column map to include processor's ghosted nodes
  // use boundary map to include slave side boundary nodes
  for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    int gid = snodecolmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    // build averaged normal at each slave node
    cnode->BuildAveragedNormal();
  }
  
  //**********************************************************************
  // contact search algorithm
  //**********************************************************************
  //lComm()->Barrier();
  //const double t_start = ds_cputime();

  if (SearchAlg()==INPAR::CONTACT::search_bfnode)          EvaluateContactSearch();
  else if (SearchAlg()==INPAR::CONTACT::search_bfele)      EvaluateContactSearchBruteForce(SearchParam());
  else if (SearchAlg()==INPAR::CONTACT::search_binarytree) EvaluateContactSearchBinarytree();
  else                                                     dserror("ERROR: Invalid contact search algorithm");

  //lComm()->Barrier();
  //const double t_end = ds_cputime()-t_start;
  //if (lComm()->MyPID()==0)
  //  cout << "Search Time (overall): " << t_end << " seconds\n";

#ifdef CONTACTFDVERTEX3D
  // define test variable for FDVertex
  vector<vector<double> > testv(0,vector<double>(6));
#endif // #ifdef CONTACTFDVERTEX3D

#ifdef CONTACTFDGP3D
  // define test variable for FDGP
  vector<vector<double> > testgps(0,vector<double>(12));
  vector<vector<double> > testgpm(0,vector<double>(12));
  vector<vector<double> > testjs(0,vector<double>(6));
  vector<vector<double> > testji(0,vector<double>(6));
#endif // #ifdef CONTACTFDGP3D

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);

    //********************************************************************
    // 1) integrate Mortar matrix D (lives on slave side only!)
    // 2) compute directional derivative of D and store into nodes
    //********************************************************************
    // we only do this on a slave element basis in the case with
    // two separate Mortar loops for D and M. If CONTACTONEMORTARLOOP
    // is applied, we combine the construction of D with the construction
    // of M and linearization of D with the linearization of M in
    // Integrator::IntegrateDerivSegment2D()!
    //********************************************************************
#ifndef CONTACTONEMORTARLOOP
    IntegrateSlave(*selement);
#endif // #ifndef CONTACTONEMORTARLOOP

    // loop over the contact candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CElement* melement = static_cast<CElement*>(ele2);

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pair)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
#ifdef CONTACTFDVERTEX3D
      IntegrateCoupling(*selement,*melement,testv);
#else
#ifdef CONTACTFDGP3D
      IntegrateCoupling(*selement,*melement,testgps,testgpm,testjs,testji);
#else
      // this is the standard case without any FD checks
      IntegrateCoupling(*selement,*melement);
#endif // #ifdef CONTACTFDGP3D
#endif // #ifdef CONTACTFDVERTEX3D
    }
  }
#ifdef CONTACTFDVERTEX3D
  // FD check of coupling vertex derivatives (3D)
  FDCheckVertex3DDeriv(testv);
#endif // #ifdef CONTACTFDVERTEX3D

#ifdef CONTACTFDGP3D
  // FD check of Gauss points (3D)
  FDCheckGP3DDeriv(testgps,testgpm,testjs,testji);
#endif // #ifdef CONTACTFDGP3D

  return;
}

/*----------------------------------------------------------------------*
 |  Contact search node-based ("brute force") (public)        popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::EvaluateContactSearch()
{
  /**********************************************************************/
  /* CONTACT SEARCH ALGORITHM:                                          */
  /* The idea of the search is to reduce the number of master / slave   */
  /* element pairs that are checked for overlap and contact by intro-   */
  /* ducing information about proximity and maybe history!              */
  /* This old version is still brute force for finding the closest      */
  /* CNode to each CNode, so it has been replaced by a more             */
  /* sophisticated approach (bounding vol. hierarchies /  binary tree). */
  /**********************************************************************/

  // loop over proc's slave nodes for closest node detection
  // use standard column map to include processor's ghosted nodes
  // use boundary map to include slave boundary nodes
  for (int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    int gid = snodecolmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);

    // find closest master node to current slave node
    double mindist = 1.0e12;
    CNode* closestnode = snode->FindClosestNode(idiscret_,mnodefullmapnobound_,mindist);
    snode->ClosestNode() = closestnode->Id();

    // proceed only if nodes are not far from each other!!!
    if (mindist<=SearchParam())
    {
      // get adjacent elements to current slave node and to closest node
      int neles = snode->NumElement();
      DRT::Element** adjslave = snode->Elements();
      int nelec = closestnode->NumElement();
      DRT::Element** adjclosest = closestnode->Elements();

      // get global element ids for closest node's adjacent elements
      std::vector<int> cids(nelec);
      for (int j=0;j<nelec;++j)
        cids[j]=adjclosest[j]->Id();

      // try to add these to slave node's adjacent elements' search list
      for (int j=0;j<neles;++j)
      {
        CElement* selement = static_cast<CElement*> (adjslave[j]);
        selement->AddSearchElements(cids);
      }
    }
  }

  // loop over all master nodes for closest node detection
  // use full overlap column map to include all nodes
  // use no boundary map to exclude slave side boundary nodes
  for (int i=0; i<mnodefullmapnobound_->NumMyElements();++i)
  {
    int gid = mnodefullmapnobound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find master node with gid %",gid);
    CNode* mnode = static_cast<CNode*>(node);

    // find closest slave node to current master node
    double mindist = 1.0e12;
    CNode* closestnode = mnode->FindClosestNode(idiscret_,snodefullmapbound_,mindist);
    mnode->ClosestNode() = closestnode->Id();

    // proceed only if nodes are not far from each other!!!
    if (mindist<=SearchParam())
    {
      // get adjacent elements to current master node and to closest node
      int nelem = mnode->NumElement();
      DRT::Element** adjmaster = mnode->Elements();
      int nelec = closestnode->NumElement();
      DRT::Element** adjclosest = closestnode->Elements();

      // get global element ids for master node's adjacent elements
      std::vector<int> mids(nelem);
      for (int j=0;j<nelem;++j)
        mids[j]=adjmaster[j]->Id();

      // try to add these to closest node's adjacent elements' search list
      for (int j=0;j<nelec;++j)
      {
        CElement* selement = static_cast<CElement*> (adjclosest[j]);
        selement->AddSearchElements(mids);
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Contact search element-based "brute force" (public)       popp 10/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateContactSearchBruteForce(const double& eps)
{
  // calculate minimal element length
  double lmin = 1.0e12;
  double enlarge = 0.0;

  // loop over all slave elements on this proc.
  for (int i=0;i<selecolmap_->NumMyElements();++i)
  {
    DRT::Element* element = idiscret_->gElement(selecolmap_->GID(0));
    if (!element) dserror("ERROR: Cannot find element with gid %\n",selecolmap_->GID(0));
    CONTACT::CElement* celement = (CElement*) element;
    if (celement->MinEdgeSize() < lmin)
      lmin= celement->MinEdgeSize();
  }

  // loop over all master elements on this proc.
  for (int i=0;i<melefullmap_->NumMyElements();++i)
  {
    DRT::Element* element = idiscret_->gElement(melefullmap_->GID(i));
    if (!element) dserror("ERROR: Cannot find element with gid %\n",melefullmap_->GID(i));
    CONTACT::CElement* celement = (CElement*) element;
    if (celement->MinEdgeSize() < lmin)
      lmin= celement->MinEdgeSize();
  }

  // compute DOP inflation length
  enlarge=eps*lmin;

  // define dopnormals
  Epetra_SerialDenseMatrix dopnormals;
  int kdop=0;

  if (dim_==2)
  {
    kdop=8;

    // setup normals for DOP
    dopnormals.Reshape(4,3);
    dopnormals(0,0)= 1; dopnormals(0,1)= 0; dopnormals(0,2)= 0;
    dopnormals(1,0)= 0; dopnormals(1,1)= 1; dopnormals(1,2)= 0;
    dopnormals(2,0)= 1; dopnormals(2,1)= 1; dopnormals(2,2)= 0;
    dopnormals(3,0)=-1; dopnormals(3,1)= 1; dopnormals(3,2)= 0;
  }
  else if (dim_==3)
  {
    kdop=18;

    // setup normals for DOP
    dopnormals.Reshape(9,3);
    dopnormals(0,0)= 1; dopnormals(0,1)= 0; dopnormals(0,2)= 0;
    dopnormals(1,0)= 0; dopnormals(1,1)= 1; dopnormals(1,2)= 0;
    dopnormals(2,0)= 0; dopnormals(2,1)= 0; dopnormals(2,2)= 1;
    dopnormals(3,0)= 1; dopnormals(3,1)= 1; dopnormals(3,2)= 0;
    dopnormals(4,0)= 1; dopnormals(4,1)= 0; dopnormals(4,2)= 1;
    dopnormals(5,0)= 0; dopnormals(5,1)= 1; dopnormals(5,2)= 1;
    dopnormals(6,0)= 1; dopnormals(6,1)= 0; dopnormals(6,2)=-1;
    dopnormals(7,0)=-1; dopnormals(7,1)= 1; dopnormals(7,2)= 0;
    dopnormals(8,0)= 0; dopnormals(8,1)=-1; dopnormals(8,2)= 1;
  }
  else
    dserror("ERROR: Problem dimension must be 2D or 3D!");

  // define slave and master slabs
  Epetra_SerialDenseMatrix sslabs(kdop/2,2);
  Epetra_SerialDenseMatrix mslabs(kdop/2,2);

  //**********************************************************************
  // perform brute-force contact search (element-based)
  //**********************************************************************
  // for every slave element
  for (int i=0; i<selecolmap_->NumMyElements();i++)
  {
    // calculate slabs
    double dcurrent = 0.0;

    //initialize slabs with first node
    int sgid=selecolmap_->GID(i);
    DRT::Element* element= idiscret_->gElement(sgid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n",sgid);
    DRT::Node** node= element->Nodes();
    CNode* cnode=static_cast<CNode*>(node[0]);
    const double* posnode = cnode->xspatial();

    // calculate slabs initialization
    for (int j=0; j<kdop/2; j++)
    {
      //= ax+by+cz=d/sqrt(aa+bb+cc)
      sslabs(j,0)=sslabs(j,1) = (dopnormals(j,0)*posnode[0]+dopnormals(j,1)*posnode[1]+dopnormals(j,2)*posnode[2])
        /sqrt((dopnormals(j,0)*dopnormals(j,0))+(dopnormals(j,1)*dopnormals(j,1))+(dopnormals(j,2)*dopnormals(j,2)));
    }

    // for int j=1, because of initialization done before
    for (int j=1;j<element->NumNode();j++)
    {
      CNode* cnode=static_cast<CNode*>(node[j]);
      posnode = cnode->xspatial();

      for(int k=0;k<kdop/2;k++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        dcurrent = (dopnormals(k,0)*posnode[0]+dopnormals(k,1)*posnode[1]+dopnormals(k,2)*posnode[2])
          /sqrt((dopnormals(k,0)*dopnormals(k,0))+(dopnormals(k,1)*dopnormals(k,1))+(dopnormals(k,2)*dopnormals(k,2)));
        if (dcurrent > sslabs(k,1))
          sslabs(k,1)=dcurrent;
        if (dcurrent < sslabs(k,0))
          sslabs(k,0)=dcurrent;
      }
    }

    // add auxiliary positions
    // (last converged positions for all slave nodes)
    for (int j=0;j<element->NumNode();j++)
    {
      //get pointer to slave node
      CNode* cnode=static_cast<CNode*>(node[j]);

      double auxpos [3] = {0.0, 0.0, 0.0};
      double scalar=0.0;
      for (int k=0; k<dim_;k++)
        scalar+=(cnode->X()[k]+cnode->uold()[k]-cnode->xspatial()[k])*cnode->n()[k];
      for (int k=0;k<dim_;k++)
        auxpos[k]= cnode->xspatial()[k]+scalar*cnode->n()[k];

      for(int j=0;j<kdop/2;j++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        dcurrent = (dopnormals(j,0)*auxpos[0]+dopnormals(j,1)*auxpos[1]+dopnormals(j,2)*auxpos[2])
          /sqrt((dopnormals(j,0)*dopnormals(j,0))+(dopnormals(j,1)*dopnormals(j,1))+(dopnormals(j,2)*dopnormals(j,2)));
        if (dcurrent > sslabs(j,1))
          sslabs(j,1)=dcurrent;
        if (dcurrent < sslabs(j,0))
          sslabs(j,0)=dcurrent;
      }
    }

    // enlarge slabs with scalar factor
    for (int j=0;j<kdop/2;j++)
    {
      sslabs(j,0)=sslabs(j,0)-enlarge;
      sslabs(j,1)=sslabs(j,1)+enlarge;
    }

    // for every master element
    for (int j=0; j<melefullmap_->NumMyElements();j++)
    {
      // calculate slabs
      double dcurrent = 0.0;

      //initialize slabs with first node
      int mgid=melefullmap_->GID(j);
      DRT::Element* element= idiscret_->gElement(mgid);
      if (!element) dserror("ERROR: Cannot find element with gid %\n",mgid);
      DRT::Node** node= element->Nodes();
      CNode* cnode=static_cast<CNode*>(node[0]);
      const double* posnode = cnode->xspatial();

      // calculate slabs initialization
      for (int k=0; k<kdop/2;k++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        mslabs(k,0)=mslabs(k,1) = (dopnormals(k,0)*posnode[0]+dopnormals(k,1)*posnode[1]+dopnormals(k,2)*posnode[2])
          /sqrt((dopnormals(k,0)*dopnormals(k,0))+(dopnormals(k,1)*dopnormals(k,1))+(dopnormals(k,2)*dopnormals(k,2)));
      }

      // for int k=1, because of initialization done before
      for (int k=1;k<element->NumNode();k++)
      {
        CNode* cnode=static_cast<CNode*>(node[k]);
        posnode = cnode->xspatial();

        for(int l=0; l<kdop/2; l++)
        {
          //= d=ax+by+cz/sqrt(aa+bb+cc)
          dcurrent = (dopnormals(l,0)*posnode[0]+dopnormals(l,1)*posnode[1]+dopnormals(l,2)*posnode[2])
            /sqrt((dopnormals(l,0)*dopnormals(l,0))+(dopnormals(l,1)*dopnormals(l,1))+(dopnormals(l,2)*dopnormals(l,2)));
          if (dcurrent > mslabs(l,1))
            mslabs(l,1)=dcurrent;
          if (dcurrent < mslabs(l,0))
            mslabs(l,0)=dcurrent;
        }
      }

      // enlarge slabs with scalar factor
      for (int k=0 ; k<kdop/2 ; k++)
      {
        mslabs(k,0)=mslabs(k,0)-enlarge;
        mslabs(k,1)=mslabs(k,1)+enlarge;
      }

      // check if slabs of current master and slave element intercept
      int nintercepts=0;
      for (int k=0;k<kdop/2;k++)
      {
        if ((sslabs(k,0)<=mslabs(k,0)&&sslabs(k,1)>=mslabs(k,0))
          ||(mslabs(k,1)>=sslabs(k,0)&&mslabs(k,0)<=sslabs(k,0))
          ||(sslabs(k,0)<=mslabs(k,0)&&sslabs(k,1)>=mslabs(k,1))
          ||(sslabs(k,0)>=mslabs(k,0)&&mslabs(k,1)>=sslabs(k,1)))
        {
          nintercepts++;
        }
      }

      //cout <<"\n"<< Comm().MyPID() << " Number of intercepts found: " << nintercepts ;

      // slabs of current master and slave element do intercept
      if (nintercepts==kdop/2)
      {
        //cout << Comm().MyPID() << "\nContact found between slave element: "
        //     << sgid <<" and master element: "<< mgid;
        DRT::Element* element= idiscret_->gElement(sgid);
        CONTACT::CElement* selement = static_cast<CONTACT::CElement*>(element);
        selement->AddSearchElements(mgid);
      }
    } // for all master elements
  } // for all slave elements

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 10/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::EvaluateContactSearchBinarytree()
{
  // get out of here if not participating in interface
  if (!lComm()) return true;

  // *********************************************************************
  // Possible versions for self contact:
  // *********************************************************************
  //
  // 1) Combined Update and Contact Search
  // -> In this case we only have to call SearchContactCombined(), which
  //    does both bottom-up update and search. Then the dynamics master/slave
  //    assignment routine UpdateMasterSlaveSets() is called.
  //
  // *********************************************************************
  if (SelfContact())
  {
    // calculate minimal element length
    binarytreeself_->SetEnlarge(false);

    // search for contact with an combined algorithm
    binarytreeself_->SearchContactCombined();

    // update master/slave sets of interface
    UpdateMasterSlaveSets();
  }

  // *********************************************************************
  // Possible versions for 2-body contact:
  // *********************************************************************
  //
  // 1) Combined Update and Contact Search
  // -> In this case we only have to call SearchCOntactCombined(), which
  //    does buth top-down update (where necessary) and search.
  //
  // 2) Separate Update and Contact Search
  // -> In this case we have to explicitly call and updating routine, i.e.
  //    UpdateTreeTopDown() or UpdateTreeBottomUp() before calling the
  //    search routine SearchContactSeparate(). Of course, the bottom-up
  //    update makes more sense here!
  //
  // *********************************************************************
  else
  {
    // calculate minimal element length
    binarytree_->SetEnlarge(false);

    // update tree in a top down way
    //binarytree_->UpdateTreeTopDown();

    // update tree in a bottom up way
    //binarytree_->UpdateTreeBottomUp();

  #ifdef CONTACTGMSHCTN
    for (int i=0;i<(int)(binarytree_->ContactMap().size());i++)
      binarytree_->ContactMap()[i].clear();
    binarytree_->ContactMap().clear();
    binarytree_->ContactMap().resize(2);
  #endif //CONTACTGMSHCTN

    // search for contact with a separate algorithm
    //binarytree_->SearchContactSeparate();

    // search for contact with an combined algorithm
    binarytree_->SearchContactCombined();
  }

	return true;
}

/*----------------------------------------------------------------------*
 |  Integrate Mortar matrix D on slave element (public)       popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateSlave(CONTACT::CElement& sele)
{
  // create an integrator instance with correct NumGP and Dim
  CONTACT::Integrator integrator(shapefcn_,sele.Shape());

  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // do the element integration (integrate and linearize D)
  int nrow = sele.NumNode();
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),nrow*Dim()));
  integrator.IntegrateDerivSlave2D3D(sele,sxia,sxib,dseg);

  // do the assembly into the slave nodes
  integrator.AssembleD(Comm(),sele,*dseg);

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlap      popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateCoupling(CONTACT::CElement& sele,
                                           CONTACT::CElement& mele)
{
  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (Dim()==2)
    CONTACT::Coupling2d coup(shapefcn_,Discret(),Dim(),sele,mele,CSegs());
  // ************************************************************** 3D ***
  else if (Dim()==3)
  {
    bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

    // ************************************************** quadratic 3D ***
    if (sele.IsQuad() || mele.IsQuad())
    {
      // only for auxiliary plane 3D version
      if (!auxplane) dserror("ERROR: Quadratic 3D contact only for AuxPlane case!");

      // build linear integration elements from quadratic CElements
      vector<RCP<CONTACT::IntElement> > sauxelements(0);
      vector<RCP<CONTACT::IntElement> > mauxelements(0);
      SplitIntElements(sele,sauxelements);
      SplitIntElements(mele,mauxelements);

      // loop over all IntElement pairs for coupling
      for (int i=0;i<(int)sauxelements.size();++i)
      {
        for (int j=0;j<(int)mauxelements.size();++j)
        {
          // create instance of coupling class
          CONTACT::Coupling3dQuad coup(shapefcn_,Discret(),Dim(),true,auxplane,
                        sele,mele,*sauxelements[i],*mauxelements[j]);
          // do coupling
          coup.EvaluateCoupling();
        }
      }
    }

    // ***************************************************** linear 3D ***
    else
    {
      // create instance of coupling class
      CONTACT::Coupling3d coup(shapefcn_,Discret(),Dim(),false,auxplane,sele,mele);
      // do coupling
      coup.EvaluateCoupling();
    }
  }
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlap      popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateCoupling(CONTACT::CElement& sele,
                                           CONTACT::CElement& mele,
                                           vector<vector<double> >& testv,
                                           bool printderiv)
{
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  if (Dim()==2)
    dserror("ERROR: FD check of coupling vertices only for 3D!");
  else if (Dim()==3)
  {
    bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");
    // here the EvaluateCoupling routine is part of the constructor!
    CONTACT::Coupling3d coup(Discret(),Dim(),false,auxplane,sele,mele,testv,printderiv);
  }
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlap      popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateCoupling(CONTACT::CElement& sele,
                                           CONTACT::CElement& mele,
                                           vector<vector<double> >& testgps,
                                           vector<vector<double> >& testgpm,
                                           vector<vector<double> >& testjs,
                                           vector<vector<double> >& testji,
                                           bool printderiv)
{
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  if (Dim()==2)
    dserror("ERROR: FD check of Gauss points and Jacobian only for 3D!");
  else if (Dim()==3)
  {
    bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");
    // here the EvaluateCoupling routine is part of the constructor!
    CONTACT::Coupling3d coup(Discret(),Dim(),false,auxplane,sele,mele,testgps,testgpm,testjs,testji,printderiv);
  }
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");

  return true;
}

/*----------------------------------------------------------------------*
 | Split CElementsinto IntElements for 3D quad. coupling      popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::SplitIntElements(CONTACT::CElement& ele,
                       vector<RCP<CONTACT::IntElement> >& auxele)
{
  // *********************************************************************
  // do splitting for given element
  // *********************************************************** quad9 ***
  if (ele.Shape()==DRT::Element::quad9)
  {
    //dserror("ERROR: Quadratic 3D contact for quad9 under construction...");

    // split into for quad4 elements
    int numnode = 4;
    DRT::Element::DiscretizationType dt = DRT::Element::quad4;

    // first integration element
    // containing parent nodes 0,4,8,7
    int nodeids[4] = {0,0,0,0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[8];
    nodeids[3] = ele.NodeIds()[7];

    vector<DRT::Node*> nodes(4);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[4];
    nodes[2] = ele.Nodes()[8];
    nodes[3] = ele.Nodes()[7];

    auxele.push_back(rcp(new IntElement(0,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));

    // second integration element
    // containing parent nodes 4,1,5,8
    nodeids[0] = ele.NodeIds()[4];
    nodeids[1] = ele.NodeIds()[1];
    nodeids[2] = ele.NodeIds()[5];
    nodeids[3] = ele.NodeIds()[8];

    nodes[0] = ele.Nodes()[4];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[5];
    nodes[3] = ele.Nodes()[8];

    auxele.push_back(rcp(new IntElement(1,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));

    // third integration element
    // containing parent nodes 8,5,2,6
    nodeids[0] = ele.NodeIds()[8];
    nodeids[1] = ele.NodeIds()[5];
    nodeids[2] = ele.NodeIds()[2];
    nodeids[3] = ele.NodeIds()[6];

    nodes[0] = ele.Nodes()[8];
    nodes[1] = ele.Nodes()[5];
    nodes[2] = ele.Nodes()[2];
    nodes[3] = ele.Nodes()[6];

    auxele.push_back(rcp(new IntElement(2,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));

    // fourth integration element
    // containing parent nodes 7,8,6,3
    nodeids[0] = ele.NodeIds()[7];
    nodeids[1] = ele.NodeIds()[8];
    nodeids[2] = ele.NodeIds()[6];
    nodeids[3] = ele.NodeIds()[3];

    nodes[0] = ele.Nodes()[7];
    nodes[1] = ele.Nodes()[8];
    nodes[2] = ele.Nodes()[6];
    nodes[3] = ele.Nodes()[3];

    auxele.push_back(rcp(new IntElement(3,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));
  }

  // *********************************************************** quad8 ***
  else if (ele.Shape()==DRT::Element::quad8)
  {
    //dserror("ERROR: Quadratic 3D contact for quad8 under construction...");

    // split into four tri3 elements and one quad4 element
    int numnodetri = 3;
    int numnodequad = 4;
    DRT::Element::DiscretizationType dttri = DRT::Element::tri3;
    DRT::Element::DiscretizationType dtquad = DRT::Element::quad4;

    // first integration element
    // containing parent nodes 0,4,7
    int nodeids[3] = {0,0,0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[7];

    vector<DRT::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[4];
    nodes[2] = ele.Nodes()[7];

    auxele.push_back(rcp(new IntElement(0,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dttri,numnodetri,nodeids,nodes,ele.IsSlave())));

    // second integration element
    // containing parent nodes 1,5,4
    nodeids[0] = ele.NodeIds()[1];
    nodeids[1] = ele.NodeIds()[5];
    nodeids[2] = ele.NodeIds()[4];

    nodes[0] = ele.Nodes()[1];
    nodes[1] = ele.Nodes()[5];
    nodes[2] = ele.Nodes()[4];

    auxele.push_back(rcp(new IntElement(1,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dttri,numnodetri,nodeids,nodes,ele.IsSlave())));

    // third integration element
    // containing parent nodes 2,6,5
    nodeids[0] = ele.NodeIds()[2];
    nodeids[1] = ele.NodeIds()[6];
    nodeids[2] = ele.NodeIds()[5];

    nodes[0] = ele.Nodes()[2];
    nodes[1] = ele.Nodes()[6];
    nodes[2] = ele.Nodes()[5];

    auxele.push_back(rcp(new IntElement(2,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dttri,numnodetri,nodeids,nodes,ele.IsSlave())));

    // fourth integration element
    // containing parent nodes 3,7,6
    nodeids[0] = ele.NodeIds()[3];
    nodeids[1] = ele.NodeIds()[7];
    nodeids[2] = ele.NodeIds()[6];

    nodes[0] = ele.Nodes()[3];
    nodes[1] = ele.Nodes()[7];
    nodes[2] = ele.Nodes()[6];

    auxele.push_back(rcp(new IntElement(3,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dttri,numnodetri,nodeids,nodes,ele.IsSlave())));

    // fifth integration element
    // containing parent nodes 4,5,6,7
    int nodeidsquad[4] = {0,0,0,0};
    nodeidsquad[0] = ele.NodeIds()[4];
    nodeidsquad[1] = ele.NodeIds()[5];
    nodeidsquad[2] = ele.NodeIds()[6];
    nodeidsquad[3] = ele.NodeIds()[7];

    vector<DRT::Node*> nodesquad(4);
    nodesquad[0] = ele.Nodes()[4];
    nodesquad[1] = ele.Nodes()[5];
    nodesquad[2] = ele.Nodes()[6];
    nodesquad[3] = ele.Nodes()[7];

    auxele.push_back(rcp(new IntElement(4,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dtquad,numnodequad,nodeidsquad,nodesquad,ele.IsSlave())));
  }

  // ************************************************************ tri6 ***
  else if (ele.Shape()==DRT::Element::tri6)
  {
    //dserror("ERROR: Quadratic 3D contact for tri6 under construction...");

    // split into for tri3 elements
    int numnode = 3;
    DRT::Element::DiscretizationType dt = DRT::Element::tri3;

    // first integration element
    // containing parent nodes 0,3,5
    int nodeids[3] = {0,0,0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[3];
    nodeids[2] = ele.NodeIds()[5];

    vector<DRT::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[3];
    nodes[2] = ele.Nodes()[5];

    auxele.push_back(rcp(new IntElement(0,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));

    // second integration element
    // containing parent nodes 3,1,4
    nodeids[0] = ele.NodeIds()[3];
    nodeids[1] = ele.NodeIds()[1];
    nodeids[2] = ele.NodeIds()[4];

    nodes[0] = ele.Nodes()[3];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[4];

    auxele.push_back(rcp(new IntElement(1,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));

    // third integration element
    // containing parent nodes 5,4,2
    nodeids[0] = ele.NodeIds()[5];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[2];

    nodes[0] = ele.Nodes()[5];
    nodes[1] = ele.Nodes()[4];
    nodes[2] = ele.Nodes()[2];

    auxele.push_back(rcp(new IntElement(2,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));

    // fourth integration element
    // containing parent nodes 4,5,3
    nodeids[0] = ele.NodeIds()[4];
    nodeids[1] = ele.NodeIds()[5];
    nodeids[2] = ele.NodeIds()[3];

    nodes[0] = ele.Nodes()[4];
    nodes[1] = ele.Nodes()[5];
    nodes[2] = ele.Nodes()[3];

    auxele.push_back(rcp(new IntElement(3,ele.Id(),ele.Type(),ele.Owner(),
        ele.Shape(),dt,numnode,nodeids,nodes,ele.IsSlave())));
  }

  // *********************************************************** quad4 ***
  else if (ele.Shape()==DRT::Element::quad4)
  {
    // 1:1 conversion to IntElement
    vector<DRT::Node*> nodes(4);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[2];
    nodes[3] = ele.Nodes()[3];

    auxele.push_back(rcp(new IntElement(0,ele.Id(),ele.Type(),ele.Owner(),
       ele.Shape(),ele.Shape(),ele.NumNode(),ele.NodeIds(),nodes,ele.IsSlave())));
  }

  // ************************************************************ tri3 ***
  else if (ele.Shape()==DRT::Element::tri3)
  {
    // 1:1 conversion to IntElement
    vector<DRT::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[2];

    auxele.push_back(rcp(new IntElement(0,ele.Id(),ele.Type(),ele.Owner(),
       ele.Shape(),ele.Shape(),ele.NumNode(),ele.NodeIds(),nodes,ele.IsSlave())));
  }

  // ********************************************************* invalid ***
  else
    dserror("ERROR: SplitIntElements called for unknown element shape!");

  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate penalty scaling factor kapp (public)            popp 11/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateKappaPenalty(CONTACT::CElement& sele)
{
  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // check for quadratic 3D case
  bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

  // ************************************************** quadratic 3D ***
  if (Dim()==3 && sele.IsQuad())
  {
    // only for auxiliary plane 3D version
    if (!auxplane) dserror("ERROR: Quadratic 3D contact only for AuxPlane case!");

    // build linear integration elements from quadratic CElements
    vector<RCP<CONTACT::IntElement> > sauxelements(0);
    SplitIntElements(sele,sauxelements);

    for (int i=0;i<(int)sauxelements.size();++i)
    {
      // do the int element integration of kappa and store into gap
#ifdef CONTACTPETROVGALERKIN
      int nrow = sauxelements[i]->NumNode();
#else
      int nrow = sele.NumNode();
#endif // #ifdef CONTACTPETROVGALERKIN
      RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

      // create an integrator instance with correct NumGP and Dim
      CONTACT::Integrator integrator(shapefcn_,sauxelements[i]->Shape());
      integrator.IntegrateKappaPenalty(sele,*(sauxelements[i]),sxia,sxib,gseg);

      // do the assembly into the slave nodes
#ifdef CONTACTPETROVGALERKIN
      integrator.AssembleG(Comm(),*(sauxelements[i]),*gseg);
#else
      integrator.AssembleG(Comm(),sele,*gseg);
#endif //#ifdef CONTACTPETROVGALERKIN
    }
  }

  else
  {
    // do the element integration of kappa and store into gap
    int nrow = sele.NumNode();
    RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));

    // create an integrator instance with correct NumGP and Dim
    CONTACT::Integrator integrator(shapefcn_,sele.Shape());
    integrator.IntegrateKappaPenalty(sele,sxia,sxib,gseg);

    // do the assembly into the slave nodes
    integrator.AssembleG(Comm(),sele,*gseg);
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate relative movement (jump) of a slave node      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateRelMov()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // parameters
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    // get some informatiom form the node
    double gap = cnode->Getg();
    int dim = cnode->NumDof();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k=0;k<3;++k)
      nz += cnode->n()[k] * cnode->lm()[k];

    vector <double> jump(dim);
    for(int dim=0;dim<Dim();dim++)
      jump[dim] = 0;

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->lmuzawa()[k]*cnode->n()[k];

    double kappa = cnode->Kappa();

    // evaluate jump (relative displacement) of this node
    // only when the node is going to be active, otherwise,
    // this value isn't needed.
    bool activeinfuture = false;

    if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_penalty)
    {
    	if (-gap >= 0) activeinfuture = true;
    }
    else if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
    		     Teuchos::getIntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")!=1)
    {
    	if (-gap >= 0) activeinfuture = true;
    }
    else if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
		         Teuchos::getIntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")==1)
    {
    	if(nz - cn*gap > 0) activeinfuture = true;
    }
    else if (Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_auglag)
    {
    	if(lmuzawan - kappa * pp * gap >= 0) activeinfuture = true;
    }
    else
    	dserror("Error in Interface::EvaluateRelMov(): Solution strategy not known!");

    if(activeinfuture==true)
    {
      vector<map<int,double> > dmap = cnode->GetD();
      vector<map<int,double> > dmapold = cnode->GetDOld();

      // check if there are entries in the old D map
      if(dmapold.size()< 1)
      	dserror("Error in Interface::EvaluateRelMov(): No old D-Map!");

      map<int,double>::iterator colcurr;
      set <int> snodes = cnode->GetSNodes();

      set<int>::iterator scurr;

      // loop over all slave nodes with an entry adjacent to this node
      for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* csnode = static_cast<CNode*>(snode);
        const int* sdofs = csnode->Dofs();

        double dik = (dmap[0])[sdofs[0]];
        double dikold = (dmapold[0])[sdofs[0]];

        map<int,double>::iterator mcurr;

        for (int dim=0;dim<csnode->NumDof();++dim)
        {
            jump[dim]-=(dik-dikold)*(csnode->xspatial()[dim]);
        }
      } //  loop over adjacent slave nodes

      vector<map<int,double> > mmap = cnode->GetM();
      vector<map<int,double> > mmapold = cnode->GetMOld();

      // check if there are entries in the old M map
      if(mmapold.size()< 1)
      	dserror("Error in Interface::EvaluateRelMov(): No old M-Map!");

      set <int> mnodes;
    	set<int>::iterator mcurr;

    	set <int> mnodescurrent = cnode->GetMNodes();
    	set <int> mnodesold = cnode->GetMNodesOld();

    	for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
    	  mnodes.insert(*mcurr);

    	for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
    	  mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        map<int,double>::iterator mcurr;

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
           jump[dim]+= (mik-mikold)*(cmnode->xspatial()[dim]);
        }
      } //  loop over master nodes

      // write it to nodes
      for(int dim=0;dim<Dim();dim++)
        cnode->jump()[dim] = jump[dim];


      // linearization of jump vector
      /*** 01  **********************************************************/
      // loop over all slave nodes (find adjacent ones to this stick node)
      for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* csnode = static_cast<CNode*>(snode);
        const int* sdofs = csnode->Dofs();

        double dik = (dmap[0])[sdofs[0]];
        double dikold=(dmapold[0])[sdofs[0]];

   	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
  	    {
          int col = csnode->Dofs()[dimrow];
          double val = -(dik-dikold);
          cnode->AddDerivJumpValue(dimrow,col,val);
  	    }
      }

      /*** 02  **********************************************************/
      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold=(mmapold[0])[mdofs[0]];

  	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
  	    {
            int col = cmnode->Dofs()[dimrow];
            double val = (mik-mikold);
            cnode->AddDerivJumpValue(dimrow,col,val);
  	    }
      }

      /*** 03 ***********************************************************/
      // we need the Lin(D-matrix) entries of this node
      map<int,map<int,double> >& ddmap = cnode->GetDerivD();
      map<int,map<int,double> >::iterator dscurr;

      // loop over all slave nodes in the DerivM-map of the stick slave node
      for (dscurr=ddmap.begin();dscurr!=ddmap.end();++dscurr)
      {
        int gid = dscurr->first;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* csnode = static_cast<CNode*>(snode);
        double* dxi = csnode->xspatial();

        // compute entry of the current stick node / slave node pair
        map<int,double>& thisdmmap = cnode->GetDerivD(gid);

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
          	double val =-colcurr->second*dxi[dimrow];
          	cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }

      /*** 04 ***********************************************************/
      // we need the Lin(M-matrix) entries of this node
      map<int,map<int,double> >& dmmap = cnode->GetDerivM();
      map<int,map<int,double> >::iterator dmcurr;

      // loop over all master nodes in the DerivM-map of the stick slave node
      for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
      {
        int gid = dmcurr->first;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        double* mxi = cmnode->xspatial();

        // compute entry of the current stick node / master node pair
        map<int,double>& thisdmmap = cnode->GetDerivM(gid);

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
          	double val =colcurr->second*mxi[dimrow];
          	cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }
    } // active nodes
  } // loop over slave nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble relative movement / jump (global)             gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleRelMov(Epetra_Vector& jumpglobal)
{
  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    int dim = cnode->NumDof();
    double* jump = cnode->jump();

    Epetra_SerialDenseVector jumpnode(dim);
    vector<int> jumpdof(dim);
    vector<int> jumpowner(dim);

    for( int k=0; k<dim; ++k )
    {
      jumpnode(k) = jump[k];
      jumpdof[k] = cnode->Dofs()[k];
      jumpowner[k] = cnode->Owner();
    }

    // do assembly
    LINALG::Assemble(jumpglobal, jumpnode, jumpdof, jumpowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate L2 Norm of tangential contact conditions     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateTangentNorm(double& cnormtan)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  cnormtan=0;

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    // get some information form node
    double* n = cnode->n();
    int dim = cnode->NumDof();

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // jump vector
    Epetra_SerialDenseMatrix jumpvec(dim,1);
    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->jump()[i];

    // jump vector
    Epetra_SerialDenseMatrix forcevec(dim,1);
    for (int i=0;i<dim;i++)
      forcevec(i,0) = cnode->lm()[i];

    // evaluate jump in tangential direction
    Epetra_SerialDenseMatrix jumptan(dim,1);
    jumptan.Multiply('N','N',1,tanplane,jumpvec,0.0);

    // norm of tangential jumps for stick nodes
    if (cnode->Active()== true and cnode->Slip()==false)
    {
     	for( int j=0;j<cnode->NumDof();++j)
        cnormtan+=jumptan(j,0)*jumptan(j,0);
    }
    else if (cnode->Active()== true and cnode->Slip()==true)
    {
      double jumptxi = 0;
      double jumpteta = 0;
      double forcen = 0;
      double forcetxi = 0;
      double forceteta = 0;

      for (int i=0;i<dim;i++)
      {
      	jumptxi+=cnode->txi()[i]*cnode->jump()[i];
      	jumpteta+=cnode->teta()[i]*cnode->jump()[i];

      	forcen+=cnode->n()[i]*cnode->lm()[i];
      	forcetxi+=cnode->txi()[i]*cnode->lm()[i];
      	forceteta+=cnode->teta()[i]*cnode->lm()[i];
      }

      //cout << "FACTOR-Direction " << (jumptxi/jumpteta)/(forcetxi/forceteta) << endl;
      //cout << "FACTOR-Magnitude" << (frcoeff*forcen)/(sqrt(forcetxi*forcetxi+forceteta*forceteta))<< endl;
    }
  } // loop over slave nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                 popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleRegNormalForces(bool& localisincontact,
                                                 bool& localactivesetchange)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    int dim = cnode->NumDof();
    double gap = cnode->Getg();

    // modified gap for zero initial gap
    // (if gap is below zero, it is explicitly set to zero)
    //double modgap = cnode->Getg();
    //if (abs(modgap) < 1.0e-10) modgap=0.0;

    double kappa = cnode->Kappa();

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->lmuzawa()[k]*cnode->n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitely to constant
    // and corresponding derivatives to zero

    for( int j=0;j<dim;++j)
      cnode->lm()[j] = i*j;

    cnode->GetDerivZ().clear();

    continue;
#endif

    //********************************************************************
    // Decision on active /  inactive nodes (regularization)
    //
    // CASE 1: Penalty approach
    // A node is activated if its weighted gap is negative or deactivated
    // if its gap is equal zero or positive.
    // -> the regularization reads: lambda_n = kappa * pp * < -gap >
    //
    // CASE 2: Augmented Lagrange approach
    // A node is activated if its Lagrange multiplier, stemming from the
    // last Uzawa Lagrange multiplier AND the current regularization is
    // negative or deactivated if its LM is equal zero or positive.
    // -> the regularization reads: lambda_n = < lmuzawa_n - kappa * pp * gap >
    //
    // As the Uzawa Lagrange multipliers are zero in the penalty approach,
    // the two cases can formally be treted identically, see below.
    // We do not need an explicit separation of cases!
    //
    //********************************************************************

    // Activate/Deactivate node and notice any change
    if( (cnode->Active() == false) && (lmuzawan - kappa * pp * gap >= 0) )
    {
        cnode->Active() = true;
        localactivesetchange = true;

        //cout << "node #" << gid << " is now active (";
        //for( int j=0; j<dim; j++)
        //  cout << " " << cnode->Dofs()[j] << " ";
        //cout << ") gap=" << gap << endl;
    }

    else if( (cnode->Active() == true) && (lmuzawan - kappa * pp * gap < 0) )
    {
        cnode->Active() = false;
        localactivesetchange = true;

        //cout << "node #" << gid << " is now inactive, gap=" << gap << endl;
    }
    //********************************************************************

    // Compute derivZ-entries with the Macauley-Bracket
    // of course, this is only done for active constraints in order
    // for linearization and r.h.s to match!
    if( cnode->Active()==true )
    {

//  	  cout << "GID " << gid << endl;
//  	  cout << "LMUZAWAN " << lmuzawan << endl;
//  	  cout << "GAP " << gap << endl;

      localisincontact = true;

      double* normal = cnode->n();

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->lm()[j] = (lmuzawan - kappa * pp * gap) * normal[j];

      // compute derivatives of lagrange multipliers and store into node

      // contribution of derivative of weighted gap
      map<int,double>& derivg = cnode->GetDerivG();
      map<int,double>::iterator gcurr;

      // contribution of derivative of normal
      vector<map<int,double> >& derivn = cnode->GetDerivN();
      map<int,double>::iterator ncurr;

      for( int j=0;j<dim;++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
          cnode->AddDerivZValue(j, gcurr->first, - kappa * pp * (gcurr->second) * normal[j]);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, - kappa * pp * gap * ncurr->second);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, + lmuzawan * ncurr->second);
      }
    }

    // be sure to remove all LM-related stuff from inactive nodes
    else
    {
      // clear lagrange multipliers
      for( int j=0;j<dim;++j) cnode->lm()[j] = 0;

      // clear derivz
      cnode->GetDerivZ().clear();

    } // Macauley-Bracket
  } // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                 gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleRegTangentForcesPenalty()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    // get some informatiom form the node
    double gap = cnode->Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->Kappa();
    double* n = cnode->n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->lmuzawa()[k]*cnode->n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);
       tanplane(0,2)=  -(n[0]*n[2]);
       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
       tanplane(1,2)=  -(n[1]*n[2]);

       tanplane(2,0)=  -(n[2]*n[0]);
       tanplane(2,1)=  -(n[2]*n[1]);
       tanplane(2,2)= 1-(n[2]*n[2]);
     }
     else if (dim==2)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);

       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
     }
     else
       dserror("Error in AssembleTangentForces: Unknown dimension.");

    // Lagrange multiplier in tangential direction
    Epetra_SerialDenseMatrix lmuzawatan(dim,1);
    lmuzawatan.Multiply('N','N',1,tanplane,lmuzawa,0.0);

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // fill vector tractionold
    vector<double> tractionold(dim);
    for (int i=0;i<dim;i++)
    	tractionold[i] = cnode->tractionold()[i];

    // Evaluate trailtraction (tractionold+temptrac in penalty case)
    vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
    	trailtraction[i]=tractionold[i]+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))
    {
    	//cout << "Node " << gid << " is stick" << endl;
    	cnode->Slip() = false;

    	// in the stick case, traction is trailtraction
    	for (int i=0;i<dim;i++)
    		cnode->traction()[i]=trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->lm()[j] = n[j]*(- kappa * ppnor * gap) + trailtraction[j];
    }
    else
    {
    	//cout << "Node " << gid << " is slip" << endl;
    	cnode->Slip() = true;

    	// in the slip case, traction is evaluated with a return map algorithm
    	for (int i=0;i<dim;i++)
    		cnode->traction()[i]=maxtantrac/magnitude*trailtraction[i];

    	// compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->lm()[j] = n[j]*(- kappa * ppnor * gap)+maxtantrac/magnitude*trailtraction[j];
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->Slip() == false)
    {
	    /***************************************** tanplane.deriv(jump) ***/
	    vector<map<int,double> >& derivjump = cnode->GetDerivJump();
	    map<int,double>::iterator colcurr;

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	    	for (int dim=0;dim<cnode->NumDof();++dim)
	    	{
	        // loop over all entries of the current derivative map
	        for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
	        {
	      		  int col = colcurr->first;
	            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
	            cnode->AddDerivZValue(dimrow,col,val);
	        }
	    	}
	    }

	    /**************************************** deriv(tanplane).jump  ***/
	    vector<map<int,double> >& derivn = cnode->GetDerivN();

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
	      {
	      	for (int dim =0;dim<cnode->NumDof();++dim)
	      	{
	      		int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->jump()[dim];
	          cnode->AddDerivZValue(dimrow,col,val);
	      	}
	      }
	    }

	    // loop over dimensions
	    for (int dim=0;dim<cnode->NumDof();++dim)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
	      {
	      	for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
	      	{
	      		int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->jump()[dim];
	          cnode->AddDerivZValue(dimrow,col,val);
	      	}
	      }
      }
    }
    // slip nodes
    else if (cnode->Active() == true && cnode->Slip()== true)
    {
	    /******************** tanplane.deriv(jump).maxtantrac/magnidude ***/

    	vector<map<int,double> >& derivjump = cnode->GetDerivJump();
	    map<int,double>::iterator colcurr;

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	    	for (int dim=0;dim<cnode->NumDof();++dim)
	    	{
	        // loop over all entries of the current derivative map
	        for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
	        {
	      		  int col = colcurr->first;
	            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
	      		  cnode->AddDerivZValue(dimrow,col,val);
	        }
	    	}
	    }

	    /******************** deriv(tanplane).jump.maxtantrac/magnitude ***/
	    vector<map<int,double> >& derivn = cnode->GetDerivN();

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
	      {
	      	for (int dim =0;dim<cnode->NumDof();++dim)
	      	{
	      		int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->jump()[dim]*maxtantrac/magnitude;
	          cnode->AddDerivZValue(dimrow,col,val);
	      	}
	      }
	    }
	    // loop over dimensions
	    for (int dim=0;dim<cnode->NumDof();++dim)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
	      {
	      	for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
	      	{
	      		int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->jump()[dim]*maxtantrac/magnitude;
	      		cnode->AddDerivZValue(dimrow,col,val);
	      	}
	      }
      }

	    /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      map<int,double>& derivg = cnode->GetDerivG();
      map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j, gcurr->first, - frcoeff*kappa * ppnor * (gcurr->second)*trailtraction[j]/magnitude);
        }
      }

	    /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      vector <double> temp(cnode->NumDof());
 	    for (int dim=0;dim<cnode->NumDof();++dim)
	      temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
	    for (int dimout=0;dimout<cnode->NumDof();++dimout)
	    {
	    	double traction = 0;
	    	for (int dim=0;dim<cnode->NumDof();++dim)
	    	  traction += tanplane(dimout,dim)*cnode->jump()[dim]*kappa*pptan;

	    	traction+= tractionold[dimout];

	 	    for (int dim=0;dim<cnode->NumDof();++dim)
	 	    {
		      // loop over all entries of the current derivative map
		      for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
		      {
		        int col = colcurr->first;
		      	double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

		        for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
		      	{
  	      	  double val1 = val*temp[dimrow];
		      		cnode->AddDerivZValue(dimrow,col,val1);
  	      	}
		      }
		    }
		  }

      // loop over dimensions
		  for (int dimout=0;dimout<cnode->NumDof();++dimout)
		  {
		    double traction = 0;
		    for (int dim=0;dim<cnode->NumDof();++dim)
		      traction += tanplane(dimout,dim)*cnode->jump()[dim]*kappa*pptan;

		    traction+=tractionold[dimout];

		    // loop over all entries of the current derivative map
  	    for (colcurr=derivn[dimout].begin();colcurr!=derivn[dimout].end();++colcurr)
	      {
		      int col = colcurr->first;

		      for (int dim=0;dim<cnode->NumDof();++dim)
		      {
		      	double val =-colcurr->second*n[dim]*cnode->jump()[dim]*traction/magnitude*pptan*kappa;
	        	for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	        	{
	         	  double val1 = val*temp[dimrow];
			     		cnode->AddDerivZValue(dimrow,col,val1);
	         	}
		      }
		    }
			}

		  // loop over dimensions
			for (int dimout=0;dimout<cnode->NumDof();++dimout)
			{
			  double traction = 0;
			  for (int dim=0;dim<cnode->NumDof();++dim)
			    traction += tanplane(dimout,dim)*cnode->jump()[dim]*kappa*pptan;

			    traction += tractionold[dimout];

		      for (int dim=0;dim<cnode->NumDof();++dim)
		      {

			    // loop over all entries of the current derivative map
	  	    for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
		      {
			      int col = colcurr->first;

			     	double val =-colcurr->second*n[dimout]*cnode->jump()[dim]*traction/magnitude*pptan*kappa;

			     	for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
				    {
		  	       double val1 = val*temp[dimrow];
				       cnode->AddDerivZValue(dimrow,col,val1);
 		      	}
          }
        }
	    }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->lm()[j] = 0;
      // clear derivz
      cnode->GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces (Aug. Lagr.)    gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleRegTangentForcesAugmented()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    // get some informatiom form the node
    double gap = cnode->Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->Kappa();
    double* n = cnode->n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->lmuzawa()[k]*cnode->n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // Lagrange multiplier in tangential direction
    Epetra_SerialDenseMatrix lmuzawatan(dim,1);
    lmuzawatan.Multiply('N','N',1,tanplane,lmuzawa,0.0);

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // Evaluate trailtraction
    vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
    	trailtraction[i]=lmuzawatan(i,0)+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))    {
    	//cout << "Node " << gid << " is stick" << endl;
    	cnode->Slip() = false;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j];
    }
    else
    {
    	//cout << "Node " << gid << " is slip" << endl;
    	cnode->Slip() = true;

     	// compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j]*maxtantrac/magnitude;
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->Slip() == false)
    {
	    /***************************************** tanplane.deriv(jump) ***/
	    vector<map<int,double> >& derivjump = cnode->GetDerivJump();
	    map<int,double>::iterator colcurr;

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	    	for (int dim=0;dim<cnode->NumDof();++dim)
	    	{
	        // loop over all entries of the current derivative map
	        for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
	        {
	      	  int col = colcurr->first;
	          double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
	          cnode->AddDerivZValue(dimrow,col,val);
	        }
	    	}
	    }

	    /******************************* deriv(tanplane).(lmuzawa+jump) ***/
	    vector<map<int,double> >& derivn = cnode->GetDerivN();

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
	      {
	      	for (int dim =0;dim<cnode->NumDof();++dim)
	      	{
	      		int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dim]*(cnode->jump()[dim]);
	          val = val - (colcurr->second)*n[dim]*(cnode->lmuzawa()[dim]);
	          cnode->AddDerivZValue(dimrow,col,val);
	      	}
	      }
	    }

	    // loop over dimensions
	    for (int dim=0;dim<cnode->NumDof();++dim)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
	      {
	      	for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
	      	{
	      		int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dimrow]*(cnode->jump()[dim]);
	          val = val-(colcurr->second)*n[dimrow]*(cnode->lmuzawa()[dim]);
	          cnode->AddDerivZValue(dimrow,col,val);
	      	}
	      }
      }
    }

    // slip nodes
    else if (cnode->Active() == true && cnode->Slip()== true)
    {
    	/***************************************** tanplane.deriv(jump) ***/
    	vector<map<int,double> >& derivjump = cnode->GetDerivJump();
	    map<int,double>::iterator colcurr;

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	    	for (int dim=0;dim<cnode->NumDof();++dim)
	    	{
	        // loop over all entries of the current derivative map
	        for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
	        {
	     		  int col = colcurr->first;
	          double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
	 	        cnode->AddDerivZValue(dimrow,col,val);
	        }
	    	}
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      vector<map<int,double> >& derivn = cnode->GetDerivN();

	    // loop over dimensions
	    for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dimrow].begin();colcurr!=derivn[dimrow].end();++colcurr)
	      {
	        for (int dim =0;dim<cnode->NumDof();++dim)
	        {
	          int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dim]*cnode->jump()[dim];
	          val = (val - (colcurr->second)*n[dim]*(cnode->lmuzawa()[dim]))*maxtantrac/magnitude;
	          cnode->AddDerivZValue(dimrow,col,val);
	        }
	      }
	    }

	    // loop over dimensions
	    for (int dim=0;dim<cnode->NumDof();++dim)
	    {
	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
	      {
	        for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
	        {
	          int col = colcurr->first;
	          double val =-pptan*kappa*(colcurr->second)*n[dimrow]*cnode->jump()[dim];
	          val = (val-(colcurr->second)*n[dimrow]*(cnode->lmuzawa()[dim]))*maxtantrac/magnitude;
	          cnode->AddDerivZValue(dimrow,col,val);
	        }
	      }
	    }

	    /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      map<int,double>& derivg = cnode->GetDerivG();
      map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
        	cnode->AddDerivZValue(j,gcurr->first,- frcoeff*kappa*ppnor*(gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( colcurr = (derivn[j]).begin(); colcurr != (derivn[j]).end(); ++colcurr )
        {
        	for( int k=0;k<cnode->NumDof();++k)
        	{
        		double val = frcoeff*(colcurr->second)*lmuzawa(j,0)*trailtraction[k]/magnitude;
        		cnode->AddDerivZValue(k,colcurr->first,val);
        	}
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      vector <double> temp(cnode->NumDof());
	      for (int dim=0;dim<cnode->NumDof();++dim)
          temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
      	double traction = 0;
      	for (int dim=0;dim<cnode->NumDof();++dim)
      	  traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
 	      {
	        // loop over all entries of the current derivative map
	        for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
	        {
	          int col = colcurr->first;
	      	  double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

	          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
	      	  {
	      	    double val1 = val*temp[dimrow];
	      	    cnode->AddDerivZValue(dimrow,col,val1);
	        	}
	        }
	      }
	    }

      // loop over dimensions
	    for (int dimout=0;dimout<cnode->NumDof();++dimout)
	    {
	      double traction = 0;
	      for (int dim=0;dim<cnode->NumDof();++dim)
    	    traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->jump()[dim]*kappa*pptan);

	      // loop over all entries of the current derivative map
	      for (colcurr=derivn[dimout].begin();colcurr!=derivn[dimout].end();++colcurr)
        {
	        int col = colcurr->first;

	        for (int dim=0;dim<cnode->NumDof();++dim)
	        {
	        	double val =-colcurr->second*n[dim]*(lmuzawa(dim,0)+cnode->jump()[dim]*pptan*kappa)*traction/magnitude;
	        	for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          	{
         	    double val1 = val*temp[dimrow];
        	  	cnode->AddDerivZValue(dimrow,col,val1);
         	  }
	        }
	      }
  		}

	    // loop over dimensions
		  for (int dimout=0;dimout<cnode->NumDof();++dimout)
		  {
		    double traction = 0;
		    for (int dim=0;dim<cnode->NumDof();++dim)
    	    traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->jump()[dim]*kappa*pptan);

	      for (int dim=0;dim<cnode->NumDof();++dim)
	      {
  		    // loop over all entries of the current derivative map
  	      for (colcurr=derivn[dim].begin();colcurr!=derivn[dim].end();++colcurr)
	        {
		        int col = colcurr->first;
		     	  double val =-colcurr->second*n[dimout]*(lmuzawa(dim,0)+cnode->jump()[dim]*pptan*kappa)*traction/magnitude;

		      	for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
			      {
	  	        double val1 = val*temp[dimrow];
		     	  	cnode->AddDerivZValue(dimrow,col,val1);
		        }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->lm()[j] = 0;
      // clear derivz
      cnode->GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble geometry-dependent lagrange multipliers (global)      popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleLM(Epetra_Vector& zglobal)
{
  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    int dim = cnode->NumDof();
    double* lm = cnode->lm();

    Epetra_SerialDenseVector lmnode(dim);
    vector<int> lmdof(dim);
    vector<int> lmowner(dim);

    for( int k=0; k<dim; ++k )
    {
      lmnode(k) = lm[k];
      lmdof[k] = cnode->Dofs()[k];
      lmowner[k] = cnode->Owner();
    }

    // do assembly
    LINALG::Assemble(zglobal, lmnode, lmdof, lmowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble derivatives of lagrange multipliers              popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleLinZ(LINALG::SparseMatrix& linzglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all slave nodes (row map)
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinZ: Node ownership inconsistency!");

    // derivz is the vector<map> we want to assemble
    vector<map<int,double> >& derivz = cnode->GetDerivZ();

    if ( (int) derivz.size()>0 )
    {
      int rowsize = cnode->NumDof();
      int colsize = (int) derivz[0].size();

      // consistency check
      for (int j=0; j<rowsize-1; ++j)
        if ((int)derivz[j].size() != (int)derivz[j+1].size())
          dserror("ERROR: AssembleLinZ: Column dim. of nodal derivz-map is inconsistent!");

      map<int,double>::iterator colcurr;

      // loop over dofs
      for ( int k=0; k<rowsize; ++k )
      {
        int row = cnode->Dofs()[k]; // row index equals global dof index of this #i node's dof k
        int l = 0;

        // loop over all directional derivative entries using the map iterator
        for( colcurr = derivz[k].begin(); colcurr != derivz[k].end(); ++colcurr )
        {
          int col = colcurr->first; // col index equals global id of directional derivative component ,l
          double val = colcurr->second;
          linzglobal.Assemble(val, row, col);
          l++;
        }

        if( l != colsize )
          dserror("ERROR: AssembleLinZ: l = %i but colsize = %i",k,colsize);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar matrices and weighted gap                 popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleDMG(LINALG::SparseMatrix& dglobal,
                                      LINALG::SparseMatrix& mglobal,
                                      Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDMG: Node ownership inconsistency!");

    /**************************************************** D-matrix ******/
    if ((cnode->GetD()).size()>0)
    {
      vector<map<int,double> > dmap = cnode->GetD();
      int rowsize = cnode->NumDof();
      int colsize = (int)dmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)dmap[j].size() != (int)dmap[j+1].size())
          dserror("ERROR: AssembleDMG: Column dim. of nodal D-map is inconsistent!");

      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = cnode->Dofs()[j];
        int k = 0;

        for (colcurr=dmap[j].begin();colcurr!=dmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;

          // do the assembly into global D matrix
          if (shapefcn_ == Interface::DualFunctions)
          {
            // check for diagonality
            if (row!=col && abs(val)>1.0e-12)
              dserror("ERROR: AssembleDMG: D-Matrix is not diagonal!");

            // create an explicitly diagonal d matrix
            if (row==col)
              dglobal.Assemble(val, row, col);
          }
          else if (shapefcn_ == Interface::StandardFunctions)
          {
            // don't check for diagonality
            // since for standard shape functions, as in general when using
            // arbitrary shape function types, this is not the case

            // create the d matrix, do not assemble zeros
            dglobal.Assemble(val, row, col);
          }

          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleDMG: k = %i but colsize = %i",k,colsize);
      }
    }

    /**************************************************** M-matrix ******/
    if ((cnode->GetM()).size()>0)
    {
      vector<map<int,double> > mmap = cnode->GetM();
      int rowsize = cnode->NumDof();
      int colsize = (int)mmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)mmap[j].size() != (int)mmap[j+1].size())
          dserror("ERROR: AssembleDMG: Column dim. of nodal M-map is inconsistent!");

      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = cnode->Dofs()[j];
        int k = 0;

        for (colcurr=mmap[j].begin();colcurr!=mmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;

          // do not assemble zeros into m matrix
          if (abs(val)>1.0e-12) mglobal.Assemble(val,row,col);
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleDMG: k = %i but colsize = %i",k,colsize);
      }
    }

    /************************************************* Mmod-matrix ******/
    if ((cnode->GetMmod()).size()>0)
    {
      vector<map<int,double> > mmap = cnode->GetMmod();
      int rowsize = cnode->NumDof();
      int colsize = (int)mmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)mmap[j].size() != (int)mmap[j+1].size())
          dserror("ERROR: AssembleDMG: Column dim. of nodal Mmod-map is inconsistent!");

      Epetra_SerialDenseMatrix Mnode(rowsize,colsize);
      vector<int> lmrow(rowsize);
      vector<int> lmcol(colsize);
      vector<int> lmrowowner(rowsize);
      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = cnode->Dofs()[j];
        int k = 0;
        lmrow[j] = row;
        lmrowowner[j] = cnode->Owner();

        for (colcurr=mmap[j].begin();colcurr!=mmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;
          lmcol[k] = col;

          Mnode(j,k)=val;
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleDMG: k = %i but colsize = %i",k,colsize);
      }

      mglobal.Assemble(-1,Mnode,lmrow,lmrowowner,lmcol);
    }

    /**************************************************** g-vector ******/
    if (cnode->Getg()!=0.0)
    {
      double gap = cnode->Getg();

      /*
      // replace integrated gap by definition of gap
      if (cnode->Active())
      {
        // check two versions of weighted gap
        double defgap = 0.0;
        double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];

        for (int j=0;j<3;++j)
          defgap-= (cnode->n()[j])*wii*(cnode->xspatial()[j]);

        vector<map<int,double> > mmap = cnode->GetM();
        map<int,double>::iterator mcurr;

        for (int m=0;m<mnodefullmap_->NumMyElements();++m)
        {
          int gid = mnodefullmap_->GID(m);
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CNode* cmnode = static_cast<CNode*>(mnode);
          const int* mdofs = cmnode->Dofs();
          bool hasentry = false;

          // look for this master node in M-map of the active slave node
          for (mcurr=mmap[0].begin();mcurr!=mmap[0].end();++mcurr)
            if ((mcurr->first)==mdofs[0])
            {
              hasentry=true;
              break;
            }

          double mik = (mmap[0])[mdofs[0]];
          double* mxi = cmnode->xspatial();

          // get out of here, if master node not adjacent or coupling very weak
          if (!hasentry || abs(mik)<1.0e-12) continue;

          for (int j=0;j<3;++j)
            defgap+= (cnode->n()[j]) * mik * mxi[j];
        }

        //cout << "SNode: " << cnode->Id() << " IntGap: " << gap << " DefGap: " << defgap << endl;
        gap = defgap;
      }
      */

      // cout << "Node ID: " << cnode->Id() << " HasProj: " << cnode->HasProj()
      //      << " IsActive: " << cnode->Active() << " Gap: " << gap << endl;

      // check if this active node has a feasible projection
      // else, one would at first think that a dserror has to be thrown!
      // (but this is not true in general, as there might indeed be an
      // active node which nervertheless has no feasible projection,
      // e.g. a slave node which is just over the edge of the master surface)
      // -> thus this check has been removed (popp 03/09)

      //if (!cnode->HasProj() && cnode->Active())
      //  dserror("ERROR: Active node ID: %i without feasible projection", cnode->Id());

      // check if this inactive node has a feasible projection
      // else, it cannot be in contact and weighted gap should be positive
      // (otherwise wrong results possible for g~ because of non-positivity
      // of dual shape functions!!!)
      // when applying a Petrov-Galerkin scheme with standard shape functions
      // for the weighted gap interpolation this problem does not exist
      // FIXME: Only the linear case considered here (popp 04/09)
#ifndef CONTACTPETROVGALERKIN
      if (!cnode->HasProj() && !cnode->Active())
      {
      	gap = 1.0e12;
      	cnode->Getg()=gap;
      }
#endif // #ifndef CONTACTPETROVGALERKIN

      Epetra_SerialDenseVector gnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);

      gnode(0) = gap;
      lm[0] = cnode->Id();
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(gglobal,gnode,lm,lmowner);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices with nodal normals / tangents           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleNT(LINALG::SparseMatrix& nglobal,
                                     LINALG::SparseMatrix& tglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleNT: Node ownership inconsistency!");

    if (Dim()==2)
    {
      // prepare assembly
      int colsize = cnode->NumDof();
      vector<int> lmrowN(1);
      vector<int> lmrowT(1);
      vector<int> lmrowownerN(1);
      vector<int> lmrowownerT(1);
      vector<int> lmcol(colsize);

      lmrowN[0] = activen_->GID(i);
      lmrowownerN[0] = cnode->Owner();
      lmrowT[0] = activet_->GID(i);
      lmrowownerT[0] = cnode->Owner();

      /**************************************************** N-matrix ******/
      Epetra_SerialDenseMatrix Nnode(1,colsize);

      // we need D diagonal entry of this node
      double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(1,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->txi()[j];
      }

      // assemble into matrix of normal vectors T
      tglobal.Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
    }

    else if (Dim()==3)
    {
      // prepare assembly
      int colsize = cnode->NumDof();
      vector<int> lmrowN(1);
      vector<int> lmrowT(2);
      vector<int> lmrowownerN(1);
      vector<int> lmrowownerT(2);
      vector<int> lmcol(colsize);

      lmrowN[0] = activen_->GID(i);
      lmrowownerN[0] = cnode->Owner();
      lmrowT[0] = activet_->GID(2*i);
      lmrowT[1] = activet_->GID(2*i+1);
      lmrowownerT[0] = cnode->Owner();
      lmrowownerT[1] = cnode->Owner();

      /**************************************************** N-matrix ******/
      Epetra_SerialDenseMatrix Nnode(1,colsize);

      // we need D diagonal entry of this node
      double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Nnode(0,j) = wii * cnode->n()[j];
      }

      // assemble into matrix of normal vectors N
      nglobal.Assemble(-1,Nnode,lmrowN,lmrowownerN,lmcol);

      /**************************************************** T-matrix ******/
      Epetra_SerialDenseMatrix Tnode(2,colsize);

      for (int j=0;j<colsize;++j)
      {
        lmcol[j] = cnode->Dofs()[j];
        Tnode(0,j) = cnode->txi()[j];
        Tnode(1,j) = cnode->teta()[j];
      }

      // assemble into matrix of normal vectors T
      tglobal.Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
    }
    else
      dserror("ERROR: Dim() must be either 2D or 3D");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives           popp 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    map<int,double>& dgmap = cnode->GetDerivG();
    map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;
      //cout << "Assemble S: " << row << " " << col << " " << val << endl;
      // do not assemble zeros into s matrix
      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix P containing tangent derivatives          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleP(LINALG::SparseMatrix& pglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleP: Node ownership inconsistency!");

    if (Dim()==2)
    {
      // prepare assembly
      vector<map<int,double> >& dtmap = cnode->GetDerivTxi();
      map<int,double>::iterator colcurr;
      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();
      int row = activet_->GID(i);

      if (mapsize==3) mapsize=2;

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleP: Column dim. of nodal DerivT-map is inconsistent!");

      // begin assembly of P-matrix
      //cout << endl << "->Assemble P for Node ID: " << cnode->Id() << endl;

      // loop over all derivative maps (=dimensions)
      for (int j=0;j<mapsize;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << row << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,row,col);
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsize);
      }
    }
    else if (Dim()==3)
    {
      // prepare assembly
      vector<map<int,double> >& dtximap = cnode->GetDerivTxi();
      vector<map<int,double> >& dtetamap = cnode->GetDerivTeta();
      map<int,double>::iterator colcurr;
      int colsizexi = (int)dtximap[0].size();
      int colsizeeta = (int)dtetamap[0].size();
      int mapsizexi = (int)dtximap.size();
      int mapsizeeta = (int)dtetamap.size();
      int rowxi = activet_->GID(2*i);
      int roweta = activet_->GID(2*i+1);

      for (int j=0;j<mapsizexi-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleS: Column dim. of nodal DerivTXi-map is inconsistent!");

      for (int j=0;j<mapsizeeta-1;++j)
        if ((int)dtetamap[j].size() != (int)dtetamap[j+1].size())
          dserror("ERROR: AssembleS: Column dim. of nodal DerivTEta-map is inconsistent!");

      // begin assembly of P-matrix
      //cout << endl << "->Assemble P for Node ID: " << cnode->Id() << endl;

      // loop over all derivative maps (=dimensions) for TXi
      for (int j=0;j<mapsizexi;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << rowxi << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,rowxi,col);
          ++k;
        }

        if (k!=colsizexi)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsizexi);
      }

      // loop over all derivative maps (=dimensions) for TEta
      for (int j=0;j<mapsizeeta;++j)
      {
        int k=0;

        // loop over all entries of the current derivative map
        for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = cnode->lm()[j]*(colcurr->second);
          //cout << "lm[" << j << "]=" << cnode->lm()[j] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble P: " << roweta << " " << col << " " << val << endl;
          // do not assemble zeros into P matrix
          if (abs(val)>1.0e-12) pglobal.Assemble(val,roweta,col);
          ++k;
        }

        if (k!=colsizeeta)
          dserror("ERROR: AssembleP: k = %i but colsize = %i",k,colsizeeta);
      }
    }
    else
      dserror("ERROR: Dim() must be either 2 or 3!");

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrices LinDM containing fc derivatives         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                                       LINALG::SparseMatrix& linmglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  /********************************************** LinDMatrix **********/
  // This is easy and can be done without communication, as the global
  // matrix lind has the same row map as the storage of the derivatives
  // of the D matrix.
  /**********************************************************************/

  // loop over all slave nodes (row map)
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i) // jtilde
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    int dim = cnode->NumDof();

    // Mortar matrix D derivatives
    map<int,map<int,double> >& dderiv = cnode->GetDerivD();

    map<int,double>::iterator colcurr;
    map<int,map<int,double> >::iterator dcurr;

    for( int j=0; j<dim; ++j ) // j
    {
      int row = cnode->Dofs()[j];

      // loop over all slave nodes dof again using the map iterator
      for ( dcurr = dderiv.begin(); dcurr != dderiv.end(); ++dcurr  ) // ktilde
      {
        // get the corresponding node
        DRT::Node* snode = idiscret_->gNode(dcurr->first);
        if (!snode)
          dserror("ERROR: Cannot find node with gid %",dcurr->first);
        CNode* csnode = static_cast<CNode*>(snode);

        double* lm = csnode->lm();

        // get dderiv row
        map<int,double>& dderivrow = dderiv[dcurr->first];

        // loop over all directional derivative entries using the map iterator
        for( colcurr = dderivrow.begin(); colcurr != dderivrow.end(); ++colcurr ) // l
        {
          int col = colcurr->first;
          double val = lm[j] * (colcurr->second);
          lindglobal.Assemble(val, row, col);
        }
      }
    }

    /** old version, using dual shape functions **
    // Mortar matrix D derivatives
    map<int,double>& dderiv = cnode->GetDerivD();
    map<int,double>::iterator colcurr;

    // current Lagrange multipliers
    double* lm = cnode->lm();

    // loop over all slave node dofs for assembly
    for (int j=0; j<numdof; ++j)
    {
      int row = cnode->Dofs()[j]; // row index equals global dof index of this #i slave node's dof j

      cout << "  j=" << j  << " row=" << row ;

      // loop over all directional derivative entries
      for (colcurr=dderiv.begin(); colcurr!=dderiv.end(); ++colcurr)
      {
        int col = colcurr->first; // col index equals global id of directional derivative component ,l

        double val = lm[j] * (colcurr->second);

        cout << " | col=" << col << ", val=" << val << " | ";

        if (abs(val)>1.0e-12)
          lindglobal.Assemble(val, row, col);
      }

      cout << " " << endl;
    }
    */
  }

  /********************************************** LinMMatrix ************/
  // This is a quite complex task and we have to do a lot of parallel
  // communication here!!! The reason for this is that the directional
  // derivative of M is stored slave-node-wise meaning that we can only
  // address the derivative of M via the slave node row map! But now
  // we have to assemble into a matrix linm which is based on the
  // master node row map!!! This is not straight forward as we CANNOT
  // do this without communication: when we have collected one slave
  // node's contribution to linm we have to assemble this portion to
  // the respective master nodes in parallel! Only the processor owning
  // the master node can do this! Therefore, we communicate the derivatives
  // of M to all procs and then explicitly call the current master proc
  // to do the assembly into the sparse matrix linm.
  /**********************************************************************/

  // loop over all slave nodes (full map)
  for (int i=0;i<snodefullmap_->NumMyElements();++i)
  {
    int gid = snodefullmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    // Mortar matrix M derivatives
    map<int,map<int,double> >& mderiv = cnode->GetDerivM();
    map<int,double>::iterator colcurr;
    map<int,map<int,double> >::iterator mcurr;

    // current Lagrange multipliers
    double* lm = cnode->lm();

    // inter-proc. communication
    // (we want to keep all procs around, although at the moment only
    // the cnode owning proc does relevant work!)
    int mastersize = 0;
    if (Comm().MyPID()==cnode->Owner())
    {
      mastersize = (int)mderiv.size();
      mcurr = mderiv.begin();
    }
    //********************************************************************
    // important notes:
    // 1) we use lComm here, as we only communicate among participating
    // procs of this interface (all others are already out of here)
    // 2) the broadcasting proc has to be given in lComm numbering, not
    // Comm numbering, which is achieved by calling the map procmap_
    //********************************************************************
    lComm()->Broadcast(&mastersize,1,procmap_[cnode->Owner()]);

    // loop over all master nodes in the DerivM-map of the current slave node
    for (int j=0;j<mastersize;++j)
    {
      int mgid = 0;
      if (Comm().MyPID()==cnode->Owner())
      {
        mgid = mcurr->first;
        ++mcurr;
      }
      lComm()->Broadcast(&mgid,1,procmap_[cnode->Owner()]);

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
      CNode* cmnode = static_cast<CNode*>(mnode);
      int mnumdof = cmnode->NumDof();

      // Mortar matrix M derivatives
      map<int,double>& thismderiv = cnode->GetDerivM()[mgid];
      int mapsize = 0;
      if (Comm().MyPID()==cnode->Owner())
        mapsize = (int)(thismderiv.size());
      lComm()->Broadcast(&mapsize,1,procmap_[cnode->Owner()]);

      // loop over all master node dofs
      for (int k=0;k<mnumdof;++k)
      {
        int row = cmnode->Dofs()[k];

        if (Comm().MyPID()==cnode->Owner())
          colcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = 0;
          double val = 0.0;
          if (Comm().MyPID()==cnode->Owner())
          {
            col = colcurr->first;
            val = lm[k] * (colcurr->second);
            ++colcurr;
          }
          lComm()->Broadcast(&col,1,procmap_[cnode->Owner()]);
          lComm()->Broadcast(&val,1,procmap_[cnode->Owner()]);

          // owner of master node has to do the assembly!!!
          if (Comm().MyPID()==cmnode->Owner())
          {
            //cout << "Assemble LinM: " << row << " " << col << " " << val << endl;
            if (abs(val)>1.0e-12) linmglobal.Assemble(-val,row,col);
          }
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (Comm().MyPID()==cnode->Owner() && colcurr!=thismderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
      }
    }

    // check for completeness of DerivM-Master-iteration
    if (Comm().MyPID()==cnode->Owner() && mcurr!=mderiv.end())
      dserror("ERROR: AssembleLinDM: Not all master entries of DerivM considered!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinStick with tangential+D+M derivatives  mgit 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
                                          LINALG::SparseMatrix& linstickDISglobal,
                                          Epetra_Vector& linstickRHSglobal)
{
	// FIXGIT: Assemble LinStick is containing a matrix for the de-
	// rivatives of the Lagrange multipliers. This is according to Hüeber.
	// Because of worse convergence, this is not implemented, but the
	// code is commented after the algorithm.

	// get out of here if not participating in interface
  if (!lComm())
    return;

  // create map of stick nodes
  RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
  RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements()==0)
    return;

  // loop over all stick nodes of the interface
  for (int i=0;i<sticknodes->NumMyElements();++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

    // prepare assembly, get information from node
    vector<map<int,double> > dtximap = cnode->GetDerivTxi();
    vector<map<int,double> > dtetamap = cnode->GetDerivTeta();

   for (int j=0;j<Dim()-1;++j)
 	  	if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

	  if (Dim()==3)
	  {
	    for (int j=0;j<Dim()-1;++j)
 	  	  if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
	  }

	  // more information from node
    double* xi = cnode->xspatial();
    double* txi = cnode->txi();
    double* teta = cnode->teta();
    double* jump = cnode->jump();

    // iterator for maps
    map<int,double>::iterator colcurr;

    // row number of entries
 	  vector<int> row (Dim()-1);
 	  if (Dim()==2)
 	  {
 	    row[0] = stickt->GID(i);
 	  }
 	  else if (Dim()==3)
 	  {
 	  	row[0] = stickt->GID(2*i);
 	  	row[1] = stickt->GID(2*i)+1;
 	  }
 	  else
 	  	dserror("ERROR: AssemblelinStick: Dimension not correct");

    // evaluation of specific components of entries to assemble
 	  double jumptxi=0;
    double jumpteta=0;
    for (int dim = 0;dim < Dim();dim++)
    {
      jumptxi += txi[dim]*jump[dim];
      jumpteta += teta[dim]*jump[dim];
    }

     // check for dimensions
	   if(Dim()==2 and (jumpteta != 0.0))
	   	dserror ("ERROR: AssembleLinSlip: jumpteta must be zero in 2D");

    // Entries on right hand side
    /************************************************ (-utxi, -uteta) ***/
    Epetra_SerialDenseVector rhsnode(Dim()-1);
    vector<int> lm(Dim()-1);
    vector<int> lmowner(Dim()-1);

    rhsnode(0) = -jumptxi;
    lm[0] = cnode->Dofs()[1];
    lmowner[0] = cnode->Owner();

    if (Dim()==3)
    {
    	rhsnode(1) = -jumpteta;
      lm[1] = cnode->Dofs()[2];
      lmowner[1] = cnode->Owner();
    }

    LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);

    // Entries from differentiation with respect to displacements
    /*** 01 **************************************** tangent.(D-Dn-1) ***/

    // we need the nodal entries of the D-matrix and the old ones
    // they must have an entry
    vector<map<int,double> > dmapold=cnode->GetDOld();
    if(dmapold.size()<1)
        dserror("ERROR: AssembleLinStick: No dmap-entries form previous time step");

    double D= (cnode->GetD()[0])[cnode->Dofs()[0]];
    double Dold= dmapold[0][cnode->Dofs()[0]];

    // loop over dimensions
    for (int dim=0;dim<cnode->NumDof();++dim)
    {
      int col = cnode->Dofs()[dim];
      double valtxi = (-1)*txi[dim]*(D-Dold);

      // do not assemble zeros into matrix
      if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

      if(Dim()==3)
      {
        double valteta = (-1)*teta[dim]*(D-Dold);
        if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
      }
    }

    /*** 02 **************************************** tangent.(M-Mn-1) ***/
    // we need the nodal entries of the M-matrix and the old ones
    vector<map<int,double> > mmap = cnode->GetM();
    vector<map<int,double> > mmapold = cnode->GetMOld();

    // map from previous time step must have an entry
    if(mmapold.size()<1)
       dserror("ERROR: AssembleLinStick: No mmap-entries form previous time step");

    // create a set of nodes including nodes according to M entries
    // from current and previous time step
    set <int> mnodes;

    // iterator
    set<int>::iterator mcurr;

    set <int> mnodescurrent = cnode->GetMNodes();
  	set <int> mnodesold = cnode->GetMNodesOld();

  	for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
  	  mnodes.insert(*mcurr);

  	for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
  	  mnodes.insert(*mcurr);

    // loop over all master nodes (find adjacent ones to this stick node)
    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
    {
      int gid = *mcurr;
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cmnode = static_cast<CNode*>(mnode);
      const int* mdofs = cmnode->Dofs();

      double mik = (mmap[0])[mdofs[0]];
      double mikold = (mmapold[0])[mdofs[0]];

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cmnode->Dofs()[dim];
        double valtxi = txi[dim]*(mik-mikold);

        // do not assemble zeros into matrix
        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

        if (Dim()==3)
        {
          double valteta = teta[dim]*(mik-mikold);
          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }
    }

    /*** 1 ******************************** deriv(tangent).(D-Dn-1).x ***/
    // loop over dimensions
    for (int j=0;j<Dim();++j)
    {
      // loop over all entries of the current derivative map (txi)
      for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
      {
        int col = colcurr->first;
        double val = (-1)*(D-Dold)*xi[j]*colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
      }

      if(Dim()==3)
      {
      	// loop over all entries of the current derivative map (teta)
        for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*(D-Dold)*xi[j]*colcurr->second;

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
        }
      }
    }

    /*** 2 ******************************** deriv(tangent).(M-Mn-1).x ***/
    // loop over all master nodes (find adjacent ones to this stick node)
    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
    {
      int gid = *mcurr;
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cmnode = static_cast<CNode*>(mnode);
      const int* mdofs = cmnode->Dofs();

      double mik = (mmap[0])[mdofs[0]];
      double mikold = (mmapold[0])[mdofs[0]];
      double* mxi = cmnode->xspatial();

      // loop over dimensions
      for (int j=0;j<Dim();++j)
      {
        // loop over all entries of the current derivative map (dtxi)
        for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = (mik-mikold)*mxi[j]*colcurr->second;
          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
        }

        if(Dim()==3)
        {
          // loop over all entries of the current derivative map (dteta)
          for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (mik-mikold)*mxi[j]*colcurr->second;
            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
          }
        }
      }
    }

    /*** 3 ************************************** tangent.Deriv(D).x  ***/
    // we need the dot product txi.x and teta.x of this node
    double txidotx = 0.0;
    double tetadotx = 0.0;
    for (int dim=0;dim<cnode->NumDof();++dim)
    {
    	txidotx += txi[dim]*xi[dim];
    	tetadotx += teta[dim]*xi[dim];
    }

    // prepare assembly
    map<int,double>& ddmap = cnode->GetDerivD()[gid];

    // loop over all entries of the current derivative map
    for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
    {
      int col = colcurr->first;
      double valtxi = (-1)*txidotx*colcurr->second;

      if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

      if(Dim()==3)
      {
      	double valteta = (-1)*tetadotx*colcurr->second;
      	if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
      }
    }

    /*** 4 *************************************** tangent.Deriv(M).x ***/
    // we need the Lin(M-matrix) entries of this node
    map<int,map<int,double> >& dmmap = cnode->GetDerivM();
    map<int,map<int,double> >::iterator dmcurr;

    // loop over all master nodes in the DerivM-map of the active slave node
    for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
    {
      int gid = dmcurr->first;
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cmnode = static_cast<CNode*>(mnode);
      double* mxi = cmnode->xspatial();

      // we need the dot product ns*xm of this node pair
      double txidotx = 0.0;
      double tetadotx = 0.0;
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
      	txidotx += txi[dim]*mxi[dim];
      	tetadotx += teta[dim]*mxi[dim];

      }
      // compute matrix entry of the current active node / master node pair
      map<int,double>& thisdmmap = cnode->GetDerivM(gid);

      // loop over all entries of the current derivative map
      for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
      {
        int col = colcurr->first;
        double val = txidotx*colcurr->second;

        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);

        if (Dim()==3)
        {
        	double val = tetadotx*colcurr->second;
        	if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
        }
      }
    }
  }

//  have a look to the beginning of the function
//  // only for coulomb friction
//  string ftype   = IParams().get<string>("friction type","none");
//  double frcoeff = IParams().get<double>("friction coefficient",0.0);
//  double cn = IParams().get<double>("semismooth cn",0.0);
//  double ct = IParams().get<double>("semismooth ct",0.0);
//
//  if (ftype == "tresca")
//    dserror ("Error: AssemblelinStick: complementary function according"
//             " to Hueber only available for Coulomb friction");
//
//  // get out of here if not participating in interface
//  if (!lComm())
//    return;
//
//  // create map of stick nodes
//  RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
//  RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);
//
//  // nothing to do if no stick nodes
//  if (sticknodes->NumMyElements()==0)
//    return;
//
//  // not yet implemented for 3D
//    if (Dim()==3)
//      dserror("ERROR: AssembleLinStick: 3D not yet implemented");
//
//  // loop over all stick nodes of the interface
//  for (int i=0;i<sticknodes->NumMyElements();++i)
//  {
//    int gid = sticknodes->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
//    CNode* cnode = static_cast<CNode*>(node);
//
//    if (cnode->Owner() != Comm().MyPID())
//      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");
//
//    // prepare assembly
//    vector<map<int,double> > dnmap = cnode->GetDerivN();
//    map<int,double>::iterator colcurr;
//
//    // calculate DerivT from DerivN
//    // only for 2D so far, in this case calculation is very easy
//    // dty =  dnx
//    // dtx = -dny
//
//    vector <map<int,double> > dtmap(Dim());
//
//    for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
//      dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));
//
//    for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
//      dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));
//
//    int colsize = (int)dtmap[0].size();
//    int mapsize = (int)dtmap.size();
//    int row = stickt->GID(i);
//    double* xi = cnode->xspatial();
//    double* txi = cnode->txi();
//    double* jump = cnode->jump();
//    double& wgap = cnode->Getg();
//
//    double utan = 0;
//    double nz = 0;
//
//    for (int dim = 0;dim < Dim();dim++)
//    {
//      utan += txi[dim]*jump[dim];
//      nz += cnode->n()[dim] * cnode->lm()[dim];
//    }
//
//    // initialization of nz if nz = wgap = 0
//    if(nz==wgap and wgap==0)
//    {
//      nz=1;
//      cout << "Warning: Initialization of nz to 1" << endl;
//    }
//
//    for (int j=0;j<mapsize-1;++j)
//      if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
//        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivT-map is inconsistent!");
//
//    // Entries on right hand side
//    /**************************************** frcoeff*cn*wgap*ct*utan ***/
//
//    Epetra_SerialDenseVector rhsnode(1);
//    vector<int> lm(1);
//    vector<int> lmowner(1);
//
//    rhsnode(0) = -frcoeff*cn*wgap*ct*utan;
//    lm[0] = cnode->Dofs()[1];
//    lmowner[0] = cnode->Owner();
//
//    LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
//
//    // Entries from differentiation with respect to lagrange multipliers
//    /*******************/
//
//    // loop over the dimension
//    for (int dim=0;dim<cnode->NumDof();++dim)
//    {
//      int col = cnode->Dofs()[dim];
//      double val = -frcoeff*cnode->n()[dim]*ct*utan;
//      // do not assemble zeros into matrix
//      if (abs(val)>1.0e-12) linstickLMglobal.Assemble(val,row,col);
//    }
//
//    // Entries from differentiation with respect to displacements
//    /************************************************** -tan.(D-Dn-1) ***/
//
//    // we need the nodal entries of the D-matrix and the old one
//    double D= (cnode->GetD()[0])[cnode->Dofs()[0]];
//    double Dold= (cnode->GetDOld()[0])[cnode->Dofs()[0]];
//
//    // loop over all derivative maps (=dimensions)
//    for (int dim=0;dim<cnode->NumDof();++dim)
//    {
//      int col = cnode->Dofs()[dim];
//      double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*txi[dim]*(D-Dold);
//     //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
//
//     // do not assemble zeros into matrix
//     if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /*************************************************** tan.(M-Mn-1) ***/
//
//    // we need the nodal entries of the M-matrix and the old one
//    vector<map<int,double> > mmap = cnode->GetM();
//    vector<map<int,double> > mmapold = cnode->GetMOld();
//
//    // create a set of nodes including nodes according to M entries
//    // from current and previous time step
//    set <int> mnodes;
//
//    for (colcurr=mmap[0].begin(); colcurr!=mmap[0].end(); colcurr++)
//      mnodes.insert((colcurr->first)/Dim());
//
//    if(mmapold.size()<1)
//    {
//      cout << "GID " << gid << endl;
//      dserror("vector too small");
//    }
//
//    for (colcurr=mmapold[0].begin(); colcurr!=mmapold[0].end(); colcurr++)
//      mnodes.insert((colcurr->first)/Dim());
//
//    set<int>::iterator mcurr;
//
//    // loop over all master nodes (find adjacent ones to this stick node)
//    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
//    {
//      int gid = *mcurr;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      CNode* cmnode = static_cast<CNode*>(mnode);
//      const int* mdofs = cmnode->Dofs();
//
//      double mik = (mmap[0])[mdofs[0]];
//      double mikold = (mmapold[0])[mdofs[0]];
//
//      // compute linstick-matrix entry of the current active node / master node pair
//      // loop over all derivative maps (=dimensions)
//      for (int dim=0;dim<cnode->NumDof();++dim)
//      {
//        int col = cmnode->Dofs()[dim];
//        double val = -frcoeff*(nz-cn*wgap)*ct*txi[dim]*(mik-mikold);
//        //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
//
//       // do not assemble zeros into matrix
//       if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//      }
//    }
//
//    /***************************************  -DerivT.(D-Dn-1).xs  ******/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all derivative maps (=dimensions)
//    for (int j=0;j<mapsize;++j)
//    {
//      int k=0;
//
//      // loop over all entries of the current derivative map
//      for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*(D-Dold)*xi[j]*colcurr->second;
//
//        // do not assemble zeros into s matrix
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//        ++k;
//      }
//
//      if (k!=colsize)
//        dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//    }
//
//    /***************************************  -DerivT.(M-Mn-1).xm  ******/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all master nodes (find adjacent ones to this stick node)
//    for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
//    {
//      int gid = *mcurr;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      CNode* cmnode = static_cast<CNode*>(mnode);
//      const int* mdofs = cmnode->Dofs();
//
//      double mik = (mmap[0])[mdofs[0]];
//      double mikold = (mmapold[0])[mdofs[0]];
//
//      double* mxi = cmnode->xspatial();
//
//      // compute linstick-matrix entry of the current active node / master node pair
//      // loop over all derivative maps (=dimensions)
//      for (int j=0;j<mapsize;++j)
//      {
//        int k=0;
//
//        // loop over all entries of the current derivative map
//        for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
//        {
//          int col = colcurr->first;
//          double val = -frcoeff*(nz-cn*wgap)*ct*(mik-mikold)*mxi[j]*colcurr->second;
//          // do not assemble zeros into matrix
//          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//          ++k;
//        }
//
//        if (k!=colsize)
//          dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//      }
//    }
//
//    /**********************************************  -T.DerivD.x  *******/
//
//    // we need the dot product n*x of this node
//    double tdotx = 0.0;
//    for (int dim=0;dim<cnode->NumDof();++dim)
//      tdotx += txi[dim]*xi[dim];
//
//    // prepare assembly
//    map<int,double>& ddmap = cnode->GetDerivD();
//
//    // loop over all entries of the current derivative map
//    for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
//    {
//      int col = colcurr->first;
//      double val = -frcoeff*(nz-cn*wgap)*ct*(-1)*tdotx*colcurr->second;
//
//      if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /***********************************************   -T.DerivM.x ******/
//
//    // we need the Lin(M-matrix) entries of this node
//    map<int,map<int,double> >& dmmap = cnode->GetDerivM();
//    map<int,map<int,double> >::iterator dmcurr;
//
//    // loop over all master nodes in the DerivM-map of the active slave node
//    for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
//    {
//      int gid = dmcurr->first;
//      DRT::Node* mnode = idiscret_->gNode(gid);
//      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
//      CNode* cmnode = static_cast<CNode*>(mnode);
//      double* mxi = cmnode->xspatial();
//
//      // we need the dot product ns*xm of this node pair
//      double tdotx = 0.0;
//      for (int dim=0;dim<cnode->NumDof();++dim)
//        tdotx += txi[dim]*mxi[dim];
//
//      // compute matrix entry of the current active node / master node pair
//      map<int,double>& thisdmmap = cnode->GetDerivM(gid);
//
//      // loop over all entries of the current derivative map
//      for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*(nz-cn*wgap)*ct*tdotx*colcurr->second;
//
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//      }
//    }
//
//    /**************************************-frcoeff*cn*ct*utan*DerivG ***/
//
//    // prepare assembly
//    map<int,double>& dgmap = cnode->GetDerivG();
//
//    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
//    {
//      int col = colcurr->first;
//      double val = +frcoeff*cn*ct*utan*colcurr->second;
//      //cout << "Assemble LinStick: " << row << " " << col << " " << val << endl;
//      // do not assemble zeros into matrix
//      if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//    }
//
//    /************************************** -frcoeff*DerivN*z+ct*utan ***/
//    // we need the nodal entries of the D-matrix and the old one
//
//    // loop over all derivative maps (=dimensions)
//    for (int j=0;j<mapsize;++j)
//    {
//      int k=0;
//
//      // loop over all entries of the current derivative map
//      for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
//      {
//        int col = colcurr->first;
//        double val = -frcoeff*cnode->lm()[j]*ct*utan*colcurr->second;
//
//        // do not assemble zeros into s matrix
//        if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row,col);
//        ++k;
//      }
//
//      if (k!=colsize)
//        dserror("ERROR: AssembleLinStick: k = %i but colsize = %i",k,colsize);
//    }
//  }
  return;
}

/*----------------------------------------------------------------------*
|  Assemble matrix LinSlip with tangential+D+M derivatives    mgit 02/09|
*----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
                                         LINALG::SparseMatrix& linslipDISglobal,
                                         Epetra_Vector& linslipRHSglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  // information from interface contact parameter list
  INPAR::CONTACT::ContactFrictionType ftype =
    Teuchos::getIntegralValue<INPAR::CONTACT::ContactFrictionType>(IParams(),"FRICTION");
  double frbound = IParams().get<double>("FRBOUND");
  double frcoeff = IParams().get<double>("FRCOEFF");
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  bool fulllin = Teuchos::getIntegralValue<int>(IParams(),"FULL_LINEARIZATION");

  // Coulomb Friction
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
#ifdef CONTACTCOMPHUEBER

  	// loop over all slip nodes of the interface
		for (int i=0;i<slipnodes_->NumMyElements();++i)
		{
		  int gid = slipnodes_->GID(i);
		  DRT::Node* node = idiscret_->gNode(gid);
		  if (!node) dserror("ERROR: Cannot find node with gid %",gid);
		  CNode* cnode = static_cast<CNode*>(node);

		  if (cnode->Owner() != Comm().MyPID())
		    dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // prepare assembly, get information from node
		  vector<map<int,double> > dnmap = cnode->GetDerivN();
	    vector<map<int,double> > dtximap = cnode->GetDerivTxi();
	    vector<map<int,double> > dtetamap = cnode->GetDerivTeta();

      // check for Dimension of derivative maps
		  for (int j=0;j<Dim()-1;++j)
	 	  	if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
	        dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

		   for (int j=0;j<Dim()-1;++j)
		 	  	if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
		        dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

			 if (Dim()==3)
			 {
			   for (int j=0;j<Dim()-1;++j)
		 	 	  if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
		        dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
			 }

      // more information from node
	    double* jump = cnode->jump();
	    double* n = cnode->n();
	    double* txi = cnode->txi();
	    double* teta = cnode->teta();
	    double* xi = cnode->xspatial();
		  double* z = cnode->lm();
		  double& wgap = cnode->Getg();

		  // iterator for maps
	    map<int,double>::iterator colcurr;

	 	  // row number of entries
	 	  vector<int> row (Dim()-1);
	 	  if (Dim()==2)
	 	  {
	 	    row[0] = slipt_->GID(i);
	 	  }
	 	  else if (Dim()==3)
	 	  {
	 	  	row[0] = slipt_->GID(2*i);
	 	  	row[1] = slipt_->GID(2*i)+1;
	 	  }
	 	  else
	 	  	dserror("ERROR: AssemblelinSlip: Dimension not correct");

	 	  // boolean variable if flag "CONTACTFRICTIONLESSFIRST" AND
	 	  // ActiveOld = true
	 	  bool friclessandfirst = false;

	 	  // evaluation of specific components of entries to assemble
	    double znor = 0;
	    double ztxi = 0;
	    double zteta = 0;
	    double jumptxi = 0;
	    double jumpteta = 0;
	    double euclidean = 0;
	    for (int i=0;i<Dim();i++)
	    {
	    	znor += n[i]*z[i];
	    	ztxi += txi[i]*z[i];
	    	zteta += teta[i]*z[i];
	    	jumptxi += txi[i]*jump[i];
	    	jumpteta += teta[i]*jump[i];
	    }

    	// evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
    	vector<double> sum1 (Dim()-1,0);
    	sum1[0] = ztxi+ct*jumptxi;
    	if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
    	if (Dim()==2) euclidean = abs(sum1[0]);
    	if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

	    // check of dimensions
	    if(Dim()==2 and (zteta != 0.0 or jumpteta != 0.0))
	    	dserror ("ERROR: AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      // check of euclidean norm
    	if (euclidean==0.0)
    		dserror ("ERROR: AssemblelinSlip: Euclidean norm is zero");

#ifdef CONTACTFRICTIONLESSFIRST

    	// in the case of frictionless contact for nodes just coming into
    	// contact, the frictionless contact condition is applied.
    	if (cnode->ActiveOld()==false)
  	  {
    		friclessandfirst=true;
    		for (int dim=0;dim<cnode->NumDof();++dim)
  	    {
    			int col = cnode->Dofs()[dim];
		    	double valtxi = txi[dim];
		    	double valteta = 0;
		    	if (Dim()==3) valteta = teta[dim];

		      if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
		      if (Dim()==3)
		        if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);

  	    }
    		if(fulllin)
  	    {
    			for (int dim=0;dim<cnode->NumDof();++dim)
          {
		        for (colcurr=dtximap[dim].begin();colcurr!=dtximap[dim].end();++colcurr)
		        {
		        	int col = colcurr->first;
		        	double valtxi = (colcurr->second)*z[dim];
		        	if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
		        }

		        if(Dim()==3)
		        {
		        	for (colcurr=dtetamap[dim].begin();colcurr!=dtetamap[dim].end();++colcurr)
	            {
		        		int col = colcurr->first;
		      	    double valteta = (colcurr->second)*z[dim];
		      	    if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
	            }
		        }
		      }
  	    }
  	  }
#endif

    	// this is not evaluated if "FRICTIONLESSFIRST" is flaged on AND the node
    	// is just coming into contact
    	if(friclessandfirst==false)
    	{
				/******************************************************************/
				// calculation of matrix entries of the linearized slip condition
				/******************************************************************/
				// 1) Entries from differentiation with respect to lagrange multipliers
				// 2) Entries on right hand side
				// 3) Entries from differentiation with respect to displacements

				// 1) Entries from differentiation with respect to lagrange multipliers
				/******************************************************************/

				// loop over the dimension
				for (int dim=0;dim<cnode->NumDof();++dim)
				{
					double valtxi = 0;
					int col = cnode->Dofs()[dim];
					double valtxi0 = euclidean*txi[dim];
					double valtxi1 = ((ztxi+ct*jumptxi)/euclidean*ztxi)*txi[dim];
					double valtxi3 = (zteta+ct*jumpteta)/euclidean*ztxi*teta[dim];
					double valtxi2 = -frcoeff*(znor-cn*wgap)*txi[dim]-frcoeff*(ztxi+ct*jumptxi)*n[dim];
					valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

					double valteta = 0;
					if (Dim()==3)
					{
						double valteta0 = euclidean*teta[dim];
						double valteta1 = ((ztxi+ct*jumptxi)/euclidean*zteta)*txi[dim];
						double valteta3 = (zteta+ct*jumpteta)/euclidean*zteta*teta[dim];
						double valteta2 = -frcoeff*(znor-cn*wgap)*teta[dim]-frcoeff*(zteta+ct*jumpteta)*n[dim];
						valteta = valteta0 + valteta1 + valteta2 + valteta3;
					}

					// do not assemble zeros into matrix
					if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
					if (Dim()==3)
						if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
				}

				// 2) Entries on right hand side
				/************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

				double valuetxi1 = -(euclidean)*ztxi+(frcoeff*(znor-cn*wgap))*(ztxi+ct*jumptxi);
				double valuetxi2 = +euclidean*ztxi;
				double valuetxi3 = (ztxi+ct*jumptxi)/euclidean*ztxi*ztxi;
				double valuetxi4 = (zteta+ct*jumpteta)/euclidean*zteta*ztxi;
				double valuetxi5 = -(frcoeff*(znor-cn*wgap))*ztxi-(frcoeff*znor)*(ztxi+ct*jumptxi);

				Epetra_SerialDenseVector rhsnode(Dim()-1);
				vector<int> lm(Dim()-1);
				vector<int> lmowner(Dim()-1);

				rhsnode(0) = (valuetxi1+valuetxi2+valuetxi3+valuetxi4+valuetxi5);
				lm[0] = cnode->Dofs()[1];
				lmowner[0] = cnode->Owner();

				if(Dim()==3)
				{
					double valueteta1 = -(euclidean)*zteta+(frcoeff*(znor-cn*wgap))*(zteta+ct*jumpteta);
					double valueteta2 = +euclidean*zteta;
					double valueteta3 = (ztxi+ct*jumptxi)/euclidean*ztxi*zteta;
					double valueteta4 = (zteta+ct*jumpteta)/euclidean*zteta*zteta;
					double valueteta5 = -(frcoeff*(znor-cn*wgap))*zteta-(frcoeff*znor)*(zteta+ct*jumpteta);

					rhsnode(1) = (valueteta1+valueteta2+valueteta3+valueteta4+valueteta5);
					lm[1] = cnode->Dofs()[2];
					lmowner[1] = cnode->Owner();
				}

				LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

				// 3) Entries from differentiation with respect to displacements
				/*** 01  *********** -Deriv(euclidean).ct.tangent.(D-Dn-1)*ztan ***/

				// we need the nodal entries of the D-matrix and the old one
				vector<map<int,double> > dmapold=cnode->GetDOld();
				if(dmapold.size()<1)
						dserror("ERROR: AssembleLinSlip: No dmap-entries form previous time step");

				double D= (cnode->GetD()[0])[cnode->Dofs()[0]];
				double Dold= dmapold[0][cnode->Dofs()[0]];

				// loop over dimensions
				for (int dim=0;dim<cnode->NumDof();++dim)
				{
					double valtxi1=0;
					double valtxi2=0;
					double valteta1=0;
					double valteta2=0;

					int col = cnode->Dofs()[dim];
					valtxi1 = (ztxi+ct*jumptxi)/euclidean*(-1)*ct*txi[dim]*(D-Dold)*ztxi;

					if(Dim()==3)
					{
						valteta1 = (ztxi+ct*jumptxi)/euclidean*(-1)*ct*txi[dim]*(D-Dold)*zteta;
						valtxi2 = (zteta+ct*jumpteta)/euclidean*(-1)*ct*teta[dim]*(D-Dold)*ztxi;
						valteta2 = (zteta+ct*jumpteta)/euclidean*(-1)*ct*teta[dim]*(D-Dold)*zteta;
					}

					// do not assemble zeros into matrix
					if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
					if (Dim()==3)
					{
						if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
						if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
						if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
					}
				 }

				/*** 02  *********** -Deriv(euclidean).ct.tangent.(M-Mn-1)*ztan ***/

				// we need the nodal entries of the M-matrix and the old one
				vector<map<int,double> > mmap = cnode->GetM();
				vector<map<int,double> > mmapold = cnode->GetMOld();

				// map from previous time step must have an entry
				if(mmapold.size()<1)
					 dserror("ERROR: AssembleLinStick: No mmap-entries form previous time step");

				// create a set of nodes including nodes according to M entries
				 // from current and previous time step
				 set <int> mnodes;

				 // iterator
				 set<int>::iterator mcurr;

				 set <int> mnodescurrent = cnode->GetMNodes();
				set <int> mnodesold = cnode->GetMNodesOld();

				for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
					mnodes.insert(*mcurr);

				for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
					mnodes.insert(*mcurr);

				// loop over all master nodes (find adjacent ones to this slip node)
				for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
				{
					int gid = *mcurr;
					DRT::Node* mnode = idiscret_->gNode(gid);
					if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
					CNode* cmnode = static_cast<CNode*>(mnode);
					const int* mdofs = cmnode->Dofs();

					double mik = (mmap[0])[mdofs[0]];
					double mikold = (mmapold[0])[mdofs[0]];

					double valtxi1=0;
					double valtxi2=0;
					double valteta1=0;
					double valteta2=0;

					// compute linstick-matrix entry of the current active node / master node pair
					// loop over all derivative maps (=dimensions)
					for (int dim=0;dim<cnode->NumDof();++dim)
					{
						int col = cmnode->Dofs()[dim];
						valtxi1 = (ztxi+ct*jumptxi)/euclidean*(+1)*ct*txi[dim]*(mik-mikold)*ztxi;

						if(Dim()==3)
						{
							valteta1 = (ztxi+ct*jumptxi)/euclidean*(+1)*ct*txi[dim]*(mik-mikold)*zteta;
							valtxi2 = (zteta+ct*jumpteta)/euclidean*(+1)*ct*teta[dim]*(mik-mikold)*ztxi;
							valteta2 = (zteta+ct*jumpteta)/euclidean*(+1)*ct*teta[dim]*(mik-mikold)*zteta;
						}

						// do not assemble zeros into matrix
						if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
						if (Dim()==3)
						{
							if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
							if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
							if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
						}
					}
				}

				/*** 03 ********************** frcoeff*znor*ct*tangent.(D-Dn-1) ***/
				// loop over dimensions
				for (int dim=0;dim<cnode->NumDof();++dim)
				{
					double valtxi=0;
					double valteta=0;
					int col = cnode->Dofs()[dim];
					valtxi = (frcoeff*(znor-cn*wgap))*ct*txi[dim]*(D-Dold);

					if (Dim()==3)
						valteta = (frcoeff*(znor-cn*wgap))*ct*teta[dim]*(D-Dold);

					// do not assemble zeros into matrix
					if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
					if (Dim()==3)
						if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
					}

				/*** 04 ********************** frcoeff*znor*ct*tangent.(M-Mn-1) ***/

				// loop over all master nodes
				for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
				{
					double valtxi=0;
					double valteta=0;

					int gid = *mcurr;
					DRT::Node* mnode = idiscret_->gNode(gid);
					if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
					CNode* cmnode = static_cast<CNode*>(mnode);
					const int* mdofs = cmnode->Dofs();

					double mik = (mmap[0])[mdofs[0]];
					double mikold = (mmapold[0])[mdofs[0]];

					// loop over all derivative maps (=dimensions)
					for (int dim=0;dim<cnode->NumDof();++dim)
					{
						int col = cmnode->Dofs()[dim];
						valtxi = (frcoeff*(znor-cn*wgap))*(-1)*ct*txi[dim]*(mik-mikold);

						if (Dim()==3)
							valteta = (frcoeff*(znor-cn*wgap))*(-1)*ct*teta[dim]*(mik-mikold);

						// do not assemble zeros into matrix
						if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
						if(Dim()==3)
							if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
					 }
				 }

				// remaining terms only in case of full linearization
				if(fulllin)
				{
					/*** 1 ********************************* euclidean.deriv(T).z ***/
					// loop over dimensions
					for (int j=0;j<Dim();++j)
					{
						// loop over all entries of the current derivative map (txi)
						for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
						{
							int col = colcurr->first;
							double val = euclidean*(colcurr->second)*z[j];

							// do not assemble zeros into s matrix
							if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
						}

						if (Dim()==3)
						{
							// loop over all entries of the current derivative map (teta)
							for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
							{
								int col = colcurr->first;
								double val = euclidean*(colcurr->second)*z[j];

								// do not assemble zeros into s matrix
								if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
							}
						}
					}

					/*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
					// loop over dimensions
					for (int j=0;j<Dim();++j)
					{
						// loop over all entries of the current derivative map (txi)
						for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
						{
							int col = colcurr->first;
							double valtxi = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*ztxi;
							double valteta = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*zteta;

						 // do not assemble zeros into matrix
							if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
							if (Dim()==3)
								if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
						}

						if(Dim()==3)
						{
							// 3D loop over all entries of the current derivative map (teta)
							for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
							{
								int col = colcurr->first;
								double valtxi = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*ztxi;
								double valteta = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*zteta;

								// do not assemble zeros into matrix
								if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
								if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
							}
						}
					}

					/*** 3 ****************** deriv(euclidean).deriv(T).jump.ztan ***/

					// loop over dimensions
					for (int j=0;j<Dim();++j)
					{
						// loop over all entries of the current derivative map (txi)
						for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
						{
							int col = colcurr->first;
							double valtxi = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
							double valteta = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

							// do not assemble zeros into s matrix
							if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
							if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
						}

						if(Dim()==3)
						{
							// loop over all entries of the current derivative map (teta)
							for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
							{
								int col = colcurr->first;
								double valtxi = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
								double valteta = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

								// do not assemble zeros into matrix
								if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
								if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
							}
						}
					}

					/*** 4 ******************* deriv(euclidean).tan.deriv(D).ztan ***/
					// we need the dot product t*x of this node
					double txidotx = 0.0;
					double tetadotx = 0.0;
					for (int dim=0;dim<cnode->NumDof();++dim)
					{
						txidotx += txi[dim]*xi[dim];
						tetadotx += teta[dim]*xi[dim];
					}

					// prepare assembly
					map<int,double>& ddmap = cnode->GetDerivD()[gid];

				// loop over all entries of the current derivative map
					for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
					{
						int col = colcurr->first;
						double valtxi1 = (-1)*(ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*ztxi;
						double valteta1 = (-1)*(ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*zteta;
						double valtxi2 = (-1)*(zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*ztxi;
						double valteta2 = (-1)*(zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*zteta;

						// do not assemble zeros into matrix
						if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
						if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
						if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
						if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
					}

					/*** 5 ******************* deriv(euclidean).tan.deriv(M).ztan ***/
					// we need the Lin(M-matrix) entries of this node
					map<int,map<int,double> >& dmmap = cnode->GetDerivM();
					map<int,map<int,double> >::iterator dmcurr;

					// loop over all master nodes in the DerivM-map of the active slave node
					for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
					{
						int gid = dmcurr->first;
						DRT::Node* mnode = idiscret_->gNode(gid);
						if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
						CNode* cmnode = static_cast<CNode*>(mnode);
						double* mxi = cmnode->xspatial();

						// we need the dot product ns*xm of this node pair
						double txidotx = 0.0;
						double tetadotx = 0.0;
						for (int dim=0;dim<cnode->NumDof();++dim)
						{
							txidotx += txi[dim]*mxi[dim];
							tetadotx += teta[dim]*mxi[dim];
						}

						// compute entry of the current active node / master node pair
						map<int,double>& thisdmmap = cnode->GetDerivM(gid);

						// loop over all entries of the current derivative map
						for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
						{
							int col = colcurr->first;
							double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*ztxi;
							double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*txidotx*colcurr->second*zteta;
							double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*ztxi;
							double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*tetadotx*colcurr->second*zteta;

							// do not assemble zeros into matrix
							if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
							if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
							if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
							if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
						}
					}

					/*** 6 **************************** (frcoeff*znor).deriv(T).z ***/
					// loop over all dimensions
					for (int j=0;j<Dim();++j)
					{
						// loop over all entries of the current derivative map (txi)
						for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
						{
							int col = colcurr->first;
							double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];

							// do not assemble zeros into matrix
							if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
						}

						if(Dim()==3)
						{
							// loop over all entries of the current derivative map (teta)
							for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
							{
								int col = colcurr->first;
								double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];

								// do not assemble zeros into matrix
								if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
							}
						}
					}

					/*** 7 ************************* (frcoeff*znor).deriv(T).jump ***/
					// loop over all dimensions
					for (int j=0;j<Dim();++j)
					{
						// loop over all entries of the current derivative map (txi)
						for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
						{
							int col = colcurr->first;
							double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];

							// do not assemble zeros into matrix
							if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
						}

						if(Dim()==3)
						{
							// loop over all entries of the current derivative map (teta)
							for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
							{
								int col = colcurr->first;
								double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];

								// do not assemble zeros into s matrix
								if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
							}
						}
					}

					/*** 8 ********************* (frcoeff*znor).ct.tan.deriv(D).x ***/
					// we need the dot product t*x of this node
					txidotx = 0.0;
					tetadotx = 0.0;

					for (int dim=0;dim<cnode->NumDof();++dim)
					{
					 txidotx += txi[dim]*xi[dim];
					 tetadotx += teta[dim]*xi[dim];
					}

					// loop over all entries of the current derivative map
					for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
					{
						int col = colcurr->first;
						double valtxi = (-1)*(-1)*(frcoeff*(znor-cn*wgap))*ct*txidotx*colcurr->second;
						double valteta = (-1)*(-1)*(frcoeff*(znor-cn*wgap))*ct*tetadotx*colcurr->second;

						// do not assemble zeros into matrix
						if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

						if (Dim()==3)
						{
						 if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
						}
					}

					/*** 9 ********************* (frcoeff*znor).ct.tan.deriv(M).x ***/

					// loop over all master nodes in the DerivM-map of the active slave node
					for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
					{
						int gid = dmcurr->first;
						DRT::Node* mnode = idiscret_->gNode(gid);
						if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
						CNode* cmnode = static_cast<CNode*>(mnode);
						double* mxi = cmnode->xspatial();

						// we need the dot product t*xm of this node pair
						double txidotx = 0.0;
						double tetadotx = 0.0;
						for (int dim=0;dim<cnode->NumDof();++dim)
						{
							txidotx += txi[dim]*mxi[dim];
							tetadotx += teta[dim]*mxi[dim];
						}

						// compute entry of the current active node / master node pair
						map<int,double>& thisdmmap = cnode->GetDerivM(gid);

						 // loop over all entries of the current derivative map
						 for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
						 {
							 int col = colcurr->first;
							 double valtxi = (-1)*(frcoeff*(znor-cn*wgap))*ct*txidotx*colcurr->second;
							 double valteta = (-1)*(frcoeff*(znor-cn*wgap))*ct*tetadotx*colcurr->second;

							 // do not assemble zeros into matrix
							 if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
							 if (Dim()==3)
							 {
								 if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
							 }
						}
					}

					/*** 10 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
					// loop over all dimensions
					for (int j=0;j<Dim();++j)
					{
						// loop over all entries of the current derivative map
						for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
						{
							int col = colcurr->first;
							double valtxi = (-1)*(ztxi+ct*jumptxi)*frcoeff*(colcurr->second)*z[j];
							double valteta = (-1)*(zteta+ct*jumpteta)*frcoeff*(colcurr->second)*z[j];

							// do not assemble zeros into s matrix
							if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
							if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
						}
					}

					/*** 11 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
					// prepare assembly
					map<int,double>& dgmap = cnode->GetDerivG();

					// loop over all entries of the current derivative map
					for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
					{
						int col = colcurr->first;
						double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);
						double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);

						// do not assemble zeros into matrix
						if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
						if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
					}
				} // if fullin
    	} // if (frictionlessandfirst == false)
    } // loop over all slip nodes of the interface
#else
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      vector<map<int,double> > dnmap = cnode->GetDerivN();

      // iterator
      map<int,double>::iterator colcurr;

      vector <map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->jump();
      double* n = cnode->n();
      double* txi = cnode->txi();
      double* xi = cnode->xspatial();
      double* z = cnode->lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double znor    = n[0]*z[0] + n[1]*z[1];
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      /******************************************************************/
      // calculation of matrix entries of the linearized slip condition
      /******************************************************************/
      // 1) Entries from differentiation with respect to lagrange multipliers
      // 2) Entries on right hand side
      // 3) Entries from differentiation with respect to displacements

      // 1) Entries from differentiation with respect to lagrange multipliers
      /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frcoff*znor).tan ***/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (prefactor*ztan+sum-frcoeff*znor)*txi[dim]-frcoeff*(ztan+ct*jumptan)*n[dim];

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = txi[dim];
#endif

        // do not assemble zeros into matrix
        if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
      }

      // 2) Entries on right hand side
      /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

      // -C and remaining terms
      double value1 = -(abs(ztan+ct*jumptan))*ztan+(frcoeff*znor)*(ztan+ct*jumptan);
      double value2 = -(frcoeff*znor)*ztan-(frcoeff*znor)*(ztan+ct*jumptan);
      double value3 = +sum*ztan+prefactor*ztan*ztan;

      Epetra_SerialDenseVector rhsnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);

      rhsnode(0) = (value1+value2+value3);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) rhsnode(0) = 0;
#endif

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements

      /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

      // we need the nodal entries of the D-matrix and the old one
      double D= (cnode->GetD()[0])[cnode->Dofs()[0]];
      double Dold= (cnode->GetDOld()[0])[cnode->Dofs()[0]];

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
       //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif


       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

      // we need the nodal entries of the M-matrix and the old one
      vector<map<int,double> > mmap = cnode->GetM();
      vector<map<int,double> > mmapold = cnode->GetMOld();

      // create a set of nodes including nodes according to M entries
      // from current and previous time step
      set <int> mnodes;

      // iterator
      set<int>::iterator mcurr;

      set <int> mnodescurrent = cnode->GetMNodes();
    	set <int> mnodesold = cnode->GetMNodesOld();

    	for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
    	  mnodes.insert(*mcurr);

    	for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
    	  mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // compute linstick-matrix entry of the current active node / master node pair
        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
          //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      /********************************* frcoeff*znor*ct*tan.(D-Dn-1) ***/

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (frcoeff*znor)*ct*txi[dim]*(D-Dold);
        //cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -frcoeff*znor*ct*tan.(M-Mn-1).xm ***/

      // loop over all master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = (frcoeff*znor)*(-1)*ct*txi[dim]*(mik-mikold);
          //cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      // remaining terms only in case of full linearization
      if(fulllin)
      {
        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

          // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        map<int,double>& ddmap = cnode->GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        map<int,map<int,double> >& dmmap = cnode->GetDerivM();
        map<int,map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CNode* cmnode = static_cast<CNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          map<int,double>& thisdmmap = cnode->GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /*********************************** -(frcoeff*znor).DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*(frcoeff*znor)*(colcurr->second)*z[j];
            //cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /****************************** (frcoeff*znor).ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*(frcoeff*znor)*ct*(colcurr->second)*jump[j];
            //cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /****************************** +(frcoeff*znor).ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*(frcoeff*znor)*ct*tdotx*colcurr->second;
           //cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /***************************** -(frcoeff*znor).ct.T.DerivM.x ***/

          // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           CNode* cmnode = static_cast<CNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           map<int,double>& thisdmmap = cnode->GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*(frcoeff*znor)*ct*tdotx*colcurr->second;
            //cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /***************************** -frcoeff*DerivN.z(ztan+ct*utan) ***/

       // loop over all derivative maps (=dimensions)
       for (int j=0;j<mapsize;++j)
       {
         int k=0;

         // loop over all entries of the current derivative map
         for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(ztan+ct*jumptan)*frcoeff*(colcurr->second)*z[j];
           //cout << "10 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val =0;
#endif

           // do not assemble zeros into s matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
           ++k;
         }

         if (k!=colsize)
           dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }
      } // if fullin
    }
#endif
  } // Coulomb friction

  // Tresca Friction
  if (ftype == INPAR::CONTACT::friction_tresca)
  {
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* cnode = static_cast<CNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      vector<map<int,double> > dnmap = cnode->GetDerivN();

      // iterator
      map<int,double>::iterator colcurr;

      vector <map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->jump();
      double* txi = cnode->txi();
      double* xi = cnode->xspatial();
      double* z = cnode->lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      /******************************************************************/
      // calculation of matrix entries of the linearized slip condition
      /******************************************************************/
      // 1) Entries from differentiation with respect to lagrange multipliers
      // 2) Entries on right hand side
      // 3) Entries from differentiation with respect to displacements

      // 1) Entries from differentiation with respect to lagrange multipliers
      /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frbound).tan ***/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = (prefactor*ztan+sum-frbound)*txi[dim];

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = txi[dim];
#endif

        // do not assemble zeros into matrix
        if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
      }

      // 2) Entries on right hand side
      /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

      // -C and remaining terms
      double value1= -(abs(ztan+ct*jumptan))*ztan+frbound*(ztan+ct*jumptan);
      double value2= -frbound*ztan+sum*ztan+prefactor*ztan*ztan;

      Epetra_SerialDenseVector rhsnode(1);
      vector<int> lm(1);
      vector<int> lmowner(1);
      rhsnode(0) = (value1+value2);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) rhsnode(0) = 0;
#endif

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements

      /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

      // we need the nodal entries of the D-matrix and the old one
      double D= (cnode->GetD()[0])[cnode->Dofs()[0]];
      double Dold= (cnode->GetDOld()[0])[cnode->Dofs()[0]];

      if (abs(Dold)<0.0001)
        dserror ("Error:No entry for Dold");

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
       //cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << endl;
#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

       // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

      // we need the nodal entries of the M-matrix and the old one
      vector<map<int,double> > mmap = cnode->GetM();
      vector<map<int,double> > mmapold = cnode->GetMOld();

      // create a set of nodes including nodes according to M entries
      // from current and previous time step
      set <int> mnodes;

      // iterator
      set<int>::iterator mcurr;

      set <int> mnodescurrent = cnode->GetMNodes();
    	set <int> mnodesold = cnode->GetMNodesOld();

    	for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
    	  mnodes.insert(*mcurr);

    	for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
    	  mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this stick node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // compute linstick-matrix entry of the current active node / master node pair
        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
          //cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      /************************************** frbound*ct*tan.(D-Dn-1) ***/

      // loop over all derivative maps (=dimensions)
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double val = frbound*ct*txi[dim]*(D-Dold);
        //cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

        // do not assemble zeros into matrix
       if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
      }

      /********************************** -frbound*ct*tan.(M-Mn-1).xm ***/

      // loop over all master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();

        double mik = (mmap[0])[mdofs[0]];
        double mikold = (mmapold[0])[mdofs[0]];

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cmnode->Dofs()[dim];
          double val = frbound*(-1)*ct*txi[dim]*(mik-mikold);
          //cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->ActiveOld()==false) val = 0;
#endif
          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }
      }

      // remaining terms only in case of full linearization
      if(fulllin)
      {
        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
            if (cnode->ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        map<int,double>& ddmap = cnode->GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        map<int,map<int,double> >& dmmap = cnode->GetDerivM();
        map<int,map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CNode* cmnode = static_cast<CNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          map<int,double>& thisdmmap = cnode->GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /****************************************** -frbound.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*(colcurr->second)*z[j];
            //cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************ -frbound.ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*ct*(colcurr->second)*jump[j];
            //cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************* +frbound.ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*frbound*ct*tdotx*colcurr->second;
           //cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /********************************  -frbound.ct.T.DerivM.x ******/

          // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           CNode* cmnode = static_cast<CNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           map<int,double>& thisdmmap = cnode->GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*frbound*ct*tdotx*colcurr->second;
            //cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->ActiveOld()==false) val = 0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }
      } // if fullin
    }
  }// Tresca friction

  return;
}


/*----------------------------------------------------------------------*
 |  initialize active set (nodes / dofs)                      popp 03/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::InitializeActiveSet()
{
  // define local variables
  int countnodes = 0;
  int countdofs = 0;
  vector<int> mynodegids(snoderowmap_->NumMyElements());
  vector<int> mydofgids(sdofrowmap_->NumMyElements());

  // define local variables for slip maps
  vector<int> myslipnodegids(0);
  vector<int> myslipdofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // This is given by the CNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (cnode->IsInitActive())
    {
      // FULL: all slave nodes are assumed to be active
      cnode->Active()=true;
      mynodegids[countnodes] = cnode->Id();
      ++countnodes;

      for (int j=0;j<numdof;++j)
      {
        mydofgids[countdofs] = cnode->Dofs()[j];
        ++countdofs;
      }
    }
  }

  // resize the temporary vectors
  mynodegids.resize(countnodes);
  mydofgids.resize(countdofs);

  // communicate countnodes and countdofs among procs
  int gcountnodes, gcountdofs;
  Comm().SumAll(&countnodes,&gcountnodes,1);
  Comm().SumAll(&countdofs,&gcountdofs,1);

  // create active node map and active dof map
  activenodes_ = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));
  activedofs_  = rcp(new Epetra_Map(gcountdofs,countdofs,&mydofgids[0],0,Comm()));

  // create an empty slip node map and slip dof map
  // for the first time step (t=0) all nodes are stick nodes
  slipnodes_   = rcp(new Epetra_Map(0,0,Comm()));
  slipdofs_    = rcp(new Epetra_Map(0,0,Comm()));

  // split active dofs into Ndofs, Tdofs and slipTdofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::BuildActiveSet()
{
  // define local variables
  int countnodes = 0;
  int countdofs = 0;
  vector<int> mynodegids(snoderowmap_->NumMyElements());
  vector<int> mydofgids(sdofrowmap_->NumMyElements());

  // define local variables
  int countslipnodes = 0;
  int countslipdofs = 0;
  vector<int> myslipnodegids(snoderowmap_->NumMyElements());
  vector<int> myslipdofgids(sdofrowmap_->NumMyElements());

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    const int numdof = cnode->NumDof();

    // add node / dofs to temporary map IF active
    if (cnode->Active())
    {
      mynodegids[countnodes] = cnode->Id();
      ++countnodes;

      for (int j=0;j<numdof;++j)
      {
        mydofgids[countdofs] = cnode->Dofs()[j];
        ++countdofs;
      }
    }

    // add node / dofs to temporary map IF slip
    if (cnode->Slip())
    {
      myslipnodegids[countslipnodes] = cnode->Id();
      ++countslipnodes;

      for (int j=0;j<numdof;++j)
      {
        myslipdofgids[countslipdofs] = cnode->Dofs()[j];
        ++countslipdofs;
       }
     }
  }

  // resize the temporary vectors
  mynodegids.resize(countnodes);
  mydofgids.resize(countdofs);
  myslipnodegids.resize(countslipnodes);
  myslipdofgids.resize(countslipdofs);

  // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
  int gcountnodes, gcountdofs, gcountslipnodes,gcountslipdofs;
  Comm().SumAll(&countnodes,&gcountnodes,1);
  Comm().SumAll(&countdofs,&gcountdofs,1);
  Comm().SumAll(&countslipnodes,&gcountslipnodes,1);
  Comm().SumAll(&countslipdofs,&gcountslipdofs,1);

  // create active node map and active dof map
  activenodes_ = rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,Comm()));
  activedofs_  = rcp(new Epetra_Map(gcountdofs,countdofs,&mydofgids[0],0,Comm()));
  slipnodes_ = rcp(new Epetra_Map(gcountslipnodes,countslipnodes,&myslipnodegids[0],0,Comm()));
  slipdofs_  = rcp(new Epetra_Map(gcountslipdofs,countslipdofs,&myslipdofgids[0],0,Comm()));

  // split active dofs into Ndofs and Tdofs and slip dofs in Tslipdofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  split active dofs into Ndofs, Tdofs and slipTdofs          popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::SplitActiveDofs()
{
  // get out of here if active set is empty
  if (activenodes_==null)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  else if (activenodes_->NumGlobalElements()==0)
  {
    activen_ = rcp(new Epetra_Map(0,0,Comm()));
    activet_ = rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countN=0;
  int countT=0;
  vector<int> myNgids(activenodes_->NumMyElements());
  vector<int> myTgids((Dim()-1)*activenodes_->NumMyElements());

  // dimension check
  double dimcheck =(activedofs_->NumGlobalElements())/(activenodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all active row nodes
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    const int numdof = cnode->NumDof();

    // add first dof to Nmap
    myNgids[countN] = cnode->Dofs()[0];
    ++countN;

    // add remaining dofs to Tmap
    for (int j=1;j<numdof;++j)
    {
      myTgids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNgids.resize(countN);
  myTgids.resize(countT);

  // communicate countN and countT among procs
  int gcountN, gcountT;
  Comm().SumAll(&countN,&gcountN,1);
  Comm().SumAll(&countT,&gcountT,1);

  // check global dimensions
  if ((gcountN+gcountT)!=activedofs_->NumGlobalElements())
    dserror("ERROR: SplitActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = rcp(new Epetra_Map(gcountN,countN,&myNgids[0],0,Comm()));
  activet_ = rcp(new Epetra_Map(gcountT,countT,&myTgids[0],0,Comm()));

  // *******************************************************************
  // EXTRACTING TANGENTIAL DOFS FROM SLIP DOFS
  // *******************************************************************

  // get out of here if slip set is empty
  if (slipnodes_==null)
  {
  	slipt_ = rcp(new Epetra_Map(0,0,Comm()));
  	return true;
  }

  if (slipnodes_->NumGlobalElements()==0)
  {
  	slipt_ = rcp(new Epetra_Map(0,0,Comm()));
  	return true;
  }

  // define local variables
  int countslipT=0;
  vector<int> myslipTgids((Dim()-1)*slipnodes_->NumMyElements());

  // dimension check
  dimcheck =(slipdofs_->NumGlobalElements())/(slipnodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slip row nodes
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    const int numdof = cnode->NumDof();

    // add dofs to slipTmap
    for (int j=1;j<numdof;++j)
    {
      myslipTgids[countslipT] = cnode->Dofs()[j];
      ++countslipT;
    }
  }

  // resize the temporary vectors
  myslipTgids.resize(countslipT);

  // communicate countslipT among procs
  int gcountslipT;
  Comm().SumAll(&countslipT,&gcountslipT,1);

  // create Tslipmap objects
  slipt_   = rcp(new Epetra_Map(gcountslipT,countslipT,&myslipTgids[0],0,Comm()));

  return true;
}

#endif  // #ifdef CCADISCRET
