/*!----------------------------------------------------------------------
\file mortar_interface.cpp
\brief One mortar coupling interface

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef PARALLEL
#include "Epetra_SerialComm.h"
#endif

#include "mortar_interface.H"
#include "mortar_integrator.H"
#include "mortar_coupling2d.H"
#include "mortar_coupling3d.H"
#include "mortar_dofset.H"
#include "mortar_binarytree.H"
#include "mortar_defines.H"
#include "../drt_contact/contactdefines.H"
#include "../drt_lib/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarInterface::MortarInterface(const int id, const Epetra_Comm& comm,
                             const int dim, const Teuchos::ParameterList& icontact) :
id_(id),
comm_(comm),
dim_(dim),
icontact_(icontact),
shapefcn_(MortarInterface::Undefined)
{
  RCP<Epetra_Comm> com = rcp(Comm().Clone());
  if (Dim()!=2 && Dim()!=3) dserror("ERROR: Mortar problem must be 2D or 3D");
  procmap_.clear();
  idiscret_ = rcp(new DRT::Discretization((string)"Mortar Interface",com));

  // overwrite shape function type
  INPAR::MORTAR::ShapeFcn shapefcn = Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(IParams(),"SHAPEFCN");
  if (shapefcn == INPAR::MORTAR::shape_dual)
   shapefcn_ = MortarInterface::DualFunctions;
  else if (shapefcn == INPAR::MORTAR::shape_standard)
    shapefcn_ = MortarInterface::StandardFunctions;

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MORTAR::MortarInterface& interface)
{
  interface.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "Mortar Interface Id " << id_ << endl;
    os << "Mortar Interface Discretization:" << endl;
  }
  os << Discret();
  return;
}

/*----------------------------------------------------------------------*
 |  finalize construction of interface (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::FillComplete()
{
  // we'd like to call idiscret_.FillComplete(true,false,false) but this
  // will assign all nodes new degrees of freedom which we don't want.
  // We would like to use the degrees of freedom that were stored in the
  // mortar nodes. To do so, we have to create and set our own
  // version of a DofSet class before we call FillComplete on the
  // interface discretization.
  // Our special dofset class will not assign new dofs but will assign the
  // dofs stored in the nodes.
  {
    RCP<MORTAR::MortarDofSet> mrtrdofset = rcp(new MORTAR::MortarDofSet());
    Discret().ReplaceDofSet(mrtrdofset);
    // do not assign dofs yet, we'll do this below after
    // shuffling around of nodes and elements (saves time)
    Discret().FillComplete(false,false,false);
  }

#ifdef MORTARBOUNDMOD
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
    MORTAR::MortarNode* node = static_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

    // candidates are slave nodes with only 1 adjacent MortarElement
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
#endif // #ifdef MORTARBOUNDMOD

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
  // column map. This way, also mortar elements will be fully ghosted on all
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
  if (SearchAlg()==INPAR::MORTAR::search_binarytree)
  {
    // create binary tree object for search and setup tree
	  binarytree_ = rcp(new MORTAR::BinaryTree(Discret(),selecolmap_,melefullmap_,Dim(),SearchParam()));
  }
	  
  return;
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                 popp 11/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveSets()
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
      bool isslave = dynamic_cast<MORTAR::MortarNode*>(Discret().gNode(gid))->IsSlave();
      bool isonbound = dynamic_cast<MORTAR::MortarNode*>(Discret().gNode(gid))->IsOnBound();
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
      bool isslave = dynamic_cast<MORTAR::MortarElement*>(Discret().gElement(gid))->IsSlave();
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
      MortarNode* mrtrnode = static_cast<MortarNode*>(node);
      bool isslave = mrtrnode->IsSlave();

      if (oldnodecolmap_->MyGID(gid))
      {
        if (isslave)
          for (int j=0;j<mrtrnode->NumDof();++j)
            sc.push_back(mrtrnode->Dofs()[j]);
        else
          for (int j=0;j<mrtrnode->NumDof();++j)
            mc.push_back(mrtrnode->Dofs()[j]);
      }

      if (isslave)
        for (int j=0;j<mrtrnode->NumDof();++j)
          scfull.push_back(mrtrnode->Dofs()[j]);
      else
        for (int j=0;j<mrtrnode->NumDof();++j)
          mcfull.push_back(mrtrnode->Dofs()[j]);

      if (!noderowmap->MyGID(gid)) continue;

      if (isslave)
        for (int j=0;j<mrtrnode->NumDof();++j)
          sr.push_back(mrtrnode->Dofs()[j]);
      else
        for (int j=0;j<mrtrnode->NumDof();++j)
          mr.push_back(mrtrnode->Dofs()[j]);
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
 |  initialize / reset mortar interface                       popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    MORTAR::MortarNode* node = static_cast<MORTAR::MortarNode*>(idiscret_->lColNode(i));

    //reset nodal normal and tangents and jumps
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
    }

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
    
    // reset feasible projection status
    node->HasProj() = false;
  }

  // loop over all elements to reset candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    MORTAR::MortarElement* element = static_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  set current and old deformation state                      popp 12/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::SetState(const string& statename, const RCP<Epetra_Vector> vec)
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
      MORTAR::MortarNode* node = static_cast<MORTAR::MortarNode*>(idiscret_->lColNode(i));
      const int numdof = node->NumDof();
      vector<double> mydisp(numdof);
      vector<int> lm(numdof);

      for (int j=0;j<numdof;++j)
        lm[j]=node->Dofs()[j];

      DRT::UTILS::ExtractMyValues(global,mydisp,lm);

      // add mydisp[2]=0 for 2D problems
      if (mydisp.size()<3)
        mydisp.resize(3);

      // set current configuration
      for (int j=0;j<3;++j)
        node->xspatial()[j]=node->X()[j]+mydisp[j];
    }

    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      MORTAR::MortarElement* element = static_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
      element->Area()=element->ComputeArea();
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate mortar coupling (public)                         popp 11/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Evaluate(bool nonlinear)
{
  // interface needs to be complete
  if (!Filled() && Comm().MyPID()==0)
      dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm()) return;

  //**********************************************************************
  // search algorithm
  //**********************************************************************
  //lComm()->Barrier();
  //const double t_start = ds_cputime();

  if (SearchAlg()==INPAR::MORTAR::search_bfnode)          EvaluateSearch();
  else if (SearchAlg()==INPAR::MORTAR::search_bfele)      EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg()==INPAR::MORTAR::search_binarytree) EvaluateSearchBinarytree();
  else                                                     dserror("ERROR: Invalid search algorithm");
  
  // loop over proc's slave nodes of the interface
  // use standard column map to include processor's ghosted nodes
  // use boundary map to include slave side boundary nodes
  for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    int gid = snodecolmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MortarNode* mrtrnode = static_cast<MortarNode*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();
  }
    
  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    MortarElement* selement = static_cast<MortarElement*>(ele1);

    //********************************************************************
    // 1) integrate Mortar matrix D (lives on slave side only!)
    // 2) compute directional derivative of D and store into nodes
    //********************************************************************
    // we only do this on a slave element basis in the case with
    // two separate Mortar loops for D and M. If MORTARONELOOP
    // is applied, we combine the construction of D with the construction
    // of M and linearization of D with the linearization of M in
    // Integrator::IntegrateDerivSegment2D()!
    //********************************************************************
#ifndef MORTARONELOOP
    IntegrateSlave(*selement,nonlinear);
#endif // #ifndef MORTARONELOOP

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      MortarElement* melement = static_cast<MortarElement*>(ele2);

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pair)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      IntegrateCoupling(*selement,*melement,nonlinear);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Search node-based ("brute force") (public)                popp 01/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::EvaluateSearch()
{
  /**********************************************************************/
  /* SEARCH ALGORITHM:                                                 */
  /* The idea of the search is to reduce the number of master / slave   */
  /* element pairs that are checked for overlap and coupling by intro-  */
  /* ducing information about proximity and maybe history!              */
  /* This old version is still brute force for finding the closest      */
  /* node to each node, so it has been replaced by a more               */
  /* sophisticated approach (bounding vol. hierarchies / binary tree).  */
  /**********************************************************************/

  // loop over proc's slave nodes for closest node detection
  // use standard column map to include processor's ghosted nodes
  // use boundary map to include slave boundary nodes
  for (int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    int gid = snodecolmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    MortarNode* snode = static_cast<MortarNode*>(node);

    // find closest master node to current slave node
    double mindist = 1.0e12;
    MortarNode* closestnode = snode->FindClosestNode(idiscret_,mnodefullmapnobound_,mindist);
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
        MortarElement* selement = static_cast<MortarElement*> (adjslave[j]);
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
    MortarNode* mnode = static_cast<MortarNode*>(node);

    // find closest slave node to current master node
    double mindist = 1.0e12;
    MortarNode* closestnode = mnode->FindClosestNode(idiscret_,snodefullmapbound_,mindist);
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
        MortarElement* selement = static_cast<MortarElement*> (adjclosest[j]);
        selement->AddSearchElements(mids);
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Search element-based "brute force" (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateSearchBruteForce(const double& eps)
{
  // calculate minimal element length
  double lmin = 1.0e12;
  double enlarge = 0.0;

  // loop over all slave elements on this proc.
  for (int i=0;i<selecolmap_->NumMyElements();++i)
  {
    DRT::Element* element = idiscret_->gElement(selecolmap_->GID(0));
    if (!element) dserror("ERROR: Cannot find element with gid %\n",selecolmap_->GID(0));
    MORTAR::MortarElement* mrtrelement = (MortarElement*) element;
    if (mrtrelement->MinEdgeSize() < lmin)
      lmin= mrtrelement->MinEdgeSize();
  }

  // loop over all master elements on this proc.
  for (int i=0;i<melefullmap_->NumMyElements();++i)
  {
    DRT::Element* element = idiscret_->gElement(melefullmap_->GID(i));
    if (!element) dserror("ERROR: Cannot find element with gid %\n",melefullmap_->GID(i));
    MORTAR::MortarElement* mrtrelement = (MortarElement*) element;
    if (mrtrelement->MinEdgeSize() < lmin)
      lmin= mrtrelement->MinEdgeSize();
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
  // perform brute-force search (element-based)
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
    MortarNode* mrtrnode=static_cast<MortarNode*>(node[0]);
    const double* posnode = mrtrnode->xspatial();

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
      MortarNode* mrtrnode=static_cast<MortarNode*>(node[j]);
      posnode = mrtrnode->xspatial();

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
      MortarNode* mrtrnode=static_cast<MortarNode*>(node[j]);

      double auxpos [3] = {0.0, 0.0, 0.0};
      double scalar=0.0;
      for (int k=0; k<dim_;k++)
        scalar+=(mrtrnode->X()[k]+mrtrnode->uold()[k]-mrtrnode->xspatial()[k])*mrtrnode->n()[k];
      for (int k=0;k<dim_;k++)
        auxpos[k]= mrtrnode->xspatial()[k]+scalar*mrtrnode->n()[k];

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
      MortarNode* mrtrnode=static_cast<MortarNode*>(node[0]);
      const double* posnode = mrtrnode->xspatial();

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
        MortarNode* mrtrnode=static_cast<MortarNode*>(node[k]);
        posnode = mrtrnode->xspatial();

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
        //cout << Comm().MyPID() << "\nCoupling found between slave element: "
        //     << sgid <<" and master element: "<< mgid;
        DRT::Element* element= idiscret_->gElement(sgid);
        MORTAR::MortarElement* selement = static_cast<MORTAR::MortarElement*>(element);
        selement->AddSearchElements(mgid);
      }
    } // for all master elements
  } // for all slave elements

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially coupling sl/ma pairs (public)      popp 10/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::EvaluateSearchBinarytree()
{
  // get out of here if not participating in interface
  if (!lComm()) return true;

  // *********************************************************************
  // Possible versions for general mortar setting:
  // *********************************************************************
  //
  // 1) Combined Update and Search
  // -> In this case we only have to call SearchCombined(), which
  //    does buth top-down update (where necessary) and search.
  //
  // 2) Separate Update and Search
  // -> In this case we have to explicitly call and updating routine, i.e.
  //    UpdateTreeTopDown() or UpdateTreeBottomUp() before calling the
  //    search routine SearchSeparate(). Of course, the bottom-up
  //    update makes more sense here!
  //
  // *********************************************************************
  
  // calculate minimal element length
  binarytree_->SetEnlarge(false);
    
  // update tree in a top down way
  //binarytree_->UpdateTreeTopDown();

  // update tree in a bottom up way
  //binarytree_->UpdateTreeBottomUp();

#ifdef CONTACTGMSHCTN
  for (int i=0;i<(int)(binarytree_->CouplingMap().size());i++)
    binarytree_->CouplingMap()[i].clear();
  binarytree_->CouplingMap().clear();
  binarytree_->CouplingMap().resize(2);
#endif //CONTACTGMSHCTN

  // search with a separate algorithm
  //binarytree_->SearchSeparate();

  // search with an combined algorithm
  binarytree_->SearchCombined();
  
	return true;
}

/*----------------------------------------------------------------------*
 |  Integrate Mortar matrix D on slave element (public)       popp 01/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::IntegrateSlave(MORTAR::MortarElement& sele,
                                             bool nonlinear)
{
  // create an integrator instance with correct NumGP and Dim
  MORTAR::Integrator integrator(shapefcn_,sele.Shape(),nonlinear);

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
bool MORTAR::MortarInterface::IntegrateCoupling(MORTAR::MortarElement& sele,
                                           MORTAR::MortarElement& mele,
                                           bool nonlinear)
{
  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (Dim()==2)
    MORTAR::Coupling2d coup(shapefcn_,Discret(),Dim(),sele,mele,nonlinear);
  // ************************************************************** 3D ***
  else if (Dim()==3)
  {
    bool auxplane = Teuchos::getIntegralValue<int>(IParams(),"COUPLING_AUXPLANE");

    // ************************************************** quadratic 3D ***
    if (sele.IsQuad() || mele.IsQuad())
    {
      // only for auxiliary plane 3D version
      if (!auxplane) dserror("ERROR: Quadratic 3D coupling only for AuxPlane case!");

      // build linear integration elements from quadratic MortarElements
      vector<RCP<MORTAR::IntElement> > sauxelements(0);
      vector<RCP<MORTAR::IntElement> > mauxelements(0);
      SplitIntElements(sele,sauxelements);
      SplitIntElements(mele,mauxelements);

      // loop over all IntElement pairs for coupling
      for (int i=0;i<(int)sauxelements.size();++i)
      {
        for (int j=0;j<(int)mauxelements.size();++j)
        {
          // create instance of coupling class
          MORTAR::Coupling3dQuad coup(shapefcn_,Discret(),Dim(),true,auxplane,
                        sele,mele,*sauxelements[i],*mauxelements[j],nonlinear);
          // do coupling
          coup.EvaluateCoupling();
        }
      }
    }

    // ***************************************************** linear 3D ***
    else
    {
      // create instance of coupling class
      MORTAR::Coupling3d coup(shapefcn_,Discret(),Dim(),false,auxplane,
                              sele,mele,nonlinear);
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
 | Split MortarElements->IntElements for 3D quad. coupling    popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::SplitIntElements(MORTAR::MortarElement& ele,
                       vector<RCP<MORTAR::IntElement> >& auxele)
{
  // *********************************************************************
  // do splitting for given element
  // *********************************************************** quad9 ***
  if (ele.Shape()==DRT::Element::quad9)
  {
    //dserror("ERROR: Quadratic 3D coupling for quad9 under construction...");

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
    //dserror("ERROR: Quadratic 3D coupling for quad8 under construction...");

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
    //dserror("ERROR: Quadratic 3D coupling for tri6 under construction...");

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
 |  Assemble geometry-dependent lagrange multipliers (global)      popp 05/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleLM(Epetra_Vector& zglobal)
{
  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    MortarNode* mrtrnode = static_cast<MortarNode*>(node);

    int dim = mrtrnode->NumDof();
    double* lm = mrtrnode->lm();

    Epetra_SerialDenseVector lmnode(dim);
    vector<int> lmdof(dim);
    vector<int> lmowner(dim);

    for( int k=0; k<dim; ++k )
    {
      lmnode(k) = lm[k];
      lmdof[k] = mrtrnode->Dofs()[k];
      lmowner[k] = mrtrnode->Owner();
    }

    // do assembly
    LINALG::Assemble(zglobal, lmnode, lmdof, lmowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar matrices                                  popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleDM(LINALG::SparseMatrix& dglobal,
                                      LINALG::SparseMatrix& mglobal)
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
    MortarNode* mrtrnode = static_cast<MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDM: Node ownership inconsistency!");

    /**************************************************** D-matrix ******/
    if ((mrtrnode->GetD()).size()>0)
    {
      vector<map<int,double> > dmap = mrtrnode->GetD();
      int rowsize = mrtrnode->NumDof();
      int colsize = (int)dmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)dmap[j].size() != (int)dmap[j+1].size())
          dserror("ERROR: AssembleDM: Column dim. of nodal D-map is inconsistent!");

      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = mrtrnode->Dofs()[j];
        int k = 0;

        for (colcurr=dmap[j].begin();colcurr!=dmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;

          // do the assembly into global D matrix
          if (shapefcn_ == MortarInterface::DualFunctions)
          {
            // check for diagonality
            if (row!=col && abs(val)>1.0e-12)
              dserror("ERROR: AssembleDM: D-Matrix is not diagonal!");

            // create an explicitly diagonal d matrix
            if (row==col)
              dglobal.Assemble(val, row, col);
          }
          else if (shapefcn_ == MortarInterface::StandardFunctions)
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
          dserror("ERROR: AssembleDM: k = %i but colsize = %i",k,colsize);
      }
    }

    /**************************************************** M-matrix ******/
    if ((mrtrnode->GetM()).size()>0)
    {
      vector<map<int,double> > mmap = mrtrnode->GetM();
      int rowsize = mrtrnode->NumDof();
      int colsize = (int)mmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)mmap[j].size() != (int)mmap[j+1].size())
          dserror("ERROR: AssembleDM: Column dim. of nodal M-map is inconsistent!");

      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = mrtrnode->Dofs()[j];
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
          dserror("ERROR: AssembleDM: k = %i but colsize = %i",k,colsize);
      }
    }

    /************************************************* Mmod-matrix ******/
    if ((mrtrnode->GetMmod()).size()>0)
    {
      vector<map<int,double> > mmap = mrtrnode->GetMmod();
      int rowsize = mrtrnode->NumDof();
      int colsize = (int)mmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)mmap[j].size() != (int)mmap[j+1].size())
          dserror("ERROR: AssembleDM: Column dim. of nodal Mmod-map is inconsistent!");

      Epetra_SerialDenseMatrix Mnode(rowsize,colsize);
      vector<int> lmrow(rowsize);
      vector<int> lmcol(colsize);
      vector<int> lmrowowner(rowsize);
      map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = mrtrnode->Dofs()[j];
        int k = 0;
        lmrow[j] = row;
        lmrowowner[j] = mrtrnode->Owner();

        for (colcurr=mmap[j].begin();colcurr!=mmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;
          lmcol[k] = col;

          Mnode(j,k)=val;
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleDM: k = %i but colsize = %i",k,colsize);
      }

      mglobal.Assemble(-1,Mnode,lmrow,lmrowowner,lmcol);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble weighted gap only                                popp 12/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleG(Epetra_Vector& gglobal)
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
    MortarNode* mrtrnode = static_cast<MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleG: Node ownership inconsistency!");

    /**************************************************** g-vector ******/
    // weighted gap
    vector<double> gap(Dim());
    
    // slave contribution (D)
    vector<map<int,double> > dmap = mrtrnode->GetD();
    map<int,double>::iterator dcurr;
    
    for (int d=0;d<snodefullmap_->NumMyElements();++d)
    {
      int gid = snodefullmap_->GID(d);
      DRT::Node* snode = idiscret_->gNode(gid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
      MortarNode* mrtrsnode = static_cast<MortarNode*>(snode);
      const int* sdofs = mrtrsnode->Dofs();
      bool hasentry = false;

      // look for this master node in D-map of the active slave node
      for (dcurr=dmap[0].begin();dcurr!=dmap[0].end();++dcurr)
        if ((dcurr->first)==sdofs[0])
        {
          hasentry=true;
          break;
        }

      double dik = (dmap[0])[sdofs[0]];

      // get out of here, if slave node not adjacent or coupling very weak
      if (!hasentry || abs(dik)<1.0e-12) continue;
      for (int j=0;j<(int)gap.size();++j) gap[j] += dik * (mrtrsnode->X()[j]);
    }

    // master contribution (M)
    vector<map<int,double> > mmap = mrtrnode->GetM();
    map<int,double>::iterator mcurr;

    for (int m=0;m<mnodefullmap_->NumMyElements();++m)
    {
      int gid = mnodefullmap_->GID(m);
      DRT::Node* mnode = idiscret_->gNode(gid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
      MortarNode* mrtrmnode = static_cast<MortarNode*>(mnode);
      const int* mdofs = mrtrmnode->Dofs();
      bool hasentry = false;

      // look for this master node in M-map of the active slave node
      for (mcurr=mmap[0].begin();mcurr!=mmap[0].end();++mcurr)
        if ((mcurr->first)==mdofs[0])
        {
          hasentry=true;
          break;
        }

      double mik = (mmap[0])[mdofs[0]];

      // get out of here, if master node not adjacent or coupling very weak
      if (!hasentry || abs(mik)<1.0e-12) continue;
      for (int j=0;j<(int)gap.size();++j) gap[j] -= mik * (mrtrmnode->X()[j]);
    }
    
    // check if constraints fulfilled in reference configuration
    //if (i==0) cout << endl;
    //printf("NODE: %d \t gap[0]: %e \t gap[1]: %e \t gap[2]: %e \n",gid,gap[0],gap[1],gap[2]);
    //fflush(stdout);
    //if (abs(gap[0])>1.0e-12 || abs(gap[1])>1.0e-12 || abs(gap[2])>1.0e-12)
    //  cout << "***WARNING*** Non-zero initial constraint condition!" << endl;
    
    // prepare assembly
    Epetra_SerialDenseVector gnode(Dim());
    vector<int> lm(Dim());
    vector<int> lmowner(Dim());

    for (int j=0; j< Dim(); ++j)
    {
      gnode(j) = gap[j];
      lm[j] = mrtrnode->Dofs()[j];
      lmowner[j] = mrtrnode->Owner();
    }
    
    // do assembly
    LINALG::Assemble(gglobal,gnode,lm,lmowner);
  }

  return;
}

#endif  // #ifdef CCADISCRET
