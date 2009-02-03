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
#include "drt_contact_coupling.H"
#include "drt_cdofset.H"
#include "contactdefines.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Interface::Interface(const int id, const Epetra_Comm& comm, const int dim) :
id_(id),
comm_(comm),
dim_(dim)
{
  RCP<Epetra_Comm> com = rcp(Comm().Clone());
  if (Dim()!=2 && Dim()!=3) dserror("ERROR: Contact problem must be 2D or 3D");
  procmap_.clear();
  idiscret_ = rcp(new DRT::Discretization((string)"Contact Interface",com));
  contactsegs_.Reshape(0,0);
  
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

  // get standard nodal column map
  RCP<Epetra_Map> oldnodecolmap = rcp (new Epetra_Map(*(Discret().NodeColMap() )));
  // get standard element column map
  RCP<Epetra_Map> oldelecolmap = rcp (new Epetra_Map(*(Discret().ElementColMap() )));

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
      if (oldnodecolmap->MyGID(gid))
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
      if (oldelecolmap->MyGID(gid))
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
      
      if (oldnodecolmap->MyGID(gid))
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

  // create binarytree object for contact search
  // binarytree_ = rcp(new CONTACT::BinaryTree(Discret(),selecolmap_,melefullmap_,Dim()));
   
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

    //reset nodal normal and tangents
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
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
        
    // reset closest node
    // (FIXME: at the moment we do not need this info. in the next
    // iteration, but it might be helpful for accelerated search!!!)
    node->ClosestNode() = -1;

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
    
    // reset nodal weighted gap
    node->Getg() = 1.0e12;

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
 
#ifdef DEBUG
  // Visualize nodal normal field
  // VisualizeGmshLight();
#endif // #ifdef DEBUG
  
  // contact search algorithm
  //lComm()->Barrier();
  //const double t_start = ds_cputime();
  EvaluateContactSearch();
  //lComm()->Barrier();
  //const double t_end = ds_cputime()-t_start;
  //if (lComm()->MyPID()==0)
  // {
  //  cout << "Brute-force search: " << t_end << " seconds\n";
  //}
  
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
    // of M in the method Integrator::AssembleM() and linearization of D
    // with the linearization of M in the method Integrator::DerivM()! 
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
      IntegrateCoupling(*selement,*melement);
    }
  }

#ifdef CONTACTFDMORTARD
  // FD check of Mortar matrix D derivatives
  FDCheckMortarDDeriv();
#endif // #ifdef CONTACTFDMORTARD
  
#ifdef CONTACTFDMORTARM
  // FD check of Mortar matrix M derivatives
  FDCheckMortarMDeriv();
#endif // #ifdef CONTACTFDMORTARM
  
  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::EvaluateContactSearch()
{
  /**********************************************************************/
  /* CONTACT SEARCH ALGORITHM:                                          */
  /* The idea of the search is to reduce the number of master / slave   */
  /* element pairs that are checked for overlap and contact by intro-   */
  /* ducing information about proximity and maybe history!              */
  /* At the moment this is still brute force for finding the closest    */
  /* CNode to each CNode, so it will have to replaced by a more         */
  /* sophisticated approach in the future (bounding vol. hierarchies?)  */
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
    if (mindist<=CONTACTCRITDIST)
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
    if (mindist<=CONTACTCRITDIST)
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
 |  Integrate Mortar matrix D on slave element (public)       popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateSlave(CONTACT::CElement& sele)
{
  // create an integrator instance with correct NumGP and Dim
  CONTACT::Integrator integrator(sele.Shape());
  
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

  // do the integration
  RCP<Epetra_SerialDenseMatrix> dseg = integrator.IntegrateD(sele,sxia,sxib);
  
  // do the linearization of D
  integrator.DerivD(sele,sxia,sxib);
  
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
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  CONTACT::Coupling coup(Discret(),sele,mele,Dim(),CSegs());
      
  return true;
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

          // check for diagonality
          if (row!=col && abs(val)>1.0e-12)
            dserror("ERROR: AssembleDMG: D-Matrix is not diagonal!");

          // create an explicitly diagonal d matrix
          if (row==col) dglobal.Assemble(val,row,col);
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
      
      // check if this node has a feasible projection
      // else, it cannot be in contact and weighted gap should be positive
      // (otherwise wrong results possible for g~ because of non-positivity
      // of dual shape functions!!!)
      //cout << "Node ID: " << cnode->Id() << " HasProj: " << cnode->HasProj() << endl;
      if (!cnode->HasProj()) gap = 1.0e12;

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
 |  Assemble matrix L / vector R for Tresca friction          mgit 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleTresca(LINALG::SparseMatrix& lglobal,
                                        Epetra_Vector& rglobal,
                                        double& frbound, double& ct)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // loop over all active slip nodes of the interface
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleTresca: Node ownership inconsistency!");
    
    // prepare assembly (only 2D so far !!!!)
    int colsize = 1;
    vector<int> lmrowT(1);
    vector<int> lmrowownerT(1);
    vector<int> lmcol(colsize);

    lmrowT[0] = slipt_->GID(i);
    lmrowownerT[0] = cnode->Owner();

    if (colsize==3)
      dserror("ERROR: AssembleTresca: 3D case not yet implemented!");

    /******************************** L-matrix and R-vector**************/
    // the L-Matrix is calculated according to Huebner, Stadler, Wohlmuth: 
    // A primal-dual active set algorithm for three dimensional contact 
    // problems with coulomb friction, 2008
    
    // we need the tangent vector and the vector of lagrange multipliers 
    double tangent[3];
    tangent[0] = cnode->txi()[0];
    tangent[1] = cnode->txi()[1];
    tangent[2] = 0.0;
    
    double z[3];
    z[0] = cnode->lm()[0];
    z[1] = cnode->lm()[1];
    z[2] =  0.0;
    
    // we get the incremental, relative displacements (incremental dis-
    // placement is equivalent to velocity in case of linear approach of 
    // time integration)
    double jump[3];
    jump[0] = cnode->jump()[0];
    jump[1] = cnode->jump()[1];
    jump[2] =  0.0;
    
//    // FIXGIT: The following the calculation of the jump is done 
      // based on quantities stored on nodes. But this is done globally 
      // in contact manager wherby following lines will be deleted soon. 
//        
//    // get nodal entries of mortar matrices D and M, in fact current and
//    // old ones
//    vector<map<int,double> > dmap = cnode->GetD();
//  	vector<map<int,double> > mmap = cnode->GetM();
//  	vector<map<int,double> > dmapold = cnode->GetDOld();
//  	vector<map<int,double> > mmapold = cnode->GetMOld();
//  	
//  	// get vectors of ids from master nodes which have entries in the M 
//  	// matrix, in fact current and old ones
//  	vector<int> mnodesvector = cnode->GetMNodes();
//  	vector<int> mnodesvectorold = cnode->GetMNodesOld();
//  	
//  	// create a merged vector which includes the ids of master nodes with
//  	// entries in M from current and last step 
//  	int mnodessize = mnodesvector.size();
//  	int mnodessizeold = mnodesvectorold.size();
//  	
//  	for(int i=0;i<mnodessizeold;++i)
//  	{
//  		bool isentry = false;
//  		for(int j=0;j<mnodessize;j++)
//  		{
//  			if (mnodesvector[j] == mnodesvectorold[i]) 
//  			{
//  			  isentry = true;
//  			  break;
//  			}
//  		}
//  		
//  		// add nodeid to vector mnodes_
//  		if (isentry == false) mnodesvector.push_back(mnodesvectorold[i]);
//  	} // end of creating merged vector
//  	
//  	// get row- and column size
//  	// mnodes now also contains old ones
//  	int rowsize1 = cnode->NumDof();
//  	mnodessize = mnodesvector.size();
//
//  	// loop over rows (degrees of freedom)
//  	for (int j=0;j<rowsize1;++j)
//  	{
//  	  // (D-Dold).xs
//  		int row1 = cnode->Dofs()[j];
//  	  jump[j] = ((dmap[j])[row1]-(dmapold[j])[row1])*(cnode->xspatial()[j]);
//  	  
//  	  // loop over according master nodes
//  	  int k;
//  	  for (k = 0;k < mnodessize;++k)
//  	  {
//  	    // get node
//  	  	int gidm = mnodesvector[k];
//  	    DRT::Node* node = idiscret_->gNode(gidm);
//        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
//        CNode* mnode = static_cast<CNode*>(node);
//        
//        // -(M-Mols).xm
//        int col1 = mnode->Dofs()[j];
//        jump[j] -= ((mmap[j])[col1]-(mmapold[j])[col1])*mnode->xspatial()[j];
//	    } // loop over master nodes
//
//  	  // Scaling with D Matrix
//  	  jump[j] = -jump[j]/(dmap[j])[row1]; 
//  	} // loop over rowsize
//  	
//  	// store it to nodes 
//  	for (int j=0;j<rowsize1;++j)
//  	{
//  		//cnode->jump()[j] = jump[j];
//  	}
  	
  	// lagrange multiplier and jump in tangential direction (projection) 
    double ztan    = tangent[0]*z[0] + tangent[1]*z[1];
    double jumptan = tangent[0]*jump[0] + tangent[1]*jump[1];

    // FIXGIT: later, we always start with stick conditions.
    // In the meantime, we set the lagrange multipliers to two
    // in case of being null!
    
    if(abs(ztan) < 0.001)
    {
      ztan = +1.2;
 	  	cout << "Warning: lagrange multiplier in tangential direction had been set to" 
    	" 100 (hard coded)" << endl;
    }
    
    // epk = gp/(abs(ztan + ct*utan))
    double temp = ztan + ct*jumptan;
    double epk = frbound/abs(temp);
    
    // Fpk = ztan*(ztan + ct*utan)/(gp*(abs(ztan + ct*utan))
    double Fpk = ztan*temp/(frbound*abs(temp));
    
//    // Modification of Robin System - not working already 
//    // First Modification: Modification of Fpk 
//    if (abs(ztan)>frbound)
//    {
//    	// Fpk
//    	Fpk = ztan*temp/(abs(ztan)*abs(temp));
//      cout << "Modification of F" << endl;   
//    }  
//    
//    // alpha
//    double alpha = ztan*temp/(abs(ztan)*abs(temp));
//      
//    // delta
//    double delta = (abs(ztan))/frbound;
//    if (delta > 1) delta = 1;
//     
//    // beta
//    double beta;
//    if (alpha < 0)
//    {
//    	beta = 1/(1-(alpha*delta));
//      cout << "Modification of beta" << " " << "GID" << " " << gid << " " << beta <<  endl;   
//    }
//    else
//    {
//    	beta = 1;
//    }	
  	
    
    // Mpk = epk(1-Fpk)
    double Mpk = epk*(1-Fpk);
    
    //hpk = epk*ct*utan + ztan*epk*Fpk 
    double hpk = epk*ct*jumptan + ztan*epk*Fpk;
    
    Epetra_SerialDenseMatrix Lnode(1,1);
    Lnode(0,0)= Mpk/(1-Mpk)*ct;
    
    // First Modification: 
    //Lnode(0,0)= ct*(1/(1-beta*Mpk)-1);
        
    lmcol[0] = cnode->Dofs()[1];
        
    //assemble into L matrix
    lglobal.Assemble(-1,Lnode,lmrowT,lmrowownerT,lmcol); 
    
    Epetra_SerialDenseVector Rnode(1);
    vector<int> lm(1);
    vector<int> lmowner(1);

    Rnode(0) = -hpk/(1-Mpk);
    
    // First Modification:
    //Rnode(0) = -hpk/(1-beta*Mpk);
        
    lm[0] = cnode->Dofs()[1];
    lmowner[0] = cnode->Owner();

    LINALG::Assemble(rglobal,Rnode,lm,lmowner);
    }
      
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing normal+D+M derivatives       popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // nothing to do if no active nodes
  if (activenodes_==null)
    return;

  // not yet implemented for 3D
  if (Dim()==3)
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
    vector<map<int,double> > dnmap = cnode->GetDerivN();
    map<int,double>::iterator colcurr;
    int colsize = (int)dnmap[0].size();
    int mapsize = (int)dnmap.size();
    int row = activen_->GID(i);
    double* xi = cnode->xspatial();
    double* n = cnode->n();
    
    // only 2D so far
    if (mapsize==3) mapsize=2;
    
    for (int j=0;j<mapsize-1;++j)
      if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
        dserror("ERROR: AssembleS: Column dim. of nodal DerivN-map is inconsistent!");
         
    /*********************************************** N_normal_slave *****/
    //cout << endl << "->Assemble N_Normal_s for Node ID: " << cnode->Id() << endl;
    
    // we need the D-matrix entry of this node
    double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
    
    // loop over all derivative maps (=dimensions)
    for (int j=0;j<mapsize;++j)
    {
      int k=0;
    
      // loop over all entries of the current derivative map
      for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
      {
        int col = colcurr->first;
        double val = wii*xi[j]*colcurr->second;
        //cout << "wii=" << wii << " xi=" << xi[0] << " yi=" << xi[1] << " deriv=" << colcurr->second << endl;
        //cout << "Assemble N_Normal_s: " << row << " " << col << " " << val << endl;
        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
        ++k;
      }

      if (k!=colsize)
        dserror("ERROR: AssembleS: k = %i but colsize = %i",k,colsize);
    }
    
    /********************************************** N_normal_master *****/
    //cout << endl << "->Assemble N_Normal_m for Node ID: " << cnode->Id() << endl;
    
    // we need the M-matrix entries of this node
    vector<map<int,double> > mmap = cnode->GetM();
    map<int,double>::iterator mcurr;
    int mcolsize = (int)mmap[0].size();
    
    // get out of here, if no M-matric entries for this node
    if ((int)mmap.size()==0) continue;
    if (mcolsize%2!=0) dserror("ERROR: AssembleS: 3D case not yet implemented!");
    
    // loop over all master nodes (find adjacent ones to this active node)
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
      
      // compute S-matrix entry of the current active node / master node pair
      // loop over all derivative maps (=dimensions)
      for (int j=0;j<mapsize;++j)
      {
        int k=0;
      
        // loop over all entries of the current derivative map
        for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = mik*mxi[j]*colcurr->second;
          //cout << "mik=" << mik << " mxi=" << mxi[0] << " myi=" << mxi[1] << " deriv=" << colcurr->second << endl;
          //cout << "Assemble N_Normal_m: " << row << " " << col << " " << val << endl;
          // do not assemble zeros into s matrix
          if (abs(val)>1.0e-12) sglobal.Assemble(-val,row,col);
          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleS: k = %i but colsize = %i",k,colsize);
      }  
    }
    
    /********************************************** N_Dmortar_slave *****/
    //cout << endl << "->Assemble N_Dmortar_s for Node ID: " << cnode->Id() << endl;
    
    // we need the dot product n*x of this node
    double ndotx = 0.0;
    for (int dim=0;dim<cnode->NumDof();++dim)
      ndotx += n[dim]*xi[dim];
    
    // prepare assembly
    map<int,double>& ddmap = cnode->GetDerivD();
        
    // loop over all entries of the current derivative map
    for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = ndotx*colcurr->second;
      //cout << "ndotx=" << ndotx <<  " deriv=" << colcurr->second << endl;
      //cout << "Assemble N_DMortar_s: " << row << " " << col << " " << val << endl;
      // do not assemble zeros into s matrix
      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }
    
    /*************************************** N_Mmortar_slave+master *****/
    //cout << endl << "->Assemble N_Mmortar for Node ID: " << cnode->Id() << endl;
    
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
      double ndotx = 0.0;
      for (int dim=0;dim<cnode->NumDof();++dim)
        ndotx += n[dim]*mxi[dim];
          
      // compute S-matrix entry of the current active node / master node pair
      map<int,double>& thisdmmap = cnode->GetDerivM(gid);
      
      // loop over all entries of the current derivative map
      for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
      {
        int col = colcurr->first;
        double val = ndotx*colcurr->second;
        //cout << "ndotx=" << ndotx <<  " deriv=" << colcurr->second << endl;
        //cout << "Assemble N_Mmortar: " << row << " " << col << " " << val << endl;
        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) sglobal.Assemble(-val,row,col);
      }
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
  
  // not yet implemented for 3D
  if (Dim()==3)
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
    
    // prepare assembly
    vector<map<int,double> > dtmap = cnode->GetDerivTxi();
    map<int,double>::iterator colcurr;
    int colsize = (int)dtmap[0].size();
    int mapsize = (int)dtmap.size();
    int row = activet_->GID(i);
        
    // only 2D so far
    if (mapsize==3) mapsize=2;
    
    for (int j=0;j<mapsize-1;++j)
      if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
        dserror("ERROR: AssembleS: Column dim. of nodal DerivT-map is inconsistent!");
         
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
  
  // not yet implemented for 3D
  if (Dim()==3)
    return;
  
  /********************************************** LinDMatrix **********/
  // This is easy and can be done without communication, as the global
  // matrix lind has the same row map as the storage of the derivatives
  // of the D matrix.
  /**********************************************************************/
  
  // loop over all slave nodes (row map)
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    int numdof = cnode->NumDof();
    
    // Mortar matrix D derivatives
    map<int,double>& dderiv = cnode->GetDerivD();
    map<int,double>::iterator colcurr;
    
    // current Lagrange multipliers
    double* lm = cnode->lm();
    
    // loop over all slave node dofs for assembly
    for (int j=0;j<numdof;++j)
    {
      int row = cnode->Dofs()[j];
      
      // loop over all directional derivative entries
      for (colcurr=dderiv.begin();colcurr!=dderiv.end();++colcurr)
      {
        int col = colcurr->first;
        double val = lm[j] * (colcurr->second);
        //cout << "Assemble LinD: " << row << " " << col << " " << val << endl;
        
        // check entries for the one mortar loop case
        // (due to tolerances and round-off errors, we might get spurious entries)
        bool assemble = true;
#ifdef CONTACTONEMORTARLOOP
        if (sdoffullmap_->LID(col) < 0)
        {
          assemble=false;
          if(abs(val)>1.0e-6) dserror("ERROR: AssembleLinDM: Invalid non-zero master entry in LinD");
        }
#endif // #ifdef CONTACTONEMORTARLOOP
        
        if (assemble && abs(val)>1.0e-12) lindglobal.Assemble(val,row,col);
      }
    }
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
  vector<int> myslipTgids(slipnodes_->NumMyElements());
   
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
