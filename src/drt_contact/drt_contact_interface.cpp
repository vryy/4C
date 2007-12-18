/*!----------------------------------------------------------------------
\file drt_contact_interface.cpp
\brief One contact interface

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_interface.H"
#include "drt_cdofset.H"
#include "../drt_lib/linalg_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Interface::Interface(const int id, const Epetra_Comm& comm) :
id_(id),
comm_(comm)
{
  RCP<Epetra_Comm> com = rcp(Comm().Clone());
  idiscret_ = rcp(new DRT::Discretization((string)"Contact Interface",com));
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
  
  // get standard nodal column map
  RCP<Epetra_Map> oldnodecolmap = rcp (new Epetra_Map(*(Discret().NodeColMap() )));
      
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
    vector<int> sc;
    vector<int> sr;
    vector<int> scfull;
    vector<int> mc;
    vector<int> mr;
    vector<int> mcfull;
    for (int i=0; i<nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      bool isslave = dynamic_cast<CONTACT::CNode*>(Discret().gNode(gid))->IsSlave();
      if (oldnodecolmap->MyGID(gid))
      {
    	  if (isslave) sc.push_back(gid);
    	  else         mc.push_back(gid);
      }  
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
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  set current deformation state                             popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::SetState(const string& statename, const RCP<Epetra_Vector> vec)
{
  idiscret_->SetState(statename, vec);
  
  //Epetra_Vector global(*idiscret_->DofColMap(),false);
  //LINALG::Export(*vec,global);
  
  // Alternative method to get vec to full overlap
  RCP<const Epetra_Vector> global = idiscret_->GetState(statename);

/*#ifdef DEBUG
  vec->Print(cout);
  comm_.Barrier();
  global->Print(cout);
#endif // #ifdef DEBUG*/
  
  // loop over all nodes to set current displacement
  // usefully overlapping column map to make disp. available on all procs
  if (statename=="displacement")
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
  	CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));
  	const int numdof = node->NumDof();
  	vector<double> mydisp(numdof);
  	vector<int> lm(numdof);
  	
  	for (int j=0;j<numdof;++j)
  		lm[j]=node->Dofs()[j];

  	DRT::Utils::ExtractMyValues(*global,mydisp,lm);
  	
  	// add mydisp[2]=0 for 2D problems
  	if (mydisp.size()<3)
  	  		mydisp.resize(3);
  	
/*#ifdef DEBUG
  	cout << "PROC " << comm_.MyPID() << " thinks that CNode " << node->Id() << ", owned by PROC " << node->Owner()
  			 << ", Dofs " << lm[0] << " " << lm[1] << " " << lm[2]
  			 << ", has the displacements): \n"
  			 << mydisp[0] << "\t\t"
  			 << mydisp[1] << "\t\t"
  			 << mydisp[2] << "\n\n";
#endif // #ifdef DEBUG*/
  	
  	// set current configuration and displacement, reset normal
  	for (int j=0;j<3;++j)
    {
  		node->u()[j]=mydisp[j];
  		node->xspatial()[j]=node->X()[j]+mydisp[j];
  		node->n()[j]=0.0;
    }
  }
  
  // loop over all elements to set current area or length
  // use fully overlapping column map to make areas available on all procs
  if (statename=="displacement")
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
  	CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
  	element->Area()=element->ComputeArea();
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                 popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Evaluate()
{
	// interface needs to be complete
	if (!Filled() && comm_.MyPID()==0)
	    dserror("ERROR: FillComplete() not called on interface %", id_);
	
	// loop over all slave nodes of the interface
	// use standard column map to include processor "border nodes" 
	for(int i=0; i<snodecolmap_->NumMyElements();++i)
	{
		int gid = snodecolmap_->GID(i);
		DRT::Node* node = idiscret_->gNode(gid);
		if (!node) dserror("ERROR: Cannot find node with gid %",gid);
		CNode* cnode = static_cast<CNode*>(node);
		
#ifdef DEBUG
		cnode->Print(cout);
		cout << endl;
#endif // #ifdef DEBUG
		
	  // build averaged normal at each slave node
	  cnode->BuildAveragedNormal();
	  
	  // find projection of each slave node to master surface
	  // NOT YET IMPLEMENTED!!!
	  // Project_SlaveToMaster();
	}
	
	// Exit 0 - Debug
	//exit(0);
  return;
}


#endif  // #ifdef CCADISCRET
