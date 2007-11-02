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
    vector<int> mc;
    vector<int> mr;
    for (int i=0; i<nodecolmap->NumMyElements(); ++i)
    {
      int gid      = nodecolmap->GID(i);
      bool isslave = dynamic_cast<CONTACT::CNode*>(Discret().gNode(gid))->IsSlave();
      if (isslave) sc.push_back(gid);
      else         mc.push_back(gid);
      if (!noderowmap->MyGID(gid)) continue;
      if (isslave) sr.push_back(gid);
      else         mr.push_back(gid);
    }
    snoderowmap_ = rcp(new Epetra_Map(-1,(int)sr.size(),&sr[0],0,Comm()));
    snodecolmap_ = rcp(new Epetra_Map(-1,(int)sc.size(),&sc[0],0,Comm()));
    mnoderowmap_ = rcp(new Epetra_Map(-1,(int)mr.size(),&mr[0],0,Comm()));
    mnodecolmap_ = rcp(new Epetra_Map(-1,(int)mc.size(),&mc[0],0,Comm()));
  }

  return;
}





#endif  // #ifdef CCADISCRET
