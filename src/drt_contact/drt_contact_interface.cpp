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
#include "../io/gmsh.H"
#include "drt_contact_projector.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Interface::Interface(const int id, const Epetra_Comm& comm) :
id_(id),
comm_(comm)
{
  RCP<Epetra_Comm> com = rcp(Comm().Clone());
  idiscret_ = rcp(new DRT::Discretization((string)"Contact Interface",com));
  contactsegs_.Reshape(0,0);
  counter_ = 0;
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
  // get standard element column map
  RCP<Epetra_Map> oldelecolmap = rcp (new Epetra_Map(*(Discret().ElementColMap() )));
      
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
  
  // do the same business for elements
  // (get row and column maps of slave and master elements seperately)
  {
  	const Epetra_Map* elerowmap = Discret().ElementRowMap();
  	const Epetra_Map* elecolmap = Discret().ElementColMap();
  	vector<int> sc;
  	vector<int> sr;
  	vector<int> scfull;
  	vector<int> mc;
  	vector<int> mr;
  	vector<int> mcfull;
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
  
  // do part of the same business for dofs
  // (get row map of slave and master dofs seperately)
  {
  	const Epetra_Map* noderowmap = Discret().NodeRowMap();
  	vector<int> sr;
  	vector<int> mr;
  	for (int i=0; i<noderowmap->NumMyElements();++i)
  	{
  		int gid = noderowmap->GID(i);
  		DRT::Node* node = Discret().gNode(gid);
  		if (!node) dserror("ERROR: Cannot find node with gid %",gid);
  		CNode* cnode = static_cast<CNode*>(node);
  	
  		if (cnode->IsSlave())
  			for (int j=0;j<cnode->NumDof();++j)
  				sr.push_back(cnode->Dofs()[j]);
  		else
  			for (int j=0;j<cnode->NumDof();++j)
  				mr.push_back(cnode->Dofs()[j]);
  	}
  	sdofrowmap_ = rcp(new Epetra_Map(-1,(int)sr.size(),&sr[0],0,Comm()));
  	mdofrowmap_ = rcp(new Epetra_Map(-1,(int)mr.size(),&mr[0],0,Comm()));
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  initiliaze / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Initialize()
{
	// loop over all nodes to reset normals, closestnode and Mortar maps
	// (use fully overlapping column map)
	for (int i=0;i<idiscret_->NumMyColNodes();++i)
	{
		CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));
		
		//reset nodal normal vector
		for (int j=0;j<3;++j)
			node->n()[j]=0.0;
		
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
		
		// reset nodal weighted gap
		node->Getg() = 0.0;
	}
	
	// loop over all elements to set current element length / area
	// and to reset contact candidates / search lists
	// (use fully overlapping column map)
	for (int i=0;i<idiscret_->NumMyColElements();++i)
	{
		CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
		element->Area()=element->ComputeArea();
		element->SearchElements().resize(0);
	}
	
	// reset matrix containing interface contact segments (gmsh)
	CSegs().Shape(0,0);
	
	return;
}

/*----------------------------------------------------------------------*
 |  set current deformation state                             popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::SetState(const string& statename, const RCP<Epetra_Vector> vec)
{
	if (statename=="displacement")
	{
		// set displacements in interface discretization
		idiscret_->SetState(statename, vec);
  
		// Get vec to full overlap
		RCP<const Epetra_Vector> global = idiscret_->GetState(statename);

		// alternative method to get vec to full overlap
		// Epetra_Vector global(*idiscret_->DofColMap(),false);
		// LINALG::Export(*vec,global);
    
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

  		DRT::UTILS::ExtractMyValues(*global,mydisp,lm);
  	
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
	
	// loop over proc's slave nodes of the interface
	// use standard column map to include processor's ghosted nodes 
	for(int i=0; i<snodecolmap_->NumMyElements();++i)
	{
		int gid = snodecolmap_->GID(i);
		DRT::Node* node = idiscret_->gNode(gid);
		if (!node) dserror("ERROR: Cannot find node with gid %",gid);
		CNode* cnode = static_cast<CNode*>(node);
/*		
#ifdef DEBUG
		cnode->Print(cout);
		cout << endl;
#endif // #ifdef DEBUG
*/		
	  // build averaged normal at each slave node
	  cnode->BuildAveragedNormal();
	}
	
#ifdef DEBUG
	comm_.Barrier();
#endif // #ifdef DEBUG
	
	// contact search algorithm
	EvaluateContactSearch();
	
	// loop over proc's slave elements of the interface for integration
	// use standard column map to include processor's ghosted elements 
	for (int i=0; i<selecolmap_->NumMyElements();++i)
	{
		int gid1 = selecolmap_->GID(i);
		DRT::Element* ele1 = idiscret_->gElement(gid1);
		if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
		CElement* selement = static_cast<CElement*>(ele1);
		
		// loop over the contact candidate master elements
		// use slave element's candidate list SearchElements !!!
		for (int j=0;j<selement->NumSearchElements();++j)
		{
			int gid2 = selement->SearchElements()[j];
			DRT::Element* ele2 = idiscret_->gElement(gid2);
			if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
			CElement* melement = static_cast<CElement*>(ele2);
			
			// check for element overlap and integrate the pair
			EvaluateOverlap_2D(*selement,*melement);
		}
	}

#ifdef DEBUG
	// Visualize node projections with gmsh
	// currently only works for 1 processor case because of writing of output
	if (comm_.NumProc()==1)
		VisualizeGmsh(CSegs());
#endif // #ifdef DEBUG
	
	// Exit 0 - Debug
	// exit(0);
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
	for (int i=0; i<snodecolmap_->NumMyElements();++i)
	{
		int gid = snodecolmap_->GID(i);
		DRT::Node* node = idiscret_->gNode(gid);
		if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
		CNode* snode = static_cast<CNode*>(node);
		
		// find closest master node to current slave node
		double mindist = 1.0e12;
		CNode* closestnode = snode->FindClosestNode(idiscret_,mnodefullmap_,mindist);
		snode->ClosestNode() = closestnode->Id();
		
		// proceed only if nodes are not far from each other!!!
		if (mindist<=CONTACT_CRITDIST)
		{
			// get adjacent elements to current slave node and to closest node
			int neles = snode->NumElement();
			DRT::Element** adj_slave = snode->Elements();
			int nelec = closestnode->NumElement();
			DRT::Element** adj_closest = closestnode->Elements();	
			
			// get global element ids for closest node's adjacent elements
			std::vector<int> cids(nelec);
			for (int j=0;j<nelec;++j)
				cids[j]=adj_closest[j]->Id();
			
			// try to add these to slave node's adjacent elements' search list
			for (int j=0;j<neles;++j)
			{
				CElement* selement = static_cast<CElement*> (adj_slave[j]);
				selement->AddSearchElements(cids);
			}
		}
	}
	
	// loop over all master nodes for closest node detection
	// use full overlap column map to include all nodes 
	for (int i=0; i<mnodefullmap_->NumMyElements();++i)
	{
		int gid = mnodefullmap_->GID(i);
		DRT::Node* node = idiscret_->gNode(gid);
		if (!node) dserror("ERROR: Cannot find master node with gid %",gid);
		CNode* mnode = static_cast<CNode*>(node);
			
		// find closest slave node to current master node
		double mindist = 1.0e12;
		CNode* closestnode = mnode->FindClosestNode(idiscret_,snodefullmap_,mindist);
		mnode->ClosestNode() = closestnode->Id();
		
		// proceed only if nodes are not far from each other!!!
		if (mindist<=CONTACT_CRITDIST)
		{
			// get adjacent elements to current master node and to closest node
			int nelem = mnode->NumElement();
			DRT::Element** adj_master = mnode->Elements();
			int nelec = closestnode->NumElement();
			DRT::Element** adj_closest = closestnode->Elements();	
			
			// get global element ids for master node's adjacent elements
			std::vector<int> mids(nelem);
			for (int j=0;j<nelem;++j)
				mids[j]=adj_master[j]->Id();
			
			// try to add these to closest node's adjacent elements' search list
			for (int j=0;j<nelec;++j)
			{
				CElement* selement = static_cast<CElement*> (adj_closest[j]);
				selement->AddSearchElements(mids);
			}				
		}
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Determine overlap and integrate sl/ma pair (public)       popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::EvaluateOverlap_2D(CONTACT::CElement& sele,
																			      CONTACT::CElement& mele)
{
	//cout << "Proc " << comm_.MyPID() << " checking pair... Slave ID: "
	//		 << sele.Id() << " Master ID: " << mele.Id() << endl;
	
	/************************************************************/
	/* There are several cases how the 2 elements can overlap.  */
	/* Handle all of them, including the ones that they don't   */
	/* overlap at all !																					*/
	/************************************************************/
	
	// create local booleans for projections of end nodes
	bool s0_hasproj = false;
	bool s1_hasproj = false;
	bool m0_hasproj = false;
	bool m1_hasproj = false;
	
	// get slave and master element nodes
	DRT::Node** mysnodes = sele.Nodes();
	if (!mysnodes)
		dserror("ERROR: EvaluateOverlap_2D: Null pointer for mysnodes!");
	DRT::Node** mymnodes = mele.Nodes();
	if (!mymnodes)
			dserror("ERROR: EvaluateOverlap_2D: Null pointer for mymnodes!");
	
	// create a 2-dimensional projector instance
	CONTACT::Projector projector(true);
		
	// project slave nodes onto master element
	vector<double> sprojxi(sele.NumNode());
	for (int i=0;i<2;++i)
	{
		CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(mysnodes[i]);
		double xi[2] = {0.0, 0.0};
		projector.Project_NodalNormal(*snode,mele,xi);
		sprojxi[i]=xi[0];
		
		// save projection if it is feasible
		// we need an expanded feasible domain in order to check pathological
		// cases due to round-off error and iteration tolerances later!
		if ((-1.0-CONTACT_PROJTOL<=sprojxi[i]) && (sprojxi[i]<=1.0+CONTACT_PROJTOL))
		{
			if (i==0) s0_hasproj=true;
			if (i==1) s1_hasproj=true;
		}		
	}
	
	// project master nodes onto slave element
	vector<double> mprojxi(mele.NumNode());
	for (int i=0;i<2;++i)
	{
		CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mymnodes[i]);
		double xi[2] = {0.0, 0.0};
		projector.Project_ElementNormal(*mnode,sele,xi);
		mprojxi[i]=xi[0];
		
		// save projection if it is feasible
		// we need an expanded feasible domain in order to check pathological
		// cases due to round-off error and iteration tolerances later!!!
		if ((-1.0-CONTACT_PROJTOL<=mprojxi[i]) && (mprojxi[i]<=1.0+CONTACT_PROJTOL))
		{
			if (i==0) m0_hasproj=true;
		  if (i==1) m1_hasproj=true;
		}		
	}
	
	/**********************************************************************/
	/* OVERLAP CASES																											*/
	/* depending on mxi and sxi overlap will be decided										*/
	/* even for 3noded CElements only the two end nodes matter in 2D!!!   */
	/**********************************************************************/
	//
	// For the non-overlapping cases, the possibility of an identical local
	// node numbering direction for both sides is taken into account!!
	// (this can happen, when elements far from each other are projected,
	// which actually should be impossible due to the CONTACT_CRITDIST
	// condition in the potential contact pair search above!
	// But you never know...)
	
	// For the overlapping cases, it is a prerequisite that the two local
	// node numbering directions are opposite!!
	// (this is the case, when the elements are sufficiently near each other,
	// which is ensured by only processing nodes that fulfill the 
	// CONTACT_CRITDIST condition above!)
	
	bool overlap = false;
	double sxia = 0.0;
	double sxib = 0.0;
	double mxia = 0.0;
	double mxib = 0.0;
	
	/* CASE 1 (NO OVERLAP):
	   no feasible projection found for any of the 4 outer element nodes  */
	
	if (!s0_hasproj && !s1_hasproj && !m0_hasproj && !m1_hasproj)
	{
		//do nothing
	}
		
	/* CASES 2-5 (NO OVERLAP):
	   feasible projection found only for 1 of the 4 outer element nodes
	   (this can happen due to the necessary projection tolerance!!!)     */
		
	else if  (s0_hasproj && !s1_hasproj && !m0_hasproj && !m1_hasproj)
	{
		if ((-1.0+CONTACT_PROJTOL<=sprojxi[0]) && (sprojxi[0]<=1.0-CONTACT_PROJTOL))
		{
			cout << "SElement Node IDs: " << mysnodes[0]->Id() << " " << mysnodes[1]->Id() << endl;
			cout << "MElement Node IDs: " << mymnodes[0]->Id() << " " << mymnodes[1]->Id() << endl;
			cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
			cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
			dserror("ERROR: EvaluateOverlap_2D: Significant overlap ignored S%i M%i!", sele.Id(), mele.Id());
		}
	}
	
	else if  (!s0_hasproj && s1_hasproj && !m0_hasproj && !m1_hasproj)
	{
		if ((-1.0+CONTACT_PROJTOL<=sprojxi[1]) && (sprojxi[1]<=1.0-CONTACT_PROJTOL))
		{
			cout << "SElement Node IDs: " << mysnodes[0]->Id() << " " << mysnodes[1]->Id() << endl;
			cout << "MElement Node IDs: " << mymnodes[0]->Id() << " " << mymnodes[1]->Id() << endl;
			cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
			cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
			dserror("ERROR: EvaluateOverlap_2D: Significant overlap ignored S%i M%i!", sele.Id(), mele.Id());
		}
	}
	
  else if  (!s0_hasproj && !s1_hasproj && m0_hasproj && !m1_hasproj)
	{
  	if ((-1.0+CONTACT_PROJTOL<=mprojxi[0]) && (mprojxi[0]<=1.0-CONTACT_PROJTOL))
  	{
  		cout << "SElement Node IDs: " << mysnodes[0]->Id() << " " << mysnodes[1]->Id() << endl;
  	  cout << "MElement Node IDs: " << mymnodes[0]->Id() << " " << mymnodes[1]->Id() << endl;
  		cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
  		cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
  		dserror("ERROR: EvaluateOverlap_2D: Significant overlap ignored S%i M%i!", sele.Id(), mele.Id());
  	}
	}
	
	else if  (!s0_hasproj && !s1_hasproj && !m0_hasproj && m1_hasproj)
	{
		if ((-1.0+CONTACT_PROJTOL<=mprojxi[1]) && (mprojxi[1]<=1.0-CONTACT_PROJTOL))
		{
			cout << "SElement Node IDs: " << mysnodes[0]->Id() << " " << mysnodes[1]->Id() << endl;
			cout << "MElement Node IDs: " << mymnodes[0]->Id() << " " << mymnodes[1]->Id() << endl;
			cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
			cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
			dserror("ERROR: EvaluateOverlap_2D: Significant overlap ignored S%i M%i!", sele.Id(), mele.Id());
		}
	}
	
	/* CASE 6 (OVERLAP):
		 feasible projection found for all 4 outer element nodes
		 (this can happen due to the necessary projection tolerance!!!)     */
	
	else if (s0_hasproj && s1_hasproj && m0_hasproj && m1_hasproj)
	{
		overlap = true;
		cout << "***WARNING***" << endl << "CONTACT::Interface::EvaluateOverlap_2D "<< endl
				 << "has detected '4 feasible projections'-case for Slave/Master pair "
				 << sele.Id() << "/" << mele.Id() << endl;
		
		// internal case 1 for global CASE 6
		// (equivalent to global CASE 7, slave fully projects onto master)
		if ((sprojxi[0]<1.0) && (sprojxi[1]>-1.0))
		{
			sxia = -1.0;
			sxib = 1.0;
			mxia = sprojxi[1];			// local node numbering always anti-clockwise!!!
			mxib = sprojxi[0];
		}
		
		// internal case 2 for global CASE 6
		// (equivalent to global CASE 8, master fully projects onto slave)
		else if ((mprojxi[0]<1.0) && (mprojxi[1]>-1.0))
		{
			mxia = -1.0;
			mxib = 1.0;
			sxia = mprojxi[1];			// local node numbering always anti-clockwise!!!
			sxib = mprojxi[0];
		}
		
		// internal case 3 for global CASE 6
	  // (equivalent to global CASE 9, both nodes no. 0 project successfully)
		else if ((sprojxi[0]<1.0) && (mprojxi[0]<1.0))
		{
			sxia = -1.0;
			sxib = mprojxi[0];			// local node numbering always anti-clockwise!!!
			mxia = -1.0;
			mxib = sprojxi[0];
		}
		
		// internal case 4 for global CASE 6
		// (equivalent to global CASE 10, both nodes no. 1 project successfully)
		else if ((sprojxi[1]>-1.0) && (mprojxi[1]>-1.0))
		{
			sxia = mprojxi[1];
			sxib = 1.0;						// local node numbering always anti-clockwise!!!
			mxia = sprojxi[1];
			mxib = 1.0;
		}
		
		// unknown internal case for global CASE 6
		else
			dserror("ERROR: EvaluateOverlap_2D: Unknown overlap case found in global case 6!");	
	}
	
	/* CASES 7-8 (OVERLAP):
		 feasible projections found for both nodes of one element, this 
	 	 means one of the two elements is projecting fully onto the other!  */
	
	else if (s0_hasproj && s1_hasproj && !m0_hasproj && !m1_hasproj)
	{
		overlap = true;
		sxia = -1.0;
		sxib = 1.0;
		mxia = sprojxi[1];			// local node numbering always anti-clockwise!!!
		mxib = sprojxi[0];
	}
	
	else if (!s0_hasproj && !s1_hasproj && m0_hasproj && m1_hasproj)
	{
		overlap = true;
		mxia = -1.0;
		mxib = 1.0;
		sxia = mprojxi[1];			// local node numbering always anti-clockwise!!!
		sxib = mprojxi[0];
	}
	
	/* CASES 9-10 (OVERLAP):
		 feasible projections found for one node of each element, due to
		 node numbering only identical local node ID pairs possible!        */
	
	else if (s0_hasproj && !s1_hasproj && m0_hasproj && !m1_hasproj)
	{
		// do the two elements really have an overlap?
		if ((sprojxi[0]>-1.0) && (mprojxi[0]>-1.0))
		{
			overlap = true;
			sxia = -1.0;
			sxib = mprojxi[0];			// local node numbering always anti-clockwise!!!
			mxia = -1.0;
			mxib = sprojxi[0];
		}
	}
	
	else if (!s0_hasproj && s1_hasproj && !m0_hasproj && m1_hasproj)
	{
		// do the two elements really have an overlap?
		if ((sprojxi[1]<1.0) && (mprojxi[1]<1.0))
		{
			overlap = true;
			sxia = mprojxi[1];
			sxib = 1.0;						// local node numbering always anti-clockwise!!!
			mxia = sprojxi[1];
			mxib = 1.0;
		}
	}
	
	/* CASES 11-14 (OVERLAP):
		 feasible projections found for 3 out of the total 4 nodes,
		 this can either lead to cases 7/8 or 9/10!                         */
	else if (s0_hasproj && s1_hasproj && m0_hasproj && !m1_hasproj)
	{
		overlap = true;
		// equivalent to global case 7
		if (mprojxi[0]>1.0)
		{
			sxia = -1.0;
			sxib = 1.0;
			mxia = sprojxi[1];		// local node numbering always anti-clockwise!!!
			mxib = sprojxi[0];
		}
		// equivalent to global case 9
		else
		{
			sxia = -1.0;
			sxib = mprojxi[0];		// local node numbering always anti-clockwise!!!
			mxia = -1.0;
			mxib = sprojxi[0];
		}
	}
	
	else if (s0_hasproj && s1_hasproj && !m0_hasproj && m1_hasproj)
	{
		overlap = true;
		// equivalent to global case 7
		if (mprojxi[1]<-1.0)
		{
			sxia = -1.0;
			sxib = 1.0;
			mxia = sprojxi[1];	// local node numbering always anti-clockwise!!!
			mxib = sprojxi[0];
		}
		// equivalent to global case 10
		else
		{
			sxia = mprojxi[1];
			sxib = 1.0;					// local node numbering always anti-clockwise!!!
			mxia = sprojxi[1];
			mxib = 1.0;
		}
	}
	
	else if (s0_hasproj && !s1_hasproj && m0_hasproj && m1_hasproj)
	{
		overlap = true;
		// equivalent to global case 8
		if (sprojxi[0]>1.0)
		{
			mxia = -1.0;
			mxib = 1.0;
			sxia = mprojxi[1];			// local node numbering always anti-clockwise!!!
			sxib = mprojxi[0];
		}
		// equivalent to global case 9
		else
		{
			sxia = -1.0;
			sxib = mprojxi[0];		// local node numbering always anti-clockwise!!!
			mxia = -1.0;
			mxib = sprojxi[0];
		}
	}
	
	else if (!s0_hasproj && s1_hasproj && m0_hasproj && m1_hasproj)
	{
		overlap = true;
		// equivalent to global case 8
		if (sprojxi[1]<-1.0)
		{
			mxia = -1.0;
			mxib = 1.0;
			sxia = mprojxi[1];	// local node numbering always anti-clockwise!!!
			sxib = mprojxi[0];
		}
		// equivalent to global case 10
		else
		{
			sxia = mprojxi[1];
			sxib = 1.0;					// local node numbering always anti-clockwise!!!
			mxia = sprojxi[1];
			mxib = 1.0;
		}
	}
		
	/* CASE DEFAULT: unknown overlap case													        */
	else
	{
		dserror("ERROR: EvaluateOverlap_2D: Unknown overlap case found!");
	}
	
	if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
		dserror("ERROR: EvaluateOverlap_2D: Determined infeasible limits!");

#ifdef DEBUG
	if (overlap)
	{
		//cout << "Found overlap!!!" << endl;
		//cout << "sxia: " << sxia << " sxib: " << sxib << endl;
		//cout << "mxia: " << mxia << " mxib: " << mxib << endl;
		
		// prepare gmsh visualization
		// currently only works for 1 processor case because of writing of output
		if (comm_.NumProc()==1)
		{
			double sxia_loc[2] = {sxia, 0.0};
			double sxib_loc[2] = {sxib, 0.0};
			double mxia_loc[2] = {mxia, 0.0};
			double mxib_loc[2] = {mxib, 0.0};
		
			double sxia_glob[3] = {0.0, 0.0, 0.0};
			double sxib_glob[3] = {0.0, 0.0, 0.0};
			double mxia_glob[3] = {0.0, 0.0, 0.0};
			double mxib_glob[3] = {0.0, 0.0, 0.0};
		
			sele.LocalToGlobal(sxia_loc,sxia_glob,true);
			sele.LocalToGlobal(sxib_loc,sxib_glob,true);
			mele.LocalToGlobal(mxia_loc,mxia_glob,true);
			mele.LocalToGlobal(mxib_loc,mxib_glob,true);
		
			Epetra_SerialDenseMatrix& segs = CSegs();
			segs.Reshape(segs.M()+1,12);
			segs(segs.M()-1,0)  = sxia_glob[0];
			segs(segs.M()-1,1)  = sxia_glob[1];
			segs(segs.M()-1,2)  = sxia_glob[2];
			segs(segs.M()-1,3)  = mxib_glob[0];
			segs(segs.M()-1,4)  = mxib_glob[1];
			segs(segs.M()-1,5)  = mxib_glob[2];
			segs(segs.M()-1,6)  = mxia_glob[0];
			segs(segs.M()-1,7)  = mxia_glob[1];
			segs(segs.M()-1,8)  = mxia_glob[2];
			segs(segs.M()-1,9)  = sxib_glob[0];
			segs(segs.M()-1,10) = sxib_glob[1];
			segs(segs.M()-1,11) = sxib_glob[2];
		}
	}
	else
	{
		//cout << "Did not find overlap!!!" << endl;
	}
#endif // #ifdef DEBUG
	
	/**********************************************************************/
	/* INTEGRATION																									    	*/
	/* depending on overlap, sxia, sxib, mxia and mxib 									*/
	/* integrate the Mortar matrices D, M and the weighted gap function   */
	/* on the overlap of the current Slave / Master CElement pair         */
	/**********************************************************************/
	
	// send non-overlapping pairs out of here
	if (!overlap)
		return true;
	
	// create a 2D integrator instance of GP order 5
	CONTACT::Integrator integrator(CONTACT_NGP,true);
	
	// do the three integrations
	RCP<Epetra_SerialDenseMatrix> D_seg = integrator.Integrate_D(sele,sxia,sxib);
	RCP<Epetra_SerialDenseMatrix> M_seg = integrator.Integrate_M(sele,sxia,sxib,mele,mxia,mxib);
	RCP<Epetra_SerialDenseVector> g_seg = integrator.Integrate_g(sele,sxia,sxib,mele,mxia,mxib);
	
	// do the three assemblies into the slave nodes
	integrator.Assemble_D(*this,sele,*D_seg);
	integrator.Assemble_M(*this,sele,mele,*M_seg);
	integrator.Assemble_g(*this,sele,*g_seg);
		
	// check for the modification of the M matrix for curved interfaces
	// (based on the paper by M. Puso / B. Wohlmuth, IJNME, 2005)
	bool modification = false;
	
	// conditions for modification
	// (1) linear shape functions for the slave elements
	// (2) 2D problems (3D unknown / unpublished ??)
	// (3) curved interface, n1!=n2
	// (4) use of dual shape functions for LM (always true in our case !!)
	if (sele.Shape()==DRT::Element::line2)
	{
		if (integrator.IsOneDimensional())
		{
			CNode* snode0 = static_cast<CNode*>(sele.Nodes()[0]);
			CNode* snode1 = static_cast<CNode*>(sele.Nodes()[1]);

			const double* n0 = snode0->n();
			const double* n1 = snode1->n();
			double delta = (n0[0]-n1[0])*(n0[0]-n1[0]) + (n0[1]-n1[1])*(n0[1]-n1[1]);
			
			if (delta>1.0e-8)
				modification = true;
		}
	}
	
	// integrate and assemble the modification, if necessary
	if (modification)
	{
		RCP<Epetra_SerialDenseMatrix> Mmod_seg = integrator.Integrate_Mmod(sele,sxia,sxib,mele,mxia,mxib);
		integrator.Assemble_Mmod(*this,sele,mele,*Mmod_seg);
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Visualize node projections with gmsh                      popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::VisualizeGmsh(const Epetra_SerialDenseMatrix& csegs)
{
	// increase counter variable by one
	counter_+=1;
	
	// construct unique filename for gmsh output
	std::ostringstream filename;
	filename << "o/gmsh_output/contact_";
	if (counter_<10)
		filename << 0 << 0 << 0 << 0;
	else if (counter_<100)
		filename << 0 << 0 << 0;
	else if (counter_<1000)
		filename << 0 << 0;
	else if (counter_<10000)
		filename << 0;
		
	filename << counter_ << ".pos";
	std::ofstream f_system(filename.str().c_str());
	
	// write output to temporary stringstream
	std::stringstream gmshfilecontent;
	gmshfilecontent << "View \" Slave and master side CElements \" {" << endl;
		
	// plot elements
	for (int i=0; i<idiscret_->NumMyColElements(); ++i)
	{
		CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
		int nnodes = element->NumNode();
		LINALG::SerialDenseMatrix coord(3,nnodes);
		double area = element->Area();
		
		// 2D linear case (2noded line elements)
		if (element->Shape()==DRT::Element::line2)
		{
			coord = element->GetNodalCoords();
			
			gmshfilecontent << "SL(" << scientific << coord(0,0) << "," << coord(1,0) << ","
							            << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
							            << coord(2,1) << ")";
			gmshfilecontent << "{" << scientific << area << "," << area << "};" << endl;
		}
		
		// 2D quadratic case (3noded line elements)
		if (element->Shape()==DRT::Element::line3)
		{
			coord = element->GetNodalCoords();
			
			gmshfilecontent << "SL2(" << scientific << coord(0,0) << "," << coord(1,0) << ","
							            << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
							            << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
							            << coord(2,2) << ")";
			gmshfilecontent << "{" << scientific << area << "," << area << "," << area << "};" << endl;
		}
	}
		
	// plot normal vectors
	for (int i=0; i<snodecolmap_->NumMyElements(); ++i)
	{
		int gid = snodecolmap_->GID(i);
		DRT::Node* node = idiscret_->gNode(gid);
		if (!node) dserror("ERROR: Cannot find node with gid %",gid);
		CNode* cnode = static_cast<CNode*>(node);
		if (!cnode) dserror("ERROR: Static Cast to CNode* failed");
		
	 	double nc[3];
	 	double nn[3];
	 	
	 	for (int j=0;j<3;++j)
	 	{
	 		nc[j]=cnode->xspatial()[j];
	 		nn[j]=100*cnode->n()[j];
	 	}
	 	
	 	gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
	 	gmshfilecontent << "{" << scientific << nn[0] << "," << nn[1] << "," << nn[2] << "};" << endl;
	}
		
	// plot contact segments (slave and master projections)
	if (csegs.M()!=0)
	{
		double sf = 100/csegs.M();
		for (int i=0; i<csegs.M(); ++i)
		{
			gmshfilecontent << "SQ(" << scientific << csegs(i,0) << "," << csegs(i,1) << ","
															 << csegs(i,2) << "," << csegs(i,3) << "," << csegs(i,4) << ","
															 << csegs(i,5) << "," << csegs(i,6) << "," << csegs(i,7) << ","
															 << csegs(i,8) << "," << csegs(i,9) << "," << csegs(i,10) << ","
															 << csegs(i,11) << ")";
			gmshfilecontent << "{" << scientific << i*sf << "," << i*sf << "," << i*sf << "," << i*sf << "};" << endl; 
		
			gmshfilecontent << "SL(" << scientific << csegs(i,0) << "," << csegs(i,1) << ","
											<< csegs(i,2) << "," << csegs(i,3) << "," << csegs(i,4) << ","
											<< csegs(i,5) << ")";
			gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl; 
		
			gmshfilecontent << "SL(" << scientific << csegs(i,6) << "," << csegs(i,7) << ","
											<< csegs(i,8) << "," << csegs(i,9) << "," << csegs(i,10) << ","
											<< csegs(i,11) << ")";
			gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
		}
	}
	
	gmshfilecontent << "};" << endl;
	
	// move everything to gmsh post-processing file and close it
	f_system << gmshfilecontent.str();
	f_system.close();
	
	return;
}

#endif  // #ifdef CCADISCRET
