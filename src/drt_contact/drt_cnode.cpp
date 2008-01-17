/*!----------------------------------------------------------------------
\file drt_cnode.cpp
\brief A class for a contact node

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_cnode.H"
#include "../drt_lib/drt_dserror.H"
#include "drt_celement.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode::CNode(int id, const double* coords, const int owner, 
                      const int numdof, const vector<int>& dofs, const bool isslave) :
DRT::Node(id,coords,owner),
isslave_(isslave),
numdof_(numdof),
dofs_(dofs),
closestnode_(-1)
{
  for (int i=0;i<3;++i)
  {
  	n()[i]=0.0;
    u()[i]=0.0;
    xspatial()[i]=X()[i];
  }
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode::CNode(const CONTACT::CNode& old) :
DRT::Node(old),
isslave_(old.isslave_),
numdof_(old.numdof_),
dofs_(old.dofs_),
closestnode_(old.closestnode_)
{
	for (int i=0;i<3;++i)
	  {
	  	n()[i]=old.n_[i];
	    u()[i]=old.u_[i];
	    xspatial()[i]=old.xspatial_[i];
	  }
	return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of CNode and return pointer to it (public)  |
 |                                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CNode* CONTACT::CNode::Clone() const
{
  CONTACT::CNode* newnode = new CONTACT::CNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::CNode& cnode)
{
  cnode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::Print(ostream& os) const
{
  // Print id and coordinates
  os << "CNode " << setw(12) << Id()
     << " Owner " << setw(4) << Owner()
     << " Coords "
     << setw(12) << X()[0] << " "
     << setw(12) << X()[1] << " "
     << setw(12) << X()[2] << " "
     << " Dofs "; 
  for (int i=0; i<(int)dofs_.size(); ++i)
    os << dofs_[i] << " ";
  if (IsSlave()) os << " Slave Side  ";
  else           os << " Master Side ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Node
  vector<char> basedata(0);
  DRT::Node::Pack(basedata);
  AddtoPack(data,basedata);
  // add isslave_
  AddtoPack(data,isslave_);
  // add numdof_
  AddtoPack(data,numdof_);
  // add dofs_
  AddtoPack(data,dofs_);
  // add xspatial_
  AddtoPack(data,xspatial_,3);
  // add n_
  AddtoPack(data,n_,3);
  // add u_
  AddtoPack(data,u_,3);
  // add closestnode_
  AddtoPack(data,closestnode_);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::Node
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Node::Unpack(basedata);
  // isslave_
  ExtractfromPack(position,data,isslave_);
  // numdof_
  ExtractfromPack(position,data,numdof_);
  // dofs_
  ExtractfromPack(position,data,dofs_);
  // xspatial_
  ExtractfromPack(position,data,xspatial_,3);
  // n_
  ExtractfromPack(position,data,n_,3);
  // u_
  ExtractfromPack(position,data,u_,3);
  // closestnode_
  ExtractfromPack(position,data,closestnode_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  Build nodal normal                                        popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::CNode::BuildAveragedNormal()
{
	int nseg = NumElement();
	DRT::Element** adj_eles = Elements();
	
	for (int i=0;i<nseg;++i)
	{
		CElement* adj_cele = static_cast<CElement*> (adj_eles[i]);
/*		
#ifdef DEBUG
		adj_cele->Print(cout);
		cout << endl;
#endif // #ifdef DEBUG	 
*/	
		// build element normal at current node
		vector<double> ele_n(3);
		adj_cele->BuildNormalAtNode(Id(),ele_n);
		double wgt = adj_cele->Area();

/*
#ifdef DEBUG
		cout << "Area for CElement " << adj_cele->Id() << " is " << wgt << endl;
#endif // #ifdef DEBUG
*/		
		// add weighted element normal to nodal normal n_
		for (int j=0;j<3;++j)
		  n()[j]+=wgt*ele_n[j];
		
		/* average normal without weighting (see Diss. HARTMANN, 2007)
		for (int j=0;j<3;++j)
		 n()[j]+=ele_n[j];*/
	}
	
	// create unit normal vector
	double length = sqrt(n()[0]*n()[0]+n()[1]*n()[1]+n()[2]*n()[2]);
	
	if (length==0.0)
	  cout << "ERROR: Nodal normal of length zero" << endl;
	else
		for (int j=0;j<3;++j)
			n()[j]/=length;

/*
#ifdef DEBUG
	cout << endl;
	cout << "Unit normal for node " << Id() << " is " << n()[0] << " "
																										<< n()[1] << " "
																										<< n()[2] << endl;
	cout << endl;
#endif // #ifdef DEBUG
*/	
  return;
}

/*----------------------------------------------------------------------*
 |  Find closest node from given node set                     popp 01/08|
 *----------------------------------------------------------------------*/
CONTACT::CNode* CONTACT::CNode::FindClosestNode(const RCP<DRT::Discretization> intdis,
  																	            const RCP<Epetra_Map> nodesearchmap,
  																	            double& mindist)
{
	CNode* closestnode = NULL;
	
	// loop over all nodes of the DRT::Discretization that are
	// included in the given Epetra_Map
	for(int i=0; i<nodesearchmap->NumMyElements();++i)
	{
		int gid = nodesearchmap->GID(i);
		DRT::Node* node = intdis->gNode(gid);
		if (!node) dserror("ERROR: FindClosestNode: Cannot find node with gid %",gid);
		CNode* cnode = static_cast<CNode*>(node);
		
		// build distance between the two nodes
		double dist = 0.0;
		const double* p1 = this->xspatial();
		const double* p2 = cnode->xspatial();
		
		for (int j=0;j<3;++j)
			dist+=(p1[j]-p2[j])*(p1[j]-p2[j]);
		dist=sqrt(dist);
		
		// new closest node found, update
		if (dist <= mindist)
		{
			mindist=dist;
			closestnode=cnode;
		}
	}
	
	if (!closestnode)
		dserror("ERROR: FindClosestNode: No closest node found at all!");
	
	return closestnode;
}

#endif  // #ifdef CCADISCRET
