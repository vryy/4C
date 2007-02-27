/*!----------------------------------------------------------------------
\file element.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_element.H"
#include "drt_node.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(int id, ElementType etype, int owner) :
ParObject(),
id_(id),
owner_(owner),
etype_(etype),
dofset_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(const DRT::Element& old) :
ParObject(old),
id_(old.id_),
owner_(old.owner_),
etype_(old.etype_),
nodeid_(old.nodeid_),
node_(old.node_),
dofset_(old.dofset_)
{
  // we do NOT want a deep copy of the condition_ as the condition
  // is only a reference in the elements anyway
  map<string,RefCountPtr<Condition> >::const_iterator fool;
  for (fool=old.condition_.begin(); fool!=old.condition_.end(); ++fool)
    SetCondition(fool->first,fool->second);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::~Element()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Element& element)
{
  element.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Element::Print(ostream& os) const
{
  os << setw(12) << Id() << " Owner " << setw(5) << Owner() << " ";
  const int nnode = NumNode();
  const int* nodes = NodeIds();
  if (nnode)
  {
    os << " Nodes ";
    for (int i=0; i<nnode; ++i) os << setw(10) << nodes[i] << " ";
  }
  // Print dofs if there are any
  if (Dof().NumDof())
    cout << Dof() << " ";
  
  // Print conditions if there are any
  int numcond = condition_.size();
  if (numcond)
  {
    os << endl << numcond << " Conditions:\n";
    map<string,RefCountPtr<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << curr->first << " ";
      os << *(curr->second) << endl;
    }
  }
  
  
  return;
}

/*----------------------------------------------------------------------*
 |  read element input dummy (public)                        mwgee 11/06|
 |  this is a base class dummy for elements that do not need            |
 |  a reading method
 *----------------------------------------------------------------------*/
bool DRT::Element::ReadElement()
{
  cout << "DRT::Element::ReadElement:\n"
       << "Base class dummy routine Element::ReadElement() called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return false;
}

/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Element::SetNodeIds(const int nnode, const int* nodes)
{
  nodeid_.resize(nnode);
  for (int i=0; i<nnode; ++i) nodeid_[i] = nodes[i];
  node_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::Pack(vector<char>& data) const
{
  data.resize(0);
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add id
  int id = Id();
  AddtoPack(data,id);
  // add owner
  int owner = Owner();
  AddtoPack(data,owner);
  // add type of element
  ElementType etype = Type();
  AddtoPack(data,etype);
  // add vector nodeid_
  AddtoPack(data,nodeid_);
  // dofset
  vector<char> dofsetpack(0);
  dofset_.Pack(dofsetpack);
  AddtoPack(data,dofsetpack);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // id_
  ExtractfromPack(position,data,id_);
  // owner_
  ExtractfromPack(position,data,owner_);
  // etype_
  ExtractfromPack(position,data,etype_);
  // nodeid_
  ExtractfromPack(position,data,nodeid_);
  // dofset_
  vector<char> dofpack(0);
  ExtractfromPack(position,data,dofpack);
  dofset_.Unpack(dofpack);
  
  // node_ is NOT communicated
  node_.resize(0);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
} 
 

/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Element::BuildNodalPointers(map<int,RefCountPtr<DRT::Node> >& nodes)
{
  int        nnode   = NumNode();
  const int* nodeids = NodeIds();
  node_.resize(nnode);
  for (int i=0; i<nnode; ++i)
  {
    map<int,RefCountPtr<DRT::Node> >::iterator curr = nodes.find(nodeids[i]);
    // this node is not on this proc
    if (curr==nodes.end()) dserror("Element %d cannot find node %d",Id(),nodeids[i]);
    else 
      node_[i] = curr->second.get();
  } 
  return true;
}

/*----------------------------------------------------------------------*
 |  Build nodal pointers                                    (protected) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
bool DRT::Element::BuildNodalPointers(DRT::Node** nodes)
{
  node_.resize(NumNode()); 
  for (int i=0; i<NumNode(); ++i) node_[i] = nodes[i];
  return true;
}


/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::GetCondition(const string& name,vector<DRT::Condition*>& out)
{
  const int num = condition_.count(name);
  out.resize(num);
  multimap<string,RefCountPtr<Condition> >::iterator startit = 
                                         condition_.lower_bound(name);
  multimap<string,RefCountPtr<Condition> >::iterator endit = 
                                         condition_.upper_bound(name);
  int count=0;
  multimap<string,RefCountPtr<Condition> >::iterator curr;
  for (curr=startit; curr!=endit; ++curr)
    out[count++] = curr->second.get();
  if (count != num) dserror("Mismatch in number of conditions found");
  return;
}

/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Element::GetCondition(const string& name)
{
  multimap<string,RefCountPtr<Condition> >::iterator curr = 
                                         condition_.find(name);
  if (curr==condition_.end()) return NULL;
  curr = condition_.lower_bound(name);
  return curr->second.get();
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(vector<int>& lm, vector<int>& lmdirich, 
                                  vector<int>& lmowner)
{
  const int numnode = NumNode();
  // count how many degrees of freedom I have
  int count=0;
  // count nodal dofs
  DRT::Node** nodes = Nodes();
  if (nodes)
    for (int i=0; i<numnode; ++i)
      count += nodes[i]->Dof().NumDof();
  // add element dofs
  count += Dof().NumDof();

  lm.resize(count);
  lmdirich.resize(count);
  lmowner.resize(count);
  for (int i=0; i<count; ++i) lmdirich[i] = 0;
  
  // fill the vector with nodal dofs
  int count2=0;
  if (nodes)
    for (int i=0; i<numnode; ++i)
    {
      DRT::Condition* dirich = nodes[i]->GetCondition("Dirichlet");
      vector<int>* flag = NULL;
      if (dirich)
      {
        if (dirich->Type()!=DRT::Condition::PointDirichlet &&
            dirich->Type()!=DRT::Condition::LineDirichlet &&
            dirich->Type()!=DRT::Condition::SurfaceDirichlet &&
            dirich->Type()!=DRT::Condition::VolumeDirichlet)
          dserror("condition with name dirichlet is not of type Dirichlet");
        flag = dirich->Get<vector<int> >("onoff");
      }
      const int owner = nodes[i]->Owner();
      for (int j=0; j<nodes[i]->Dof().NumDof(); ++j)
      {
        if (flag)
          if ((*flag)[j]) 
            lmdirich[count2] = 1;
        lmowner[count2] = owner;
        lm[count2++]    = nodes[i]->Dof()[j];
      }
    }

  // fill the vector with element dofs
  vector<int>* flag = NULL;
  DRT::Condition* dirich = GetCondition("Dirichlet");
  if (dirich)
  {
    if (dirich->Type()!=DRT::Condition::PointDirichlet &&
        dirich->Type()!=DRT::Condition::LineDirichlet &&
        dirich->Type()!=DRT::Condition::SurfaceDirichlet &&
        dirich->Type()!=DRT::Condition::VolumeDirichlet)
      dserror("condition with name dirichlet is not of type Dirichlet");
    flag = dirich->Get<vector<int> >("onoff");
  }
  const int owner = Owner();
  for (int j=0; j<Dof().NumDof(); ++j)
  {
    if (flag)
      if ((*flag)[j]) 
        lmdirich[count2] = 1;
    lmowner[count2] = owner;
    lm[count2++]    = Dof()[j];
  }
    
  if (count2!=count) dserror("Mismatch in no. of dofs");
  
  return;
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Element::LocationVector(vector<int>& lm, vector<int>& lmowner)
{
  const int numnode = NumNode();
  // count how many degrees of freedom I have
  int count=0;
  // count nodal dofs
  DRT::Node** nodes = Nodes();
  if (nodes)
    for (int i=0; i<numnode; ++i)
      count += nodes[i]->Dof().NumDof();
  // add element dofs
  count += Dof().NumDof();

  lm.resize(count);
  lmowner.resize(count);
  
  // fill the vector with nodal dofs
  int count2=0;
  if (nodes)
    for (int i=0; i<numnode; ++i)
    {
      const int owner = nodes[i]->Owner();
      for (int j=0; j<nodes[i]->Dof().NumDof(); ++j)
      {
        lmowner[count2] = owner;
        lm[count2++]    = nodes[i]->Dof()[j];
      }
    }

  // fill the vector with element dofs
  const int owner = Owner();
  for (int j=0; j<Dof().NumDof(); ++j)
  {
    lmowner[count2] = owner;
    lm[count2++]    = Dof()[j];
  }
    
  if (count2!=count) dserror("Mismatch in no. of dofs");
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Element::Evaluate(ParameterList& params,
                           DRT::Discretization&      discretization,
                           vector<int>&              lm,
                           Epetra_SerialDenseMatrix& elemat1,
                           Epetra_SerialDenseMatrix& elemat2,
                           Epetra_SerialDenseVector& elevec1,
                           Epetra_SerialDenseVector& elevec2,
                           Epetra_SerialDenseVector& elevec3)
{
  cout << "DRT::Element::Evaluate:\n"
       << "Base class dummy routine DRT::Element::Evaluate(...) called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return -1;
}

/*----------------------------------------------------------------------*
 |  evaluate Neumann BC dummy (public)                       mwgee 01/07|
 *----------------------------------------------------------------------*/
int DRT::Element::EvaluateNeumann(ParameterList& params, 
                                  DRT::Discretization&      discretization,
                                  DRT::Condition&           condition,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseVector& elevec1)
{
  cout << "DRT::Element::EvaluateNeumann:\n"
       << "Base class dummy routine DRT::Element::EvaluateNeumann(...) called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return -1;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
