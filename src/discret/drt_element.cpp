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
  // we want a true deep copy of the condition_
  map<string,RefCountPtr<Condition> >::const_iterator fool;
  for (fool=old.condition_.begin(); fool!=old.condition_.end(); ++fool)
    condition_[fool->first] = rcp(new Condition(*(fool->second)));

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
  os << setw(12) << Id() << " Owner " << setw(5) << Owner();
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
    map<string,RefCountPtr<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << "Conditions : " << curr->first << " ";
      os << *(curr->second);
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
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Element::Pack(int& size) const
{
  const int sizeint            = sizeof(int);
  const int sizetype           = sizeof(enum ElementType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // get data and size of dofset_
  int dofsetsize=0;
  const char* dofsetpack = dofset_.Pack(dofsetsize);
  
  
  // get the size of all conditions
  int condsize=0;
  map<string,RefCountPtr<Condition> >::const_iterator curr;
  for (curr=condition_.begin(); curr != condition_.end(); ++curr)
  {
    condsize += SizeString(curr->first);
    int tmp=0;
    const char* condpack = curr->second->Pack(tmp);
    condsize += tmp;
    delete [] condpack;
  }

  // calculate size of vector
  size = 
  sizeint                +   // holds size itself
  sizeint                +   // type of this instance of ParObject, see top of ParObject.H
  sizeint                +   // holds Id()
  sizeint                +   // holds Owner()
  sizetype               +   // holds type of element
  SizeVector(nodeid_)    +   // nodeid_
  dofsetsize             +   // dofset_
  sizeint                +   // no. of objects in condition_
  condsize               +   // condition_
  0;                         // continue to add data here...

  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;

  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Element;
  AddtoPack(position,data,type);
  // add Id()
  int id = Id();
  AddtoPack(position,data,id);
  // add Id()
  int owner = Owner();
  AddtoPack(position,data,owner);
  // add type of element
  ElementType etype = Type();
  AddtoPack(position,data,etype);
  // add vector nodeid_
  AddVectortoPack(position,data,nodeid_);
  // dofset_
  AddtoPack(position,data,dofsetpack,dofsetsize);
  delete [] dofsetpack;
  // condition_
  int num = condition_.size(); // no. of objects
  AddtoPack(position,data,num);
  for (curr=condition_.begin(); curr != condition_.end(); ++curr)
  {
    AddStringtoPack(position,data,curr->first);
    int tmp=0;
    const char* condpack = curr->second->Pack(tmp);
    AddtoPack(position,data,condpack,tmp);
    delete [] condpack;
  }
  // continue to add stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Element::Unpack(const char* data)
{
  int position = 0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_Element) dserror("Wrong instance type in data");
  // extract id
  ExtractfromPack(position,data,id_);
  // extract owner_
  ExtractfromPack(position,data,owner_);
  // extract type
  ExtractfromPack(position,data,etype_);
  // extract nodeid_
  ExtractVectorfromPack(position,data,nodeid_);
  // dofset_
  dofset_.Unpack(&data[position]);
  int dofsetsize = dofset_.SizePack(&data[position]);
  position += dofsetsize;
  // extract condition_
  int num=0;
  ExtractfromPack(position,data,num);
  for (int i=0; i<num; ++i)
  {
    string name;
    ExtractStringfromPack(position,data,name);
    RefCountPtr<Condition> cond = rcp(new Condition());
    cond->Unpack(&data[position]);
    int size = cond->SizePack(&data[position]);
    position += size;
    SetCondition(name,cond);
  }
  // node_ is NOT communicated
  node_.resize(0);

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
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
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Condition* DRT::Element::GetCondition(const string& name)
{
  map<string,RefCountPtr<Condition> >::const_iterator curr = condition_.find(name);
  if (curr != condition_.end()) return curr->second.get();
  else                          return NULL;
  return NULL;
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
        dsassert(dirich->Type()==DRT::Condition::condition_Dirichlet,"condition with name Dirichlet is not of type condition_Dirichlet");
        flag = dirich->GetVector<int>("onoff");
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
    dsassert(dirich->Type()==DRT::Condition::condition_Dirichlet,"condition with name Dirichlet is not of type condition_Dirichlet");
    flag = dirich->GetVector<int>("onoff");
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






#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
