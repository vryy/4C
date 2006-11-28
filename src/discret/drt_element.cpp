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
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::Element(int id, ElementType etype, int owner) :
ParObject(),
id_(id),
owner_(owner),
etype_(etype)
{
  nodeid_.resize(0);
  node_.resize(0);
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
node_(old.node_)
{

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Element::~Element()
{
  nodeid_.clear();
  node_.clear();
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
  os << Id() << " Owner " << Owner();
  const int nnode = NumNode();
  const int* nodes = NodeIds();
  if (nnode)
  {
    os << " Nodes ";
    for (int i=0; i<nnode; ++i) os << nodes[i] << " ";
  }
  // Print conditions if there are any
  int numcond = condition_.size();
  if (numcond)
  {
    os << endl; // end the previous line
    map<string,RefCountPtr<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << "Condition : " << curr->first << endl;
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
  const int sizeint    = sizeof(int);
  const int sizetype   = sizeof(enum ElementType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

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
const DRT::Condition* DRT::Element::GetCondition(const string& name) const
{
  map<string,RefCountPtr<Condition> >::const_iterator curr = condition_.find(name);
  if (curr != condition_.end()) return curr->second.get();
  else                          return NULL;
  return NULL;
}








#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
