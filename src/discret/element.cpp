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

#include "element.H"
#include "dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::Element(int id, ElementType etype) :
ParObject(),
id_(id),
etype_(etype)
{
  nodeid_.resize(0);
  node_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::Element(const CCADISCRETIZATION::Element& old) :
ParObject(old),
id_(old.id_),
etype_(old.etype_),
nodeid_(old.nodeid_),
node_(old.node_)
{

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::~Element()
{
  nodeid_.clear();
  node_.clear();
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CCADISCRETIZATION::Element& element)
{
  element.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  read element input dummy (public)                        mwgee 11/06|
 |  this is a base class dummy for elements that do not need            |
 |  a reading method
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Element::ReadElement()
{
  cout << "CCADISCRETIZATION::Element::ReadElement:\n"
       << "Base class dummy routine Element::ReadElement() called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return false;
}

/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Element::SetNodeIds(const int nnode, const int* nodes)
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
const char* CCADISCRETIZATION::Element::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizetype   = sizeof(enum ElementType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // calculate size of vector
  size = 
  sizeint                +   // holds size itself
  sizeint                +   // holds Id()
  sizetype               +   // holds type of element
  sizeint                +   // holds size of vector nodeid_
  nodeid_.size()*sizeint +   // holds vector nodeid_
  0;                         // continue to add data here...

  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;

  // add size
  AddtoPack(position,data,size);
  // add Id()
  int id = Id();
  AddtoPack(position,data,id);
  // add type of element
  ElementType etype = Type();
  AddtoPack(position,data,etype);
  // add vector nodeid_
  AddVectortoPack(position,data,nodeid_);
  // continue to add stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Element::Unpack(const char* data)
{
  int position = 0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // extract id
  ExtractfromPack(position,data,id_);
  // extract type
  ExtractfromPack(position,data,etype_);
  // extract nodeid_
  ExtractVectorfromPack(position,data,nodeid_);
  // node_ is NOT communicated
  node_.resize(0);

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  Build nodal pointers                                       (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Element::BuildNodalPointers(map<int,RefCountPtr<CCADISCRETIZATION::Node> >& nodes)
{
  int        nnode   = NumNode();
  const int* nodeids = NodeIds();
  node_.resize(nnode);
  for (int i=0; i<nnode; ++i)
  {
    map<int,RefCountPtr<CCADISCRETIZATION::Node> >::iterator curr = nodes.find(nodeids[i]);
    if (curr==nodes.end()) return false;
    node_[i] = curr->second.get();
  } 
  return true;
}










#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
