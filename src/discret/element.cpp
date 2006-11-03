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



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::Element(int id, ElementType etype) :
ParObject(),
id_(id),
etype_(etype)
{
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
node_(old.node_)
{

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element::~Element()
{
  node_.clear();
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CCADISCRETIZATION::Element& element)
{
  element.Print(); 
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
void CCADISCRETIZATION::Element::SetNodes(const int nnode, const int* nodes)
{
  node_.resize(nnode);
  for (int i=0; i<nnode; ++i) node_[i] = nodes[i];
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
  sizeint               +   // holds size itself
  sizeint               +   // holds Id()
  sizetype              +   // holds type of element
  sizeint               +   // holds size of vector node_
  node_.size()*sizeint  +   // holds vector node_
  0;                        // continue to add data here...

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
  // add vector node_
  AddVectortoPack(position,data,node_);
  // continue to add stuff here

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
  // extract node_
  ExtractVectorfromPack(position,data,node_);

  if (position != size)
  {
    cout << "CCADISCRETIZATION::Element::Unpack:\n"
         << "Mismatch in size of data " << size << " <-> " << position << endl
         << __FILE__ << ":" << __LINE__ << endl;
    exit(EXIT_FAILURE);
  }

  return true;
}












#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
