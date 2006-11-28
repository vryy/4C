/*!----------------------------------------------------------------------
\file designnode.cpp
\brief A node that is part of a CAD design description

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_designnode.H"
#include "drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DesignNode::DesignNode(int id, const double* coords, const int owner) :
Node(id,coords,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DesignNode::DesignNode(const DRT::DesignNode& old) :
Node(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DesignNode::~DesignNode()
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public)   |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::DesignNode* DRT::DesignNode::Clone() const
{
  DRT::DesignNode* newnode = new DRT::DesignNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::DesignNode::Print(ostream& os) const
{
  Node::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::DesignNode::Pack(int& size) const
{
  const int sizeint    = sizeof(int);

  int basesize=0;
  const char* basedata = Node::Pack(basesize);
  
  size = 
  sizeint  + // holds size 
  sizeint  + // type of this instance of ParObject, see top of ParObject.H
  basesize + // holds basedata
  0;         // continue to add stuff here
  
  char* data = new char[size];
  int position=0;
  
  // add size
  AddtoPack(position,data,size);    
  // ParObject type
  int type = ParObject_DesignNode;
  AddtoPack(position,data,type);
  // add basedata
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;  

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::DesignNode::Unpack(const char* data)
{
  int position=0;
  
  // extract size
  int size=0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_DesignNode) dserror("Wrong instance type in data");
  
  // extract base class
  int basesize = SizePack(&data[position]);
  Node::Unpack(&data[position]);
  position += basesize;
  
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
