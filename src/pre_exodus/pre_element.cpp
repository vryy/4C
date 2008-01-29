/*----------------------------------------------------------------------*/
/*!
\file pre_node.H

\brief preprocessor reader for exodusII format 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with the preprocess node
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS

#include "pre_element.H"

/*!
\class PreElement

\brief PreElement is an element class for preprocess handling of elements. It is 
adopted to DRT::Elements

\author maf (frenzel@lnm.mw.tum.de)
*/


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::PreElement::PreElement(int id, ElementBlock::Shape distype)
{
  id_ =id;
  distype_ = distype;
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::PreElement::~PreElement()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::PreElement::PreElement(const PreElement& old)
{
  id_ = old.id_;
  distype_ = old.distype_;
  nodeid_ = old.nodeid_;
  node_ = old.node_;
  return;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                     maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::PreElement::Print(ostream& os) const
{
  os << "Id: " << Id() << " Shape: " ;
  os << ShapeToString(distype_) << endl;
  const int nnode = NumNode();
  const int* nodes = NodeIds();
  if (nnode)
  {
    os << " Nodes ";
    for (int i=0; i<nnode; ++i) os << setw(10) << nodes[i] << " ";
  }
  os << endl;
}

/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                       maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::PreElement::SetNodeIds(const int nnode, const int* nodes)
{
  nodeid_.resize(nnode);
  for (int i=0; i<nnode; ++i) nodeid_[i] = nodes[i];
  node_.resize(0);
  return;
}


#endif //D_EXODUS
