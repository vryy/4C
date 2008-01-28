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

#include "pre_node.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::PreNode::PreNode(int id, const double* coords)
{
  id_ =id;
  for (int i=0; i<3; ++i) x_[i] = coords[i];
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::PreNode::~PreNode()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/08|
 *----------------------------------------------------------------------*/
EXODUS::PreNode::PreNode(const PreNode& old)
{
  id_ = old.id_;
  for (int i=0; i<3; ++i) x_[i] = old.x_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  print node    (public)                                     maf 01/08|
 *----------------------------------------------------------------------*/
void EXODUS::PreNode::Print(ostream& os, bool storeid) const
{
  if (storeid) os << "Id: " << StoreId();
  os << " ExoID: " << ExoId() << " Coords: "
     << X()[0] << "," << X()[1] << "," << X()[2] << endl;
}

#endif //D_EXODUS
