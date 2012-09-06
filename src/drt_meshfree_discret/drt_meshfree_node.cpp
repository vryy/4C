/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_node.cpp
 *
 * \brief A derived class of (meshfree) nodes with additional pointer to nodes
 *
 * A class of meshfree nodes derived from the generic DRT:Node and used in
 * meshfree discretisations. Additional to the element pointer it also has an
 * pointer to DRT::Nodes to assign connectivity. The meshfree node should be
 * used as geometric points defining integration cells and thus be dof-free.
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
 * \date January, 2012
 *
\*--------------------------------------------------------------------------*/

#include "drt_meshfree_node.H"
#include "../drt_lib/drt_dserror.H"

/*==========================================================================*\
 *                                                                          *
 * class MeshfreeNodeType                                                   *
 *                                                                          *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 | self-instantiation as parallel object type (?)                 nis Mar12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNodeType DRT::MESHFREE::MeshfreeNodeType::instance_;

/*--------------------------------------------------------------------------*
 | creates meshfree node                                 (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
DRT::ParObject* DRT::MESHFREE::MeshfreeNodeType::Create( const std::vector<char> & data )
{
  double dummycoord[3] = {999.,999.,999.};
  DRT::MESHFREE::MeshfreeNode* object = new DRT::MESHFREE::MeshfreeNode(-1,dummycoord,-1);
  object->Unpack(data);
  return object;
}


/*==========================================================================*\
 *                                                                          *
 * class MeshfreeNode                                                       *
 *                                                                          *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode::MeshfreeNode(int id, const double* coords, const int owner) :
DRT::Node(id, coords, owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode::MeshfreeNode(const DRT::MESHFREE::MeshfreeNode& old) :
Node(old)
{
  return;
}

/*--------------------------------------------------------------------------*
 | Deep copy this instance of MeshfreeNode and return pointer to it         |
 |                                                       (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode* DRT::MESHFREE::MeshfreeNode::Clone() const
{
  DRT::MESHFREE::MeshfreeNode* newnode = new DRT::MESHFREE::MeshfreeNode(*this);
  return newnode;
}

/*--------------------------------------------------------------------------*
 |  dtor                                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode::~MeshfreeNode()
{
  return;
}

/*--------------------------------------------------------------------------*
 | Delete a node from the meshfree node                  (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeNode::DeleteNodePtr(int gid)
{
  for (unsigned int i = 0; i<nodeid_.size(); i++){
    if (nodeid_[i]==gid){
      nodeid_.erase(nodeid_.begin()+i);
      node_.erase(node_.begin()+i);
      return;
    }
  }
  dserror("Connectivity issues: No node with specified gid to delete in meshfree node.");
  return;
}

/*----------------------------------------------------------------------*
 |  print this meshfree node                          (public) nis Mar12|
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeNode::Print(ostream& os) const
{
  // Print id and coordinates
  os << "MeshfreeNode " << setw(12) << Id()
     << " Owner " << setw(4) << Owner()
     << " Coords "
     << setw(12) << X()[0] << " "
     << setw(12) << X()[1] << " "
     << setw(12) << X()[2] << " ";
}
/*----------------------------------------------------------------------*
 |  << operator for meshfree node                     (public) nis Mar12|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::MESHFREE::MeshfreeNode& meshfreenode)
{
  meshfreenode.Print(os);
  return os;
}
