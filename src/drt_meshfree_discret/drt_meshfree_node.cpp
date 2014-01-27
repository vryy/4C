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
  MeshfreeNode* object = new MeshfreeNode(-1,dummycoord,-1);
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
  DRT::Node(id, coords, owner),
  nodeid_(0),
  node_(0),
  facedim_(4)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeNode::MeshfreeNode(const DRT::MESHFREE::MeshfreeNode& old) :
  Node(old),
  nodeid_(old.nodeid_),
  node_(old.node_),
  facedim_(old.facedim_)
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

/*----------------------------------------------------------------------*
 |  print this meshfree node                          (public) nis Mar12|
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeNode::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Meshfree";
  DRT::Node::Print(os);
  return;
}
/*----------------------------------------------------------------------*
 |  << operator for meshfree node                     (public) nis Mar12|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::MeshfreeNode& meshfreenode)
{
  meshfreenode.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 | Delete a node from the meshfree node               (protected) nis Jan12 |
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

/*--------------------------------------------------------------------------*
 | Set vectors of ids and pointers to my nodes        (protected) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeNode::SetMyNodeTopology(
  const std::vector<int>& nodeid,
  const std::map<int,DRT::MESHFREE::MeshfreeNode*>& nodemap)
{
  const int nnode = nodeid.size();
  nodeid_.resize(nnode);
  node_.resize(nnode);
  std::vector<int>::const_iterator constit;
  std::vector<int>::iterator it = nodeid_.begin();
  std::vector<DRT::Node*>::iterator ptrit = node_.begin();
  for(constit=nodeid.begin(); constit!=nodeid.end(); ++constit)
  {
    *it = *constit;
    *ptrit = nodemap.at(*constit);
    ++it;
    ++ptrit;
  }
}

