/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_node.H
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

#ifndef DRT_MESHFREE_NODE_H
#define DRT_MESHFREE_NODE_H

#include "../drt_lib/drt_node.H"

namespace DRT
{

namespace MESHFREE
{
/*==========================================================================*/
// forward declaration
/*==========================================================================*/
class MeshfreeDiscretization;

/*--------------------------------------------------------------------------*/
/*!
 * \brief The parallel object type class of meshfree nodes, creating the same
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeNodeType : public DRT::ParObjectType
{
public:

  //!< name of specific parallel object type
  std::string Name() const { return "MeshfreeNodeType"; }

  //!< returning instance of specific parallel object type
  static MeshfreeNodeType & Instance() { return instance_; };

  //!< create parallel object of this parallel object type
  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  //!< static instance of this parallel object type (for self-instantiation)
  static MeshfreeNodeType instance_;
};

/*--------------------------------------------------------------------------*/
/*!
 * \brief A derived class of (meshfree) nodes with additional pointer to nodes
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeNode : public DRT::Node
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of the control point
  */
  friend class Discretization;

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *
   * \param id     (in): A globally unique node id
   * \param coords (in): vector of nodal coordinates, length 3
   * \param owner  (in): Owner of this node.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeNode(int id, const double* coords, const int owner);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a MeshfreeNode
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeNode(const DRT::MESHFREE::MeshfreeNode& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy the derived class and return pointer to it
   *
   * This method is sort of a copy constructor for a class derived from
   * DRT::Node. It allows to copy construct the derived class without
   * knowing what it actually is using the base class Node.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual DRT::MESHFREE::MeshfreeNode* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeNode();

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of nodes adjacent to this node
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int NumNode() const { return nodeid_.size(); }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return id's of nodes adjacent to this meshfree node
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  const int* NodeIds() const
  {
    if (nodeid_.size()) return &nodeid_[0];
    else                return NULL;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return ptr to vector of nodes ptrs
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual DRT::Node** Nodes()
  { if (!nodeid_.empty()) return &node_[0]; else return NULL; }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return const ptr to vector of const nodes ptrs
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual const DRT::Node*const* Nodes() const
  { if (!nodeid_.empty()) return (const DRT::Node*const*)(&node_[0]); else return NULL; }


  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return eucledian distance to a given point
   *
   * \param x (in): point to which distance is computed
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual double DistToPoint (const double* x) const
  {
    return sqrt(pow(x_[0]-x[0],2.) + pow(x_[1]-x[1],2.) + pow(x_[2]-x[2],2.));
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add a node to the meshfree node
   *
   * \param nodeptr (in): pointer to DRT::Node to be added
   *
   * \warning It is your own responsibility to make sure that there will not
   *          be any double entries!
   *
   * Adds entry at the end of nodeid_ and node_ pointers
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AddNodePtr(DRT::Node*& nodeptr)
  {
    nodeid_.push_back(nodeptr->Id());
    node_.push_back(nodeptr);
    return;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add a node to the meshfree node
   *
   * \param gid     (in): global id of DRT::Node to be added
   * \param nodeptr (in): pointer to DRT::Node to be added
   *
   * \warning It is your own responsibility to make sure that there will not
   *          be any double entries! (This would be disastrous!)
   *
   * Adds entry at the end of nodeid_ and node_ pointers
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void AddNodePtr(const int& gid, DRT::Node*& nodeptr)
  {
    nodeid_.push_back(gid);
    node_.push_back(nodeptr);
    return;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete a node from the meshfree node
   *
   * \param gid (in): global id of node to be deleted
   *
   * Searches for position of node with specified gid and deletes entry in
   * vectors node_ and nodeid_
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void DeleteNodePtr(int gid);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this meshfree node
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Print(ostream& os) const ;

protected:

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Clear vector of pointers to my nodes
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void ClearMyNodeTopology()
  {
    nodeid_.clear();
    node_.clear();
    return;
  }

protected:
  // internal data:

  vector<int>        nodeid_;  //!< List of gids of my nodes, length NumNode()
  vector<DRT::Node*> node_;    //!< pointers to adjacent nodes

}; // end MeshfreeNode

} // end namespace MESHFREE
} // end namespace DRT

// << operator
ostream& operator << (ostream& os, const DRT::MESHFREE::MeshfreeNode& meshfreenode);

#endif  // #ifndef DRT_MESHFREE_NODE_H
