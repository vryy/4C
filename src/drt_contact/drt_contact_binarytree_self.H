/*!----------------------------------------------------------------------
\file drt_contact_binarytree_self.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_BINARYTREE_SELF_H
#define DRT_CONTACT_BINARYTREE_SELF_H

#include "../drt_lib/drt_discret.H"
#include "drt_celement.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

//! @name Enums and Friends
enum BinaryTreeSelfNodeType
{
  SELF_INNER,        ///< indicates an inner node (has children)
  SELF_LEAF,         ///< indicates a leaf node (no further children)
  SELF_NO_ELEMENTS,  ///< inidicates that there are no elements on this (root) treenode
  SELF_UNDEFINED     ///< indiactes an undefined tree node
};

//@}

/*!
\brief A class representing one tree node of the binary tree for self contact

Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

\author popp (popp@lnm.mw.tum.de)
*/
class BinaryTreeSelfNode
{
public:
  
  /*!
  \brief Constructor of a tree node
  
  \param type           type of BinaryTreeSelfNode
  \param discret        contact interface discretization 
  \param parent         points to parent tree node
  \param elelist        list of all elements in BinaryTreeSelfNode
  \param dopnormals     reference to DOP normals
  \param samplevectors  reference to sample vectors
  \param kdop           refernece to no. of vertices
  \param dim            dimension of problem    
  \param layer          current layer of treenode
  \param treenodes      references to treenodes storage scheme 
  
  */
  BinaryTreeSelfNode(BinaryTreeSelfNodeType type, DRT::Discretization& discret, 
                     RCP<BinaryTreeSelfNode> parent, vector<int> elelist,
                     const Epetra_SerialDenseMatrix& dopnormals,
                     const Epetra_SerialDenseMatrix& samplevectors,
                     const int& kdop,
                     const int& dim,const int& nvectors, const int layer,
                     vector<vector<RCP<BinaryTreeSelfNode> > > & treenodes);
  
  /*!
  \brief Destructor
  
  */
  virtual ~BinaryTreeSelfNode() {};
  
  //! @name Evaluation methods
  
  /*!
  \brief Calculate slabs of dop in initialization
  
  */
  void CalculateSlabsDop(bool isinit);
  
  /*!
  \brief Update slabs of current treenode in bottom up way
  
  */
  void UpdateSlabsBottomUp(double & eps);
  
  /*!
  \brief Calculate the logical array of qualified sample vectors for leaf nodes
  
  */
  void CalculateQualifiedVectors();
  
  /*!
  \brief Update the logical array of qualified sample vectors for non-leaf nodes
  
  */
  void UpdateQualifiedVectorsBottomUp();
  
  /*!
  \brief Return logical array of qualified sample vectors
  
  */
  virtual inline vector<bool> QualifiedVectors() { return qualifiedvectors_; }
  
  /*!
  \brief Set logical array of qualified sample vectors
  
  */
  virtual inline void SetQualifiedVectors(vector<bool> qualifiedvectors)
  { qualifiedvectors_ = qualifiedvectors; }
  
  /*!
  \brief Enlarge geometry of a Treenode by an offset, dependent on size
  
  */    
  void EnlargeGeometry(double& eps);
  
  /*!
  \brief Print type of TreeNode to std::cout
  
  */
  void PrintType();
  
  /*!
  \brief Print slabs to std::cout
  
  */
  void PrintSlabs();
  
  /*!
  \brief Set slabs of current treenode
  
  */
  void SetSlabs(Epetra_SerialDenseMatrix& newslabs);
  
  //@}  
  
  //! @name Access methods
  
  /*!
  \brief Get communicator
  
  */
  virtual const Epetra_Comm& Comm() const { return idiscret_.Comm(); }
  
  
  /*!
  \brief Complete tree by filling treenode storage scheme
  
  */
  void CompleteTree(int layer, double& enlarge);
  
  /*!
    
  /brief Get discretization of the interface
  
  */
  virtual inline DRT::Discretization& Discret() const { return idiscret_; }
  
  /*!
  \brief Return pointer to type of treenode
  
  */
  virtual inline BinaryTreeSelfNodeType Type() const { return type_; }
  
  /*!
  \brief Set type of treenode
  
  */
  virtual inline void SetType(BinaryTreeSelfNodeType type) { type_ = type; }
  
  /*!
  \brief Return pointer to adjacent treenodes
  
  */
  virtual inline vector<RCP<BinaryTreeSelfNode> >  AdjacentTreenodes() { return adjacentTreenodes_; }
  
  /*!
  \brief set adjacent Treefnodes
  
  */
  virtual inline void SetAdjacentTnodes(vector<RCP<BinaryTreeSelfNode> > adjTnodes)
  { adjacentTreenodes_ = adjTnodes; }
  
  /*!
  \brief Return list of endnodes
  
  */
  virtual inline vector<int> Endnodes() { return endnodes_; }
  
  /*!
  \brief Set list of endnodes
  
  */
  virtual inline void SetEndnodes(vector<int> endnodes) { endnodes_= endnodes; }
   
  /*!
  \brief Update list of endnodes with endnodes of children
  
  */
  virtual void UpdateEndnodes() ;
  
  /*!
  \brief Return pointer to right child
  
  */
  virtual inline RCP<BinaryTreeSelfNode> Rightchild() const { return rightchild_; }
      
  /*!
  \brief Return pointer to left child
  
  */
  virtual inline RCP<BinaryTreeSelfNode> Leftchild() const { return leftchild_; }
  
  /*!
  \brief set children of a Binary Tree Node
  
  */
  void SetChildren(RCP<BinaryTreeSelfNode> leftchild,RCP<BinaryTreeSelfNode> rightchild);
  
  /*!
  \brief Return pointer to parent
  
  */
  virtual inline RCP<BinaryTreeSelfNode> Parent() const { return parent_; }
  
  /*!
  \brief set parent of treenode
  
  */
  virtual inline void SetParent(RCP<BinaryTreeSelfNode> parent) { parent_=parent; }
  
  /*!
  \brief Return pointer to slabs of DOP
  
  */
  virtual inline Epetra_SerialDenseMatrix Slabs() const { return slabs_; }
  
  /*!
  \brief Return pointer to nomals of DOP
  
  */
  virtual inline Epetra_SerialDenseMatrix Dopnormals() const { return dopnormals_; }
  
  /*!
  \brief Return no. of vertices
  
  */
  virtual inline int Kdop() const { return kdop_; }
  
  /*!
  \brief Return layer of current TreeNode
  
  */
  virtual inline int Layer(){ return layer_; }
  
  /*!
  \brief Set layer of current TreeNode
  
  */
  virtual inline void SetLayer(int layer) { layer_ = layer; }
  
  /*!
  \brief Return pointer to element list of TreeNode
  
  */
  virtual inline vector<int> Elelist() const { return elelist_; }
  
  /*!
  \brief Return dim of Problem
  
  */
  virtual inline const int& Dim() { return dim_; }
  
  //@}
    
private:
  
  //! type of BinaryTreeSelfNode
  BinaryTreeSelfNodeType type_;
  
  //! contact interface discretization
  DRT::Discretization& idiscret_;
  
  //! pointer to the parent BinaryTreeSelfNode
  RCP<BinaryTreeSelfNode> parent_; 

  //! pointer to the left child TreeNode
  RCP<BinaryTreeSelfNode> leftchild_;
  
  //! pointer to the right child TreeNode
  RCP<BinaryTreeSelfNode> rightchild_;
  
  //! geometry slabs of TreeNode, saved as Min|Max
  Epetra_SerialDenseMatrix slabs_;
  
  //! logical array of qualified sample vectors of current tree node
  vector<bool> qualifiedvectors_;
  
  //! vector with global IDs of end nodes of a surface (2D), -1 if there are no end-nodes
  vector<int> endnodes_ ;
       
  //! list with all elements of BinaryTreeSelfNode
  vector<int> elelist_;
  
  //! vector pointers to adjacent treenodes on the same layer
  vector<RCP<BinaryTreeSelfNode> >  adjacentTreenodes_;
  
  //! reference to DOP normals
  const Epetra_SerialDenseMatrix& dopnormals_;
  
  //! reference to sample vectors
  const Epetra_SerialDenseMatrix& samplevectors_;
  
  //! reference to no. of vertices
  const int& kdop_;
  
  //! reference to dim. of problem
  const int& dim_;
  
  //! reference to number of sample vectors
  const int& nvectors_;
  
  //! layer of TreeNode in Tree (0=rootnode!)
  int layer_;
 
  //! reference to storage scheme of all treenodes, sorted by layer
  vector<vector<RCP<BinaryTreeSelfNode> > > &  treenodes_;
  
  // relational operators for binary tree nodes
  
  //! operator <
  friend inline bool operator< (const RCP<BinaryTreeSelfNode> node1, const RCP<BinaryTreeSelfNode> node2) 
  {
    if (node1->elelist_.size() < node2->elelist_.size())
      return true;
    else if (node1->elelist_.size() == node2->elelist_.size() and node1->elelist_[0] < node2->elelist_[0])
      return true;
    else 
      return false;
  }
  
  //! operator >
  friend inline bool operator> (const RCP<BinaryTreeSelfNode> node1, const RCP<BinaryTreeSelfNode> node2) 
  {
    return operator< (node2,node1);
  }
  
  //! operator <=
  friend inline bool operator<= (const RCP<BinaryTreeSelfNode> node1, const RCP<BinaryTreeSelfNode> node2) 
  {
    return !operator>(node1,node2);
  }
  
  //! operator >=
  friend inline bool operator>= (const RCP<BinaryTreeSelfNode> node1, const RCP<BinaryTreeSelfNode> node2) 
  {
    return !operator< (node1,node2);
  }
  
  //! operator ==
  friend inline bool operator== (const RCP<BinaryTreeSelfNode> node1, const RCP<BinaryTreeSelfNode> node2)
  {
    if (node1->elelist_.size() != node2->elelist_.size())
      return false;
    else if (node1->elelist_[0] == node2->elelist_[0])
      return true;
    else 
      return false;
  }
  
  //! operator !=
  friend inline bool operator!= (const RCP<BinaryTreeSelfNode> node1, const RCP<BinaryTreeSelfNode> node2) 
  {
    return !operator== (node1,node2);
  }

}; // class BinaryTreeSelfNode


/*!
\brief A class representing one edge of the dual graph for self contact search

Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

\author popp (popp@lnm.mw.tum.de)
*/
class DualEdge
{ 
public:
  
  //relational operators for dual edges
  
  //! operator ==
  friend inline bool operator== (const RCP<DualEdge> edge1, const RCP<DualEdge> edge2)
  {
    if ((edge1->node1_ == edge2->node1_) and (edge1->node2_ == edge2->node2_))
      return true;
    else if ((edge1->node2_ == edge2->node1_) and (edge1->node1_ == edge2->node2_))
      return true;
    else
      return false;
  }
  
  //! operator !=
  friend inline bool operator!= (const RCP<DualEdge> edge1, const RCP<DualEdge> edge2) 
  {
    return !operator== (edge1,edge2);
  }
  
  //! operator <
  friend inline bool operator< (const RCP<DualEdge> edge1, const RCP<DualEdge> edge2) 
  {
    if (edge1->costs_ < edge2->costs_)
      return true;
    else if (edge1->costs_ > edge2->costs_)
      return false;
    else if (edge1 != edge2)
    {
      if (edge1->GreaterNode() < edge2->GreaterNode()) 
        return true;
      else if (edge1->GreaterNode() == edge2->GreaterNode())
      {
        if (edge1->LesserNode() < edge2->LesserNode())
          return true;
        else 
          return false;
      }
      else
        return false;
    }
    else 
      return false;
  }
  
  //! operator >
  friend inline bool operator> (const RCP<DualEdge> edge1, const RCP<DualEdge> edge2) 
  {
    return operator< (edge2,edge1);
  }

  //! operator <=
  friend inline bool operator<= (const RCP<DualEdge> edge1, const RCP<DualEdge> edge2) 
  {
    return !operator> (edge1,edge2);
  }
  
  //! operator >=
  friend inline bool operator>= (const RCP<DualEdge> edge1, const RCP<DualEdge> edge2) 
  {
    return !operator< (edge1,edge2);
  }

  /*!
  \brief Constructor of a dual edge

  */
  DualEdge (RCP<BinaryTreeSelfNode> node1_, RCP<BinaryTreeSelfNode> node2_,const int& dim);
 
  /*!
  \brief Destructor

  */
  virtual ~DualEdge() {};
 
  /*!
  \brief Calculate the cost function of a dual edge

  */
  void CalculateCosts();
 
  /*!
  \brief Return costs

  */
  virtual inline double Costs() const { return costs_; }
  
  /*!
  \brief Return first node of dual edge

  */
  virtual inline RCP<BinaryTreeSelfNode> GetNode1() const { return node1_; }
  
  /*!
  \brief Return second node of dual edge

  */
  virtual inline RCP<BinaryTreeSelfNode> GetNode2() const { return node2_; }
  
     
  /*!
  \brief return greater node of dual edge

  */
  virtual inline RCP<BinaryTreeSelfNode> GreaterNode()
  {
    if (node1_ > node2_)       return node1_;
    else if (node2_ > node1_)  return node2_;
    else                       return node1_;
  }
  
  virtual inline RCP<BinaryTreeSelfNode> LesserNode()
  {
    if (node1_ > node2_)       return node2_;
    else if (node2_ > node1_)  return node1_;
    else                       return node1_;
  }
  
  /*!
  \brief Return common tree node of two dual edges

  */
  virtual inline RCP<BinaryTreeSelfNode> CommonNode(RCP<DualEdge> treenode)
  {
    RCP<BinaryTreeSelfNode> node1 = treenode->GetNode1();
    RCP<BinaryTreeSelfNode> node2 = treenode->GetNode2();
    
    if (GetNode1() == node1 or GetNode2() == node1)
      return node1;
    else if (GetNode1() == node2 or GetNode2() == node2)
       return node2;
    else
      return null;
  }
  
private:
 
  //! first node of dual edge
  RCP<BinaryTreeSelfNode> node1_;
  
  //! second node of dual edge
  RCP<BinaryTreeSelfNode> node2_;
  
  //! cost function value fo dual edge
  double costs_;
 
  //! reference to dim. of problem
  const int& dim_;

}; // class DualEdge


/*!
\brief A class for performing self contact search in 2D / 3D based
       on a binary search tree and dual graphs
       
Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

\author popp (popp@lnm.mw.tum.de)
*/

class BinaryTreeSelf
{
public:
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  For now, we only consider the serial case (1 processor) here!!!
  
  \param discret (in):    The contact interface discretization
  \param elements (in):   All elements on self contact interface (fully overlapping map)
  \param dim (in):        The problem dimension 
  
  */
  BinaryTreeSelf(DRT::Discretization& discret,
                 RCP<Epetra_Map> elements,
                 int dim, double eps);

  /*!
  \brief Destructor

  */
  virtual ~BinaryTreeSelf() {};
  
  //! @name Access methods
  
  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return idiscret_.Comm(); }
  
  /*!
  \brief Get discretization of the interface

  */
  virtual inline DRT::Discretization& Discret() const { return idiscret_; }
  
  /*!
  \brief Get matrix of DOP normals

  */
  virtual inline const Epetra_SerialDenseMatrix& DopNormals() { return dopnormals_; }
  
  /*!
  \brief Get matrix of sample vectors

  */
  virtual inline const Epetra_SerialDenseMatrix& SampleVectors() { return samplevectors_; }
  
  /*!
  \brief Get number of vertices of DOP

  */
  virtual inline const int& Kdop() { return kdop_; }
  
  /*!
  \brief Return dim of Problem

  */
  virtual inline const int& Dim() { return dim_; }
  
  /*!
  \brief Return no. of sample vectors

  */
  virtual inline const int& Nvectors() { return nvectors_; }
  
  /*!
  \brief Return reference to storage scheme of all treenodes

  */
  virtual inline vector<vector<RCP<BinaryTreeSelfNode> > >& Treenodes() { return treenodes_; }
  
  /*!
  \brief Set the vector of adjacent tree nodes for leaf-nodes in the lowest layer

  */
  void CalculateAdjacentLeaves();
  
  /*!
  \brief Calculate the vector of adjacent tree nodes of inner tree nodes

  */
  void CalculateAdjacentTnodes();
 
  /*!
  \brief Return factor "eps" to enlarge dops

  */
  virtual inline double & Eps() { return eps_; }
  
  /*!
  \brief Return factor "enlarge" to enlarge dops

  */
  virtual inline double & Enlarge() { return enlarge_; }
    
  //@}
  
  //! @name Evaluation methods
 
  /*!
  \brief Calculate minimal element length / inflation factor "enlarge"

  */
  void SetEnlarge(bool isinit);
   
  /*!
  \brief Initialize Tree in a bottom up way based on dual graph

  */ 
  void InitializeTreeBottomUp(map<RCP<DualEdge>,vector <RCP<DualEdge> > > *dualGraph);

  /*!
  \brief Evaluate Binary search tree for combined search and update

  */
  void SearchContactCombined();
  
  /*!
  \brief Evaluate Binary search tree for self contact search 

  */ 
  void SearchSelfContact(RCP<BinaryTreeSelfNode> treenode);
  
  /*!
  \brief Find contact of adjacent surfaces

  */ 
  void EvaluateContactAndAdjacency(RCP<BinaryTreeSelfNode> treenode1,
                                   RCP<BinaryTreeSelfNode> treenode2,
                                   bool isadjacent);
  
  /*!
  \brief Test for adjacency (2D)

  */ 
  bool TestAdjacent2D (RCP<BinaryTreeSelfNode> treenode1,
                       RCP<BinaryTreeSelfNode> treenode2);
   
  /*!
  \brief Test for adjacency (2D)
      
  */  
  bool TestAdjacent3D (RCP<BinaryTreeSelfNode> treenode1,
                       RCP<BinaryTreeSelfNode> treenode2);
   
  /*!
  \brief Master/Slave sorting for selfcontact
      
  */ 
  void MasterSlaveSorting(int eleID,bool isslave);
   
  /*!
  \brief Update normals and qualified sample vectors of the whole tree
      
  */ 
  void UpdateNormals();
   
  /*!
  \brief Update normals and qualified sample vectors of one treenode
      
  */ 
  void UpdateQualifiedVectors(RCP<BinaryTreeSelfNode> treenode);

  //@}
      
protected:
    
  // don't want = operator and cctor
  BinaryTreeSelf operator = (const BinaryTreeSelf& old);
  BinaryTreeSelf(const BinaryTreeSelf& old);

  DRT::Discretization&                       idiscret_;       // contact interface discretization
  RCP<Epetra_Map>                            elements_;       // all contact elements on surface (full map)
  
  const int                                  dim_;            // problem dimension (2D or 3D)
  int                                        kdop_;           // sets k for DOP (8 for 2D, 18 for 3D)   
  int                                        nvectors_;       // defining number of sample vectors
  Epetra_SerialDenseMatrix                   dopnormals_;     // defining vertices of DOP
  Epetra_SerialDenseMatrix                   samplevectors_;  // defining sample vectors
  double                                     eps_;            // epsilon for enlarging dops (of user)
  double                                     enlarge_;        // to enlarge dops
  
  RCP<BinaryTreeSelfNode>                    root_;           // root treenode
  vector<vector<RCP<BinaryTreeSelfNode> > >  treenodes_;      // storage of all treenodes, sorted by layers
  map<int,vector<int> >                      contactpairs_;   // contacting elements (key: GID, argument: GIDs of partners)
  map<int,vector<RCP<BinaryTreeSelfNode> > > adjacencymatrix_;// map of adjacent elements, sorted by global id (only needed in 3D)
  map<int, RCP<BinaryTreeSelfNode> >         leafsmap_;       // map of all leaf nodes, sorted by global id
     
}; // class BinaryTreeSelf
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_BINARYTREE_H
#endif  // #ifdef CCADISCRET
