/*!----------------------------------------------------------------------
\file selfcontact_binarytree.H

\brief Search tree for self-contact problems

\level 2

\maintainer Alexander Popp

*-----------------------------------------------------------------------*/
#ifndef SELFCONTACT_BINARYTREE_SELF_H
#define SELFCONTACT_BINARYTREE_SELF_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_mortar/abstract_binarytree.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CONTACT
{
  // forward declarations

  //! @name Enums and Friends
  enum SelfBinaryTreeNodeType
  {
    SELFCO_INNER,        ///< indicates an inner node (has children)
    SELFCO_LEAF,         ///< indicates a leaf node (no further children)
    SELFCO_NO_ELEMENTS,  ///< inidicates that there are no elements on this (root) treenode
    SELFCO_UNDEFINED     ///< indiactes an undefined tree node
  };

  //@}

  /*!
  \brief A class representing one tree node of the binary tree for self contact

  Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

  \author popp (popp@lnm.mw.tum.de)
  */
  class SelfBinaryTreeNode
  {
   public:
    /*!
    \brief Constructor of a tree node

    \param type           type of SelfBinaryTreeNode
    \param discret        contact interface discretization
    \param parent         points to parent tree node
    \param elelist        list of all elements in SelfBinaryTreeNode
    \param dopnormals     reference to DOP normals
    \param samplevectors  reference to sample vectors
    \param kdop           refernece to no. of vertices
    \param dim            dimension of problem
    \param layer          current layer of treenode
    \param treenodes      references to treenodes storage scheme

    */
    SelfBinaryTreeNode(SelfBinaryTreeNodeType type, DRT::Discretization& discret,
        Teuchos::RCP<SelfBinaryTreeNode> parent, std::vector<int> elelist,
        const Epetra_SerialDenseMatrix& dopnormals, const Epetra_SerialDenseMatrix& samplevectors,
        const int& kdop, const int& dim, const int& nvectors, const int layer,
        std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& treenodes);

    /*!
    \brief Destructor

    */
    virtual ~SelfBinaryTreeNode(){};

    //! @name Evaluation methods

    /*!
    \brief Calculate slabs of dop

    */
    void CalculateSlabsDop();

    /*!
    \brief Update slabs of current treenode in bottom up way

    */
    void UpdateSlabsBottomUp(double& eps);

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
    virtual std::vector<bool> QualifiedVectors() { return qualifiedvectors_; }

    /*!
    \brief Set logical array of qualified sample vectors

    */
    virtual void SetQualifiedVectors(std::vector<bool> qualifiedvectors)
    {
      qualifiedvectors_ = qualifiedvectors;
    }

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
    virtual const Epetra_Comm& Comm() const;


    /*!
    \brief Complete tree by filling treenode storage scheme

    */
    void CompleteTree(int layer, double& enlarge);

    /*!

    /brief Get discretization of the interface

    */
    virtual DRT::Discretization& Discret() const { return idiscret_; }

    /*!
    \brief Return pointer to type of treenode

    */
    virtual SelfBinaryTreeNodeType Type() const { return type_; }

    /*!
    \brief Set type of treenode

    */
    virtual void SetType(SelfBinaryTreeNodeType type) { type_ = type; }

    /*!
    \brief Return pointer to adjacent treenodes

    */
    virtual std::vector<Teuchos::RCP<SelfBinaryTreeNode>> AdjacentTreenodes()
    {
      return adjacentTreenodes_;
    }

    /*!
    \brief set adjacent Treefnodes

    */
    virtual void SetAdjacentTnodes(std::vector<Teuchos::RCP<SelfBinaryTreeNode>> adjTnodes)
    {
      adjacentTreenodes_ = adjTnodes;
    }

    /*!
    \brief Return list of endnodes

    */
    virtual std::vector<int> Endnodes() { return endnodes_; }

    /*!
    \brief Set list of endnodes

    */
    virtual void SetEndnodes(std::vector<int> endnodes) { endnodes_ = endnodes; }

    /*!
    \brief Update list of endnodes with endnodes of children

    */
    virtual void UpdateEndnodes();

    /*!
    \brief Return pointer to right child

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> Rightchild() const { return rightchild_; }

    /*!
    \brief Return pointer to left child

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> Leftchild() const { return leftchild_; }

    /*!
    \brief set children of a Binary Tree Node

    */
    void SetChildren(
        Teuchos::RCP<SelfBinaryTreeNode> leftchild, Teuchos::RCP<SelfBinaryTreeNode> rightchild);

    /*!
    \brief Return pointer to parent

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> Parent() const { return parent_; }

    /*!
    \brief set parent of treenode

    */
    virtual void SetParent(Teuchos::RCP<SelfBinaryTreeNode> parent) { parent_ = parent; }

    /*!
    \brief Return pointer to slabs of DOP

    */
    virtual Epetra_SerialDenseMatrix Slabs() const { return slabs_; }

    /*!
    \brief Return pointer to nomals of DOP

    */
    virtual Epetra_SerialDenseMatrix Dopnormals() const { return dopnormals_; }

    /*!
    \brief Return no. of vertices

    */
    virtual int Kdop() const { return kdop_; }

    /*!
    \brief Return layer of current TreeNode

    */
    virtual int Layer() { return layer_; }

    /*!
    \brief Set layer of current TreeNode

    */
    virtual void SetLayer(int layer) { layer_ = layer; }

    /*!
    \brief Return pointer to element list of TreeNode

    */
    virtual std::vector<int> Elelist() const { return elelist_; }

    /*!
    \brief Return dim of Problem

    */
    virtual const int& Dim() { return dim_; }

    //@}

   private:
    // don't want = operator and cctor
    SelfBinaryTreeNode operator=(const SelfBinaryTreeNode& old);
    SelfBinaryTreeNode(const SelfBinaryTreeNode& old);

    //! type of SelfBinaryTreeNode
    SelfBinaryTreeNodeType type_;

    //! contact interface discretization
    DRT::Discretization& idiscret_;

    //! pointer to the parent SelfBinaryTreeNode
    Teuchos::RCP<SelfBinaryTreeNode> parent_;

    //! pointer to the left child TreeNode
    Teuchos::RCP<SelfBinaryTreeNode> leftchild_;

    //! pointer to the right child TreeNode
    Teuchos::RCP<SelfBinaryTreeNode> rightchild_;

    //! geometry slabs of TreeNode, saved as Min|Max
    Epetra_SerialDenseMatrix slabs_;

    //! logical array of qualified sample vectors of current tree node
    std::vector<bool> qualifiedvectors_;

    //! vector with global IDs of end nodes of a surface (2D), -1 if there are no end-nodes
    std::vector<int> endnodes_;

    //! list with all elements of SelfBinaryTreeNode
    std::vector<int> elelist_;

    //! vector pointers to adjacent treenodes on the same layer
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> adjacentTreenodes_;

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
    std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& treenodes_;

    // relational operators for binary tree nodes

    //! operator <
    friend bool operator<(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      if (node1->elelist_.size() < node2->elelist_.size())
        return true;
      else if (node1->elelist_.size() == node2->elelist_.size() and
               node1->elelist_[0] < node2->elelist_[0])
        return true;
      else
        return false;
    }

    //! operator >
    friend bool operator>(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      return operator<(node2, node1);
    }

    //! operator <=
    friend bool operator<=(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      return !operator>(node1, node2);
    }

    //! operator >=
    friend bool operator>=(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      return !operator<(node1, node2);
    }

    //! operator ==
    friend bool operator==(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      if (node1->elelist_.size() != node2->elelist_.size())
        return false;
      else if (node1->elelist_[0] == node2->elelist_[0])
        return true;
      else
        return false;
    }

    //! operator !=
    friend bool operator!=(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      return !operator==(node1, node2);
    }

  };  // class SelfBinaryTreeNode


  /*!
  \brief A class representing one edge of the dual graph for self contact search

  Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

  \author popp (popp@lnm.mw.tum.de)
  */
  class SelfDualEdge
  {
   public:
    // relational operators for dual edges

    //! operator ==
    friend bool operator==(
        const Teuchos::RCP<SelfDualEdge> edge1, const Teuchos::RCP<SelfDualEdge> edge2)
    {
      if ((edge1->node1_ == edge2->node1_) and (edge1->node2_ == edge2->node2_))
        return true;
      else if ((edge1->node2_ == edge2->node1_) and (edge1->node1_ == edge2->node2_))
        return true;
      else
        return false;
    }

    //! operator !=
    friend bool operator!=(
        const Teuchos::RCP<SelfDualEdge> edge1, const Teuchos::RCP<SelfDualEdge> edge2)
    {
      return !operator==(edge1, edge2);
    }

    //! operator <
    friend bool operator<(
        const Teuchos::RCP<SelfDualEdge> edge1, const Teuchos::RCP<SelfDualEdge> edge2)
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
    friend bool operator>(
        const Teuchos::RCP<SelfDualEdge> edge1, const Teuchos::RCP<SelfDualEdge> edge2)
    {
      return operator<(edge2, edge1);
    }

    //! operator <=
    friend bool operator<=(
        const Teuchos::RCP<SelfDualEdge> edge1, const Teuchos::RCP<SelfDualEdge> edge2)
    {
      return !operator>(edge1, edge2);
    }

    //! operator >=
    friend bool operator>=(
        const Teuchos::RCP<SelfDualEdge> edge1, const Teuchos::RCP<SelfDualEdge> edge2)
    {
      return !operator<(edge1, edge2);
    }

    /*!
    \brief Constructor of a dual edge

    */
    SelfDualEdge(Teuchos::RCP<SelfBinaryTreeNode> node1_, Teuchos::RCP<SelfBinaryTreeNode> node2_,
        const int& dim);

    /*!
    \brief Destructor

    */
    virtual ~SelfDualEdge(){};

    /*!
    \brief Return costs

    */
    virtual double Costs() const { return costs_; }

    /*!
    \brief Return first node of dual edge

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> GetNode1() const { return node1_; }

    /*!
    \brief Return second node of dual edge

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> GetNode2() const { return node2_; }

    /*!
    \brief Return common tree node of two dual edges

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> CommonNode(Teuchos::RCP<SelfDualEdge> treenode)
    {
      Teuchos::RCP<SelfBinaryTreeNode> node1 = treenode->GetNode1();
      Teuchos::RCP<SelfBinaryTreeNode> node2 = treenode->GetNode2();

      if (GetNode1() == node1 or GetNode2() == node1)
        return node1;
      else if (GetNode1() == node2 or GetNode2() == node2)
        return node2;
      else
        return Teuchos::null;
    }

   private:
    /*!
    \brief Calculate the cost function of a dual edge

    */
    void CalculateCosts();

    /*!
    \brief return greater node of dual edge

    */
    virtual Teuchos::RCP<SelfBinaryTreeNode> GreaterNode()
    {
      if (node1_ > node2_)
        return node1_;
      else if (node2_ > node1_)
        return node2_;
      else
        return node1_;
    }

    virtual Teuchos::RCP<SelfBinaryTreeNode> LesserNode()
    {
      if (node1_ > node2_)
        return node2_;
      else if (node2_ > node1_)
        return node1_;
      else
        return node1_;
    }

    // don't want = operator and cctor
    SelfDualEdge operator=(const SelfDualEdge& old);
    SelfDualEdge(const SelfDualEdge& old);

    //! first node of dual edge
    Teuchos::RCP<SelfBinaryTreeNode> node1_;

    //! second node of dual edge
    Teuchos::RCP<SelfBinaryTreeNode> node2_;

    //! cost function value fo dual edge
    double costs_;

    //! reference to dim. of problem
    const int& dim_;

  };  // class SelfDualEdge


  /*!
  \brief A class for performing self contact search in 2D / 3D based
         on a binary search tree and dual graphs

  Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

  \author popp (popp@lnm.mw.tum.de)
  */

  class SelfBinaryTree : public MORTAR::AbstractBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    Constructs an instance of this class.<br>
    For now, we only consider the serial case (1 processor) here!!!

    \param discret (in):    The contact interface discretization
    \param lcomm (in):      Local interface communicator
    \param elements (in):   All elements on self contact interface (fully overlapping map)
    \param dim (in):        The problem dimension

    */
    SelfBinaryTree(DRT::Discretization& discret, const Epetra_Comm* lcomm,
        Teuchos::RCP<Epetra_Map> elements, int dim, double eps);

    /*!
    \brief Destructor

    */
    virtual ~SelfBinaryTree(){};

    //! @name Evaluation methods

    /*!
    \brief Evaluate search self binarytree:
           call SetEnlarge() and SearchContactSeparate()

    */
    void EvaluateSearch() final;

    /*!
    \brief Initialize the binary tree

    */
    void Init() final;

    //@}

   private:
    //! @name Access methods

    /*!
    \brief Get communicator

    */
    virtual const Epetra_Comm& Comm() const;

    /*!
    \brief Get discretization of the interface

    */
    virtual DRT::Discretization& Discret() const { return idiscret_; }

    /*!
    \brief Get local interface communicator (only contains participating processors)

    */
    virtual const Epetra_Comm* lComm() const { return lcomm_; }

    /*!
    \brief Get matrix of DOP normals

    */
    virtual const Epetra_SerialDenseMatrix& DopNormals() { return dopnormals_; }

    /*!
    \brief Get matrix of sample vectors

    */
    virtual const Epetra_SerialDenseMatrix& SampleVectors() { return samplevectors_; }

    /*!
    \brief Get number of vertices of DOP

    */
    virtual const int& Kdop() { return kdop_; }

    /*!
    \brief Return dim of Problem

    */
    virtual const int& Dim() { return dim_; }

    /*!
    \brief Return no. of sample vectors

    */
    virtual const int& Nvectors() { return nvectors_; }

    /*!
    \brief Return reference to storage scheme of all treenodes

    */
    virtual std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& Treenodes()
    {
      return treenodes_;
    }

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
    virtual double& Eps() { return eps_; }

    /*!
    \brief Return factor "enlarge" to enlarge dops

    */
    virtual double& Enlarge() { return enlarge_; }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Initialize Tree in a bottom up way based on dual graph

    */
    void InitializeTreeBottomUp(
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualGraph);

    /*!
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    void SetEnlarge();

    /*!
    \brief Evaluate Binary search tree for separate update (bottom-up,
           on whole interface) and search

    This is the default version, because it outperforms the combined
    algorithm for large self contact problems.

    */
    void SearchContactSeparate();

    /*!
    \brief Evaluate Binary search tree for combined update (top-down,
           only where necessary) and search

    */
    void SearchContactCombined();

    /*!
    \brief Evaluate Binary search tree for self contact search

    */
    void SearchSelfContactSeparate(Teuchos::RCP<SelfBinaryTreeNode> treenode);
    void SearchSelfContactCombined(Teuchos::RCP<SelfBinaryTreeNode> treenode);

    /*!
    \brief Evaluate Binary search tree for contact search between separate roots
           (this is more or less identical to two-body contact search)

    */
    void SearchRootContactSeparate(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);
    void SearchRootContactCombined(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);

    /*!
    \brief Find contact of adjacent surfaces

    */
    void EvaluateContactAndAdjacency(Teuchos::RCP<SelfBinaryTreeNode> treenode1,
        Teuchos::RCP<SelfBinaryTreeNode> treenode2, bool isadjacent);

    /*!
    \brief Test for adjacency (2D)

    */
    bool TestAdjacent2D(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);

    /*!
    \brief Test for adjacency (2D)

    */
    bool TestAdjacent3D(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);

    /*!
    \brief Master/Slave sorting for selfcontact

    */
    void MasterSlaveSorting(int eleID, bool isslave);

    /*!
    \brief Update normals and qualified sample vectors of the whole tree

    */
    void UpdateNormals();


    // don't want = operator and cctor
    SelfBinaryTree operator=(const SelfBinaryTree& old);
    SelfBinaryTree(const SelfBinaryTree& old);

    DRT::Discretization& idiscret_;      // contact interface discretization
    const Epetra_Comm* lcomm_;           // local interface communicator
    Teuchos::RCP<Epetra_Map> elements_;  // all contact elements on surface (full map)

    const int dim_;                           // problem dimension (2D or 3D)
    int kdop_;                                // sets k for DOP (8 for 2D, 18 for 3D)
    int nvectors_;                            // defining number of sample vectors
    Epetra_SerialDenseMatrix dopnormals_;     // defining vertices of DOP
    Epetra_SerialDenseMatrix samplevectors_;  // defining sample vectors
    double eps_;                              // epsilon for enlarging dops (of user)
    double enlarge_;                          // to enlarge dops

    // root treenodes
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> roots_;

    // storage of all treenodes, sorted by layers
    std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>> treenodes_;

    // storage of all treenodes, sorted by layers
    std::map<int, std::vector<int>> contactpairs_;

    // map of adjacent elements, sorted by global id (only needed in 3D)
    std::map<int, std::vector<Teuchos::RCP<SelfBinaryTreeNode>>> adjacencymatrix_;

    // map of all leaf nodes, sorted by global id
    std::map<int, Teuchos::RCP<SelfBinaryTreeNode>> leafsmap_;

  };  // class SelfBinaryTree
}  // namespace CONTACT


#endif  // #ifndef DSELFCONTACT_BINARYTREE_H
