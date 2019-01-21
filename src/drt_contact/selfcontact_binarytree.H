/*!----------------------------------------------------------------------
\file selfcontact_binarytree.H

\brief Search tree for self-contact problems

\level 2

\maintainer Alexander Popp

*-----------------------------------------------------------------------*/
#ifndef SELFCONTACT_BINARYTREE_SELF_H
#define SELFCONTACT_BINARYTREE_SELF_H

#include "Epetra_Map.h"
#include "../drt_mortar/base_binarytree.H"

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
    SELFCO_NO_ELEMENTS,  ///< indicates that there are no elements on this (root) treenode
    SELFCO_UNDEFINED     ///< indicates an undefined tree node
  };

  //@}

  /*!
  \brief A class representing one tree node of the binary tree for self contact

  Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

  \author popp (popp@lnm.mw.tum.de)
  */
  class SelfBinaryTreeNode : public MORTAR::BaseBinaryTreeNode
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
    \param kdop           reference to no. of vertices
    \param dim            dimension of problem
    \param layer          current layer of tree node
    \param treenodes      references to tree nodes storage scheme

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
    \brief Update slabs of current tree node in bottom up way

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
    std::vector<bool> QualifiedVectors() const { return qualifiedvectors_; }

    /*!
    \brief Set logical array of qualified sample vectors

    */
    void SetQualifiedVectors(std::vector<bool> qualifiedvectors)
    {
      qualifiedvectors_ = qualifiedvectors;
    }

    /*!
    \brief Print type of tree node to std::cout

    */
    virtual void PrintType() final;

    /*!
    \brief Set slabs of current treenode

    */
    void SetSlabs(Epetra_SerialDenseMatrix& newslabs);
    //@}

    //! @name Access methods

    /*!
    \brief Get communicator

    */
    const Epetra_Comm& Comm() const;


    /*!
    \brief Complete tree by filling tree node storage scheme

    */
    void CompleteTree(int layer, double& enlarge);

    /*!
    \brief Return pointer to type of treenode

    */
    SelfBinaryTreeNodeType Type() const { return type_; }

    /*!
    \brief Set type of tree node

    */
    void SetType(SelfBinaryTreeNodeType type) { type_ = type; }

    /*!
    \brief Return pointer to adjacent tree nodes

    */
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> AdjacentTreenodes() { return adjacentTreenodes_; }

    /*!
    \brief set adjacent tree nodes

    */
    void SetAdjacentTnodes(std::vector<Teuchos::RCP<SelfBinaryTreeNode>> adjTnodes)
    {
      adjacentTreenodes_ = adjTnodes;
    }

    /*!
    \brief Return list of endnodes

    */
    std::vector<int> Endnodes() const { return endnodes_; }

    /*!
    \brief Set list of endnodes

    */
    void SetEndnodes(std::vector<int> endnodes) { endnodes_ = endnodes; }

    /*!
    \brief Update list of endnodes with endnodes of children

    */
    void UpdateEndnodes();

    /*!
    \brief Return pointer to right child

    */
    Teuchos::RCP<SelfBinaryTreeNode> Rightchild() const { return rightchild_; }

    /*!
    \brief Return pointer to left child

    */
    Teuchos::RCP<SelfBinaryTreeNode> Leftchild() const { return leftchild_; }

    /*!
    \brief set children of a Binary Tree Node

    */
    void SetChildren(
        Teuchos::RCP<SelfBinaryTreeNode> leftchild, Teuchos::RCP<SelfBinaryTreeNode> rightchild);

    /*!
    \brief Return pointer to parent

    */
    Teuchos::RCP<SelfBinaryTreeNode> Parent() const { return parent_; }

    /*!
    \brief set parent of treenode

    */
    void SetParent(Teuchos::RCP<SelfBinaryTreeNode> parent) { parent_ = parent; }
    //@}

   private:
    // don't want = operator and cctor
    SelfBinaryTreeNode operator=(const SelfBinaryTreeNode& old);
    SelfBinaryTreeNode(const SelfBinaryTreeNode& old);

    //! type of SelfBinaryTreeNode
    SelfBinaryTreeNodeType type_;

    // the pointers to the parent as well as to the left and right child are not moved to the
    // BaseBinaryTreeNode as this would require a lot of dynamic casting and thereby complicating
    // the readability of the code
    //! pointer to the parent SelfBinaryTreeNode
    Teuchos::RCP<SelfBinaryTreeNode> parent_;

    //! pointer to the left child TreeNode
    Teuchos::RCP<SelfBinaryTreeNode> leftchild_;

    //! pointer to the right child TreeNode
    Teuchos::RCP<SelfBinaryTreeNode> rightchild_;

    //! logical array of qualified sample vectors of current tree node
    std::vector<bool> qualifiedvectors_;

    //! vector with global IDs of end nodes of a surface (2D), -1 if there are no end-nodes
    std::vector<int> endnodes_;

    //! vector pointers to adjacent treenodes on the same layer
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> adjacentTreenodes_;

    //! reference to sample vectors
    const Epetra_SerialDenseMatrix& samplevectors_;

    //! reference to number of sample vectors
    const int& nvectors_;

    //! reference to storage scheme of all treenodes, sorted by layer
    std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& treenodes_;

    // relational operators for binary tree nodes

    //! operator <
    friend bool operator<(
        const Teuchos::RCP<SelfBinaryTreeNode> node1, const Teuchos::RCP<SelfBinaryTreeNode> node2)
    {
      if (node1->Elelist().size() < node2->Elelist().size())
        return true;
      else if (node1->Elelist().size() == node2->Elelist().size() and
               node1->Elelist()[0] < node2->Elelist()[0])
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
      if (node1->Elelist().size() != node2->Elelist().size())
        return false;
      else if (node1->Elelist()[0] == node2->Elelist()[0])
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
    double Costs() const { return costs_; }

    /*!
    \brief Return first node of dual edge

    */
    Teuchos::RCP<SelfBinaryTreeNode> GetNode1() const { return node1_; }

    /*!
    \brief Return second node of dual edge

    */
    Teuchos::RCP<SelfBinaryTreeNode> GetNode2() const { return node2_; }

    /*!
    \brief Return common tree node of two dual edges

    */
    Teuchos::RCP<SelfBinaryTreeNode> CommonNode(Teuchos::RCP<SelfDualEdge> treenode)
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
    Teuchos::RCP<SelfBinaryTreeNode> GreaterNode()
    {
      if (node1_ > node2_)
        return node1_;
      else if (node2_ > node1_)
        return node2_;
      else
        return node1_;
    }

    Teuchos::RCP<SelfBinaryTreeNode> LesserNode()
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

  class SelfBinaryTree : public MORTAR::BaseBinaryTree
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
    \brief Evaluate search self binary tree:
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
    const Epetra_Comm& Comm() const;

    /*!
    \brief Get local interface communicator (only contains participating processors)

    */
    const Epetra_Comm* lComm() const { return lcomm_; }

    /*!
    \brief Get matrix of sample vectors

    */
    const Epetra_SerialDenseMatrix& SampleVectors() const { return samplevectors_; }

    /*!
    \brief Return no. of sample vectors

    */
    const int& Nvectors() const { return nvectors_; }

    /*!
    \brief Return reference to storage scheme of all tree nodes

    */
    std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& Treenodes() { return treenodes_; }

    /*!
    \brief Set the vector of adjacent tree nodes for leaf-nodes in the lowest layer

    */
    void CalculateAdjacentLeaves();

    /*!
    \brief Calculate the vector of adjacent tree nodes of inner tree nodes

    */
    void CalculateAdjacentTnodes();
    //@}

    //! @name Evaluation methods

    /*!
    \brief Calculate the dual graph

    \param [out] dualgraph:  construction of binary tree is based on this data
    \param [in]  elelist:    list (gids) of all contact elements of the surface

     */
    void CalculateDualGraph(
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualGraph,
        const std::vector<int>& elelist);

    /*!
    \brief Calculate number of slabs intersections of treenode 1 and 2

    \param [in] treenode1:  self binary tree node
    \param [in] treenode2:  self binary tree node

     */
    int CalculateSlabsIntercepts(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);

    /*!
    \brief Initialize Tree in a bottom up way based on dual graph

    \param [in] dualgraph:  construction of binary tree is based on this data

    */
    void InitializeTreeBottomUp(
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualGraph);

    /*!
    \brief Initialize internal variables

     */
    void InitInternalVariables() final;

    /*!
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    void SetEnlarge() final;

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

    const Epetra_Comm* lcomm_;           // local interface communicator
    Teuchos::RCP<Epetra_Map> elements_;  // all contact elements on surface (full map)

    int nvectors_;                            // defining number of sample vectors
    Epetra_SerialDenseMatrix samplevectors_;  // defining sample vectors

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
