/*!----------------------------------------------------------------------
\file selfcontact_binarytree.H

\brief Search tree for self-contact problems

\level 2

\maintainer Alexander Popp

*-----------------------------------------------------------------------*/
#ifndef SELFCONTACT_BINARYTREE_SELF_H
#define SELFCONTACT_BINARYTREE_SELF_H

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "../drt_lib/drt_element.H"
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
    \param nonsmoothsurf  non-smooth self contact surface
    \param treenodes      references to tree nodes storage scheme

    */
    SelfBinaryTreeNode(SelfBinaryTreeNodeType type, DRT::Discretization& discret,
        Teuchos::RCP<SelfBinaryTreeNode> parent, std::vector<int> elelist,
        const Epetra_SerialDenseMatrix& dopnormals, const Epetra_SerialDenseMatrix& samplevectors,
        const int& kdop, const int& dim, const int& nvectors, const int layer,
        const bool nonsmoothsurf,
        std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& treenodes);

    /*!
    \brief Destructor

    */
    virtual ~SelfBinaryTreeNode(){};

    //! @name Evaluation methods

    /*!
    \brief Update slabs of current tree node in bottom up way

    */
    void UpdateSlabsBottomUp(double& eps) final;

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
    \brief Print type of tree node to std::cout

    */
    void PrintType() final;
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
    \brief set parent of tree node

    */
    void SetParent(Teuchos::RCP<SelfBinaryTreeNode> parent) { parent_ = parent; }

    /*!
    \brief Return owner of current tree node

    */
    int Owner() const { return owner_; }

    /*!
    \brief set owner of tree node

    */
    void SetOwner(int treenodeowner) { owner_ = treenodeowner; }

    /*!
    \brief set owner of parent according to owner of children

    */
    void SetParentOwner(int leftchildowner, int rightchildowner);
    //@}

   private:
    // don't want = operator and cctor
    SelfBinaryTreeNode operator=(const SelfBinaryTreeNode& old) = delete;
    SelfBinaryTreeNode(const SelfBinaryTreeNode& old) = delete;

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

    //! owner of the tree node
    int owner_;

    //! tree node is part of a non-smooth self contact surface
    const bool nonsmoothsurf_;

    //! reference to storage scheme of all tree nodes, sorted by layer
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
    SelfDualEdge operator=(const SelfDualEdge& old) = delete;
    SelfDualEdge(const SelfDualEdge& old) = delete;

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
    \param iparams (in):    interface specific parameter list
    \param elements (in):   All elements on self contact interface (fully overlapping map)
    \param dim (in):        The problem dimension

    */
    SelfBinaryTree(DRT::Discretization& discret, const Epetra_Comm* lcomm,
        const Teuchos::ParameterList& iparams, Teuchos::RCP<Epetra_Map> elements, int dim,
        double eps);

    /*!
    \brief Destructor

    */
    virtual ~SelfBinaryTree(){};

    //! @name Evaluation methods

    /*!
    \brief Evaluate search self binary tree: call SetEnlarge() and SearchContact()

    */
    void EvaluateSearch() final;

    /*!
    \brief Initialize the self binary tree

    */
    void Init() override;

   protected:
    /*!
    \brief Initialize the leaf nodes and related map

    \param elelist:  gids of all contact elements of current surface

    */
    void InitLeafNodesAndMap(std::vector<int>& elelist);
    //@}

    //! @name Access methods

    /*!
    \brief Get write access to the adjacency matrix

    */
    std::map<int, std::vector<Teuchos::RCP<SelfBinaryTreeNode>>>& SetAdjacencymatrix()
    {
      return adjacencymatrix_;
    }

    /*!
    \brief Get communicator

    */
    const Epetra_Comm& Comm() const;

    /*!
    \brief Get access to the contact pairs

    */
    std::map<int, std::vector<int>> ContactPairs() const { return contactpairs_; }

    /*!
    \brief Get write access to the contact pairs

    */
    std::map<int, std::vector<int>>& SetContactPairs() { return contactpairs_; }

    /*!
    \brief Get local interface communicator (only contains participating processors)

    */
    const Epetra_Comm* lComm() const { return lcomm_; }

    /*!
    \brief Get map of leaf nodes

    */
    std::map<int, Teuchos::RCP<SelfBinaryTreeNode>> Leafsmap() const { return leafsmap_; }

    /*!
    \brief Get write access to map of leaf nodes

    */
    std::map<int, Teuchos::RCP<SelfBinaryTreeNode>>& SetLeafsmap() { return leafsmap_; }

    /*!
    \brief Return no. of sample vectors

    */
    const int& Nvectors() const { return nvectors_; }

    /*!
    \brief Get root nodes

    */
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> Roots() const { return roots_; }

    /*!
    \brief Get write access to root nodes

    */
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>>& SetRoots() { return roots_; }

    /*!
    \brief Get matrix of sample vectors

    */
    const Epetra_SerialDenseMatrix& SampleVectors() const { return samplevectors_; }

    /*!
    \brief Return reference to storage scheme of all tree nodes

    */
    std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>> Treenodes() const
    {
      return treenodes_;
    }
    //@}

    //! @name Evaluation methods

    /*!
    \brief Decide whether tree nodes need to be added to contact pairs

    \param [in]  treenode1:       first tree node
    \param [in]  treenode2:       second tree node

    */
    virtual void AddTreeNodesToContactPairs(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);
    /*!
    \brief Set the vector of adjacent tree nodes for leaf-nodes in the lowest layer

    */
    void CalculateAdjacentLeaves();

    /*!
    \brief Calculate the vector of adjacent tree nodes of inner tree nodes

    */
    void CalculateAdjacentTnodes();

    /*!
    \brief Calculate the adjacent tree nodes and adjacent dual edges of current element

    \param [in]  possadjids:   vector of global IDs of possible adjacent elements
    \param [in]  gid:          global id of current element
    \param [in]  adjElementk:  k-th adjacent element of current element
    \param [in]  node1:        leaf node with global id "gid"
    \param [out] adjtreenodes: vector of adjacent tree nodes (elements) of current element
    \param [out] adjdualedges: vector of adjacent dual edges containing current element

     */
    void CalculateAdjacentTreeNodesAndDualEdges(std::vector<int>& possadjids, const int gid,
        DRT::Element* adjElementk, Teuchos::RCP<SelfBinaryTreeNode>& node1,
        std::vector<Teuchos::RCP<SelfBinaryTreeNode>>& adjtreenodes,
        std::vector<Teuchos::RCP<SelfDualEdge>>& adjdualedges);

    /*!
    \brief Get the mortar element specific number of first order nodes

    \param [in] element:  element of which number of first order nodes shall be determined

    */
    int GetEleSpecificNumNodes(DRT::Element* element);

    /*!
    \brief Get the (contracted) node that combines the nodes of the contracted edge

    \param [in]      contractedEdge:  dual edge that is contracted
    \param [in,out]  contractedNode:  node that consists of both nodes of contracted edge

    */
    virtual void GetContractedNode(Teuchos::RCP<SelfDualEdge>& contractedEdge,
        Teuchos::RCP<SelfBinaryTreeNode>& contractedNode);

    /*!
    \brief Initialize internal variables

     */
    void InitInternalVariables() final;

    /*!
    \brief Master/Slave sorting for self contact

    */
    void MasterSlaveSorting(int eleID, bool isslave);

    /*!
    \brief Evaluate Binary search tree for self contact search

    */
    void SearchSelfContact(Teuchos::RCP<SelfBinaryTreeNode> treenode);

    /*!
    \brief Evaluate Binary search tree for contact search between separate roots
           (this is more or less identical to two-body contact search)

    */
    void SearchRootContact(
        Teuchos::RCP<SelfBinaryTreeNode> treenode1, Teuchos::RCP<SelfBinaryTreeNode> treenode2);

    /*!
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    void SetEnlarge() final;

    /*!
    \brief Update the dual graph and determine the root nodes

    \param [in] contractedEdge: contracted dual edge
    \param [in] adjEdges:       vector of all adjacent edges of the contracted edge
    \param [in] newnode:        node that represents the contracted self dual edge
    \param [in/out] dualgraph:  construction of binary tree is based on this data

    */
    void UpdateDualGraph(Teuchos::RCP<SelfDualEdge>& contractedEdge,
        std::vector<Teuchos::RCP<SelfDualEdge>>& adjEdges,
        Teuchos::RCP<SelfBinaryTreeNode>& newNode,
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualgraph);

    /*!
    \brief Update normals and qualified sample vectors of the whole tree

    */
    void UpdateNormals();

   private:
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
    \brief Evaluate binary search tree

    \note Search and update is carried out in a separate way. There has also been a combined
    approach, but this has been removed as part of GitLab Issue 181, as it is outperformed by the
    separate approach for large self contact problems!

    */
    virtual void SearchContact();

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
    //@}

    //! @name Print methods for debug and development purposes

    /*!
    \brief Plot the adjacency matrix

     */
    void PlotAdjacencyMatrix();

    /*!
    \brief Plot the dual graph

     */
    void PlotDualGraph(
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>> dualgraph);

    /*!
    \brief Plot root nodes and the self binary tree

     */
    void PlotRootsAndTree();
    //@}

    // don't want = operator and cctor
    SelfBinaryTree operator=(const SelfBinaryTree& old) = delete;
    SelfBinaryTree(const SelfBinaryTree& old) = delete;

    //! Local interface communicator
    const Epetra_Comm* lcomm_;

    //! All contact elements on surface (full map)
    Teuchos::RCP<Epetra_Map> elements_;

    //! Interface-specific parameter list
    const Teuchos::ParameterList& iparams_;

    /*!
    \brief Defining number of sample vectors

    \todo What are sample vectors?
    */
    int nvectors_;

    /*!
    \brief Defining sample vectors

    \todo What are sample vectors?
    */
    Epetra_SerialDenseMatrix samplevectors_;

    //! root treenodes
    std::vector<Teuchos::RCP<SelfBinaryTreeNode>> roots_;

    //! storage of all treenodes, sorted by layers
    std::vector<std::vector<Teuchos::RCP<SelfBinaryTreeNode>>> treenodes_;

    //! storage of all treenodes, sorted by layers
    std::map<int, std::vector<int>> contactpairs_;

    //! map of adjacent elements, sorted by global id (only needed in 3D)
    std::map<int, std::vector<Teuchos::RCP<SelfBinaryTreeNode>>> adjacencymatrix_;

    //! map of all leaf nodes, sorted by global id
    std::map<int, Teuchos::RCP<SelfBinaryTreeNode>> leafsmap_;

  };  // class SelfBinaryTree
}  // namespace CONTACT


#endif  // #ifndef DSELFCONTACT_BINARYTREE_H
