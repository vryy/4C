/*!----------------------------------------------------------------------
\file mortar_binarytree.H

\brief A class for performing mortar search in 2D/3D based on binarytrees

\level 1

\maintainer Alexander Popp

*----------------------------------------------------------------------*/
#ifndef MORTAR_BINARYTREE_H
#define MORTAR_BINARYTREE_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "base_binarytree.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  // forward declarations

  //! @name Enums and Friends

  /// Type of binary tree node
  enum BinaryTreeNodeType
  {
    SLAVE_INNER,        ///< indicates a slave inner node (has children)
    SLAVE_LEAF,         ///< indicates a slave leaf node (no further children)
    MASTER_INNER,       ///< indicates a master inner node (has children)
    MASTER_LEAF,        ///< indicates a master leaf node (no further children)
    NOSLAVE_ELEMENTS,   ///< indicates that there are no slave elements on this (root) treenode
    NOMASTER_ELEMENTS,  ///< indicates that there are no master elements on this (root) treenode
    UNDEFINED           ///< indicates an undefined tree node
  };

  //@}

  /*!
  \brief A class representing one tree node of the binary search tree

  Refer also to the Semesterarbeit of Thomas Eberl, 2009

  \author popp (popp@lnm.mw.tum.de)
  */
  class BinaryTreeNode : public BaseBinaryTreeNode
  {
   public:
    /*!
    \brief constructor of a tree node

    \param type        type of BinaryTreeNode
    \param discret     interface discretization
    \param parent      points to parent tree node
    \param elelist     list of all elements in BinaryTreeNode
    \param dopnormals  reference to DOP normals
    \param kdop        reference to no. of vertices
    \param dim         dimension of problem
    \param useauxpos   bool whether auxiliary position is used when computing dops
    \param layer       current layer of tree node
    \param ...map      references to maps

    */
    BinaryTreeNode(BinaryTreeNodeType type, DRT::Discretization& discret,
        Teuchos::RCP<BinaryTreeNode> parent, std::vector<int> elelist,
        const Epetra_SerialDenseMatrix& dopnormals, const int& kdop, const int& dim,
        const bool& useauxpos, const int layer,
        std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& streenodesmap,
        std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& mtreenodesmap,
        std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& sleafsmap,
        std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& mleafsmap);

    /*!
    \brief Destructor

    */
    virtual ~BinaryTreeNode(){};

    //! @name Evaluation methods

    /*!
    \brief Update slabs of current treenode in bottom up way

    */
    void UpdateSlabsBottomUp(double& eps);

    /*!
    \brief Initialize Tree

    */
    void InitializeTree(double& eps);

    /*!
    \brief Divide a TreeNode into two child nodes

    */
    void DivideTreeNode();

    /*!
    \brief Print type of tree node to std::cout

    */
    virtual void PrintType() final;

    /*!
    \brief Set slabs of current treenode

    */
    void SetSlabs(Epetra_SerialDenseMatrix& newslabs);

    //@}

    //! @name Access and modification methods

    /*!
    \brief Get communicator

    */
    const Epetra_Comm& Comm() const;

    /*!
    \brief Return type of treenode

    */
    BinaryTreeNodeType Type() const { return type_; }

    /*!
    \brief Set type of treenode

    */
    void SetType(MORTAR::BinaryTreeNodeType type) { type_ = type; }

    /*!
    \brief Return pointer to right child

    */
    Teuchos::RCP<BinaryTreeNode> Rightchild() const { return rightchild_; }

    /*!
    \brief Return pointer to left child

    */
    Teuchos::RCP<BinaryTreeNode> Leftchild() const { return leftchild_; }

    /*!
    \brief Return pointer to parent

    */
    Teuchos::RCP<BinaryTreeNode> Parent() const { return parent_; }

    /*!
    \brief Set children of a treenode

    */
    void SetChildren(
        Teuchos::RCP<BinaryTreeNode> leftchild, Teuchos::RCP<BinaryTreeNode> rightchild);

    /*!
    \brief Set parent of treenode

    */
    void SetParent(Teuchos::RCP<BinaryTreeNode> parent) { parent_ = parent; }

    /*!
    \brief Return different maps of tree nodes

    */
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Streenodesmap()
    {
      return streenodesmap_;
    }
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Mtreenodesmap()
    {
      return mtreenodesmap_;
    }
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Sleafsmap() { return sleafsmap_; }
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Mleafsmap() { return mleafsmap_; }

    //@}

   private:
    // don't want = operator and cctor
    BinaryTreeNode operator=(const BinaryTreeNode& old);
    BinaryTreeNode(const BinaryTreeNode& old);

    //! type of BinaryTreeNode
    MORTAR::BinaryTreeNodeType type_;

    // the pointers to the parent as well as to the left and right child are not moved to the
    // BaseBinaryTreeNode as this would require a lot of dynamic casting and thereby complicating
    // the readability of the code
    //! pointer to the parent BinaryTreeNode
    Teuchos::RCP<BinaryTreeNode> parent_;

    //! pointer to the left child TreeNode
    Teuchos::RCP<BinaryTreeNode> leftchild_;

    //! pointer to the right child TreeNode
    Teuchos::RCP<BinaryTreeNode> rightchild_;

    //! reference to map of all slave treenodes, sorted by layer
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& streenodesmap_;

    //! reference to map of all master treenodes, sorted by layer
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& mtreenodesmap_;

    //! reference to map of all slave leaf treenodes
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& sleafsmap_;

    //! reference to map of all master leaf treenodes
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& mleafsmap_;

  };  // class BinaryTreeNode


  /*!
  \brief A class for performing search in 2D/3D based on binary trees

  Refer also to the Semesterarbeit of Thomas Eberl, 2009

  \author popp (popp@lnm.mw.tum.de)
  */
  class BinaryTree : public BaseBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    Constructs an instance of this class.<br>

    \param discret (in):     The interface discretization
    \param selements (in):   All slave elements (column map)
    \param melements (in):   All master elements (fully overlapping map)
    \param dim (in):         The problem dimension
    \param eps (in):
    \param useauxpos (in):

    */
    BinaryTree(DRT::Discretization& discret, Teuchos::RCP<Epetra_Map> selements,
        Teuchos::RCP<Epetra_Map> melements, int dim, double eps, bool useauxpos);

    /*!
    \brief Destructor

    */
    virtual ~BinaryTree(){};

    //! @name Query methods

    /*!
    \brief Evaluate search tree to get corresponding master elements for the slave elements

    */
    void EvaluateSearch() final;

    /*!
    \brief Initialize the binary tree

    */
    void Init() final;

   private:
    /*!
    \brief clear found search elements

    */
    void InitSearchElements();

    /*!
    \brief Print full tree

    */
    void PrintTree(Teuchos::RCP<BinaryTreeNode> treenode);

    /*!
    \brief Print full tree out of map of treenodes

    */
    void PrintTreeOfMap(std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& treenodesmap);

    //@}

    //! @name Access methods

    /*!
    \brief Get communicator

    */
    const Epetra_Comm& Comm() const;

    /*!
    \brief Return reference to slave treenodesmap

    */
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Streenodesmap()
    {
      return streenodesmap_;
    }

    /*!
    \brief Return reference to master treenodesmap

    */
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Mtreenodesmap()
    {
      return mtreenodesmap_;
    }

    /*!
    \brief Return reference to coupling treenodesmap

    */
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& CouplingMap() { return couplingmap_; }

    /*!
    \brief Return pointer to sroot-treenode

    */
    Teuchos::RCP<BinaryTreeNode>& Sroot() { return sroot_; }

    /*!
    \brief Return pointer to mroot-treenode

    */
    Teuchos::RCP<BinaryTreeNode>& Mroot() { return mroot_; }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Initialize internal variables

     */
    void InitInternalVariables() final;

    /*!
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    void SetEnlarge() override;

    /*!
    \brief Update master and slave tree in a topdown way

    */
    void UpdateTreeTopDown()
    {
      EvaluateUpdateTreeTopDown(sroot_);
      EvaluateUpdateTreeTopDown(mroot_);
      return;
    }

    /*!
    \brief Evaluate update of master and slave tree in a topdown way

    */
    void EvaluateUpdateTreeTopDown(Teuchos::RCP<BinaryTreeNode> treenode);

    /*!
    \brief Updates master and slave tree in a bottom up way

    */
    void UpdateTreeBottomUp()
    {
      EvaluateUpdateTreeBottomUp(streenodesmap_);
      EvaluateUpdateTreeBottomUp(mtreenodesmap_);
      return;
    }

    /*!
    \brief Evaluate update of master and slave tree in a bottom up way

    */
    void EvaluateUpdateTreeBottomUp(
        std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& treenodesmap);

    /*!
    \brief Evaluate Binary search tree for separate update and search

    */
    void SearchSeparate()
    {
      EvaluateSearchSeparate(sroot_, mroot_);
      return;
    }

    /*!
    \brief Evaluate Binary search tree for separate update and search

    */
    void EvaluateSearchSeparate(
        Teuchos::RCP<BinaryTreeNode> streenode, Teuchos::RCP<BinaryTreeNode> mtreenode);

    /*!
    \brief Evaluate Binary search tree for combined search and update

    */
    void SearchCombined()
    {  // if proc has no slave elements --> return;
      if (sroot_->Type() == NOSLAVE_ELEMENTS) return;
      sroot_->CalculateSlabsDop();
      sroot_->EnlargeGeometry(Enlarge());
      mroot_->CalculateSlabsDop();
      mroot_->EnlargeGeometry(Enlarge());
      EvaluateSearchCombined(sroot_, mroot_);
      return;
    }

    /*!
    \brief Evaluate Binary search tree for combined search and update

    */
    void EvaluateSearchCombined(
        Teuchos::RCP<BinaryTreeNode> streenode, Teuchos::RCP<BinaryTreeNode> mtreenode);


    // don't want = operator and cctor
    BinaryTree operator=(const BinaryTree& old);
    BinaryTree(const BinaryTree& old);

    Teuchos::RCP<Epetra_Map> selements_;  // all slave elements on surface (column map)
    Teuchos::RCP<Epetra_Map> melements_;  // all master elements on surface (full map)
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>
        streenodesmap_;  // map of all slave treenodes, sorted by layers
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>
        mtreenodesmap_;  // map of all master treenodes, sorted by layers
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>
        couplingmap_;  // map of all treenodes, that possibly couple, st/mt
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>
        sleafsmap_;  // map of all slave leaf treenodes, [0]=leftchild,[1]=rightchild
    std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>
        mleafsmap_;  // map of all master leaf treenodes, [0]=leftchild,[1]=rightlchild
    Teuchos::RCP<BinaryTreeNode> sroot_;  // slave root treenode
    Teuchos::RCP<BinaryTreeNode> mroot_;  // master root treenode
    bool useauxpos_;  // bool whether auxiliary position is used when computing dops

  };  // class BinaryTree
}  // namespace MORTAR

#endif  // #ifndef MORTAR_BINARYTREE_H
