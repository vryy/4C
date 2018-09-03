/*!----------------------------------------------------------------------
\file mortar_binarytree.H

\brief A class for performing mortar search in 2D/3D based on binarytrees

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_BINARYTREE_H
#define MORTAR_BINARYTREE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"


// forward declarations
namespace DRT
{
  class Discretization;
}

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
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
  class BinaryTreeNode
  {
   public:
    /*!
    \brief constructor of a tree node

    \param type        type of BinaryTreeNode
    \param discret     interface discretization
    \param parent      points to parent tree node
    \param elelist     list of all elements in BinaryTreeNode
    \param dopnormals  reference to DOP normals
    \param kdop        refernece to no. of vertices
    \param dim         dimension of problem
    \param useauxpos   bool whether auxiliary position is used when computing dops
    \param layer       current layer of treenode
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
    \brief Calculate slabs of dop

    */
    void CalculateSlabsDop();

    /*!
    \brief Update slabs of current treenode in bottom up way

    */
    void UpdateSlabsBottomUp(double& eps);

    /*!
    \brief Initialize Tree

    */
    void InitializeTree(double& eps);

    /*!
    \brief Enlarge geometry of a Treenode by an offset, dependent on size

    */
    void EnlargeGeometry(double& eps);

    /*!
    \brief Divide a TreeNode into two child nodes

    */
    void DivideTreeNode();

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

    //! @name Access and modification methods

    /*!
    \brief Get communicator

    */
    virtual const Epetra_Comm& Comm() const;

    /*!
    \brief Get discretization of the interface

    */
    virtual DRT::Discretization& Discret() const { return idiscret_; }

    /*!
    \brief Return type of treenode

    */
    virtual BinaryTreeNodeType Type() const { return type_; }

    /*!
    \brief Set type of treenode

    */
    virtual void SetType(MORTAR::BinaryTreeNodeType type) { type_ = type; }

    /*!
    \brief Return pointer to right child

    */
    virtual Teuchos::RCP<BinaryTreeNode> Rightchild() const { return rightchild_; }

    /*!
    \brief Return pointer to left child

    */
    virtual Teuchos::RCP<BinaryTreeNode> Leftchild() const { return leftchild_; }

    /*!
    \brief Return pointer to parent

    */
    virtual Teuchos::RCP<BinaryTreeNode> Parent() const { return parent_; }

    /*!
    \brief Set children of a treenode

    */
    void SetChildren(
        Teuchos::RCP<BinaryTreeNode> leftchild, Teuchos::RCP<BinaryTreeNode> rightchild);

    /*!
    \brief Set parent of treenode

    */
    virtual void SetParent(Teuchos::RCP<BinaryTreeNode> parent) { parent_ = parent; }

    /*!
    \brief Return slabs of DOP

    */
    virtual Epetra_SerialDenseMatrix& Slabs() { return slabs_; }

    /*!
    \brief Return number of k-DOP vertices

    */
    virtual int Kdop() const { return kdop_; }

    /*!
    \brief Return layer of current TreeNode

    */
    virtual int Layer() const { return layer_; }

    /*!
    \brief Set layer of current TreeNode

    */
    virtual void SetLayer(int layer) { layer_ = layer; }

    /*!
    \brief Return reference to element list of TreeNode

    */
    virtual std::vector<int>& Elelist() { return elelist_; }

    /*!
    \brief Return dim of Problem

    */
    virtual const int& Dim() { return dim_; }

    /*!
    \brief Return different maps of tree nodes

    */
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Streenodesmap()
    {
      return streenodesmap_;
    }
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Mtreenodesmap()
    {
      return mtreenodesmap_;
    }
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Sleafsmap()
    {
      return sleafsmap_;
    }
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Mleafsmap()
    {
      return mleafsmap_;
    }

    //@}

    //! @name Visualization methods (GMSH)

    /*!
    \brief Print slabs of DOP to file for GMSH output
    \param filename     filename to which treenodes are plotted

    */
    void PrintDopsForGmsh(std::string filename);

    /*!
    \brief Plot a point in GMSH to given file

    */
    void PlotGmshPoint(std::string filename, double* position0, int nr);

    /*!
    \brief Plot a quadrangle in GMSH to given file

    */
    void PlotGmshQuadrangle(std::string filename, double* position0, double* position1,
        double* position2, double* position3);

    /*!
    \brief Plot a triangle in GMSH to given file

    */
    void PlotGmshTriangle(
        std::string filename, double* position0, double* position1, double* position2);

    //@}

   private:
    // don't want = operator and cctor
    BinaryTreeNode operator=(const BinaryTreeNode& old);
    BinaryTreeNode(const BinaryTreeNode& old);

    //! type of BinaryTreeNode
    MORTAR::BinaryTreeNodeType type_;

    //! interface discretization
    DRT::Discretization& idiscret_;

    //! pointer to the parent BinaryTreeNode
    Teuchos::RCP<BinaryTreeNode> parent_;

    //! pointer to the left child TreeNode
    Teuchos::RCP<BinaryTreeNode> leftchild_;

    //! pointer to the right child TreeNode
    Teuchos::RCP<BinaryTreeNode> rightchild_;

    //! geometry slabs of TreeNode, saved as Min|Max
    Epetra_SerialDenseMatrix slabs_;

    //! list with all elements of BinaryTreeNode
    std::vector<int> elelist_;

    //! reference to DOP normals
    const Epetra_SerialDenseMatrix& dopnormals_;

    //! reference to no. of vertices
    const int& kdop_;

    //! reference to dim. of problem
    const int& dim_;

    //! bool whether auxiliary position is used when computing dops
    const bool& useauxpos_;

    //! layer of TreeNode in Tree (0=rootnode!)
    int layer_;

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
  class BinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    Constructs an instance of this class.<br>

    \param discret (in):     The interface discretization
    \param selements (in):   All slave elements (column map)
    \param melements (in):   All master elements (fully overlapping map)
    \param dim (in):         The problem dimension

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
    void EvaluateSearch();

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
    virtual const Epetra_Comm& Comm() const;

    /*!
    \brief Get discretization of the interface

    */
    virtual DRT::Discretization& Discret() const { return idiscret_; }

    /*!
    \brief Get matrix of DOP normals

    */
    virtual const Epetra_SerialDenseMatrix& DopNormals() { return dopnormals_; }

    /*!
    \brief Get number of vertices of DOP

    */
    virtual const int& Kdop() { return kdop_; }

    /*!
    \brief Return dim of Problem

    */
    virtual const int& Dim() { return dim_; }

    /*!
    \brief Return reference to slave treenodesmap

    */
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Streenodesmap()
    {
      return streenodesmap_;
    }

    /*!
    \brief Return reference to master treenodesmap

    */
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& Mtreenodesmap()
    {
      return mtreenodesmap_;
    }

    /*!
    \brief Return reference to coupling treenodesmap

    */
    virtual std::vector<std::vector<Teuchos::RCP<BinaryTreeNode>>>& CouplingMap()
    {
      return couplingmap_;
    }

    /*!
    \brief Return pointer to sroot-treenode

    */
    virtual Teuchos::RCP<BinaryTreeNode>& Sroot() { return sroot_; }

    /*!
    \brief Return pointer to mroot-treenode

    */
    virtual Teuchos::RCP<BinaryTreeNode>& Mroot() { return mroot_; }

    /*!
    \brief Return minimal element length

    */
    virtual double& MinLengthEle() { return minlengthele_; }

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
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    void SetEnlarge();

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
      if (sroot_->Type() == 4) return;
      sroot_->CalculateSlabsDop();
      sroot_->EnlargeGeometry(enlarge_);
      mroot_->CalculateSlabsDop();
      mroot_->EnlargeGeometry(enlarge_);
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

    DRT::Discretization& idiscret_;       // mortar interface discretization
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
    double minlengthele_;                 // minimum of length of fe-elements
    const int dim_;                       // problem dimension (2D or 3D)
    double eps_;                          // epsilon for enlarging dops (of user)
    bool useauxpos_;  // bool whether auxiliary position is used when computing dops
    double enlarge_;  // to enlarge dops
    Epetra_SerialDenseMatrix dopnormals_;  // defining vertices of DOP
    int kdop_;                             // sets k for DOP (8 for 2D, 18 for 3D)

  };  // class BinaryTree
}  // namespace MORTAR

#endif  // #ifndef MORTAR_BINARYTREE_H
