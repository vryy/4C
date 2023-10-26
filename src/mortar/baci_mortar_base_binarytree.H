/*----------------------------------------------------------------------*/
/*! \file
\brief A base class for binary trees and binary tree nodes providing common functionality

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MORTAR_BASE_BINARYTREE_H
#define BACI_MORTAR_BASE_BINARYTREE_H

#include "baci_mortar_abstract_binarytree.H"

#include <Teuchos_RCP.hpp>

#include <vector>

// forward declarations
namespace DRT
{
  class Discretization;
}
namespace CONTACT
{
  class SelfBinaryTree;
  class SelfDualEdge;
  class UnbiasedSelfBinaryTree;
}  // namespace CONTACT

namespace MORTAR
{
  /*!
  \brief A base class for binary tree nodes

  */
  class BaseBinaryTreeNode : public AbstractBinaryTreeNode
  {
    // these classes need access to protected methods of this class and are therefore defined friend
    // classes here
    friend class BinaryTree;
    friend class CONTACT::SelfBinaryTree;
    friend class CONTACT::SelfDualEdge;
    friend class CONTACT::UnbiasedSelfBinaryTree;

   public:
    /*!
    \brief Standard constructor of a base binary tree node

    \param [in] discret:     interface discretization
    \param [in] elelist:     list of all elements in BaseBinaryTreeNode
    \param [in] dopnormals:  reference to DOP normals
    \param [in] kdop:        reference to no. of vertices
    \param [in] dim:         dimension of the problem
    \param [in] useauxpos:   bool whether auxiliary position is used when computing dops
    \param [in] layer:       current layer of tree node

    */
    BaseBinaryTreeNode(DRT::Discretization& discret, std::vector<int> elelist,
        const CORE::LINALG::SerialDenseMatrix& dopnormals, const int& kdop, const int& dim,
        const bool& useauxpos, const int layer);

    /*!
    \brief Destructor

    */
    ~BaseBinaryTreeNode() override{};

    //! @name Evaluation methods

    /*!
    \brief Calculate slabs of dop

    */
    void CalculateSlabsDop() override;

    /*!
    \brief Update slabs of current tree node in bottom up way

    */
    void UpdateSlabsBottomUp(double& eps) override = 0;

    /*!
    \brief Enlarge geometry of a tree node by an offset, dependent on size

    */
    void EnlargeGeometry(double& enlarge) override;
    //@}

    //! @name Print methods
    // please note: these methods do not get called and are therefore untested. However, they might
    // be of interest for debug and development purposes
    /*!
    \brief Print type of tree node to std::cout

    */
    virtual void PrintType() = 0;

    /*!
    \brief Print slabs to std::cout

    */
    void PrintSlabs();
    //@}

    //! @name Visualization methods (GMSH)
    // please note: these methods only get called if "MORTARGMSHTN" of mortar_defines.H is
    // activated. However, this might be interesting functionality for debug and development
    // purposes

    /*!
    \brief Print slabs of DOP to file for GMSH output
    \param filename     filename to which tree nodes are plotted

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

   protected:
    //! @name Access and modification methods
    /*!
    \brief Return dim of Problem

    */
    const int& Dim() const { return dim_; }

    /*!
    \brief Get discretization of the interface

    */
    DRT::Discretization& Discret() const { return idiscret_; }

    /*!
    \brief Return pointer to normals of DOP

    */
    const CORE::LINALG::SerialDenseMatrix& Dopnormals() const { return dopnormals_; }

    /*!
    \brief Return pointer to element list of tree node

    */
    std::vector<int> Elelist() const { return elelist_; }

    /*!
    \brief Return no. of vertices

    */
    const int& Kdop() const { return kdop_; }

    /*!
    \brief Return layer of current tree node

    */
    int Layer() const { return layer_; }

    /*!
    \brief Set layer of current tree node

    */
    void SetLayer(int layer) { layer_ = layer; }

    /*!
    \brief Return pointer to slabs of DOP

    */
    CORE::LINALG::SerialDenseMatrix& Slabs() { return slabs_; }

    /*!
    \brief Return bool indicating whether auxiliary position is used when computing dops

    */
    const bool& UseAuxPos() const { return useauxpos_; }
    //@}

   private:
    //! dimension of the problem
    const int dim_;
    //! reference to DOP normals
    const CORE::LINALG::SerialDenseMatrix& dopnormals_;
    //! list containing the gids of all elements of the tree node
    std::vector<int> elelist_;
    //! interface discretization
    DRT::Discretization& idiscret_;
    //! number of vertices
    const int kdop_;
    //! layer of tree node in tree (0=root node!)
    int layer_;
    //! geometry slabs of tree node, saved as Min|Max
    CORE::LINALG::SerialDenseMatrix slabs_;
    //! auxiliary position is used when computing dops
    const bool useauxpos_;
  };  // class BaseBinaryTreeNode

  /*!
  \brief A base class for binary trees

  */
  class BaseBinaryTree : public AbstractBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    \param [in] discret: interface discretization
    \param [in] dim:     dimension of the problem
    \param [in] eps:     factor used to enlarge dops
    */
    BaseBinaryTree(DRT::Discretization& discret, int dim, double eps);

    /*!
    \brief Destructor

    */
    ~BaseBinaryTree() override{};

    //! @name Evaluation methods

    /*!
    \brief Evaluate search tree

    */
    void EvaluateSearch() override = 0;

    /*!
    \brief Initialize the base binary tree

    */
    void Init() override;

    /*!
    \brief Calculate minimal element length / inflation factor "enlarge"

    */
    virtual void SetEnlarge() = 0;
    //@}

   protected:
    //! @name Access and modification methods
    /*!
    \brief Return dim of the problem

    */
    const int& Dim() const { return dim_; }

    /*!
    \brief Get discretization of the interface

    */
    DRT::Discretization& Discret() const { return idiscret_; }

    /*!
    \brief Get matrix of DOP normals

    */
    const CORE::LINALG::SerialDenseMatrix& DopNormals() const { return dopnormals_; }

    /*!
    \brief Return factor "enlarge" to enlarge dops

    */
    double& Enlarge() { return enlarge_; }

    /*!
    \brief Return factor "eps" to set "enlarge"

    */
    const double& Eps() const { return eps_; }

    /*!
    \brief Get number of vertices of DOP

    */
    const int& Kdop() const { return kdop_; }
    //@}

   private:
    /*!
    \brief Initialize internal variables

     */
    virtual void InitInternalVariables() = 0;

    //! interface discretization
    DRT::Discretization& idiscret_;
    //! problem dimension (2D or 3D)
    const int dim_;
    //! normals of DOP
    CORE::LINALG::SerialDenseMatrix dopnormals_;
    //! needed to enlarge dops
    double enlarge_;
    //! epsilon for enlarging dops (of user)
    const double eps_;
    //! set k for DOP (8 for 2D, 18 for 3D)
    int kdop_;
  };  // class BaseBinaryTree
}  // namespace MORTAR

#endif  // MORTAR_BASE_BINARYTREE_H
