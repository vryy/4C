/*----------------------------------------------------------------------*/
/*!
\file base_binarytree.H

\brief A base class for binary trees and binary tree nodes providing common functionality

\level 1

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

#ifndef BASE_BINARYTREE_H
#define BASE_BINARYTREE_H

#include "abstract_binarytree.H"
#include "Teuchos_RCP.hpp"

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

  \author schmidt (schmidt@lnm.mw.tum.de)
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
        const Epetra_SerialDenseMatrix& dopnormals, const int& kdop, const int& dim,
        const bool& useauxpos, const int layer);

    /*!
    \brief Destructor

    */
    virtual ~BaseBinaryTreeNode(){};

    //! @name Evaluation methods

    /*!
    \brief Calculate slabs of dop

    */
    virtual void CalculateSlabsDop();

    /*!
    \brief Update slabs of current tree node in bottom up way

    */
    virtual void UpdateSlabsBottomUp(double& eps) = 0;

    /*!
    \brief Enlarge geometry of a tree node by an offset, dependent on size

    */
    virtual void EnlargeGeometry(double& eps);
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
    const Epetra_SerialDenseMatrix& Dopnormals() const { return dopnormals_; }

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
    Epetra_SerialDenseMatrix& Slabs() { return slabs_; }

    /*!
    \brief Return bool indicating whether auxiliary position is used when computing dops

    */
    const bool& UseAuxPos() const { return useauxpos_; }
    //@}

   private:
    const int& dim_;                              // reference to dim. of problem
    const Epetra_SerialDenseMatrix& dopnormals_;  // reference to DOP normals
    std::vector<int> elelist_;                    // list with all elements of tree node
    DRT::Discretization& idiscret_;               // interface discretization
    const int& kdop_;                             // reference to no. of vertices
    int layer_;                                   // layer of tree node in tree (0=root node!)
    Epetra_SerialDenseMatrix slabs_;              // geometry slabs of TreeNode, saved as Min|Max
    const bool& useauxpos_;  // bool whether auxiliary position is used when computing dops

  };  // class BaseBinaryTreeNode

  /*!
  \brief A base class for binary trees

  \author schmidt (schmidt@lnm.mw.tum.de)
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
    virtual ~BaseBinaryTree(){};

    //! @name Evaluation methods

    /*!
    \brief Evaluate search tree

    */
    virtual void EvaluateSearch() = 0;

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
    const Epetra_SerialDenseMatrix& DopNormals() const { return dopnormals_; }

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

    // internal members
    DRT::Discretization& idiscret_;        // interface discretization
    const int dim_;                        // problem dimension (2D or 3D)
    Epetra_SerialDenseMatrix dopnormals_;  // defining vertices of DOP
    double enlarge_;                       // needed to enlarge dops
    const double eps_;                     // epsilon for enlarging dops (of user)
    int kdop_;                             // set k for DOP (8 for 2D, 18 for 3D)

  };  // class BaseBinaryTree
}  // namespace MORTAR

#endif  // #ifndef BASE_BINARYTREE_H
