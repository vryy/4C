/*-----------------------------------------------------------------------*/
/*! \file
\brief An abstract interface for binary trees

\level 1

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_ABSTRACT_BINARYTREE_HPP
#define FOUR_C_MORTAR_ABSTRACT_BINARYTREE_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MORTAR
{
  /*!
  \brief An abstract interface for binary trees

  */
  class AbstractBinaryTreeNode
  {
   public:
    /*!
    \brief Standard constructor

    */
    AbstractBinaryTreeNode() { return; };

    /*!
    \brief Destructor

    */
    virtual ~AbstractBinaryTreeNode() = default;

    //! @name Evaluation methods

    /*!
    \brief Calculate slabs of dop

    */
    virtual void CalculateSlabsDop() = 0;

    /*!
    \brief Update slabs of current treenode in bottom up way

    */
    virtual void UpdateSlabsBottomUp(double& eps) = 0;

    /*!
    \brief Enlarge geometry of a Treenode by an offset, dependent on size

    */
    virtual void EnlargeGeometry(double& eps) = 0;
  };  // class AbstractBinaryTreeNode


  /*!
  \brief An abstract interface for binary trees

  */
  class AbstractBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    */
    AbstractBinaryTree() { return; };

    /*!
    \brief Destructor

    */
    virtual ~AbstractBinaryTree() = default;

    //! @name Query methods

    /*!
    \brief Evaluate search tree to get corresponding master elements for the slave elements

    */
    virtual void evaluate_search() = 0;

    /*!
    \brief initialize the binary tree

    */
    virtual void Init() = 0;
  };  // class AbstractBinaryTree
}  // namespace MORTAR

FOUR_C_NAMESPACE_CLOSE

#endif
