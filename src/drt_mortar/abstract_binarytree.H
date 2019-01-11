/*!----------------------------------------------------------------------
\file abstract_binarytree.H

\brief An abstract interface for binary trees

\level 1

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#ifndef ABSTRACT_BINARYTREE_H
#define ABSTRACT_BINARYTREE_H

namespace MORTAR
{
  /*!
  \brief An abstract interface for binary trees

  \author schmidt (schmidt@lnm.mw.tum.de)
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
    virtual ~AbstractBinaryTree(){};

    //! @name Query methods

    /*!
    \brief Evaluate search tree to get corresponding master elements for the slave elements

    */
    virtual void EvaluateSearch() = 0;

    /*!
    \brief initialize the binary tree

    */
    virtual void Init() = 0;
  };  // class AbstractBinaryTree
}  // namespace MORTAR

#endif  // #ifndef ABSTRACT_BINARYTREE_H
