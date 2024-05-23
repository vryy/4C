/*-----------------------------------------------------------------------*/
/*! \file
\brief Search tree for unbiased self-contact problems

\level 2

*/
/*-----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_SELFCONTACT_BINARYTREE_UNBIASED_HPP
#define FOUR_C_CONTACT_SELFCONTACT_BINARYTREE_UNBIASED_HPP

#include "4C_config.hpp"

#include "4C_contact_selfcontact_binarytree.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
  \brief A class for performing unbiased self contact search in 2D / 3D based on a binary search
  tree and processor specific dual graphs to enhance parallel search

  */

  class UnbiasedSelfBinaryTree : public SelfBinaryTree
  {
   public:
    /*!
    \brief Standard constructor

    \param discret (in):    The contact interface discretization
    \param iparams (in):    interface specific parameter list
    \param elements (in):   All elements on self contact interface (fully overlapping map)
    \param dim (in):        The problem dimension

    */
    UnbiasedSelfBinaryTree(DRT::Discretization& discret, const Teuchos::ParameterList& iparams,
        Teuchos::RCP<Epetra_Map> elements, int dim, double eps);


    /*!
    \brief Initialize the unbiased self binary tree

    */
    void Init() final;

   private:
    //! @name Evaluation methods

    /*!
    \brief Decide whether tree nodes need to be added to contact pairs

    \param [in]  treenode1:  first tree node
    \param [in]  treenode2:  second tree node

    */
    void add_tree_nodes_to_contact_pairs(Teuchos::RCP<SelfBinaryTreeNode> treenode1,
        Teuchos::RCP<SelfBinaryTreeNode> treenode2) final;

    /*!
    \brief Calculate the processor specific dual graph

    \param [out] dualgraph:  construction of binary tree is based on this data
    \param [in]  elelist:    list (gids) of all contact elements of the surface
    \param [in]  p:          number of current (not necessarily calling) processor

    */
    void calculate_proc_specific_dual_graph(
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>* dualGraph,
        const std::vector<int>& elelist, const int p);

    /*!
    \brief Use list of possible contact pairs to define the search elements

    */
    void define_search_elements();

    /*!
    \brief Get the mortar element specific number of first order nodes

    \param [in]      contractedEdge:  dual edge that is contracted
    \param [in,out]  contractedNode:  node that consists of both nodes of contracted edge

    */
    void GetContractedNode(Teuchos::RCP<SelfDualEdge>& contractedEdge,
        Teuchos::RCP<SelfBinaryTreeNode>& contractedNode) final;

    /*!
    \brief Initialize Tree in a bottom up way based on dual graph

    \param [in] procdualgraph:  processor specific dual graph, construction of binary tree is based
    on this data

    */
    void initialize_tree_bottom_up(std::map<int,
        std::map<Teuchos::RCP<SelfDualEdge>, std::vector<Teuchos::RCP<SelfDualEdge>>>>*
            procdualGraph);

    /*!
    \brief Checks roughly whether self contact of two elements shall be evaluated (3D)

    This method checks if the normal at the slave element center and the vector connecting this
    slave element center and the center of the element it is projected to (master element) point in
    the same direction. All these vectors are evaluated in the reference configuration to check the
    initial state. Only if both vectors point in the same direction integration shall be performed.

    */
    bool RoughCheckRefConfig(int ele1gid, int ele2gid);

    /*!
    \brief Evaluate unbiased binary search tree

    */
    void SearchContact() final;

    /*!
    \brief Communicate the Search Elements to all processors

    */
    void communicate_search_elements_all_procs();
    //@}

    // don't want = operator and cctor
    UnbiasedSelfBinaryTree operator=(const SelfBinaryTree& old) = delete;
    UnbiasedSelfBinaryTree(const SelfBinaryTree& old) = delete;

    //! use two half pass approach
    const bool two_half_pass_;
    //! perform reference configuration check for non-smooth self contact
    const bool check_nonsmooth_selfcontactsurface_;
    //! the contact pairs are communicated to all processors
    const bool searchele_all_proc_;
  };  // class UnbiasedSelfBinaryTree
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
