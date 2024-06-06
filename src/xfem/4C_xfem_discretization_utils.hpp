/*----------------------------------------------------------------------*/
/*! \file

\brief Basic discretization-related tools used in XFEM routines

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_XFEM_DISCRETIZATION_UTILS_HPP
#define FOUR_C_XFEM_DISCRETIZATION_UTILS_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_lib_discret_faces.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <string>
#include <vector>


FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace CORE::Nodes
{
  class Node;
}

namespace CORE::Elements
{
  class Element;
}

namespace XFEM
{
  namespace UTILS
  {
    void PrintDiscretizationToStream(Teuchos::RCP<DRT::Discretization> dis,
        const std::string& disname, bool elements, bool elecol, bool nodes, bool nodecol,
        bool faces, bool facecol, std::ostream& s,
        std::map<int, CORE::LINALG::Matrix<3, 1>>* curr_pos = nullptr);

    class XFEMDiscretizationBuilder
    {
     public:
      /// constructor
      XFEMDiscretizationBuilder(){/* should stay empty! */};

      void setup_xfem_discretization(const Teuchos::ParameterList& xgen_params,
          Teuchos::RCP<DRT::Discretization> dis, int numdof = 4) const;

      //! setup xfem discretization and embedded discretization
      void setup_xfem_discretization(const Teuchos::ParameterList& xgen_params,
          Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<DRT::Discretization> embedded_dis,
          const std::string& embedded_cond_name, int numdof = 4) const;

      /*! \brief Setup xfem discretization and embedded discretization
       *  by using a boundary condition vector
       *
       *  Split the given discretization into a XFEM discretization, which
       *  is located next to the (enriched) conditioned interfaces and a standard
       *  discretization
       *
       *           ___boundary cond___    _____ enriched and conditioned boundary
       *          /                   \  /      interface nodes (o)
       *                                /
       *          o---o---o---o---o---o
       *         /   /   /   /   /   /|       enriched element row (xFem discret.)
       *        o---o---o---o---o---o +   <== (enriched (o) and std. nodes (+))
       *        | 0 | 1 | 2 | 3 | 4 |/|
       *        +---+---+---+---+---+ +       standard element row (std. discret.)
       *        | 5 | 6 | 7 | 8 | 9 |/    <== (only std. nodes (+))
       *        +---+---+---+---+---+
       *
       *                 __
       *                |  |
       *               _|  |_
       *               \    /
       *                \  /
       *                 \/
       *
       *  We get one new cut xFem discretization, which is connected to the
       *  conditioned boundary interface (o)
       *
       *          o---o---o---o---o---o
       *         /   /   /   /   /   /|
       *        o---o---o---o---o---o +
       *        | 0 | 1 | 2 | 3 | 4 |/   <== xstruct_dis_ptr
       *        +---+---+---+---+---+
       *
       *  and the remaining standard discretization (+)
       *
       *          +---+---+---+---+---+
       *         /   /   /   /   /   /|
       *        +---+---+---+---+---+ +
       *        | 5 | 6 | 7 | 8 | 9 |/   <== struct_dis_ptr_
       *        +---+---+---+---+---+
       *
       *  The two discretizations share the same node ID's at the coupling interface,
       *  but differ in the global degrees of freedom ID's!
       *
       *  \author hiermeier
       *  \date 06/16 */
      int setup_xfem_discretization(const Teuchos::ParameterList& xgen_params,
          Teuchos::RCP<DRT::Discretization> src_dis, Teuchos::RCP<DRT::Discretization> target_dis,
          const std::vector<CORE::Conditions::Condition*>& boundary_conds) const;

     private:
      //! split a discretization into two by removing conditioned nodes
      //! in source and adding to target
      void split_discretization_by_condition(
          Teuchos::RCP<DRT::Discretization> sourcedis,  //< initially contains all
          Teuchos::RCP<DRT::Discretization> targetdis,  //< initially empty
          std::vector<CORE::Conditions::Condition*>&
              conditions,  //< conditioned nodes to be shifted to target
          const std::vector<std::string>& conditions_to_copy  //< conditions to copy to target
      ) const;

      /*! split the discretization by removing the given elements and nodes in
       *  the source discretization and adding them to the target discretization */
      void split_discretization(Teuchos::RCP<DRT::Discretization> sourcedis,
          Teuchos::RCP<DRT::Discretization> targetdis,
          const std::map<int, CORE::Nodes::Node*>& sourcenodes,
          const std::map<int, CORE::Nodes::Node*>& sourcegnodes,
          const std::map<int, Teuchos::RCP<CORE::Elements::Element>>& sourceelements,
          const std::vector<std::string>& conditions_to_copy) const;

      //! re-partitioning of newly created discretizations (e.g. split by condition)
      void redistribute(Teuchos::RCP<DRT::Discretization> dis, std::vector<int>& noderowvec,
          std::vector<int>& nodecolvec) const;

      /*! \brief Split a discretization into two parts by removing elements near boundary
       *         conditions
       *
       *  Split a volume source discretization into one part which is directly
       *  connected to the boundary condition face elements and the other. Currently
       *  not tested.
       *
       *  \date 06/16
       *  \author hiermeier  */
      void split_discretization_by_boundary_condition(
          const Teuchos::RCP<DRT::Discretization>& sourcedis,
          const Teuchos::RCP<DRT::Discretization>& targetdis,
          const std::vector<CORE::Conditions::Condition*>& boundary_conds,
          const std::vector<std::string>& conditions_to_copy) const;

      /** \brief remove conditions which are no longer part of the splitted
       *         partial discretizations, respectively
       *
       *  \author  hiermeier \date 10/16 */
      Teuchos::RCP<CORE::Conditions::Condition> split_condition(
          const CORE::Conditions::Condition* src_cond, const std::vector<int>& nodecolvec,
          const Epetra_Comm& comm) const;
    };  // class XFEMDiscretizationBuilder
  }     // namespace UTILS

  class DiscretizationXWall : public DRT::DiscretizationFaces
  {
   public:
    /*!
    \brief Standard Constructor

    \param name: name of this discretization
    \param comm: Epetra comm object associated with this discretization
    \param n_dim: number of space dimensions of this discretization
    */
    DiscretizationXWall(const std::string name, Teuchos::RCP<Epetra_Comm> comm, unsigned int n_dim);



    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to assign_degrees_of_freedom()))

    Additional input nodal dof set: If the node contains more than one set of dofs, which can be
    evaluated, the number of the set needs to be given. Currently only the case for XFEM.

    \param dof (in)         : vector of dof gids (to be filled)
    \param nds (in)         : number of dofset
    \param nodaldofset (in) : number of nodal dofset
    \param node (in)        : the node
    \param element (in)     : the element (optionally)
    */
    void Dof(std::vector<int>& dof, const CORE::Nodes::Node* node, unsigned nds,
        unsigned nodaldofset, const CORE::Elements::Element* element = nullptr) const override
    {
      if (nds > 1) FOUR_C_THROW("xwall discretization can only handle one dofset at the moment");

      FOUR_C_ASSERT(nds < dofsets_.size(), "undefined dof set");
      FOUR_C_ASSERT(havedof_, "no dofs assigned");

      std::vector<int> totaldof;
      dofsets_[nds]->Dof(totaldof, node, nodaldofset);

      if (element == nullptr && element->Shape() == CORE::FE::CellType::hex8)
        FOUR_C_THROW("element required for location vector of hex8 element");

      int size;
      if (element != nullptr)
        size = std::min((int)totaldof.size(), element->NumDofPerNode(*node));
      else
        size = (int)totaldof.size();
      // only take the first dofs that have a meaning for all elements at this node
      for (int i = 0; i < size; i++) dof.push_back(totaldof.at(i));


      return;
    }
  };  // class DiscretizationXWall
}  // namespace XFEM


FOUR_C_NAMESPACE_CLOSE

#endif
