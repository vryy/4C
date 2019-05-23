/*----------------------------------------------------------------------*/
/*!

\brief Basic discretization-related tools used in XFEM routines

\level 1

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_DISCRETIZATION_UTILS_H_
#define XFEM_DISCRETIZATION_UTILS_H_

#include <vector>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_fixedsizematrix.H"

class Epetra_Comm;
namespace DRT
{
  class Discretization;
  class DiscretizationInterface;
  class Condition;
  class Element;
  class Node;
}  // namespace DRT

namespace XFEM
{
  namespace UTILS
  {
    void PrintDiscretizationToStream(Teuchos::RCP<DRT::Discretization> dis,
        const std::string& disname, bool elements, bool elecol, bool nodes, bool nodecol,
        bool faces, bool facecol, std::ostream& s,
        std::map<int, LINALG::Matrix<3, 1>>* curr_pos = NULL);

    class XFEMDiscretizationBuilder
    {
     public:
      /// constructor
      XFEMDiscretizationBuilder(){/* should stay empty! */};

      void SetupXFEMDiscretization(const Teuchos::ParameterList& xgen_params,
          Teuchos::RCP<DRT::Discretization> dis, int numdof = 4) const;

      //! setup xfem discretization and embedded discretization
      void SetupXFEMDiscretization(const Teuchos::ParameterList& xgen_params,
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
      int SetupXFEMDiscretization(const Teuchos::ParameterList& xgen_params,
          Teuchos::RCP<DRT::DiscretizationInterface> src_dis,
          Teuchos::RCP<DRT::DiscretizationInterface> target_dis,
          const std::vector<DRT::Condition*>& boundary_conds) const;
      int SetupXFEMDiscretization(const Teuchos::ParameterList& xgen_params,
          Teuchos::RCP<DRT::Discretization> src_dis, Teuchos::RCP<DRT::Discretization> embedded_dis,
          const std::vector<DRT::Condition*>& boundary_conds) const;

     private:
      //! split a discretization into two by removing conditioned nodes
      //! in source and adding to target
      void SplitDiscretizationByCondition(
          Teuchos::RCP<DRT::Discretization> sourcedis,  //< initially contains all
          Teuchos::RCP<DRT::Discretization> targetdis,  //< initially empty
          std::vector<DRT::Condition*>& conditions,  //< conditioned nodes to be shifted to target
          const std::vector<std::string>& conditions_to_copy  //< conditions to copy to target
          ) const;

      /*! split the discretization by removing the given elements and nodes in
       *  the source discretization and adding them to the target discretization */
      void SplitDiscretization(Teuchos::RCP<DRT::Discretization> sourcedis,
          Teuchos::RCP<DRT::Discretization> targetdis, const std::map<int, DRT::Node*>& sourcenodes,
          const std::map<int, DRT::Node*>& sourcegnodes,
          const std::map<int, Teuchos::RCP<DRT::Element>>& sourceelements,
          const std::vector<std::string>& conditions_to_copy) const;

      //! re-partitioning of newly created discretizations (e.g. split by condition)
      void Redistribute(Teuchos::RCP<DRT::Discretization> dis, std::vector<int>& noderowvec,
          std::vector<int>& nodecolvec) const;

      /*! \brief Split a discretization into two parts by removing elements near boundary
       *         conditions
       *
       *  Split a volume source discretization into one part which is directly
       *  connected to the boundary condition face elements and the other. Currently
       *  only tested for the XCONTACT implementation.
       *
       *  \date 06/16
       *  \author hiermeier  */
      void SplitDiscretizationByBoundaryCondition(
          const Teuchos::RCP<DRT::Discretization>& sourcedis,
          const Teuchos::RCP<DRT::Discretization>& targetdis,
          const std::vector<DRT::Condition*>& boundary_conds,
          const std::vector<std::string>& conditions_to_copy) const;

      /** \brief remove conditions which are no longer part of the splitted
       *         partial discretizations, respectively
       *
       *  \author  hiermeier \date 10/16 */
      Teuchos::RCP<DRT::Condition> SplitCondition(const DRT::Condition* src_cond,
          const std::vector<int>& nodecolvec, const Epetra_Comm& comm) const;
    };  // class XFEMDiscretizationBuilder
  }     // namespace UTILS
}  // namespace XFEM


#endif /* XFEM_DISCRETIZATION_UTILS_H_ */
