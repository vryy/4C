/*----------------------------------------------------------------------*/
/*! \file

\brief utils functions for conditions

\level 1


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_LIB_CONDITION_UTILS_HPP
#define FOUR_C_LIB_CONDITION_UTILS_HPP

#include "baci_config.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class MapExtractor;
}

namespace DRT
{
  class Condition;
  class Discretization;
  class Element;
  class Node;

  namespace UTILS
  {
    // forward declaration
    class ConditionSelector;

    /// std unary function version of Epetra_Map::MyGID()
    struct MyGID : public std::unary_function<int, bool>
    {
      const Epetra_Map* emap_;
      MyGID(const Epetra_Map* emap) : emap_(emap) {}
      bool operator()(int gid) const { return emap_->MyGID(gid); }
    };

    /// find all local nodes from discretization marked with condition
    /*!
      Loop all conditions of the given discretization, find the ones with the
      specified name and collect the locally owned node ids in the supplied
      set. The nodes vector is unique and ordered on output.

      \param dis : (in) Discretization
      \param condname : (in) name of condition in question
      \param nodes : (out) empty set on input, filled with nodal gids of local nodes
    */
    void FindConditionedNodes(
        const DRT::Discretization& dis, const std::string& condname, std::vector<int>& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis, const std::string& condname,
        std::map<int, DRT::Node*>& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(
        const DRT::Discretization& dis, const std::string& condname, std::set<int>& nodeset);

    /// find all local nodes from discretization marked with condition
    /*!
      Loop all conditions of the given discretization, find the ones with the
      specified name and collect the locally owned node ids in the suppied
      set. The nodes vector is unique and ordered on output.

      \param dis : (in) Discretization
      \param conds : (in) conditions in question
      \param nodes : (out) empty set on input, filled with nodal gids of local nodes
    */
    void FindConditionedNodes(const DRT::Discretization& dis,
        const std::vector<DRT::Condition*>& conds, std::vector<int>& nodes);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis,
        const std::vector<DRT::Condition*>& conds, std::map<int, DRT::Node*>& nodes);

    /// find all local nodes from discretization marked with condition and
    /// put them into a map indexed by Id of the condition
    void FindConditionedNodes(const DRT::Discretization& dis,
        const std::vector<DRT::Condition*>& conds, std::map<int, std::map<int, DRT::Node*>>& nodes);

    /// find all local nodes from discretization marked with condition and
    /// put them into a vector indexed by Id of the condition
    void FindConditionedNodes(const DRT::Discretization& dis,
        const std::vector<DRT::Condition*>& conds,
        std::map<int, Teuchos::RCP<std::vector<int>>>& nodes, bool use_coupling_id = true);

    /// find all local nodes from discretization marked with condition
    void FindConditionedNodes(const DRT::Discretization& dis,
        const std::vector<DRT::Condition*>& conds, std::set<int>& nodeset);


    /// collect all local nodes and elements in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param elements unique map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis, std::map<int, DRT::Node*>& nodes,
        std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::string& condname);

    /// collect all nodes (in- and excluding 'ghosts') and
    /// elements (including ghosts) in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param ghostnodes overlapping map of nodes
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis, std::map<int, DRT::Node*>& nodes,
        std::map<int, DRT::Node*>& gnodes, std::map<int, Teuchos::RCP<DRT::Element>>& elements,
        const std::vector<DRT::Condition*>& conds);

    /// collect all elements in a condition including ghosts
    /*!
      \param elements overlapping map of elements
      \param vector containing condition pointers
     */
    void FindConditionObjects(std::map<int, Teuchos::RCP<DRT::Element>>& elements,
        const std::vector<DRT::Condition*>& conds);

    /// collect all nodes (in- and excluding 'ghosts') and
    /// elements (including ghosts) in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param ghostnodes overlapping map of nodes
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis, std::map<int, DRT::Node*>& nodes,
        std::map<int, DRT::Node*>& gnodes, std::map<int, Teuchos::RCP<DRT::Element>>& elements,
        const std::string& condname);

    /// collect all nodes (in- and excluding 'ghosts') and
    /// elements (including ghosts) in a condition
    /*!
      \param dis discretization
      \param nodes unique map of nodes
      \param ghostnodes overlapping map of nodes
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
        std::map<int, std::map<int, DRT::Node*>>& nodes,
        std::map<int, std::map<int, DRT::Node*>>& gnodes,
        std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>& elements,
        const std::string& condname);

    /// collect all elements in a condition including ghosts
    /*!
      \param dis discretization
      \param elements overlapping map of elements
      \param condname name of condition
     */
    void FindConditionObjects(const DRT::Discretization& dis,
        std::map<int, Teuchos::RCP<DRT::Element>>& elements, const std::string& condname,
        const int label = -1);


    /// Create element map to given condition
    /*!
      \pre Geometry must be available in condition

      \param dis : (in) Discretization
      \param condname : (in) name of condition in question
      \param colmap : (in) whether we want to create a (ghosted) column map
    */
    Teuchos::RCP<Epetra_Map> ConditionElementMap(
        const DRT::Discretization& dis, const std::string& condname, bool colmap);

    /// Find all conditions with given name that all nodes of the element have in common
    /*!
      \param ele (in) the element
      \param condname (in) name of the condition to look for
      \param condition (out) all conditions that cover all element nodes
    */
    void FindElementConditions(const DRT::Element* ele, const std::string& condname,
        std::vector<DRT::Condition*>& condition);

    /// row map with nodes from condition
    Teuchos::RCP<Epetra_Map> ConditionNodeRowMap(
        const DRT::Discretization& dis, const std::string& condname);

    /// col map with nodes from condition
    Teuchos::RCP<Epetra_Map> ConditionNodeColMap(
        const DRT::Discretization& dis, const std::string& condname);

    /// create the set of column element gids that have conditioned nodes
    /*!
      \note These are not elements from the condition geometry. Rather the
      gids of actual discretization elements are listed.
     */
    Teuchos::RCP<std::set<int>> ConditionedElementMap(
        const DRT::Discretization& dis, const std::string& condname);

    /// Writes debug output for volumetric coupling
    void WriteBoundarySurfacesVolumeCoupling(
        std::map<std::vector<int>, Teuchos::RCP<DRT::Element>>
            surfmap,  // map containing the surfaces of the boundary condition
        int condID,   // condition ID
        int numproc,  // number of procs
        int mypid     // current proc
    );

    /*!
     * \brief This method checks whether handed in conditions are defined on the same set of nodes
     *
     * @param[in] condition1  first condition to be tested
     * @param[in] condition2  second condition to be tested
     * @param[in] mustmatch   both conditions must match
     * @return flag indicating if both conditions are defined on the same set of nodes
     */
    bool HaveSameNodes(const DRT::Condition* const condition1,
        const DRT::Condition* const condition2, bool mustmatch);

  }  // namespace UTILS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
