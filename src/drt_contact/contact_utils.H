/*----------------------------------------------------------------------------*/
/**
\file contact_utils.H

\brief Contains a summary of contact utility functions

\maintainer Michael Hiermeier

\date Jun 15, 2016

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_UTILS_H_
#define SRC_DRT_CONTACT_CONTACT_UTILS_H_

#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Condition;
  class DiscretizationInterface;
  class Node;
  class Element;
  class Condition;
}  // namespace DRT
namespace LINALG
{
  class SerialDenseMatrix;
}  // namespace LINALG
namespace CONTACT
{
  namespace UTILS
  {
    /// Get the solid to solid contact conditions
    int GetContactConditions(std::vector<DRT::Condition*>& cond_group,
        const std::vector<DRT::Condition*> beamandsolidcontactconditions,
        const bool& throw_error = true);

    /// Find the solid to solid contact conditions and combine them to contact condition groups
    int GetContactConditionGroups(std::vector<std::vector<DRT::Condition*>>& ccond_grps,
        const DRT::DiscretizationInterface& discret_wrapper, const bool& throw_error = true);

    /// Combine the solid to solid contact conditions to contact condition groups
    void GetContactConditionGroups(std::vector<std::vector<DRT::Condition*>>& ccond_grps,
        const std::vector<DRT::Condition*>& cconds);

    /// Gather information which side is master and which side is slave
    void GetMasterSlaveSideInfo(std::vector<bool>& isslave, std::vector<bool>& isself,
        const std::vector<DRT::Condition*> cond_grp);

    /// write conservation data to an output file
    void WriteConservationDataToFile(const int mypid, const int interface_id, const int nln_iter,
        const LINALG::SerialDenseMatrix& conservation_data, const std::string& ofile_path,
        const std::string& prefix);

    /** \brief Detect DBC slave nodes and elements
     *
     *  Check all slave contact conditions. If the optional condition tag
     *  "RemoveDBCSlaveNodes" can be found in the slave condition line,
     *  all slave nodes and adjacent elements are added to the corresponding
     *  sets.
     *
     *  A possible condition line can look like
     *  E 7 - 1 Slave Inactive FrCoeffOrBound 0.0 AdhesionBound 0.0 Solidcontact RemoveDBCSlaveNodes
     *
     *  \author hiermeier \date 01/18 */
    class DbcHandler
    {
     public:
      /// remove constructor and destructor
      DbcHandler() = delete;
      ~DbcHandler() = delete;

      /** \brief Detect all slave nodes and elements which hold Dbc information
       *
       *  \param(in)  str_discret: structural discretization
       *  \param(in)  ccond_grps:  contact condition groups
       *  \param(out) dbc_slave_nodes: set containing all slave nodes which hold
       *                               DBC information
       *  \param(out) dbc_slave_eles: set containing all slave elements which
       *                              contain at least one DBC slave node
       *
       *  \author hiermeier \date 01/18 */
      static void DetectDbcSlaveNodesAndElements(const DRT::DiscretizationInterface& str_discret,
          const std::vector<std::vector<DRT::Condition*>>& ccond_grps,
          std::set<const DRT::Node*>& dbc_slave_nodes,
          std::set<const DRT::Element*>& dbc_slave_eles);

     private:
      static void DetectDbcSlaveNodes(std::map<const DRT::Node*, int>& dbc_slave_node_map,
          const DRT::DiscretizationInterface& str_discret,
          const std::vector<const DRT::Condition*>& sl_conds);

      static void DetectDbcSlaveElements(std::set<const DRT::Element*>& dbc_slave_nodes,
          const std::map<const DRT::Node*, int>& dbc_slave_node_map,
          const std::vector<const DRT::Condition*>& sl_conds);
    };  // class DbcHandler

  }  // namespace UTILS
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_CONTACT_UTILS_H_ */
