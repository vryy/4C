/*-----------------------------------------------------------------------*/
/*! \file
\brief Contains a summary of contact utility functions


\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_CONTACT_UTILS_H
#define BACI_CONTACT_UTILS_H

#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

// forward declarations
namespace DRT
{
  class Condition;
  class Discretization;
  class Node;
  class Element;
  class Condition;
}  // namespace DRT

namespace CORE::LINALG
{
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace CONTACT
{
  namespace UTILS
  {
    /// Get the solid to solid contact conditions
    int GetContactConditions(std::vector<DRT::Condition*>& contact_conditions,
        const std::vector<DRT::Condition*>& beamandsolidcontactconditions,
        const bool& throw_error = true);

    /// Find the solid to solid contact conditions and combine them to contact condition groups
    int GetContactConditionGroups(std::vector<std::vector<DRT::Condition*>>& ccond_grps,
        const DRT::Discretization& discret_wrapper, const bool& throw_error = true);

    /// Combine the solid to solid contact conditions to contact condition groups
    void GetContactConditionGroups(std::vector<std::vector<DRT::Condition*>>& ccond_grps,
        const std::vector<DRT::Condition*>& cconds);

    /// Gather information which side is master and which side is slave
    void GetMasterSlaveSideInfo(std::vector<bool>& isslave, std::vector<bool>& isself,
        const std::vector<DRT::Condition*>& cond_grp);

    /**
     * \brief Gather information on initialization (Active/Inactive)
     *
     * \param [in,out]  Two_half_pass: two half pass approach applied for current condition group
     * \param [in,out]  Check_nonsmooth_selfcontactsurface: reference configuration check for
     *                  non-smooth self contact shall be performed for current condition group
     * \param [in,out]  Searchele_AllProc: Search elements on all processors
     * \param [in,out]  isactive:  condition is set active
     * \param [in]      isslave:   condition is defined as slave side
     * \param [in]      isself:    condition is self contact condition
     * \param [in]      cond_grp: current contact condition group (i.e. conditions with same ID)
     *
     * \author cschmidt \date 11/18 */
    void GetInitializationInfo(bool& Two_half_pass, bool& Check_nonsmooth_selfcontactsurface,
        bool& Searchele_AllProc, std::vector<bool>& isactive, std::vector<bool>& isslave,
        std::vector<bool>& isself, const std::vector<DRT::Condition*>& cond_grp);

    /// write conservation data to an output file
    void WriteConservationDataToFile(const int mypid, const int interface_id, const int nln_iter,
        const CORE::LINALG::SerialDenseMatrix& conservation_data, const std::string& ofile_path,
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
      static void DetectDbcSlaveNodesAndElements(const DRT::Discretization& str_discret,
          const std::vector<std::vector<DRT::Condition*>>& ccond_grps,
          std::set<const DRT::Node*>& dbc_slave_nodes,
          std::set<const DRT::Element*>& dbc_slave_eles);

     private:
      static void DetectDbcSlaveNodes(std::map<const DRT::Node*, int>& dbc_slave_node_map,
          const DRT::Discretization& str_discret,
          const std::vector<const DRT::Condition*>& sl_conds);

      static void DetectDbcSlaveElements(std::set<const DRT::Element*>& dbc_slave_eles,
          const std::map<const DRT::Node*, int>& dbc_slave_nodes,
          const std::vector<const DRT::Condition*>& sl_conds);
    };  // class DbcHandler

  }  // namespace UTILS
}  // namespace CONTACT


#endif  // CONTACT_UTILS_H
