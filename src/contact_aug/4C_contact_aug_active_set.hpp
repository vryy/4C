/*----------------------------------------------------------------------------*/
/*! \file
\brief Identify the correct active set.


\level 2
*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_ACTIVE_SET_HPP
#define FOUR_C_CONTACT_AUG_ACTIVE_SET_HPP

#include "4C_config.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class ParamsInterface;
  namespace Aug
  {
    class Strategy;

    /// \brief class to handle the active set identification
    /**
     *  \author hiermeier \date 03/18 */
    class ActiveSet
    {
      enum class Status : int
      {
        unevaluated = -1,
        unchanged = 0,
        changed = 1
      };

      static std::string status2_string(enum Status status)
      {
        switch (status)
        {
          case Status::unevaluated:
            return "Status::unevaluated";
          case Status::unchanged:
            return "Status::unchanged";
          case Status::changed:
            return "Status::changed";
          default:
            return "Unknown active set status";
        }
      }

     public:
      /// constructor
      ActiveSet(Strategy& strategy) : strategy_(strategy){};

      /** @brief Update the active set status
       *
       *  @param(in) cparams  Contact parameter interface
       *
       *  @author hiermeier @date 03/18 */
      void Update(const CONTACT::ParamsInterface& cparams);

      /** @brief print the active set status
       *
       *  @param(out) os modified screen/file output stream
       *
       *  @author hiermeier @date 03/18 */
      void Print(std::ostream& os) const;

     private:
      /// get the global update status
      Status update_status(const CONTACT::ParamsInterface& cparams) const;

      /// @brief set the initial active status
      /** This routine is only executed in the very first step or directly after
       *  a restart. Furthermore, the active set of an interface is only changed
       *  if no natural change has been detected.
       *
       *  \note If a TangDisConstFext predictor is used, it is possible that the
       *  first execution does not change anything, since the contact is not
       *  initiated by a displacement change but rather by a Neumann load change.
       *  However, if the displacement change does not naturally set any node
       *  active, this scenario should be correctly treated.
       *
       *  @param(in) cparams  Contact parameter interface
       *  @param(in) istatus  Interface active set status vector
       *
       *  \author hiermeier \date 03/18 */
      Status update_initial_status(
          const CONTACT::ParamsInterface& cparams, const std::vector<enum Status>& istatus) const;

      /// update related node and dof maps
      void update_maps(const CONTACT::ParamsInterface& cparams);

      /// Merge the active set status of different contact interfaces
      Status merge(const std::vector<Status>& istatus) const;

      /// Perform a sanity check
      void sanity_check(const CONTACT::ParamsInterface& cparams, const enum Status gstatus) const;

      /// @brief Check if the active set update shall be skipped
      /**
       *  @return TRUE if the active set update shall be skipped.
       *
       *  @author hiermeier @date 03/18 */
      bool skip_update() const;

      /// post update routine
      void post_update(const CONTACT::ParamsInterface& cparams, const enum Status gstatus);

     private:
      /// reference to the underlying augmented contact strategy
      Strategy& strategy_;
    };
  }  // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
