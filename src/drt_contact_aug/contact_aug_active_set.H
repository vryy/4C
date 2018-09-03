/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_active_set.H

\brief Identify the correct active set.

\maintainer Michael Hiermeier

\date Mar 14, 2018

\level 2
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_ACTIVE_SET_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_ACTIVE_SET_H_

#include <string>
#include <vector>

namespace CONTACT
{
  class ParamsInterface;
  namespace AUG
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

      static std::string Status2String(enum Status status)
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

      /// destructor
      ~ActiveSet() = default;

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
      Status UpdateStatus(const CONTACT::ParamsInterface& cparams) const;

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
      Status UpdateInitialStatus(
          const CONTACT::ParamsInterface& cparams, const std::vector<enum Status>& istatus) const;

      /// update related node and dof maps
      void UpdateMaps(const CONTACT::ParamsInterface& cparams);

      /// Merge the active set status of different contact interfaces
      Status Merge(const std::vector<Status>& istatus) const;

      /// Perform a sanity check
      void SanityCheck(const CONTACT::ParamsInterface& cparams, const enum Status gstatus) const;

      /// @brief Check if the active set update shall be skipped
      /**
       *  @return TRUE if the active set update shall be skipped.
       *
       *  @author hiermeier @date 03/18 */
      bool SkipUpdate() const;

      /// post update routine
      void PostUpdate(const CONTACT::ParamsInterface& cparams, const enum Status gstatus);

     private:
      /// reference to the underlying augmented contact strategy
      Strategy& strategy_;
    };
  }  // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_ACTIVE_SET_H_ */
