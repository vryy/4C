/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_timemonitor.H

\brief Fast time monitor. E. g. to measure the element evaluation times.

\maintainer Michael Hiermeier

\date May 23, 2018

\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_TIMEMONITOR_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_TIMEMONITOR_H_

#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include <Epetra_Comm.h>

namespace CONTACT
{
  namespace AUG
  {
    enum class TimeID : unsigned
    {
      Deriv2nd_NonUnitSlaveElementNormal = 0,
      Deriv2nd_UnitSlaveElementNormal = 1,
      Deriv2nd_Jacobian = 2,
      Deriv1st_NonUnitSlaveElementNormal = 3,
      Deriv1st_MXiGP = 4,
      INCOMPLETE_Add_Jac_Deriv2nd_GapN = 5,
      INCOMPLETE_Add_Deriv1st_GapN_Deriv1st_Jac = 6,
      /*--- MAX_VALUE: Must stay at the end ------------*/
      MAX_TIME_ID
    };

    enum class GlobalTimeID : unsigned
    {
      IntegrateDerivSlaveElement = 0,
      IntegrateDerivEle2D = 1,
      IntegrateDerivSegment2D = 2,
      IntegrateDerivEle3D = 3,
      IntegrateDerivCell3DAuxPlane = 4,
      /*--- MAX_VALUE: Must stay at the end ------------*/
      MAX_TIME_ID
    };

    /** \brief Basic enum --> string converter
     *
     *  See the specializations for more sophisticated solutions.
     *  \author hiermeier \date 05/18 */
    template <typename enum_class>
    inline std::string TimeID2Str(enum_class id)
    {
      static_assert(std::is_same<unsigned, typename std::underlying_type<enum_class>::type>::value,
          "The template ENUM_CLASS must use UNSIGNED INT as underlying type!");

      std::ostringstream oss;
      oss << "ENUM-#" << static_cast<unsigned>(id);
      return oss.str();
    };

    /** \brief Specialization of the enum --> string converter
     *  \author hiermeier \date 05/18 */
    template <>
    inline std::string TimeID2Str(TimeID id)
    {
      switch (id)
      {
        case TimeID::Deriv2nd_NonUnitSlaveElementNormal:
          return "Deriv2nd_NonUnitSlaveElementNormal";
        case TimeID::Deriv2nd_UnitSlaveElementNormal:
          return "Deriv2nd_UnitSlaveElementNormal";
        case TimeID::Deriv2nd_Jacobian:
          return "Deriv2nd_Jacobian";
        case TimeID::Deriv1st_NonUnitSlaveElementNormal:
          return "Deriv1st_NonUnitSlaveElementNormal";
        case TimeID::Deriv1st_MXiGP:
          return "Deriv1st_MXiGP";
        case TimeID::INCOMPLETE_Add_Jac_Deriv2nd_GapN:
          return "INCOMPLETE_Add_Jac_Deriv2nd_GapN";
        case TimeID::INCOMPLETE_Add_Deriv1st_GapN_Deriv1st_Jac:
          return "INCOMPLETE_Add_Deriv1st_GapN_Deriv1st_Jac";
        default:
          return "UNKNOWN";
      }
    }

    /** \brief Execution time measuring with minimal overhead
     *  \author hiermeier \date 5/18 */
    template <typename enum_class>
    class TimeMonitor
    {
     public:
      TimeMonitor();
      ~TimeMonitor() = default;

      /// set communicator
      void setComm(const Epetra_Comm* comm)
      {
        /* The communicator must be copied, since the ownership must lie in this
         * object (the time monitor), since the time monitor is maybe deleted
         * after the communicator object.                             hiermeier */
        comm_ = std::shared_ptr<const Epetra_Comm>(comm->Clone());
      }

      /// start the timer for the corresponding enumerator
      void start(const enum_class id);

      /// stop the timer for the corresponding enumerator
      void stop(const enum_class id);

      /// write the result overview to the provided stream (screen or file)
      void write(std::ostream& os);

      /// reset the time monitor
      void reset();

      /** get the accumulated time over all enumerators of the ENUM_CLASS on
       *  this processor */
      double getMyTotalTime() const;

      /** \brief return the time interval between the very last start and stop call
       *
       *  \note \"Last\" means in this context, really the last start and stop
       *  call on the TimeMonitor object, independently of the enumerator.
       *
       *  \author hiermeier \date 05/18 */
      inline double getLastTimeIncr() const { return last_incr_; };

     private:
      std::vector<std::pair<double, double>> timings_;
      double last_incr_ = 0.0;
      std::shared_ptr<const Epetra_Comm> comm_;
    };

    typedef TimeMonitor<GlobalTimeID> GlobalTimeMonitor;
  }  // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_TIMEMONITOR_H_ */
