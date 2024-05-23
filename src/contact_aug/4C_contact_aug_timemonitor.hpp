/*----------------------------------------------------------------------------*/
/*! \file
\brief Fast time monitor. E. g. to measure the element evaluation times.


\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_TIMEMONITOR_HPP
#define FOUR_C_CONTACT_AUG_TIMEMONITOR_HPP

#include "4C_config.hpp"

#include <Epetra_Comm.h>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
  {
    /** enumerator list for the inner CONTACT::AUG::Integrator functions
     *  which can be found in the respective policy */
    enum class TimeID : unsigned
    {
      deriv2nd_non_unit_slave_element_normal,
      deriv2nd_unit_slave_element_normal,
      Deriv2nd_Jacobian,
      deriv1st_non_unit_slave_element_normal,
      Deriv1st_MXiGP,
      INCOMPLETE_Add_Jac_Deriv2nd_GapN,
      INCOMPLETE_Add_Deriv1st_GapN_Deriv1st_Jac,
      /*--- MAX_VALUE: Must stay at the end ------------*/
      MAX_TIME_ID
    };

    /** enumerator list of the \"global\" CONTACT::AUG::Integrator functions
     *  or other element specific time consuming functions. */
    enum class GlobalTimeID : unsigned
    {
      integrate_deriv_slave_element,
      IntegrateDerivEle2D,
      integrate_deriv_segment2_d,
      IntegrateDerivEle3D,
      integrate_deriv_cell3_d_aux_plane,
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
        case TimeID::deriv2nd_non_unit_slave_element_normal:
          return "deriv2nd_non_unit_slave_element_normal";
        case TimeID::deriv2nd_unit_slave_element_normal:
          return "deriv2nd_unit_slave_element_normal";
        case TimeID::Deriv2nd_Jacobian:
          return "Deriv2nd_Jacobian";
        case TimeID::deriv1st_non_unit_slave_element_normal:
          return "deriv1st_non_unit_slave_element_normal";
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
      /// constructor
      TimeMonitor();

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
      /** The pairs hold: first -> last start wall time, second -> accumulated
       *  time since the last reset call. The vector entries correspond to the
       *  different involved enumerators */
      std::vector<std::pair<double, double>> timings_;

      /// time increment between the very last start and stop call
      double last_incr_ = 0.0;

      /// pointer to the Epetra communciator object (has the ownership)
      std::shared_ptr<const Epetra_Comm> comm_;
    };

    typedef TimeMonitor<GlobalTimeID> GlobalTimeMonitor;
  }  // namespace AUG
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
