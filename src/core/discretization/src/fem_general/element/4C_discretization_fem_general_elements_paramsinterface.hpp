/*-----------------------------------------------------------*/
/*! \file

\brief collection of a serial of possible element actions to pass information
       down to elements

\level 3


*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_ELEMENTS_PARAMSINTERFACE_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_ELEMENTS_PARAMSINTERFACE_HPP

#include "4C_config.hpp"

#include "4C_legacy_enum_definitions_element_actions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::UTILS
{
  class FunctionManager;
}

namespace Core::Elements
{

  /*! \brief Parameter interface for the element <--> time integrator data exchange
   *
   *  Pure virtual interface class. This class is supposed to replace the current
   *  tasks of the Teuchos::ParameterList.
   *  Please consider to derive a special interface class, if you need special parameters inside
   *  of your element. Keep the Evaluate call untouched and cast the interface object to the
   *  desired specification when and where you need it.
   *
   *  ToDo Currently we set the interface in the elements via the Teuchos::ParameterList.
   *  Theoretically, the Teuchos::ParameterList can be replaced by the interface itself!
   *
   *  \date 03/2016
   *  \author hiermeier */
  class ParamsInterface
  {
   public:
    //! destructor
    virtual ~ParamsInterface() = default;

    //! @name Access general control parameters
    //! @{
    //! get the desired action type
    virtual enum ActionType GetActionType() const = 0;

    //! get the current total time for the evaluate call
    virtual double GetTotalTime() const = 0;

    //! get the current time step
    virtual double GetDeltaTime() const = 0;

    //! get function manager
    virtual const Core::UTILS::FunctionManager* get_function_manager() const = 0;
    //! @}
  };  // class ParamsInterface
}  // namespace Core::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
