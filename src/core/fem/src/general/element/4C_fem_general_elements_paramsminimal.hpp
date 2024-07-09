/*----------------------------------------------------------------------*/
/*! \file

\brief Minimal implementation of the parameter interface for the element
<--> time integrator data exchange

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_ELEMENTS_PARAMSMINIMAL_HPP
#define FOUR_C_FEM_GENERAL_ELEMENTS_PARAMSMINIMAL_HPP

#include "4C_config.hpp"

#include "4C_fem_general_elements_paramsinterface.hpp"
#include "4C_utils_function_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  /*!
     \brief Minimal implementation of the parameter interface for the element <--> time integrator
     data exchange
    */
  class ParamsMinimal : public ParamsInterface
  {
   public:
    //! constructor
    ParamsMinimal() : ele_action_(none), total_time_(-1.0), delta_time_(-1.0){};

    enum ActionType get_action_type() const override { return ele_action_; };

    double get_total_time() const override { return total_time_; };

    double get_delta_time() const override { return delta_time_; };

    const Core::UTILS::FunctionManager* get_function_manager() const override
    {
      return function_manager_;
    }

    /*! @name set routines which are used to set the parameters of the data container
     *
     *  These functions are not allowed to be called by the elements! */
    //! @{
    //! set the action type
    inline void set_action_type(const enum Core::Elements::ActionType& actiontype)
    {
      ele_action_ = actiontype;
    }

    //! set the total time for the evaluation call
    inline void set_total_time(const double& total_time) { total_time_ = total_time; }

    //! set the current time step for the evaluation call
    inline void set_delta_time(const double& dt) { delta_time_ = dt; }

    //! store function manager
    void set_function_manager(const Core::UTILS::FunctionManager& function_manager)
    {
      function_manager_ = &function_manager;
    }

    //! @}

   private:
    //! @name General element control parameters
    //! @{
    //! Current action type
    enum ActionType ele_action_;

    //! total time for the evaluation
    double total_time_;

    //! current time step for the evaluation
    double delta_time_;

    //! function manager
    const Core::UTILS::FunctionManager* function_manager_;
    //! @}
  };  // class ParamsMinimal

}  // namespace Core::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
