/*----------------------------------------------------------------------*/
/*!
 \file drt_elements_paramsminimal.H

\brief Minimal implementation of the parameter interface for the element <--> time integrator data exchange

\level 3

\maintainer Michael Hiermeier
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_ELEMENTS_PARAMSMINIMAL_H_
#define SRC_DRT_LIB_DRT_ELEMENTS_PARAMSMINIMAL_H_

#include "drt_elements_paramsinterface.H"

namespace DRT {
namespace ELEMENTS {

  /*!
     \brief Minimal implementation of the parameter interface for the element <--> time integrator data exchange
    */
  class ParamsMinimal : public ParamsInterface
  {
  public:
    //! constructor
    ParamsMinimal()
        : ele_action_(none),
          total_time_(-1.0),
          delta_time_(-1.0)
        {};

    //! destructor
    virtual ~ParamsMinimal() {};

    //! @name Access general control parameters
    //! @{
    //! get the desired action type
    virtual enum ActionType GetActionType() const {return ele_action_;};

    //! get the current total time for the evaluate call
    virtual double GetTotalTime() const {return total_time_;};

    //! get the current time step
    virtual double GetDeltaTime() const {return delta_time_;};
    //! @}

    /*! @name set routines which are used to set the parameters of the data container
     *
     *  These functions are not allowed to be called by the elements! */
    //! @{
    //! set the action type
    inline void SetActionType(const enum DRT::ELEMENTS::ActionType& actiontype)
    {
      ele_action_ = actiontype;
    }

    //! set the total time for the evaluation call
    inline void SetTotalTime(const double& total_time)
    {
      total_time_ = total_time;
    }

    //! set the current time step for the evaluation call
    inline void SetDeltaTime(const double& dt)
    {
      delta_time_ = dt;
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
    //! @}
  };  // class ParamsMinimal

} // namespace ELEMENTS
} // namespace DRT


#endif /* SRC_DRT_LIB_DRT_ELEMENTS_PARAMSMINIMAL_H_ */
