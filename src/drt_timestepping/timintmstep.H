/*----------------------------------------------------------------------*/
/*!
\file timintmstep.H

\brief Multi-step functionalities for time integration

<pre>
\level 0
\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef TIMINTMSTEP_H
#define TIMINTMSTEP_H

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/* A collection of time integration utilities */
namespace TIMINT
{
  /*====================================================================*/
  /*!
   * \brief This is the base object for holding multi-step solution quantities,
   *        e.g. displacements \f$D_n\f$, \f$D_{n-1}\f$, \f$D_{n-2}\f$,
   *        \f$\ldots\f$, or time points \f$t_n\f$, \f$t_{n-1}\f$,
   *        \f$t_{n-2}\f$, \f$\ldots\f$, etc.
   *
   * Multi-step quantities occur --surprisingly indeed-- in multi-step time
   * integrators. However, they are needed in single-step integrators
   * as well, if the auxiliary scheme is multi-step. The idea is to be
   * able to resize the multi-step quantities dynamically and have
   * intrinsic update mechanism such that a single-step integrator
   * can provide the data on which a multi-step auxiliary integrator
   * can work.
   *
   * \author bborn
   * \date 07/08
   */
  template <typename STATE>
  class TimIntMStepBase
  {
   public:
    //! @name Access functions
    //@{

    //! Access state object by time step index
    STATE& operator[](const int step  //!< inquiry step
    )
    {
      if (not StepExists(step)) dserror("Step %d is not admissible", step);
      return state_[IndexByStep(step)];
    }

    //! Access state object by time step index (read-only)
    const STATE& operator[](const int step  //!< inquiry step
        ) const
    {
      if (not StepExists(step)) dserror("Step %d is not admissible", step);
      return state_[IndexByStep(step)];
    }

    //! Access state object by time step index and wrap into Teuchos::RCP
    Teuchos::RCP<STATE> operator()(const int step  //!< inquiry step
    )
    {
      if (not StepExists(step)) dserror("Step %d is not admissible", step);
      return Teuchos::rcp<STATE>(&(state_[IndexByStep(step)]), false);
    }

    //@}

    //! @name Query functions
    //@{

    //! step index of most past and future steps
    std::pair<int, int> GetSteps() { return std::make_pair(steppast_, stepfuture_); }

    //! Number of steps
    int NumSteps() { return steps_; }
    //@}

   protected:
    //! @name Allocation
    //@{

    //! Dummy constructor
    TimIntMStepBase() : steppast_(0), stepfuture_(0), steps_(0), state_() { return; }

    /*! Another dummy constructor which sets vector limits and its size
     *  but does not allocate anything
     */
    TimIntMStepBase(const int steppast,  //!< lower index bound
        const int stepfuture             //!< higher index bound, >= lower bound
        )
        : steppast_(steppast),
          stepfuture_(stepfuture),
          steps_(stepfuture >= steppast ? stepfuture - steppast + 1 : 0),
          state_()
    {
      // verify a positive #steps_
      if (steps_ <= 0) dserror("Past step must be lower or equal to future step");

      return;
    }

    //@}

    //! @name Verification
    //@{

    //! Check sanity prior resize
    bool ResizeSane(const int steppast,  //!< lower index bound
        const int stepfuture             //!< higher index bound, >= lower bound
    )
    {
      bool sane = true;

      sane = sane and (steppast <= stepfuture);
      if (not sane) dserror("Past step must be lower than future step");

      sane = sane and (stepfuture == stepfuture_);
      if (not sane) dserror("Future step cannot be changed");

      return sane;
    }

    //! Determine whether step lies in given bounds
    bool StepExists(const int step  //!< inquired step
        ) const
    {
      return ((step >= steppast_) and (step <= stepfuture_));
    }

    //! determine whether index lies in given bounds
    bool IndexExists(const int index  //!< inquired index
        ) const
    {
      return ((index >= 0) and (index < (int)state_.size()));
    }

    //! Map step number to vector index
    int StepByIndex(const int index  //!< vector index
    )
    {
      return index + steppast_;
    }

    //! Map vector index to step
    int IndexByStep(const int step  //!< step number
        ) const
    {
      int index = step - steppast_;
      dsassert(IndexExists(index), "step is not permissible!");
      return (unsigned)index;
    }

    //@}

    //! @name Steps
    //@{

    int steppast_;    ///< lowest number
    int stepfuture_;  ///< highest number
    int steps_;       ///< total

    //@}

    //! Multi-step quantity, it's stored vectorially
    std::vector<STATE> state_;
  };

  /*====================================================================*/
  /*!
   * \brief General derived class for multi-step quantities of
   *        simple type.
   *        This is going to be used with \c double and \c int
   */
  template <typename STATE>
  class TimIntMStep : public TimIntMStepBase<STATE>
  {
   protected:
    //! template base class
    typedef TimIntMStepBase<STATE> MStepBase;

   public:
    //! @name Life
    //@{

    //! Constructor
    TimIntMStep(const int steppast,  //!< lower index bound
        const int stepfuture,        //!< higher index bound, >= lower bound
        const STATE init             //!< initialise to init
        )
        : MStepBase(steppast, stepfuture)
    {
      // allocate the vectors themselves
      for (int index = 0; index < MStepBase::steps_; ++index)
      {
        MStepBase::state_.push_back(STATE(init));
      }

      return;
    }

    //! Resize
    void Resize(const int steppast,  //!< lower index bound
        const int stepfuture,        //!< higher index bound, >= lower bound
        const STATE init             //!< initialise to init
    )
    {
      // check this
      {
        bool sane = MStepBase::ResizeSane(steppast, stepfuture);
        if (not sane) dserror("Sanity check not passed.");
      }

      // add states for steps in past
      if (steppast < MStepBase::steppast_)
      {
        for (int past = MStepBase::steppast_; past > steppast; --past)
        {
          MStepBase::state_.insert(MStepBase::state_.begin(), STATE(init));
        }
        MStepBase::steppast_ = steppast;
        MStepBase::steps_ = MStepBase::stepfuture_ - MStepBase::steppast_ + 1;
      }

      return;
    }

    //@}

    //! @name Actions
    //@{

    //! Set entry at #step to #value
    void SetStep(const int step, const STATE value)
    {
      MStepBase::state_.at(MStepBase::IndexByStep(step)) = value;
    }

    /*! \brief Update multi-step state
     *
     *  i.e. state_{n} := state_{n+1},
     *       state_{n-1} := state_{n},
     *       state_{n} := state_{n+1},
     *       etc.
     */
    void UpdateSteps(const STATE staten  //!< state_{n+1}
    )
    {
      for (int ind = 0; ind < MStepBase::steps_ - 1; ++ind)
      {
        MStepBase::state_[ind] = MStepBase::state_[ind + 1];
      }
      MStepBase::state_[MStepBase::steps_ - 1] = staten;

      return;
    }

    //@}
  };

  /*====================================================================*/
  /*! \brief Specialize general #TimIntMStepBase object
   *         for \c Epetra_Vector as needed for state vectors
   *         like displacements, velocities and accelerations
   */
  template <>
  class TimIntMStep<Epetra_Vector> : public TimIntMStepBase<Epetra_Vector>
  {
   protected:
    //! base template class
    typedef TimIntMStepBase<Epetra_Vector> MStepBase;

   public:
    //! @name Life
    //@{

    //! Dummy constructor
    TimIntMStep() : MStepBase::TimIntMStepBase() { ; }

    //! Constructor
    TimIntMStep(const int steppast,   //!< lower index bound
        const int stepfuture,         //!< higher index bound, >= lower bound
        const Epetra_Map* dofrowmap,  //!< vector layout from discretization
        const bool inittozero         //!< initialise to zero, if true
        )
        : MStepBase(steppast, stepfuture)
    {
      // allocate the vectors themselves
      for (int index = 0; index < steps_; ++index)
      {
        state_.push_back(Epetra_Vector(*dofrowmap, inittozero));
      }

      return;
    }

    /*! \brief Resize
     *
     *  State vectors are added and placed according to their
     *  indices #steppast to #stepfuture
     */
    void Resize(const int steppast,   //!< lower index bound
        const int stepfuture,         //!< higher index bound, >= lower bound
        const Epetra_Map* dofrowmap,  //!< vector layout from discretization
        const bool inittozero         //!< initialise to zero, if true
    )
    {
      // check this
      {
        bool sane = MStepBase::ResizeSane(steppast, stepfuture);
        if (not sane) dserror("Sanity check not passed.");
      }

      // add states for steps in past
      if (steppast < steppast_)
      {
        for (int past = MStepBase::steppast_; past > steppast; --past)
        {
          MStepBase::state_.insert(
              MStepBase::state_.begin(), Epetra_Vector(*dofrowmap, inittozero));
        }
        MStepBase::steppast_ = steppast;
        MStepBase::steps_ = MStepBase::stepfuture_ - MStepBase::steppast_ + 1;
      }

      return;
    }

    /*! \brief Replace maps and initialize to zero
     *
     *  State vectors cleared and rebuild with given map
     *  take care that underlying discret_ contains the same maps
     */
    void ReplaceMaps(const Epetra_Map* dofrowmap  //!< new vector layout
    )
    {
      state_.clear();
      // allocate the vectors themselves
      for (int index = 0; index < steps_; ++index)
      {
        state_.push_back(Epetra_Vector(*dofrowmap, true));
      }

      return;
    }

    //@}

    //! @name Actions
    //@{

    /*! \brief Update multi-step state
     *
     *  i.e. state_{n} := state_{n+1},
     *       state_{n-1} := state_{n},
     *       state_{n} := state_{n+1},
     *       etc.
     */
    void UpdateSteps(const Epetra_Vector& staten  //!< state_{n+1}
    )
    {
      for (int ind = 0; ind < MStepBase::steps_ - 1; ++ind)
      {
        (MStepBase::state_[ind]).Update(1.0, (MStepBase::state_[ind + 1]), 0.0);
      }
      (MStepBase::state_[steps_ - 1]).Update(1.0, staten, 0.0);

      return;
    }

    //@}

  };  // class TimIntMStep<Epetra_Vector>
}  // namespace TIMINT

/*----------------------------------------------------------------------*/
#endif  // #ifndef TIMINTMSTEP_H
