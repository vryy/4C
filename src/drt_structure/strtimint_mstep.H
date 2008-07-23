/*----------------------------------------------------------------------*/
/*!
\file strtimint_mstep.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_MSTEP_H
#define STRTIMINT_MSTEP_H

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief This is the base object for holding multi-step quantities,
   *        e.g. displacements \f$D_n\f$, \f$D_{n-1}\f$, \f$D_{n-2}\f$,
   *        \f$\ldots\f$, or time points \f$t_n\f$, \f$t_{n-1}\f$,
   *        \f$t_{n-2}\f$, \f$\ldots\f$, etc.
   *
   * Multi-step quantities occur surprisingly indeed in multi-step time
   * integrators. However, they are needed in single-step integrators
   * as well, if the auxiliar scheme is multi-step. The idea is to be
   * able to resize the multi-step quantities dynamically and have 
   * intrinsic update mechanism such that a single-step integrator
   * can provide the data on which a multi-step auxiliar integrator
   * can work.
   *
   * \author bborn
   * \date 07/08
   */
  template <typename STATE> class TimIntMStepBase
  {

  public:

    //! @name Access functions
    //@{

    //! Access state object by time step
    STATE& operator[]
    (
      const int step  //!< inquiry step
    )
    {
      if (not StepExists(step)) dserror("Step %d is not admissible", step);
      return state_[IndexByStep(step)];
    }

    //! Access state object by time step index as RCP
    Teuchos::RCP<STATE> operator()
    (
      const int step  //!< inquire step
    )
    { 
      if (not StepExists(step)) dserror("Step %d is not admissible", step);
      return Teuchos::rcp<STATE>(&(state_[IndexByStep(step)]),false);
    }

    //@}

  protected:

    //! @name Allocation
    //@{

    //! Dummy constructor
    TimIntMStepBase()
    : steppast_(0),
      stepfuture_(0),
      steps_(0),
      state_()
    {
      // that's it
      return;
    }

    //! Another dummy constructor which sets vector limits and its size
    //! but does not allocate anything
    TimIntMStepBase
    (
      const int steppast,  //!< lower index bound
      const int stepfuture  //!< higher index bound, >= lower bound
    )
    : steppast_(steppast),
      stepfuture_(stepfuture),
      steps_(stepfuture>=steppast ? stepfuture-steppast+1 : 0),
      state_()
    {
      // verify a positive #steps_
      dsassert(steps_>0, "Past step must be lower or equal to future step");

      // bye
      return;
    }

    //@}

    //! @name Verfication
    //@{

    //! Check sanity prior resize
    bool ResizeSane
    (
      const int steppast,
      const int stepfuture
    )
    {
      bool sane = true;

      sane = sane and (steppast <= stepfuture);
      if (not sane)
        dserror("Past step must be lower than future step");
      
      sane = sane and (stepfuture == stepfuture_);
      if (not sane)
        dserror("Future step cannot be changed");

      // deliver it
      return sane;
    }

    //! Determine wether step lies in given bounds
    bool StepExists
    (
      const int step  //!< inquired step
    )
    {
      return ( (step >= steppast_) and (step <= stepfuture_) );
    }  

    //! determine wether index lies in given bounds
    bool IndexExists
    (
      const int index //!< inquired index
    )
    {
      return ( (index >= 0) and (index < (int) state_.size()) );
    }

    //! Map step number to vector index
    int StepByIndex
    (
      const int index  //!< vector index
    )
    { 
      return index + steppast_; 
    }

    //! Map vector index to step
    int IndexByStep
    (
      const int step  //!< step number
    )
    {
      int index = step - steppast_;
      dsassert(IndexExists(index), "step is not permissible!");
      return (unsigned) index;
    }

    //@}

    //! @name Steps
    //@{

    //! lowest number
    int steppast_;
    //! highest number
    int stepfuture_;
    //! total
    int steps_;

    //@}

    //! multiple state vectors, i.e. displ., veloc., accel.
    std::vector<STATE> state_;

  };

  /*====================================================================*/
  /*!
   * \brief General derived class for multi-step quantities of
   *        simple type. This going to be used with \c double and \c int
   */
  template <typename STATE> class TimIntMStep
  : public TimIntMStepBase<STATE>
  {

  protected :

    //! template base class
    typedef TimIntMStepBase<STATE> TBase;

  public :

    //! @name Life
    //@{

    //! Constructor
    TimIntMStep
    (
      const int steppast,  //!< lower index bound
      const int stepfuture,  //!< higher index bound, >= lower bound
      const STATE init  //!< initialse to init
    )
    : TBase(steppast, stepfuture)
    {
      // allocate the vectors themselves
      for (int index=0; index<TBase::steps_; ++index)
      {
        TBase::state_.push_back(STATE(init));
      }
  
      // good bye
      return;
    }

    //! Resize
    void Resize
    (
      const int steppast,  //!< lower index bound
      const int stepfuture,  //!< higher index bound, >= lower bound
      const STATE init  //!< initialise to init
    )
    {
      // check this
      bool sane = TBase::ResizeSane(steppast, stepfuture);
      if (not sane) dserror("Sanity check not passed.");
      
      // add states for steps in past
      if (steppast < TBase::steppast_)
      {
        for (int past=TBase::steppast_; past>steppast; --past)
        {
          TBase::state_.insert(TBase::state_.begin(), STATE(init));
        }
        TBase::steppast_ = steppast;
        TBase::steps_ = TBase::stepfuture_ - TBase::steppast_ + 1;
      }

      // farewell
      return;
    }

    //@}

    //! @name Actions
    //@{

    //! Set entrie at #step to #value
    void SetStep
    (
      const int step,
      const STATE value
    )
    {
      TBase::state_.at(TBase::IndexByStep(step)) = value;
    }

    //! Update multi-step state
    //! i.e. state_{n} := state_{n+1}, 
    //!      state_{n-1} := state_{n},
    //!      state_{n} := state_{n+1},
    //!      etc.
    void UpdateSteps
    (
      const STATE staten  //!< state_{n+1}
    )
    {
      for (int ind=0; ind<TBase::steps_-1; ++ind)
      {
        TBase::state_[ind] = TBase::state_[ind+1];
      }
      TBase::state_[TBase::steps_-1] = staten;
      
      // ciao
      return;
    }

    //@}

  };

  /*====================================================================*/
  /*! \brief Explicitly declare #double instantiated class 
   *         for multi-step qunatities of \c double typewhich need 
   *         to be generated by template
   */
  template class TimIntMStep<double>;

  /*====================================================================*/
  /*! \brief Explicitly declare #int instantiated class 
   *         for multi-step qunatities of \c int type which need 
   *         to be generated by template
   */
  template class TimIntMStep<int>;

  /*====================================================================*/
  /*! \brief Specialise general #TimIntMStepBase object 
   *         for \c Epetra_Vector as needed for state vectors
   *         -- displacements, velocities and accelerations
   *         (alright, accelerations do not truely belong to the phase 
   *         state) -- in structural dynamics
   */
  template <> class TimIntMStep<Epetra_Vector>
  : public TimIntMStepBase<Epetra_Vector>
  {

  protected :
    
    //! base template class
    typedef TimIntMStepBase<Epetra_Vector> TBase;

  public :

    //! @name Life
    //@{

    //! Dummy constructor 
    TimIntMStep()
    : TBase::TimIntMStepBase()
    { ; }

    //! Constructor
    TimIntMStep
    (
      const int steppast,  //!< lower index bound
      const int stepfuture,  //!< higher index bound, >= lower bound
      const Epetra_Map* dofrowmap,  //!< vector layout from discretization
      const bool inittozero  //!< initialse to zero, if true
    )
    : TBase(steppast, stepfuture)
    {
      // allocate the vectors themselves
      for (int index=0; index<steps_; ++index)
      {
        state_.push_back(Epetra_Vector(*dofrowmap, inittozero));
      }
  
      // good bye
      return;
    }

    //! Resize
    //!
    //! State vectors are added and placed according to their
    //! indices #steppast to #stepfuture
    void Resize
    (
      const int steppast,  //!< lower index bound
      const int stepfuture,  //!< higher index bound, >= lower bound
      const Epetra_Map* dofrowmap,  //!< vector layout from discretization
      const bool inittozero  //!< initialse to zero, if true
    )
    { 
      // check this
      bool sane = TBase::ResizeSane(steppast, stepfuture);
      if (not sane) dserror("Sanity check not passed.");

      // add states for steps in past
      if (steppast < steppast_)
      {
        for (int past=steppast_; past>steppast; --past)
        {
          state_.insert(state_.begin(), Epetra_Vector(*dofrowmap, inittozero));
        }
        steppast_ = steppast;
        steps_ = stepfuture_ - steppast_ + 1;
      }

      // farewell
      return;
    }

    //@}

    //! @name Actions
    //@{

    //! Update multi-step state
    //! i.e. state_{n} := state_{n+1}, 
    //!      state_{n-1} := state_{n},
    //!      state_{n} := state_{n+1},
    //!      etc.
    void UpdateSteps
    (
      const Epetra_Vector& staten  //!< state_{n+1}
    )
    {
      for (int ind=0; ind<steps_-1; ++ind)
      {
        (state_[ind]).Update(1.0, (state_[ind+1]), 0.0);
      }
      (state_[steps_-1]).Update(1.0, staten, 0.0);
      
      // ciao
      return;
    }

    //@}

  };  // class TimIntMStep<Epetra_Vector>
  
}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_MSTEP_H
#endif  // #ifdef CCADISCRET
