/*----------------------------------------------------------------------*/
/*!
\file strutimint_state.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_STATE_H
#define STRUTIMINT_STATE_H

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Object holds multiple displacements, velocities and
   *    accelerations vectors (or similiar), of structural time 
   *        integration scheme, as they occur due to multiple step or stage 
   *        character of the scheme.
   *        (Alright, accelerations do not belong truely to the phase state. 
   *        Moreover, the state is also constituted by the elementwise stored 
   *        hybrid/mixed quantities.)
   *        The main purpose arises from time-step size adaptivity
   *        in which the auxiliar scheme might be a multi-step scheme
   *        and thus needs to access older states, which are not directly
   *        needed by the -- for instance single-step -- marching method.
   * \author bborn
   * \date 07/08
   */
  class StruTimIntState
  {

  public:

    //! Dummy constructor
    StruTimIntState();

    //! Constructor
    StruTimIntState
    (
      const int indexlow,  //!< lower index bound
      const int indexhigh,  //!< higher index bound, >= lower bound
      const Epetra_Map* dofrowmap,  //!< vector layout from discretization
      const bool inittozero  //!< initialse to zero, if true
    );

    //! Resize
    void Resize
    (
      const int indexlow,  //!< lower index bound
      const int indexhigh,  //!< higher index bound, >= lower bound
      const bool inittozero  //!< initialse to zero, if true
    );

    //! @name Actions
    //@{

    //! Update multi-step state
    //! i.e. state_{n} := state_{n+1}, 
    //!      state_{n-1} := state_{n},
    //!      state_{n} := state_{n+1},
    //!      etc.
    void UpdateSteps
    (
      const Teuchos::RCP<Epetra_Vector> staten  //!< state_{n+1}
    );

    //@}

    //! @name Access functions
    //@{

    //! Access state vectors by time step index
    Teuchos::RCP<Epetra_Vector> operator()
    (
      const int step  //!< inquire step
    )
    { 
      if (not StepExists(step)) dserror("Step %d is not admissible", step);
      return state_[IndexByStep(step)];
    }

    //! Access state vectors by default time step
    Teuchos::RCP<Epetra_Vector> operator()()
    { 
      return operator()(0);
    }

    //@}
  

  protected:

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

    //! @name Steps
    //@{

    //! lowest number
    int steppast_;
    //! highest number
    int stepfuture_;
    //! total
    int steps_;

    //@}

    //! vector layout from discretization
    const Epetra_Map* dofrowmap_;

    //! multiple state vectors, i.e. displ., veloc., accel.
    std::vector< Teuchos::RCP<Epetra_Vector> > state_;

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_STATE_H
#endif  // #ifdef CCADISCRET
