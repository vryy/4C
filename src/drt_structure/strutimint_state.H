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

/*======================================================================*/
/*!
\brief Objects holds phase states, i.e. displacements, velocities and
       accelerations, of structural time integration scheme.
       (Alright, accelerations do not belong truely to the phase state. 
       Moreover, the state is also constituted by the elementwise stored 
       hybrid/mixed quantities.)
       The main purpose arises from time-step size adaptivity
       in which the auxiliar scheme ought to access the states
       originallt computed by the marching method.
\author bborn
\date 07/08
*/
class StruTimIntState
{

public:

  //! constructor
  StruTimIntState
  (
    const int indexlow,  //!< lower index bound
    const int indexhigh,  //!< higher index bound, >= lower bound
    const Epetra_Map* dofrowmap,  //!< vector layout from discretization
    const bool inittozero  //!< initialse to zero, if true
  );

  //! Update states
  //! ... state_{n-1} := state_{n}, state_{n} := state_{n+1}, etc
  void UpdateStep();

  //! @name Access state vectors by index
  //@{

  //! Displacements at index
  Teuchos::RCP<Epetra_Vector> Dis
  (
    const int step  //!< inquire step
  )
  { 
    dsassert(StepExists(step), "Step is not admissible");
    return dis_[IndexByStep(step)];
  }

  //! Velocitis at index 
  Teuchos::RCP<Epetra_Vector> Vel
  (
    const int step  //!< inquire step
  )
  { 
    dsassert(StepExists(step), "Step is not admissible");
    return vel_[IndexByStep(step)]; 
  }

  //! Accelerations at index
  Teuchos::RCP<Epetra_Vector> Acc
  (
    const int step  //!< inquired step
  ) 
  { 
    dsassert(StepExists(step), "Step is not admissible");
    return acc_[IndexByStep(step)];
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
    return ( (index >= 0) and (index < (int) dis_.size()) );
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

  //! Lowest index
  int steppast_;
  int stepfuture_;
  int steps_;

  //! vector layout from discretization
  const Epetra_Map* dofrowmap_;

  //! displacement states
  std::vector<Teuchos::RCP<Epetra_Vector> > dis_;

  //! velcocity states
  std::vector<Teuchos::RCP<Epetra_Vector> > vel_;

  //! acceleration states
  std::vector<Teuchos::RCP<Epetra_Vector> > acc_;

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_STATE_H
#endif  // #ifdef CCADISCRET
