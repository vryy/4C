/*----------------------------------------------------------------------*/
/*!
\file strutimint_state.cpp

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "strutimint_state.H"


/*----------------------------------------------------------------------*/
/* Dummy constructor */
STR::StruTimIntState::StruTimIntState
(
)
: steppast_(0),
  stepfuture_(0),
  steps_(0),
  dofrowmap_(NULL),
  state_(Teuchos::null)
{
  // that's it
  return;
}

/*----------------------------------------------------------------------*/
/* constructor */
STR::StruTimIntState::StruTimIntState
(
  const int steppast,
  const int stepfuture,
  const Epetra_Map* dofrowmap,
  const bool inittozero
)
: steppast_(steppast),
  stepfuture_(stepfuture),
  steps_(stepfuture>=steppast ? stepfuture-steppast+1 : 0),
  dofrowmap_(dofrowmap),
  state_()
{
  // verify a positive #length_
  dsassert(steps_>0, "Past step must be lower or equal to future step");

  // allocate place for steps_-times state_
  state_.resize(steps_);

  // allocate the vectors themselves
  for (int index=0; index<steps_; ++index)
  {
    state_[index] = LINALG::CreateVector(*dofrowmap_, inittozero);
  }
  
  // good bye
  return;
}

/*----------------------------------------------------------------------*/
/* Resize */
void STR::StruTimIntState::Resize
(
  const int steppast,
  const int stepfuture,
  const bool inittozero
)
{
  // check this
  dsassert(steppast <= stepfuture, 
           "Past step must be lower than future step");
  
  // verify this
  dsassert(stepfuture == stepfuture_,
           "Future step cannot be changed");

  // add states for steps in past
  if (steppast < steppast_)
  {
    for (int past=steppast_; past>steppast; --past)
    {
      state_.insert(state_.begin(),
                  LINALG::CreateVector(*dofrowmap_, inittozero));
    }
    steppast_ = steppast;
    steps_ = stepfuture_ - steppast_ + 1;
  }

  // farewell
  return;
}

/*----------------------------------------------------------------------*/
/* Update steps */
void STR::StruTimIntState::UpdateSteps
(
  const Teuchos::RCP<Epetra_Vector> staten
)
{
  for (int ind=0; ind<steps_-1; ++ind)
  {
    state_[ind]->Update(1.0, *(state_[ind+1]), 0.0);
  }
  state_[steps_-1]->Update(1.0, *staten, 0.0);

  // ciao
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
