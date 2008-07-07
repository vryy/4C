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
/* constructor */
StruTimIntState::StruTimIntState
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
  dis_(),
  vel_(),
  acc_()
{
  // verify a positive #length_
  dsassert(steps_>0, "Past step must be lower or equal to future step");

  // allocate place for displacement, velocitiy, etc vectors
  dis_.resize(steps_);
  vel_.resize(steps_);
  acc_.resize(steps_);

  // allocate the displacement vectors themselves
  for (int index=0; index<steps_; ++index)
  {
    dis_[index] = LINALG::CreateVector(*dofrowmap_, inittozero);
    vel_[index] = LINALG::CreateVector(*dofrowmap_, inittozero);
    acc_[index] = LINALG::CreateVector(*dofrowmap_, inittozero);
  }
  
  // good bye
  return;
}

/*----------------------------------------------------------------------*/
/* Resize */
void StruTimIntState::Resize
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
      dis_.insert(dis_.begin(),
                  LINALG::CreateVector(*dofrowmap_, inittozero));
      vel_.insert(vel_.begin(),
                  LINALG::CreateVector(*dofrowmap_, inittozero));
      acc_.insert(acc_.begin(),
                  LINALG::CreateVector(*dofrowmap_, inittozero));
    }
    steppast_ = steppast;
  }

  // farewell
  return;
}

/*----------------------------------------------------------------------*/
/* Update states */
void StruTimIntState::UpdateStep()
{
  for (int ind=0; ind<steps_-1; ++ind)
  {
    dis_[ind]->Update(1.0, *(dis_[ind+1]), 0.0);
    vel_[ind]->Update(1.0, *(vel_[ind+1]), 0.0);
    acc_[ind]->Update(1.0, *(acc_[ind+1]), 0.0);
  }

  // ciao
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
