#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef RESULTTEST

#include "fsi_ale_resulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
}


FSI::AleResultTest::AleResultTest(AleLinear& ale)
  : ale_(ale)
{
}

void FSI::AleResultTest::TestNode(RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one ale discretization supported for testing");

  if (ale_.discret_->HaveGlobalNode(res->node))
  {
    DRT::Node* actnode = ale_.discret_->gNode(res->node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != ale_.discret_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& dispnpmap = ale_.dispnp_->Map();

    string position = res->position;
    if (position=="dispx")
    {
      result = (*ale_.dispnp_)[dispnpmap.LID(ale_.discret_->Dof(actnode,0))];
    }
    else if (position=="dispy")
    {
      result = (*ale_.dispnp_)[dispnpmap.LID(ale_.discret_->Dof(actnode,1))];
    }
    else if (position=="dispz")
    {
      result = (*ale_.dispnp_)[dispnpmap.LID(ale_.discret_->Dof(actnode,2))];
    }
    else
    {
      dserror("position '%s' not supported in ale testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}

bool FSI::AleResultTest::Match(RESULTDESCR* res)
{
  return res->field==ale;
}

#endif
#endif
#endif
