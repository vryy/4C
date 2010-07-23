/*----------------------------------------------------------------------*/
/*!
\file ale_resulttest.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "ale_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "ale_lin.H"
#include "ale_springs.H"
#include "ale_laplace.H"
#include "ale_springs_fixed_ref.H"

#ifdef PARALLEL
#include <mpi.h>
#endif


ALE::AleResultTest::AleResultTest(ALE::AleLinear& ale)
: aledis_(ale.Discretization()),
    dispnp_(ale.Disp())
{
}

ALE::AleResultTest::AleResultTest(ALE::AleSprings& ale)
: aledis_(ale.Discretization()),
    dispnp_(ale.Disp())
{
}

ALE::AleResultTest::AleResultTest(ALE::AleLaplace& ale)
: aledis_(ale.Discretization()),
    dispnp_(ale.Disp())
{
}

ALE::AleResultTest::AleResultTest(ALE::AleSpringsFixedRef& ale)
: aledis_(ale.Discretization()),
    dispnp_(ale.Disp())
{
}

void ALE::AleResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one ale discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  if (aledis_->HaveGlobalNode(node))
  {
    DRT::Node* actnode = aledis_->gNode(node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != aledis_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& dispnpmap = dispnp_->Map();

    std::string position;
    res.ExtractString("POSITION",position);
    if (position=="dispx")
    {
      result = (*dispnp_)[dispnpmap.LID(aledis_->Dof(actnode,0))];
    }
    else if (position=="dispy")
    {
      result = (*dispnp_)[dispnpmap.LID(aledis_->Dof(actnode,1))];
    }
    else if (position=="dispz")
    {
      result = (*dispnp_)[dispnpmap.LID(aledis_->Dof(actnode,2))];
    }
    else
    {
      dserror("position '%s' not supported in ale testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}

bool ALE::AleResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("ALE");
}

#endif
