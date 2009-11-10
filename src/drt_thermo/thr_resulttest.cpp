/*----------------------------------------------------------------------*/
/*!
\file thr_resulttest.cpp

\brief tesing of structure calculation results

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>

#include "thr_resulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
THR::ResultTest::ResultTest(TimInt& tintegrator)
{
  temp_ = tintegrator.Temp();
  rate_ = tintegrator.Rate();
  thrdisc_ = tintegrator.Discretization();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void THR::ResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one thermal discretization supported for testing, but we have %g", dis);

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  // this implementation does not allow testing of heatfluxes
  if (thrdisc_->HaveGlobalNode(node))
  {
    const DRT::Node* actnode = thrdisc_->gNode(node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != thrdisc_->Comm().MyPID())
      return;

    // verbose output
    //cout << "TESTING THERMAL RESULTS with THR::ResultTest::TestNode(..)" << endl;

    std::string position;
    res.ExtractString("POSITION",position);
    bool unknownpos = true;  // make sure the result value string can be handled
    double result = 0.0;  // will hold the actual result of run

    // test temperature
    if (temp_ != Teuchos::null)
    {
      const Epetra_BlockMap& tempmap = temp_->Map();

      if (position == "temp")
      {
        unknownpos = false;
        result = (*temp_)[tempmap.LID(thrdisc_->Dof(actnode,0))];
      }
    }

    // test temperature rates
    if (rate_ != Teuchos::null)
    {
      const Epetra_BlockMap& ratemap = rate_->Map();

      if (position == "rate")
      {
        unknownpos = false;
        result = (*rate_)[ratemap.LID(thrdisc_->Dof(actnode,0))];
      }
    }

    // catch position strings, which are not handled by thermo result test
    if (unknownpos)
      dserror("position '%s' not supported in structure testing", position.c_str());

    // compare values
    const int err = CompareValues(result, res);
    nerr += err;
    test_count++;

    // verbose output
    cout.precision(16);
    cout << "RESULT "  << test_count
         << " IS " << std::scientific << result
         << " AND " << ((err==0) ? "OKAY" : "INCORRECT")
         << endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool THR::ResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  // res.field is a enum of type _FIELDTYP and can be found in headers/enums.h
  return (res.HaveNamed("THERMAL"));
}


#endif /* CCADISCRET */
