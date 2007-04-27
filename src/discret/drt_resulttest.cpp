
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef RESULTTEST

#include "drt_resulttest.H"

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
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*!
 * \brief An array that contains the expected results.
 *
 * An array that contains the expected results.
 * Along with the real values the position in the mesh ist stored. The
 * data is read from the input file.
 *
 * \author uk
 * \date 06/04
 */
struct _RESULTDESCR      *resultdescr;
}

using namespace DRT;

ResultTest::ResultTest()
{
}

ResultTest::~ResultTest()
{
}

void ResultTest::TestElement(RESULTDESCR* res, int& nerr, int& test_count)
{
  dserror("no element test availabe");
}

void ResultTest::TestNode(RESULTDESCR* res, int& nerr, int& test_count)
{
  dserror("no node test availabe");
}

void ResultTest::TestSpecial(RESULTDESCR* res, int& nerr, int& test_count)
{
  dserror("no special case test availabe");
}


/*----------------------------------------------------------------------*/
/*!
 \brief Compare \a actresult with \a givenresult and return 0 if they are
 considered to be equal.

 Compare \a actresult with \a givenresult and return 0 if they are
 considered to be equal.

 \param err        the file where to document both values
 \param res        the describtion of the expected result including name and tolerance

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
int ResultTest::CompareValues(double actresult, RESULTDESCR *res)
{
  FILE *err = allfiles.out_err;
  int ret = 0;
  double givenresult = res->value;

  fprintf(err,"actual = %24.16f, given = %24.16f, diff = %24.16f\n",
          actresult, givenresult, actresult-givenresult);
  if (!(FABS(FABS(actresult-givenresult)-FABS(actresult-givenresult)) < res->tolerance) )
  {
    printf("RESULTCHECK: %s is NAN!\n", res->name);
    ret = 1;
  }
  else if (FABS(actresult-givenresult) > res->tolerance)
  {
    printf("RESULTCHECK: %s not correct. actresult=%15.9e, givenresult=%15.9e\n",
           res->name, actresult, givenresult);
    ret = 1;
  }
  return ret;
}


///////////////////////////////////////////////////////////////////


ResultTestManager::ResultTestManager(const Epetra_Comm& comm)
  : comm_(comm)
{
}

ResultTestManager::~ResultTestManager()
{
}


void ResultTestManager::AddFieldTest(Teuchos::RefCountPtr<ResultTest> test)
{
  fieldtest_.push_back(test);
}


void ResultTestManager::TestAll()
{
  FILE *err = allfiles.out_err;
  INT nerr = 0;
  INT test_count = 0;

  if (comm_.MyPID()==0)
    cout << "\n" GRAY_LIGHT "Checking results ..." END_COLOR "\n";

  for (int i=0; i<genprob.numresults; ++i)
  {
    RESULTDESCR* res = &(resultdescr[i]);

    for (unsigned j=0; j<fieldtest_.size(); ++j)
    {
      if (fieldtest_[j]->Match(res))
      {
        if (res->element != -1)
          fieldtest_[j]->TestElement(res,nerr,test_count);
        else if (res->node != -1)
          fieldtest_[j]->TestNode(res,nerr,test_count);
        else
          fieldtest_[j]->TestSpecial(res,nerr,test_count);
      }
    }
  }

  int numerr;
  comm_.SumAll(&nerr,&numerr,1);

  if (numerr > 0)
  {
    dserror("Result check failed with %d errors out of %d tests", numerr, genprob.numresults);
  }
  else
    fprintf(err,"===========================================\n");

  /* test_count == -1 means we had a special test routine. It's thus
   * illegal to use both a special routine and single tests. But who
   * wants that? */
  if (test_count > -1)
  {
    int count;
    comm_.SumAll(&test_count,&count,1);

    /* It's indeed possible to count more tests than expected if you
     * happen to test values of a boundary element. We don't want this
     * dserror to go off in that case. */
    if (count < genprob.numresults)
    {
      dserror("expected %d tests but performed %d", genprob.numresults, count);
    }
  }

  cout << "\n" MAGENTA_LIGHT "OK" END_COLOR "\n";
}


#endif
#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */
