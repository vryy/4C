/*----------------------------------------------------------------------*/
/*!
\file drt_resulttest.H

\brief general result test framework

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_RESULTTEST_H
#define DRT_RESULTTEST_H

#ifdef CCADISCRET

#include <Epetra_Comm.h>
#include <Teuchos_RefCountPtr.hpp>

#include <vector>

struct _RESULTDESCR;

namespace DRT
{

  /*!
    \class ResultTest

    \brief Pure virtual base class of all field test classes

    The idea is to have a subclass of this for every algorithm class
    that needs result testing. The Match method needs to be defined to
    state if a particular result value belongs to that field and needs
    to be checked here. And then there are testing methods for element
    tests, nodal tests and special cases (like beltrami fluid
    flow). These methods provide dummy (dserror) implementations and
    have to be redefined in subclasses to actually do the testing.

    The basis of the testing mechanism is the old C RESULTDESCR
    structure that contains one line of the RESULT DESCRIPTION section
    of our input files. Such a line states which internal baci
    variable is to be checked against what value with which
    tolerance. Pretty simple actually. But the meaning of such a line
    is algorithm specific. That's why there are subclasses needed.

    \author u.kue
  */
  class ResultTest
  {
  public:

    /// not yet doxumented
    ResultTest();

    /// not yet doxumented
    virtual ~ResultTest();

    /// perform element value test
    virtual void TestElement(_RESULTDESCR* res, int& nerr, int& test_count);

    /// perform nodal value test
    virtual void TestNode(_RESULTDESCR* res, int& nerr, int& test_count);

    /// perform special case test
    virtual void TestSpecial(_RESULTDESCR* res, int& nerr, int& test_count);

    /// tell whether this field test matches to a given line
    virtual bool Match(_RESULTDESCR* res) = 0;

  protected:

    /// compare a calculated value with the expected one
    virtual int CompareValues(double actresult, _RESULTDESCR *res);
  };

  /*!
    \class ResultTestManager

    \brief Manager class of result test framework

    You have to create one object of this class to test the results of
    your calculation. For each field involved you will want to add a
    specific field test class (derived from ResultTest). Afterwards
    just start testing...

    \author u.kue
  */
  class ResultTestManager
  {
  public:

    /// not yet doxumented
    ResultTestManager(const Epetra_Comm& comm);

    /// not yet doxumented
    virtual ~ResultTestManager();

    /// add field specific result test object
    void AddFieldTest(Teuchos::RefCountPtr<ResultTest> test);

    /// do all tests of all fields including appropiate output
    void TestAll();

  private:

    /// our communicator
    const Epetra_Comm& comm_;

    /// set of field specific result test objects
    std::vector<Teuchos::RefCountPtr<ResultTest> > fieldtest_;
  };

}

#endif /* CCADISCRET       */

#endif
