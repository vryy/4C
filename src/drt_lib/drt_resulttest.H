/*----------------------------------------------------------------------*/
/*!

\brief general result test framework

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_RESULTTEST_H
#define DRT_RESULTTEST_H


#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

#include <vector>


namespace DRT
{
  namespace INPUT
  {
    class DatFileReader;
    class LineDefinition;
    class Lines;
  }  // namespace INPUT

  /*!
    \brief Pure virtual base class of all field test classes

    The idea is to have a subclass of this for every algorithm class
    that needs result testing. The Match method needs to be defined to
    state if a particular result value belongs to that field and needs
    to be checked here. And then there are testing methods for element
    tests, nodal tests and special cases (like beltrami fluid
    flow). These methods provide dummy (dserror) implementations and
    have to be redefined in subclasses to actually do the testing.

    \author u.kue
  */
  class ResultTest
  {
   public:
    /// not yet documented
    ResultTest(const std::string name = "NONE");

    /// not yet documented
    virtual ~ResultTest();

    /// perform element value test
    virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// perform nodal value test
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// perform special case test
    virtual void TestSpecial(
        DRT::INPUT::LineDefinition& res, int& nerr, int& test_count, int& unevaluated_test_count);

    /// perform special case test
    virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// tell whether this field test matches to a given line
    virtual bool Match(DRT::INPUT::LineDefinition& res);

   protected:
    //! compare a calculated value with the expected one
    //!
    //! There is a difference between node/element based results and special results.
    //! Node/element based results have to be compared at a specific node/element.
    //! Special results are not attached to a specific node/element, but to the
    //! overall algorithm.
    virtual int CompareValues(double actresult, std::string type, DRT::INPUT::LineDefinition& res);

   private:
    /// specific name of a field test
    const std::string myname_;
  };

  /*!
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
    /// add field specific result test object
    void AddFieldTest(Teuchos::RCP<ResultTest> test);

    /// do all tests of all fields including appropiate output
    void TestAll(const Epetra_Comm& comm);

    /// Definition of valid input section
    Teuchos::RCP<DRT::INPUT::Lines> ValidResultLines();

    /// read the baci input file and set up all curves
    void ReadInput(DRT::INPUT::DatFileReader& reader);

   private:
    /// set of field specific result test objects
    std::vector<Teuchos::RCP<ResultTest>> fieldtest_;

    /// expected results
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> results_;
  };

}  // namespace DRT

void PrintResultDescrDatHeader();

#endif
