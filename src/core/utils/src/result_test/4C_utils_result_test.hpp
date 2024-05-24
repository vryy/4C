/*----------------------------------------------------------------------*/
/*! \file

\brief general result test framework

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_RESULT_TEST_HPP
#define FOUR_C_UTILS_RESULT_TEST_HPP


#include "4C_config.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace INPUT
{
  class DatFileReader;
  class LineDefinition;
  class Lines;
}  // namespace INPUT

namespace CORE::UTILS
{
  /*!
    \brief Base class of all field test classes

    The idea is to have a subclass of this for every algorithm class
    that needs result testing. The Match method needs to be defined to
    state if a particular result value belongs to that field and needs
    to be checked here. And then there are testing methods for element
    tests, nodal tests and special cases (like beltrami fluid
    flow). These methods provide dummy (FOUR_C_THROW) implementations and
    have to be redefined in subclasses to actually do the testing.

    \author u.kue
  */
  class ResultTest
  {
   public:
    /// not yet documented
    explicit ResultTest(std::string name = "NONE");

    /**
     * Virtual destructor.
     */
    virtual ~ResultTest() = default;

    /// perform element value test
    virtual void TestElement(INPUT::LineDefinition& res, int& nerr, int& test_count);

    /*!
     * @brief  perform nodal value test
     *
     * @param[in] res         input file line containing result test specification
     * @param[in] nerr        number of failed result tests
     * @param[in] test_count  number of result tests
     */
    virtual void test_node(INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// perform special case test
    virtual void TestSpecial(
        INPUT::LineDefinition& res, int& nerr, int& test_count, int& unevaluated_test_count);

    /*!
     * @brief  perform special case test
     *
     * @param[in] res         input file line containing result test specification
     * @param[in] nerr        number of failed result tests
     * @param[in] test_count  number of result tests
     */
    virtual void TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count);

    /// tell whether this field test matches to a given line
    virtual bool Match(INPUT::LineDefinition& res);

   protected:
    //! compare a calculated value with the expected one
    //!
    //! There is a difference between node/element based results and special results.
    //! Node/element based results have to be compared at a specific node/element.
    //! Special results are not attached to a specific node/element, but to the
    //! overall algorithm.
    virtual int CompareValues(double actresult, std::string type, INPUT::LineDefinition& res);

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

    /// Store the parsed @p results.
    void SetParsedLines(std::vector<INPUT::LineDefinition> results);

   private:
    /// set of field specific result test objects
    std::vector<Teuchos::RCP<ResultTest>> fieldtest_;

    /// expected results
    std::vector<INPUT::LineDefinition> results_;
  };

}  // namespace CORE::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
