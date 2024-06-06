/*--------------------------------------------------------------------------*/
/*! \file

\brief testing of lubrication calculation results

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_LUBRICATION_RESULTTEST_HPP
#define FOUR_C_LUBRICATION_RESULTTEST_HPP


#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::Nodes
{
  class Node;
}

namespace LUBRICATION
{
  // forward declaration
  class TimIntImpl;

  /*!
    \brief lubrication specific result test class

    \author wirtz
    \date 11/15
  */
  class ResultTest : public Core::UTILS::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ResultTest(Teuchos::RCP<TimIntImpl> lubrication);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "pre"
     */
    void test_node(Input::LineDefinition& res, int& nerr, int& test_count) override;

    //! test special quantity not associated with a particular element or node
    void TestSpecial(Input::LineDefinition& res, int& nerr, int& test_count) override;

   protected:
    //! get nodal result to be tested
    double result_node(const std::string quantity,  //! name of quantity to be tested
        Core::Nodes::Node* node                     //! node carrying the result to be tested
    ) const;

    //! get special result to be tested
    virtual double result_special(const std::string quantity  //! name of quantity to be tested
    ) const;

   private:
    /// Teuchos::RCP to lubrication discretization
    Teuchos::RCP<Discret::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// number of iterations in last newton iteration
    int mynumiter_;
  };
}  // namespace LUBRICATION


FOUR_C_NAMESPACE_CLOSE

#endif
