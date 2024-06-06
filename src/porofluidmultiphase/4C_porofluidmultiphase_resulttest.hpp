/*----------------------------------------------------------------------*/
/*! \file
 \brief result test for multiphase porous fluid

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_RESULTTEST_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_RESULTTEST_HPP



#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::Nodes
{
  class Node;
}

namespace Core::Elements
{
  class Element;
}

namespace POROFLUIDMULTIPHASE
{
  // forward declaration
  class TimIntImpl;

  /*!
    \brief POROFLUIDMULTIPHASE specific result test class

    \author vuong
    \date 08/16
  */
  class ResultTest : public Core::UTILS::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ResultTest(TimIntImpl& porotimint);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "pre"
     */
    void test_node(Input::LineDefinition& res, int& nerr, int& test_count) override;

    /// our version of element value tests
    void TestElement(Input::LineDefinition& res, int& nerr, int& test_count) override;

    //! test special quantity not associated with a particular element or node
    void TestSpecial(Input::LineDefinition& res, int& nerr, int& test_count) override;

   protected:
    //! get nodal result to be tested
    double result_node(const std::string quantity,  //! name of quantity to be tested
        Core::Nodes::Node* node                     //! node carrying the result to be tested
    ) const;

    //! get element result to be tested
    double result_element(const std::string quantity,  //! name of quantity to be tested
        const Core::Elements::Element* element         //! element carrying the result to be tested
    ) const;

    //! get special result to be tested
    virtual double result_special(const std::string quantity  //! name of quantity to be tested
    ) const;

   private:
    //! time integrator
    const TimIntImpl& porotimint_;
  };
}  // namespace POROFLUIDMULTIPHASE


FOUR_C_NAMESPACE_CLOSE

#endif
