/*----------------------------------------------------------------------*/
/*! \file

\brief testing of scalar transport calculation results

\level 1


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_RESULTTEST_HPP
#define FOUR_C_SCATRA_RESULTTEST_HPP

#include "baci_config.hpp"

#include "baci_lib_resulttest.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT

namespace SCATRA
{
  // forward declaration
  class ScaTraTimIntImpl;

  /*!
    \brief scalar-transport specific result test class

    \author gjb
    \date 07/08
  */
  class ScaTraResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ScaTraResultTest(Teuchos::RCP<ScaTraTimIntImpl> scatratimint);


    /// our version of nodal value tests
    /*!
      Possible position flags is only "phi"
     */
    void TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    //! test special quantity not associated with a particular element or node
    void TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   protected:
    //! get nodal result to be tested
    virtual double ResultNode(const std::string quantity,  //! name of quantity to be tested
        DRT::Node* node                                    //! node carrying the result to be tested
    ) const;

    //! get special result to be tested
    virtual double ResultSpecial(const std::string quantity  //! name of quantity to be tested
    ) const;

    //! time integrator
    const Teuchos::RCP<const ScaTraTimIntImpl> scatratimint_;

   private:
  };
}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif
