/*--------------------------------------------------------------------------*/
/*! \file

\brief testing of lubrication calculation results

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef BACI_LUBRICATION_RESULTTEST_H
#define BACI_LUBRICATION_RESULTTEST_H


#include "baci_lib_resulttest.H"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT

namespace LUBRICATION
{
  // forward declaration
  class TimIntImpl;

  /*!
    \brief lubrication specific result test class

    \author wirtz
    \date 11/15
  */
  class ResultTest : public DRT::ResultTest
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
    void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    //! test special quantity not associated with a particular element or node
    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   protected:
    //! get nodal result to be tested
    double ResultNode(const std::string quantity,  //! name of quantity to be tested
        DRT::Node* node                            //! node carrying the result to be tested
    ) const;

    //! get special result to be tested
    virtual double ResultSpecial(const std::string quantity  //! name of quantity to be tested
    ) const;

   private:
    /// Teuchos::RCP to lubrication discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// number of iterations in last newton iteration
    int mynumiter_;
  };
}  // namespace LUBRICATION


BACI_NAMESPACE_CLOSE

#endif  // LUBRICATION_RESULTTEST_H
