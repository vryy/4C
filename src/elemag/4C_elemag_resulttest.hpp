/*----------------------------------------------------------------------*/
/*! \file

\brief testing of electromagnetic calculation results

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_RESULTTEST_HPP
#define FOUR_C_ELEMAG_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace ELEMAG
{
  class ElemagTimeInt;

  class ElemagResultTest : public CORE::UTILS::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    ElemagResultTest(ElemagTimeInt& elemagalgo);


    /// nodal value tests
    /*!
      Possible position flags is only "pressure"
     */
    void test_node(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to elemagstical discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// Error vector
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> error_;

  };  // class ElemagResultTest

}  // namespace ELEMAG

FOUR_C_NAMESPACE_CLOSE

#endif