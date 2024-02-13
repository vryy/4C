/*----------------------------------------------------------------------*/
/*! \file

\brief testing of electromagnetic calculation results

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ELEMAG_RESULTTEST_HPP
#define BACI_ELEMAG_RESULTTEST_HPP

#include "baci_config.hpp"

#include "baci_lib_resulttest.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace ELEMAG
{
  class ElemagTimeInt;

  class ElemagResultTest : public DRT::ResultTest
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
    void TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to elemagstical discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// Error vector
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> error_;

  };  // class ElemagResultTest

}  // namespace ELEMAG

BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_RESULTTEST_H