/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for HDG problems

\level 3


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_RESULTTEST_HDG_HPP
#define FOUR_C_SCATRA_RESULTTEST_HDG_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensevector.hpp"
#include "4C_scatra_resulttest.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SCATRA
{
  // forward declaration
  class TimIntHDG;

  // class implementation
  class HDGResultTest : public ScaTraResultTest
  {
   public:
    //! constructor
    HDGResultTest(Teuchos::RCP<ScaTraTimIntImpl> timint);

   private:
    //! get nodal result to be tested
    double ResultNode(const std::string quantity,  //! name of quantity to be tested
        DRT::Node* node                            //! node carrying the result to be tested
    ) const override;

    //! time integrator
    Teuchos::RCP<const TimIntHDG> scatratiminthdg_;

    Teuchos::RCP<CORE::LINALG::SerialDenseVector> errors_;

  };  // class HDGResultTest : public ScaTraResultTest
}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif
