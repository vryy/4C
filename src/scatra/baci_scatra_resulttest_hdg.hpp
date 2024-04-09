/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for HDG problems

\level 3


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_RESULTTEST_HDG_HPP
#define FOUR_C_SCATRA_RESULTTEST_HDG_HPP

#include "baci_config.hpp"

#include "baci_linalg_serialdensevector.hpp"
#include "baci_scatra_resulttest.hpp"

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif
