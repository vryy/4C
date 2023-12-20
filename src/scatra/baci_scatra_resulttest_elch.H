/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for electrochemistry problems


\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SCATRA_RESULTTEST_ELCH_H
#define BACI_SCATRA_RESULTTEST_ELCH_H

#include "baci_config.H"

#include "baci_scatra_resulttest.H"
#include "baci_scatra_timint_elch.H"

BACI_NAMESPACE_OPEN

namespace SCATRA
{

  // class implementation
  class ElchResultTest : public ScaTraResultTest
  {
   public:
    //! constructor
    ElchResultTest(Teuchos::RCP<ScaTraTimIntElch> elchtimint);

   private:
    //! return pointer to elch time integrator after cast
    Teuchos::RCP<const SCATRA::ScaTraTimIntElch> ElchTimInt() const
    {
      return Teuchos::rcp_dynamic_cast<const SCATRA::ScaTraTimIntElch>(scatratimint_);
    };

    //! get special result to be tested
    double ResultSpecial(const std::string quantity  //! name of quantity to be tested
    ) const override;
  };  // class ElchResultTest : public ScaTraResultTest
}  // namespace SCATRA
BACI_NAMESPACE_CLOSE

#endif
