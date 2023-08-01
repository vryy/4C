/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for electrochemistry problems


\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SCATRA_RESULTTEST_ELCH_H
#define BACI_SCATRA_RESULTTEST_ELCH_H

#include "baci_scatra_resulttest.H"
#include "baci_scatra_timint_elch.H"

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
    ) const;
  };  // class ElchResultTest : public ScaTraResultTest
}  // namespace SCATRA
#endif
