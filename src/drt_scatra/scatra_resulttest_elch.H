/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for electrochemistry problems


\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_RESULTTEST_ELCH_H
#define SCATRA_RESULTTEST_ELCH_H

#include "scatra_resulttest.H"

namespace SCATRA
{
  // forward declaration
  class ScaTraTimIntElch;

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
