/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for HDG problems

\level 3


*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SCATRA_RESULTTEST_HDG_H
#define BACI_SCATRA_RESULTTEST_HDG_H

#include "baci_linalg_serialdensevector.H"
#include "baci_scatra_resulttest.H"

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
    virtual double ResultNode(const std::string quantity,  //! name of quantity to be tested
        DRT::Node* node                                    //! node carrying the result to be tested
    ) const;

    //! time integrator
    Teuchos::RCP<const TimIntHDG> scatratiminthdg_;

    Teuchos::RCP<CORE::LINALG::SerialDenseVector> errors_;

  };  // class HDGResultTest : public ScaTraResultTest
}  // namespace SCATRA
#endif
