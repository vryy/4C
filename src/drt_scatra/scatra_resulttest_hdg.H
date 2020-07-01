/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for HDG problems

\level 3


*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_RESULTTEST_HDG_H
#define SCATRA_RESULTTEST_HDG_H

#include "scatra_resulttest.H"

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

  };  // class HDGResultTest : public ScaTraResultTest
}  // namespace SCATRA
#endif
