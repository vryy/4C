/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest_elch.H

\brief result tests for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
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

      //! get result to be tested
      const double GetResult(
          const std::string   quantity,   //! name of quantity to be tested
          DRT::Node*          node        //! node carrying the result to be tested
          ) const;

      //! time integrator
      const Teuchos::RCP<ScaTraTimIntElch>& elchtimint_;
  }; // class ElchResultTest : public ScaTraResultTest
} // namespace SCATRA
#endif
