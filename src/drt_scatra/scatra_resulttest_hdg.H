/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest_hdg.H

\brief result tests for HDG problems

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15264
</pre>
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
      virtual double ResultNode(
          const std::string   quantity,   //! name of quantity to be tested
          DRT::Node*          node        //! node carrying the result to be tested
          ) const;

      //! time integrator
      Teuchos::RCP<const TimIntHDG> scatratiminthdg_;

  }; // class HDGResultTest : public ScaTraResultTest
} // namespace SCATRA
#endif
