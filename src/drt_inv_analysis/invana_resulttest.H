/*----------------------------------------------------------------------*/
/*!
\file invana_resulttest.H

\brief tesing of inverse analysis results
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

*/


/*----------------------------------------------------------------------*/

#ifndef INVANA_RESULTTEST_H
#define INVANA_RESULTTEST_H


#include "../drt_lib/drt_resulttest.H"
#include "Epetra_MultiVector.h"

namespace DRT
{
  class Discretization;
}

namespace INVANA
{
  // forward declarations
  class InvanaControl;
  class MatParManager;


  /*!
  \brief Invana specific result test class

  \author keh
  */
  class InvanaResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief not documented yet
    */
    InvanaResultTest(InvanaControl& ia);

    virtual ~InvanaResultTest(){};

    /*!
    Possible flag is the name of the material parameter as used by the material, e.g AAAstopro:
    YOUNG or BETA
    */
    virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) { return; }

   private:
    //! pointer to discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! pointer to unknown vector with values
    Teuchos::RCP<Epetra_MultiVector> mysol_;

    //! pointer to manager of the material parameters to be optimized
    Teuchos::RCP<MatParManager> matman_;

    //! the inverse analysis algorithm object
    InvanaControl& ia_;
  };
}  // namespace INVANA


#endif
