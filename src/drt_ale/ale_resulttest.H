/*----------------------------------------------------------------------*/
/*!
\file ale_resulttest.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ALE_RESULTTEST_H
#define ALE_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

namespace DRT
{
  class Discretization;
}

namespace ALE {

// forward declarations
class Ale;

/// Result test subclass for linear ale algorithm
class AleResultTest : public DRT::ResultTest
{
public:
  AleResultTest(ALE::Ale& ale);

  /// our version of nodal value tests
  /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
   */
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  /// pointer to ALE discretization
  Teuchos::RCP<DRT::Discretization> aledis_;

  /// pointer to ALE displacement result vector
  Teuchos::RCP<Epetra_Vector>       dispnp_;

};

}

#endif
