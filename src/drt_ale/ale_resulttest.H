/*----------------------------------------------------------------------------*/
/*!
\file ale_resulttest.H

\brief resutl test for ALE field

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE_RESULTTEST_H
#define ALE_RESULTTEST_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

#include "../drt_adapter/ad_ale.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace DRT
{
  class Discretization;
}

namespace ALE
{
  class Ale;
}

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ALE
{
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
    Teuchos::RCP<const DRT::Discretization> aledis_;

    /// pointer to ALE displacement result vector
    Teuchos::RCP<const Epetra_Vector> dispnp_;
  };

}  // namespace ALE

#endif /* ALE_RESULTTEST_H */
