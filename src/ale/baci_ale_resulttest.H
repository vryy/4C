/*----------------------------------------------------------------------------*/
/*! \file

\brief Result tests for pure ALE problems

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef BACI_ALE_RESULTTEST_H
#define BACI_ALE_RESULTTEST_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_adapter_ale.H"
#include "baci_lib_resulttest.H"

#include <Epetra_Vector.h>

BACI_NAMESPACE_OPEN

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
    void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /// pointer to ALE discretization
    Teuchos::RCP<const DRT::Discretization> aledis_;

    /// pointer to ALE displacement result vector
    Teuchos::RCP<const Epetra_Vector> dispnp_;
  };

}  // namespace ALE

BACI_NAMESPACE_CLOSE

#endif  // ALE_RESULTTEST_H
