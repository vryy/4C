/*----------------------------------------------------------------------*/
/*!
\file acou_inv_resulttest.H

\brief testing of photoacoustic reconstruction results

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_INV_RESULTTEST_H
#define ACOU_INV_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"
#include <Epetra_Vector.h>

// forward declarations
namespace DRT
{
class Discretization;
}

namespace ACOU
{
class PatImageReconstruction;

class AcouInvResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief constructor
  */
  AcouInvResultTest(PatImageReconstruction& invalgo);


  /// our version of nodal value tests
  /*!
    Possible position flags is only "pressure"
   */
  virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  /// Teuchos::RCP to scatra discretization
  Teuchos::RCP<DRT::Discretization> dis_;

  /// Teuchos::RCP to solution vector
  Teuchos::RCP<Epetra_MultiVector>       mysol_;

}; // class AcouInvResultTest



} // namespace ACOU

#endif /* ACOU_INV_RESULTTEST_H */
