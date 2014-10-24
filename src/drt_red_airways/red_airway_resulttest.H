/*----------------------------------------------------------------------*/
/*!
\file red_airway_resulttest.H

\brief testing of Red_Airway calculation results

<pre>
Maintainer: Christian Roth
            roth@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef RED_AIRWAY_RESULTTEST_H
#define RED_AIRWAY_RESULTTEST_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_resulttest.H"
#include "airwayimplicitintegration.H"

namespace DRT
{
  class Discretization;
}

namespace AIRWAY{

// forward declaration
class RedAirwayImplicitTimeInt;

/*!
  \brief red_airways specific result test class

  \author Mahmoud Ismail
  \date 04/11
*/
class RedAirwayResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief constructor
  */
  RedAirwayResultTest(RedAirwayImplicitTimeInt& red_airways);


  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  /// Teuchos::RCP to scalar transport discretization
  Teuchos::RCP<DRT::Discretization> dis_;
  /// Teuchos::RCP to solution vector
  Teuchos::RCP<Epetra_Vector>       mysol_;

};

} // namespace AIRWAY

#endif //RED_AIRWAY_RESULTTEST_H

