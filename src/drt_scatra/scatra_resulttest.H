/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.H

\brief tesing of scalar transport calculation results

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_RESULTTEST_H
#define SCATRA_RESULTTEST_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_resulttest.H"

namespace DRT
{
  class Discretization;
}

namespace SCATRA{

// forward declaration
class ScaTraTimIntImpl;

/*!
  \brief scalar-transport specific result test class

  \author gjb
  \date 07/08
*/
class ScaTraResultTest : public DRT::ResultTest
{
public:

  /*!
  \brief constructor
  */
  ScaTraResultTest(ScaTraTimIntImpl& scatra);


  /// our version of nodal value tests
  /*!
    Possible position flags is only "phi"
   */
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  /// RCP to scalar transport discretization
  Teuchos::RCP<DRT::Discretization> dis_;
  /// RCP to solution vector
  Teuchos::RCP<Epetra_Vector>       mysol_;
  /// RCP to mass/heat flux vector
  Teuchos::RCP<Epetra_MultiVector>  myflux_;
  /// RCP to scatra structure biofilm growth vector
  Teuchos::RCP<Epetra_MultiVector>  mystrgrowth_;
  /// RCP to scatra fluid biofilm growth vector
  Teuchos::RCP<Epetra_MultiVector>  myfldgrowth_;
};

} // namespace SCATRA

#endif
