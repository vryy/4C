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
  ScaTraResultTest(Teuchos::RCP<ScaTraTimIntImpl> scatra);


  /// our version of nodal value tests
  /*!
    Possible position flags is only "phi"
   */
  virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

  virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

private:

  /// Teuchos::RCP to scalar transport discretization
  Teuchos::RCP<DRT::Discretization> dis_;
  /// Teuchos::RCP to solution vector
  Teuchos::RCP<Epetra_Vector>       mysol_;
  /// Teuchos::RCP to solution at nodes (meshfree)
  Teuchos::RCP<Epetra_Vector>       myvan_;
  /// Teuchos::RCP to mass/heat flux vector
  Teuchos::RCP<Epetra_MultiVector>  myflux_;
  /// Teuchos::RCP to electrode kinetics vector
  Teuchos::RCP<std::vector<double> >  myelectkin_;
  /// Teuchos::RCP to scatra structure biofilm growth vector
  Teuchos::RCP<Epetra_MultiVector>  mystrgrowth_;
  /// Teuchos::RCP to scatra fluid biofilm growth vector
  Teuchos::RCP<Epetra_MultiVector>  myfldgrowth_;
  /// number of iterations in last newton iteration
  int mynumiter_;
};

} // namespace SCATRA

#endif
