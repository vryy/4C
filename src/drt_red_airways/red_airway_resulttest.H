/*---------------------------------------------------------------------*/
/*!

\brief Testing of calculation results for reduced elements

\maintainer Carolin Geitner

\level 3

*/
/*---------------------------------------------------------------------*/

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

namespace AIRWAY
{
  // Forward declaration
  class RedAirwayImplicitTimeInt;

  /*!
    \brief red_airways specific result test class
  */
  class RedAirwayResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    RedAirwayResultTest(RedAirwayImplicitTimeInt& red_airways);


    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    /*!
      Element test routine, our version of element value tests
    */
    virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   private:
    /// Teuchos::RCP to scalar transport discretization
    Teuchos::RCP<DRT::Discretization> dis_;
    /// Teuchos::RCP to nodal solution vector
    Teuchos::RCP<Epetra_Vector> mysol_;
    /// Teuchos::RCP to element solution vector
    Teuchos::RCP<Epetra_Vector> myelemsol_;
  };

}  // namespace AIRWAY

#endif  // RED_AIRWAY_RESULTTEST_H
