/*---------------------------------------------------------------------*/
/*! \file

\brief Testing of calculation results for reduced elements


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
    /// Teuchos::RCP to nodal solution vector containing pressure
    Teuchos::RCP<Epetra_Vector> mynodesol_pressure_;
    /// Teuchos::RCP to nodal solution vector containing flow in
    Teuchos::RCP<Epetra_Vector> mynodesol_flow_in_;
    /// Teuchos::RCP to nodal solution vector containing flow out
    Teuchos::RCP<Epetra_Vector> mynodesol_flow_out_;

    /// Teuchos::RCP to element solution vector containing external pressure of element
    Teuchos::RCP<Epetra_Vector> myelemsol_pressure_external_;
    /// Teuchos::RCP to element solution vector containing acinus volume
    Teuchos::RCP<Epetra_Vector> myelemsol_acinivol_;
    /// Teuchos::RCP to element solution vector containing airway volume
    Teuchos::RCP<Epetra_Vector> myelemsol_airwayvol_;
    /// Teuchos::RCP to element solution vector containing open status of airway
    Teuchos::RCP<Epetra_Vector> myelemsol_open_;
    /// Teuchos::RCP to element solution vector containing opening trajectory of airway
    Teuchos::RCP<Epetra_Vector> myelemsol_openingTrajectory_;
  };

}  // namespace AIRWAY

#endif  // RED_AIRWAY_RESULTTEST_H
