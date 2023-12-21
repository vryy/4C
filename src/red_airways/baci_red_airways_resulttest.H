/*---------------------------------------------------------------------*/
/*! \file

\brief Testing of calculation results for reduced elements


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_RED_AIRWAYS_RESULTTEST_H
#define BACI_RED_AIRWAYS_RESULTTEST_H

#include "baci_config.H"

#include "baci_lib_resulttest.H"
#include "baci_red_airways_implicitintegration.H"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

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


    void TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    /*!
      Element test routine, our version of element value tests
    */
    void TestElement(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

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

BACI_NAMESPACE_CLOSE

#endif  // RED_AIRWAYS_RESULTTEST_H
