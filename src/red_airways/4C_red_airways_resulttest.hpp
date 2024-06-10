/*---------------------------------------------------------------------*/
/*! \file

\brief Testing of calculation results for reduced elements


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_RED_AIRWAYS_RESULTTEST_HPP
#define FOUR_C_RED_AIRWAYS_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_red_airways_implicitintegration.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Airway
{
  // Forward declaration
  class RedAirwayImplicitTimeInt;

  /*!
    \brief red_airways specific result test class
  */
  class RedAirwayResultTest : public Core::UTILS::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    RedAirwayResultTest(RedAirwayImplicitTimeInt& red_airways);


    void test_node(Input::LineDefinition& res, int& nerr, int& test_count) override;

    /*!
      Element test routine, our version of element value tests
    */
    void TestElement(Input::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to scalar transport discretization
    Teuchos::RCP<Core::FE::Discretization> dis_;
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
    Teuchos::RCP<Epetra_Vector> myelemsol_opening_trajectory_;
  };

}  // namespace Airway

FOUR_C_NAMESPACE_CLOSE

#endif
