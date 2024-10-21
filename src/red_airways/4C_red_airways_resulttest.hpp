#ifndef FOUR_C_RED_AIRWAYS_RESULTTEST_HPP
#define FOUR_C_RED_AIRWAYS_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_red_airways_implicitintegration.hpp"
#include "4C_utils_result_test.hpp"

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
  class RedAirwayResultTest : public Core::Utils::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    RedAirwayResultTest(RedAirwayImplicitTimeInt& red_airways);


    void test_node(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

    /*!
      Element test routine, our version of element value tests
    */
    void test_element(
        const Core::IO::InputParameterContainer& container, int& nerr, int& test_count) override;

   private:
    /// Teuchos::RCP to scalar transport discretization
    Teuchos::RCP<Core::FE::Discretization> dis_;
    /// Teuchos::RCP to nodal solution vector containing pressure
    Teuchos::RCP<Core::LinAlg::Vector<double>> mynodesol_pressure_;
    /// Teuchos::RCP to nodal solution vector containing flow in
    Teuchos::RCP<Core::LinAlg::Vector<double>> mynodesol_flow_in_;
    /// Teuchos::RCP to nodal solution vector containing flow out
    Teuchos::RCP<Core::LinAlg::Vector<double>> mynodesol_flow_out_;

    /// Teuchos::RCP to element solution vector containing external pressure of element
    Teuchos::RCP<Core::LinAlg::Vector<double>> myelemsol_pressure_external_;
    /// Teuchos::RCP to element solution vector containing acinus volume
    Teuchos::RCP<Core::LinAlg::Vector<double>> myelemsol_acinivol_;
    /// Teuchos::RCP to element solution vector containing airway volume
    Teuchos::RCP<Core::LinAlg::Vector<double>> myelemsol_airwayvol_;
    /// Teuchos::RCP to element solution vector containing open status of airway
    Teuchos::RCP<Core::LinAlg::Vector<double>> myelemsol_open_;
    /// Teuchos::RCP to element solution vector containing opening trajectory of airway
    Teuchos::RCP<Core::LinAlg::Vector<double>> myelemsol_opening_trajectory_;
  };

}  // namespace Airway

FOUR_C_NAMESPACE_CLOSE

#endif
