// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_RED_AIRWAYS_RESULTTEST_HPP
#define FOUR_C_RED_AIRWAYS_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_red_airways_implicitintegration.hpp"
#include "4C_utils_result_test.hpp"

#include <memory>

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
    /// std::shared_ptr to scalar transport discretization
    std::shared_ptr<Core::FE::Discretization> dis_;
    /// std::shared_ptr to nodal solution vector containing pressure
    std::shared_ptr<Core::LinAlg::Vector<double>> mynodesol_pressure_;
    /// std::shared_ptr to nodal solution vector containing flow in
    std::shared_ptr<Core::LinAlg::Vector<double>> mynodesol_flow_in_;
    /// std::shared_ptr to nodal solution vector containing flow out
    std::shared_ptr<Core::LinAlg::Vector<double>> mynodesol_flow_out_;

    /// std::shared_ptr to element solution vector containing external pressure of element
    std::shared_ptr<Core::LinAlg::Vector<double>> myelemsol_pressure_external_;
    /// std::shared_ptr to element solution vector containing acinus volume
    std::shared_ptr<Core::LinAlg::Vector<double>> myelemsol_acinivol_;
    /// std::shared_ptr to element solution vector containing airway volume
    std::shared_ptr<Core::LinAlg::Vector<double>> myelemsol_airwayvol_;
    /// std::shared_ptr to element solution vector containing open status of airway
    std::shared_ptr<Core::LinAlg::Vector<double>> myelemsol_open_;
    /// std::shared_ptr to element solution vector containing opening trajectory of airway
    std::shared_ptr<Core::LinAlg::Vector<double>> myelemsol_opening_trajectory_;
  };

}  // namespace Airway

FOUR_C_NAMESPACE_CLOSE

#endif
