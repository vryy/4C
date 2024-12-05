// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TURBULENCE_TURBULENT_FLOW_ALGORITHM_HPP
#define FOUR_C_FLUID_TURBULENCE_TURBULENT_FLOW_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_fluid_result_test.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class FluidDiscretExtractor;

  class TurbulentFlowAlgorithm
  {
   public:
    /// constructor
    explicit TurbulentFlowAlgorithm(MPI_Comm comm, const Teuchos::ParameterList& fdyn);

    /// virtual destructor
    virtual ~TurbulentFlowAlgorithm() = default;

    /// do time loop
    void time_loop();

    /// read restart
    /// only during inflow generation
    void read_restart(const int restart);

    /// do result check
    std::shared_ptr<Core::Utils::ResultTest> do_result_check()
    {
      return fluidalgo_->fluid_field()->create_field_test();
    };

   private:
    /// method to transfer inflow velocity from inflow discret to compete discret
    void transfer_inflow_velocity();

    /// discretization of the compete domain
    std::shared_ptr<Core::FE::Discretization> fluiddis_;
    /// discretization of the separate part
    std::shared_ptr<Core::FE::Discretization> inflowdis_;
    /// object for a redistributed evaluation of of the separated part
    std::shared_ptr<FluidDiscretExtractor> inflowgenerator_;
    /// instance of fluid algorithm
    std::shared_ptr<Adapter::FluidBaseAlgorithm> fluidalgo_;
    /// instance of fluid inflow algorithm
    std::shared_ptr<Adapter::FluidBaseAlgorithm> inflowfluidalgo_;
    /// number of time steps
    int step_;
    /// number of development steps
    int numtimesteps_;
    /// velocity/pressure at time n+1 to be transferred to the complete fluid field
    std::shared_ptr<Core::LinAlg::Vector<double>> velnp_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
