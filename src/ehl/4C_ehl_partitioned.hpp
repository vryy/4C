// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_EHL_PARTITIONED_HPP
#define FOUR_C_EHL_PARTITIONED_HPP


#include "4C_config.hpp"

#include "4C_ehl_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace EHL
{
  class Partitioned : public Base
  {
   public:
    /// setup EHL algorithm
    explicit Partitioned(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string lubrication_disname);  // Problem builder

    /// setup of single fields (if needed)
    void setup_system() override{};

    /// time loop of coupled problem
    void timeloop() override;

   protected:
    /// prepare time step of single fields
    void prepare_time_step() override;

    //! perform iteration loop between fields
    void outer_loop();

    //! update time step and print to screen
    void update_and_output();

    //! convergence check of outer loop
    bool convergence_check(int itnum);

    /// do one inner iteration loop step of the structure
    void do_struct_step();

    /// do one inner iteration loop step of the lubrication
    void do_lubrication_step();

    //! pressure increment of the outer loop
    std::shared_ptr<Core::LinAlg::Vector<double>> preincnp_;
    //! displacement increment of the outer loop
    std::shared_ptr<Core::LinAlg::Vector<double>> dispincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;
  };
}  // namespace EHL

FOUR_C_NAMESPACE_CLOSE

#endif
