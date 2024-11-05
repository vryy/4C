// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_HPP
#define FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_cardiac_monodomain.hpp"
#include "4C_scatra_timint_hdg.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class TimIntCardiacMonodomainHDG : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntHDG
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainHDG(std::shared_ptr<Core::FE::Discretization> dis,
        std::shared_ptr<Core::LinAlg::Solver> solver,
        std::shared_ptr<Teuchos::ParameterList> params,
        std::shared_ptr<Teuchos::ParameterList> sctratimintparams,
        std::shared_ptr<Teuchos::ParameterList> extraparams,
        std::shared_ptr<Core::IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void update() override;

    void collect_runtime_output_data() override;

   protected:
    void element_material_time_update() override;

    void write_restart() const override;

    void collect_problem_specific_runtime_output_data(
        std::shared_ptr<Core::LinAlg::Vector<double>> interpolatedPhi) override;

    //! adapt material
    void pack_material() override;

    //! adapt material
    void unpack_material() override;

    //! project material field
    void project_material() override;

    //! read restart
    void read_restart(
        const int step, std::shared_ptr<Core::IO::InputControl> input = nullptr) override;

   private:
    //! activation time
    std::shared_ptr<Core::LinAlg::Vector<double>> activation_time_interpol_;

    //! element data
    std::shared_ptr<std::vector<char>> data_;


  };  // class TimIntCardiacMonodomainHDG

}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
