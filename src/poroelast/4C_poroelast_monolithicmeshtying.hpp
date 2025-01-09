// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_MONOLITHICMESHTYING_HPP
#define FOUR_C_POROELAST_MONOLITHICMESHTYING_HPP


#include "4C_config.hpp"

#include "4C_poroelast_monolithic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class CouplingPoroMortar;
}

namespace PoroElast
{
  class MonolithicMeshtying : public Monolithic
  {
   public:
    explicit MonolithicMeshtying(MPI_Comm comm, const Teuchos::ParameterList& timeparams,
        std::shared_ptr<Core::LinAlg::MapExtractor> porosity_splitter);

    //! Setup the monolithic system
    void setup_system() override;

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>>
                      iterinc,  //!< increment between iteration i and i+1
        bool firstiter = false) override;

    //! use monolithic update and set old meshtying quantities at the end of a timestep
    void update() override;

    //! Recover Lagrange Multiplier after Newton step
    void recover_lagrange_multiplier_after_newton_step(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iterinc) override;

    //! build meshtying specific norms where meshtying constraint residuals are evaluated separately
    void build_convergence_norms() override;

    //! extractor to split fluid RHS vector for convergence check
    //! should be named fluidvelocityactiverowdofmap_
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> fluid_vel_active_dof_extractor() const
    {
      return fvelactiverowdofmap_;
    }

    //! setup meshtying activedof extractors
    void setup_extractor();

    //! decide convergence with additional evaluation of meshtying constraint residuals
    bool converged() override;

    //! setup solver with additional residual tolerances for meshtying
    bool setup_solver() override;

    //! contains header to print_newton_iter with meshtying solver tolerance
    void print_newton_iter_header_stream(std::ostringstream& oss) override;

    //! contains text to print_newton_iter with meshtying residuals
    void print_newton_iter_text_stream(std::ostringstream& oss) override;

   private:
    //! nonlinear mortar adapter used to evaluate meshtying
    std::shared_ptr<Adapter::CouplingPoroMortar> mortar_adapter_;

    //! fluid velocity dof row map split in active dofs and the rest (no pressures)
    std::shared_ptr<Core::LinAlg::MultiMapExtractor>
        fvelactiverowdofmap_;  //!< should be named fluidvelocityactiverowdofmap_, but kept shorter

    double normrhsfactiven_;  //!< norm of coupling part of residual forces (fluid )

    double tolfres_ncoup_;  //!< residuum tolerance for porofluid normal coupling condition
  };
}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
