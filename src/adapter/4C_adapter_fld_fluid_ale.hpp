// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_FLUID_ALE_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_ALE_HPP

/*----------------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_adapter_fld_moving_boundary.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_Comm.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace Adapter
{
  class AleFluidWrapper;
}  // namespace Adapter

namespace FSI
{
  class InterfaceCorrector;
}

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace Adapter
{
  /// fluid on ale
  class FluidAle : public FluidMovingBoundary
  {
   public:
    FluidAle(const Teuchos::ParameterList& prbdyn, std::string condname);

    /// fluid field
    const std::shared_ptr<Adapter::Fluid>& fluid_field() override { return fluid_; }

    /// ale field
    const std::shared_ptr<Adapter::AleFluidWrapper>& ale_field() const { return ale_; }

    /// discretization
    std::shared_ptr<Core::FE::Discretization> discretization() override;

    /// fluid interface
    std::shared_ptr<FLD::Utils::MapExtractor> const& interface() const override
    {
      return fluid_->interface();
    }

    /// Prepare a single time step
    void prepare_time_step() override;

    /// Update to go from time step \f$t_n\f$ to \f$t_{n+1}\f$
    void update() override;

    /// Output current state of simulation
    void output() override;

    /// Read resatart data
    double read_restart(int step  ///< step number to restart from
        ) override;

    void nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel) override;

    virtual void nonlinear_solve_vol_coupl(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel,
        std::shared_ptr<FSI::InterfaceCorrector> icorrector);

    void apply_interface_values(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> ivel) override;

    std::shared_ptr<Core::LinAlg::Vector<double>> relaxation_solve(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp, double dt) override;

    std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_forces() override;
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_velnp() override;
    std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_veln() override;

    int itemax() const override { return fluid_->itemax(); }
    void set_itemax(int itemax) override { fluid_->set_itemax(itemax); }

    std::shared_ptr<Core::LinAlg::Vector<double>> integrate_interface_shape() override;

    std::shared_ptr<Core::Utils::ResultTest> create_field_test() override;

   protected:
    //! @name Transfer helpers
    //@{

    /// field transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv  ///< ALE vector (to be converted)
    ) const;

    /// field transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv  ///< ALE vector (to be converted)
    ) const;

    /// interface transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_ale(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv  ///< Fluid vector (to be converted)
    ) const;

    /// interface transform
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_ale(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv  ///< Fluid vector (to be converted)
    ) const;

    //@}

    /// coupling of fluid and ale (whole field)
    std::shared_ptr<Coupling::Adapter::CouplingBase> coupfa_;

    /// coupling of fluid and ale (interface or volume...)
    std::shared_ptr<Coupling::Adapter::CouplingBase> icoupfa_;

    /// coupling of fluid and ale for the ale update condition
    std::shared_ptr<Coupling::Adapter::Coupling> aucoupfa_;

   private:
    /// problem-specific Fluid-wrapper
    std::shared_ptr<Adapter::Fluid> fluid_;

    /// problem-specific ALE-wrapper
    std::shared_ptr<Adapter::AleFluidWrapper> ale_;

    /// problem specific time parameter list
    const Teuchos::ParameterList& timeparams_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
