/*----------------------------------------------------------------------------*/
/*! \file

\brief Solver for fluid field on a moving ALE mesh

\level 1

*/
/*----------------------------------------------------------------------------*/

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
#include <Teuchos_RCP.hpp>

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
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() override { return fluid_; }

    /// ale field
    const Teuchos::RCP<Adapter::AleFluidWrapper>& ale_field() const { return ale_; }

    /// discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() override;

    /// fluid interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& interface() const override
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

    void nonlinear_solve(Teuchos::RCP<Core::LinAlg::Vector<double>> idisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> ivel) override;

    virtual void nonlinear_solve_vol_coupl(Teuchos::RCP<Core::LinAlg::Vector<double>> idisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> ivel,
        Teuchos::RCP<FSI::InterfaceCorrector> icorrector);

    void apply_interface_values(Teuchos::RCP<Core::LinAlg::Vector<double>> idisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> ivel) override;

    Teuchos::RCP<Core::LinAlg::Vector<double>> relaxation_solve(
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp, double dt) override;

    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_forces() override;
    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_velnp() override;
    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_veln() override;

    int itemax() const override { return fluid_->itemax(); }
    void set_itemax(int itemax) override { fluid_->set_itemax(itemax); }

    Teuchos::RCP<Core::LinAlg::Vector<double>> integrate_interface_shape() override;

    Teuchos::RCP<Core::UTILS::ResultTest> create_field_test() override;

   protected:
    //! @name Transfer helpers
    //@{

    /// field transform
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        Teuchos::RCP<Core::LinAlg::Vector<double>> iv  ///< ALE vector (to be converted)
    ) const;

    /// field transform
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> ale_to_fluid_field(
        Teuchos::RCP<const Core::LinAlg::Vector<double>> iv  ///< ALE vector (to be converted)
    ) const;

    /// interface transform
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> fluid_to_ale(
        Teuchos::RCP<Core::LinAlg::Vector<double>> iv  ///< Fluid vector (to be converted)
    ) const;

    /// interface transform
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> fluid_to_ale(
        Teuchos::RCP<const Core::LinAlg::Vector<double>> iv  ///< Fluid vector (to be converted)
    ) const;

    //@}

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<Coupling::Adapter::CouplingBase> coupfa_;

    /// coupling of fluid and ale (interface or volume...)
    Teuchos::RCP<Coupling::Adapter::CouplingBase> icoupfa_;

    /// coupling of fluid and ale for the ale update condition
    Teuchos::RCP<Coupling::Adapter::Coupling> aucoupfa_;

   private:
    /// problem-specific Fluid-wrapper
    Teuchos::RCP<Adapter::Fluid> fluid_;

    /// problem-specific ALE-wrapper
    Teuchos::RCP<Adapter::AleFluidWrapper> ale_;

    /// problem specific time parameter list
    const Teuchos::ParameterList& timeparams_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
