// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_FLUID_ALE_XFEM_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_ALE_XFEM_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_ale.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{

  /// fluid with moving interfaces implemented by the XFEM
  class FluidAleXFEM : public FluidAle
  {
   public:
    /// constructor
    explicit FluidAleXFEM(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<Core::FE::Discretization> boundary_discretization();

    /// communication object at the struct interface
    virtual Teuchos::RCP<FLD::Utils::MapExtractor> const& struct_interface();

    //@}

    /*========================================================================*/
    //! @name Solver calls
    /*========================================================================*/

    /// nonlinear solve
    void nonlinear_solve(Teuchos::RCP<Core::LinAlg::Vector<double>> idisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> ivel) override;

    /// relaxation solve
    Teuchos::RCP<Core::LinAlg::Vector<double>> relaxation_solve(
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp, double dt) override;
    //@}

    /*========================================================================*/
    //! @name Extract interface forces
    /*========================================================================*/

    /// After the fluid solve we need the forces at the FSI interface.
    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_forces() override;
    //@}

    /*========================================================================*/
    //! @name extract helpers
    /*========================================================================*/

    /// extract the interface velocity at time t^(n+1)
    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_velnp() override;

    /// extract the interface velocity at time t^n
    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_veln() override;
    //@}
    //@}
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
