// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_FLUID_IMMERSED_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_IMMERSED_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_adapter_fld_moving_boundary.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /// fluid with moving interfaces
  class FluidImmersed : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FluidImmersed(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() override { return fluid_->fluid_field(); }

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() override;

    /// communication object at the interface
    Teuchos::RCP<FLD::Utils::MapExtractor> const& interface() const override;

    //@}

    /*========================================================================*/
    //! @name Time step helpers
    /*========================================================================*/

    /// start new time step
    void prepare_time_step() override;

    /// update at time step end
    void update() override;

    /// output results
    void output() override;

    /// read restart information for given time step
    double read_restart(int step) override;

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

    /*========================================================================*/
    //! @name Number of Newton iterations
    /*========================================================================*/

    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver
    /// Not used for IMMERSED yet !

    /// get the maximum number of iterations from the fluid field
    int itemax() const override { return fluidadapter_->itemax(); }

    /// set the maximum number of iterations for the fluid field
    void set_itemax(int itemax) override { fluid_field()->set_itemax(itemax); }

    /// add dirichlet conditions to dirichlet condmap before next fluid solve
    virtual void add_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    /// remove dirichlet conditions from dirichlet condmap after last fluid solve
    virtual void remove_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    Teuchos::RCP<Core::LinAlg::Vector<double>> integrate_interface_shape() override;

    /// create the testing of fields
    Teuchos::RCP<Core::Utils::ResultTest> create_field_test() override;


   private:
    /// fluid base algorithm object
    Teuchos::RCP<Adapter::FluidBaseAlgorithm> fluid_;
    Teuchos::RCP<Adapter::Fluid> fluidadapter_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
