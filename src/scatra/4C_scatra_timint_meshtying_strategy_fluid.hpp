// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_HPP
#define FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_HPP

#include "4C_config.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_scatra_timint_meshtying_strategy_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace FLD
{
  class Meshtying;
}

namespace ScaTra
{
  /*!
  \brief Fluid-fluid meshtying strategy for standard scalar transport problems

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the fluid-fluid meshtying strategy for standard
  scalar transport problems.

  */

  class MeshtyingStrategyFluid : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyFluid(ScaTra::ScaTraTimIntImpl* scatratimint);

    //! return global map of degrees of freedom
    const Epetra_Map& dof_row_map() const override;

    //! compute meshtying residual terms and their linearizations
    void evaluate_meshtying() override;

    //! include Dirichlet conditions into condensation
    void include_dirichlet_in_condensation() const override;

    //! initialize meshtying objects
    void init_meshtying() override;

    bool system_matrix_initialization_needed() const override { return true; }

    std::shared_ptr<Core::LinAlg::SparseOperator> init_system_matrix() const override;

    std::shared_ptr<Core::LinAlg::MultiMapExtractor> interface_maps() const override
    {
      FOUR_C_THROW("InterfaceMaps() is not implemented in MeshtyingStrategyFluid.");
      return nullptr;
    }

    //! setup meshtying objects
    void setup_meshtying() override;

    //! solve resulting linear system of equations
    void solve(const std::shared_ptr<Core::LinAlg::Solver>& solver,         //!< solver
        const std::shared_ptr<Core::LinAlg::SparseOperator>& systemmatrix,  //!< system matrix
        const std::shared_ptr<Core::LinAlg::Vector<double>>& increment,     //!< increment vector
        const std::shared_ptr<Core::LinAlg::Vector<double>>& residual,      //!< residual vector
        const std::shared_ptr<Core::LinAlg::Vector<double>>& phinp,  //!< state vector at time n+1
        const int iteration,  //!< number of current Newton-Raphson iteration
        Core::LinAlg::SolverParams& solver_params) const override;

    //! return linear solver for global system of linear equations
    const Core::LinAlg::Solver& solver() const override;

   protected:
    //! instantiate strategy for Newton-Raphson convergence check
    void init_conv_check_strategy() override;

    //! fluid-fluid meshtying algorithm for internal interface
    std::shared_ptr<FLD::Meshtying> meshtying_;

    //! type of fluid-fluid meshtying
    enum Inpar::FLUID::MeshTying type_;

   private:
    //! copy constructor
    MeshtyingStrategyFluid(const MeshtyingStrategyFluid& old);
  };  // class MeshtyingStrategyFluid
}  // namespace ScaTra
FOUR_C_NAMESPACE_CLOSE

#endif
