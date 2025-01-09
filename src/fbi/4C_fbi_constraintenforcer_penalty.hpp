// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_CONSTRAINTENFORCER_PENALTY_HPP
#define FOUR_C_FBI_CONSTRAINTENFORCER_PENALTY_HPP

#include "4C_config.hpp"

#include "4C_fbi_constraintenforcer.hpp"
#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace Adapter
namespace Core::Geo
{
  class SearchTree;
}

namespace FBI
{
  class FBIGeometryCoupler;
}

namespace Core::LinAlg
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace Core::LinAlg
namespace Adapter
{
  class ConstraintEnforcerFactory;

  /**
   * \brief Implements the constraint enforcement technique of a penalty approach (Mortar and GPTS)
   * (for fluid-beam interaction)
   *
   * This class decides which information to pass to the actual (partitioned) algorithm.
   * It assembles the stiffness and force contributions of a penalty constraintenforcement approach
   * coupling the velocities of to non-matching meshes
   *
   */
  class FBIPenaltyConstraintenforcer : public FBIConstraintenforcer
  {
    friend ConstraintEnforcerFactory;

   public:
    /**
     * \brief Sets up the constraint enforcer
     *
     *\param[in] structure wrapper for the structure solver
     *\param[in] fluid moving boundary wrapper for the fluid solver
     */
    void setup(std::shared_ptr<Adapter::FSIStructureWrapper> structure,
        std::shared_ptr<Adapter::FluidMovingBoundary> fluid) override;

    /// Handle fbi specific output
    void output(double time, int step) override;

   protected:
    /** \brief You will have to use the Adapter::ConstraintEnforcerFactory
     *
     * \param[in] bridge an object managing the pair contributins
     * \param[in] geometrycoupler an object managing the search, parallel communication, etc.
     */

    FBIPenaltyConstraintenforcer(std::shared_ptr<Adapter::FBIConstraintBridge> bridge,
        std::shared_ptr<FBI::FBIGeometryCoupler> geometrycoupler)
        : FBIConstraintenforcer(bridge, geometrycoupler) {};

    /**
     * \brief Computes the contributions to the stiffness matrix of the structure field.
     *
     *
     * \returns coupling contributions to the structure system matrix
     */
    std::shared_ptr<const Core::LinAlg::SparseMatrix> assemble_structure_coupling_matrix()
        const override;

    /**
     * \brief Computes the contributions to the stiffness matrix of the fluid field.
     *
     * \returns coupling contributions to the fluid system matrix
     */
    std::shared_ptr<const Core::LinAlg::SparseOperator> assemble_fluid_coupling_matrix()
        const override;

    /**
     * \brief Computes the contributions to the rhs of the structure field.
     *
     * \returns coupling contributions to the structure residual
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> assemble_structure_coupling_residual()
        const override;

    /**
     * \brief Computes the contributions to the rhs of the slave field.
     *
     * \returns coupling contributions to the fluid residual
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> assemble_fluid_coupling_residual() const override;

    /// Interface to do preparations to solve the fluid
    void prepare_fluid_solve() override;

    /// Output the constraint violation
    virtual void print_violation(double time, int step);

   private:
    FBIPenaltyConstraintenforcer() = delete;
  };
}  // namespace Adapter
FOUR_C_NAMESPACE_CLOSE

#endif
