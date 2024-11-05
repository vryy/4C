// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_DIRICHLETNEUMANNSLIDEALE_HPP
#define FOUR_C_FSI_DIRICHLETNEUMANNSLIDEALE_HPP

#include "4C_config.hpp"

#include "4C_fsi_dirichletneumann.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  namespace Utils
  {
    class SlideAleUtils;
  }
  /**
   *  \brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
   *   with sliding ALE-structure interfaces
   *
   * This class implements the abstract interface FSI::DirichletNeumann
   * for the algorithm class of Dirichlet-Neumann partitioned FSI problems.
   * Specifically, this class enables tangential sliding between the solid and the ALE mesh, while
   * the tangential relative motion between solid and fluid particles adheres to a perfect stick
   * condition.
   *
   * FluidOp() takes a interface displacement, applies it to the ale
   * field, solves the ale field, calculates the interface velocities,
   * applies them to the fluid field, solves the fluid field on the
   * newly deformed fluid mesh and returns the interface forces.
   *
   * StructOp() takes interface forces, applies them to the structural
   * field, solves the field and returns the interface displacements.
   *
   * Furthermore this class contains a Remeshing() method containing the
   * the computation of rotation free ALE displacement values, remeshing
   * of the fluid field and reevaluation of the Mortar interface.
   */
  class DirichletNeumannSlideale : public DirichletNeumann
  {
    friend class DirichletNeumannFactory;

   protected:
    /**
     *  \brief almost empty ctor within saving of important things
     *
     * You will have to use the FSI::DirichletNeumannFactory to create an instance of this class
     *
     * @param[in] comm Communicator
     */
    explicit DirichletNeumannSlideale(const Epetra_Comm& comm);

   public:
    /// setup this object
    void setup() override;

    /**
     *  \brief Perform remeshing to account for mesh sliding at the interface
     *
     * This encompasses computation of rotation free ALE displacement values, remeshing
     * of the fluid field and reevaluation of the Mortar interface terms.
     */
    void remeshing() override;

   protected:
    /** \brief interface fluid operator
     *
     * Solve the fluid field problem.  Since the fluid field is the Dirichlet partition, the
     * interface displacement is prescribed as a Dirichlet boundary condition.
     *
     * \param[in] idisp interface displacement
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface force
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> fluid_op(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp, const FillType fillFlag) final;

    /** \brief interface structural operator
     *
     * Solve the structure field problem.  Since the structure field is the Neumann partition, the
     * interface forces are prescribed as a Neumann boundary condition.
     *
     * \param[in] iforce interface force
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface displacement
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> struct_op(
        std::shared_ptr<Core::LinAlg::Vector<double>> iforce, const FillType fillFlag) final;

    /// predictor
    std::shared_ptr<Core::LinAlg::Vector<double>> initial_guess() override;

   private:
    //! Sliding Ale helper class
    std::shared_ptr<FSI::Utils::SlideAleUtils> slideale_;

    //! Displacement of slave side of the interface
    std::shared_ptr<Core::LinAlg::Vector<double>> islave_;

    //! Slave displacement on master side at every time step begin
    std::shared_ptr<Core::LinAlg::Vector<double>> ft_stemp_;
  };

}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
