// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_HPP
#define FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_fsi_monolithic_nonox.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FluidFluidFSI;
}  // namespace Adapter

namespace FSI
{
  /// Monolithic fluid-fluid FSI algorithm (fluidsplit)
  /*!
   * Split of the fluid matrix, while the structure matrix remains
   * unmodified / condensation of interface fluid velocities.
   */
  class FluidFluidMonolithicFluidSplitNoNOX : public MonolithicNoNOX
  {
    friend class FSI::FSIResultTest;

   public:
    /// constructor
    explicit FluidFluidMonolithicFluidSplitNoNOX(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// initialize, read parameters and create merged DOF row map
    void setup_system() override;

   protected:
    //! @name Apply current field state to system

    //! setup composed right hand side from field solvers
    //!
    //! The RHS consists of three contributions from:
    //! 1) the single fields residuals
    //! 2) the Lagrange multiplier field lambda_
    //! 3) terms in the first nonlinear iteration
    void setup_rhs(Core::LinAlg::Vector<double>& f,  ///< empty rhs vector (to be filled)
        bool firstcall = false  ///< indicates whether this is the first nonlinear iteration or not
        ) override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix() override;
    //@}

    /// create merged map of DOF in the final system from all fields
    void create_combined_dof_row_map() override;

    /// initial guess for subsequent fields
    void initial_guess(std::shared_ptr<Core::LinAlg::Vector<double>> ig) override;

    /// create merged Dirichlet map from single field maps
    std::shared_ptr<Epetra_Map> combined_dbc_map() override;

    /// Newton loop
    void newton() override;

    //! Extract the three field vectors from a given composed vector
    //!
    //! In analogy to NOX, x is step increment \f$\Delta x\f$
    //! that brings us from \f$t^{n}\f$ to \f$t^{n+1}\f$:
    //! \f$x^{n+1} = x^{n} + \Delta x\f$
    //!
    //! Iteration increments, that are needed internally in the single fields,
    //! have to be computed somewhere else.
    //!
    //! \param x  (i) composed vector that contains all field vectors
    //! \param sx (o) structural displacements
    //! \param fx (o) fluid velocities and pressure
    //! \param ax (o) ale displacements
    void extract_field_vectors(std::shared_ptr<const Core::LinAlg::Vector<double>> x,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& fx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& ax) override;

    /// compute the Lagrange multiplier (FSI stresses) for the current time step
    void recover_lagrange_multiplier() override;

    /// compute the residual and incremental norms required for convergence check
    void build_convergence_norms() override;

    /// read restart data
    void read_restart(int step) override;

    /// output of fluid, structure & ALE-quantities and Lagrange multiplier
    void output() override;

    /*!
     * In case of a change in the fluid DOF row maps during the Newton loop (full Newton approach),
     * reset vectors accordingly.
     * \author kruse
     * \date 05/14
     */
    void handle_fluid_dof_map_change_in_newton() override;

    /*!
     * Determine a change in fluid DOF map
     * \param (in) : DOF map of fluid increment vector
     * \return : true, in case of a mismatch between map of increment vector
     * and inner fluid DOF map after evaluation
     * \author kruse
     * \date 05/14
     */
    bool has_fluid_dof_map_changed(const Epetra_BlockMap& fluidincrementmap) override;

   private:
    /// build block vector from field vectors
    void setup_vector(Core::LinAlg::Vector<double>& f, const Core::LinAlg::Vector<double>& sv,
        const Core::LinAlg::Vector<double>& fv, const Core::LinAlg::Vector<double>& av,
        double fluidscale);

    /// access type-cast pointer to problem-specific fluid-wrapper
    const std::shared_ptr<Adapter::FluidFluidFSI>& fluid_field() { return MonolithicNoNOX::fluid_; }

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    std::shared_ptr<Coupling::Adapter::MatrixRowColTransform> fggtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> fgitransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> figtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> aigtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> fmiitransform_;
    std::shared_ptr<Coupling::Adapter::MatrixRowColTransform> fmgitransform_;
    std::shared_ptr<Coupling::Adapter::MatrixRowColTransform> fmggtransform_;

    //@}

    /// @name Recovery of Lagrange multiplier at the end of each time step

    //@{
    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! fluid) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> lambda_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fggprev_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fgiprev_;

    //! block \f$F^G_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fmggprev_;

    //! block \f$F^G_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fmgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fgicur_;

    //! block \f$F^G_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fmggcur_;

    //! block \f$F^G_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fmgicur_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX
    //! iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> ddialeinc_;
    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> ddginc_;
    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration
    //! \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> duiinc_;
    //@}
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
