/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi (fluidsplit)
using XFEM

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_HPP
#define FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_HPP

#include "4C_config.hpp"

#include "4C_fsi_monolithic_nonox.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Coupling;
  class FluidFluidFSI;
}  // namespace Adapter

namespace Core::LinAlg
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace Core::LinAlg

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
    void SetupSystem() override;

   protected:
    //! @name Apply current field state to system

    //! setup composed right hand side from field solvers
    //!
    //! The RHS consists of three contributions from:
    //! 1) the single fields residuals
    //! 2) the Lagrange multiplier field lambda_
    //! 3) terms in the first nonlinear iteration
    void setup_rhs(Epetra_Vector& f,  ///< empty rhs vector (to be filled)
        bool firstcall = false  ///< indicates whether this is the first nonlinear iteration or not
        ) override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix() override;
    //@}

    /// create merged map of DOF in the final system from all fields
    void create_combined_dof_row_map() override;

    /// initial guess for subsequent fields
    void initial_guess(Teuchos::RCP<Epetra_Vector> ig) override;

    /// create merged Dirichlet map from single field maps
    Teuchos::RCP<Epetra_Map> combined_dbc_map() override;

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
    void extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax) override;

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
    void setup_vector(Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv,
        Teuchos::RCP<const Epetra_Vector> fv, Teuchos::RCP<const Epetra_Vector> av,
        double fluidscale);

    /// access type-cast pointer to problem-specific fluid-wrapper
    const Teuchos::RCP<Adapter::FluidFluidFSI>& fluid_field() { return MonolithicNoNOX::fluid_; }

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<Core::LinAlg::MatrixRowColTransform> fggtransform_;
    Teuchos::RCP<Core::LinAlg::MatrixRowTransform> fgitransform_;
    Teuchos::RCP<Core::LinAlg::MatrixColTransform> figtransform_;
    Teuchos::RCP<Core::LinAlg::MatrixColTransform> aigtransform_;
    Teuchos::RCP<Core::LinAlg::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<Core::LinAlg::MatrixRowColTransform> fmgitransform_;
    Teuchos::RCP<Core::LinAlg::MatrixRowColTransform> fmggtransform_;

    //@}

    /// @name Recovery of Lagrange multiplier at the end of each time step

    //@{
    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! fluid) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fggprev_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fgiprev_;

    //! block \f$F^G_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmggprev_;

    //! block \f$F^G_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fgicur_;

    //! block \f$F^G_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmggcur_;

    //! block \f$F^G_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmgicur_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddialeinc_;
    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;
    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;
    //@}
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
