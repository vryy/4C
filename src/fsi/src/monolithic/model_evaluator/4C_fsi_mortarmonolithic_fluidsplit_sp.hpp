// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_MORTARMONOLITHIC_FLUIDSPLIT_SP_HPP
#define FOUR_C_FSI_MORTARMONOLITHIC_FLUIDSPLIT_SP_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_inpar_fsi.hpp"

class Epetra_Comm;
namespace NOX
{
  namespace Epetra
  {
    class Group;
  }

  namespace StatusTest
  {
    class Combo;
  }
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Coupling::Adapter
{
  class Coupling;
  class CouplingMortar;
}  // namespace Coupling::Adapter

namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace FSI
{
  class OverlappingBlockMatrix;

  namespace Utils
  {
    class SlideAleUtils;
  }  // namespace Utils
}  // namespace FSI

namespace FSI
{
  class MortarMonolithicFluidSplitSaddlePoint : public BlockMonolithic
  {
    friend class FSI::FSIResultTest;

   public:
    explicit MortarMonolithicFluidSplitSaddlePoint(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    void setup_system() final;

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) final;

    //@}

    /// the composed system matrix
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> system_matrix() const override;

    /// read restart
    void read_restart(int step) final;

    //! @name Time Adaptivity
    //@{

    /*! \brief Select \f$\Delta t_{min}\f$ of all proposed time step sizes
     *         based on error estimation
     *
     *  Depending on the chosen method (fluid or structure split), only 3 of the
     *  6 available norms are useful. Each of these three norms delivers a new
     *  time step size. Select the minimum of these three as the new time step
     *  size.
     */
    double select_dt_error_based() const final;

    /*! \brief Check whether time step is accepted or not
     *
     *  In case that the local truncation error is small enough, the time step
     *  is accepted.
     */
    bool set_accepted() const final;

    //@}

   protected:
    void create_system_matrix();

    Teuchos::RCP<::NOX::StatusTest::Combo> create_status_test(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp) final;

    void update() final;

    void output() final;

    /// Write Lagrange multiplier
    void output_lambda() final;

    /*!
    @copydoc FSI::Monolithic::extract_field_vectors

    Since this is a saddle-point formulation, we also need to extract the Lagrange multipliers
    #lag_mult_ from the monolithic solution vector.
    */
    void extract_field_vectors(std::shared_ptr<const Core::LinAlg::Vector<double>> x,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& fx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& ax,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& lagx);


    /*!
    @copydoc FSI::Monolithic::Evaluate

    Since this is a saddle-point formulation, we also need to evaluate the Lagrange multipliers
    #lag_mult_.
    */
    void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>>
            step_increment  ///< increment between time step n and n+1
        ) override;

    void initial_guess(std::shared_ptr<Core::LinAlg::Vector<double>> initial_guess) final;

   private:
    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid, ALE and Lagrange multipliers to an global FSI
     *  DOF row map.
     */
    void create_combined_dof_row_map() final;

    /*! \brief Create the DOF row map for lagrange
     *
     *  Create the DOF row map for lagrange multiplier based on the fluid interface field and
     *  last GID of the ALE field DOF row map
     */
    virtual void create_lagrange_multiplier_dof_row_map();

    virtual void combine_field_vectors(
        Core::LinAlg::Vector<double>& f,  ///< composed vector containing all field vectors
        std::shared_ptr<const Core::LinAlg::Vector<double>> solid_vector,  ///< structural DOFs
        std::shared_ptr<const Core::LinAlg::Vector<double>> fluid_vector,  ///< fluid DOFs
        std::shared_ptr<const Core::LinAlg::Vector<double>> ale_vector,    ///< ale DOFs
        std::shared_ptr<const Core::LinAlg::Vector<double>>
            lag_mult_vector,  /// < lagrange multiplier
        bool fullvectors);

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a
     *  FSI-global condition map and other map.
     */
    void setup_dbc_map_extractor() final;

    /// setup RHS contributions based on single field residuals
    void setup_rhs_residual(Core::LinAlg::Vector<double>& f) final;

    /// setup RHS contributions based on the Lagrange multiplier field
    void setup_rhs_lambda(Core::LinAlg::Vector<double>& f) final;

    /// setup RHS contributions based on terms for first nonlinear iteration
    void setup_rhs_firstiter(Core::LinAlg::Vector<double>& f) final;

    //! Create #lag_mult_
    virtual void set_lag_mult();

    //! Set #notsetup_ = true after redistribution
    void set_not_setup() override { notsetup_ = true; }

    //! @name Methods for infnorm-scaling of the system
    //!@{

    /// apply infnorm scaling to linear block system
    void scale_system(
        Core::LinAlg::BlockSparseMatrixBase& mat, Core::LinAlg::Vector<double>& b) override;

    /// undo infnorm scaling from scaled solution
    void unscale_solution(Core::LinAlg::BlockSparseMatrixBase& mat, Core::LinAlg::Vector<double>& x,
        Core::LinAlg::Vector<double>& b) override;

    //!@}

    /*! block system matrix
     *  System matrix has a 6x6-block structure corresponding to the vector of unknowns
     *
     *  \f$\Delta x^T = [\Delta d_I^{S,n+1}~\Delta d_\Gamma^{S,n+1}~\Delta u_I^{F,n+1}~\Delta
     * u_\Gamma^{F,n+1}~\Delta d_I^{G,n+1}~\Lambda^{n+1}]\f$.
     *
     * As for the code, we have a 4x4 system.
     */
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> systemmatrix_;

    /// communicator
    const Epetra_Comm& comm_;

    /// @name Matrix block transform objects to handle row and column map exchange for matrix blocks

    /// Coupling of structure and fluid at the interface
    std::shared_ptr<Coupling::Adapter::CouplingMortar> coupling_solid_fluid_mortar_;

    /// Helper variable for the transformation of aleunknowns onto the slave side
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> ale_inner_interf_transform_;

    /// Helper variable for the transformation of fluid unknowns onto the slave side
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> fluid_mesh_inner_inner_transform_;

    ///@}

    /// @name infnorm scaling
    //!@{

    std::shared_ptr<Core::LinAlg::Vector<double>> srowsum_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scolsum_;
    std::shared_ptr<Core::LinAlg::Vector<double>> arowsum_;
    std::shared_ptr<Core::LinAlg::Vector<double>> acolsum_;

    //!@}

    /// additional ale residual to avoid incremental ale errors
    std::shared_ptr<Core::LinAlg::Vector<double>> aleresidual_;

    //! DOF map of Lagrange multiplier unknowns
    std::shared_ptr<const Epetra_Map> lag_mult_dof_map_;

    //! Lagrange multiplier
    std::shared_ptr<Core::LinAlg::Vector<double>> lag_mult_;

    //! Lagrange multiplier from previous time step
    std::shared_ptr<Core::LinAlg::Vector<double>> lag_mult_old_;

    //! Flag to indicate if Setup has not been called yet
    bool notsetup_;

  };  // class MortarMonolithicFluidSplitSaddlePoint
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif