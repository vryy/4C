// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FPSI_MONOLITHIC_PLAIN_HPP
#define FOUR_C_FPSI_MONOLITHIC_PLAIN_HPP

#include "4C_config.hpp"

#include "4C_fpsi_monolithic.hpp"
#include "4C_inpar_fpsi.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Adapter
{
  class Structure;
  class Fluid;
}  // namespace Adapter


namespace FSI
{
  class MonolithicFluidSplit;
}  // namespace FSI

namespace Core::LinAlg
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace Core::LinAlg

namespace FPSI
{
  class MonolithicPlain : public Monolithic
  {
   public:
    explicit MonolithicPlain(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    /*! do the setup for the monolithic system


    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix


    */
    void setup_system() override;

    //! @name Apply current field state to system


    //@}

    /// the composed system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> system_matrix() const override
    {
      return systemmatrix_;
    }

    /// setup composed right hand side from field solvers
    void setup_rhs(bool firstcall = false) override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

    /// Recover the Lagrange multiplier at the interface   mayr.mt (03/2012)
    void recover_lagrange_multiplier() override;


   protected:
    //! set full monolithic dof row map
    /*!
    A subclass calls this method (from its constructor) and thereby
    defines the number of blocks, their maps and the block order. The block
    maps must be row maps by themselves and must not contain identical GIDs.
    */
    void set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// extract the three field vectors from a given composed vector
    void extract_field_vectors(Teuchos::RCP<const Core::LinAlg::Vector<double>>
                                   x,  ///< composed vector that contains all field vectors
        Teuchos::RCP<const Core::LinAlg::Vector<double>>& sx,   ///< structural displacements
        Teuchos::RCP<const Core::LinAlg::Vector<double>>& pfx,  ///< fluid velocities and pressure
        Teuchos::RCP<const Core::LinAlg::Vector<double>>& fx,   ///< fluid velocities and pressure
        Teuchos::RCP<const Core::LinAlg::Vector<double>>& ax,   ///< ale displacements
        bool firstiter_                                         ///< firstiteration?
        ) override;

   private:
    /// build block vector from field vectors
    void setup_vector(
        Core::LinAlg::Vector<double>& f,  ///< composed vector that contains all field vectors
        const Core::LinAlg::Vector<double>& sv,   ///< structural dofs
        const Core::LinAlg::Vector<double>& pfv,  ///< poro fluid velocities and pressure
        const Core::LinAlg::Vector<double>& fv,   ///< fluid velocities and pressure
        const Core::LinAlg::Vector<double>& av,   ///< ale displacements
        double fluidscale                         ///< residual scaling for fluid
    );

    void setup_rhs_lambda(
        Core::LinAlg::Vector<double>& f  ///< composed vector that contains all field vectors
    );

    void setup_rhs_first_iter(
        Core::LinAlg::Vector<double>& f  ///< composed vector that contains all field vectors
    );

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform>
        fggtransform_;  // g_fsi/g_fsi || F->S/F->S transform (FSI/FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform>
        fggtransform2_;  // g_fis/g_fpsi || F->S/F->S transform (FSI/FPSI)
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform>
        fmgitransform_;  // g_fsi/i transform || F->S/F->A transform (FSI/FluidVolume)

    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform>
        fgitransform1_;  // g_fsi || F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform>
        fgitransform2_;  // g_fsi || F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform>
        cfgtransform_;  // full fluid_field || F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform>
        cfptransform_;  // full fluid_field || F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform>
        cfptransform2_;  // full fluid_field || F->S transform (FSI)

    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        figtransform1_;  // for Row/Col-Map for Full - fluid_field & F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        figtransform2_;  // for Row/Col-Map for Full - fluid_field & F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        figtransform3_;  // for Row/Col-Map for Full - fluid_field & F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        figtransform4_;  // for Row/Col-Map for Full - fluid_field & F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        aigtransform_;  // for Row/Col-Map for Full - ale_field & F->S transform (FSI)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        aigtransform2_;  // for Row/Col-Map for Full - ale_field & F->S transform (FPSI)

    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        couplingcoltransform_;  // for Row/Col-Map for Full - fluid_field & F->A transform
                                // (FluidVolume)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform>
        couplingcoltransformfs_;  // for Row/Col-Map for Full - fluid_field & F->S transform (FPSI)
    ///@}


    /// @name infnorm scaling

    Teuchos::RCP<Core::LinAlg::Vector<double>> srowsum_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scolsum_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> arowsum_;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acolsum_;

    //@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Core::LinAlg::Vector<double>> aleresidual_;


    /// @name Recovery of Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! fluid) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{F,n+1}\f$ onto the fluid at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::Vector<double>> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{F,n+1}\f$ onto the fluid at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::Vector<double>> fgprev_;

    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current NOX iteration \f$i+1\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> ddginc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> duiinc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Core::LinAlg::Vector<double>> disgprev_;

    //! inner velocity solution of fluid at previous NOX iteration
    Teuchos::RCP<const Core::LinAlg::Vector<double>> soliprev_;

    //! interface velocity solution of the fluid at previous NOX iteration
    Teuchos::RCP<const Core::LinAlg::Vector<double>> solgprev_;

    //! inner ALE displacement solution at previous NOX iteration
    Teuchos::RCP<const Core::LinAlg::Vector<double>> solialeprev_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> ddialeinc_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fgicur_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fggprev_;

    //! block \f$F_{\Gamma I,i+1}^G\f$ of fluid shape derivatives matrix at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmgicur_;

    //! block \f$F_{\Gamma I,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration
    //! \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}^G\f$ of fluid shape derivatives matrix at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmggcur_;

    //! block \f$F_{\Gamma\Gamma,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration
    //! \f$i\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> fmggprev_;

    //@}

  };  // class monolithic_plain
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
