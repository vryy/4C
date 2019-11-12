/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FPSI problem with matching grids using a monolithic scheme
       in its plain form without any fancy splitting and condensation stuff

\level 3

\maintainer  Johannes Kremheller
*/

/*----------------------------------------------------------------------*/

#ifndef FPSI_MONOLITHICPLAIN_H
#define FPSI_MONOLITHICPLAIN_H

#include "fpsi_monolithic.H"
#include "../drt_inpar/inpar_fpsi.H"


namespace ADAPTER
{
  class Structure;
  class Fluid;
}  // namespace ADAPTER


namespace FSI
{
  class MonolithicFluidSplit;
}  // namespace FSI

namespace LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace LINALG

namespace FPSI
{
  class Monolithic_Plain : public Monolithic
  {
   public:
    explicit Monolithic_Plain(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    /*! do the setup for the monolithic system


    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix


    */
    virtual void SetupSystem();

    //! @name Apply current field state to system


    //@}

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    /// Recover the Lagrange multiplier at the interface   mayr.mt (03/2012)
    virtual void RecoverLagrangeMultiplier();


   protected:
    //! set full monolithic dof row map
    /*!
    A subclass calls this method (from its constructor) and thereby
    defines the number of blocks, their maps and the block order. The block
    maps must be row maps by themselves and must not contain identical GIDs.
    */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// extract the three field vectors from a given composed vector
    void ExtractFieldVectors(
        Teuchos::RCP<const Epetra_Vector> x,    ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector>& sx,  ///< structural displacements
        Teuchos::RCP<const Epetra_Vector>& pfx,  ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector>& fx,   ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector>& ax,   ///< ale displacements
        bool firstiter_                          ///< firstiteration?
    );

   private:
    /// build block vector from field vectors
    void SetupVector(Epetra_Vector& f,          ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector> sv,   ///< structural dofs
        Teuchos::RCP<const Epetra_Vector> pfv,  ///< poro fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector> av,   ///< ale displacements
        double fluidscale                       ///< residual scaling for fluid
    );

    void SetupRHSLambda(Epetra_Vector& f  ///< composed vector that contains all field vectors
    );

    void SetupRHSFirstIter(Epetra_Vector& f  ///< composed vector that contains all field vectors
    );

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<LINALG::MatrixRowColTransform>
        fggtransform_;  // g_fsi/g_fsi || F->S/F->S transform (FSI/FSI)
    Teuchos::RCP<LINALG::MatrixRowColTransform>
        fggtransform2_;  // g_fis/g_fpsi || F->S/F->S transform (FSI/FPSI)
    Teuchos::RCP<LINALG::MatrixRowColTransform>
        fmgitransform_;  // g_fsi/i transform || F->S/F->A transform (FSI/FluidVolume)

    Teuchos::RCP<LINALG::MatrixRowTransform> fgitransform1_;  // g_fsi || F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixRowTransform> fgitransform2_;  // g_fsi || F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixRowTransform>
        Cfgtransform_;  // full FluidField || F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixRowTransform>
        Cfptransform_;  // full FluidField || F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixRowTransform>
        Cfptransform2_;  // full FluidField || F->S transform (FSI)

    Teuchos::RCP<LINALG::MatrixColTransform>
        figtransform1_;  // for Row/Col-Map for Full - FluidField & F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixColTransform>
        figtransform2_;  // for Row/Col-Map for Full - FluidField & F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixColTransform>
        figtransform3_;  // for Row/Col-Map for Full - FluidField & F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixColTransform>
        figtransform4_;  // for Row/Col-Map for Full - FluidField & F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixColTransform>
        aigtransform_;  // for Row/Col-Map for Full - AleField & F->S transform (FSI)
    Teuchos::RCP<LINALG::MatrixColTransform>
        aigtransform2_;  // for Row/Col-Map for Full - AleField & F->S transform (FPSI)

    Teuchos::RCP<LINALG::MatrixColTransform>
        couplingcoltransform_;  // for Row/Col-Map for Full - FluidField & F->A transform
                                // (FluidVolume)
    Teuchos::RCP<LINALG::MatrixColTransform>
        couplingcoltransformfs_;  // for Row/Col-Map for Full - FluidField & F->S transform (FPSI)
    ///@}


    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;


    /// @name Recovery of Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! fluid) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{F,n+1}\f$ onto the fluid at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{F,n+1}\f$ onto the fluid at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Epetra_Vector> fgprev_;

    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> disgprev_;

    //! inner velocity solution of fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> soliprev_;

    //! interface velocity solution of the fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solgprev_;

    //! inner ALE displacement solution at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solialeprev_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddialeinc_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgicur_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggprev_;

    //! block \f$F_{\Gamma I,i+1}^G\f$ of fluid shape derivatives matrix at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgicur_;

    //! block \f$F_{\Gamma I,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration
    //! \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}^G\f$ of fluid shape derivatives matrix at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggcur_;

    //! block \f$F_{\Gamma\Gamma,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration
    //! \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggprev_;

    //@}

  };  // class monolithic_plain
}  // namespace FPSI

#endif
