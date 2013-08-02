/*!----------------------------------------------------------------------
\file fsi_monolithic_xfem.cpp
\brief Control routine for monolithic XFSI using XFEM

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHIC_XFEM_H
#define FSI_MONOLITHIC_XFEM_H

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_fsi_xfem/fsi_xfem_algorithm.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;

    class DebugWriter;
  }
}

namespace FLD
{
  namespace UTILS
  {
    class FluidXFluidMapExtractor;
  }
}

namespace FSI
{

  /// monolithic XFSI algorithm
  /*!

//TODO:
    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

    \author Shadan Shahmiri
    \date  11/2011
  */
  class MonolithicXFEM : public AlgorithmXFEM
  {
    friend class FSIResultTest;

  public:
    explicit MonolithicXFEM(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix

    </pre>

    */
    virtual void SetupSystem();

    void SetupNewSystem();

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall);

    /// setup composed system matrix from field solvers
    //void SetupSystemMatrixold(LINALG::BlockSparseMatrixBase& mat);

    virtual void SetupSystemMatrix();

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);


    //! @name Access methods for subclasses

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
    */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<LINALG::MultiMapExtractor> Extractor() const { return blockrowdofmap_; }

    //@}

    // newton for xfem
    virtual void Newton();

  protected:



    void LinearSolve();

    //! create linear solver (setup of parameter lists, etc...)
    void CreateLinearSolver();


    virtual void Update();




    //! time loop of the monolithic system
    virtual void Timeloop();

    //! prepare the time step for fluid and structure
    virtual void PrepareTimeStep();


    //! extract the two field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx);


    virtual void ReadRestart(int step);

    virtual void Output();



    virtual void RecoverLagrangeMultiplier();

    virtual void BuildCovergenceNorms();

//    void SetDefaultParameters(const Teuchos::ParameterList& fsidyn,
//                              Teuchos::ParameterList& list);


    //! @name Print Newton Iteration info

    //! print to screen information about residual forces and displacements
    void PrintNewtonIter();

    //! contains text to PrintNewtonIter
    void PrintNewtonIterText();

    //! contains header to PrintNewtonIter
    void PrintNewtonIterHeader();

    //! print statistics of converged Newton-Raphson iteration
    //void PrintNewtonConv();

    //@}


    //! structural displacement increment of interface DOFs due to predictor or inhomogeneous DBCs
    Teuchos::RCP<Epetra_Vector> ddgpred_;

    //! @name General solver parameters
    //@{

    bool solveradapttol_;  //!< adapt solver tolerance
    double solveradaptolbetter_;  //!< tolerance to which is adpated ????
    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver

    bool merge_fsi_blockmatrix_; //!< merged blockmatrix, used for solver setup!

    //@}


  private:

    const Teuchos::ParameterList& fsidyn_;  //!< FSI dynamic parameter list


    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     double fluidscale);

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;




    //! outer iteration loop when active fluid dofsets change, calls inner Newton-Raphson iterations within each outer iteration
    void Solve();

    //! inner iteration loop return if converged or not (in case of changing fluid maps)
    bool NewtonFull();

    //! Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc ///< increment between time step n and n+1
        );

    //! check convergence of Newton iteration (public)
    bool Converged();

    //! @name Output

//    //! print to screen information about residual forces and displacements
//    //! \author lw (originally in STR) \date 12/07
//    void PrintNewtonIter();
//
//    //! contains text to PrintNewtonIter
//    //! \author lw (originally in STR) \date 12/07
//    void PrintNewtonIterText(
//      FILE* ofile  //!< output file handle
//      );
//
//    //! contains header to PrintNewtonIter
//    //! \author lw (originally) \date 12/07
//    void PrintNewtonIterHeader(
//      FILE* ofile  //!< output file handle
//      );
//
//    //! print statistics of converged Newton-Raphson iteration
//    void PrintNewtonConv();



    /// block system matrix
    //Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;


    //! dofrowmap splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> blockrowdofmap_;



    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

//    Teuchos::RCP<UTILS::MatrixRowColTransform> sggtransform_;
//    Teuchos::RCP<UTILS::MatrixRowTransform> sgitransform_;
//    Teuchos::RCP<UTILS::MatrixColTransform> sigtransform_;
//    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;
//
//    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;
//    Teuchos::RCP<UTILS::MatrixColTransform> fmgitransform_;
//
//    Teuchos::RCP<UTILS::MatrixColTransform> fsaigtransform_;
//    Teuchos::RCP<UTILS::MatrixColTransform> fsmgitransform_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;
//
//    Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor> xfluidfluidsplitter_;
//
//    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;
//
//    int currentstep_;
//
//    int relaxing_ale_every_;
//
//    bool relaxing_ale_;

    /// @name Some quantities to recover the Langrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the structure)
    //! evalutated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    // lambda lives at the slave side (here at stucture)
    //Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the structure at current NOX iteration \f$i+1\f$
    //Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{S,n+1}\f$ onto the structure at previous NOX iteration \f$i\f$
    //Teuchos::RCP<const Epetra_Vector> fgpre_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddiinc_;

    //! inner displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solipre_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solgpre_;

    //! block \f$S_{\Gamma I,i+1}\f$ of structural matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sgicur_;

    //! block \f$S_{\Gamma I,i}\f$ of structural matrix at previous NOX iteration \f$i\f$
    //Teuchos::RCP<const LINALG::SparseMatrix> sgipre_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sggcur_;

    //! block \f$S_{\Gamma\Gamma,i}\f$ of structural matrix at previous NOX iteration \f$i\f$
    //Teuchos::RCP<const LINALG::SparseMatrix> sggpre_;


    bool firstcall_;

    // sum of increments, increment w.r.t the old timestep t^n
    Teuchos::RCP<Epetra_Vector> x_sum_;

    /// global Newton increment = iteration increment Delta x = x^n+1_i+1 - x^n+1_i
    Teuchos::RCP<Epetra_Vector> iterinc_;

    /// global residual vector
    Teuchos::RCP<Epetra_Vector> rhs_;

    /// global vector for combined fluid and structure system filled with zeros used for DBCs
    Teuchos::RCP<Epetra_Vector> zeros_;

    /// output stream
    Teuchos::RCP<std::ofstream> log_;

    /// @name special debugging output

    Teuchos::RCP<UTILS::DebugWriter> sdbg_;
    Teuchos::RCP<UTILS::DebugWriter> fdbg_;


    //--------------------------------------------------------------------------//

    double tolrhs_;

    int ns_;  //!< length of structural dofs
    int nf_;  //!< length of fluid dofs
    int nfv_; //!< length of fluid velocity dofs
    int nfp_; //!< length of fluid pressure dofs
    int nall_;

    double normrhs_;  //!< norm of residual forces
    double norminc_;  //!< norm of solution increment

    //L2-NORMS
    //--------------------------------------------------------------------------//
    double normstrrhsL2_;     //!< norm of structural residual
    double normflvelrhsL2_;   //!< norm of inner fluid velocity residual
    double normflpresrhsL2_;  //!< norm of fluid pressure residual

    //--------------------------------------------------------------------------//
    double normstrincL2_;     //!< norm of inner structural increment
    double normflvelincL2_;   //!< norm of inner fluid velocity residual forces
    double normflpresincL2_;  //!< norm of fluid pressure residual forces
    //--------------------------------------------------------------------------//

    //Inf-NORMS
    //--------------------------------------------------------------------------//
    double normstrrhsInf_;     //!< norm of structural residual
    double normflvelrhsInf_;   //!< norm of inner fluid velocity residual
    double normflpresrhsInf_;  //!< norm of fluid pressure residual

    //--------------------------------------------------------------------------//
    double normstrincInf_;     //!< norm of inner structural increment
    double normflvelincInf_;   //!< norm of inner fluid velocity residual forces
    double normflpresincInf_;  //!< norm of fluid pressure residual forces

    //--------------------------------------------------------------------------//
    int iter_;                //!< iteration step
    int itermax_;     //!< maximally permitted iterations
    int itermin_;     //!< minimally requested iteration

    int iter_outer_;          //!< iteration counter for outer Newton loop
    int itermax_outer_;       //!< maximally permitted iterations

    //@}

    //! @name Convergence criterion and convergence tolerances

    enum INPAR::FSI::ConvNorm normtypeinc_;  //!< convergence check for increment
    enum INPAR::FSI::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum INPAR::FSI::BinaryOp combincfres_;  //!< binary operator to combine temperatures and forces

    double tolinc_;   //!< tolerance residual temperatures
    double tolfres_;  //!< tolerance force residual

    double TOL_DIS_RES_L2_;
    double TOL_DIS_RES_INF_;
    double TOL_DIS_INC_L2_;
    double TOL_DIS_INC_INF_;
    double TOL_PRE_RES_L2_;
    double TOL_PRE_RES_INF_;
    double TOL_PRE_INC_L2_;
    double TOL_PRE_INC_INF_;
    double TOL_VEL_RES_L2_;
    double TOL_VEL_RES_INF_;
    double TOL_VEL_INC_L2_;
    double TOL_VEL_INC_INF_;
    //@}



  };
}

#endif
