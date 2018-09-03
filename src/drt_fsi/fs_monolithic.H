

#ifndef FS_MONOLITHIC_H
#define FS_MONOLITHIC_H

#include "fsi_overlapprec.H"

#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_ale/ale.H"

#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
//#define BLOCKMATRIXMERGE

namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}  // namespace LINALG

namespace ADAPTER
{
  class Coupling;
  class AleFluidWrapper;
}  // namespace ADAPTER

namespace NOX
{
  namespace FSI
  {
    class AdaptiveNewtonNormF;
  }
}  // namespace NOX

namespace FSI
{
  // forward declarations

  namespace UTILS
  {
    class DebugWriter;
    class MatrixColTransform;
  }  // namespace UTILS

  /// Base class for Freesurface block preconditioning matrices
  class BlockPreconditioningMatrixFS
      : public LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
  {
   public:
    BlockPreconditioningMatrixFS(const LINALG::MultiMapExtractor& maps, ::ADAPTER::Fluid& fluid,
        ::ADAPTER::AleFluidWrapper& ale, int symmetric, double omega = 1.0, int iterations = 1,
        double fomega = 1.0, int fiterations = 0, FILE* err = NULL);

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();

   protected:
    /// (symmetric) Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    /// merge block matrix for direct solve
    void MergeSolve(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// Richardson iteration on one block using the given flags
    static void LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
        const LINALG::SparseMatrix& innerOp, Teuchos::RCP<Epetra_Vector> x,
        Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx, int iterations,
        double omega, FILE* err, const Epetra_Comm& comm);

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<LINALG::Preconditioner> fluidsolver_;
    Teuchos::RCP<LINALG::Preconditioner> alesolver_;

    Teuchos::RCP<LINALG::Preconditioner> constalesolver_;

    //@}

    /// Symmetric block GS preconditioner in monolithic FSI or ordinary GS
    int symmetric_;

    /// \name Richardson iteration
    //@{

    double omega_;
    int iterations_;
    double fomega_;
    int fiterations_;

    //@}

    /// log file
    FILE* err_;

#ifdef BLOCKMATRIXMERGE
    /// debug merged sparse
    Teuchos::RCP<LINALG::SparseMatrix> sparse_;
#endif
  };


  /// special version of block matrix that includes the FSI block preconditioner
  /*!
    The normal block matrix is enhanced by a ApplyInverse() method that does
    the Gauss-Seidel block preconditioning explicitly for FSI block matrices.
  */
  class OverlappingBlockMatrixFS : public BlockPreconditioningMatrixFS
  {
   public:
    /// construction
    OverlappingBlockMatrixFS(const LINALG::MultiMapExtractor& maps, ::ADAPTER::Fluid& fluid,
        ::ADAPTER::AleFluidWrapper& ale, bool structuresplit, int symmetric, double omega = 1.0,
        int iterations = 1, double fomega = 1.0, int fiterations = 0, FILE* err = NULL);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

    /// setup of block preconditioners
    virtual void SetupPreconditioner();

   protected:
    /// symmetric Gauss-Seidel block preconditioner
    virtual void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;

    ::ADAPTER::Fluid& fluid_;
    ::ADAPTER::AleFluidWrapper& ale_;
  };



  /// monolithic FSI algorithm base for Freesurface Problem
  /*!

  content

    \author --
    \date --
   */
  class MonolithicBaseFS : public ::ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit MonolithicBaseFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// virtual destructor to support polymorph destruction
    virtual ~MonolithicBaseFS();

    /// read restart data
    void ReadRestart(int step);

    /// access to Fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }

    /// access to ale field
    const Teuchos::RCP<ADAPTER::AleFluidWrapper>& AleField() { return ale_; }

   protected:
    //! @name Time loop building blocks

    /// start a new time step
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// write output
    virtual void Output();

    //@}

    //! @name Transfer helpers

    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    ADAPTER::Coupling& FluidAleCoupling();

    const ADAPTER::Coupling& FluidAleCoupling() const;

   private:
    /// coupling of fluid and ale
    Teuchos::RCP<ADAPTER::Coupling> coupfa_;

    /// underlying fluid of the FS problem
    Teuchos::RCP<::ADAPTER::Fluid> fluid_;

    /// underlying ale of the FS problem
    Teuchos::RCP<::ADAPTER::AleFluidWrapper> ale_;
  };


  /// base class of all monolithic FSI algorithms for Freesurface Problem
  /*!

  content

    \author --
    \date --
   */
  class MonolithicMainFS : public MonolithicBaseFS,
                           public NOX::Epetra::Interface::Required,
                           public NOX::Epetra::Interface::Jacobian,
                           public NOX::Epetra::Interface::Preconditioner,
                           public NOX::Direction::UserDefinedFactory
  {
   public:
    explicit MonolithicMainFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// outer level FSI time loop
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //! @name NOX methods

    /// compute FSI residual
    bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

    /// preconditioner
    bool computePreconditioner(
        const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams = 0);

    //@}

    /// create my own direction object
    /*!
      MonolithicMainFS is a (inherits from)
      NOX::Direction::UserDefinedFactory. This is an implementation
      detail. This way we can construct a specialized direction object at a
      place where we know about the status tests. This is the whole point
      here. Our specialized direction is of the type NOX::FSI::Newton, the
      normal Newton direction enhanced with adaptive tolerance control for the
      internal linear (iterative) solver.
     */
    virtual Teuchos::RCP<NOX::Direction::Generic> buildDirection(
        const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const;

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall = false) = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() = 0;

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc  ///< increment between time step n and n+1
    );

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) = 0;

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) {}

   protected:
    /// setup solver for global block system
    virtual Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem(
        Teuchos::ParameterList& nlParams, NOX::Epetra::Vector& noxSoln,
        Teuchos::RCP<NOX::Utils> utils) = 0;

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp) = 0;

    /// extract the two field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax) = 0;

    //! @name Access methods for subclasses

    /// output utility
    Teuchos::RCP<NOX::Utils> Utils() const { return utils_; }

    /// full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return blockrowdofmap_.FullMap(); }

    /// set full monolithic dof row map
    /*!
      A subclass calls this method (from its constructor) and thereby
      defines the number of blocks, their maps and the block order. The block
      maps must be row maps by themselves and must not contain identical GIDs.
     */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// extractor to communicate between full monolithic map and block maps
    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //@}

    /// flags passed to NOX
    Teuchos::ParameterList& NOXParameterList() { return noxparameterlist_; }

    /// setup list with default parameters
    void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

    /// add a status test to be used for adaptive linear solver convergence
    void AddStatusTest(Teuchos::RCP<NOX::FSI::AdaptiveNewtonNormF> test)
    {
      statustests_.push_back(test);
    }

   private:
    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    /// output utilities
    Teuchos::RCP<NOX::Utils> utils_;

    /// flags passed to NOX
    Teuchos::ParameterList noxparameterlist_;

    /// keep the status tests available so we can connect them with our
    /// adaptive Newton direction
    std::vector<Teuchos::RCP<NOX::FSI::AdaptiveNewtonNormF>> statustests_;

    /// @name special debugging output

    Teuchos::RCP<UTILS::DebugWriter> fdbg_;

    //@}
  };



  /// Monolithic FSI with block system matrix for Freesurface Problem
  class BlockMonolithicFS : public MonolithicMainFS
  {
   public:
    explicit BlockMonolithicFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! @name NOX methods

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

    /// preconditioner
    bool computePreconditioner(
        const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams = 0);

    //@}

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() { SetupSystemMatrix(*SystemMatrix()); }

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) = 0;

    //@}

    /// the composed system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const = 0;

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) { ScaleSystem(*SystemMatrix(), b); }

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b)
    {
      UnscaleSolution(*SystemMatrix(), x, b);
    }

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(
        LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
    {
    }

    //@}

   protected:
    /// start a new time step
    virtual void PrepareTimeStep();

   private:
    /*! \brief Counter of iterations to reuse the block matrix preconditioner
     *
     *  Rebuild preconditioner as soon as this counter is zero.
     *
     *  \note We enforce rebuilding the preconditioner at the beginning of
     *  every time step.
     */
    int precondreusecount_;
  };



  /// monolithic Freesurface algorithm
  /*!
  adapted from MonolithicStructureSplit

    \sa
    \author
    \date
   */
  class MonolithicFS : public BlockMonolithicFS
  {
   public:
    explicit MonolithicFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(
        LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

   protected:
    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem(Teuchos::ParameterList& nlParams,
        NOX::Epetra::Vector& noxSoln, Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp);

    /// extract the two field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax);


   private:
    /// build block vector from field vectors
    void SetupVector(Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> fv,
        Teuchos::RCP<const Epetra_Vector> av);

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrixFS> systemmatrix_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<ADAPTER::Coupling> icoupfa_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;

    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmgitransform_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;
  };
}  // namespace FSI


namespace NOX
{
  namespace FSI
  {
    /// Special NOX group that always sets Jacobian and RHS at the same time.
    class GroupFS : public NOX::Epetra::Group
    {
     public:
      GroupFS(::FSI::MonolithicMainFS& mfsi, Teuchos::ParameterList& printParams,
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x,
          const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys);

      /// fetch the known Jacobian and RHS from the field solvers
      void CaptureSystemState();

      virtual NOX::Abstract::Group::ReturnType computeF();

      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p);

     private:
      ::FSI::MonolithicMainFS& mfsi_;
    };
  }  // namespace FSI
}  // namespace NOX



#endif
