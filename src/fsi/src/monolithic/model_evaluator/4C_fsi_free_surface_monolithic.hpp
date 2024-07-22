/*----------------------------------------------------------------------*/
/*! \file

 \brief Free-surface flow as monolithic problem


 \level 3
 */

#ifndef FOUR_C_FSI_FREE_SURFACE_MONOLITHIC_HPP
#define FOUR_C_FSI_FREE_SURFACE_MONOLITHIC_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_ale.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fsi_overlapprec.hpp"

#include <NOX.H>
#include <NOX_Epetra.H>
#include <Teuchos_TimeMonitor.hpp>

// debug flag to merge the MFSI block matrix to one sparse matrix
// and use the fluid solver to solve for it
// #define BLOCKMATRIXMERGE

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseMatrix;
  class MapExtractor;
  class MatrixColTransform;
}  // namespace Core::LinAlg

namespace Adapter
{
  class Coupling;
  class AleFluidWrapper;
}  // namespace Adapter

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
  }  // namespace UTILS

  /// Base class for Freesurface block preconditioning matrices
  class BlockPreconditioningMatrixFS
      : public Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>
  {
   public:
    BlockPreconditioningMatrixFS(const Core::LinAlg::MultiMapExtractor& maps, Adapter::Fluid& fluid,
        Adapter::AleFluidWrapper& ale, int symmetric, double omega = 1.0, int iterations = 1,
        double fomega = 1.0, int fiterations = 0, FILE* err = nullptr);

    /** \name Mathematical functions */
    //@{

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    //@}

    /// setup of block preconditioners
    void setup_preconditioner() override;

   protected:
    /// (symmetric) Gauss-Seidel block preconditioner
    virtual void sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    /// merge block matrix for direct solve
    void merge_solve(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /// Richardson iteration on one block using the given flags
    static void local_block_richardson(Teuchos::RCP<Core::LinAlg::Preconditioner> solver,
        const Core::LinAlg::SparseMatrix& innerOp, Teuchos::RCP<Epetra_Vector> x,
        Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx, int iterations,
        double omega, FILE* err, const Epetra_Comm& comm);

    /** \name Field solver objects */
    //@{

    Teuchos::RCP<Core::LinAlg::Preconditioner> fluidsolver_;
    Teuchos::RCP<Core::LinAlg::Preconditioner> alesolver_;

    Teuchos::RCP<Core::LinAlg::Preconditioner> constalesolver_;

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
    Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse_;
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
    OverlappingBlockMatrixFS(const Core::LinAlg::MultiMapExtractor& maps, Adapter::Fluid& fluid,
        Adapter::AleFluidWrapper& ale, bool structuresplit, int symmetric, double omega = 1.0,
        int iterations = 1, double fomega = 1.0, int fiterations = 0, FILE* err = nullptr);

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    const char* Label() const override;

    //@}

    /// setup of block preconditioners
    void setup_preconditioner() override;

   protected:
    /// symmetric Gauss-Seidel block preconditioner
    void sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// split is in structural matrix, interface equations belong to fluid block
    bool structuresplit_;

    Adapter::Fluid& fluid_;
    Adapter::AleFluidWrapper& ale_;
  };



  /// monolithic FSI algorithm base for Freesurface Problem
  /*!

  content

    \author --
    \date --
   */
  class MonolithicBaseFS : public Adapter::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit MonolithicBaseFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);


    /// read restart data
    void read_restart(int step) override;

    /// access to Fluid field
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() { return fluid_; }

    /// access to ale field
    const Teuchos::RCP<Adapter::AleFluidWrapper>& ale_field() { return ale_; }

   protected:
    //! @name Time loop building blocks

    /// start a new time step
    void prepare_time_step() override;

    /// take current results for converged and save for next time step
    void update() override;

    /// write output
    void output() override;

    //@}

    //! @name Transfer helpers

    virtual Teuchos::RCP<Epetra_Vector> ale_to_fluid(Teuchos::RCP<Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> ale_to_fluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    Core::Adapter::Coupling& fluid_ale_coupling();

    const Core::Adapter::Coupling& fluid_ale_coupling() const;

   private:
    /// coupling of fluid and ale
    Teuchos::RCP<Core::Adapter::Coupling> coupfa_;

    /// underlying fluid of the FS problem
    Teuchos::RCP<Adapter::Fluid> fluid_;

    /// underlying ale of the FS problem
    Teuchos::RCP<Adapter::AleFluidWrapper> ale_;
  };


  /// base class of all monolithic FSI algorithms for Freesurface Problem
  /*!

  content

    \author --
    \date --
   */
  class MonolithicMainFS : public MonolithicBaseFS,
                           public ::NOX::Epetra::Interface::Required,
                           public ::NOX::Epetra::Interface::Jacobian,
                           public ::NOX::Epetra::Interface::Preconditioner,
                           public ::NOX::Direction::UserDefinedFactory
  {
   public:
    explicit MonolithicMainFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// outer level FSI time loop
    void timeloop(const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface);

    //! @name NOX methods

    /// compute FSI residual
    bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override;

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
        Teuchos::ParameterList* precParams = nullptr) override;

    //@}

    /// create my own direction object
    /*!
      MonolithicMainFS is a (inherits from)
      ::NOX::Direction::UserDefinedFactory. This is an implementation
      detail. This way we can construct a specialized direction object at a
      place where we know about the status tests. This is the whole point
      here. Our specialized direction is of the type NOX::FSI::Newton, the
      normal Newton direction enhanced with adaptive tolerance control for the
      internal linear (iterative) solver.
     */
    Teuchos::RCP<::NOX::Direction::Generic> buildDirection(
        const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) const override;

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void setup_rhs(Epetra_Vector& f, bool firstcall = false) = 0;

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix() = 0;

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void evaluate(
        Teuchos::RCP<const Epetra_Vector> step_increment  ///< increment between time step n and n+1
    );

    /// Extract initial guess from fields
    virtual void initial_guess(Teuchos::RCP<Epetra_Vector> ig) = 0;

    /// apply infnorm scaling to linear block system
    virtual void scale_system(Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void unscale_solution(Epetra_Vector& x, Epetra_Vector& b) {}

   protected:
    /// setup solver for global block system
    virtual Teuchos::RCP<::NOX::Epetra::LinearSystem> create_linear_system(
        Teuchos::ParameterList& nlParams, ::NOX::Epetra::Vector& noxSoln,
        Teuchos::RCP<::NOX::Utils> utils) = 0;

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<::NOX::StatusTest::Combo> create_status_test(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp) = 0;

    /// extract the two field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax) = 0;

    //! @name Access methods for subclasses

    /// output utility
    Teuchos::RCP<::NOX::Utils> utils() const { return utils_; }

    /// full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() const { return blockrowdofmap_.full_map(); }

    /// set full monolithic dof row map
    /*!
      A subclass calls this method (from its constructor) and thereby
      defines the number of blocks, their maps and the block order. The block
      maps must be row maps by themselves and must not contain identical GIDs.
     */
    void set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps);

    /// extractor to communicate between full monolithic map and block maps
    const Core::LinAlg::MultiMapExtractor& extractor() const { return blockrowdofmap_; }

    //@}

    /// flags passed to NOX
    Teuchos::ParameterList& nox_parameter_list() { return noxparameterlist_; }

    /// setup list with default parameters
    void set_default_parameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

    /// add a status test to be used for adaptive linear solver convergence
    void add_status_test(Teuchos::RCP<NOX::FSI::AdaptiveNewtonNormF> test)
    {
      statustests_.push_back(test);
    }

   private:
    /// dof row map splitted in (field) blocks
    Core::LinAlg::MultiMapExtractor blockrowdofmap_;

    /// output utilities
    Teuchos::RCP<::NOX::Utils> utils_;

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
    bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) override;

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
        Teuchos::ParameterList* precParams = nullptr) override;

    //@}

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    void setup_system_matrix() override { setup_system_matrix(*system_matrix()); }

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) = 0;

    //@}

    /// the composed system matrix
    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> system_matrix() const = 0;

    /// apply infnorm scaling to linear block system
    void scale_system(Epetra_Vector& b) override { scale_system(*system_matrix(), b); }

    /// undo infnorm scaling from scaled solution
    void unscale_solution(Epetra_Vector& x, Epetra_Vector& b) override
    {
      unscale_solution(*system_matrix(), x, b);
    }

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void scale_system(Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b) {}

    /// undo infnorm scaling from scaled solution
    virtual void unscale_solution(
        Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
    {
    }

    //@}

   protected:
    /// start a new time step
    void prepare_time_step() override;

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
    void setup_rhs(Epetra_Vector& f, bool firstcall = false) override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

    //@}

    /// Extract initial guess from fields
    void initial_guess(Teuchos::RCP<Epetra_Vector> ig) override;

    /// the composed system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> system_matrix() const override
    {
      return systemmatrix_;
    }

    /// apply infnorm scaling to linear block system
    void scale_system(Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b) override;

    /// undo infnorm scaling from scaled solution
    void unscale_solution(
        Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b) override;

   protected:
    /// setup solver for global block system
    Teuchos::RCP<::NOX::Epetra::LinearSystem> create_linear_system(Teuchos::ParameterList& nlParams,
        ::NOX::Epetra::Vector& noxSoln, Teuchos::RCP<::NOX::Utils> utils) override;

    /// setup of NOX convergence tests
    Teuchos::RCP<::NOX::StatusTest::Combo> create_status_test(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp) override;

    /// extract the two field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    void extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax) override;


   private:
    /// build block vector from field vectors
    void setup_vector(Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> fv,
        Teuchos::RCP<const Epetra_Vector> av);

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrixFS> systemmatrix_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<Core::Adapter::Coupling> icoupfa_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<Core::LinAlg::MatrixColTransform> aigtransform_;

    Teuchos::RCP<Core::LinAlg::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<Core::LinAlg::MatrixColTransform> fmgitransform_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    Inpar::FSI::LinearBlockSolver linearsolverstrategy_;
  };
}  // namespace FSI


namespace NOX
{
  namespace FSI
  {
    /// Special NOX group that always sets Jacobian and RHS at the same time.
    class GroupFS : public ::NOX::Epetra::Group
    {
     public:
      GroupFS(FourC::FSI::MonolithicMainFS& mfsi, Teuchos::ParameterList& printParams,
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
          const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys);

      /// fetch the known Jacobian and RHS from the field solvers
      void capture_system_state();

      ::NOX::Abstract::Group::ReturnType computeF() override;

      ::NOX::Abstract::Group::ReturnType computeJacobian() override;

      ::NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p) override;

     private:
      FourC::FSI::MonolithicMainFS& mfsi_;
    };
  }  // namespace FSI
}  // namespace NOX



FOUR_C_NAMESPACE_CLOSE

#endif
