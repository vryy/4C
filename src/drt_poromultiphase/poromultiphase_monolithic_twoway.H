/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled monolithic solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_MONOLITHIC_TWOWAY_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_MONOLITHIC_TWOWAY_H_

#include "poromultiphase_monolithic.H"
#include <Epetra_Time.h>

#include "../drt_inpar/inpar_poromultiphase.H"

// forward declarations
namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class Solver;
}  // namespace LINALG

namespace DRT
{
  namespace UTILS
  {
    class LocsysManager;
  }
}  // namespace DRT

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseMonolithicTwoWay : public PoroMultiPhaseMonolithic
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseMonolithicTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~PoroMultiPhaseMonolithicTwoWay(){};

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs);

    /// setup
    virtual void SetupSystem();

    /// time step of coupled problem
    virtual void TimeStep();

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const LINALG::MultiMapExtractor> Extractor() const { return blockrowdofmap_; }

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    void Evaluate(Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx,
        const bool firstcall);

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables) public for access from monolithic scatra problem
    void UpdateFieldsAfterConvergence(
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx);

    // access to monolithic rhs vector
    Teuchos::RCP<const Epetra_Vector> RHS() const { return rhs_; }

    // access to monolithic block system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const { return systemmatrix_; }

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const { return combinedDBCMap_; };

   protected:
    //! Newton output to screen
    virtual void NewtonOutput();

    //! Newton error check after loop
    virtual void NewtonErrorCheck();

    //! build the combined dirichletbcmap
    virtual void BuildCombinedDBCMap();

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    virtual void SetupRHS();

    virtual Teuchos::RCP<Epetra_Vector> SetupStructurePartofRHS();

    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(Epetra_Vector& f,         //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> fv   //!< vector containing only fluid dofs
    );

    //! extract the field vectors from a given composed vector.
    /*!
     x is the sum of all increments up to this point.
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations)
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() { SetupSystemMatrix(*systemmatrix_); }

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    /// setup composed system matrix from field solvers
    virtual void SetupMaps();

    // Setup solver for monolithic system
    virtual bool SetupSolver();

    //! build the block null spaces
    virtual void BuildBlockNullSpaces(Teuchos::RCP<LINALG::Solver>& solver);

    //! Evaluate mechanical-fluid system matrix
    virtual void ApplyStrCouplMatrix(
        Teuchos::RCP<LINALG::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
    );

    //! Evaluate fluid-mechanical system matrix
    virtual void ApplyFluidCouplMatrix(
        Teuchos::RCP<LINALG::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
    );

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> iterinc);

    //! return structure fluid coupling sparse matrix
    Teuchos::RCP<LINALG::SparseMatrix> StructFluidCouplingMatrix();

    //! return fluid structure coupling sparse matrix
    Teuchos::RCP<LINALG::SparseMatrix> FluidStructCouplingMatrix();

    //! Solve the linear system of equations
    void LinearSolve();

    //! equilibrate the system of equations
    void EquilibrateSystem();

    //! equilibrate matrix rows
    void EquilibrateMatrixRows(LINALG::SparseMatrix& matrix,  //!< matrix
        const Teuchos::RCP<Epetra_Vector>&
            invrowsums  //!< sums of absolute values of row entries in matrix
        ) const;

    //! compute inverse sums of absolute values of matrix row entries
    void ComputeInvRowSums(const LINALG::SparseMatrix& matrix,  //!< matrix
        const Teuchos::RCP<Epetra_Vector>&
            invrowsums  //!< inverse sums of absolute values of row entries in matrix
        ) const;

    //! Create the linear solver
    virtual void CreateLinearSolver(
        const Teuchos::ParameterList& solverparams, const int solvertype);

    //! Setup Newton-Raphson
    void SetupNewton();

    //! Print Header to screen
    void PrintHeader();

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables)
    void UpdateFieldsAfterConvergence();

    // build norms for convergence check
    void BuildConvergenceNorms();

    void PoroFDCheck();

    // check for convergence
    bool Converged();

    /// Print user output that structure field is disabled
    void PrintStructureDisabledInfo();


    //! linear solver
    bool directsolve_;
    //! convergence tolerance for increments
    double ittolinc_;
    //! convergence tolerance for residuals
    double ittolres_;
    //! maximally permitted iterations
    int itmax_;
    //! minimally necessary iterations
    int itmin_;
    //! current iteration step
    int itnum_;
    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

    Teuchos::RCP<Epetra_Vector> iterinc_;  //!< increment between Newton steps k and k+1
    //!< \f$\Delta{x}^{<k>}_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> rhs_;  //!< rhs of Poroelasticity system

    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver
    double solveradaptolbetter_;           //!< tolerance to which is adpated ?
    bool solveradapttol_;                  //!< adapt solver tolerance


    //@}

    //! @name Global matrixes

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

    //! structure-fluid coupling matrix
    Teuchos::RCP<LINALG::SparseOperator> k_sf_;
    //! fluid-structure coupling matrix
    Teuchos::RCP<LINALG::SparseOperator> k_fs_;

    //@}

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> blockrowdofmap_;

    //! flag for equilibration of global system of equations
    INPAR::POROMULTIPHASE::EquilibrationMethods equilibration_;

    //! inverse sums of absolute values of row entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invrowsums_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<Epetra_Map> combinedDBCMap_;

    double tolinc_;   //!< tolerance residual increment
    double tolfres_;  //!< tolerance force residual

    double tolinc_struct_;   //!< tolerance residual increment for structure displacements
    double tolfres_struct_;  //!< tolerance force residual for structure displacements

    double tolinc_fluid_;   //!< tolerance residual increment for fluid
    double tolfres_fluid_;  //!< tolerance force residual for fluid

    double normrhs_;  //!< norm of residual forces
    double norminc_;  //!< norm of residual unknowns

    double normrhsfluid_;  //!< norm of residual forces (fluid )
    double normincfluid_;  //!< norm of residual unknowns (fluid )

    double normrhsstruct_;  //!< norm of residual forces (structure)
    double normincstruct_;  //!< norm of residual unknowns (structure)

    enum INPAR::POROMULTIPHASE::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum INPAR::POROMULTIPHASE::VectorNorm vectornorminc_;   //!< type of norm for increments

    Epetra_Time timernewton_;  //!< timer for measurement of solution time of newton iterations
    double dtsolve_;           //!< linear solver time
    double dtele_;             //!< time for element evaluation + build-up of system matrix

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

    //! flag for finite difference check
    INPAR::POROMULTIPHASE::FDCheck fdcheck_;

  };  // PoroMultiPhasePartitioned

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseMonolithicTwoWayArteryCoupling : public PoroMultiPhaseMonolithicTwoWay
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseMonolithicTwoWayArteryCoupling(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~PoroMultiPhaseMonolithicTwoWayArteryCoupling(){};

    //! extract the field vectors from a given composed vector.
    /*!
     x is the sum of all increments up to this point.
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations)
                   and artery pressure
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx);

   protected:
    /// setup composed system matrix from field solvers
    virtual void SetupMaps();

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    /// setup global rhs
    virtual void SetupRHS();

    //! build the combined dirichletbcmap
    virtual void BuildCombinedDBCMap();

    //! Create the linear solver
    virtual void CreateLinearSolver(
        const Teuchos::ParameterList& solverparams, const int solvertype);

    //! build the block null spaces
    virtual void BuildArteryBlockNullSpace(
        Teuchos::RCP<LINALG::Solver>& solver, const int& arteryblocknum);

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_artporo_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> blockrowdofmap_artporo_;
  };


}  // namespace POROMULTIPHASE

#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_MONOLITHIC_TWOWAY_H_ */
