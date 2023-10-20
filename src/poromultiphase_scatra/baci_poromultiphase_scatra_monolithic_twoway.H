/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled monolithic algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_SCATRA_MONOLITHIC_TWOWAY_H
#define BACI_POROMULTIPHASE_SCATRA_MONOLITHIC_TWOWAY_H

#include "baci_inpar_solver.H"
#include "baci_poromultiphase_scatra_monolithic.H"
#include "baci_poromultiphase_scatra_utils.H"

#include <Teuchos_Time.hpp>

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
  class SparseOperator;
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class Solver;
  class Equilibration;
  enum class EquilibrationMethod;
}  // namespace CORE::LINALG

namespace POROMULTIPHASESCATRA
{
  //! monolithic coupling algorithm of poromultiphasescatra framework
  class PoroMultiPhaseScaTraMonolithicTwoWay : public PoroMultiPhaseScaTraMonolithic
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraMonolithicTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    /// setup
    void SetupSystem() override;

    /// setup solver (only needed for poromultiphase monolithic coupling)
    void SetupSolver() override;

    /// time step of coupled problem
    void TimeStep() override;

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> Extractor() const
    {
      return blockrowdofmap_;
    }

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const { return combinedDBCMap_; };


   protected:
    //! Setup Newton-Raphson iteration
    void SetupNewton();

    //! Setup full map
    virtual void SetupMaps();

    //! Setup monolithic rhs-vector
    virtual void SetupRHS();

    //! build the combined dirichletbcmap
    virtual void BuildCombinedDBCMap();

    //! build the block null spaces
    virtual void BuildBlockNullSpaces();

    //! create the linear solver
    void CreateLinearSolver(
        const Teuchos::ParameterList& solverparams, const INPAR::SOLVER::SolverType solvertype);

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> iterinc);

    //! extract the field vectors from a given composed vector x.
    /*!
     \param x   (i) composed vector that contains all field vectors
     \param stx (o) structural vector (e.g. displacements)
     \param flx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations)
     and pressures of artery network
     \param scx (o) scatra vector (primary variables of scatra field, i.e. mass fraction)
      and mass fractions in 1D artery network
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& stx, Teuchos::RCP<const Epetra_Vector>& flx,
        Teuchos::RCP<const Epetra_Vector>& scx);

    //! extract only the 3D field vectors from a given composed vector x.
    /*!
     \param x   (i) composed vector that contains all field vectors
     \param stx (o) structural vector (e.g. displacements)
     \param flx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations)
     of 3D field
     \param scx (o) scatra vector (primary variables of scatra field, i.e. mass fraction)
     of 3D field
     */
    void Extract3DFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& stx, Teuchos::RCP<const Epetra_Vector>& flx,
        Teuchos::RCP<const Epetra_Vector>& scx);

    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(Epetra_Vector& f,  //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector>
            pv,  //!< vector containing structural + fluid dofs, i.e. poro dofs
        Teuchos::RCP<const Epetra_Vector> sv  //!< vector containing only scatra dofs
    );

    //! setup monolithic system matrix
    virtual void SetupSystemMatrix();

    //! print header
    void PrintHeader();

    //! solve linear system of equations
    void LinearSolve();

    //! convergence check
    bool Converged();

    //! build norms
    virtual void BuildConvergenceNorms();

    //! output
    void NewtonOutput();

    //! check for convergence
    void NewtonErrorCheck();

    //! update the single fields after convergence
    void UpdateFieldsAfterConvergence();

    //! update the scatra field
    virtual void UpdateScatra(Teuchos::RCP<const Epetra_Vector> scatrainc);

    //! return structure fluid coupling sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> PoroFluidScatraCouplingMatrix();

    //! return scatra structure coupling sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ScatraStructCouplingMatrix();

    //! return scatra fluid coupling sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ScatraPoroFluidCouplingMatrix();

    //! evaluate scatra field
    virtual void EvaluateScatra();

    //! evaluate porofluid-scatra coupling sparse matrix
    void ApplyPoroFluidScatraCouplMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> k_pfs);

    //! evaluate scatra-structure coupling sparse matrix
    void ApplyScatraStructCouplMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> k_sps);

    //! evaluate scatra-porofluid coupling sparse matrix
    void ApplyScatraPoroFluidCouplMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> k_spf);

    // update the single fields after convergence
    void PrintStructureDisabledInfo();

    //! FD-Check
    void PoroMultiPhaseScaTraFDCheck();

    //! convergence tolerance (increment)
    double ittolinc_;
    //! convergence tolerance (residual)
    double ittolres_;
    //! maximally permitted iterations
    int itmax_;
    //! minimally necessary iterations
    int itmin_;
    //! current iteration step
    int itnum_;

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockrowdofmap_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<CORE::LINALG::Equilibration> equilibration_;

    //! equilibration method applied to system matrix
    CORE::LINALG::EquilibrationMethod equilibration_method_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<Epetra_Map> combinedDBCMap_;

    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

    Teuchos::RCP<Epetra_Vector> iterinc_;  //!< increment between Newton steps k and k+1
    //!< \f$\Delta{x}^{<k>}_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> rhs_;  //!< rhs of struct-fluid-scatra system

    Teuchos::RCP<CORE::LINALG::Solver> solver_;  //!< linear algebraic solver
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ?
    bool solveradapttol_;                        //!< adapt solver tolerance

    // do we solve the structure?
    bool solve_structure_;

    // for building blocks
    int struct_offset_;

    //! block systemmatrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> systemmatrix_;

    //! structure-scatra coupling matrix --> we do not have it (yet)
    // Teuchos::RCP<CORE::LINALG::SparseMatrix> k_pss_;

    //! fluid-scatra coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_pfs_;

    //! scatra-structure coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_sps_;
    //! scatra-fluid coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_spf_;

    double tolinc_;   //!< tolerance residual increment
    double tolfres_;  //!< tolerance force residual

    double tolinc_struct_;   //!< tolerance residual increment for structure displacements
    double tolfres_struct_;  //!< tolerance force residual for structure displacements

    double tolinc_fluid_;   //!< tolerance residual increment for fluid
    double tolfres_fluid_;  //!< tolerance force residual for fluid

    double tolinc_scatra_;   //!< tolerance residual increment for scatra
    double tolfres_scatra_;  //!< tolerance force residual for scatra

    double normrhs_;  //!< norm of residual forces

    double normrhsfluid_;  //!< norm of residual forces (fluid )
    double normincfluid_;  //!< norm of residual unknowns (fluid )

    double normrhsstruct_;  //!< norm of residual forces (structure)
    double normincstruct_;  //!< norm of residual unknowns (structure)

    double normrhsscatra_;  //!< norm of residual forces (scatra)
    double normincscatra_;  //!< norm of residual unknowns (scatra)

    double normrhsart_;       //!< norm of residual (artery)
    double normincart_;       //!< norm of residual unknowns (artery)
    double arterypressnorm_;  //!< norm of artery pressure

    double normrhsartsca_;  //!< norm of residual (artery-scatra)
    double normincartsca_;  //!< norm of residual unknowns (artery-scatra)
    double arteryscanorm_;  //!< norm of artery scatra mass-fractions

    double maxinc_;  //!< maximum increment
    double maxres_;  //!< maximum residual

    enum INPAR::POROMULTIPHASESCATRA::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum INPAR::POROMULTIPHASESCATRA::VectorNorm vectornorminc_;   //!< type of norm for increments

    Teuchos::Time timernewton_;  //!< timer for measurement of solution time of newton iterations
    double dtsolve_;             //!< linear solver time
    double dtele_;               //!< time for element evaluation + build-up of system matrix

    //! flag for finite difference check
    INPAR::POROMULTIPHASESCATRA::FDCheck fdcheck_;


  };  // PoroMultiPhaseScatraMonolithic

  //! monolithic coupling algorithm of poromultiphasescatra framework coupled with 1D artery network
  class PoroMultiPhaseScaTraMonolithicTwoWayArteryCoupling
      : public PoroMultiPhaseScaTraMonolithicTwoWay
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraMonolithicTwoWayArteryCoupling(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

   private:
    // Setup full map
    void SetupMaps() override;

    // update the scatra field
    void UpdateScatra(Teuchos::RCP<const Epetra_Vector> scatrainc) override;

    //! extract the field vectors from a given composed vector x.
    /*!
     \param x   (i) composed vector that contains all field vectors
     \param stx (o) structural vector (e.g. displacements)
     \param flx (o) fluid vector (primary variables of fluid field, i.e. pressures or saturations)
     and pressures of artery network
     \param scx (o) scatra vector (primary variables of scatra field, i.e. mass fraction)
      and mass fractions in 1D artery network
     */
    void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& stx, Teuchos::RCP<const Epetra_Vector>& flx,
        Teuchos::RCP<const Epetra_Vector>& scx) override;

    //! setup monolithic system matrix
    void SetupSystemMatrix() override;

    // Setup monolithic rhs-vector
    void SetupRHS() override;

    //! evaluate scatra field
    void EvaluateScatra() override;

    //! build the combined dirichletbcmap
    void BuildCombinedDBCMap() override;

    /// setup
    void SetupSystem() override;

    //! return arteryscatra-artery coupling sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ArteryScatraArteryCouplingMatrix();

    //! evaluate arteryscatra-artery coupling sparse matrix
    void ApplyArteryScatraArteryCouplMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> k_asa);

    //! build the block null spaces
    void BuildBlockNullSpaces() override;

    //! build norms for convergence check
    void BuildConvergenceNorms() override;

    //! dof row map (not splitted), only artery and porofluid
    Teuchos::RCP<Epetra_Map> fullmap_artporo_;

    //! dof row map splitted in (field) blocks, only artery and porofluid
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockrowdofmap_artporo_;

    //! dof row map (not splitted), only artery and artery-scatra
    Teuchos::RCP<Epetra_Map> fullmap_artscatra_;

    //! dof row map splitted in (field) blocks, only artery and artery-scatra
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockrowdofmap_artscatra_;

    //! artscatra-artery coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_asa_;

    //! flag if nodal coupling active or not
    bool nodal_coupl_inactive_;

  };  // PoroMultiPhaseScaTraMonolithicTwoWayArteryCoupling


}  // namespace POROMULTIPHASESCATRA



#endif  // POROMULTIPHASE_SCATRA_MONOLITHIC_TWOWAY_H
