/*----------------------------------------------------------------------*/
/*! \file

 \brief General framework for monolithic fpsi solution schemes

\level 3

\maintainer  Johannes Kremheller
 */

/*----------------------------------------------------------------------*/

#ifndef FPSI_MONOLITHIC_H
#define FPSI_MONOLITHIC_H

// FPSI includes
#include "fpsi.H"
#include "fpsi_coupling.H"

// FSI includes
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_fsi/fsi_monolithic.H"

// POROELAST includes
#include "../drt_poroelast/poro_base.H"

// ALE includes
#include "../drt_ale/ale.H"
#include "../drt_adapter/ad_ale_fpsi.H"

// ADAPTER includes
#include "../drt_adapter/ad_fld_fluid_fpsi.H"

// LINALG includes
#include "../linalg/linalg_mapextractor.H"

// INPAR includes
#include "../drt_inpar/inpar_fpsi.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}  // namespace LINALG

namespace POROELAST
{
  class Monolithic;
}

namespace ADAPTER
{
  class Coupling;
  class AleFpsiWrapper;
}  // namespace ADAPTER

namespace FSI
{
  namespace UTILS
  {
    class DebugWriter;
    class MonolithicDebugWriter;
    class MatrixColTransform;
  }  // namespace UTILS
}  // namespace FSI
/*----------------------------------------------------------------------*/

namespace FPSI
{
  class MonolithicBase : public FPSI_Base
  {
   public:
    //! ctor
    explicit MonolithicBase(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    //! dtor
    virtual ~MonolithicBase();

    //! read restart data
    virtual void ReadRestart(int step);

    //! start a new time step
    virtual void PrepareTimeStep();

    //! take current results for converged and save for next time step
    virtual void Update();

    //! calculate stresses, strains, energies
    virtual void PrepareOutput();

    //! Output routine accounting for Lagrange multiplier at the interface
    virtual void Output();

    //! @name access sub-fields
    const Teuchos::RCP<POROELAST::Monolithic>& PoroField() { return poroelast_subproblem_; };
    const Teuchos::RCP<ADAPTER::FluidFPSI>& FluidField() { return fluid_subproblem_; };
    const Teuchos::RCP<ADAPTER::AleFpsiWrapper>& AleField() { return ale_; };

    //@}

    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap;
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap;

   protected:
    //! underlying poroelast problem
    Teuchos::RCP<POROELAST::Monolithic> poroelast_subproblem_;
    //! underlying fluid of the FSI problem
    Teuchos::RCP<ADAPTER::FluidFPSI> fluid_subproblem_;
    //! underlying ale of the FSI problem
    Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale_;

    //! flag defines if FSI Interface exists for this problem
    bool FSI_Interface_exists_;

   public:
    //! FPSI coupling object (does the interface evaluations)
    Teuchos::RCP<FPSI::FPSICoupling>& FPSICoupl() { return fpsicoupl_; }

    //! @name Access General Couplings
    ADAPTER::Coupling& FluidAleCoupling() { return *coupfa_; }

    const ADAPTER::Coupling& FluidAleCoupling() const { return *coupfa_; }

    // Couplings for FSI
    ADAPTER::Coupling& StructureFluidCoupling_FSI() { return *coupsf_fsi_; }
    ADAPTER::Coupling& StructureAleCoupling_FSI() { return *coupsa_fsi_; }
    ADAPTER::Coupling& InterfaceFluidAleCoupling_FSI() { return *icoupfa_fsi_; }

    const ADAPTER::Coupling& StructureFluidCoupling_FSI() const { return *coupsf_fsi_; }
    const ADAPTER::Coupling& StructureAleCoupling_FSI() const { return *coupsa_fsi_; }
    const ADAPTER::Coupling& InterfaceFluidAleCoupling_FSI() const { return *icoupfa_fsi_; }

    //@}

   protected:
    //! @name Transfer helpers
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    virtual Teuchos::RCP<Epetra_Vector> StructToFluid_FSI(
        Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct_FSI(
        Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToAle_FSI(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct_FSI(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle_FSI(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid_FSI(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidInterface_FSI(
        Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

   private:
    //! FPSI - COUPLING
    //! coupling of fluid and ale in the entire fluid volume
    Teuchos::RCP<ADAPTER::Coupling> coupfa_;

    //! FSI - COUPLING
    //! coupling of structure and fluid at the interface
    Teuchos::RCP<ADAPTER::Coupling> coupsf_fsi_;
    //! coupling of structure and ale at the interface
    Teuchos::RCP<ADAPTER::Coupling> coupsa_fsi_;
    //! coupling of fluid and ale in the entire fluid volume
    Teuchos::RCP<ADAPTER::Coupling> coupfa_fsi_;
    //! coupling of all interface fluid and ale dofs
    Teuchos::RCP<ADAPTER::Coupling> icoupfa_fsi_;
    //! coupling of FPSI+FSI interface overlapping dofs of structure and freefluid
    Teuchos::RCP<ADAPTER::Coupling> iffcoupsf_fsi_;

    //! FPSI Coupling Object
    Teuchos::RCP<FPSI::FPSICoupling> fpsicoupl_;
  };
  // MonolithicBase

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //<<<<<<<<<<<<<<<<<<<<<<  MonolithicBase -> Monolithic  >>>>>>>>>>>>>>>>>>>>>
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  class Monolithic : public MonolithicBase
  {
    friend class FSI::UTILS::MonolithicDebugWriter;

   public:
    //! ctor
    explicit Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    //! setup fpsi system
    virtual void SetupSystem();

    //! setup fsi part of the system
    virtual void SetupSystem_FSI();

    //! perform time loop
    void Timeloop();

    //! prepare time loop
    void PrepareTimeloop();

    //! solve one time step
    virtual void TimeStep();

    //! perform result test
    void TestResults(const Epetra_Comm& comm);

    //! build RHS vector from sub fields
    virtual void SetupRHS(bool firstcall = false) = 0;

    //! build system matrix form sub fields + coupling
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) = 0;
    //! build system matrix form sub fields + coupling
    virtual void SetupSystemMatrix() { SetupSystemMatrix(*SystemMatrix()); }

    //! access system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const = 0;

    /// setup solver
    virtual void SetupSolver();

    /// Recover the Lagrange multiplier at the interface   mayr.mt (03/2012)
    virtual void RecoverLagrangeMultiplier()
    {
      dserror("RecoverLagrangeMultiplier: Not Implemented in Base Class!");
    }

    /// Extract specific columns from Sparse Matrix
    void ExtractColumnsfromSparse(Teuchos::RCP<Epetra_CrsMatrix> src,  ///< source Matrix
        const Teuchos::RCP<const Epetra_Map>&
            colmap,  ///< map with column gids to be extracted! (gid which are not in the source
                     ///< Matrix will be ignored!)
        Teuchos::RCP<Epetra_CrsMatrix> dst);  ///< destination Matrix (will be filled!)

    //! Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc);  ///< increment between time step n and n+1

    //! setup of newton scheme
    void SetupNewton();

    //! finite difference check for fpsi systemmatrix
    void FPSIFDCheck();

    //! solve linear system
    void LinearSolve();

    //! solve using line search method
    void LineSearch(Teuchos::RCP<LINALG::SparseMatrix>& sparse);

    //! create linear solver (setup of parameter lists, etc...)
    void CreateLinearSolver();

    //! build convergence norms after solve
    void BuildConvergenceNorms();

    //! print header and results of newton iteration to screen
    void PrintNewtonIter();

    //! print header of newton iteration
    void PrintNewtonIterHeader(FILE* ofile);

    //! print results of newton iteration
    void PrintNewtonIterText(FILE* ofile);

    //! perform convergence check
    bool Converged();

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const { return blockrowdofmap_.FullMap(); }

    //! map of all dofs on Dirichlet-Boundary
    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    //! extractor to communicate between full monolithic map and block maps
    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //! set conductivity (for fps3i)
    void SetConductivity(double conduct);

    //! external acces to rhs vector (used by xfpsi)
    Teuchos::RCP<Epetra_Vector>& RHS() { return rhs_; }  // TodoAge: will be removed again!

   protected:
    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;
    //! dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;
    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;
    //! increment between Newton steps k and k+1
    Teuchos::RCP<Epetra_Vector> iterinc_;
    Teuchos::RCP<Epetra_Vector> iterincold_;
    //! zero vector of full length
    Teuchos::RCP<Epetra_Vector> zeros_;
    //! linear algebraic solver
    Teuchos::RCP<LINALG::Solver> solver_;
    //! rhs of FPSI system
    Teuchos::RCP<Epetra_Vector> rhs_;
    Teuchos::RCP<Epetra_Vector> rhsold_;

    Teuchos::RCP<const Epetra_Vector> meshdispold_;

    Teuchos::RCP<Epetra_Vector> porointerfacedisplacementsold_;

    //! adapt solver tolerancePoroField()->SystemSparseMatrix()
    bool solveradapttol_;
    int linesearch_;
    double linesearch_counter;
    double solveradaptolbetter_;

    bool active_FD_check_;  // indicates if evaluate() is called from FDCheck (firstiter should not
                            // be added anymore!!!)

    /// extract the three field vectors from a given composed vector
    /*
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) poroelast dofs
     \param fx (o) free fluid velocities and pressure
     \param ax (o) ale displacements
     \param firstiter_ (i) firstiteration? - how to evaluate FSI-velocities
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& pfx,
        Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax,
        bool firstiter_) = 0;

    /// setup list with default parameters
    void SetDefaultParameters(const Teuchos::ParameterList& fpsidynparams);

    //! block ids of the monolithic system
    int porofluid_block_;
    int structure_block_;
    int fluid_block_;
    int ale_i_block_;

   private:
    //! flag for direct solver of linear system
    bool directsolve_;

    enum INPAR::FPSI::ConvergenceNorm normtypeinc_;
    enum INPAR::FPSI::ConvergenceNorm normtypefres_;
    enum INPAR::FPSI::BinaryOp combinedconvergence_;

    double toleranceiterinc_;
    double toleranceresidualforces_;
    std::vector<double>
        toleranceresidualforceslist_;  // order of fields: porofluidvelocity, porofluidpressure,
                                       // porostructure, fluidvelocity, fluidpressure, ale
    std::vector<double>
        toleranceiterinclist_;  // order of fields: porofluidvelocity, porofluidpressure,
                                // porostructure, fluidvelocity, fluidpressure, ale

    int maximumiterations_;
    int minimumiterations_;
    double normofrhs_;
    double normofrhsold_;
    double normofiterinc_;
    double normofiterincold_;

    double normrhsfluidvelocity_;
    double normrhsfluidpressure_;
    double normrhsporofluidvelocity_;
    double normrhsporofluidpressure_;
    double normrhsporointerface_;
    double normrhsfluidinterface_;
    double normrhsporostruct_;
    double normrhsfluid_;
    double normrhsale_;

    double normofiterincporostruct_;
    double normofiterincporofluid_;
    double normofiterincfluid_;
    double normofiterincporofluidvelocity_;
    double normofiterincporofluidpressure_;
    double normofiterincfluidvelocity_;
    double normofiterincfluidpressure_;
    double normofiterincale_;
    double normofiterincfluidinterface_;
    double normofiterincporointerface_;

    double sqrtnall_;  //!< squareroot of lenght of all dofs
    double sqrtnfv_;   //!< squareroot of length of fluid velocity dofs
    double sqrtnfp_;   //!< squareroot of length of fluid pressure dofs
    double sqrtnpfv_;  //!< squareroot of length of porofluid velocity dofs
    double sqrtnpfp_;  //!< squareroot of length of porofluid pressure dofs
    double sqrtnps_;   //!< squareroot of length of porostruct dofs
    double sqrtna_;    //!< squareroot of length of ale dofs

    double norm1_alldof_;  //!< sum of absolute values of all dofs
    double norm1_fv_;      //!< sum of absolute fluid velocity values
    double norm1_fp_;      //!< sum of absolute fluid pressure values
    double norm1_pfv_;     //!< sum of absolute poro fluid velocity values
    double norm1_pfp_;     //!< sum of absolute poro fluid pressure values
    double norm1_ps_;      //!< sum of absolute poro structural displacements values
    double norm1_a_;       //!< sum of absolute ale displacements values

    //! iteration step
    int iter_;

    int printscreen_;    ///> print infos to standard out every printscreen_ steps
    bool printiter_;     ///> print intermediate iterations during solution
    bool printerrfile_;  ///>print infos to error file #errfile_
    FILE* errfile_;      ///> error file handle
    //! timer for solution technique
    Epetra_Time timer_;

    bool isfirsttimestep_;

    //! hydraulic conductivity (needed for coupling in case of probtype fps3i)
    double conductivity_;

   protected:
    bool islinesearch_;
    //!  flag is true if this is the first Newton iteration, false otherwise
    bool firstcall_;
  };
  // class Monolithic

}  // namespace FPSI

#endif  // FPSI_MONOLITHIC_H
