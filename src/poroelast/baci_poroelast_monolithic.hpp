/*----------------------------------------------------------------------*/
/*! \file

 \brief  Basis of all monolithic poroelasticity algorithms

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_MONOLITHIC_HPP
#define FOUR_C_POROELAST_MONOLITHIC_HPP

#include "baci_config.hpp"

#include "baci_inpar_poroelast.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_poroelast_base.hpp"
#include "baci_poroelast_utils.hpp"

namespace Teuchos
{
  class Time;
}

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class SparseMatrix;
  class SparseOperator;

  class BlockSparseMatrixBase;
  class Solver;

  class Equilibration;
  enum class EquilibrationMethod;
}  // namespace CORE::LINALG

namespace POROELAST
{
  //! base class of all monolithic Poroelasticity algorithms
  class Monolithic : public PoroBase
  {
   public:
    //! create using a Epetra_Comm
    Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter);

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     ReadMesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic Poroelasticity system
    void SetupSystem() override;

    //! setup composed right hand side from field solvers
    void SetupRHS(bool firstcall = false) override;

    //! start a new time step
    void PrepareTimeStep() override;

    //! setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() { SetupSystemMatrix(*systemmatrix_); }

    //! setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat);

    //! setup equilibration of system matrix
    void SetupEquilibration();

    //! setup newton solver
    virtual void SetupNewton();


    //! build the combined dirichletbcmap
    void BuildCombinedDBCMap() override;

    //! @name Access methods for subclasses

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> Extractor() const override
    {
      return blockrowdofmap_;
    }

    //!@}

    //! @name Access methods

    //! composed system matrix
    // remove this method!
    // this method merges the block matrix when called.
    // As this is very expensive this,this method is not meant to be used any more.
    // Use BlockSystemMatrix() instead and assemble the blocks separately, if necessary.
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override;

    //! block system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      return systemmatrix_;
    }

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() override;

    //! dof row map of Structure field
    Teuchos::RCP<const Epetra_Map> DofRowMapStructure() override;

    //! dof row map of Fluid field
    Teuchos::RCP<const Epetra_Map> DofRowMapFluid() override;

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const override { return combinedDBCMap_; }

    //! right hand side vector
    Teuchos::RCP<const Epetra_Vector> RHS() override { return rhs_; }

    //! zero all entries in iterinc vector
    void ClearPoroIterinc();

    //! replaces the iterinc with poroinc
    void UpdatePoroIterinc(Teuchos::RCP<const Epetra_Vector> poroinc);

    //! iter_ += 1
    void IncrementPoroIter();

    //! FluidField()->SystemMatrix()->RangeMap()
    const Epetra_Map& FluidRangeMap();

    //! FluidField()->SystemMatrix()->DomainMap()
    const Epetra_Map& FluidDomainMap();

    //! StructureField()->SystemMatrix()->DomainMap()
    const Epetra_Map& StructureDomainMap();

    //!@}

    //! solve linear system
    void LinearSolve();

    //! create linear solver (setup of parameter lists, etc...)
    void CreateLinearSolver();

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override;

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc (with structural and fluid
    //! increment separately)
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> s_iterinc, Teuchos::RCP<const Epetra_Vector> f_iterinc);

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    //! and assemble systemmatrix and rhs-vector
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc,  //!< increment between iteration i and i+1
        bool firstiter) override;

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    //! and assemble systemmatrix and rhs-vector
    void Evaluate(Teuchos::RCP<const Epetra_Vector>
                      s_iterinc,  //!< structural increment between iteration i and i+1
        Teuchos::RCP<const Epetra_Vector>
            f_iterinc,  //!< fluid increment between iteration i and i+1
        bool firstiter) override;

    //! evaluate fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void EvaluateFields(Teuchos::RCP<const Epetra_Vector> iterinc);

    //! evaluate fields seperately at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void EvaluateFields(
        Teuchos::RCP<const Epetra_Vector> s_iterinc, Teuchos::RCP<const Epetra_Vector> f_iterinc);

    //! extract initial guess from fields
    //! returns \f$\Delta x_{n+1}^{<k>}\f$
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    //! is convergence reached of iterative solution technique?
    //! keep your fingers crossed...
    virtual bool Converged();

    //! inner newton iteration
    void Solve() override;

    //! perform one time step (setup + solve + output)
    void DoTimeStep() override;

    //! @name Output

    //! print to screen information about residual forces and displacements
    virtual void PrintNewtonIter();

    //! contains text to PrintNewtonIter
    virtual void PrintNewtonIterText(FILE* ofile  //!< output file handle
    );

    //! contains text to PrintNewtonIter
    virtual void PrintNewtonIterTextStream(std::ostringstream& oss);

    //! contains header to PrintNewtonIter
    virtual void PrintNewtonIterHeader(FILE* ofile  //!< output file handle
    );

    //! contains header to PrintNewtonIter
    virtual void PrintNewtonIterHeaderStream(std::ostringstream& oss);

    //! print statistics of converged Newton-Raphson iteration
    void PrintNewtonConv();

    //!@}

    //! finite difference check of stiffness matrix
    [[maybe_unused]] void PoroFDCheck();

    //! Evaluate no penetration condition
    void EvaluateCondition(Teuchos::RCP<CORE::LINALG::SparseOperator> Sysmat,
        POROELAST::coupltype coupltype = POROELAST::fluidfluid);

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplierAfterTimeStep() {}

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each
    //! iteration step (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplierAfterNewtonStep(
        Teuchos::RCP<const Epetra_Vector> iterinc);

    //! Setup solver for monolithic system
    bool SetupSolver() override;

    //! read restart data
    void ReadRestart(const int step) override;

   protected:
    //! Aitken
    void Aitken();

    //! Aitken Reset
    [[maybe_unused]] void AitkenReset();

    //! @name Apply current field state to system

    //! Evaluate mechanical-fluid system matrix
    virtual void ApplyStrCouplMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
    );

    //! Evaluate fluid-mechanical system matrix
    virtual void ApplyFluidCouplMatrix(
        Teuchos::RCP<CORE::LINALG::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
    );

    //!@}

    //! convergence check for Newton solver
    virtual void BuildConvergenceNorms();

    //! extract the field vectors from a given composed vector. Different for fluid and structure
    //! split
    /*!
     x is the sum of all increments up to this point.
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (e.g. velocities and pressure)
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        bool firstcall = false);

    //! @name General purpose algorithm members
    //!@{

    bool solveradapttol_;                        //!< adapt solver tolerance
    double solveradaptolbetter_;                 //!< tolerance to which is adpated ????
    Teuchos::RCP<CORE::LINALG::Solver> solver_;  //!< linear algebraic solver

    //!@}

    //! @name Printing and output
    //!@{

    int printscreen_;  //!< print infos to standard out every printscreen_ steps
    bool printiter_;   //!< print intermediate iterations during solution

    //!@}

    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

    Teuchos::RCP<Epetra_Vector> rhs_;  //!< rhs of Poroelasticity system

    //!@}

    enum INPAR::STR::DynamicType strmethodname_;  //!< enum for STR time integration

    //! @name Global matrixes

    //! block systemmatrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> systemmatrix_;

    //! structure-fluid coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_sf_;
    //! fluid-structure coupling matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> k_fs_;

    //!@}

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockrowdofmap_;

    //! dirichlet map of monolithic system
    Teuchos::RCP<Epetra_Map> combinedDBCMap_;

    //! return structure fluid coupling sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> StructFluidCouplingMatrix();

    //! return fluid structure coupling sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> FluidStructCouplingMatrix();

    //! return structure fluid coupling block sparse matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> StructFluidCouplingBlockMatrix();

    //! return fluid structure coupling block sparse matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> FluidStructCouplingBlockMatrix();


    //! @name poro-contact

    //! apply current velocity of fluid  to ContactMangager if contact problem
    void SetPoroContactStates();

    //! assemble relevant matrixes for porocontact and meshtying
    void EvalPoroMortar();

    //! flag activation poro contact no penetration condition
    bool no_penetration_;

    //!@}

    //! build block vector from field vectors, e.g. rhs, increment vector
    virtual void SetupVector(Epetra_Vector& f,  //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv,   //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> fv    //!< vector containing only fluid dofs
    );

    //! @name Iterative solution technique

    enum INPAR::POROELAST::ConvNorm normtypeinc_;   //!< convergence check for residual temperatures
    enum INPAR::POROELAST::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum INPAR::POROELAST::BinaryOp
        combincfres_;  //!< binary operator to combine temperatures and forces
    enum INPAR::POROELAST::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum INPAR::POROELAST::VectorNorm vectornorminc_;   //!< type of norm for increments

    double tolinc_;   //!< tolerance residual increment
    double tolfres_;  //!< tolerance force residual

    double tolinc_struct_;   //!< tolerance residual increment for structure displacements
    double tolfres_struct_;  //!< tolerance force residual for structure displacements

    double tolinc_velocity_;   //!< tolerance residual increment for fluid velocity field
    double tolfres_velocity_;  //!< tolerance force residual for fluid velocity field

    double tolinc_pressure_;   //!< tolerance residual increment for fluid pressure field
    double tolfres_pressure_;  //!< tolerance force residual for fluid pressure field

    double tolinc_porosity_;   //!< tolerance residual increment for porosity field
    double tolfres_porosity_;  //!< tolerance force residual for porosity field

    int itermax_;     //!< maximally permitted iterations
    int itermin_;     //!< minimally requested iteration
    double normrhs_;  //!< norm of residual forces
    double norminc_;  //!< norm of residual unknowns

    double normrhsfluidvel_;   //!< norm of residual forces (fluid velocity)
    double normincfluidvel_;   //!< norm of residual unknowns (fluid velocity)
    double normrhsfluidpres_;  //!< norm of residual forces (fluid pressure)
    double normincfluidpres_;  //!< norm of residual unknowns (fluid pressure)
    double normrhsfluid_;      //!< norm of residual forces (fluid )
    double normincfluid_;      //!< norm of residual unknowns (fluid )

    double normrhsstruct_;  //!< norm of residual forces (structure)
    double normincstruct_;  //!< norm of residual unknowns (structure)

    double normrhsporo_;  //!< norm of residual forces (porosity)
    double normincporo_;  //!< norm of residual unknowns (porosity)

    Teuchos::RCP<Teuchos::Time> timer_;  //!< timer for solution technique

    int iter_;  //!< iteration step

    //!@}

    //! @name Various global forces

    Teuchos::RCP<Epetra_Vector> iterinc_;  //!< increment between Newton steps k and k+1
    //!< \f$\Delta{x}^{<k>}_{n+1}\f$

    //!@}

    //! flag for direct solver
    bool directsolve_;

    //! @name Aitken relaxation

    //! difference of last two solutions
    // del = r^{i+1}_{n+1} = d^{i+1}_{n+1} - d^i_{n+1}
    Teuchos::RCP<Epetra_Vector> del_;
    //! difference of difference of last two pair of solutions
    // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
    Teuchos::RCP<Epetra_Vector> delhist_;
    //! Aitken factor
    double mu_;
    //!@}

    //! @name matrix equilibration

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<CORE::LINALG::Equilibration> equilibration_;

    //! equilibration method applied to system matrix
    CORE::LINALG::EquilibrationMethod equilibration_method_;
    //!@}

    //!@}
  };

}  // namespace POROELAST


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // POROELAST_MONOLITHIC_H
