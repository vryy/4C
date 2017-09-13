/*----------------------------------------------------------------------*/
/*!
\file tsi_monolithic.H

\maintainer Alexander Seitz

\level 1

\brief Basis of all monolithic TSI algorithms that perform a coupling between
       the structure field equation and temperature field equations

*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 11/10 |
 *----------------------------------------------------------------------*/
#ifndef TSI_MONOLITHIC_H
#define TSI_MONOLITHIC_H


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/10 |
 *----------------------------------------------------------------------*/
#include "tsi_algorithm.H"
#include "../drt_inpar/inpar_tsi.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_thermo.H"

#include <Epetra_Time.h>

/*----------------------------------------------------------------------*
 |                                                           dano 11/10 |
 *----------------------------------------------------------------------*/
  // forward declarations
namespace LINALG
{
  class SparseMatrix;
  class MultiMapExtractor;

  class BlockSparseMatrixBase;
  class Solver;
}

namespace DRT
{
  namespace UTILS
  {
    class LocsysManager;
  }
}

namespace ADAPTER
{
 class MortarVolCoupl;
}

//! TSI: Thermo-Structure Interaction
namespace TSI
{
  namespace UTILS
  {
    // forward declaration of clone strategy
    class ThermoStructureCloneStrategy;
  }

  //! monolithic TSI algorithm
  //!
  //!  Base class of TSI algorithms. Derives from StructureBaseAlgorithm and
  //!  ThermoBaseAlgorithm with temperature field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, thermo. In that order.
  //!
  //!  \note There is the Algorithm class for general purpose TSI algorithms.
  //!  This simplifies the monolithic implementation.
  //!
  //!  \author u.kue
  //!  \date 02/08
  class Monolithic : public Algorithm
  {
  public:

    explicit Monolithic(
      const Epetra_Comm& comm,
      const Teuchos::ParameterList& sdynparams
      );

    //! virtual destructor to support polymorph destruction
    virtual ~Monolithic();


    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling
    2.) get maps for all blocks in the system (and for the whole system as well)
        create combined map
    3.) create system matrix

    </pre>

    \note We want to do this setup after reading the restart information, not
    directly in the constructor. This is necessary since during restart (if
    ReadMesh is called), the dofmaps for the blocks might get invalid.
    */
    //! Setup the monolithic TSI system
    virtual void SetupSystem();

    /// non-linear solve, i.e. (multiple) corrector
    virtual void Solve();

    //! outer level TSI time loop
    virtual void TimeLoop();

    //! read restart data
    virtual void ReadRestart(
      int step  //!< step number where the calculation is continued
      );

    //! @name Apply current field state to system

    //! setup composed right hand side from field solvers
    void SetupRHS();

    //! setup composed system matrix from field solvers
    void SetupSystemMatrix();

    //! composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const
    {
      return systemmatrix_;
    }

    //! solve linear TSI system
    void LinearSolve();

    //! create linear solver (setup of parameter lists, etc...)
    void CreateLinearSolver();

    //! Evaluate mechanical-thermal system matrix
    void ApplyStrCouplMatrix(
      Teuchos::RCP<LINALG::SparseMatrix> k_st  //!< mechanical-thermal stiffness matrix
      );

    //! Evaluate thermal-mechanical system matrix
    void ApplyThrCouplMatrix(
      Teuchos::RCP<LINALG::SparseMatrix> k_ts  //!< thermal-mechanical tangent matrix
      );

    //! Evaluate thermal-mechanical system matrix for geonln + heat convection BC
    void ApplyThrCouplMatrix_ConvBC(
      Teuchos::RCP<LINALG::SparseMatrix> k_ts  //!< thermal-mechanical tangent matrix
      );

    //@}

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
      Teuchos::RCP<Epetra_Vector> stepinc //!< increment between time step n and n+1
      );

    //! extract initial guess from fields
    //! returns \f$\Delta x_{n+1}^{<k>}\f$
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    //! is convergence reached of iterative solution technique?
    //! keep your fingers crossed...
    //! \author lw (originally in STR) \date 12/07
    bool Converged();

    //! outer iteration loop
    void NewtonFull();

    //! apply DBC to all blocks
    void ApplyDBC();

    //! do pseudo-transient continuation nonlinear iteration
    //!
    //! Pseudo-transient continuation is a variant of a full newton which has a
    //! larger convergence radius than newton and is therefore more stable
    //! and/or can do larger time steps
    //!
    //! originally by mwgee for structural analysis \date 03/12
    void PTC();

    //! @name Output

    //! print to screen information about residual forces and displacements
    //! \author lw (originally in STR) \date 12/07
    void PrintNewtonIter();

    //! contains text to PrintNewtonIter
    //! \author lw (originally in STR) \date 12/07
    void PrintNewtonIterText(
      FILE* ofile  //!< output file handle
      );

    //! contains header to PrintNewtonIter
    //! \author lw (originally) \date 12/07
    void PrintNewtonIterHeader(
      FILE* ofile  //!< output file handle
      );

    //! print statistics of converged Newton-Raphson iteration
    void PrintNewtonConv();

    //! Determine norm of force residual
    double CalculateVectorNorm(
      const enum INPAR::TSI::VectorNorm norm,  //!< norm to use
      const Teuchos::RCP<Epetra_Vector> vect  //!< the vector of interest
      );

    //@}

    //! apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    //! undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(
      LINALG::BlockSparseMatrixBase& mat,
      Epetra_Vector& x,
      Epetra_Vector& b
      );

  protected:

    //! @name Time loop building blocks

    //! start a new time step
    virtual void PrepareTimeStep();

    //! calculate stresses, strains, energies
    virtual void PrepareOutput();
    //@}

    //! convergence check for Newton solver
    bool ConvergenceCheck(int itnum, int itmax, double ittol);

    //! extract the three field vectors from a given composed vector
    /*!
      x is the sum of all increments up to this point.
      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural vector (e.g. displacements)
      \param tx (o) thermal vector (e.g. temperatures)
      */
    virtual void ExtractFieldVectors(
      Teuchos::RCP<Epetra_Vector> x,
      Teuchos::RCP<Epetra_Vector>& sx,
      Teuchos::RCP<Epetra_Vector>& tx
      );

    //! @name Access methods for subclasses

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
    */
    void SetDofRowMaps();

    //! combined DBC map
    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    //! extractor to communicate between full monolithic map and block maps
    Teuchos::RCP<LINALG::MultiMapExtractor> Extractor() const { return blockrowdofmap_; }

    //! setup list with default parameters
    void SetDefaultParameters();

    //! recover structural and thermal Lagrange multipliers
    // this takes into account the dependence on off diagonal blocks
    void RecoverStructThermLM();

    //@}

    //! @name General purpose algorithm members
    //@{

    bool solveradapttol_;  //!< adapt solver tolerance
    double solveradaptolbetter_;  //!< tolerance to which is adpated ????
    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver

    //@}

    //! @name Printing and output
    //@{

    bool printiter_;  //!< print intermediate iterations during solution
    bool printerrfile_;  //!< print infos to error file #errfile_
    FILE* errfile_;  //!< error file handle

    //! calculate nodal values (displacements, temperatures, reaction forces) at
    //! specific nodes used for validation of implementation with literature
    //! here: validation of thermoplasticity with e.g. Simo and Miehe (1992)
    void CalculateNeckingTSIResults();

    //@}

    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}

    //! enum for STR time integartion
    enum INPAR::STR::DynamicType strmethodname_;


    //! apply temperature state on structure discretization
    virtual void ApplyThermoCouplingState(Teuchos::RCP<const Epetra_Vector> temp,
                                  Teuchos::RCP<const Epetra_Vector> temp_res=Teuchos::null);

    //! apply structural displacements and velocities on thermo discretization
    virtual void ApplyStructCouplingState(Teuchos::RCP<const Epetra_Vector> disp,
                                  Teuchos::RCP<const Epetra_Vector> vel);

  private:

    //! if just rho_inf is specified for genAlpha, the other parameters in the global parameter
    //! list need to be adapted accordingly
    void FixTimeIntegrationParams();

    const Teuchos::ParameterList& tsidyn_;  //!< TSI dynamic parameter list
    const Teuchos::ParameterList& tsidynmono_;  //!< monolithic TSI dynamic parameter list

    //! dofrowmap splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> blockrowdofmap_;

    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(
      Epetra_Vector& f,  //!< vector of length of all dofs
      Teuchos::RCP<const Epetra_Vector> sv,  //!< vector containing only structural dofs
      Teuchos::RCP<const Epetra_Vector> tv  //!< vector containing only thermal dofs
      );

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

    //! off diagonal matrixes
    Teuchos::RCP<LINALG::SparseMatrix> k_st_;
    Teuchos::RCP<LINALG::SparseMatrix> k_ts_;

    bool merge_tsi_blockmatrix_;  //!< bool whether TSI block matrix is merged

    //! @name iterative solution technique

    enum INPAR::TSI::NlnSolTech soltech_;  //!< kind of iteration technique or
                                           //!< nonlinear solution technique

    enum INPAR::TSI::ConvNorm normtypeinc_;  //!< convergence check for increments
    enum INPAR::TSI::ConvNorm normtyperhs_;  //!< convergence check for residual forces
    enum INPAR::STR::ConvNorm normtypedisi_;  //!< convergence check for residual displacements
    enum INPAR::STR::ConvNorm normtypestrrhs_;  //!< convergence check for residual forces
    enum INPAR::THR::ConvNorm normtypetempi_;  //!< convergence check for residual temperatures
    enum INPAR::THR::ConvNorm normtypethrrhs_;  //!< convergence check for residual thermal forces

    enum INPAR::TSI::BinaryOp combincrhs_;  //!< binary operator to combine increments and forces

    enum INPAR::TSI::VectorNorm iternorm_;  //!< vector norm to check TSI values with
    enum INPAR::TSI::VectorNorm iternormstr_;  //!< vector norm to check structural values with
    enum INPAR::TSI::VectorNorm iternormthr_;  //!< vector norm to check thermal values with

    double tolinc_;  //!< tolerance for increment
    double tolrhs_;  //!< tolerance for rhs
    double toldisi_;  //!< tolerance for displacement increments
    double tolstrrhs_;  //!< tolerance for structural rhs
    double toltempi_;  //!< tolerance for temperature increments
    double tolthrrhs_;  //!< tolerance for thermal rhs

    double normrhs_;  //!< norm of residual forces
    double normrhsiter0_; //!< norm of residual force of 1st iteration
    double norminc_;  //!< norm of residual unknowns
    double norminciter0_;  //!< norm of residual unknowns of 1st iteration
    double normdisi_;  //!< norm of residual displacements
    double normdisiiter0_;  //!< norm of residual displacements of 1st iteration
    double normstrrhs_;  //!< norm of structural residual forces
    double normstrrhsiter0_;  //!< norm of structural residual forces of 1st iteration
    double normtempi_;  //!< norm of residual temperatures
    double normtempiiter0_;  //!< norm of residual temperatures of 1st iteration
    double normthrrhs_;  //!< norm of thermal residual forces
    double normthrrhsiter0_;  //!< norm of thermal residual forces of 1st iteration

    int iter_;  //!< iteration step
    int itermax_;  //!< maximally permitted iterations
    int itermin_;  //!< minimally requested iteration

    const Teuchos::ParameterList& sdyn_;  //!< structural dynamic parameter list

    Epetra_Time timernewton_; //!< timer for measurement of solution time of newton iterations
    double dtsolve_; //!< linear solver time
    double dtcmt_; //!< contact evaluation time
    //@}

    //! @name Pseudo-transient continuation parameters

    double ptcdt_;  //!< pseudo time step size for PTC
    double dti_;  //!< scaling factor for PTC (initially 1/ptcdt_, then adapted)

    //@}

    //! @name Various global forces

    //! rhs of TSI system
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! increment between Newton steps k and k+1 \f$\Delta{x}^{<k>}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> iterinc_;

    //! global velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> vel_;

    //! Dirichlet BCs with local co-ordinate system
    Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

    //@}

    //! @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;  //!< sum of absolute values of the rows of the structural block
    Teuchos::RCP<Epetra_Vector> scolsum_;  //!< sum of absolute values of the column of the structural block
    Teuchos::RCP<Epetra_Vector> trowsum_;  //!< sum of absolute values of the rows of the thermal block
    Teuchos::RCP<Epetra_Vector> tcolsum_;  //!< sum of absolute values of the column of the thermal block

    //@}

  }; // Monolithic

} // namespace TSI


/*----------------------------------------------------------------------*/
#endif  // TSI_MONOLITHIC_H
