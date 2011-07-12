/*----------------------------------------------------------------------*/
/*!
\file tsi_monolithic.H

\brief Basis of all monolithic TSI algorithms that perform a coupling between
       the structure field equation and temperature field equations

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 11/10 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TSI_MONOLITHIC_H
#define TSI_MONOLITHIC_H

/*----------------------------------------------------------------------*
 | headers                                                   dano 11/10 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_inpar/inpar_tsi.H"

/*----------------------------------------------------------------------*
 |                                                           dano 11/10 |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;

  class BlockSparseMatrixBase;
  class Solver;
}

/// TSI: Thermo-Structure Interaction
namespace TSI
{
  /// monolithic TSI algorithm base
  //!
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
  class MonolithicBase : public ADAPTER::AlgorithmBase,
                         public ADAPTER::StructureBaseAlgorithm,
                         public ADAPTER::ThermoBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit MonolithicBase(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~MonolithicBase();

    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

  protected:

    //! @name Time loop building blocks

    /// start a new time step
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// write output
    virtual void Output();

    //@}

  }; // MonolithicBase


  /// base class of all monolithic TSI algorithms
  /*!

    Monolithic TSI is a Netwon solver on a block matrix with field blocks.

    \author u.kue
    \date 02/08
   */
  class Monolithic : public MonolithicBase
  {

  public:
    explicit Monolithic(
      Epetra_Comm& comm,
      const Teuchos::ParameterList& sdynparams
      );

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
    // Setup the monolithic TSI system
    virtual void SetupSystem();

    /// outer level TSI time loop
    void TimeLoop(
      const Teuchos::ParameterList& sdynparams  //!< structural dynamic input parameters
      );

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    // implemented in FSI in the "split level"
    void SetupRHS();

    /// setup composed system matrix from field solvers
    void SetupSystemMatrix(const Teuchos::ParameterList& sdynparams);

    /// composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const
    {
      return systemmatrix_;
    }

    /// solve linear TSI system
    void LinearSolve();

    //! Evaluate mechanical-thermal system matrix
    void ApplyStrCouplMatrix(
      Teuchos::RCP<LINALG::SparseMatrix> k_st,  //!< mechanical-thermal stiffness matrix
      const Teuchos::ParameterList& sdynparams
      );

    //! Evaluate thermal-mechanical system matrix
    void ApplyThrCouplMatrix(
      Teuchos::RCP<LINALG::SparseMatrix> k_ts,  //!< thermal-mechanical tangent matrix
      const Teuchos::ParameterList& sdynparams
      );

    //@}

    /// evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc //!< increment between time step n and n+1
      );

    /// extract initial guess from fields
    //! returns \f$\Delta x_{n+1}^{<k>}\f$
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    //! is convergence reached of iterative solution technique?
    //! keep your fingers crossed...
    //! \author lw  \date 12/07
    bool Converged();

    /// outer iteration loop
    void NewtonFull(const Teuchos::ParameterList& sdyn);

    //! @name Output

    //! print to screen information about residual forces and displacements
    //! \author lw (originally) \date 12/07
    void PrintNewtonIter();

    //! contains text to PrintNewtonIter
    //! \author lw (originally) \date 12/07
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

    //@}

  protected:

    //! @name Access methods

    //! velocity calculation given the displacements (like in FSI)
    Teuchos::RCP<Epetra_Vector> CalcVelocity(
      Teuchos::RCP<const Epetra_Vector> dispnp
      );

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
      Teuchos::RCP<const Epetra_Vector> x,
      Teuchos::RCP<const Epetra_Vector>& sx,
      Teuchos::RCP<const Epetra_Vector>& tx
      );

    //! @name Access methods for subclasses

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const
    {
      return blockrowdofmap_.FullMap();
    }

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
    */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    //! combined DBC map
    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    //! extractor to communicate between full monolithic map and block maps
    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }
    
    //! apply contact to off diagonal block k_st
    /*!
    Within the condensation process of Lagrange multipliers, the structural 
    LM are expressed in terms of the displacements, and, in the TSI case,
    also in terms of the temperatures. This results in modification in the
    off diagonal block k_st
    */
    void ApplyStructContact(Teuchos::RCP<LINALG::SparseMatrix>& k_st);

    //! apply contact to off diagonal block k_ts
    /*!
    The contribution to the energy equation form thermal contact and the 
    thermal contact condition both depend on displacement dofs (deformed 
    configuration). This leads to additional entries in k_ts. Additionally,
    the condensation process of thermal LM leads to contributions to k_ts.
    */
    void ApplyThermContact(Teuchos::RCP<LINALG::SparseMatrix>& k_ts );
    
    //! entries from linearized mortar matrices D and M in energy
    //! equation (with respect to displacements)
    void AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                       LINALG::SparseMatrix& linmglobal);
    
    //! entries from linearized thermal contact condition
    //! (with respect to displacements)
    void AssembleThermContCondition(LINALG::SparseMatrix& lindisglobal);

    //! recover structural and thermal Lagrange multipliers
    // this calls the separate fields
    // plus takes into account the dependence on off diagonal blocks
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

    bool printscreen_;  //!< print infos to standard out
    bool printiter_;  //!< print intermediate iterations during solution
    bool printerrfile_;  //!< print infos to error file #errfile_
    FILE* errfile_;  //!< error file handle

    //@}

    //! @name Global vectors
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    //@}

    //! enum for STR time integartion
    enum INPAR::STR::DynamicType strmethodname_;


  private:

    //! dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(
      Epetra_Vector &f,  //!< vector of length of all dofs
      Teuchos::RCP<const Epetra_Vector> sv,  //!< vector containing only structural dofs
      Teuchos::RCP<const Epetra_Vector> tv  //!< vector containing only thermal dofs
      );

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

    //! @name Iterative solution technique

    enum INPAR::TSI::ConvNorm normtypeinc_;  //!< convergence check for residual temperatures
    enum INPAR::TSI::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum INPAR::TSI::BinaryOp combincfres_;  //!< binary operator to combine temperatures and forces

    double tolinc_;  //!< tolerance residual temperatures
    double tolfres_;  //!< tolerance force residual
    int iter_;  //!< iteration step
    int itermax_;  //!< maximally permitted iterations
    int itermin_;  //!< minimally requested iteration
    double normrhs_;  //!< norm of residual forces
    double norminc_;  //!< norm of residual unknowns

    //@}

    //! @name Various global forces

    Teuchos::RCP<Epetra_Vector> rhs_;  //!< rhs of TSI system
    Teuchos::RCP<Epetra_Vector> freact_;  //!< reaction force
    Teuchos::RCP<Epetra_Vector> iterinc_;  //!< increment between Newton steps k and k+1
                                           //!< \f$\Delta{x}^{<k>}_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> veln_;  //!< global velocities
                                        //!< \f${V}_{n+1}\f$
                                        //!< at \f$t_{n+1}\f$
    
    // contact manager
    RCP<MORTAR::ManagerBase> cmtman_;
    
    // thermo contact manager
    Teuchos::RCP<THR::ThermoContactMan> thermcontman_;
    
    RCP<Epetra_Vector> zstruct_ ;
    RCP<Epetra_Vector> ztherm_;
    RCP<LINALG::SparseMatrix> kst_;  // stiffness block K_ss (needed for LM)
    RCP<LINALG::SparseMatrix> invd_;  // stiffness block K_ss (needed for LM)
    RCP<LINALG::SparseMatrix> kts_;  // stiffness block K_ss (needed for LM)
    RCP<LINALG::SparseMatrix> invdtherm_;  // stiffness block K_ss (needed for LM)
    
    
    //@}

  }; // Monolithic

} // namespace TSI


/*----------------------------------------------------------------------*/
#endif  // TSI_MONOLITHIC_H
#endif  // CCADISCRET
