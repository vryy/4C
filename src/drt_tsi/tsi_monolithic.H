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
    explicit Monolithic(Epetra_Comm& comm);

    ///
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
    void TimeLoop();

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    // implemented in FSI in the "split level"
    void SetupRHS(Epetra_Vector& f,bool firstcall=false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(
      LINALG::BlockSparseMatrixBase& mat  //!< block matrix
      );

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> stepinc //!< increment between time step n and n+1
      );

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    //! Predict constant displacement/temperature, displcacement/temperature rate,
    //! i.e. the initial guess is equal to the last converged step
    //! except Dirichlet BCs
    void PredictConstValueRate();

  protected:

    /// Outer Iteration loop
    void NewtonFull();

    /// convergence check for Newton solver
    bool ConvergenceCheck(int itnum, int itmax, double ittol);

    /// extract the three field vectors from a given composed vector
    /*!
      x is the sum of all increments up to this point.
      \param x  (i) composed vector that contains all field vectors
      \param sx (o) displacements
      \param tx (o) temperatures
      */
    virtual void ExtractFieldVectors(
      Teuchos::RCP<const Epetra_Vector> x,
      Teuchos::RCP<const Epetra_Vector>& sx,
      Teuchos::RCP<const Epetra_Vector>& tx
      );

    //! @name Access methods for subclasses

    /// full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const
    {
      return blockrowdofmap_.FullMap();
    }

    /// set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
    */
    void SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

    /// extractor to communicate between full monolithic map and block maps
    const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }

    //@}

    //! @name General purpose algorithm members
    //@{

    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver

    //@}

  private:

    /// dof row map splitted in (field) blocks
    LINALG::MultiMapExtractor blockrowdofmap_;

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,  ///< full RHS vector of TSI
      Teuchos::RCP<const Epetra_Vector> sv,  ///< structural RHS vector
      Teuchos::RCP<const Epetra_Vector> tv  ///< thermal RHS vector
      );

    // TODO: type of systemmatrix is correct??
    /// block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > systemmatrix_;
    //    Teuchos::RCP<TSIBlockMatrix> systemmatrix_;

//    //! @name Iterative solution technique
//    int iter_;  //!< iteration step
//    double normfres_;  //!< norm of residual forces
//    double normtempi_;  //!< norm of residual unknowns

    //@}

    /// needed for ConvergenceCheck()
    /// maximum iteration steps
    int itmax_;
    /// temperature tolerance
    double ittol_;
    /// temperature increment of the outer loop
    Teuchos::RCP<Epetra_Vector> tempincnp_;  //!< \f$\Delta{T}^{<k>}_{n+1}\f$
    /// displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;  //!< \f$\Delta{d}^{<k>}_{n+1}\f$

  }; // Monolithic


  /// Base class for all TSI block preconditioning matrices
  class TSIBlockMatrix : public LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
  {
  public:

    /// construction
    TSIBlockMatrix(
      const LINALG::MultiMapExtractor& maps,
      ADAPTER::Structure& structure,
      ADAPTER::Thermo& thermo
      );

    /// setup of block matrix
    virtual void SetupTSIBlockMatrix();

    /** \name Attribute access functions */
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

  protected:

    /** \name Field solver objects (implemented in class BlockPreconditioningMatrix)*/
    //@{

    Teuchos::RCP<LINALG::Solver> structuresolver_;
    Teuchos::RCP<LINALG::Solver> thermosolver_;

    //@}

  };  // TSIBlockMatrix

} // namespace TSI

/*----------------------------------------------------------------------*/
#endif  // TSI_MONOLITHIC_H
#endif  // CCADISCRET
