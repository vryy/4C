/*-----------------------------------------------------------*/
/*!
\file str_timint_explicit.H

\maintainer Philipp Farah

\date Aug 12, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_TIMINT_EXPLICIT_H_
#define STR_TIMINT_EXPLICIT_H_


#include "str_timint_base.H"

namespace STR {
namespace TIMINT {
  /** \brief Explicit time integration strategy
   *
   * \author Michael Hiermeier */
  class Explicit : public Base
  {
  public:
    //! constructor
    Explicit();

    //! destructor
    virtual ~Explicit() {};

    virtual void Setup();

    virtual int Integrate();

    virtual int IntegrateStep();

    //! Prepare time step
    virtual void PrepareTimeStep();


    //!  Evaluate routine for coupled problems with monolithic approach
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc);

    /*! Do the nonlinear solve, i.e. (multiple) corrector,
     *  for the time step. All boundary conditions have
     *  been set. */
    virtual INPAR::STR::ConvergenceStatus Solve();

    /* Linear structure solve with just an interface load */
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear();

    //! prepare partiton step
    virtual void PreparePartitionStep();

    //! Update routine for coupled problems with monolithic approach with time adaptivity
    virtual void Update(double endtime);

    //! print summary after step
    virtual void PrintStep();

    //! Get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo();

    //! Access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat();

    /// @name Access linear system of equation via adapter (implicit only!)
    ///@{
    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess();

    /// right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> GetF() const;

    /// Return reaction forces at \f$t_{n+1}\f$ (read and write)
    virtual Teuchos::RCP<Epetra_Vector> Freact();

    //! Return stiffness,
    //! i.e. force residual differentiated by displacements
    //!      (structural block only)
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// Return stiffness,
    /// i.e. force residual differentiated by displacements
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    ///! FixMe switch structure field to block matrix in fsi simulations
    virtual void UseBlockMatrix(
        Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);

    ///@}
  };
} // namespace TimInt
} // namespace STR

#endif /* STR_TIMINT_EXPLICIT_H_ */
