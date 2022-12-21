/*-----------------------------------------------------------*/
/*! \file

\brief explicit structural time integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_TIMINT_EXPLICIT_H_
#define STR_TIMINT_EXPLICIT_H_


#include "str_timint_base.H"

namespace STR
{
  namespace TIMINT
  {
    /** \brief Explicit time integration strategy
     *
     * \author Michael Hiermeier */
    class Explicit : public Base
    {
     public:
      //! constructor
      Explicit();

      //! destructor
      virtual ~Explicit(){};

      void Setup() override;

      int Integrate() override;

      int IntegrateStep() override;

      //! Prepare time step
      void PrepareTimeStep() override;

      //!  Update state incrementally for coupled problems with monolithic approach
      void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      //!  Evaluate routine for coupled problems with monolithic approach
      void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override;

      /// set the state of the nox group and the global state data container
      /// see class \ref ADAPTER::StructureNew for detailed documentation.
      void SetState(const Teuchos::RCP<Epetra_Vector>& x) override;

      /*! Do the nonlinear solve, i.e. (multiple) corrector,
       *  for the time step. All boundary conditions have
       *  been set. */
      INPAR::STR::ConvergenceStatus Solve() override;

      //! prepare partition step
      void PreparePartitionStep() override;

      //! Update routine for coupled problems with monolithic approach with time adaptivity
      void Update(double endtime) override;

      //! print summary after step
      void PrintStep() override;

      //! Get type of thickness scaling for thin shell structures
      INPAR::STR::STC_Scale GetSTCAlgo() override;

      //! Access to scaling matrix for STC
      Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat() override;

      /// @name XFEM/XSTR related methods
      ///@{

      /// There is nothing to destroy in the explicit case
      bool DestroyNoxState() override { return true; };

      ///@}

      /// @name Access linear system of equation via adapter (implicit only!)
      ///@{
      /// initial guess of Newton's method
      Teuchos::RCP<const Epetra_Vector> InitialGuess() override;

      /// right-hand-side of Newton's method
      Teuchos::RCP<const Epetra_Vector> GetF() const override;

      /// Return reaction forces at \f$t_{n+1}\f$ (read and write)
      Teuchos::RCP<Epetra_Vector> Freact() override;

      //! Return stiffness,
      //! i.e. force residual differentiated by displacements
      //!      (structural block only)
      Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() override;

      /// Return stiffness,
      /// i.e. force residual differentiated by displacements
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override;

      ///! FixMe switch structure field to block matrix in fsi simulations
      void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps) override;

      ///@}
    };
  }  // namespace TIMINT
}  // namespace STR

#endif /* STR_TIMINT_EXPLICIT_H_ */
