/*-----------------------------------------------------------*/
/*! \file

\brief This class summarizes the functionality which all
       implicit time integration strategies share and have in
       common.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_IMPLICITBASE_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_IMPLICITBASE_H_

#include "str_timint_base.H"

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX
namespace STR
{
  namespace TIMINT
  {
    /** \brief Abstract class for all implicit based time integration strategies
     *
     *  This class is supposed to work as a connector between or a wrapper of the
     *  different implicit time integration strategies. It summarizes the functionality
     *  which all of the different implicit strategies share.
     *
     *  \author Michael Hiermeier */
    class ImplicitBase : public Base
    {
     public:
      /// constructor
      ImplicitBase();

      /// destructor
      virtual ~ImplicitBase(){};

      /// Get type of thickness scaling for thin shell structures (derived)
      virtual INPAR::STR::STC_Scale GetSTCAlgo();

      /// Get stc matrix (derived)
      virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat();

      /// Update routine for coupled problems with monolithic approach with time adaptivity
      virtual void Update(double endtime);

      /// @name Access linear system of equation via adapter (implicit only!)
      /// @{
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
      virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);
      /// @}

      /// @name XFEM/XSTR related methods
      ///@{

      /// destroy the nonlinear solver state and all internal stored state variables
      virtual bool DestroyNoxState();

      ///@}

      //! print summary after step
      virtual void PrintStep();

      /// @name Pure virtual routines (have to be implemented in the derived implicit classes)
      /// @{
      /// Take the time and integrate (time loop) (pure virtual)
      virtual int Integrate() = 0;
      /// integrate the current step (pure virtual)
      virtual int IntegrateStep() = 0;
      /// nonlinear solve (pure virtual)
      virtual INPAR::STR::ConvergenceStatus Solve() = 0;
      /// start new time step
      virtual void PrepareTimeStep() = 0;
      /// specific method for iterative staggered partitioned TSI
      virtual void PreparePartitionStep() = 0;
      /// set the state of the nox group and the global state data container
      /// see class \ref ADAPTER::StructureNew for detailed documentation.
      virtual void SetState(const Teuchos::RCP<Epetra_Vector>& x) = 0;
      /// @}

     protected:
      /// Returns the current solution group (pure virtual)
      virtual const NOX::Abstract::Group& GetSolutionGroup() const = 0;

      //! Returns the current solution group ptr
      virtual Teuchos::RCP<NOX::Abstract::Group> SolutionGroupPtr() = 0;
    };
  }  // namespace TIMINT
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_IMPLICITBASE_H_ */
