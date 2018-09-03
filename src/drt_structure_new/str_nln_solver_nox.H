/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_nox.H

\brief Structural non-linear %NOX::NLN solver.

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_NOX_H_
#define STR_NLN_SOLVER_NOX_H_

#include "str_nln_solver_generic.H"  // base class
#include "NOX_StatusTest_Generic.H"

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
  namespace Epetra
  {
    class LinearSystem;
  }  // namespace Epetra
  namespace NLN
  {
    class GlobalData;
    class Problem;
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER
  }      // namespace NLN
}  // namespace NOX

namespace LINALG
{
  class Solver;
}  // namespace LINALG

namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      class Nox : public Generic
      {
       public:
        //! constructor
        Nox();

        //! destructor
        virtual ~Nox(){};

        //! derived from the base class
        virtual void Setup();

        //! derived from the base class
        virtual void Reset();

        //! derived from the base class
        virtual INPAR::STR::ConvergenceStatus Solve();

        //! returns the outer status test object pointer
        const NOX::StatusTest::Generic& GetOStatusTest() const
        {
          CheckInitSetup();
          if (ostatus_.is_null()) dserror("The outer status test object is not defined!");
          return *ostatus_;
        }

        //! returns the outer status test object pointer
        Teuchos::RCP<const NOX::NLN::INNER::StatusTest::Generic> GetIStatusPtr() const
        {
          CheckInitSetup();
          return istatus_;
        }

        //! get number of nonlinear iterations (derived)
        virtual int GetNumNlnIterations() const;

       protected:
        //! Reset the non-linear solver parameters and variables
        void ResetParams();

        //! Convert the final nox status into a structural status
        enum INPAR::STR::ConvergenceStatus ConvertFinalStatus(
            const NOX::StatusTest::StatusType& finalstatus) const;

       private:
        //! @name Variables which stay constant after Init() and Setup() call
        //!@{
        //! pointer to the nox nln global data container
        Teuchos::RCP<NOX::NLN::GlobalData> nlnglobaldata_;
        //!@}

        //! @name variables which are reset in each Solve() call
        //!@{
        //! NOX non-linear problem class
        //! The main task is to manage the non-linear solver creation
        Teuchos::RCP<NOX::NLN::Problem> problem_;

        //! linear system class
        Teuchos::RCP<NOX::Epetra::LinearSystem> linsys_;

        //! outer status test
        Teuchos::RCP<NOX::StatusTest::Generic> ostatus_;

        //! inner status test
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> istatus_;

        //! NOX non-linear solver
        Teuchos::RCP<NOX::Solver::Generic> nlnsolver_;
        //!@}
      };  // class Nox
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


#endif /* STR_NLN_SOLVER_NOX_H_ */
