/*-----------------------------------------------------------*/
/*! \file

\brief Structural non-linear %NOX::NLN solver.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_NOX_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_NOX_HPP

#include "4C_config.hpp"

#include "4C_structure_new_nln_solver_generic.hpp"  // base class

#include <NOX_StatusTest_Generic.H>

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
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
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

namespace CORE::LINALG
{
  class Solver;
}  // namespace CORE::LINALG

namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      /*! Interface to NOX as nonlinear solver in structural dynamics
       *
       */
      class Nox : public Generic
      {
       public:
        //! constructor
        Nox();


        //! derived from the base class
        void Setup() override;

        //! derived from the base class
        void Reset() override;

        //! derived from the base class
        INPAR::STR::ConvergenceStatus Solve() override;

        //! returns the outer status test object pointer
        const ::NOX::StatusTest::Generic& GetOStatusTest() const
        {
          CheckInitSetup();
          FOUR_C_ASSERT(!ostatus_.is_null(), "The outer status test object is not defined!");
          return *ostatus_;
        }

        //! returns the outer status test object pointer
        Teuchos::RCP<const NOX::NLN::INNER::StatusTest::Generic> GetIStatusPtr() const
        {
          CheckInitSetup();
          return istatus_;
        }

        //! get number of nonlinear iterations (derived)
        int GetNumNlnIterations() const override;

       protected:
        //! Reset the non-linear solver parameters and variables
        virtual void ResetParams();

        //! Convert the final nox status into a structural status
        enum INPAR::STR::ConvergenceStatus ConvertFinalStatus(
            const ::NOX::StatusTest::StatusType& finalstatus) const;

       protected:
        //! pointer to the nox nln global data container
        Teuchos::RCP<NOX::NLN::GlobalData> nlnglobaldata_;

        //! NOX non-linear solver
        Teuchos::RCP<::NOX::Solver::Generic> nlnsolver_;

       private:
        //! @name Variables which stay constant after Init() and Setup() call
        //!@{

        //!@}

        //! @name variables which are reset in each Solve() call
        //!@{

        /*! \brief NOX non-linear problem class
         *
         *  The main task is to manage the non-linear solver creation
         */
        Teuchos::RCP<NOX::NLN::Problem> problem_;

        //! linear system class
        Teuchos::RCP<::NOX::Epetra::LinearSystem> linsys_;

        //! outer status test
        Teuchos::RCP<::NOX::StatusTest::Generic> ostatus_;

        //! inner status test
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> istatus_;

        //!@}
      };  // class Nox
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
