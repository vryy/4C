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
  namespace Nln
  {
    class GlobalData;
    class Problem;
    namespace Inner
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace Inner
  }      // namespace Nln
}  // namespace NOX

namespace Core::LinAlg
{
  class Solver;
}  // namespace Core::LinAlg

namespace STR
{
  namespace Nln
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
        void setup() override;

        //! derived from the base class
        void reset() override;

        //! derived from the base class
        Inpar::STR::ConvergenceStatus Solve() override;

        //! returns the outer status test object pointer
        const ::NOX::StatusTest::Generic& get_outer_status_test() const
        {
          check_init_setup();
          FOUR_C_ASSERT(!ostatus_.is_null(), "The outer status test object is not defined!");
          return *ostatus_;
        }

        //! returns the outer status test object pointer
        Teuchos::RCP<const NOX::Nln::Inner::StatusTest::Generic> get_inner_status_ptr() const
        {
          check_init_setup();
          return istatus_;
        }

        //! get number of nonlinear iterations (derived)
        int get_num_nln_iterations() const override;

       protected:
        //! Reset the non-linear solver parameters and variables
        virtual void reset_params();

        //! Convert the final nox status into a structural status
        enum Inpar::STR::ConvergenceStatus convert_final_status(
            const ::NOX::StatusTest::StatusType& finalstatus) const;

       protected:
        //! pointer to the nox nln global data container
        Teuchos::RCP<NOX::Nln::GlobalData> nlnglobaldata_;

        //! NOX non-linear solver
        Teuchos::RCP<::NOX::Solver::Generic> nlnsolver_;

       private:
        //! @name Variables which stay constant after init() and setup() call
        //!@{

        //!@}

        //! @name variables which are reset in each Solve() call
        //!@{

        /*! \brief NOX non-linear problem class
         *
         *  The main task is to manage the non-linear solver creation
         */
        Teuchos::RCP<NOX::Nln::Problem> problem_;

        //! linear system class
        Teuchos::RCP<::NOX::Epetra::LinearSystem> linsys_;

        //! outer status test
        Teuchos::RCP<::NOX::StatusTest::Generic> ostatus_;

        //! inner status test
        Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> istatus_;

        //!@}
      };  // class Nox
    }     // namespace SOLVER
  }       // namespace Nln
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
