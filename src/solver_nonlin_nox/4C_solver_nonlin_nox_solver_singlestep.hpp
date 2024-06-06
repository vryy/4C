/*-----------------------------------------------------------*/
/*! \file

\brief Extension of Trilinos's ::NOX::Solver::SingleStep to include 4C's specific inner/outer
tests

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SOLVER_SINGLESTEP_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SOLVER_SINGLESTEP_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Solver_SingleStep.H>  // base class

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest
    namespace Inner
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace Inner
    namespace Solver
    {
      class SingleStep : public ::NOX::Solver::SingleStep
      {
       public:
        //! Constructor
        /*!
          See reset(::NOX::Abstract::Group&, ::NOX::StatusTest::Generic&, Teuchos::ParameterList&)
          for description
         */
        SingleStep(const Teuchos::RCP<::NOX::Abstract::Group>& grp,
            const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<Teuchos::ParameterList>& params);


        [[nodiscard]] ::NOX::StatusTest::StatusType getStatus() const override;

        //! Returns the ::NOX::Utils object
        [[nodiscard]] const ::NOX::Utils& GetUtils() const;

       protected:
        //! initialize additional variables after base class initialization
        void init(const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests);

        void printUpdate() override;
      };  // class SingleStep
    }     // namespace Solver
  }       // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
