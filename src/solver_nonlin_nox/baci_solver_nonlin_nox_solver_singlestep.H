/*-----------------------------------------------------------*/
/*! \file

\brief Extension of Trilinos's ::NOX::Solver::SingleStep to include Baci's specific inner/outer
tests

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_SOLVER_SINGLESTEP_H
#define BACI_SOLVER_NONLIN_NOX_SOLVER_SINGLESTEP_H

#include "baci_config.H"

#include "baci_solver_nonlin_nox_forward_decl.H"

#include <NOX_Solver_SingleStep.H>  // base class

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace StatusTest
    {
      enum QuantityType : int;
    }  // namespace StatusTest
    namespace INNER
    {
      namespace StatusTest
      {
        class Generic;
      }  // namespace StatusTest
    }    // namespace INNER
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
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<Teuchos::ParameterList>& params);


        [[nodiscard]] ::NOX::StatusTest::StatusType getStatus() const override;

        //! Returns the ::NOX::Utils object
        [[nodiscard]] const ::NOX::Utils& GetUtils() const;

       protected:
        //! initialize additional variables after base class initialization
        void init(const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests);

        void printUpdate() override;
      };  // class SingleStep
    }     // namespace Solver
  }       // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif /* BACI_SOLVER_NONLIN_NOX_SOLVER_SINGLESTEP_H */
