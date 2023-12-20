/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a %::NOX::Epetra::Group
       to use with nonlinear singlestep solver.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_SINGLESTEP_GROUP_H
#define BACI_SOLVER_NONLIN_NOX_SINGLESTEP_GROUP_H

#include "baci_config.H"

#include "baci_solver_nonlin_nox_constraint_interface_required.H"
#include "baci_solver_nonlin_nox_enum_lists.H"
#include "baci_solver_nonlin_nox_group.H"

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace SINGLESTEP
    {
      class Group : public NOX::NLN::Group
      {
       public:
        //! Standard constructor
        Group(Teuchos::ParameterList& printParams,  //!< printing parameters
            Teuchos::ParameterList& grpOptionParams,
            const Teuchos::RCP<::NOX::Epetra::Interface::Required>&
                i,                           //!< basically the NOXified time integrator
            const ::NOX::Epetra::Vector& x,  //!< current solution vector
            const Teuchos::RCP<::NOX::Epetra::LinearSystem>&
                linSys  //!< linear system, matrix and RHS etc.
        );

        /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
          valid shared linear system. */
        Group(const NOX::NLN::SINGLESTEP::Group& source, ::NOX::CopyType type = ::NOX::DeepCopy);

        //! generate a clone of the given object concerning the given \c CopyType
        Teuchos::RCP<::NOX::Abstract::Group> clone(::NOX::CopyType type) const override;

        //! assign operator
        ::NOX::Abstract::Group& operator=(const ::NOX::Epetra::Group& source) override;

        //! compute/update the current state variables
        void computeX(
            const NOX::NLN::SINGLESTEP::Group& grp, const ::NOX::Epetra::Vector& d, double step);
        void computeX(const ::NOX::Abstract::Group& grp, const ::NOX::Abstract::Vector& d,
            double step) override;

       private:
        //! Throw an NOX_error
        void throwError(const std::string& functionName, const std::string& errorMsg) const;
      };
    }  // namespace SINGLESTEP
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLVER_NONLIN_NOX_SINGLESTEP_GROUP_H
