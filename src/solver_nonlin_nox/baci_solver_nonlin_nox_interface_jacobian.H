/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN extension of the %NOX::Epetra jacobian
       interface.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_INTERFACE_JACOBIAN_H
#define BACI_SOLVER_NONLIN_NOX_INTERFACE_JACOBIAN_H

#include "baci_solver_nonlin_nox_enum_lists.H"

#include <NOX_Epetra_Interface_Jacobian.H>  // base class
#include <NOX_Epetra_Interface_Required.H>
#include <Teuchos_RCPDecl.hpp>

// forward declaration
namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
  namespace NLN
  {
    namespace Interface
    {
      class Jacobian : public NOX::Epetra::Interface::Jacobian
      {
       public:
        //! Constructor.
        Jacobian(){};

        //! Destructor.
        ~Jacobian() override{};

        /*! \brief Compute RHS and Jacobian at once.
         *
         *  \return TRUE if computation was successful. */
        virtual bool computeFandJacobian(
            const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac) = 0;

        /*! \brief Compute the correction system of given type.
         *
         *  \return TRUE if computation was successful. */
        virtual bool computeCorrectionSystem(const enum CorrectionType type,
            const NOX::Abstract::Group& grp, const Epetra_Vector& x, Epetra_Vector& rhs,
            Epetra_Operator& jac)
        {
          return false;
        };

        virtual Teuchos::RCP<CORE::LINALG::SparseMatrix>
        CalcJacobianContributionsFromElementLevelForPTC() = 0;
      };
    }  // end namespace Interface
  }    // end namespace NLN
}  // end namespace NOX


#endif  // SOLVER_NONLIN_NOX_INTERFACE_JACOBIAN_H
