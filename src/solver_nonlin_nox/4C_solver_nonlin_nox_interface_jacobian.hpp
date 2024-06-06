/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN extension of the %NOX::Epetra jacobian
       interface.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_INTERFACE_JACOBIAN_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_INTERFACE_JACOBIAN_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Epetra_Interface_Jacobian.H>  // base class
#include <NOX_Epetra_Interface_Required.H>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace NOX
{
  namespace Nln
  {
    namespace Interface
    {
      class Jacobian : public ::NOX::Epetra::Interface::Jacobian
      {
       public:
        //! Constructor.
        Jacobian(){};

        /*! \brief Compute RHS and Jacobian at once.
         *
         *  \return TRUE if computation was successful. */
        virtual bool computeFandJacobian(
            const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac) = 0;

        /*! \brief Compute the correction system of given type.
         *
         *  \return TRUE if computation was successful. */
        virtual bool compute_correction_system(const enum CorrectionType type,
            const ::NOX::Abstract::Group& grp, const Epetra_Vector& x, Epetra_Vector& rhs,
            Epetra_Operator& jac)
        {
          return false;
        };

        virtual Teuchos::RCP<Core::LinAlg::SparseMatrix>
        calc_jacobian_contributions_from_element_level_for_ptc() = 0;
      };
    }  // end namespace Interface
  }    // end namespace Nln
}  // end namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
