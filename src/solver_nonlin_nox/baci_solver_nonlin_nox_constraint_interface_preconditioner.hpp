/*---------------------------------------------------------------------*/
/*! \file

\brief Preconditioner interface for constrained problems.



\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_CONSTRAINT_INTERFACE_PRECONDITIONER_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_CONSTRAINT_INTERFACE_PRECONDITIONER_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_enum_lists.hpp"
#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Epetra_Interface_Preconditioner.H>
#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace CORE::LINALG
{
  class Solver;
}

namespace NOX
{
  namespace NLN
  {
    namespace CONSTRAINT
    {
      namespace Interface
      {
        class Preconditioner;
      }  // namespace Interface
      // typedef
      typedef std::map<NOX::NLN::SolutionType, Teuchos::RCP<Interface::Preconditioner>>
          PrecInterfaceMap;

      namespace Interface
      {
        class Preconditioner : public ::NOX::Epetra::Interface::Preconditioner
        {
         public:
          /*! \brief Is the (CURRENT) system to solve a saddle-point system?
           *
           *  This check is supposed to return TRUE, only if the CURRENT system
           *  of equations is a saddle-point system. So in the case of inequality
           *  constraints, there is the possibility, that all constraints are
           *  inactive. In such a case the CURRENT system has no saddle-point shape
           *  and the function should return FALSE.
           *  Nevertheless, this may change during one of the following iterations!
           *
           *  \author Michael Hiermeier
           *  \date 04/2016 */
          virtual bool IsSaddlePointSystem() const = 0;

          /*! \brief Is the (CURRENT) system to solve a condensed system?
           *
           *  This check is supposed to return TRUE, only if the CURRENT system
           *  of equations involves any condensed quantities. So in the case of
           *  inequality constraints, there is the possibility, that all constraints
           *  are inactive. In such a case the CURRENT system needs no condensation
           *  and the function should return FALSE.
           *  Nevertheless, this may change during one of the following iterations!
           *
           *  \author Michael Hiermeier
           *  \date 04/2016 */
          virtual bool IsCondensedSystem() const = 0;

          //! Get necessary maps for the preconditioner.
          virtual void FillMapsForPreconditioner(
              std::vector<Teuchos::RCP<Epetra_Map>>& maps) const = 0;

          //! Get the corresponding linear solver (optional)
          virtual CORE::LINALG::Solver* GetLinearSolver() const { return nullptr; };
        };
      }  // namespace Interface
    }    // namespace CONSTRAINT
  }      // namespace NLN
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
