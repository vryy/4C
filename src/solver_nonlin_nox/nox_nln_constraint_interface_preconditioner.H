/*---------------------------------------------------------------------*/
/*!
\file nox_nln_constraint_interface_preconditioner.H

\brief Preconditioner interface for constrained problems.

\maintainer Michael Hiermeier

\date Apr 22, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_CONSTRAINT_INTERFACE_PRECONDITIONER_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_CONSTRAINT_INTERFACE_PRECONDITIONER_H_


#include <NOX_Epetra_Interface_Preconditioner.H>
#include "nox_nln_enum_lists.H"
#include <Teuchos_RCP.hpp>

// forward declarations...
class Epetra_Map;

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
        class Preconditioner : public NOX::Epetra::Interface::Preconditioner
        {
         public:
          //! destructor
          virtual ~Preconditioner(){};

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
        };
      }  // namespace Interface
    }    // namespace CONSTRAINT
  }      // namespace NLN
}  // namespace NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_CONSTRAINT_INTERFACE_PRECONDITIONER_H_ */
