/*-----------------------------------------------------------*/
/*!
\file nox_nln_interface_jacobian.H

\brief %NOX::NLN extension of the %NOX::Epetra jacobian
       interface.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INTERFACE_JACOBIAN_H_
#define NOX_NLN_INTERFACE_JACOBIAN_H_

#include <NOX_Epetra_Interface_Jacobian.H>  // base class
#include <NOX_Epetra_Interface_Required.H>

// forward declaration
namespace LINALG
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
        virtual ~Jacobian(){};

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

        virtual Teuchos::RCP<LINALG::SparseMatrix>
        CalcJacobianContributionsFromElementLevelForPTC() = 0;
      };
    }  // end namespace Interface
  }    // end namespace NLN
}  // end namespace NOX


#endif /* NOX_NLN_INTERFACE_JACOBIAN_H_ */
