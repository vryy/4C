/*-----------------------------------------------------------*/
/*!
\file nox_nln_interface_jacobian.H

\brief %NOX::NLN extension of the %NOX::Epetra jacobian
       interface.

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_INTERFACE_JACOBIAN_H_
#define NOX_NLN_INTERFACE_JACOBIAN_H_

#include <NOX_Epetra_Interface_Jacobian.H>    // base class
#include <NOX_Epetra_Interface_Required.H>

namespace NOX {
namespace NLN {
namespace Interface{

class Jacobian : public NOX::Epetra::Interface::Jacobian
{
public:

  //! Constructor.
  Jacobian() {};

  //! Destructor.
  virtual ~Jacobian() {};

  /*! \brief Compute RHS and Jacobian at once.
   *
   * Returns true if computation was successful. */
  virtual bool computeFandJacobian(
      const Epetra_Vector& x,
      Epetra_Vector& rhs,
      Epetra_Operator& jac) = 0;
};
} // end namespace Interface
} // end namespace NLN
} // end namespace NOX


#endif /* NOX_NLN_INTERFACE_JACOBIAN_H_ */
