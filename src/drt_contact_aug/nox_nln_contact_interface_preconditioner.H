/*-----------------------------------------------------------*/
/*!
\file nox_nln_contact_interface_preconditioner.H

\maintainer Michael Hiermeier

\date Jul 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONTACT_INTERFACE_PRECONDITIONER_H_
#define NOX_NLN_CONTACT_INTERFACE_PRECONDITIONER_H_

#include "NOX_Epetra_Interface_Preconditioner.H"

// forward declaration
class Epetra_Vector;
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}
namespace MORTAR {
  class StrategyBase;
}

namespace NOX {
namespace NLN {
namespace CONTACT {
namespace Interface {
  class Preconditioner : public NOX::Epetra::Interface::Preconditioner
  {
  public:
    //! Constructor
    Preconditioner() {};

    //! Destructor
    virtual ~Preconditioner() {};

    //! derived and unused at the moment
    virtual bool computePreconditioner(const Epetra_Vector& x,
                       Epetra_Operator& M,
                       Teuchos::ParameterList* precParams = 0) {return false;};

    //! return bool indicating if contact is defined
    virtual bool HaveContact() = 0;

    //! return bool indicating if meshtying is defined
    virtual bool HaveMeshtying() = 0;

    virtual MORTAR::StrategyBase& GetStrategy() = 0;
  };
} // namespace Interface
} // namespace CONTACT
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_CONTACT_INTERFACE_PRECONDITIONER_H_ */
