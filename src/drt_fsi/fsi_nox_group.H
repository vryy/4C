#ifdef CCADISCRET

#ifndef FSI_NOX_GROUP_H
#define FSI_NOX_GROUP_H

#include <NOX_Epetra_Group.H>
#include "fsi_monolithic.H"

namespace FSI
{

  /// Special NOX group that always sets Jacobian and RHS at the same time.
  class NOXGroup : public NOX::Epetra::Group
  {
  public:

    NOXGroup(Monolithic& mfsi,
             Teuchos::ParameterList& printParams,
	     const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
	     const NOX::Epetra::Vector& x,
	     const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys);

    virtual NOX::Abstract::Group::ReturnType computeF();

    virtual NOX::Abstract::Group::ReturnType computeJacobian();

  private:

    Monolithic& mfsi_;
  };
}

#endif
#endif
