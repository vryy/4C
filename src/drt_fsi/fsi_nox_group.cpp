#ifdef CCADISCRET

#include "fsi_nox_group.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::NOXGroup::NOXGroup(Monolithic& mfsi,
                        Teuchos::ParameterList& printParams,
                        const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
                        const NOX::Epetra::Vector& x,
                        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys)
  : NOX::Epetra::Group(printParams,i,x,linSys),
    mfsi_(mfsi)
{
  // we know we already have the first linear system calculated

  mfsi_.SetupRHS(RHSVector.getEpetraVector());
  mfsi_.SetupSystemMatrix(*mfsi_.SystemMatrix());

  isValidJacobian = true;
  isValidRHS = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType FSI::NOXGroup::computeF()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeF();
  if (ret==NOX::Abstract::Group::Ok)
  {
    mfsi_.SetupSystemMatrix(*mfsi_.SystemMatrix());
    isValidJacobian = true;
  }
  return ret;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType FSI::NOXGroup::computeJacobian()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeJacobian();
  if (ret==NOX::Abstract::Group::Ok)
  {
    mfsi_.SetupRHS(RHSVector.getEpetraVector());
    isValidRHS = true;
  }
  return ret;
}


#endif
