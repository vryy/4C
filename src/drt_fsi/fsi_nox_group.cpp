/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of NOX::Group for FSI

\level 1

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------*/

#include "fsi_nox_group.H"
#include "fsi_monolithicinterface.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::Group::Group(::FSI::MonolithicInterface& mfsi, Teuchos::ParameterList& printParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys)
    : NOX::Epetra::Group(printParams, i, x, linSys), mfsi_(mfsi)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::Group::CaptureSystemState()
{
  // we know we already have the first linear system calculated

  mfsi_.SetupRHS(RHSVector.getEpetraVector(), true);
  mfsi_.SetupSystemMatrix();

  sharedLinearSystem.getObject(this);
  isValidJacobian = true;
  isValidRHS = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::FSI::Group::computeF()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeF();
  if (ret == NOX::Abstract::Group::Ok)
  {
    if (not isValidJacobian)
    {
      mfsi_.SetupSystemMatrix();
      sharedLinearSystem.getObject(this);
      isValidJacobian = true;
    }
  }
  return ret;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::FSI::Group::computeJacobian()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeJacobian();
  if (ret == NOX::Abstract::Group::Ok)
  {
    if (not isValidRHS)
    {
      mfsi_.SetupRHS(RHSVector.getEpetraVector());
      isValidRHS = true;
    }
  }
  return ret;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::FSI::Group::computeNewton(Teuchos::ParameterList& p)
{
  mfsi_.ScaleSystem(RHSVector.getEpetraVector());
  NOX::Abstract::Group::ReturnType status = NOX::Epetra::Group::computeNewton(p);
  mfsi_.UnscaleSolution(NewtonVector.getEpetraVector(), RHSVector.getEpetraVector());

  // check return value of computeNewton call
  if (status == NOX::Abstract::Group::NotConverged || status == NOX::Abstract::Group::Failed)
    dserror("NOX::FSI::Group::computeNewton: linear solver not converged...");

  return status;
}
