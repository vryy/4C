/*----------------------------------------------------------------------*/
/*!
\file strtimint_noxgroup.cpp
\brief Declarations of NOX group for non-linear solution techniques
       used within implicit structural time integration

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_noxgroup.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::STR::Group::Group
(
  ::STR::TimIntImpl& sti,
  Teuchos::ParameterList& printParams,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
  const NOX::Epetra::Vector& x,
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys
)
: NOX::Epetra::Group(printParams, i, x, linSys),
  sti_(sti)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::STR::Group::computeF()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeF();
  if (ret == NOX::Abstract::Group::Ok)
  {
    if (not isValidJacobian)
    {
      sharedLinearSystem.getObject(this);
      isValidJacobian = true;
    }
  }
  return ret;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::STR::Group::computeJacobian()
{
  NOX::Abstract::Group::ReturnType ret
    = NOX::Epetra::Group::computeJacobian();
  if (ret == NOX::Abstract::Group::Ok)
  {
    if (not isValidRHS)
    {
      isValidRHS = true;
    }
  }
  return ret;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::STR::Group::computeNewton
(
  Teuchos::ParameterList& p
)
{
  NOX::Abstract::Group::ReturnType status
    = NOX::Epetra::Group::computeNewton(p);
  return status;
}

/*----------------------------------------------------------------------*/
