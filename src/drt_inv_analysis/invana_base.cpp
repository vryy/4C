/*----------------------------------------------------------------------*/
/*!
\file invana_base.cpp

\brief Base class for the inverse analysis


<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "invana_base.H"
#include "matpar_manager.H"
#include "objective_funct.H"
#include "regularization_base.H"
#include "optimizer_base.H"
#include "invana_resulttest.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*/
INVANA::InvanaBase::InvanaBase():
discret_(Teuchos::null),
objfunct_(Teuchos::null),
matman_(Teuchos::null),
regman_(Teuchos::null),
isinit_(false)
{;}

/*----------------------------------------------------------------------*/
void INVANA::InvanaBase::Init(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<INVANA::ObjectiveFunct> objfunct,
    Teuchos::RCP<INVANA::MatParManager> matman,
    Teuchos::RCP<INVANA::RegularizationBase> regman)
{
  discret_=discret;
  objfunct_=objfunct;
  matman_=matman;
  regman_=regman;

  isinit_=true;
}

const Epetra_MultiVector& INVANA::InvanaBase::InitialGuess()
{
  return matman_->InitialParams();
}

const Epetra_Comm& INVANA::InvanaBase::Comm()
{
  return discret_->Comm();
}

Teuchos::RCP<Epetra_Map> INVANA::InvanaBase::VectorRowLayout()
{
  return matman_->ParamLayoutMapUnique();
}
