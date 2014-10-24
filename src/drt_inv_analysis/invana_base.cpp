/*----------------------------------------------------------------------*/
/*!
 * \file invana_base.cpp

<pre>
Maintainer: Sebastian Kehl
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
/* standard constructor                                      keh 09/14  */
/*----------------------------------------------------------------------*/
INVANA::InvanaBase::InvanaBase():
discret_(Teuchos::null),
objfunct_(Teuchos::null),
matman_(Teuchos::null),
regman_(Teuchos::null),
optimizer_(Teuchos::null),
isinit_(false)
{;}

void INVANA::InvanaBase::Init(Teuchos::RCP<DRT::Discretization> discret,
                                   Teuchos::RCP<INVANA::ObjectiveFunct> objfunct,
                                   Teuchos::RCP<INVANA::MatParManager> matman,
                                   Teuchos::RCP<INVANA::RegularizationBase> regman,
                                   Teuchos::RCP<INVANA::OptimizerBase> optimizer,
                                   Teuchos::RCP<INVANA::InvanaBase> optprob)
{
  discret_=discret;
  objfunct_=objfunct;
  matman_=matman;
  regman_=regman;
  optimizer_=optimizer;

  optimizer_->Init(VectorRowLayout(),VectorColLayout(),NumVectors(),optprob.create_weak());
  optimizer_->Setup();

  isinit_=true;
}

const Epetra_MultiVector& INVANA::InvanaBase::InitialGuess()
{
  return matman_->InitialParams();
}

void INVANA::InvanaBase::Solve(int restart)
{
  if (!isinit_) dserror("InvanaBase is not inititialzed. Call Init() first");

  if (restart) optimizer_->ReadRestart(restart);
  optimizer_->Integrate();

  return;
}

const Epetra_Comm& INVANA::InvanaBase::Comm()
{
  return discret_->Comm();
}

Teuchos::RCP<Epetra_Map> INVANA::InvanaBase::VectorRowLayout()
{
  return matman_->ParamLayoutMapUnique();
}

Teuchos::RCP<Epetra_Map> INVANA::InvanaBase::VectorColLayout()
{
  return matman_->ParamLayoutMap();
}

double INVANA::InvanaBase::NumVectors()
{
  return matman_->NumVectors();
}

Teuchos::RCP<DRT::ResultTest> INVANA::InvanaBase::CreateFieldTest()
{
  return Teuchos::rcp(new InvanaResultTest(*this));
}
