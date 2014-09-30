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
/* standard constructor                                      keh 10/13  */
/*----------------------------------------------------------------------*/
STR::INVANA::InvanaBase::InvanaBase():
discret_(Teuchos::null),
objfunct_(Teuchos::null),
matman_(Teuchos::null),
regman_(Teuchos::null),
optimizer_(Teuchos::null),
isinit_(false)
{;}

void STR::INVANA::InvanaBase::Init(Teuchos::RCP<DRT::Discretization> discret,
                                   Teuchos::RCP<STR::INVANA::ObjectiveFunct> objfunct,
                                   Teuchos::RCP<STR::INVANA::MatParManager> matman,
                                   Teuchos::RCP<STR::INVANA::RegularizationBase> regman,
                                   Teuchos::RCP<STR::INVANA::OptimizerBase> optimizer,
                                   Teuchos::RCP<STR::INVANA::InvanaBase> optprob)
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

const Epetra_MultiVector& STR::INVANA::InvanaBase::InitialGuess()
{
  return matman_->InitialParams();
}

void STR::INVANA::InvanaBase::Solve(int restart)
{
  if (!isinit_) dserror("InvanaBase is not inititialzed. Call Init() first");

  if (restart) optimizer_->ReadRestart(restart);
  optimizer_->Integrate();

  return;
}

const Epetra_Comm& STR::INVANA::InvanaBase::Comm()
{
  return discret_->Comm();
}

Teuchos::RCP<Epetra_Map> STR::INVANA::InvanaBase::VectorRowLayout()
{
  return matman_->ParamLayoutMapUnique();
}

Teuchos::RCP<Epetra_Map> STR::INVANA::InvanaBase::VectorColLayout()
{
  return matman_->ParamLayoutMap();
}

double STR::INVANA::InvanaBase::NumVectors()
{
  return matman_->NumVectors();
}

/*----------------------------------------------------------------------*/
/* Creates the field test                                               */
Teuchos::RCP<DRT::ResultTest> STR::INVANA::InvanaBase::CreateFieldTest()
{
  return Teuchos::rcp(new InvanaResultTest(*this));
}
