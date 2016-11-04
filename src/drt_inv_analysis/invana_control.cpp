/*----------------------------------------------------------------------------*/
/*!
\file invana_control.cpp

\brief Control object to handle solution of the inverse analysis

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/

/*----------------------------------------------------------------------------*/
/* headers */
#include "invana_control.H"

// TEUCHOS
#include "Teuchos_ParameterList.hpp"

// INVANA
#include "invana_base.H"
#include "invana_factory.H"
#include "invana_resulttest.H"
#include "invana_writer.H"

// INTERNAL OPTIMIZER
#include "optimizer_factory.H"
#include "optimizer_base.H"

// BACI
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_timestepping/timintmstep.H"


/*----------------------------------------------------------------------------*/
INVANA::InvanaControl::InvanaControl() :
invprob_(Teuchos::null),
invanaopt_(Teuchos::null),
input_(Teuchos::null),
x_(Teuchos::null),
f_(Teuchos::null),
val_(0.0)
{
  return;
}

/*----------------------------------------------------------------------------*/
void INVANA::InvanaControl::Init(const Teuchos::ParameterList& invp)
{
  // create an instance of an optimization problem
  INVANA::InvanaFactory invfac;
  invprob_ = invfac.Create(invp);

  // optimization algorithm
  INVANA::OptimizerFactory optimizerfac;
  invanaopt_ = optimizerfac.Create(invp);

  invanaopt_->Init(invprob_);
  invanaopt_->Setup();

  return;
}

int INVANA::InvanaControl::Solve(int restart)
{
  invanasolve(restart);

  return 0;
}

/*----------------------------------------------------------------------*/
void INVANA::InvanaControl::invanasolve(int restart)
{

  // solve
  if (restart) InvanaOpti()->ReadRestart(restart);
  InvanaOpti()->Integrate();

  // store solution
  f_=Teuchos::rcp(new Epetra_MultiVector(InvanaOpti()->GetGradientView()));
  val_=InvanaOpti()->GetObjFunctValView();
  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> INVANA::InvanaControl::CreateFieldTest()
{
  return Teuchos::rcp(new INVANA::InvanaResultTest(*this));
}
