/*----------------------------------------------------------------------------*/
/*! \file
\brief Control object to handle solution of the inverse analysis

\level 3

*/

/*----------------------------------------------------------------------------*/
/* headers */
#include "baci_inv_analysis_control.H"

// TEUCHOS
#include <Teuchos_ParameterList.hpp>

// INVANA
#include "baci_inv_analysis_base.H"
#include "baci_inv_analysis_factory.H"
#include "baci_inv_analysis_resulttest.H"
#include "baci_inv_analysis_writer.H"

// INTERNAL OPTIMIZER
#include "baci_inv_analysis_optimizer_factory.H"
#include "baci_inv_analysis_optimizer_base.H"

// BACI
#include "baci_utils_exceptions.H"
#include "baci_inpar_statinvanalysis.H"
#include "baci_lib_globalproblem.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_timestepping_mstep.H"


/*----------------------------------------------------------------------------*/
INVANA::InvanaControl::InvanaControl()
    : invprob_(Teuchos::null),
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
  f_ = Teuchos::rcp(new Epetra_MultiVector(InvanaOpti()->GetGradientView()));
  val_ = InvanaOpti()->GetObjFunctValView();
  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> INVANA::InvanaControl::CreateFieldTest()
{
  return Teuchos::rcp(new INVANA::InvanaResultTest(*this));
}
