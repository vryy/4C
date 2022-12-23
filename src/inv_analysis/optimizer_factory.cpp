/*----------------------------------------------------------------------*/
/*! \file
\brief Factory for optimzation algorithm within the inverse analysis

\level 3


*/
/*----------------------------------------------------------------------*/
#include "optimizer_factory.H"
#include "optimizer_lbfgs.H"
#include "optimizer_smc.H"
#include "optimizer_bruteforce.H"
#include "optimizer_mh.H"
#include "optimizer_smc_predict.H"
#include "invana_base.H"
#include "dserror.H"
#include "inpar_statinvanalysis.H"

#include "Teuchos_ParameterList.hpp"

/*----------------------------------------------------------------------*/
INVANA::OptimizerFactory::OptimizerFactory() { ; }

Teuchos::RCP<INVANA::OptimizerBase> INVANA::OptimizerFactory::Create(
    const Teuchos::ParameterList& invp)
{
  Teuchos::RCP<OptimizerBase> opti = Teuchos::null;

  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(invp, "STAT_INV_ANALYSIS"))
  {
    case INPAR::INVANA::stat_inv_lbfgs:
    {
      opti = Teuchos::rcp(new INVANA::OptimizerLBFGS(invp));
    }
    break;
    case INPAR::INVANA::stat_inv_bruteforce:
    {
      opti = Teuchos::rcp(new INVANA::OptimizerBruteForce(invp));
    }
    break;
    case INPAR::INVANA::stat_inv_smc:
    {
      opti = Teuchos::rcp(new INVANA::OptimizerSMC(invp));
    }
    break;
    case INPAR::INVANA::stat_inv_mh:
    {
      opti = Teuchos::rcp(new INVANA::OptimizerMH(invp));
    }
    break;
    case INPAR::INVANA::stat_inv_prediction:
    {
      opti = Teuchos::rcp(new INVANA::PredictionSMC(invp));
    }
    break;
    default:
      dserror("Unknown type of statistical inverse analysis");
      break;
  }

  return opti;
}
