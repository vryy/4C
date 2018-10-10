/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_factory.cpp
 * \brief Factory for optimzation algorithm within the inverse analysis

\level 3

\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
*/
/*----------------------------------------------------------------------*/
#include "optimizer_factory.H"
#include "optimizer_lbfgs.H"
#include "optimizer_smc.H"
#include "optimizer_bruteforce.H"
#include "optimizer_mh.H"
#include "optimizer_smc_predict.H"
#include "invana_base.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_statinvanalysis.H"

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
#if __cplusplus >= 201103L
      opti = Teuchos::rcp(new INVANA::OptimizerSMC(invp));
#else
      dserror("Compile with >=C++11.");
#endif
    }
    break;
    case INPAR::INVANA::stat_inv_mh:
    {
#if __cplusplus >= 201103L
      opti = Teuchos::rcp(new INVANA::OptimizerMH(invp));
#else
      dserror("Compile with >=C++11.");
#endif
    }
    break;
    case INPAR::INVANA::stat_inv_prediction:
    {
#if __cplusplus >= 201103L
      opti = Teuchos::rcp(new INVANA::PredictionSMC(invp));
#else
      dserror("Compile with >=C++11.");
#endif
    }
    break;
    default:
      dserror("Unknown type of statistical inverse analysis");
      break;
  }

  return opti;
}
