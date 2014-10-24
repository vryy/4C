/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_factory.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "optimizer_factory.H"
#include "optimizer_lbfgs.H"
#include "optimizer_graddesc.H"
#include "optimizer_mc.H"
#include "invana_base.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_statinvanalysis.H"

#include "Teuchos_ParameterList.hpp"


/*----------------------------------------------------------------------*/
/* standard constructor                                      keh 08/14  */
/*----------------------------------------------------------------------*/
INVANA::OptimizerFactory::OptimizerFactory()
{;}

Teuchos::RCP<INVANA::OptimizerBase> INVANA::OptimizerFactory::Create(const Teuchos::ParameterList& invp)
{

  Teuchos::RCP<OptimizerBase> opti=Teuchos::null;

  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(invp,"STAT_INV_ANALYSIS"))
  {

    case INPAR::INVANA::stat_inv_graddesc:
    {
      opti = Teuchos::rcp(new INVANA::OptimizerGradDesc(invp));
    }
    break;

    case INPAR::INVANA::stat_inv_lbfgs:
    {
      opti = Teuchos::rcp(new INVANA::OptimizerLBFGS(invp));
    }
    break;
    case INPAR::INVANA::stat_inv_mc:
    {
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
{
        opti = Teuchos::rcp(new INVANA::OptimizerMC(invp));
}
#else
        dserror("Install new Boost version");
#endif
    }
    break;
    default:
      dserror("Unknown type of statistical inverse analysis");
    break;
  }

  return opti;

}
