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
STR::INVANA::OptimizerFactory::OptimizerFactory()
{;}

Teuchos::RCP<STR::INVANA::OptimizerBase> STR::INVANA::OptimizerFactory::Create(const Teuchos::ParameterList& invp)
{

  Teuchos::RCP<OptimizerBase> opti=Teuchos::null;

  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvAnalysisType>(invp,"STAT_INV_ANALYSIS"))
  {

    case INPAR::STR::stat_inv_graddesc:
    {
      opti = Teuchos::rcp(new STR::INVANA::OptimizerGradDesc(invp));
    }
    break;

    case INPAR::STR::stat_inv_lbfgs:
    {
      opti = Teuchos::rcp(new STR::INVANA::OptimizerLBFGS(invp));
    }
    break;
    case INPAR::STR::stat_inv_mc:
    {
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
{
        opti = Teuchos::rcp(new STR::INVANA::OptimizerMC(invp));
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
