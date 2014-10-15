/*----------------------------------------------------------------------*/
/*!
 * \file invana_factory.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "invana_factory.H"
#include "invana_base.H"
#include "invana_auglagr.H"
#include "matpar_manager.H"
#include "matpar_manager_uniform.H"
#include "matpar_manager_elementwise.H"
#include "objective_funct.H"
#include "objective_funct_disp.H"
#include "objective_funct_surfcurr.H"
#include "regularization_base.H"
#include "regularization_tikhonov.H"
#include "regularization_totalvariation.H"
#include "optimizer_factory.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*/
/* standard constructor                                      keh 08/14  */
/*----------------------------------------------------------------------*/
STR::INVANA::InvanaFactory::InvanaFactory()
{;}

Teuchos::RCP<STR::INVANA::InvanaBase> STR::INVANA::InvanaFactory::Create(Teuchos::RCP<DRT::Discretization> discret,const Teuchos::ParameterList& invp)
{

  Teuchos::RCP<InvanaBase> optprob=Teuchos::null;

  // similarity measure
  Teuchos::RCP<STR::INVANA::ObjectiveFunct> objfunct=Teuchos::null;
  switch (DRT::INPUT::IntegralValue<INPAR::STR::StatInvObjFunctType>(invp,"OBJECTIVEFUNCT"))
  {
    case INPAR::STR::stat_inv_obj_disp:
    {
      objfunct = Teuchos::rcp(new STR::INVANA::ObjectiveFunctDisp(discret));
    }
    break;
    case INPAR::STR::stat_inv_obj_surfcurr:
    {
      objfunct = Teuchos::rcp(new STR::INVANA::ObjectiveFunctSurfCurrRepresentation(discret));
    }
    break;
    case INPAR::STR::stat_inv_obj_none:
    {
      dserror("choose some type of objective function");
    }
    break;
  }

  // parameterization
  Teuchos::RCP<STR::INVANA::MatParManager> matman = Teuchos::null;
  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvMatParametrization>(invp,"PARAMETRIZATION"))
  {
    case INPAR::STR::stat_inv_mp_elementwise:
    {
      matman = Teuchos::rcp(new STR::INVANA::MatParManagerPerElement(discret));
    }
    break;
    case INPAR::STR::stat_inv_mp_uniform:
    {
      matman = Teuchos::rcp(new STR::INVANA::MatParManagerUniform(discret));
    }
    break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
      break;
  }

  // regularization!
  Teuchos::RCP<STR::INVANA::RegularizationBase> regman = Teuchos::null;
  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvRegularization>(invp,"REGULARIZATION"))
  {
    case INPAR::STR::stat_inv_reg_none:
      break;
    case INPAR::STR::stat_inv_reg_tikhonov:
    {
      regman = Teuchos::rcp(new STR::INVANA::RegularizationTikhonov(invp));
    }
    break;
    case INPAR::STR::stat_inv_reg_totalvariation:
    {
      regman = Teuchos::rcp(new STR::INVANA::RegularizationTotalVariation(invp));
    }
    break;
  }
  if (regman!=Teuchos::null)
  {
    regman->Init(discret,matman->GetConnectivityData());
    regman->Setup();
  }

  // optimization algorithm
  STR::INVANA::OptimizerFactory optimizerfac;
  Teuchos::RCP<STR::INVANA::OptimizerBase> opti = optimizerfac.Create(invp);

  // in here comes the switch for various optimization problems. So far only augmented lagrange functional
  optprob = Teuchos::rcp(new STR::INVANA::InvanaAugLagr());

  // here some checks of valid combinations of objective functions, parametrizations and
  // regularizations should come

  // initialize optimization problem
  optprob->Init(discret,objfunct,matman,regman,opti,optprob);
  optprob->Setup();

  return optprob;

}
