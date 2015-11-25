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
#include "regularization_tvdtikh.H"
#include "optimizer_factory.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*/
/* standard constructor                                      keh 08/14  */
/*----------------------------------------------------------------------*/
INVANA::InvanaFactory::InvanaFactory()
{;}

Teuchos::RCP<INVANA::InvanaBase> INVANA::InvanaFactory::Create(Teuchos::RCP<DRT::Discretization> discret,const Teuchos::ParameterList& invp)
{

  Teuchos::RCP<InvanaBase> optprob=Teuchos::null;

  // similarity measure
  Teuchos::RCP<INVANA::ObjectiveFunct> objfunct=Teuchos::null;
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvObjFunctType>(invp,"OBJECTIVEFUNCT"))
  {
    case INPAR::INVANA::stat_inv_obj_disp:
    {
      objfunct = Teuchos::rcp(new INVANA::ObjectiveFunctDisp(discret));
    }
    break;
    case INPAR::INVANA::stat_inv_obj_surfcurr:
    {
#if defined(HAVE_Kokkos)
      objfunct = Teuchos::rcp(new INVANA::SurfCurrentGroup(discret));
#else
      dserror("You need Kokkos for Surface Current based objective functions");
#endif
    }
    break;
    case INPAR::INVANA::stat_inv_obj_none:
    {
      dserror("choose some type of objective function");
    }
    break;
  }

  // parameterization
  Teuchos::RCP<INVANA::MatParManager> matman = Teuchos::null;
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(invp,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
    {
      matman = Teuchos::rcp(new INVANA::MatParManagerPerElement(discret));
    }
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
    {
      matman = Teuchos::rcp(new INVANA::MatParManagerUniform(discret));
    }
    break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
      break;
  }
  matman->Init(invp);
  matman->Setup();

  // regularization!
  Teuchos::RCP<INVANA::RegularizationBase> regman = Teuchos::null;
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp,"REGULARIZATION"))
  {
    case INPAR::INVANA::stat_inv_reg_none:
      break;
    case INPAR::INVANA::stat_inv_reg_tikhonov:
    {
      regman = Teuchos::rcp(new INVANA::RegularizationTikhonov());
    }
    break;
    case INPAR::INVANA::stat_inv_reg_totalvariation:
    {
      regman = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
    }
    break;
    case INPAR::INVANA::stat_inv_reg_tvdtikh:
    {
      regman = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
    }
    break;
  }
  if (regman!=Teuchos::null)
  {
    regman->Init(discret,matman->GetConnectivityData());
    regman->Setup(invp);
  }

  // optimization algorithm
  INVANA::OptimizerFactory optimizerfac;
  Teuchos::RCP<INVANA::OptimizerBase> opti = optimizerfac.Create(invp);

  // in here comes the switch for various optimization problems. So far only augmented lagrange functional
  optprob = Teuchos::rcp(new INVANA::InvanaAugLagr());

  // here some checks of valid combinations of objective functions, parametrizations and
  // regularizations should come

  // initialize optimization problem
  optprob->Init(discret,objfunct,matman,regman,opti,optprob);
  optprob->Setup();

  return optprob;

}
