/*----------------------------------------------------------------------*/
/*!
\file invana_factory.cpp

\brief Factory for the inverse analysis
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "invana_factory.H"

// Invana
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

// Input
#include "../drt_inpar/inpar_statinvanalysis.H"

//Baci
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

// Teuchos
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
INVANA::InvanaFactory::InvanaFactory()
{;}

Teuchos::RCP<INVANA::InvanaBase> INVANA::InvanaFactory::Create(
    const Teuchos::ParameterList& invp)
{
  // get the discretizatio to work with
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // check for a proper status of the discretization
  if (!actdis->Filled() or !actdis->HaveDofs())
    actdis->FillComplete();

  // initialize the basic inverse problem
  Teuchos::RCP<InvanaBase> optprob=Teuchos::null;

  // similarity measure
  Teuchos::RCP<INVANA::ObjectiveFunct> objfunct=Teuchos::null;
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvObjFunctType>(invp,"OBJECTIVEFUNCT"))
  {
    case INPAR::INVANA::stat_inv_obj_disp:
    {
      objfunct = Teuchos::rcp(new INVANA::ObjectiveFunctDisp(actdis));
    }
    break;
    case INPAR::INVANA::stat_inv_obj_surfcurr:
    {
#if defined( HAVE_Kokkos )
      objfunct = Teuchos::rcp(new INVANA::SurfCurrentGroup(actdis));
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
      matman = Teuchos::rcp(new INVANA::MatParManagerPerElement(actdis));
    }
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
    {
      matman = Teuchos::rcp(new INVANA::MatParManagerUniform(actdis));
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
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(
      invp,"REGULARIZATION"))
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
  }
  if (regman!=Teuchos::null)
  {
    regman->Init(matman->GetConnectivityData());
    regman->Setup(invp);
  }

  // In here comes the switch for various optimization problems.
  // So far only augmented lagrange functional
  optprob = Teuchos::rcp(new INVANA::InvanaAugLagr());

  // Some checks of valid combinations of objective functions,
  // parametrizations and regularizations should come here.

  // initialize optimization problem
  optprob->Init(actdis,objfunct,matman,regman);
  optprob->Setup();

  return optprob;

}
