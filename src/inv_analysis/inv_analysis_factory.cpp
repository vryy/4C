/*----------------------------------------------------------------------*/
/*! \file
\brief Factory for the inverse analysis

\level 3

*/
/*----------------------------------------------------------------------*/
#include "inv_analysis_factory.H"

// Invana
#include "inv_analysis_base.H"
#include "inv_analysis_auglagr.H"
#include "inv_analysis_matpar_manager.H"
#include "inv_analysis_matpar_manager_uniform.H"
#include "inv_analysis_matpar_manager_elementwise.H"
#include "inv_analysis_matpar_manager_patchwise.H"
#include "inv_analysis_matpar_manager_tvsvd.H"
#include "inv_analysis_objective_funct.H"
#include "inv_analysis_objective_funct_disp.H"
#include "inv_analysis_objective_funct_surfcurr.H"
#include "inv_analysis_regularization_base.H"
#include "inv_analysis_regularization_tikhonov.H"
#include "inv_analysis_regularization_totalvariation.H"
#include "inv_analysis_initial_guess.H"

// Input
#include "inpar_statinvanalysis.H"

// Baci
#include "utils_exceptions.H"
#include "lib_discret.H"
#include "lib_globalproblem.H"
#include "comm_utils.H"

// Teuchos
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
INVANA::InvanaFactory::InvanaFactory() { ; }

Teuchos::RCP<INVANA::InvanaBase> INVANA::InvanaFactory::Create(const Teuchos::ParameterList& invp)
{
  // get the discretizatio to work with
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // check for a proper status of the discretization
  if (!actdis->Filled() or !actdis->HaveDofs()) actdis->FillComplete();

  // initialize the basic inverse problem
  Teuchos::RCP<InvanaBase> optprob = Teuchos::null;

  // similarity measure
  Teuchos::RCP<INVANA::ObjectiveFunct> objfunct = Teuchos::null;
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvObjFunctType>(invp, "OBJECTIVEFUNCT"))
  {
    case INPAR::INVANA::stat_inv_obj_disp:
    {
      objfunct = Teuchos::rcp(new INVANA::ObjectiveFunctDisp(actdis));
    }
    break;
    case INPAR::INVANA::stat_inv_obj_surfcurr:
    {
      objfunct = Teuchos::rcp(new INVANA::SurfCurrentGroup(actdis));
    }
    break;
    case INPAR::INVANA::stat_inv_obj_none:
    {
      dserror("choose some type of objective function");
    }
    break;
  }

  // parametrization
  Teuchos::RCP<INVANA::MatParManager> matman = Teuchos::null;
  switch (
      DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(invp, "PARAMETRIZATION"))
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
    case INPAR::INVANA::stat_inv_mp_patchwise:
    {
      matman = Teuchos::rcp(new INVANA::MatParManagerPerPatch(actdis));
    }
    break;
    case INPAR::INVANA::stat_inv_mp_tvsvd:
    {
      matman = Teuchos::rcp(new INVANA::MatParManagerTVSVD(actdis));
    }
    break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
      break;
  }
  matman->Init(invp, objfunct);
  matman->Setup();

  // initial guess
  Teuchos::RCP<InitialGuess> guess = Teuchos::rcp(new InitialGuess(invp));
  guess->Compute(actdis, matman, objfunct);

  // regularization
  Teuchos::RCP<INVANA::RegularizationBase> regman = Teuchos::null;
  switch (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp, "REGULARIZATION"))
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
  if (regman != Teuchos::null)
  {
    regman->Init(invp, matman->GetConnectivityData(), guess);
    regman->Setup(invp);
  }

  // In here comes the switch for various optimization problems.
  // So far only augmented lagrange functional
  optprob = Teuchos::rcp(new INVANA::InvanaAugLagr());

  // Some checks of valid combinations of objective functions,
  // parametrizations and regularizations should come here.

  // initialize optimization problem
  optprob->Init(actdis, objfunct, matman, regman, guess);
  optprob->Setup();

  return optprob;
}
