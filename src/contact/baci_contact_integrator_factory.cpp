/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired integrator object.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_integrator_factory.H"

#include "baci_contact_aug_integrator.H"
#include "baci_contact_ehl_integrator.H"
#include "baci_contact_nitsche_integrator.H"
#include "baci_contact_nitsche_integrator_fpi.H"
#include "baci_contact_nitsche_integrator_fsi.H"
#include "baci_contact_nitsche_integrator_poro.H"
#include "baci_contact_nitsche_integrator_ssi.H"
#include "baci_contact_nitsche_integrator_ssi_elch.H"
#include "baci_contact_nitsche_integrator_tsi.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoIntegrator> CONTACT::INTEGRATOR::Factory::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const DRT::Element::DiscretizationType& slave_type, const Epetra_Comm& comm) const
{
  Teuchos::RCP<CONTACT::CoIntegrator> integrator = Teuchos::null;
  switch (sol_type)
  {
    case INPAR::CONTACT::solution_augmented:
    case INPAR::CONTACT::solution_std_lagrange:
    case INPAR::CONTACT::solution_steepest_ascent:
    case INPAR::CONTACT::solution_steepest_ascent_sp:
    case INPAR::CONTACT::solution_combo:
    {
      integrator = Teuchos::rcp<CONTACT::CoIntegrator>(
          new CONTACT::AUG::IntegrationWrapper(mortar_params, slave_type, comm));
      break;
    }
    case INPAR::CONTACT::solution_nitsche:
    {
      if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::tsi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitscheTsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::ssi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitscheSsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::ssi_elch)
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitscheSsiElch(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::poroelast ||
               mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::poroscatra)
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitschePoro(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::fsi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitscheFsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::fpi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitscheFpi(mortar_params, slave_type, comm));
      }
      else
      {
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitsche(mortar_params, slave_type, comm));
      }
      break;
    }
    case INPAR::CONTACT::solution_penalty:
    case INPAR::CONTACT::solution_multiscale:
    {
      if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar_params, "ALGORITHM") ==
          INPAR::MORTAR::algorithm_gpts)
        integrator =
            Teuchos::rcp(new CONTACT::CoIntegratorNitsche(mortar_params, slave_type, comm));
      else
        integrator = Teuchos::rcp(new CONTACT::CoIntegrator(mortar_params, slave_type, comm));
      break;
    }
    case INPAR::CONTACT::solution_lagmult:
    case INPAR::CONTACT::solution_uzawa:
    {
      integrator = Teuchos::rcp(new CONTACT::CoIntegrator(mortar_params, slave_type, comm));
      break;
    }
    case INPAR::CONTACT::solution_ehl:
    {
      integrator = Teuchos::rcp(new CONTACT::CoIntegratorEhl(mortar_params, slave_type, comm));

      break;
    }
    default:
    {
      dserror("Unsupported solving strategy! (stype = %s | %d)",
          INPAR::CONTACT::SolvingStrategy2String(sol_type).c_str(), sol_type);
      exit(EXIT_FAILURE);
    }
  }  // end switch

  return integrator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoIntegrator> CONTACT::INTEGRATOR::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const DRT::Element::DiscretizationType& slave_type, const Epetra_Comm& comm)
{
  Factory factory;
  return factory.BuildIntegrator(sol_type, mortar_params, slave_type, comm);
}
