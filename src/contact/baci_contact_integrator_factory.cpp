/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired integrator object.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "baci_contact_integrator_factory.hpp"

#include "baci_contact_aug_integrator.hpp"
#include "baci_contact_ehl_integrator.hpp"
#include "baci_contact_nitsche_integrator.hpp"
#include "baci_contact_nitsche_integrator_fpi.hpp"
#include "baci_contact_nitsche_integrator_fsi.hpp"
#include "baci_contact_nitsche_integrator_poro.hpp"
#include "baci_contact_nitsche_integrator_ssi.hpp"
#include "baci_contact_nitsche_integrator_ssi_elch.hpp"
#include "baci_contact_nitsche_integrator_tsi.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Integrator> CONTACT::INTEGRATOR::Factory::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const CORE::FE::CellType& slave_type, const Epetra_Comm& comm) const
{
  Teuchos::RCP<CONTACT::Integrator> integrator = Teuchos::null;
  switch (sol_type)
  {
    case INPAR::CONTACT::solution_augmented:
    case INPAR::CONTACT::solution_std_lagrange:
    case INPAR::CONTACT::solution_steepest_ascent:
    case INPAR::CONTACT::solution_steepest_ascent_sp:
    case INPAR::CONTACT::solution_combo:
    {
      integrator = Teuchos::rcp<CONTACT::Integrator>(
          new CONTACT::AUG::IntegrationWrapper(mortar_params, slave_type, comm));
      break;
    }
    case INPAR::CONTACT::solution_nitsche:
    {
      if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::tsi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheTsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::ssi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheSsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::ssi_elch)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheSsiElch(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::poroelast ||
               mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::poroscatra)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitschePoro(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::fsi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheFsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == INPAR::CONTACT::fpi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheFpi(mortar_params, slave_type, comm));
      }
      else
      {
        integrator = Teuchos::rcp(new CONTACT::IntegratorNitsche(mortar_params, slave_type, comm));
      }
      break;
    }
    case INPAR::CONTACT::solution_penalty:
    case INPAR::CONTACT::solution_multiscale:
    {
      if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(mortar_params, "ALGORITHM") ==
          INPAR::MORTAR::algorithm_gpts)
        integrator = Teuchos::rcp(new CONTACT::IntegratorNitsche(mortar_params, slave_type, comm));
      else
        integrator = Teuchos::rcp(new CONTACT::Integrator(mortar_params, slave_type, comm));
      break;
    }
    case INPAR::CONTACT::solution_lagmult:
    case INPAR::CONTACT::solution_uzawa:
    {
      integrator = Teuchos::rcp(new CONTACT::Integrator(mortar_params, slave_type, comm));
      break;
    }
    case INPAR::CONTACT::solution_ehl:
    {
      integrator = Teuchos::rcp(new CONTACT::IntegratorEhl(mortar_params, slave_type, comm));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported solving strategy! (stype = %s | %d)",
          INPAR::CONTACT::SolvingStrategy2String(sol_type).c_str(), sol_type);
      exit(EXIT_FAILURE);
    }
  }  // end switch

  return integrator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Integrator> CONTACT::INTEGRATOR::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const CORE::FE::CellType& slave_type, const Epetra_Comm& comm)
{
  Factory factory;
  return factory.BuildIntegrator(sol_type, mortar_params, slave_type, comm);
}

FOUR_C_NAMESPACE_CLOSE
