/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired integrator object.

\level 2


*/
/*---------------------------------------------------------------------*/

#include "4C_contact_integrator_factory.hpp"

#include "4C_contact_aug_integrator.hpp"
#include "4C_contact_ehl_integrator.hpp"
#include "4C_contact_nitsche_integrator.hpp"
#include "4C_contact_nitsche_integrator_fpi.hpp"
#include "4C_contact_nitsche_integrator_fsi.hpp"
#include "4C_contact_nitsche_integrator_poro.hpp"
#include "4C_contact_nitsche_integrator_ssi.hpp"
#include "4C_contact_nitsche_integrator_ssi_elch.hpp"
#include "4C_contact_nitsche_integrator_tsi.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Integrator> CONTACT::INTEGRATOR::Factory::BuildIntegrator(
    const Inpar::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const Core::FE::CellType& slave_type, const Epetra_Comm& comm) const
{
  Teuchos::RCP<CONTACT::Integrator> integrator = Teuchos::null;
  switch (sol_type)
  {
    case Inpar::CONTACT::solution_augmented:
    case Inpar::CONTACT::solution_std_lagrange:
    case Inpar::CONTACT::solution_steepest_ascent:
    case Inpar::CONTACT::solution_steepest_ascent_sp:
    case Inpar::CONTACT::solution_combo:
    {
      integrator = Teuchos::rcp<CONTACT::Integrator>(
          new CONTACT::Aug::IntegrationWrapper(mortar_params, slave_type, comm));
      break;
    }
    case Inpar::CONTACT::solution_nitsche:
    {
      if (mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::tsi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheTsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::ssi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheSsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::ssi_elch)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheSsiElch(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::poroelast ||
               mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::poroscatra)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitschePoro(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::fsi)
      {
        integrator =
            Teuchos::rcp(new CONTACT::IntegratorNitscheFsi(mortar_params, slave_type, comm));
      }
      else if (mortar_params.get<int>("PROBTYPE") == Inpar::CONTACT::fpi)
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
    case Inpar::CONTACT::solution_penalty:
    case Inpar::CONTACT::solution_multiscale:
    {
      if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(mortar_params, "ALGORITHM") ==
          Inpar::Mortar::algorithm_gpts)
        integrator = Teuchos::rcp(new CONTACT::IntegratorNitsche(mortar_params, slave_type, comm));
      else
        integrator = Teuchos::rcp(new CONTACT::Integrator(mortar_params, slave_type, comm));
      break;
    }
    case Inpar::CONTACT::solution_lagmult:
    case Inpar::CONTACT::solution_uzawa:
    {
      integrator = Teuchos::rcp(new CONTACT::Integrator(mortar_params, slave_type, comm));
      break;
    }
    case Inpar::CONTACT::solution_ehl:
    {
      integrator = Teuchos::rcp(new CONTACT::IntegratorEhl(mortar_params, slave_type, comm));

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported solving strategy! (stype = %s | %d)",
          Inpar::CONTACT::SolvingStrategy2String(sol_type).c_str(), sol_type);
      exit(EXIT_FAILURE);
    }
  }  // end switch

  return integrator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Integrator> CONTACT::INTEGRATOR::BuildIntegrator(
    const Inpar::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const Core::FE::CellType& slave_type, const Epetra_Comm& comm)
{
  Factory factory;
  return factory.BuildIntegrator(sol_type, mortar_params, slave_type, comm);
}

FOUR_C_NAMESPACE_CLOSE
