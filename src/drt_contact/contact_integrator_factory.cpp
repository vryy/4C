/*---------------------------------------------------------------------*/
/*!
\file contact_integrator_factory.cpp

\brief Factory to create the desired integrator object.

\level 2

\maintainer Michael Hiermeier

\date 04/2016

*/
/*---------------------------------------------------------------------*/

#include "contact_integrator_factory.H"

// supported contact integrators
#include "contact_integrator.H"
#include "../drt_contact_aug/contact_augmented_integrator.H"
#include "contact_nitsche_integrator.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::INTEGRATOR::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoIntegrator> CONTACT::INTEGRATOR::Factory::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& stype,
    Teuchos::ParameterList& p_mortar,
    const DRT::Element::DiscretizationType& eletype,
    const Epetra_Comm& comm) const
{
  Teuchos::RCP<CONTACT::CoIntegrator> integrator = Teuchos::null;
  switch (stype)
  {
    case INPAR::CONTACT::solution_augmented:
    {
      integrator = Teuchos::rcp(new CONTACT::AugmentedIntegrator(
          p_mortar,eletype,comm));
      break;
    }
    case INPAR::CONTACT::solution_nitsche:
    {
      integrator = Teuchos::rcp(new CONTACT::CoIntegratorNitsche(
          p_mortar,eletype,comm));
      break;
    }
    case INPAR::CONTACT::solution_penalty:
    {
      if(DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(p_mortar, "ALGORITHM")
          == INPAR::MORTAR::algorithm_gpts)
        integrator = Teuchos::rcp(new CONTACT::CoIntegratorNitsche(
            p_mortar,eletype,comm));
      else
        integrator = Teuchos::rcp(new CONTACT::CoIntegrator(
                  p_mortar,eletype,comm));
      break;
    }
    default:
    {
      integrator = Teuchos::rcp(new CONTACT::CoIntegrator(
          p_mortar,eletype,comm));
      break;
    }
  } // end switch

  return integrator;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoIntegrator> CONTACT::INTEGRATOR::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& stype,
    Teuchos::ParameterList& p_mortar,
    const DRT::Element::DiscretizationType& eletype,
    const Epetra_Comm& comm)
{
  Factory factory;
  return factory.BuildIntegrator(stype,p_mortar,eletype,comm);
}
