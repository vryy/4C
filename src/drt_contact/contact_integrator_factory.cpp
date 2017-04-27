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
#include "../drt_contact_aug/contact_augmented_integrator.H"
#include "contact_nitsche_integrator.H"
#include "contact_nitsche_integrator_tsi.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::INTEGRATOR::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoIntegrator> CONTACT::INTEGRATOR::Factory::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& sol_type,
    Teuchos::ParameterList& p_mortar,
    const DRT::Element::DiscretizationType& slave_type,
    const Epetra_Comm& comm) const
{
  Teuchos::RCP<CONTACT::CoIntegrator> integrator = Teuchos::null;
  switch (sol_type)
  {
    case INPAR::CONTACT::solution_augmented:
    case INPAR::CONTACT::solution_std_lagrange:
    case INPAR::CONTACT::solution_steepest_ascent:
    case INPAR::CONTACT::solution_combo:
    {
      integrator = Teuchos::rcp<CONTACT::CoIntegrator>(
          new CONTACT::AUG::IntegrationWrapper( p_mortar,slave_type,comm) );
      break;
    }
    case INPAR::CONTACT::solution_xcontact:
    {
      integrator = Teuchos::rcp<CONTACT::CoIntegrator>( BuildXIntegrator(
          p_mortar,slave_type,comm) );
      break;
    }
    case INPAR::CONTACT::solution_nitsche:
    {
      if (p_mortar.get<int>("PROBTYPE")==INPAR::CONTACT::tsi)
        integrator = Teuchos::rcp( new CONTACT::CoIntegratorNitscheTsi(
            p_mortar,slave_type,comm) );
      else
        integrator = Teuchos::rcp( new CONTACT::CoIntegratorNitsche(
            p_mortar,slave_type,comm) );
      break;
    }
    case INPAR::CONTACT::solution_penalty:
    {
      if(DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>( p_mortar, "ALGORITHM" )
          == INPAR::MORTAR::algorithm_gpts)
        integrator = Teuchos::rcp( new CONTACT::CoIntegratorNitsche(
            p_mortar,slave_type,comm ) );
      else
        integrator = Teuchos::rcp( new CONTACT::CoIntegrator(
                  p_mortar,slave_type,comm ) );
      break;
    }
    case INPAR::CONTACT::solution_lagmult:
    case INPAR::CONTACT::solution_uzawa:
    {
      integrator = Teuchos::rcp(new CONTACT::CoIntegrator(
          p_mortar,slave_type,comm));
      break;
    }
    default:
    {
      dserror( "Unsupported solving strategy! (stype = %s | %d)",
          INPAR::CONTACT::SolvingStrategy2String( sol_type ).c_str(), sol_type );
      exit( EXIT_FAILURE );
    }
  } // end switch

  return integrator;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoIntegrator* CONTACT::INTEGRATOR::Factory::BuildXIntegrator(
    Teuchos::ParameterList& p_mortar,
    const DRT::Element::DiscretizationType& slave_type,
    const Epetra_Comm& comm) const
{
  CONTACT::CoIntegrator* xintegrator = NULL;
  switch (slave_type)
  {
    case DRT::Element::line2:
    {
      xintegrator = BuildConcreteXIntegrator<DRT::Element::line2>(p_mortar,comm);
      break;
    }
    default:
    {
      dserror("Unsupported slave element type! (eletype = %d | %s)",
          slave_type,DRT::DistypeToString(slave_type).c_str());
      break;
    }
  }
  return xintegrator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoIntegrator> CONTACT::INTEGRATOR::BuildIntegrator(
    const INPAR::CONTACT::SolvingStrategy& sol_type,
    Teuchos::ParameterList& p_mortar,
    const DRT::Element::DiscretizationType& slave_type,
    const Epetra_Comm& comm)
{
  Factory factory;
  return factory.BuildIntegrator(sol_type,p_mortar,slave_type,comm);
}
