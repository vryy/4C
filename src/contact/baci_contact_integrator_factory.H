/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired integrator object.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_CONTACT_INTEGRATOR_FACTORY_H
#define BACI_CONTACT_INTEGRATOR_FACTORY_H

#include "baci_inpar_contact.H"
#include "baci_lib_element.H"

// forward declaration
class Epetra_Comm;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace CONTACT
{
  class CoIntegrator;
  namespace INTEGRATOR
  {
    class Factory
    {
     public:
      /*! \brief Build the desired contact integrator
       *
       * \date 04/16
       * \author hiermeier */
      Teuchos::RCP<CONTACT::CoIntegrator> BuildIntegrator(
          const INPAR::CONTACT::SolvingStrategy& stype, Teuchos::ParameterList& p_mortar,
          const DRT::Element::DiscretizationType& slave_type, const Epetra_Comm& comm) const;
    };  // class Factory

    // non-member function, please call this one from outside!
    Teuchos::RCP<CONTACT::CoIntegrator> BuildIntegrator(
        const INPAR::CONTACT::SolvingStrategy& stype, Teuchos::ParameterList& p_mortar,
        const DRT::Element::DiscretizationType& eletype, const Epetra_Comm& comm);
  }  // namespace INTEGRATOR
}  // namespace CONTACT


#endif  // CONTACT_INTEGRATOR_FACTORY_H
