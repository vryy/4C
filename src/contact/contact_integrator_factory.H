/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired integrator object.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef SRC_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_
#define SRC_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_

#include "inpar_contact.H"
#include "element.H"

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


#endif /* SRC_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_ */
