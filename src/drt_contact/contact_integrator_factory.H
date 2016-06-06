/*---------------------------------------------------------------------*/
/*!
\file contact_integrator_factory.H

\brief Factory to create the desired integrator object.

\level 2

\maintainer Michael Hiermeier

\date 04/2016

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_
#define SRC_DRT_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_

#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_element.H"

// forward declaration
class Epetra_Comm;
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace CONTACT {
  class CoIntegrator;
namespace INTEGRATOR {
  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    /*! \brief Build the desired contact integrator
     *
     * \date 04/16
     * \author hiermeier */
    Teuchos::RCP<CONTACT::CoIntegrator> BuildIntegrator(
        const INPAR::CONTACT::SolvingStrategy& stype,
        Teuchos::ParameterList& p_mortar,
        const DRT::Element::DiscretizationType& eletype,
        const Epetra_Comm& comm) const;
  }; // class Factory

  // non-member function, please call this one from outside!
  Teuchos::RCP<CONTACT::CoIntegrator> BuildIntegrator(
      const INPAR::CONTACT::SolvingStrategy& stype,
      Teuchos::ParameterList& p_mortar,
      const DRT::Element::DiscretizationType& eletype,
      const Epetra_Comm& comm);
} // namespace INTEGRATOR
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_ */
