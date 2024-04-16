/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired integrator object.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_INTEGRATOR_FACTORY_HPP
#define FOUR_C_CONTACT_INTEGRATOR_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_inpar_contact.hpp"
#include "baci_lib_element.hpp"

// forward declaration
class Epetra_Comm;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  class Integrator;
  namespace INTEGRATOR
  {
    class Factory
    {
     public:
      /*! \brief Build the desired contact integrator
       *
       * \date 04/16
       * \author hiermeier */
      Teuchos::RCP<CONTACT::Integrator> BuildIntegrator(
          const INPAR::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
          const CORE::FE::CellType& slave_type, const Epetra_Comm& comm) const;
    };  // class Factory

    // non-member function, please call this one from outside!
    Teuchos::RCP<CONTACT::Integrator> BuildIntegrator(
        const INPAR::CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
        const CORE::FE::CellType& slave_type, const Epetra_Comm& comm);
  }  // namespace INTEGRATOR
}  // namespace CONTACT


BACI_NAMESPACE_CLOSE

#endif
