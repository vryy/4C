/*----------------------------------------------------------------------*/
/*! \file

\brief Implements a container class containing contact constitutivelaw parameters from the input
file as well as a quick access wrapper for those parameters

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


CONTACT::CONSTITUTIVELAW::Parameter::Parameter(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container>
        coconstlawdata  ///< read and validate contactconstitutivelaw data (of 'slow' access)
    )
    : offset_(*coconstlawdata->Get<double>("Offset")){};
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Container::Container(
    const int id, const INPAR::CONTACT::ConstitutiveLawType type, const std::string name)
    : INPAR::InputParameterContainer(), id_(id), type_(type), name_(name), params_(Teuchos::null)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Container::Print(std::ostream& os) const
{
  os << "ContactConstitutiveLaw " << Id() << " " << Name() << " :: ";

  INPAR::InputParameterContainer::Print(os);
}

FOUR_C_NAMESPACE_CLOSE
