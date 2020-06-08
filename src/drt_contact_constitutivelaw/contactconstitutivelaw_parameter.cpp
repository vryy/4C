/*----------------------------------------------------------------------*/
/*! \file

\brief Implements a container class containing contact constitutivelaw parameters from the input
file as well as a quick access wrapper for those parameters

\level 1

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/

#include "contactconstitutivelaw_parameter.H"

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_dserror.H"


CONTACT::CONSTITUTIVELAW::Parameter::Parameter(
    const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container>
        coconstlawdata  ///< read and validate contactconstitutivelaw data (of 'slow' access)
    )
    : offset_(coconstlawdata->GetDouble("Offset")){};
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Container::Container(
    const int id, const INPAR::CONTACT::ConstitutiveLawType type, const std::string name)
    : DRT::Container(), id_(id), type_(type), name_(name), params_(Teuchos::null)
{
  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Container::Print(std::ostream& os) const
{
  os << "ContactConstitutiveLaw " << Id() << " " << Name() << " :: ";

  DRT::Container::Print(os);

  return;
}
