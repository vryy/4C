// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


CONTACT::CONSTITUTIVELAW::Parameter::Parameter(
    const std::shared_ptr<const CONTACT::CONSTITUTIVELAW::Container>
        coconstlawdata  ///< read and validate contactconstitutivelaw data (of 'slow' access)
    )
    : offset_(coconstlawdata->get<double>("Offset")) {};
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Container::Container(
    const int id, const Inpar::CONTACT::ConstitutiveLawType type, const std::string name)
    : Core::IO::InputParameterContainer(), id_(id), type_(type), name_(name), params_(nullptr)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Container::print(std::ostream& os) const
{
  os << "ContactConstitutiveLaw " << id() << " " << name() << " :: ";

  Core::IO::InputParameterContainer::print(os);
}

FOUR_C_NAMESPACE_CLOSE
