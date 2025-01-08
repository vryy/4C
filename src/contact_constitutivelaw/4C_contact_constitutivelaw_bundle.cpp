// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_constitutivelaw_bundle.hpp"

#include "4C_contact_constitutivelaw_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*/
CONTACT::CONSTITUTIVELAW::Bundle::Bundle() : readfromproblem_(0) {}
/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Bundle::insert(
    int matid, std::shared_ptr<CONTACT::CONSTITUTIVELAW::Container> mat)
{
  map_.insert(std::pair<int, std::shared_ptr<CONTACT::CONSTITUTIVELAW::Container>>(matid, mat));
}

/*----------------------------------------------------------------------*/
int CONTACT::CONSTITUTIVELAW::Bundle::find(const int id) const
{
  if (map_.find(id) == map_.end())
    return -1;
  else
    return map_.find(id)->first;
}

/*----------------------------------------------------------------------*/
void CONTACT::CONSTITUTIVELAW::Bundle::make_parameters()
{
  for (std::map<int, std::shared_ptr<CONTACT::CONSTITUTIVELAW::Container>>::iterator m =
           map_.begin();
      m != map_.end(); ++m)
  {
    int lawid = m->first;

    // 1st try
    {
      // indirectly add quick access parameter members
      std::shared_ptr<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> law =
          CONTACT::CONSTITUTIVELAW::ConstitutiveLaw::factory(lawid);
      // check if allocation was successful
      std::shared_ptr<CONTACT::CONSTITUTIVELAW::Container> lawpar = m->second;
      if (law != nullptr) continue;
    }
    FOUR_C_THROW(
        "Allocation of quick access parameters failed for contact constitutivelaw %d", lawid);
  }
}

/*----------------------------------------------------------------------*/
std::shared_ptr<CONTACT::CONSTITUTIVELAW::Container> CONTACT::CONSTITUTIVELAW::Bundle::by_id(
    const int id) const
{
  std::map<int, std::shared_ptr<CONTACT::CONSTITUTIVELAW::Container>>::const_iterator m =
      map_.find(id);

  if (map_.size() == 0) FOUR_C_THROW("No contact constitutivelaws available, num=%d", id);

  if (m == map_.end())
    FOUR_C_THROW("Contact Constitutive Law 'Law %d' could not be found", id);
  else
    return m->second;

  // catch up
  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
