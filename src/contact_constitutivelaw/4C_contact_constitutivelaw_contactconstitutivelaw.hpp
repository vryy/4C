// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_CONSTITUTIVELAW_CONTACTCONSTITUTIVELAW_HPP
#define FOUR_C_CONTACT_CONSTITUTIVELAW_CONTACTCONSTITUTIVELAW_HPP



#include "4C_config.hpp"

#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"
#include "4C_contact_node.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class Parameter;

    /**
     * \brief The ConstitutiveLaw class provides a framework in order to relate the contact gap to
     * the contact pressure using information like micro roughness for contact problems.
     */
    class ConstitutiveLaw
    {
     public:
      /// return type of this constitutive law
      virtual Inpar::CONTACT::ConstitutiveLawType get_constitutive_law_type() const = 0;

      /// Return quick accessible Contact Constitutive Law parameter data
      virtual CONTACT::CONSTITUTIVELAW::Parameter* parameter() const = 0;

      virtual double evaluate(double gap, CONTACT::Node* cnode) = 0;
      virtual double evaluate_deriv(double gap, CONTACT::Node* cnode) = 0;

      /* \brief create Contact ConstitutiveLaw object given the id of the constitutive law in the
       * input file
       */
      static Teuchos::RCP<ConstitutiveLaw> factory(const int id);

      /* \brief create Contact Constitutivelaw object given input information for the constitutive
       * law \param[in] Container holding the Coefficients for the Contact ConstitutiveLaw
       */
      static Teuchos::RCP<ConstitutiveLaw> factory(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> contactconstitutivelawdata);

      virtual ~ConstitutiveLaw() = default;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
