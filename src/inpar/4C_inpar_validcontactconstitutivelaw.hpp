/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid contact constitutive laws for input

\level 1


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definitions */


#ifndef FOUR_C_INPAR_VALIDCONTACTCONSTITUTIVELAW_HPP
#define FOUR_C_INPAR_VALIDCONTACTCONSTITUTIVELAW_HPP

/*----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <iostream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class LawDefinition;
  }
}  // namespace CONTACT

namespace INPUT
{
  /// construct list with all contact constitutive laws and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>>
  ValidContactConstitutiveLaws();

  /** \brief print all known contact constitutive law sections without contents
   *
   * \param[in] contactconstitutivlawlist list of contact constitutive law definitions
   */
  void PrintEmptyContactConstitutiveLawDefinitions(std::ostream& stream,
      std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>&
          contactconstitutivlawlist);

}  // namespace INPUT

/// print empty contact constitutivelaw sections
void PrintContactConstitutiveLawDatHeader();


FOUR_C_NAMESPACE_CLOSE

#endif
