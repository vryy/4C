/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid contact constitutive laws for input

\level 1


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definitions */


#ifndef BACI_INPAR_VALIDCONTACTCONSTITUTIVELAW_H
#define BACI_INPAR_VALIDCONTACTCONSTITUTIVELAW_H

/*----------------------------------------------------------------------*/
#include "baci_config.H"

#include <Teuchos_RCP.hpp>

#include <iostream>
#include <vector>

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class LawDefinition;
  }
}  // namespace CONTACT

namespace DRT
{
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
}  // namespace DRT

/// print empty contact constitutivelaw sections
void PrintContactConstitutiveLawDatHeader();


BACI_NAMESPACE_CLOSE

#endif
