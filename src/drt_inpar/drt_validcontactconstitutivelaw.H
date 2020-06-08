/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid contact constitutive laws for input

\level 1

\maintainer Nora Hagmeyer

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definitions */


#ifndef DRT_VALIDCONSTITUTIVELAWS_H
#define DRT_VALIDCONSTITUTIVELAWS_H

/*----------------------------------------------------------------------*/
#include <vector>
#include <iostream>

#include <Teuchos_RCP.hpp>

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
            contactconstitutivlawlist,
        bool color = true);

  }  // namespace INPUT
}  // namespace DRT

/// print empty contact constitutivelaw sections
void PrintContactConstitutiveLawDatHeader();


#endif
