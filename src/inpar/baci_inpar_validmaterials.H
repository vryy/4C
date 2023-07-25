/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid materials for input

\level 1


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definitions */


#ifndef BACI_INPAR_VALIDMATERIALS_H
#define BACI_INPAR_VALIDMATERIALS_H

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>
#include <string>
#include <iostream>

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>


namespace DRT
{
  namespace INPUT
  {
    class MaterialDefinition;

    /// construct list with all materials and documentation
    Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>> ValidMaterials();

    /// print all known material sections without contents
    void PrintEmptyMaterialDefinitions(
        std::ostream& stream, std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& condlist);
  }  // namespace INPUT
}  // namespace DRT

/// print empty material sections
void PrintMaterialDatHeader();


#endif
