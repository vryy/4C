/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid materials for input

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definitions */


#ifndef DRT_VALIDMATERIALS_H
#define DRT_VALIDMATERIALS_H

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
    void PrintEmptyMaterialDefinitions(std::ostream& stream,
        std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>& condlist, bool color = true);

  }  // namespace INPUT
}  // namespace DRT

/// print empty material sections
void PrintMaterialDatHeader();


#endif
