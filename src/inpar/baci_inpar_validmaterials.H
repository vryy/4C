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
#include "baci_config.H"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

BACI_NAMESPACE_OPEN



namespace INPUT
{
  class MaterialDefinition;

  /// construct list with all materials and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>> ValidMaterials();

  /// print all known material sections without contents
  void PrintEmptyMaterialDefinitions(
      std::ostream& stream, std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>& matlist);
}  // namespace INPUT

/// print empty material sections
void PrintMaterialDatHeader();


BACI_NAMESPACE_CLOSE

#endif
