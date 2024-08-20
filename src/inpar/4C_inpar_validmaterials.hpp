/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid materials for input

\level 1


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definitions */


#ifndef FOUR_C_INPAR_VALIDMATERIALS_HPP
#define FOUR_C_INPAR_VALIDMATERIALS_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class MaterialDefinition;
}

namespace Input
{
  /// construct list with all materials and documentation
  Teuchos::RCP<std::vector<Teuchos::RCP<Mat::MaterialDefinition>>> valid_materials();

  /// print all known material sections without contents
  void print_empty_material_definitions(
      std::ostream& stream, std::vector<Teuchos::RCP<Mat::MaterialDefinition>>& matlist);
}  // namespace Input

/// print empty material sections
void print_material_dat_header();


FOUR_C_NAMESPACE_CLOSE

#endif
