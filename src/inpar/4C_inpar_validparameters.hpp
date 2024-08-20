/*----------------------------------------------------------------------*/
/*! \file
\brief Setup of the list of valid input parameters

\level 1

*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_INPAR_VALIDPARAMETERS_HPP
#define FOUR_C_INPAR_VALIDPARAMETERS_HPP

#include "4C_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class Pstream;
}


namespace Input
{
  /// construct list with all parameters and documentation
  Teuchos::RCP<const Teuchos::ParameterList> valid_parameters();

  /// print flag sections of dat file with given list
  void print_dat_header(std::ostream& stream, const Teuchos::ParameterList& list,
      std::string parentname = "", bool comment = true);

  /**
   * Return true if the @p list contains any parameter that has whitespace in the key name.
   *
   * @note This is needed for the NOX parameters whose keywords and value have white spaces and
   * thus '=' are inserted to distinguish them.
   */
  bool need_to_print_equal_sign(const Teuchos::ParameterList& list);

}  // namespace Input


/*! print list of valid parameters with documentation */
void print_valid_parameters();

/*! print help message */
void print_help_message();

/*! print flag sections of dat file with default flags */
void print_default_dat_header();


FOUR_C_NAMESPACE_CLOSE

#endif
