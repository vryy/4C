#ifndef FOUR_C_INPAR_VALIDPARAMETERS_HPP
#define FOUR_C_INPAR_VALIDPARAMETERS_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

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
  /**
   * Construct a `Teuchos::ParameterList` with all parameters and their documentation.
   *
   * @return A `Teuchos::RCP` to a constant `Teuchos::ParameterList` containing all valid
   * parameters.
   */
  Teuchos::RCP<const Teuchos::ParameterList> valid_parameters();

  /**
   * Print the values of parameter entries from a dat file based on the provided parameter list.
   *
   * This function prints the header of the dat file, including parameter names and their values.
   * It processes sublists and parameters, formatting them appropriately.
   *
   * @param stream      The output stream to which the dat header will be printed.
   * @param list        The parameter list containing the parameters and sublists to be printed.
   * @param parentname  The name of the parent section to which the current section belongs
   *                    (default is an empty string).
   * @param comment     A flag indicating whether to print comments (default is true).
   */
  void print_dat_header(std::ostream& stream, const Teuchos::ParameterList& list,
      std::string parentname = "", bool comment = true);

  /**
   * Print the documentation for a parameter entry.
   *
   * This helper function prints the documentation string associated with a parameter entry.
   * It is used internally by `print_dat_header` to include comments in the output.
   *
   * @param stream The output stream to which the documentation will be printed.
   * @param entry The parameter entry whose documentation will be printed.
   */
  void print_documentation(std::ostream& stream, const Teuchos::ParameterEntry& entry);

  /**
   * Print a sublist from the parameter list.
   *
   * This helper function prints the section header for a sublist and recursively
   * calls `print_dat_header` to print the contents of the sublist.
   * It is used internally by `print_dat_header` to handle nested lists.
   *
   * @param stream The output stream to which the sublist will be printed.
   * @param parentname The name of the parent section to which the sublist belongs.
   * @param name The name of the sublist to be printed.
   * @param list The parameter list containing the sublist to be printed.
   * @param comment A flag indicating whether to print comments (default is true).
   */
  void print_sublist(std::ostream& stream, const std::string& parentname, const std::string& name,
      const Teuchos::ParameterList& list, bool comment);

  /**
   * Print the details of a parameter entry.
   *
   * This helper function prints the name and value of a parameter entry. If the parameter
   * value is not printable, it catches the exception and prints the demangled type name instead.
   * It is used internally by `print_dat_header`.
   *
   * @param stream The output stream to which the parameter details will be printed.
   * @param entry The parameter entry whose details will be printed.
   * @param name The name of the parameter.
   * @param list The parameter list containing the parameter entry.
   * @param comment A flag indicating whether to print comments (default is true).
   */
  void print_parameter(std::ostream& stream, const Teuchos::ParameterEntry& entry,
      const std::string& name, const Teuchos::ParameterList& list, bool comment);

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
