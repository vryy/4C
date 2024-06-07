/*----------------------------------------------------------------------*/
/*! \file

\brief Central storage of element input line definitions

\level 0


*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FEM_GENERAL_ELEMENT_DEFINITION_HPP
#define FOUR_C_FEM_GENERAL_ELEMENT_DEFINITION_HPP

#include "4C_config.hpp"

#include "4C_io_linedefinition.hpp"

#include <Teuchos_RCP.hpp>

#include <iostream>
#include <map>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Core::Elements
{
  /// Collection of valid element dat file line definitions
  /*!
    The actual definition is done by each element's type class.
   */
  class ElementDefinition
  {
   public:
    /// Setup of
    void setup_valid_element_lines();

    /// print all valid element definitions to stream
    void print_element_dat_header_to_stream(std::ostream& stream);

    /// Print section header to stream
    void PrintSectionHeader(std::ostream& stream, std::string name);

    /// Print valid input lines for given element type
    void PrintElementLines(std::ostream& stream, std::string name);

    /// return line definitions for given element type
    Input::LineDefinition* ElementLines(std::string name, std::string distype);

   private:
    /// input line definitions per element type
    std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_;
  };

}  // namespace Core::Elements


void PrintElementDatHeader();


FOUR_C_NAMESPACE_CLOSE

#endif
