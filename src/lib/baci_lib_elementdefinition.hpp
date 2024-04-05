/*----------------------------------------------------------------------*/
/*! \file

\brief Central storage of element input line definitions

\level 0


*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_LIB_ELEMENTDEFINITION_HPP
#define FOUR_C_LIB_ELEMENTDEFINITION_HPP

#include "baci_config.hpp"

#include "baci_io_linedefinition.hpp"

#include <Teuchos_RCP.hpp>

#include <iostream>
#include <map>
#include <string>
#include <vector>

BACI_NAMESPACE_OPEN


namespace INPUT
{
  class Lines;

  /// Collection of valid element dat file line definitions
  /*!
    The actual definition is done by each element's type class.
   */
  class ElementDefinition
  {
   public:
    /// Setup of
    void SetupValidElementLines();

    /// print all valid element definitions to stream
    void PrintElementDatHeaderToStream(std::ostream& stream);

    /// Print section header to stream
    void PrintSectionHeader(std::ostream& stream, std::string name);

    /// Print valid input lines for given element type
    void PrintElementLines(std::ostream& stream, std::string name);

    /// return line definitions for given element type
    INPUT::LineDefinition* ElementLines(std::string name, std::string distype);

   private:
    /// input line definitions per element type
    std::map<std::string, std::map<std::string, LineDefinition>> definitions_;
  };

}  // namespace INPUT


void PrintElementDatHeader();


BACI_NAMESPACE_CLOSE

#endif
