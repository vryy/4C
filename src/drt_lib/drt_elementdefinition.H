#ifdef CCADISCRET

#ifndef DRT_ELEMENTDEFINITION_H
#define DRT_ELEMENTDEFINITION_H

#ifdef __cplusplus

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>

#include "drt_linedefinition.H"

namespace DRT
{
namespace INPUT
{

  class Lines;

  /// Collection of valid element dat file line definitions
  class ElementDefinition
  {
  public:

    /// Setup of
    void SetupValidElementLines();

    /// Print section header to stream
    void PrintSectionHeader(std::ostream& stream, std::string name, bool color=false);

    /// Print valid input lines for given element type
    void PrintElementLines(std::ostream& stream, std::string name);

    /// return line definitions for given element type
    DRT::INPUT::LineDefinition* ElementLines(std::string name, std::string distype);

  private:

    /// @name Elementspecific element lines

    void SetupFluid3Lines();
    void SetupSolidH8Lines();

    //@}

    /// input line definitions per element type
    std::map<std::string,std::map<std::string,LineDefinition> > definitions_;
  };

}
}

#endif

#ifdef __cplusplus
extern "C"
{
#endif

  void PrintElementDatHeader();

#ifdef __cplusplus
}
#endif

#endif
#endif
