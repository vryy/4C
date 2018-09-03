/*----------------------------------------------------------------------*/
/*!
\file drt_elementdefinition.H

\brief Central storage of element input line definitions

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_ELEMENTDEFINITION_H
#define DRT_ELEMENTDEFINITION_H

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
      void PrintSectionHeader(std::ostream& stream, std::string name, bool color = false);

      /// Print valid input lines for given element type
      void PrintElementLines(std::ostream& stream, std::string name);

      /// return line definitions for given element type
      DRT::INPUT::LineDefinition* ElementLines(std::string name, std::string distype);

     private:
      /// input line definitions per element type
      std::map<std::string, std::map<std::string, LineDefinition>> definitions_;
    };

  }  // namespace INPUT
}  // namespace DRT


void PrintElementDatHeader();


#endif
