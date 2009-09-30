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

    /// print all valid element defintions to stream
    void PrintElementDatHeaderToStream(std::ostream& stream);

    /// Print section header to stream
    void PrintSectionHeader(std::ostream& stream, std::string name, bool color=false);

    /// Print valid input lines for given element type
    void PrintElementLines(std::ostream& stream, std::string name);

    /// return line definitions for given element type
    DRT::INPUT::LineDefinition* ElementLines(std::string name, std::string distype);

  private:

    /// @name Elementspecific element lines

    void SetupArtLines();
    void SetupBeam2Lines();
    void SetupBeam2rLines();
    void SetupBeam3Lines();
    void SetupConstrele2Lines();
    void SetupConstrele3Lines();
    void SetupPtet4Lines();
    void SetupShell8Lines();
    void SetupSolid3Lines();
    void SetupSolidh20Lines();
    void SetupSolidh27Lines();
    void SetupSolidh8Lines();
    void SetupSolidh8p1j1Lines();
    void SetupSolidsh8Lines();
    void SetupSolidsh8p8Lines();
    void SetupSolidshw6Lines();
    void SetupSolidt10Lines();
    void SetupSolidt4Lines();
    void SetupSolidw6Lines();
    void SetupTorsion2Lines();
    void SetupTorsion3Lines();
    void SetupTruss2Lines();
    void SetupTruss3Lines();
    void SetupWallLines();

    void SetupCombust3Lines();
    void SetupCondif2Lines();
    void SetupCondif3Lines();
    void SetupFluid2Lines();
    void SetupFluid3Lines();
    void SetupTranspLines();
    void SetupXdiff3Lines();
    void SetupXfluid3Lines();

    void SetupAle2Lines();
    void SetupAle3Lines();

    void SetupThermoLines();

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
