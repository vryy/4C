/*----------------------------------------------------------------------*/
/*!
\file tsi_utils.H

\brief utility functions for tsi problems

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 11/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TSI_UTILS_H
#define TSI_UTILS_H

#include <set>
#include <string>
#include <functional>
#include <iostream>

/*----------------------------------------------------------------------*
 |  headers                                                  dano 11/09 |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace TSI
{
/// TSI::UTILS: Random stuff that might be helpful when dealing with TSI problems
namespace UTILS
{
  //! \brief implementation of special clone strategy for automatic generation
  //!        of thermo discretization from a given structure discretization
  class ThermoStructureCloneStrategy
  {
  public:

    /// constructor
    explicit ThermoStructureCloneStrategy(){}
    /// destructor
    virtual ~ThermoStructureCloneStrategy() {}

  protected:

    /// determine element type string and whether element is copied or not
    bool DetermineEleType(
      DRT::Element* actele,
      const bool ismyele,
      vector<string>& eletype
      );

    /// set element-specific data (material etc.)
    void SetElementData(
      Teuchos::RCP<DRT::Element> newele,
      DRT::Element* oldele,
      const int matid,
      const bool isnurbs
      );

    /// returns conditions names to be copied (source and target name)
    map<string,string> ConditionsToCopy();

    /// check for correct material
    void CheckMaterialType(const int matid);

  private:

  }; // class ThermoStructureCloneStrategy

} // namespace UTILS

} // namespace TSI

/*----------------------------------------------------------------------*/
#endif // TSI_UTILS_H
#endif // CCADISCRET
