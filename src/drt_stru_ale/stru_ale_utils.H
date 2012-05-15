/*----------------------------------------------------------------------*/
/*!
\file stru_ale_utils.H

\brief utility functions for structural ale problems

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*
 | definitions                                               mgit 04/11 |
 *----------------------------------------------------------------------*/

#ifndef STRU_ALE_UTILS_H
#define STRU_ALE_UTILS_H

#include <set>
#include <string>
#include <functional>
#include <iostream>

/*----------------------------------------------------------------------*
 |  headers                                                  mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 |                                                           mgit 04/11 |
 *----------------------------------------------------------------------*/
namespace STRU_ALE
{
namespace UTILS
{
  //! \brief implementation of special clone strategy for automatic generation
  //!        of thermo discretization from a given structure discretization
  class AleStructureCloneStrategy
  {
  public:

    /// constructor
    explicit AleStructureCloneStrategy(){}
    /// destructor
    virtual ~AleStructureCloneStrategy() {}

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

  }; // class AleStructureCloneStrategy

} // namespace UTILS

} // namespace STRU_ALE

/*----------------------------------------------------------------------*/
#endif // TSI_UTILS_H
