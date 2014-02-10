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
#ifndef TSI_UTILS_H
#define TSI_UTILS_H


/*----------------------------------------------------------------------*
 | header                                                    dano 11/09 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Element;
  class Discretization;
}

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace TSI
{
//! TSI::UTILS: Random stuff that might be helpful when dealing with TSI problems
namespace UTILS
{
  //! \brief implementation of special clone strategy for automatic generation
  //!        of thermo discretization from a given structure discretization
  class ThermoStructureCloneStrategy
  {
  public:

    //! constructor
    explicit ThermoStructureCloneStrategy(){}
    //! destructor
    virtual ~ThermoStructureCloneStrategy() {}

  protected:

    //! determine element type std::string and whether element is copied or not
    bool DetermineEleType(
      DRT::Element* actele,
      const bool ismyele,
      std::vector<std::string>& eletype
      );

    //! set element-specific data (material etc.)
    void SetElementData(
      Teuchos::RCP<DRT::Element> newele,
      DRT::Element* oldele,
      const int matid,
      const bool isnurbs
      );

    //! returns conditions names to be copied (source and target name)
    std::map<std::string,std::string> ConditionsToCopy();

    //! check for correct material
    void CheckMaterialType(const int matid);

  private:

  }; // class ThermoStructureCloneStrategy

  //! setup TSI, clone the structural discretization
  void SetupTSI(const Epetra_Comm& comm);


  void SetMaterialPointersMatchingGrid(
      Teuchos::RCP<const DRT::Discretization> sourcedis,
      Teuchos::RCP<const DRT::Discretization> targetdis);

} // namespace UTILS

  //! prints the BACI tsi-logo on the screen
  void printlogo();

} // namespace TSI


/*----------------------------------------------------------------------*/
#endif // TSI_UTILS_H
