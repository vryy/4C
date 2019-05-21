/*----------------------------------------------------------------------*/
/*!

\brief mesh clone strategy for scalar transport problems

\level 1

\maintainer Anh-Tu Vuong
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_UTILS_CLONESTRATEGY_H
#define SCATRA_UTILS_CLONESTRATEGY_H

#include <Teuchos_RCP.hpp>


namespace DRT
{
  class Element;
}

namespace SCATRA
{
  /*!
  \brief implementation of special clone strategy for automatic generation
         of scatra from a given fluid discretization

   */
  class ScatraFluidCloneStrategy
  {
   public:
    /// constructor
    explicit ScatraFluidCloneStrategy() {}
    /// destructor
    virtual ~ScatraFluidCloneStrategy() {}

    /// returns conditions names to be copied (source and target name)
    virtual std::map<std::string, std::string> ConditionsToCopy();

   protected:
    /// determine element type std::string and whether element is copied or not
    virtual bool DetermineEleType(
        DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

    /// set element-specific data (material etc.)
    void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid,
        const bool isnurbs);

    /// check for correct material
    void CheckMaterialType(const int matid);

   private:
  };  // class ScatraFluidCloneStrategy

  /*!
  \brief implementation of special clone strategy for automatic generation
         of scatra from a given fluid discretization

   */
  class ScatraReactionCloneStrategy
  {
   public:
    /// constructor
    explicit ScatraReactionCloneStrategy() {}
    /// destructor
    virtual ~ScatraReactionCloneStrategy() {}

    /// returns conditions names to be copied (source and target name)
    virtual std::map<std::string, std::string> ConditionsToCopy();

   protected:
    /// determine element type std::string and whether element is copied or not
    virtual bool DetermineEleType(
        DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

    /// set element-specific data (material etc.)
    void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid,
        const bool isnurbs);

    /// check for correct material
    void CheckMaterialType(const int matid);

   private:
  };  // class ScatraFluidCloneStrategy

}  // namespace SCATRA

#endif
