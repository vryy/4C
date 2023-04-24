/*----------------------------------------------------------------------*/
/*! \file
\brief strategy for cloning scatra discretization from structure discretization

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef SSI_CLONESTRATEGY_H
#define SSI_CLONESTRATEGY_H

#include <Teuchos_RCP.hpp>
#include "inpar_scatra.H"

// forward declaration
namespace DRT
{
  class Element;
}

namespace SSI
{
  /*!
  \brief strategy for cloning scatra discretization from structure discretization

  For some scatra-structure interaction problems, the scatra discretization is obtained through
  cloning from the structure discretization.

  \date 09/17
  */

  class ScatraStructureCloneStrategy
  {
   public:
    explicit ScatraStructureCloneStrategy() = default;

    virtual ~ScatraStructureCloneStrategy() = default;

    //! return map with original names of conditions to be cloned as key values, and final names of
    //! cloned conditions as mapped values
    virtual std::map<std::string, std::string> ConditionsToCopy() const;

    //! get impltype of scatra element from structure element
    //!
    //! \param ele     element whose SCATRA::ImplType shall be determined
    //! \return        impltype of the scatra element
    virtual INPAR::SCATRA::ImplType GetImplType(DRT::Element* ele);

   protected:
    //! check material of cloned element
    //!
    //! \param matid     material of cloned element
    void CheckMaterialType(const int matid);

    //! decide whether element should be cloned or not, and if so, determine type of cloned element
    //!
    //! \param actele       current element on source discretization
    //! \param ismyele      ownership flag
    //! \param eletype      vector storing types of cloned elements
    //! \return
    virtual bool DetermineEleType(
        DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

    //! provide cloned element with element specific data (material etc.)
    //!
    //! \param newele    current cloned element on target discretization
    //! \param oldele    current element on source discretization
    //! \param matid     material of cloned element
    //! \param isnurbs   nurbs flag
    virtual void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele,
        const int matid, const bool isnurbs);
  };

  class ScatraStructureCloneStrategyManifold : public ScatraStructureCloneStrategy
  {
   public:
    std::map<std::string, std::string> ConditionsToCopy() const override;

    void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid,
        const bool isnurbs) override;
  };

}  // namespace SSI
#endif
