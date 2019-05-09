/*----------------------------------------------------------------------*/
/*!
\brief strategy for cloning scatra discretization from structure discretization

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/
#ifndef SSI_CLONESTRATEGY_H
#define SSI_CLONESTRATEGY_H

#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_scatra.H"

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

  \author schmidt (schmidt@lnm.mw.tum.de)
  \date 09/17
  */

  class ScatraStructureCloneStrategy
  {
   public:
    //! constructor
    explicit ScatraStructureCloneStrategy() { return; };

    //! destructor
    virtual ~ScatraStructureCloneStrategy() { return; };

    //! return map with original names of conditions to be cloned as key values, and final names of
    //! cloned conditions as mapped values
    virtual std::map<std::string, std::string> ConditionsToCopy();

    //! return SCATRA::ImplType of the element
    virtual INPAR::SCATRA::ImplType GetImplType(
        DRT::Element* ele  //! element whose SCATRA::ImplType shall be determined
    );

   protected:
    //! check material of cloned element
    void CheckMaterialType(const int matid  //! material of cloned element
    );

    //! decide whether element should be cloned or not, and if so, determine type of cloned element
    virtual bool DetermineEleType(
        DRT::Element* actele,              //! current element on source discretization
        const bool ismyele,                //! ownership flag
        std::vector<std::string>& eletype  //! vector storing types of cloned elements
    );

    //! provide cloned element with element specific data (material etc.)
    virtual void SetElementData(
        Teuchos::RCP<DRT::Element> newele,  //! current cloned element on target discretization
        DRT::Element* oldele,               //! current element on source discretization
        const int matid,                    //! material of cloned element
        const bool isnurbs                  //! nurbs flag
    );
  };  // class ScatraStructureCloneStrategy
}  // namespace SSI
#endif
