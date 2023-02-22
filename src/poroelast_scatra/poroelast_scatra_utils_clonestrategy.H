/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid for scatra elements

\level 2

 *----------------------------------------------------------------------*/


#ifndef POROELASTSCATRA_CLONESTRATEGY_H_
#define POROELASTSCATRA_CLONESTRATEGY_H_

#include "poroelast_utils_clonestrategy.H"

#include <Teuchos_RCP.hpp>
#include "inpar_scatra.H"

namespace DRT
{
  class Element;
}

namespace POROELASTSCATRA
{
  namespace UTILS
  {

    class PoroelastCloneStrategyforScatraElements : public POROELAST::UTILS::PoroelastCloneStrategy
    {
     public:
      //! constructor
      explicit PoroelastCloneStrategyforScatraElements() = default;
      //! destructor
      virtual ~PoroelastCloneStrategyforScatraElements() = default;
      //! determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);
    };

    //! \brief implementation of special clone strategy for automatic generation
    //!        of scatra discretization from a given structure discretization for porous media
    class PoroScatraCloneStrategy
    {
     public:
      //! constructor
      explicit PoroScatraCloneStrategy() = default;
      //! destructor
      virtual ~PoroScatraCloneStrategy() = default;

      //! return SCATRA::ImplType of the element
      INPAR::SCATRA::ImplType GetImplType(
          DRT::Element* ele  //! element whose SCATRA::ImplType shall be determined
      );

     protected:
      //! determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

      //! check for correct material
      void CheckMaterialType(const int matid);

      //! provide cloned element with element specific data (material etc.)
      void SetElementData(
          Teuchos::RCP<DRT::Element> newele,  //! current cloned element on target discretization
          DRT::Element* oldele,               //! current element on source discretization
          const int matid,                    //! material of cloned element
          const bool isnurbs                  //! nurbs flag
      );

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy() const;
    };
  }  // namespace UTILS
}  // namespace POROELASTSCATRA

#endif /* POROELASTSCATRA_CLONESTRATEGY_H_ */
