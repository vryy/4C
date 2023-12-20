/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid

\level 2

 *----------------------------------------------------------------------*/


#ifndef BACI_POROELAST_UTILS_CLONESTRATEGY_H
#define BACI_POROELAST_UTILS_CLONESTRATEGY_H



#include "baci_config.H"

#include <Teuchos_RCPDecl.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Element;
}

namespace POROELAST
{
  namespace UTILS
  {
    //! \brief implementation of special clone strategy for automatic generation
    //!        of fluid discretization from a given structure discretization
    class PoroelastCloneStrategy
    {
     public:
      //! constructor
      explicit PoroelastCloneStrategy() = default;
      //! destructor
      virtual ~PoroelastCloneStrategy() = default;

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy() const;


      //! determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

     protected:
      //! set element-specific data (material etc.)
      virtual void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele,
          const int matid, const bool isnurbs);

      //! check for correct material
      void CheckMaterialType(const int matid);

      //! set anisotropic permeability directions onto fluid element
      void SetAnisotropicPermeabilityDirectionsOntoFluid(
          Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele);

      //! set anisotropic permeability nodal coefficients onto fluid element
      void SetAnisotropicPermeabilityNodalCoeffsOntoFluid(
          Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele);
    };

  }  // namespace UTILS
}  // namespace POROELAST

BACI_NAMESPACE_CLOSE

#endif  // POROELAST_UTILS_CLONESTRATEGY_H
