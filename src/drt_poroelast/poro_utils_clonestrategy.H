/*----------------------------------------------------------------------*/
/*! \file

 \brief strategy to clone porofluid from porous solid

\level 2

 *----------------------------------------------------------------------*/


#ifndef PORO_CLONESTRATEGY_H_
#define PORO_CLONESTRATEGY_H_

#include "../drt_ssi/ssi_clonestrategy.H"

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
      std::map<std::string, std::string> ConditionsToCopy();

     protected:
      //! determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);
      //! set element-specific data (material etc.)
      virtual void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele,
          const int matid, const bool isnurbs);

      //! check for correct material
      void CheckMaterialType(const int matid);

      //! set anisotropic permeability directions onto fluid element
      void SetAnisotropicPermeabilityDirectionsOntoFluid(
          Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele);
    };

    //! \brief implementation of special clone strategy for automatic generation
    //!        of scatra discretization from a given structure discretization for porous media
    class PoroScatraCloneStrategy : public SSI::ScatraStructureCloneStrategy
    {
     public:
      using my = SSI::ScatraStructureCloneStrategy;

      //! constructor
      explicit PoroScatraCloneStrategy() = default;

      //! return SCATRA::ImplType of the element
      INPAR::SCATRA::ImplType GetImplType(
          DRT::Element* ele  //! element whose SCATRA::ImplType shall be determined
          ) override;

     protected:
      //! determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype) override;

      //! provide cloned element with element specific data (material etc.)
      void SetElementData(
          Teuchos::RCP<DRT::Element> newele,  //! current cloned element on target discretization
          DRT::Element* oldele,               //! current element on source discretization
          const int matid,                    //! material of cloned element
          const bool isnurbs                  //! nurbs flag
          ) override;

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy() override;
    };
  }  // namespace UTILS
}  // namespace POROELAST

#endif /* PORO_CLONESTRATEGY_H_ */
