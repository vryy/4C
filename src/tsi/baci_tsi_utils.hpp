/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for tsi problems

\level 2
*/

#ifndef FOUR_C_TSI_UTILS_HPP
#define FOUR_C_TSI_UTILS_HPP


#include "baci_config.hpp"

#include "baci_coupling_volmortar_utils.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Element;
  class Discretization;
}  // namespace DRT

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace TSI
{
  namespace UTILS
  {
    //! \brief implementation of special clone strategy for automatic generation
    //!        of thermo discretization from a given structure discretization
    class ThermoStructureCloneStrategy
    {
     public:
      //! constructor
      explicit ThermoStructureCloneStrategy() {}
      //! destructor
      virtual ~ThermoStructureCloneStrategy() = default;

     protected:
      //! determine element type std::string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

      //! set element-specific data (material etc.)
      void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid,
          const bool isnurbs);

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy() const;

      //! check for correct material
      void CheckMaterialType(const int matid);

     private:
    };  // class ThermoStructureCloneStrategy

    //! setup TSI, clone the structural discretization
    void SetupTSI(const Epetra_Comm& comm);


    void SetMaterialPointersMatchingGrid(Teuchos::RCP<const DRT::Discretization> sourcedis,
        Teuchos::RCP<const DRT::Discretization> targetdis);

    //! strategy for material assignment for non matching meshes with TSI

    /// Helper class for assigning materials for volumetric coupling of non conforming meshes (TSI)
    /*!
     When coupling two overlapping discretizations, most often one discretization needs access
     to the corresponding element/material on the other side. For conforming meshes this is straight
     forward as there is one unique element on the other side and therefore one unique material,
     which can be accessed. However, for non conforming meshes there are potentially several
     elements overlapping. Therefore, some rule for assigning materials is needed. This class is
     meant to do that. It gets the element to which it shall assign a material and a vector of IDs
     of the overlapping elements of the other discretization.

     In case of TSI, we also need the kinematic type of the structural element to be known in the
     thermo element, which is why there is a special strategy for TSI. Note that this is not yet
     working for inhomogeneous material properties.

     \author vuong 10/14
     */
    class TSIMaterialStrategy : public CORE::VOLMORTAR::UTILS::DefaultMaterialStrategy
    {
     public:
      //! constructor
      TSIMaterialStrategy(){};

      //! assignment of thermo material to structure material
      void AssignMaterial2To1(const CORE::VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele1,
          const std::vector<int>& ids_2, Teuchos::RCP<DRT::Discretization> dis1,
          Teuchos::RCP<DRT::Discretization> dis2) override;

      //! assignment of structure material to thermo material
      void AssignMaterial1To2(const CORE::VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele2,
          const std::vector<int>& ids_1, Teuchos::RCP<DRT::Discretization> dis1,
          Teuchos::RCP<DRT::Discretization> dis2) override;
    };

  }  // namespace UTILS

  //! prints the BACI tsi-logo on the screen
  void printlogo();

}  // namespace TSI


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // TSI_UTILS_H
