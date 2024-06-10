/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for tsi problems

\level 2
*/

#ifndef FOUR_C_TSI_UTILS_HPP
#define FOUR_C_TSI_UTILS_HPP


#include "4C_config.hpp"

#include "4C_coupling_volmortar_utils.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

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
      bool determine_ele_type(
          Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

      //! set element-specific data (material etc.)
      void set_element_data(Teuchos::RCP<Core::Elements::Element> newele,
          Core::Elements::Element* oldele, const int matid, const bool isnurbs);

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> conditions_to_copy() const;

      //! check for correct material
      void check_material_type(const int matid);

     private:
    };  // class ThermoStructureCloneStrategy

    //! setup TSI, clone the structural discretization
    void SetupTSI(const Epetra_Comm& comm);


    void SetMaterialPointersMatchingGrid(Teuchos::RCP<const Core::FE::Discretization> sourcedis,
        Teuchos::RCP<const Core::FE::Discretization> targetdis);

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
    class TSIMaterialStrategy : public Core::VolMortar::UTILS::DefaultMaterialStrategy
    {
     public:
      //! constructor
      TSIMaterialStrategy(){};

      //! assignment of thermo material to structure material
      void AssignMaterial2To1(const Core::VolMortar::VolMortarCoupl* volmortar,
          Core::Elements::Element* ele1, const std::vector<int>& ids_2,
          Teuchos::RCP<Core::FE::Discretization> dis1,
          Teuchos::RCP<Core::FE::Discretization> dis2) override;

      //! assignment of structure material to thermo material
      void AssignMaterial1To2(const Core::VolMortar::VolMortarCoupl* volmortar,
          Core::Elements::Element* ele2, const std::vector<int>& ids_1,
          Teuchos::RCP<Core::FE::Discretization> dis1,
          Teuchos::RCP<Core::FE::Discretization> dis2) override;
    };

  }  // namespace UTILS

  //! prints the 4C tsi-logo on the screen
  void printlogo();

}  // namespace TSI


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
