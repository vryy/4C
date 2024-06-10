/*----------------------------------------------------------------------*/
/*! \file

 \brief utils functions related to volmortar

\level 2

 *----------------------------------------------------------------------*/

#include "4C_coupling_volmortar_utils.hpp"

#include "4C_coupling_volmortar.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element_center.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  assign material to discretization A                       vuong 09/14|
 *----------------------------------------------------------------------*/
void Core::VolMortar::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(
    const Core::VolMortar::VolMortarCoupl* volmortar, Core::Elements::Element* ele1,
    const std::vector<int>& ids_2, Teuchos::RCP<Core::FE::Discretization> dis1,
    Teuchos::RCP<Core::FE::Discretization> dis2)
{
  if (ele1 == nullptr) FOUR_C_THROW("ERROR: Got nullptr pointer for AssignMaterial for element!");

  // if no corresponding element found -> leave
  if (ids_2.empty()) return;

  // default strategy: take material of element with closest center in reference coordinates
  Core::Elements::Element* ele2 = nullptr;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords1 = Core::FE::element_center_refe_coords(*ele1);

    for (unsigned i = 0; i < ids_2.size(); ++i)
    {
      Core::Elements::Element* actele2 = dis2->gElement(ids_2[i]);
      std::vector<double> centercoords2 = Core::FE::element_center_refe_coords(*actele2);

      Core::LinAlg::Matrix<3, 1> diffcoords(true);

      for (int j = 0; j < 3; ++j) diffcoords(j, 0) = centercoords1[j] - centercoords2[j];

      if (diffcoords.Norm2() - mindistance < 1e-16)
      {
        mindistance = diffcoords.Norm2();
        ele2 = actele2;
      }
    }
  }

  // assign additional material to element A
  ele1->AddMaterial(ele2->Material());

  // done
  return;
};

/*----------------------------------------------------------------------*
 |  assign material to discretization B                       vuong 09/14|
 *----------------------------------------------------------------------*/
void Core::VolMortar::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(
    const Core::VolMortar::VolMortarCoupl* volmortar, Core::Elements::Element* ele2,
    const std::vector<int>& ids_1, Teuchos::RCP<Core::FE::Discretization> dis1,
    Teuchos::RCP<Core::FE::Discretization> dis2)
{
  if (ele2 == nullptr) FOUR_C_THROW("ERROR: Got nullptr pointer for AssignMaterial for element!");

  // if no corresponding element found -> leave
  if (ids_1.empty()) return;

  // default strategy: take material of element with closest center in reference coordinates
  Core::Elements::Element* ele1 = nullptr;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords2 = Core::FE::element_center_refe_coords(*ele2);

    for (unsigned i = 0; i < ids_1.size(); ++i)
    {
      Core::Elements::Element* actele1 = dis1->gElement(ids_1[i]);
      std::vector<double> centercoords1 = Core::FE::element_center_refe_coords(*actele1);

      Core::LinAlg::Matrix<3, 1> diffcoords(true);

      for (int j = 0; j < 3; ++j) diffcoords(j, 0) = centercoords1[j] - centercoords2[j];

      if (diffcoords.Norm2() - mindistance < 1e-16)
      {
        mindistance = diffcoords.Norm2();
        ele1 = actele1;
      }
    }
  }

  // assign additional material to element B
  ele2->AddMaterial(ele1->Material());

  // done
  return;
};

FOUR_C_NAMESPACE_CLOSE
