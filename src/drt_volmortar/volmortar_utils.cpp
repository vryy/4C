/*----------------------------------------------------------------------*/
/*!
 \file volmortar_utils.cpp

 \brief utils functions related to volmortar

\level 2

\maintainer Alexander Popp
 *----------------------------------------------------------------------*/

#include "volmortar_coupling.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "volmortar_utils.H"

/*----------------------------------------------------------------------*
 |  assign material to discretization A                       vuong 09/14|
 *----------------------------------------------------------------------*/
void VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(
    const VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele1, const std::vector<int>& ids_2,
    Teuchos::RCP<DRT::Discretization> dis1, Teuchos::RCP<DRT::Discretization> dis2)
{
  if (ele1 == NULL) dserror("ERROR: Got NULL pointer for AssignMaterial for element!");

  // if no corresponding element found -> leave
  if (ids_2.empty()) return;

  // default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele2 = NULL;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords1 = DRT::UTILS::ElementCenterRefeCoords(ele1);

    for (unsigned i = 0; i < ids_2.size(); ++i)
    {
      DRT::Element* actele2 = dis2->gElement(ids_2[i]);
      std::vector<double> centercoords2 = DRT::UTILS::ElementCenterRefeCoords(actele2);

      LINALG::Matrix<3, 1> diffcoords(true);

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
void VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(
    const VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele2, const std::vector<int>& ids_1,
    Teuchos::RCP<DRT::Discretization> dis1, Teuchos::RCP<DRT::Discretization> dis2)
{
  if (ele2 == NULL) dserror("ERROR: Got NULL pointer for AssignMaterial for element!");

  // if no corresponding element found -> leave
  if (ids_1.empty()) return;

  // default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele1 = NULL;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords2 = DRT::UTILS::ElementCenterRefeCoords(ele2);

    for (unsigned i = 0; i < ids_1.size(); ++i)
    {
      DRT::Element* actele1 = dis1->gElement(ids_1[i]);
      std::vector<double> centercoords1 = DRT::UTILS::ElementCenterRefeCoords(actele1);

      LINALG::Matrix<3, 1> diffcoords(true);

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
