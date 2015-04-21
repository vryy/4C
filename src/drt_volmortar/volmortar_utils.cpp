/*----------------------------------------------------------------------*/
/*!
 \file volmortar_utils.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "volmortar_coupling.H"

#include "../drt_lib/drt_discret.H"

#include "volmortar_utils.H"

/*----------------------------------------------------------------------*
 |  assign material to discretization A                       vuong 09/14|
 *----------------------------------------------------------------------*/
void VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* ele1,
    const std::vector<int>& ids_2,
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  if(ele1==NULL)
    dserror("ERROR: Got NULL pointer for AssignMaterial for element!");

  //if no corresponding element found -> leave
  if(ids_2.empty())
    return;

  //default strategy: take only material of first element found
  DRT::Element* ele2 = dis2->gElement(ids_2[0]);

  //assign additional material to element A
  ele1->AddMaterial(ele2->Material());

  //done
  return;
};

/*----------------------------------------------------------------------*
 |  assign material to discretization B                       vuong 09/14|
 *----------------------------------------------------------------------*/
void VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* ele2,
    const std::vector<int>& ids_1,
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  if(ele2==NULL)
    dserror("ERROR: Got NULL pointer for AssignMaterial for element!");

  //if no corresponding element found -> leave
  if(ids_1.empty())
    return;

  //default strategy: take only material of first element found
  DRT::Element* ele1 = dis1->gElement(ids_1[0]);

  //assign additional material to element B
  ele2->AddMaterial(ele1->Material());

  //done
  return;
};
