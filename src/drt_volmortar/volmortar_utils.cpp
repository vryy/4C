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
 *----------------------------------------------------------------------*/
void VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterialBToA(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* Aele,
    const std::vector<int>& Bids,
    Teuchos::RCP<DRT::Discretization> disA,
    Teuchos::RCP<DRT::Discretization> disB)
{
  if(Aele==NULL)
    dserror("got NULL pointer for AssignMaterial for element!");

  //if no corresponding element found -> leave
  if(Bids.empty())
    return;

  //default strategy: take only material of first element found
  DRT::Element* Bele = disB->gElement(Bids[0]);

  Aele->AddMaterial(Bele->Material());

  //done
  return;
};

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterialAToB(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* Bele,
    const std::vector<int>& Aids,
    Teuchos::RCP<DRT::Discretization> disA,
    Teuchos::RCP<DRT::Discretization> disB)
{
  if(Bele==NULL)
    dserror("got NULL pointer for AssignMaterial for element!");

  //if no corresponding element found -> leave
  if(Aids.empty())
    return;

  //default strategy: take only material of first element found
  DRT::Element* Aele = disA->gElement(Aids[0]);

  Bele->AddMaterial(Aele->Material());

  //done
  return;
};
